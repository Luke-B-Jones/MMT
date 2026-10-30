#!/usr/bin/env python3
"""MATAFILER4 / GeneCat layout for MMT.

  index   map + MATAFILER4 run + GeneCat MAGvsGC + bin files + GeneCat representative FASTAs -> tables in --out
  assign  one sample's CRAM (samtools view -h text on stdin) -> which read lies on which MAG / MGS (tsv) + read depth
          of every MGS-MAG contig in 1 kb bins (tsv)
  split   a by-strand BAM -> one unaligned BAM per group (MGS or MAG: 2nd column of --reads)

Names and rules follow GeneCat's own code (Mods/Binning.pm getRepresentBins + createBinCtgs, Mods/GenoMetaAss.pm
readMapS), so MMT finds exactly the contigs GeneCat used:
  MAG name             "<sample or alias>__<bin>", split at the LAST "__"
  co-assembly alias    "<SmplID>M<n>" = the n-th sample of an AssmblGrps group in map order -> SmplID
  bins of a sample     <assembly dir>/Binning/<binner>/<SmplID>   (lines: contig TAB bin)
  assembly dir         first line of <MF4 run>/<SmplID>/assemblies/metag/assembly.txt (shared by a co-assembly group)
  representative       starred MAG of the MGS, else its first non-canopy MAG (MAGvsGC row order; rows need not be grouped)
  its FASTA            <GeneCat>/Bin_<binner>/Genomes/MGS_ctg/<MGS>.ctgs.<MAG>, ".gz" stripped, ".fna" added unless
                       it ends in .fa/.fna/.fasta, then ".gz"
Every representative's contigs from the bin file must equal GeneCat's FASTA, else the index stops (the rules would be
wrong). units.tsv = one row per assembly x MGS with a MAG there, with checksums of the contig lists, so a result is
redone only when the contigs behind it change (not when unrelated lines of the map or MAGvsGC change). Reads: MATAFILER4 keeps every aligned read in <MF4 run>/<SmplID>/mapping/<SmplID>-smd.cram; alignments that
failed its filter are flagged 0x200 (QC fail), duplicates 0x400.
movies.tsv = the read groups of every raw BAM (movie, chemistry): the batch key for kinetics.
abundance.tsv = sample x MGS pairs with abundance > 0 in GeneCat's MGS.matL0.txt (presence, also where no MAG was binned).
"""
import argparse, collections, glob, gzip, hashlib, os, re, sys

NAME_RE = re.compile(r"^(.+)__(.+)$")    # greedy first group = split at the last "__" (as GeneCat)


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def stop(msg):
    sys.exit("ERROR: " + msg)


# ============================================================================ index
def read_map(path):
    """-> samples [(SmplID, raw_bam, problem)], alias -> SmplID"""
    mapdir = os.path.dirname(os.path.abspath(path))
    header, dirpath, samples, alias, group_n = None, "", [], {}, collections.Counter()
    for raw in open(path):
        line = raw.rstrip("\r\n")
        if not line.strip():
            continue
        if line.startswith("#"):
            m = re.match(r"#DirPath\s+(\S+)", line)
            if m:
                dirpath = m.group(1) if os.path.isabs(m.group(1)) else os.path.join(mapdir, m.group(1))
            elif line[1:].split("\t")[0] == "SmplID":
                header = line[1:].split("\t")
            continue
        if header is None:
            stop(f"map {path}: no '#SmplID' header line before the first sample")
        f = line.split("\t")
        row = dict(zip(header, f + [""] * (len(header) - len(f))))
        s = row["SmplID"].strip()
        grp = row.get("AssmblGrps", "").strip()
        if grp:                                   # GenoMetaAss: altCurSmp = SmplID . "M" . count within the group
            group_n[grp] += 1
            alias[f"{s}M{group_n[grp]}"] = s
        raw_bam, problem = "", ""
        prefix, rpath = row.get("SmplPrefix", "").strip(), row.get("Path", "").strip()
        if prefix:
            base = os.path.join(dirpath, prefix)
            cands = [base] if base.endswith(".bam") else [base + ".bam"]
            hits = [c for c in cands if os.path.isfile(c)] or sorted(glob.glob(base + "*.bam"))
        elif rpath:
            hits = sorted(glob.glob(os.path.join(dirpath, rpath, "*.bam")))
            hits = [h for h in hits if h.endswith(".hifi_reads.bam")] or hits
        else:
            hits = []
        if len(hits) == 1:
            raw_bam = hits[0]
        else:
            problem = "no_raw_bam" if not hits else "several_raw_bams"
        samples.append((s, raw_bam, problem))
    return samples, alias


def first_line(path):
    with open(path) as fh:
        return fh.readline().strip()


MOVIE_COLS = ("rg_id", "movie", "platform", "read_type", "binding_kit", "sequencing_kit", "basecaller", "frame_rate")


def movies(raw):
    """read groups of a raw BAM header: [(rg_id, movie, platform, read_type, binding_kit, sequencing_kit, basecaller,
    frame_rate)] ; None when the header cannot be read"""
    try:
        import pysam
        with pysam.AlignmentFile(raw, check_sq=False) as f:
            rgs = f.header.to_dict().get("RG", [])
    except Exception:
        return None
    out = []
    for rg in rgs:
        ds = dict(kv.split("=", 1) for kv in str(rg.get("DS", "")).split(";") if "=" in kv)
        out.append((rg.get("ID", ""), rg.get("PU", ""), rg.get("PM", ""), ds.get("READTYPE", ""), ds.get("BINDINGKIT", ""),
                    ds.get("SEQUENCINGKIT", ""), ds.get("BASECALLERVERSION", ""), ds.get("FRAMERATEHZ", "")))
    return out


def rep_fasta_name(mgs, mag):   # Binning.pm createBinCtgs
    name = re.sub(r"\.gz$", "", f"{mgs}.ctgs.{mag}", flags=re.I)
    if not re.search(r"\.(fa|fna|fasta)$", name, flags=re.I):
        name += ".fna"
    return name + ".gz"


def fasta_lengths(path):
    lens, name = collections.OrderedDict(), None
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lens[name] = 0
            elif name is not None:
                lens[name] += len(line.strip())
    return lens


def read_abundance(path, samples, mgs_ids, report):
    """GeneCat's MGS abundance matrix -> [(sample, MGS, value)] for values > 0. Rows MGS x columns samples or the
    transpose (whichever matches more names) ; a header one field short (R row names) is padded."""
    if not path or not os.path.isfile(path):
        report(f"GeneCat abundance: not found ({path}) - presence comes from MAGs only")
        return []
    rows = [l.rstrip("\r\n").split("\t") for l in opener(path) if l.strip()]
    if len(rows) < 2:
        report(f"GeneCat abundance: {path} is empty"); return []
    hd, body = rows[0], rows[1:]
    if len(hd) == len(body[0]) - 1:
        hd = [""] + hd
    S, M = set(samples), set(mgs_ids)
    ms = sum(h in S for h in hd[1:]) * sum(r[0] in M for r in body)        # rows MGS, columns samples
    sm = sum(h in M for h in hd[1:]) * sum(r[0] in S for r in body)        # rows samples, columns MGS
    if not ms and not sm:
        report(f"GeneCat abundance: {path}: no row / column names match the MGS and samples - not used")
        return []
    out = []
    for r in body:
        for h, v in zip(hd[1:], r[1:]):
            smp, mgs = (h, r[0]) if ms >= sm else (r[0], h)
            if smp in S and mgs in M:
                try:
                    if float(v) > 0: out.append((smp, mgs, v))
                except ValueError:
                    pass
    report(f"GeneCat abundance: {path} ({'MGS x samples' if ms >= sm else 'samples x MGS'}): "
           f"{len({m for _, m, _ in out})} MGS present in {len({s for s, _, _ in out})} samples ({len(out)} pairs > 0)")
    return out


def index(a):
    os.makedirs(a.out, exist_ok=True)
    rep_lines, fatal, seen = [], [], []     # seen: files read whose change rebuilds the index (assembly.txt, bin files)
    report = lambda s: rep_lines.append(s)

    # ---- samples
    maprows, alias = read_map(a.map)
    sample = {}
    for s, raw, problem in maprows:
        sd = os.path.join(a.mf4_run, s)
        cram = os.path.join(sd, "mapping", f"{s}-smd.cram")
        at = os.path.join(sd, "assemblies", "metag", "assembly.txt")
        seen.append(at)
        asm = os.path.normpath(first_line(at)) if os.path.isfile(at) and first_line(at) else ""
        why = [problem] if problem else []
        if not os.path.isfile(cram):
            why.append("no_cram")
        if not asm or not os.path.isdir(asm):
            why.append("no_assembly")
        sample[s] = dict(raw=raw, cram=cram, asm=asm, status=",".join(why) or "ok")
    n_ok = sum(v["status"] == "ok" for v in sample.values())
    report(f"samples in map: {len(sample)}   usable: {n_ok}   co-assembly aliases: {len(alias)}")
    for s, v in sample.items():
        if v["status"] != "ok":
            report(f"  sample {s}: {v['status']}")
    if n_ok == 0:
        fatal.append("no usable sample (raw BAM + CRAM + assembly)")
    # read groups (movie, chemistry) of every raw BAM: the batch key of kinetics
    mov, unread = [], []
    for s, v in sample.items():
        if not v["raw"]:
            continue
        seen.append(v["raw"])
        m = movies(v["raw"])
        if m is None:
            unread.append(s)
        else:
            mov += [(s,) + r for r in m]
    report(f"raw BAM read groups: {len(mov)} in {len(sample) - len(unread)} BAMs"
           + (f"   header not readable: {len(unread)} (e.g. {unread[0]})" if unread else ""))

    # ---- MAGvsGC (GeneCat getRepresentBins)
    mvg = os.path.join(a.genecat, f"Bin_{a.binner}", "MAGvsGC.txt.gz")
    rep, members, qual = {}, collections.defaultdict(list), {}
    with gzip.open(mvg, "rt") as fh:
        hd = fh.readline().rstrip("\n").split("\t")
        if "Representative4MGS" not in hd:
            stop(f"no Representative4MGS column in {mvg}")
        ir = hd.index("Representative4MGS")
        ic = hd.index("Completeness") if "Completeness" in hd else None
        ix = hd.index("Contamination") if "Contamination" in hd else None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            mag, mgs = f[0], f[1]
            canopy = mag.startswith("Cano__")
            if not canopy and ((ir < len(f) and f[ir] == "*") or mgs not in rep):
                rep[mgs] = mag
            if mgs in ("", "?") or canopy:
                continue
            members[mgs].append(mag)
            qual[mag] = (f[ic] if ic is not None and ic < len(f) else "NA",
                         f[ix] if ix is not None and ix < len(f) else "NA")
    rep = {m: r for m, r in rep.items() if m in members}
    report(f"MGS in MAGvsGC: {len(members)}   member MAGs (non-canopy): {sum(map(len, members.values()))}")
    abund = read_abundance(a.abundance, list(sample), list(members), report)
    if a.abundance and os.path.isfile(a.abundance):
        seen.append(a.abundance)

    # ---- MAG -> SmplID -> assembly -> bin file -> contigs
    mag_info, by_binfile = {}, collections.defaultdict(dict)
    for mgs, mags in members.items():
        for mag in mags:
            info = dict(MGS=mgs, sample="", asm="", contigs=[], status="ok")
            mag_info[mag] = info
            m = NAME_RE.match(mag)
            if not m:
                info["status"] = "bad_name"; continue
            s = alias.get(m.group(1), m.group(1))
            if s not in sample:
                info["status"] = "sample_not_in_map"; continue
            info["sample"], info["asm"] = s, sample[s]["asm"]
            if not info["asm"]:
                info["status"] = "no_assembly"; continue
            binfile = os.path.join(info["asm"], "Binning", a.binner, s)
            by_binfile[binfile].setdefault(m.group(2), []).append(mag)
    for binfile, bins in by_binfile.items():
        seen.append(binfile)
        if not os.path.isfile(binfile):
            for mags in bins.values():
                for mag in mags:
                    mag_info[mag]["status"] = "no_bin_file"
            continue
        with open(binfile) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2 or f[0] == "Sequence ID" or f[1] not in bins:
                    continue
                for mag in bins[f[1]]:
                    mag_info[mag]["contigs"].append(f[0])
        for mags in bins.values():
            for mag in mags:
                if not mag_info[mag]["contigs"]:
                    mag_info[mag]["status"] = "bin_not_in_bin_file"
    st = collections.Counter(i["status"] for i in mag_info.values())
    report("member MAGs: " + "  ".join(f"{k} {v}" for k, v in sorted(st.items())))
    for mag, i in list(mag_info.items()):
        if i["status"] != "ok" and st[i["status"]] <= 20:
            report(f"  MAG {mag} ({i['MGS']}): {i['status']}")

    # ---- representatives: GeneCat's FASTA must hold exactly the bin-file contigs
    gdir = os.path.join(a.genecat, f"Bin_{a.binner}", "Genomes", "MGS_ctg")
    mgs_rows, rep_ctg_rows, owner, differ, dup = [], [], {}, [], []
    for mgs in sorted(members, key=lambda x: (len(x), x)):
        r = rep[mgs]
        ri = mag_info.get(r)
        fa = os.path.join(gdir, rep_fasta_name(mgs, r))
        status, lens = "ok", {}
        if ri is None or ri["status"] != "ok":
            status = "rep_" + (ri["status"] if ri else "missing")
        elif not os.path.isfile(fa):
            status = "no_rep_fasta"
        else:
            lens = fasta_lengths(fa)
            if set(lens) != set(ri["contigs"]):
                status = "rep_contigs_differ"
                differ.append(f"{mgs} ({r}): FASTA {len(lens)} contigs, bin file {len(ri['contigs'])}")
        if status == "ok":
            for c, L in lens.items():
                if c in owner:
                    dup.append(f"{c} in {owner[c]} and {mgs}")
                owner[c] = mgs
                rep_ctg_rows.append((mgs, c, L))
        n_ok = sum(mag_info[m]["status"] == "ok" for m in members[mgs])
        mgs_rows.append((mgs, status, r, ri["sample"] if ri else "", ri["asm"] if ri else "", fa if status == "ok" else "",
                         len(lens), sum(lens.values()), len(members[mgs]), n_ok))
    st = collections.Counter(r[1] for r in mgs_rows)
    report("MGS: " + "  ".join(f"{k} {v}" for k, v in sorted(st.items())))
    for r in mgs_rows:
        if r[1] != "ok" and st[r[1]] <= 20:
            report(f"  MGS {r[0]} (representative {r[2]}): {r[1]}")
    if differ:
        fatal.append(f"{len(differ)} representative(s) differ from GeneCat's own FASTA - MMT's naming rules do not "
                     "match this GeneCat version: " + "; ".join(differ[:5]))
    if dup:
        fatal.append(f"{len(dup)} contig name(s) occur in two representatives: " + "; ".join(dup[:5]))

    # ---- write
    with open(os.path.join(a.out, "samples.tsv"), "w") as o:
        o.write("sample\tstatus\traw_bam\tcram\tassembly\n")
        for s, v in sample.items():
            o.write(f"{s}\t{v['status']}\t{v['raw']}\t{v['cram']}\t{v['asm']}\n")
    with open(os.path.join(a.out, "mgs.tsv"), "w") as o:
        o.write("MGS\tstatus\trep_MAG\trep_sample\trep_assembly\trep_fasta\trep_contigs\trep_bp\tn_MAGs\tn_MAGs_ok\n")
        for r in mgs_rows:
            o.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(a.out, "mags.tsv"), "w") as o:
        o.write("MAG\tMGS\tis_rep\tstatus\tsample\tassembly\tn_contigs\tcompleteness\tcontamination\n")
        for mag, i in mag_info.items():
            c, x = qual[mag]
            o.write(f"{mag}\t{i['MGS']}\t{int(rep[i['MGS']] == mag)}\t{i['status']}\t{i['sample']}\t{i['asm']}\t"
                    f"{len(i['contigs'])}\t{c}\t{x}\n")
    with gzip.open(os.path.join(a.out, "contigs.tsv.gz"), "wt", compresslevel=3) as o:
        o.write("assembly\tcontig\tMAG\tMGS\tis_rep\n")
        for mag, i in mag_info.items():
            if i["status"] == "ok":
                isrep = int(rep[i["MGS"]] == mag)
                for c in i["contigs"]:
                    o.write(f"{i['asm']}\t{c}\t{mag}\t{i['MGS']}\t{isrep}\n")
    # assembly x MGS units ; checksums of what a result depends on (sorted "MAG contig" lines)
    md5 = lambda lines: hashlib.md5("\n".join(sorted(lines)).encode()).hexdigest()[:16]
    unit, asm_lines = collections.defaultdict(lambda: {"mags": [], "all": [], "rep": []}), collections.defaultdict(list)
    ok_mgs = {r[0] for r in mgs_rows if r[1] == "ok"}
    for mag, i in mag_info.items():
        if i["status"] != "ok":
            continue
        isrep = rep[i["MGS"]] == mag
        for c in i["contigs"]:
            asm_lines[i["asm"]].append(f"{c}\t{mag}\t{i['MGS']}\t{int(isrep)}")
        if i["MGS"] not in ok_mgs:
            continue
        u = unit[(i["asm"], i["MGS"])]
        u["mags"].append(mag)
        u["all"] += [f"{mag}\t{c}" for c in i["contigs"]]
        if isrep:
            u["rep"] = [f"{mag}\t{c}" for c in i["contigs"]]
    with open(os.path.join(a.out, "units.tsv"), "w") as o:
        o.write("assembly\tMGS\tMAGs\trep_here\tmd5_all\tmd5_rep\n")
        for (asm, mgs), u in sorted(unit.items()):
            o.write(f"{asm}\t{mgs}\t{','.join(sorted(u['mags']))}\t{int(bool(u['rep']))}\t{md5(u['all'])}\t"
                    f"{md5(u['rep']) if u['rep'] else ''}\n")
    with open(os.path.join(a.out, "assemblies.tsv"), "w") as o:
        o.write("assembly\tmd5_contigs\n")
        for asm in sorted(asm_lines):
            o.write(f"{asm}\t{md5(asm_lines[asm])}\n")
    with gzip.open(os.path.join(a.out, "rep_contigs.tsv.gz"), "wt", compresslevel=3) as o:
        o.write("MGS\tcontig\tlength\n")
        for r in rep_ctg_rows:
            o.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(a.out, "abundance.tsv"), "w") as o:
        o.write("sample\tMGS\tabundance\n")
        for r in abund:
            o.write("\t".join(r) + "\n")
    with open(os.path.join(a.out, "movies.tsv"), "w") as o:
        o.write("sample\t" + "\t".join(MOVIE_COLS) + "\n")
        for r in mov:
            o.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(a.out, "inputs.tsv"), "w") as o:      # same size:mtime format as common.sh stamp()
        for f in sorted(set(seen)):
            st = os.stat(f) if os.path.exists(f) else None
            o.write(f"{f}\t{f'{st.st_size}:{int(st.st_mtime)}' if st else 'missing'}\n")
    with open(os.path.join(a.out, "report.txt"), "w") as o:
        o.write("\n".join(rep_lines + [f"FATAL: {x}" for x in fatal]) + "\n")
    print("\n".join(rep_lines))
    if fatal:
        stop(" | ".join(fatal))


# ============================================================================ assign
RE_SPAN = re.compile(r"(\d+)[MDN=X]")
BIN = 1000


def add_depth(d, p0, span, clen):
    """reference bp of an alignment [p0, p0 + span) (0-based) into 1 kb bins d"""
    p1 = min(p0 + span, clen)                              # end exclusive
    if span <= 0 or p1 <= p0:
        return
    b0, b1 = p0 // BIN, (p1 - 1) // BIN
    if b0 == b1:
        d[b0] += p1 - p0
    else:
        d[b0] += BIN * (b0 + 1) - p0
        for b in range(b0 + 1, b1):
            d[b] += BIN
        d[b1] += p1 - BIN * b1


def assign(a):
    """every primary, non-QC-fail alignment on an MGS-MAG contig -> one row (duplicates kept, flagged: consumers
    decide) ; depth = reference bp of the primary and supplementary, non-QC-fail, non-duplicate alignments per 1 kb bin
    of every MGS-MAG contig (zero bins included: contig lengths from the CRAM header)"""
    c2 = {}
    with gzip.open(a.contigs, "rt") as fh:
        next(fh)
        for line in fh:
            asm, ctg, mag, mgs, isrep = line.rstrip("\n").split("\t")
            if asm == a.assembly:
                c2[ctg] = (mgs, mag, isrep)
    n, seen, clen, depth = collections.Counter(), set(), {}, {}
    with open(a.out, "w") as o:
        o.write("read\tMGS\tMAG\tis_rep\tcontig\tpos\tmapq\tdup\n")
        for line in sys.stdin:
            if line[0] == "@":
                if line.startswith("@SQ"):
                    t = dict(x.split(":", 1) for x in line.rstrip("\n").split("\t")[1:] if ":" in x)
                    if t.get("SN") in c2:
                        clen[t["SN"]] = int(t["LN"])
                        depth[t["SN"]] = [0] * ((int(t["LN"]) + BIN - 1) // BIN)
                continue
            f = line.split("\t", 7)
            q, flag, ctg, pos, mapq, cig = f[0], int(f[1]), f[2], f[3], f[4], f[5]
            if flag & 0x104:                       # unmapped / secondary (normally removed upstream)
                continue
            if flag & 0x800:                       # supplementary: depth only
                d = depth.get(ctg)
                if d is not None and not flag & 0x600:
                    n["depth_supplementary_alignments"] += 1
                    add_depth(d, int(pos) - 1, sum(int(x) for x in RE_SPAN.findall(cig)), clen[ctg])
                continue
            n["cram_primary_alignments"] += 1
            if flag & 0x200:
                n["skipped_qcfail_0x200"] += 1
                continue
            dup = 1 if flag & 0x400 else 0
            hit = c2.get(ctg)
            if hit is None:
                n["not_on_an_MGS_MAG"] += 1
                continue
            if q in seen:
                n["skipped_repeated_read"] += 1
                continue
            seen.add(q)
            n["assigned_reads"] += 1
            n["assigned_reads_on_representative"] += hit[2] == "1"
            n["assigned_duplicates"] += dup
            o.write(f"{q}\t{hit[0]}\t{hit[1]}\t{hit[2]}\t{ctg}\t{pos}\t{mapq}\t{dup}\n")
            d = depth.get(ctg)
            if dup or d is None:
                continue
            add_depth(d, int(pos) - 1, sum(int(x) for x in RE_SPAN.findall(cig)), clen[ctg])
    n["mgs_contigs_in_assembly"] = len(c2)
    n["mgs_contigs_in_cram_header"] = len(clen)
    with open(a.stats, "w") as s:
        for k in ("cram_primary_alignments", "skipped_qcfail_0x200", "skipped_repeated_read", "not_on_an_MGS_MAG",
                  "assigned_reads", "assigned_reads_on_representative", "assigned_duplicates", "mgs_contigs_in_assembly",
                  "mgs_contigs_in_cram_header", "depth_supplementary_alignments"):
            s.write(f"{k}\t{n[k]}\n")
    with open(a.depth, "w") as o:
        o.write("MGS\tMAG\tcontig\tc_start\tc_end\tbp\n")
        for c in sorted(depth):
            for b, v in enumerate(depth[c]):
                o.write(f"{c2[c][0]}\t{c2[c][1]}\t{c}\t{b * BIN + 1}\t{min(clen[c], (b + 1) * BIN)}\t{v}\n")
    if n["cram_primary_alignments"] == 0:
        stop("no primary alignment read from the CRAM")
    if c2 and n["assigned_reads"] == 0 and n["not_on_an_MGS_MAG"] > 0:
        stop(f"none of {n['not_on_an_MGS_MAG']} aligned reads lies on one of the {len(c2)} MGS contigs of this assembly - "
             "do the contig names in the CRAM match the bin files?")


# ============================================================================ split
def split(a):
    import pysam
    want = {}
    with open(a.reads) as fh:
        next(fh)
        for line in fh:
            r, m = line.rstrip("\n").split("\t")[:2]
            want[r] = m
    groups = sorted(set(want.values()))
    os.makedirs(a.outdir, exist_ok=True)
    counts, mols = collections.Counter(), collections.defaultdict(set)
    for b in range(0, len(groups), a.max_open):          # stay below the open-file limit: one pass per batch
        batch = set(groups[b:b + a.max_open])
        with pysam.AlignmentFile(a.bam, "rb", check_sq=False, threads=a.threads) as f:
            outs = {m: pysam.AlignmentFile(os.path.join(a.outdir, f".{m}.bam"), "wb", template=f) for m in batch}
            for r in f:
                base = r.query_name.rsplit("/", 1)[0]       # <movie>/<zmw>/ccs/fwd -> <movie>/<zmw>/ccs
                m = want.get(base)
                if m in batch:
                    outs[m].write(r)
                    counts[m] += 1
                    mols[m].add(base)
            for o in outs.values():
                o.close()
    with open(os.path.join(a.outdir, "counts.tsv"), "w") as o:
        o.write("group\tmolecules\tpseudo_reads\n")
        for m in groups:
            if counts[m]:
                os.replace(os.path.join(a.outdir, f".{m}.bam"), os.path.join(a.outdir, f"{m}.bam"))
                o.write(f"{m}\t{len(mols[m])}\t{counts[m]}\n")
            else:
                os.remove(os.path.join(a.outdir, f".{m}.bam"))
    if want and not counts:
        stop("no by-strand read matched the selected reads - check the bystrandify read names (<read>/fwd|rev)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("index")
    p.add_argument("--mf4-run", required=True); p.add_argument("--genecat", required=True)
    p.add_argument("--map", required=True); p.add_argument("--binner", default="SB"); p.add_argument("--out", required=True)
    p.add_argument("--abundance", default="")
    p = sp.add_parser("assign")
    p.add_argument("--contigs", required=True); p.add_argument("--assembly", required=True)
    p.add_argument("--out", required=True); p.add_argument("--stats", required=True); p.add_argument("--depth", required=True)
    p = sp.add_parser("split")
    p.add_argument("--bam", required=True); p.add_argument("--reads", required=True)
    p.add_argument("--outdir", required=True); p.add_argument("--max-open", type=int, default=400)
    p.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    {"index": index, "assign": assign, "split": split}[a.cmd](a)


if __name__ == "__main__":
    main()
