#!/usr/bin/env python3
"""mmt export-mot: a finished MATAFILER4 run + its GeneCat MGS stage -> a Methyl-o-Type (mot) input table + library.

  python3 motexport.py --run MF4_RUN --genecat GC --table OUT.tsv [--map MAP] [--binner TAG] [--min-mapq N]
                       [--exclude-qcfail] [--abundance FILE] [--threads N] [--samtools PATH]

The table is mot 0.7's input table (mot help input): sample genome bam reads ref mgs abundance role, one row per
sample x MAG x raw BAM; sample = the SmplID, genome = the MAG, mgs = its MGS, abundance = GeneCat's MGS abundance.

What is used (the rules of GeneCat's and MATAFILER4's own code; see docs/export-mot.md):
  maps         --map, else the maps GeneCat was run with (<GC>/LOGandSUB/GCmaps.ori), else GeneCat's copies
               (LOGandSUB/inmap.txt -> map.N.txt); samples = map rows with a folder in the run
  final MAGs   every non-canopy row of <GC>/Bin_<binner>/MAGvsGC.txt.gz (MGS "?" included: binned, quality-passed,
               no MGS) ; contigs from the bin file <assembly>/Binning/<binner>/<SmplID> (mf4.bin_contigs)
  sequence     <assembly>/scaffolds.fasta.filt(.gz) = the binning input and MATAFILER4's mapping reference ;
               circular = listed in <assembly>/scaffolds.fasta.circ ; checked against GeneCat's MGS_ctg FASTAs
  raw BAMs     <run>/<S>/input_raw.txt (MATAFILER4's record of the files it read), else the map (MF4 rule) ;
               every molecule is attributed to a BAM by its movie (BAM read groups + first records)
  reads        <run>/<S>/mapping/<S>-smd.cram: the sample's reads aligned to its own (co-)assembly; with
               -retainAllCRAM it keeps the alignments MATAFILER4 rejected (0x200 QC fail, 0x400 duplicate)

Assignment: competitive within the (co-)assembly a sample was mapped to - not across the cohort (each CRAM saw only
its own assembly). Per molecule (movie/zmw) of a CRAM, from its primary + supplementary records:
  unique                all on contigs of one final MAG, primary MAPQ >= --min-mapq  -> selected for mot
  ambiguous_low_mapq    one MAG, primary MAPQ below the threshold (a near-equal placement elsewhere in the assembly)
  split_across_mags     pieces on two or more MAGs
  split_mag_unbinned    pieces on one MAG and on contigs of no final MAG
  contig_in_two_mags    on a contig two final MAGs claim
  reference_mismatch    on a contig whose CRAM header length / M5 differs from the exported sequence
  no_primary            only supplementary records
Molecules only on contigs of no final MAG are counted, not listed. Secondary records never decide (MAPQ already
reflects competition) but are recorded. QC-failed / duplicate-flagged molecules are kept with their flags and selected
when unique (mot filters alignments itself) unless --exclude-qcfail.

Exit status: 0 COMPLETE, 3 PARTIAL (rows written for everything that validated; the rest listed), 1 FAILED.
"""
import argparse, collections, concurrent.futures as cf, datetime, gzip, hashlib, io, json, os, re, shutil, subprocess
import sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mf4                                                                     # noqa: E402

METHOD = 1          # bump when per-sample outputs change -> every sample is redone on the next run
STATUS = ["unique", "ambiguous_low_mapq", "split_across_mags", "split_mag_unbinned", "contig_in_two_mags",
          "reference_mismatch", "no_primary"]
BAMST = ["ok", "movie_not_in_bams", "movie_in_two_bams", "not_pacbio_name", "bam_not_mot_ready"]
MOVIE_RE = re.compile(r"^m\d+[A-Za-z]?_\d{6}_\d{6}(_s\d+)?$")   # PacBio movie (RS II, Sequel, Revio)
SAFE = re.compile(r"[^A-Za-z0-9._-]")
MOT_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,99}$")       # mot's sample / genome name rule (src/inputs.cpp)
MOT_COLUMNS = ["sample", "genome", "bam", "reads", "ref", "mgs", "abundance", "role"]   # mot 0.7 input table
OFF, BAD, SHARED = -1, -2, -3                                      # contig codes besides a MAG index


def log(msg):
    print(f"[{time.strftime('%F %T')}] {msg}", flush=True)


def md5s(*parts):
    return hashlib.md5("|".join(map(str, parts)).encode()).hexdigest()[:16]


def stamp(p):
    try:
        st = os.stat(p)
        return f"{st.st_size}:{int(st.st_mtime)}"
    except OSError:
        return "missing"


def safe(name):
    return SAFE.sub("_", name)


def mot_name(name):      # a mot sample / genome name (it names an output folder): safe characters, <= 100 (cut + checksum)
    s = safe(name)
    return s if len(s) <= 100 else f"{s[:91]}_{md5s(name)[:8]}"


def num(x):              # a finite number, else NaN (every comparison with NaN is False)
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v if v == v and v not in (float("inf"), float("-inf")) else float("nan")


def stem(mag):           # SemiBin MAG names end in .fa.gz: file and sample names without it
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", mag, flags=re.I) or mag


def write_atomic(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "w") as o:
        o.write(text)
    os.replace(tmp, path)


def tsv(path, head, rows):
    write_atomic(path, "\t".join(head) + "\n" + "".join("\t".join("" if x is None else str(x) for x in r) + "\n"
                                                         for r in rows))


def gz_lines(path, lines):          # deterministic gzip (mtime 0): a rerun gives the same bytes
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, compresslevel=6) as o:
        o.write(("\n".join(lines) + "\n").encode())
    os.replace(tmp, path)


def kits(chem):                     # mot's same_chemistry: binding kit + sequencing kit
    return tuple(chem.split("/")[:2])


# ============================================================================ discovery
def find_binner(gc, explicit, default):
    have = sorted(d[4:] for d in os.listdir(gc)
                  if d.startswith("Bin_") and os.path.isfile(os.path.join(gc, d, "MAGvsGC.txt.gz")))
    if explicit:
        if explicit not in have:
            mf4.stop(f"--binner {explicit}: no {gc}/Bin_{explicit}/MAGvsGC.txt.gz (found: {', '.join(have) or 'none'})")
        return explicit
    if len(have) == 1:
        return have[0]
    if not have:
        mf4.stop(f"no Bin_<binner>/MAGvsGC.txt.gz in {gc} - is this the GeneCat folder and has its MGS stage finished?")
    if default in have:
        return default
    mf4.stop(f"several binners in {gc} ({', '.join(have)}): choose one with --binner")


def find_maps(gc, override):
    if override:
        return [os.path.abspath(override)], "--map"
    lg = os.path.join(gc, "LOGandSUB")
    ori = os.path.join(lg, "GCmaps.ori")
    if os.path.isfile(ori):
        maps = [m for m in open(ori).read().strip().split(",") if m]
        if maps and all(os.path.isfile(m) for m in maps):
            return maps, f"the maps GeneCat was run with ({ori})"
    inmap = os.path.join(lg, "inmap.txt")
    maps = []
    if os.path.isfile(inmap):
        maps = [m if os.path.isabs(m) else os.path.join(lg, m) for m in open(inmap).read().split()]
    if not maps and os.path.isdir(lg):                               # GeneCat's own fallback: LOGandSUB/map.N.txt
        maps = sorted((os.path.join(lg, f) for f in os.listdir(lg) if re.fullmatch(r"map\.\d+\.txt", f)),
                      key=lambda p: int(p.rsplit(".", 2)[-2]))
    maps = [m for m in maps if os.path.isfile(m)]
    if not maps:
        mf4.stop(f"no map: give --map FILE (GeneCat's {lg}/GCmaps.ori and inmap.txt name none that exists)")
    return maps, f"GeneCat's copies of the maps ({lg})"


def bam_probe(path):
    """-> dict(ok, movies, chemistry, read_types, kinetics, error) from the header and the first 200 records"""
    import pysam
    r = dict(ok=False, movies=[], chemistry="", read_types="", kinetics=False, error="")
    try:
        with pysam.AlignmentFile(path, "rb", check_sq=False) as f:
            rgs = f.header.to_dict().get("RG", [])
            mv, chem, rt = [], set(), set()
            for rg in rgs:
                ds = dict(kv.split("=", 1) for kv in str(rg.get("DS", "")).split(";") if "=" in kv)
                if rg.get("PU") and rg["PU"] not in mv:
                    mv.append(rg["PU"])
                chem.add(f"{ds.get('BINDINGKIT', '')}/{ds.get('SEQUENCINGKIT', '')}")
                rt.add(ds.get("READTYPE", ""))
            for i, a in enumerate(f.fetch(until_eof=True)):
                m = a.query_name.split("/", 1)[0]
                if "/" in a.query_name and m not in mv:
                    mv.append(m)
                r["kinetics"] = r["kinetics"] or a.has_tag("fi") or a.has_tag("ri") or a.has_tag("ip")
                if i >= 199:
                    break
        r.update(ok=True, movies=mv, chemistry=",".join(sorted(chem)), read_types=",".join(sorted(rt)))
    except Exception as e:                                        # not a BAM, truncated, unreadable
        r["error"] = f"{type(e).__name__}: {str(e)[:100]}"
    return r


def choose_bams(s, run, map_bams, probe):
    """the sample's raw BAMs: MATAFILER4's record (input_raw.txt) when every BAM there opens, else the map ->
    (paths, provenance, notes, other inputs)"""
    notes, others = [], []
    raw = os.path.join(run, s, "input_raw.txt")
    rec = []
    if os.path.isfile(raw):
        for p in re.split(r"[;\n]", open(raw).read()):
            p = p.strip()
            if p:
                (rec if p.endswith(".bam") else others).append(p)
    good = lambda ps: [p for p in ps if os.path.isfile(p) and probe(p)["ok"]]
    if rec and good(rec) == rec:
        if {os.path.realpath(p) for p in rec} != {os.path.realpath(p) for p in map_bams}:
            notes.append(f"the map names other BAMs ({len(map_bams)}) than MATAFILER4 read ({len(rec)}): used "
                         "MATAFILER4's (input_raw.txt)")
        return rec, "input_raw.txt", notes, others
    if rec:
        notes.append(f"input_raw.txt names {len(rec) - len(good(rec))} BAM(s) that are missing or do not open "
                     f"(e.g. {next(p for p in rec if p not in good(rec))})")
    if good(map_bams):
        if rec:
            notes.append("used the BAMs the map names now (relocated?) - checked by movie against the CRAM")
        return good(map_bams), "map", notes, others
    if good(rec):
        return good(rec), "input_raw.txt (partly)", notes, others
    if map_bams:
        notes.append(f"the map's BAM(s) do not open: {map_bams[0]}")
    return [], "none", notes, others


def discover(a):
    rep_lines, fatal = [], []
    report = rep_lines.append
    run, gc = os.path.abspath(a.run), os.path.abspath(a.genecat)
    binner = find_binner(gc, a.binner, a.binner_default)
    maps, map_src = find_maps(gc, a.map)
    report(f"binner Bin_{binner} ; maps: {', '.join(maps)} ({map_src})")

    # ---- samples of this run
    rows, alias, not_in_run = [], {}, []
    for m in maps:
        r, al = mf4.scan_map(m)
        rows += [dict(x, map=m) for x in r]
        alias.update(al)
    ids = collections.Counter(r["sample"] for r in rows)
    if any(k > 1 for k in ids.values()):
        fatal.append("sample(s) more than once in the maps: " + ", ".join(s for s, k in ids.items() if k > 1)[:300])
    probes = {}

    def probe(p):
        if p not in probes:
            probes[p] = bam_probe(p)
        return probes[p]
    samples = collections.OrderedDict()
    for r in rows:
        s = r["sample"]
        sd = os.path.join(run, s)
        if not os.path.isdir(sd):
            not_in_run.append(s); continue
        at = os.path.join(sd, "assemblies", "metag", "assembly.txt")
        asm = os.path.normpath(mf4.first_line(at)) if os.path.isfile(at) and mf4.first_line(at) else ""
        bams, prov, notes, others = choose_bams(s, run, r["bams"], probe)
        if r["problem"]:
            notes.append(r["problem"])
        cram = os.path.join(sd, "mapping", f"{s}-smd.cram")
        samples[s] = dict(sample=s, group=r["group"], asm=asm if asm and os.path.isdir(asm) else "", cram=cram,
                          has_cram=os.path.isfile(cram), bams=bams, provenance=prov, notes=notes,
                          other_inputs=others, map=r["map"])
    if not samples:
        mf4.stop(f"no sample of the maps has a folder in {run} - is -i the MATAFILER4 run folder (#OutPath/#RunID)?")
    report(f"samples: {len(samples)} in the run" + (f" ; {len(not_in_run)} map rows without a run folder "
                                                     f"(e.g. {not_in_run[0]})" if not_in_run else ""))
    bams = collections.OrderedDict()
    for s, v in samples.items():
        for p in v["bams"]:
            pr = probe(p)
            b = bams.setdefault(os.path.realpath(p), dict(path=os.path.abspath(p), samples=[], **pr))
            b["samples"].append(s)
    for b in bams.values():
        b["mot_ready"] = b["ok"] and b["kinetics"] and set(b["read_types"].split(",")) <= {"CCS"}
        b["why"] = ("" if b["mot_ready"] else b["error"] or ("no kinetics (fi/ri/ip) in the first records"
                    if not b["kinetics"] else f"READTYPE {b['read_types']} (mot needs CCS / HiFi)"))

    # ---- final MAGs (MAGvsGC, non-canopy) and their contigs
    mvg = os.path.join(gc, f"Bin_{binner}", "MAGvsGC.txt.gz")
    mags, qual, rep, seen_mag = collections.OrderedDict(), {}, {}, collections.Counter()
    with gzip.open(mvg, "rt") as fh:
        hd = fh.readline().rstrip("\n").split("\t")
        ir = hd.index("Representative4MGS") if "Representative4MGS" in hd else None
        ic = hd.index("Completeness") if "Completeness" in hd else None
        ix = hd.index("Contamination") if "Contamination" in hd else None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or f[0].startswith("Cano__"):
                continue
            mag, mgs = f[0], f[1] if f[1] not in ("?", "") else ""
            seen_mag[mag] += 1
            if mgs and ((ir is not None and ir < len(f) and f[ir] == "*") or mgs not in rep):
                rep[mgs] = mag
            mags[mag] = mgs
            qual[mag] = (f[ic] if ic is not None and ic < len(f) else "NA", f[ix] if ix is not None and ix < len(f) else "NA")
    dup = [m for m, k in seen_mag.items() if k > 1]
    if dup:
        report(f"MAGvsGC lists {len(dup)} MAG(s) more than once (one MAG each kept): {', '.join(dup[:5])}")
    info, _ = mf4.bin_contigs(mags, alias, {s: v["asm"] for s, v in samples.items()}, binner)
    for m, i in info.items():
        i.update(is_rep=int(bool(i["MGS"]) and rep.get(i["MGS"]) == m), completeness=qual[m][0],
                 contamination=qual[m][1], stem=stem(m))
    # contigs claimed by two MAGs of one assembly
    owner = collections.defaultdict(list)
    for m, i in info.items():
        if i["status"] == "ok":
            for c in i["contigs"]:
                owner[(i["asm"], c)].append(m)
    shared = {k: v for k, v in owner.items() if len(v) > 1}
    for (asm, c), ms in shared.items():
        for m in ms:
            info[m].setdefault("shared", []).append(c)
    st = collections.Counter(i["status"] for i in info.values())
    report(f"final MAGs (MAGvsGC, non-canopy): {len(info)}  " + "  ".join(f"{k} {v}" for k, v in sorted(st.items())))
    if shared:
        report(f"WARNING: {len(shared)} contig(s) claimed by two MAGs of one assembly (reads there: contig_in_two_mags)")
    # unique file stems
    by_stem = collections.defaultdict(list)
    for m, i in info.items():
        by_stem[safe(i["stem"])].append(m)
    for st_, ms in by_stem.items():
        for m in ms:
            info[m]["file"] = st_ if len(ms) == 1 else f"{st_}_{md5s(m)[:8]}"
    return dict(run=run, gc=gc, binner=binner, maps=maps, map_src=map_src, samples=samples, bams=bams, mags=info,
                alias=alias, rep=rep, report=rep_lines, fatal=fatal, not_in_run=not_in_run)


# ============================================================================ genomes
def read_fasta(path, want=None):
    """yield (name, sequence) ; with want: only those names (the rest is skipped without joining lines)"""
    name, buf = None, []
    with mf4.opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (want is None or name in want):
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            elif want is None or name in want:
                buf.append(line.strip())
        if name is not None and (want is None or name in want):
            yield name, "".join(buf)


def asm_fasta(asm):
    for f in ("scaffolds.fasta.filt", "scaffolds.fasta.filt.gz"):
        if os.path.isfile(os.path.join(asm, f)) and os.path.getsize(os.path.join(asm, f)):
            return os.path.join(asm, f)
    return ""


def genome_job(job):
    """one assembly: every final MAG's contigs -> <lib>/genomes/<file>.fa.gz (bgzip + .fai + .gzi) ; contig table"""
    import pysam
    asm, mags, lib, gcdir, fp = job["asm"], job["mags"], job["lib"], job["gcdir"], job["fp"]
    gdir = os.path.join(lib, "genomes")
    keep = os.path.join(gdir, ".fp", md5s(asm))
    if os.path.isfile(keep + ".fp") and open(keep + ".fp").read() == fp and all(
            os.path.isfile(os.path.join(gdir, m["file"] + ".fa.gz.fai")) for m in mags.values()):
        return json.load(open(keep + ".json"))
    fa = asm_fasta(asm)
    res = dict(asm=asm, fasta=fa, contigs=[], mags={})
    want = {c: m for m, v in mags.items() for c in v["contigs"]}
    if not fa:
        for m in mags:
            res["mags"][m] = dict(status="no_assembly_fasta")
        return res
    circ = {}
    cf_ = os.path.join(asm, "scaffolds.fasta.circ")
    if os.path.isfile(cf_):
        for n, sq in read_fasta(cf_, set(want)):
            circ[n] = hashlib.md5(sq.upper().encode()).hexdigest()
    circ_known = os.path.isfile(cf_)
    tmpd = os.path.join(gdir, f".{md5s(asm)}.{os.getpid()}.tmp")
    shutil.rmtree(tmpd, ignore_errors=True); os.makedirs(tmpd)
    outs, got = {}, collections.defaultdict(dict)
    for n, sq in read_fasta(fa, set(want)):
        m = want[n]
        if n in got[m]:
            continue                                               # a duplicate record: the first one counts
        if m not in outs:
            outs[m] = open(os.path.join(tmpd, mags[m]["file"] + ".fa"), "w")
        outs[m].write(f">{n}\n" + "\n".join(sq[i:i + 80] for i in range(0, len(sq), 80)) + "\n")
        h = hashlib.md5(sq.upper().encode()).hexdigest()
        got[m][n] = (len(sq), h)
    for o in outs.values():
        o.close()
    for m, v in mags.items():
        miss = [c for c in v["contigs"] if c not in got[m]]
        r = dict(status="ok", n_contigs=len(got[m]), bp=sum(L for L, _ in got[m].values()), missing=miss[:5],
                 n_missing=len(miss), n_circular=sum(c in circ for c in got[m]), circ_differs=[],
                 md5=hashlib.md5("".join(f"{c}:{h}\n" for c, (_, h) in sorted(got[m].items())).encode()).hexdigest())
        r["circ_differs"] = [c for c in got[m] if c in circ and circ[c] != got[m][c][1]][:5]
        if miss:
            r["status"] = "contigs_missing_from_assembly"
        # GeneCat's own FASTA of this MAG (representatives at least): same contigs, same sequence
        gfa = os.path.join(gcdir, mf4.rep_fasta_name(v["MGS"], m)) if v["MGS"] else ""
        r["genecat_fasta"] = gfa if gfa and os.path.isfile(gfa) else ""
        if r["genecat_fasta"] and r["status"] == "ok":
            g = {n: hashlib.md5(sq.upper().encode()).hexdigest() for n, sq in read_fasta(gfa)}
            if set(g) != set(v["contigs"]):
                r["status"] = "contigs_differ_from_genecat"
            elif any(g[c] != got[m][c][1] for c in g):
                r["status"] = "sequence_differs_from_genecat"
        if r["status"] == "ok":
            src = os.path.join(tmpd, v["file"] + ".fa")
            dst = os.path.join(gdir, v["file"] + ".fa.gz")
            for x in (dst, dst + ".fai", dst + ".gzi"):
                if os.path.exists(x):
                    os.remove(x)
            pysam.tabix_compress(src, dst, force=True)               # BGZF: mot reads it with its .fai
            pysam.faidx(dst)
            for c, (L, h) in got[m].items():
                res["contigs"].append((m, c, L, h, (1 if c in circ else 0) if circ_known else "NA"))
        res["mags"][m] = r
    shutil.rmtree(tmpd, ignore_errors=True)
    os.makedirs(os.path.dirname(keep), exist_ok=True)
    write_atomic(keep + ".json", json.dumps(res))
    write_atomic(keep + ".fp", fp)
    return res


# ============================================================================ CRAM -> molecules
def sample_job(job):
    """one sample's CRAM -> molecule catalog (parquet), selections per MAG x BAM (gz), stats (json)"""
    import array
    import numpy as np, pyarrow as pa, pyarrow.parquet as pq
    S, lib = job["sample"], job["lib"]
    work = os.path.join(lib, ".work"); os.makedirs(work, exist_ok=True)
    fpf, stf = os.path.join(work, f"{safe(S)}.fp"), os.path.join(work, f"{safe(S)}.stats.json")
    if os.path.isfile(fpf) and open(fpf).read() == job["fp"] and os.path.isfile(stf):
        st = json.load(open(stf)); st["reused"] = True
        return st
    t0 = time.time()
    mags = job["mags"]                              # [MAG] in this assembly (index = code)
    ref = job["ref"]                                # contig -> (code, length, md5) for contigs of final MAGs
    movie_bam = job["movie_bam"]                    # movie -> BAM index, -1 = in two BAMs
    st = dict(sample=S, cram=job["cram"], method=METHOD, records=0, primary=0, supplementary=0, secondary=0,
              qcfail_records=0, dup_records=0, non_pacbio_names=0, non_pacbio_example="", header_contigs=0,
              mag_contigs=len(ref), mag_contigs_not_in_cram=0, mag_contigs_mismatch=0, mismatch_examples=[],
              cram_has_m5=False)
    env = dict(os.environ, REF_PATH="/nonexistent/mmt", REF_CACHE="/nonexistent/mmt")
    p = subprocess.Popen([job["samtools"], "view", "-h", "--no-PG", "-F", "0x4", "--input-fmt-option",
                          "required_fields=31", job["cram"]], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         text=True, bufsize=1 << 20, env=env)     # required_fields 31 = 0x1f: newer htslib rejects hex
    code, ctg_names, cidx = {}, [], {}              # CRAM contig -> MAG code / OFF / BAD / SHARED ; -> index
    movies, mv_idx = [], {}
    K, C, Q, F, P, I = (array.array(t) for t in ("q", "i", "B", "H", "I", "i"))
    for line in p.stdout:
        if line[0] == "@":
            if line.startswith("@SQ"):
                t = dict(x.split(":", 1) for x in line.rstrip("\n").split("\t")[1:] if ":" in x)
                n = t.get("SN", "")
                cidx[n] = len(ctg_names); ctg_names.append(n)
                st["header_contigs"] += 1
                st["cram_has_m5"] = st["cram_has_m5"] or "M5" in t
                if n in ref:
                    c, L, h = ref[n]                    # L None: a MAG that was not exported (no sequence to compare)
                    ok = L is None or (str(L) == t.get("LN") and ("M5" not in t or t["M5"].lower() == h))
                    if c >= 0 and not ok:
                        st["mag_contigs_mismatch"] += 1
                        if len(st["mismatch_examples"]) < 5:
                            st["mismatch_examples"].append(f"{n}: CRAM LN {t.get('LN')} M5 {t.get('M5', '-')[:8]} "
                                                           f"vs FASTA {L} {h[:8]}")
                    code[n] = c if ok else BAD
                else:
                    code[n] = OFF
            continue
        f = line.split("\t", 6)
        q, flag = f[0], int(f[1])
        st["records"] += 1
        if flag & 0x100:
            st["secondary"] += 1
        elif flag & 0x800:
            st["supplementary"] += 1
        else:
            st["primary"] += 1
        st["qcfail_records"] += bool(flag & 0x200); st["dup_records"] += bool(flag & 0x400)
        i = q.find("/")
        j = q.find("/", i + 1)
        zs = q[i + 1:j] if j > 0 else q[i + 1:]
        if i <= 0 or not zs.isdigit():
            st["non_pacbio_names"] += 1
            st["non_pacbio_example"] = st["non_pacbio_example"] or q[:80]
            continue
        mv = q[:i]
        mi = mv_idx.get(mv)
        if mi is None:
            mi = mv_idx[mv] = len(movies); movies.append(mv)
        K.append(mi << 32 | int(zs)); C.append(code.get(f[2], OFF)); Q.append(min(int(f[4]), 255))
        F.append(flag); P.append(int(f[3])); I.append(cidx.get(f[2], -1))
    err = p.stderr.read(); rc = p.wait()
    if rc != 0:
        raise RuntimeError(f"{S}: samtools view {job['cram']} failed (exit {rc}): {err.strip()[:300]}")
    st["mag_contigs_not_in_cram"] = sum(1 for n in ref if n not in code)
    K, C, Q, F, P, I = (np.frombuffer(x, dtype=d) if len(x) else np.zeros(0, d) for x, d in
                        ((K, np.int64), (C, np.int32), (Q, np.uint8), (F, np.uint16), (P, np.uint32), (I, np.int32)))
    sec = (F & 0x100) > 0
    prim = (F & 0x900) == 0
    o = np.lexsort((~prim, K))                       # by molecule, primary records first
    K, C, Q, F, P, I, sec, prim = K[o], C[o], Q[o], F[o], P[o], I[o], sec[o], prim[o]
    uk, start, cnt = np.unique(K, return_index=True, return_counts=True)
    n = uk.size
    red = lambda x, fn=np.add: fn.reduceat(x, start) if n else np.zeros(0, x.dtype)
    dec = ~sec                                        # primary + supplementary decide
    n_prim, n_sup, n_sec = red(prim.astype(np.int32)), red((dec & ~prim).astype(np.int32)), red(sec.astype(np.int32))
    any_ = lambda m: red((m & dec).astype(np.int32)) > 0
    has_bad, has_shared, has_off = any_(C == BAD), any_(C == SHARED), any_(C == OFF)
    BIG = np.int32(1 << 30)
    cmin = red(np.where(dec & (C >= 0), C, BIG), np.minimum)
    cmax = red(np.where(dec & (C >= 0), C, np.int32(-1)), np.maximum)
    has_mag = cmax >= 0
    mapq = red(np.where(prim, Q, np.uint8(255)), np.minimum)
    qcf = red((prim & ((F & 0x200) > 0)).astype(np.int32)) > 0
    dup = red((prim & ((F & 0x400) > 0)).astype(np.int32)) > 0
    sec_mag = red(np.where(sec & (C >= 0), C, np.int32(-1)), np.maximum)
    status = np.full(n, 0, np.int8)                  # unique
    status[(n_prim > 0) & (mapq < job["min_mapq"])] = 1
    status[n_prim == 0] = 6
    status[has_mag & has_off] = 3
    status[has_mag & (cmin != cmax)] = 2
    status[has_shared] = 4
    status[has_bad] = 5
    on = has_mag | has_bad | has_shared                # on a final MAG (secondary-only contacts count as off-MAG)
    st["molecules"] = int(n)
    st["molecules_off_mag"] = int((~on).sum())
    mol_mag = np.where(has_mag, np.where(n_prim > 0, C[start], cmin), -1)   # the primary's MAG (else the first MAG)
    mol_mag = np.where(mol_mag >= 0, mol_mag, np.where(has_mag, cmin, -1))
    keep = np.flatnonzero(on)
    mv_of = (uk >> 32).astype(np.int64); zmw = (uk & 0xffffffff).astype(np.uint32)
    bam_of = np.array([movie_bam.get(m, -2) for m in movies], np.int32)[mv_of] if n else np.zeros(0, np.int32)
    pacbio_mv = np.array([bool(MOVIE_RE.match(m)) or m in movie_bam for m in movies], bool)[mv_of] if n else \
        np.zeros(0, bool)
    bst = np.where(bam_of >= 0, 0, np.where(bam_of == -1, 2, np.where(pacbio_mv, 1, 3))).astype(np.int8)
    ready = np.array(job["bam_ready"] + [False], bool)[np.where(bam_of >= 0, bam_of, len(job["bam_ready"]))]
    bst[(bst == 0) & ~ready] = 4
    mag_ok = np.array([job["mag_ok"][m] for m in mags] + [False], bool)[np.where(mol_mag >= 0, mol_mag, len(mags))]
    sel = (status == 0) & (bst == 0) & mag_ok & (~qcf if job["exclude_qcfail"] else True)
    # other MAGs of split / secondary contacts (rare: a python loop over those molecules only)
    other = np.full(n, "", object); secm = np.full(n, "", object)
    for k in np.flatnonzero((status == 2) | ((sec_mag >= 0) & (sec_mag != mol_mag))):
        s0, e0 = start[k], start[k] + cnt[k]
        cs = C[s0:e0]
        if status[k] == 2:
            other[k] = ",".join(sorted({mags[c] for c in cs[dec[s0:e0] & (cs >= 0)]} - {mags[mol_mag[k]]}))
        sm = {mags[c] for c in cs[sec[s0:e0] & (cs >= 0)]} - ({mags[mol_mag[k]]} if mol_mag[k] >= 0 else set())
        secm[k] = ",".join(sorted(sm))
    # ---- selections: one gz per MAG x BAM, molecules sorted by movie name then ZMW
    adir = os.path.join(lib, "assignments")
    tmpd = os.path.join(adir, f".{safe(S)}.{os.getpid()}.tmp")
    shutil.rmtree(tmpd, ignore_errors=True); os.makedirs(tmpd)
    mrank = np.argsort(np.argsort(np.array(movies, object))) if movies else np.zeros(0, int)
    sel_counts = collections.Counter()
    si = np.flatnonzero(sel)
    for mg in np.unique(mol_mag[si]) if si.size else []:
        for b in np.unique(bam_of[si][mol_mag[si] == mg]):
            idx = si[(mol_mag[si] == mg) & (bam_of[si] == b)]
            idx = idx[np.lexsort((zmw[idx], mrank[mv_of[idx]]))]
            d = os.path.join(tmpd, job["mag_file"][mags[mg]])
            os.makedirs(d, exist_ok=True)
            gz_lines(os.path.join(d, job["bam_label"][b] + ".reads.txt.gz"),
                     [f"{movies[m]}/{z}" for m, z in zip(mv_of[idx], zmw[idx])])
            sel_counts[(mags[mg], int(b))] = int(idx.size)
    final = os.path.join(adir, safe(S))
    old = final + f".{os.getpid()}.old"
    if os.path.isdir(final):
        os.replace(final, old)
    os.replace(tmpd, final)
    shutil.rmtree(old, ignore_errors=True)
    # ---- catalog: every molecule touching a final MAG
    k = keep
    dic = lambda idx, vals: pa.DictionaryArray.from_arrays(pa.array(idx, pa.int32()), pa.array(vals, pa.string()))
    mags_x = mags + [""]
    bam_paths = job["bam_paths"] + ["", ""]
    t = pa.table({
        "sample": dic(np.zeros(k.size, np.int32), [S]),
        "MAG": dic(np.where(mol_mag[k] >= 0, mol_mag[k], len(mags)), mags_x),
        "movie": dic(mv_of[k], movies), "zmw": pa.array(zmw[k], pa.uint32()),
        "source_bam": dic(np.where(bam_of[k] >= 0, bam_of[k], len(bam_paths) - 1), bam_paths),
        "status": dic(status[k].astype(np.int32), STATUS), "bam_status": dic(bst[k].astype(np.int32), BAMST),
        "selected": pa.array(sel[k]), "mapq": pa.array(np.where(mapq[k] == 255, 0, mapq[k]), pa.uint8()),
        "primary_mapped": pa.array(n_prim[k] > 0), "qc_fail": pa.array(qcf[k]), "duplicate": pa.array(dup[k]),
        "n_primary": pa.array(n_prim[k], pa.uint16()), "n_supplementary": pa.array(n_sup[k], pa.uint16()),
        "n_secondary": pa.array(n_sec[k], pa.uint16()),
        "contig": dic(np.where(I[start[k]] >= 0, I[start[k]], 0), ctg_names or [""]),
        "pos": pa.array(P[start[k]], pa.uint32()),
        "other_MAGs": pa.array(other[k].tolist(), pa.string()), "secondary_MAGs": pa.array(secm[k].tolist(), pa.string()),
    })
    cdir = os.path.join(lib, "catalog", "molecules"); os.makedirs(cdir, exist_ok=True)
    tp = os.path.join(cdir, f".{safe(S)}.{os.getpid()}.tmp")
    pq.write_table(t, tp, compression="zstd")
    os.replace(tp, os.path.join(cdir, f"{safe(S)}.parquet"))
    # ---- stats
    gi = np.where(mol_mag[k] >= 0, mol_mag[k], len(mags))       # MAG of each listed molecule ("" = none)
    cnt_ = lambda mask: np.bincount(gi[mask], minlength=len(mags) + 1)
    metrics = {"on_mag": np.ones(k.size, bool), "qc_fail": qcf[k], "duplicate": dup[k], "selected": sel[k],
               "selected_qc_fail": sel[k] & qcf[k], "selected_duplicate": sel[k] & dup[k]}
    metrics.update({x: status[k] == i for i, x in enumerate(STATUS)})
    metrics.update({"bam_" + x: bst[k] == i for i, x in enumerate(BAMST)})
    counts = {key: cnt_(msk) for key, msk in metrics.items()}
    per_mag = {}
    for g, m in enumerate(mags + [""]):
        if counts["on_mag"][g]:
            per_mag[m] = {key: int(v[g]) for key, v in counts.items() if v[g]}
    st.update(per_mag=per_mag,
              selections=[[m, b, c] for (m, b), c in sorted(sel_counts.items())],
              movies=movies[:50], n_movies=len(movies),
              movies_unresolved=sorted({movies[m] for m in mv_of[k][bst[k] == 1]})[:20],
              movies_ambiguous=sorted({movies[m] for m in mv_of[k][bst[k] == 2]})[:20],
              seconds=round(time.time() - t0, 1), reused=False)
    for nm in STATUS + ["off_mag"]:
        st["molecules_" + nm] = int((status[k] == STATUS.index(nm)).sum()) if nm != "off_mag" else st["molecules_off_mag"]
    st["molecules_on_mags"] = int(k.size)
    st["molecules_selected"] = int(sel.sum())
    write_atomic(stf, json.dumps(st))
    write_atomic(fpf, job["fp"])
    return st


# ============================================================================ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True); ap.add_argument("--genecat", required=True)
    ap.add_argument("--table", required=True); ap.add_argument("--map", default="")
    ap.add_argument("--binner", default=""); ap.add_argument("--binner-default", default="SB")
    ap.add_argument("--min-mapq", type=int, default=10); ap.add_argument("--exclude-qcfail", action="store_true")
    ap.add_argument("--threads", type=int, default=4); ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--abundance", default="")                   # default <GC>/Bin_<binner>/Annotation/Abundance/MGS.matL0.txt
    ap.add_argument("--version", default="")
    a = ap.parse_args()
    t_start = time.time()
    table = os.path.abspath(a.table)
    lib = re.sub(r"\.tsv$", "", table) + ".motlib"
    os.makedirs(lib, exist_ok=True)
    D = discover(a)
    for line in D["report"]:
        log(line)
    if D["fatal"]:
        mf4.stop(" | ".join(D["fatal"]))
    # GeneCat's MGS abundance (MGS x samples): the table's abundance column
    D["abundance_path"] = os.path.abspath(a.abundance) if a.abundance else \
        os.path.join(D["gc"], f"Bin_{D['binner']}", "Annotation", "Abundance", "MGS.matL0.txt")
    if a.abundance and not os.path.isfile(D["abundance_path"]):
        mf4.stop(f"--abundance {a.abundance}: not found")
    ab_lines = []
    D["abundance"] = {(s_, g_): v_ for s_, g_, v_ in mf4.read_abundance(
        D["abundance_path"], list(D["samples"]), sorted({i["MGS"] for i in D["mags"].values() if i["MGS"]}),
        ab_lines.append)}
    for line in ab_lines:
        log(line)
    if not os.path.isfile(D["abundance_path"]):
        D["abundance_path"] = ""
    samples, bams, mags = D["samples"], D["bams"], D["mags"]
    unresolved = []                                              # (sample, MAG, BAM, reason)

    # ---- genomes, one job per assembly
    gcdir = os.path.join(D["gc"], f"Bin_{D['binner']}", "Genomes", "MGS_ctg")
    by_asm = collections.defaultdict(dict)
    for m, i in mags.items():
        if i["status"] == "ok":
            by_asm[i["asm"]][m] = dict(contigs=i["contigs"], MGS=i["MGS"], file=i["file"])
    jobs = []
    for asm, ms in by_asm.items():
        fa = asm_fasta(asm)
        fp = md5s(METHOD, fa, stamp(fa), stamp(os.path.join(asm, "scaffolds.fasta.circ")),
                  md5s(*sorted(f"{m}\t{c}" for m, v in ms.items() for c in v["contigs"])),
                  md5s(*sorted(f"{m}\t{v['file']}" for m, v in ms.items())),
                  md5s(*(stamp(os.path.join(gcdir, mf4.rep_fasta_name(v["MGS"], m))) for m, v in sorted(ms.items())
                         if v["MGS"])))
        jobs.append(dict(asm=asm, mags=ms, lib=lib, gcdir=gcdir, fp=fp))
    log(f"genomes: {sum(map(len, by_asm.values()))} MAG FASTAs from {len(by_asm)} assemblies")
    contig_rows, asm_fp = [], {}
    with cf.ProcessPoolExecutor(max(1, min(a.threads, len(jobs) or 1))) as ex:
        for job, res in zip(jobs, ex.map(genome_job, jobs)):
            asm_fp[job["asm"]] = job["fp"]
            for m, r in res["mags"].items():
                mags[m].update(genome=r, fasta=os.path.join(lib, "genomes", mags[m]["file"] + ".fa.gz")
                               if r["status"] == "ok" else "")
                if r["status"] != "ok":
                    mags[m]["status"] = r["status"]
            contig_rows += res["contigs"]
    for m, i in mags.items():
        i.setdefault("genome", {}); i.setdefault("fasta", "")
    genecat_checked = sum(bool(i["genome"].get("genecat_fasta")) for i in mags.values())
    bad_rules = [m for m, i in mags.items() if i["status"] == "contigs_differ_from_genecat"]
    if bad_rules:
        mf4.stop(f"{len(bad_rules)} MAG(s) whose bin-file contigs differ from GeneCat's own FASTA - the naming rules do "
                 f"not match this GeneCat version: {', '.join(bad_rules[:5])}")
    st = collections.Counter(i["status"] for i in mags.values())
    log("MAGs: " + "  ".join(f"{k} {v}" for k, v in sorted(st.items()))
        + f" ; {genecat_checked} checked base by base against GeneCat's FASTA")

    # ---- samples, one job per CRAM
    ready_bam = {p: b["mot_ready"] for p, b in bams.items()}
    sjobs = []
    asm_mags = collections.defaultdict(list)
    for m, i in mags.items():
        if i["asm"] and i.get("contigs"):
            asm_mags[i["asm"]].append(m)
    ctg_ref = {}
    for m, c, L, h, _ in contig_rows:
        ctg_ref[(mags[m]["asm"], c)] = (m, L, h)
    for s, v in samples.items():
        if not v["has_cram"]:
            v["status"] = "no_cram"; continue
        if not v["asm"]:
            v["status"] = "no_assembly"; continue
        ml = sorted(asm_mags.get(v["asm"], []))
        if not ml:
            v["status"] = "no_final_MAG_in_assembly"; continue
        ref = {}
        for idx, m in enumerate(ml):
            for c in mags[m]["contigs"]:
                if (v["asm"], c) in ctg_ref:
                    ref[c] = (SHARED if c in mags[m].get("shared", []) else idx,) + ctg_ref[(v["asm"], c)][1:]
                else:          # a MAG not exported: its reads are still classified, never selected
                    ref.setdefault(c, (SHARED if c in mags[m].get("shared", []) else idx, None, None))
        bl = [os.path.realpath(p) for p in v["bams"]]
        movie_bam = {}
        for bi, p in enumerate(bl):
            for mv in bams[p]["movies"]:
                movie_bam[mv] = bi if movie_bam.get(mv, bi) == bi else -1
        labels, used = [], set()
        for p in bl:
            lb = safe(re.sub(r"\.bam$", "", os.path.basename(p)))
            lb = lb if lb not in used else f"{lb}_{md5s(p)[:6]}"
            used.add(lb); labels.append(lb)
        fp = md5s(METHOD, v["cram"], stamp(v["cram"]), asm_fp.get(v["asm"], ""), a.min_mapq, a.exclude_qcfail,
                  *[f"{p}:{stamp(p)}:{ready_bam[p]}:{','.join(bams[p]['movies'])}" for p in bl],
                  *[f"{m}:{mags[m]['file']}:{mags[m]['status']}" for m in ml])
        sjobs.append(dict(sample=s, lib=lib, cram=v["cram"], fp=fp, mags=ml, ref=ref, movie_bam=movie_bam,
                          bam_ready=[ready_bam[p] for p in bl], bam_paths=[bams[p]["path"] for p in bl],
                          bam_label=labels, mag_ok={m: mags[m]["status"] == "ok" for m in ml},
                          mag_file={m: mags[m]["file"] for m in ml}, min_mapq=a.min_mapq,
                          exclude_qcfail=a.exclude_qcfail, samtools=a.samtools))
    log(f"reads: {len(sjobs)} CRAMs to read ({sum(1 for v in samples.values() if v.get('status'))} samples without "
        f"a usable CRAM / assembly / MAG)")
    stats = {}
    with cf.ProcessPoolExecutor(max(1, min(a.threads, len(sjobs) or 1))) as ex:
        futs = {ex.submit(sample_job, j): j for j in sjobs}
        for fu in cf.as_completed(futs):
            j = futs[fu]
            try:
                stats[j["sample"]] = fu.result()
                s_ = stats[j["sample"]]
                log(f"  {j['sample']}: {s_['molecules_on_mags']} molecules on MAGs, {s_['molecules_selected']} selected"
                    + (" (reused)" if s_.get("reused") else f" ({s_['seconds']} s)"))
            except Exception as e:
                samples[j["sample"]]["status"] = f"cram_failed: {str(e)[:200]}"
                log(f"  {j['sample']}: FAILED {e}")
    for s, v in samples.items():
        v.setdefault("status", "ok")
    finish(a, D, lib, table, sjobs, stats, contig_rows, unresolved, t_start)


def finish(a, D, lib, table, sjobs, stats, contig_rows, unresolved, t_start):
    samples, bams, mags = D["samples"], D["bams"], D["mags"]
    jobs = {j["sample"]: j for j in sjobs}
    rows, amag = [], []                                         # mot rows ; assignment summary rows
    ab = D.get("abundance", {})                                 # (SmplID, MGS) -> GeneCat's value (> 0)
    for s, v in samples.items():
        if v["status"] != "ok":
            if v["status"] != "no_final_MAG_in_assembly":           # an assembly without final MAGs: nothing to export
                ms = sorted(m for m, i in mags.items() if v["asm"] and i["asm"] == v["asm"])
                unresolved += [(s, m, "", v["status"]) for m in ms] or [(s, "", "", v["status"])]
            continue
        st, j = stats[s], jobs[s]
        sel = collections.defaultdict(dict)
        for m, b, c in st["selections"]:
            sel[m][b] = c
        # mot fits one adapter per sample (all its genomes, one technical state): one chemistry per sample
        chems = {kits(bams[os.path.realpath(j["bam_paths"][b])]["chemistry"]) for m in j["mags"] for b in sel.get(m, {})}
        for m in j["mags"]:
            c = st["per_mag"].get(m, {})
            out = sel.get(m, {})
            why = ("" if mags[m]["status"] == "ok" else f"MAG {mags[m]['status']}") or \
                  ("several chemistries in one sample (mot fits one adapter per sample; give each chemistry's BAMs "
                   "a sample of their own)" if len(chems) > 1 else "")
            mgs = mags[m]["MGS"]
            for b, n in sorted(out.items()):
                bam = j["bam_paths"][b]
                reads = os.path.join(lib, "assignments", safe(s), mags[m]["file"], j["bam_label"][b] + ".reads.txt.gz")
                r = dict(sample=mot_name(s), genome=mot_name(mags[m]["file"]), bam=bam, reads=reads,
                         ref=mags[m]["fasta"], mgs=mgs or ".", abundance=ab.get((s, mgs), ".") if mgs else ".",
                         role="native", _s=s, _m=m, _n=n)
                if why:
                    unresolved.append((s, m, bam, why))
                else:
                    rows.append(r)
            for b, bam in enumerate(j["bam_paths"] or [""]):
                bm = os.path.realpath(bam) if bam else ""
                amag.append((s, m, mags[m]["MGS"], bam, c.get("on_mag", 0), c.get("unique", 0),
                             c.get("ambiguous_low_mapq", 0), c.get("split_across_mags", 0), c.get("split_mag_unbinned", 0),
                             c.get("contig_in_two_mags", 0), c.get("reference_mismatch", 0), c.get("no_primary", 0),
                             c.get("qc_fail", 0), c.get("duplicate", 0), out.get(b, 0) if bam else 0,
                             int(bool(out.get(b)) and not why),
                             why or ("" if out.get(b) else "no selected molecule" if bm and bams[bm]["mot_ready"]
                                     else (bams[bm]["why"] if bm else "no raw BAM"))))
        # molecules on MAGs whose BAM could not be named
        unres = sum(c.get("bam_movie_not_in_bams", 0) + c.get("bam_movie_in_two_bams", 0) for c in st["per_mag"].values())
        if unres:
            unresolved.append((s, "", "", f"{unres} molecule(s) on MAGs whose movie is in none / two of the sample's "
                                          f"BAMs (movies {', '.join(st['movies_unresolved'] + st['movies_ambiguous'])[:200]})"))
        if st["mag_contigs_mismatch"] or st["mag_contigs_not_in_cram"]:
            unresolved.append((s, "", "", f"CRAM reference differs from the exported sequence: {st['mag_contigs_mismatch']} "
                                          f"MAG contig(s) with another length/M5, {st['mag_contigs_not_in_cram']} not in "
                                          f"the CRAM ({'; '.join(st['mismatch_examples'][:2])})"))
        for p in v["bams"]:
            b = bams[os.path.realpath(p)]
            if not b["mot_ready"]:
                n_on = st["molecules_on_mags"]
                if st["non_pacbio_names"] < st["records"] and n_on:
                    unresolved.append((s, "", b["path"], f"raw BAM not usable by mot: {b['why']}"))
        if not v["bams"] and st["molecules_on_mags"] and st["non_pacbio_names"] < st["records"]:
            unresolved.append((s, "", "", "PacBio reads in the CRAM but no raw BAM found (map / input_raw.txt)"))
    # MAGs that could not be exported (a MAG of a sample outside this run is not this run's to export)
    for m, i in mags.items():
        if i["status"] not in ("ok", "sample_not_in_map"):
            unresolved.append((i["sample"], m, "", f"MAG not exported: {i['status']}"))
    # ---- one BAM in two samples: their selections must be disjoint (mot refuses a molecule in two samples)
    by_bam = collections.defaultdict(list)
    for r in rows:
        by_bam[os.path.realpath(r["bam"])].append(r)
    drop = set()
    for bam, rs in by_bam.items():
        if len({r["_s"] for r in rs}) < 2:
            continue
        sets = {id(r): set(gzip.open(r["reads"], "rt").read().split()) for r in rs}
        for x in range(len(rs)):
            for y in range(x + 1, len(rs)):
                if rs[x]["_s"] != rs[y]["_s"] and sets[id(rs[x])] & sets[id(rs[y])]:
                    drop.update((id(rs[x]), id(rs[y])))
    for r in rows:
        if id(r) in drop:
            unresolved.append((r["_s"], r["_m"], r["bam"], "BAM shared with another sample and the same molecules "
                                                           "selected for both (MATAFILER4 mapped them twice)"))
    rows = [r for r in rows if id(r) not in drop]
    rows.sort(key=lambda r: (r["sample"], r["genome"], r["bam"]))

    # ---- validation of what is written
    checks = []
    def check(name, ok, detail=""):
        checks.append(dict(check=name, ok=bool(ok), detail=detail))
    names = collections.Counter((r["sample"], r["genome"], r["bam"]) for r in rows)
    bad = [x for r in rows for x in (r["sample"], r["genome"]) if not MOT_NAME.match(x)]
    check("mot sample and genome names valid (mot rule: letters, digits . _ -, <= 100)", not bad, ", ".join(bad)[:300])
    check("one row per sample x genome x BAM", all(k == 1 for k in names.values()))
    refs, srcs, gsrc, chem_of = (collections.defaultdict(set) for _ in range(4))
    for r in rows:
        refs[(r["sample"], r["genome"])].add((r["ref"], r["mgs"], r["abundance"]))
        srcs[r["sample"]].add(r["_s"])
        gsrc[r["genome"]].add(r["_m"])
        chem_of[r["sample"]].add(kits(bams[os.path.realpath(r["bam"])]["chemistry"]))
    check("one reference, mgs and abundance per sample x genome", all(len(v) == 1 for v in refs.values()))
    check("one MATAFILER4 sample per mot sample name, one MAG per genome name",
          all(len(v) == 1 for v in srcs.values()) and all(len(v) == 1 for v in gsrc.values()))
    check("one chemistry per mot sample (mot's rule)", all(len(v) == 1 for v in chem_of.values()))
    check("abundance a number >= 0 or '.'", all(r["abundance"] == "." or num(r["abundance"]) >= 0 for r in rows),
          ", ".join(sorted({r["abundance"] for r in rows if r["abundance"] != "." and not num(r["abundance"]) >= 0}))[:300])
    # every selection file, one BAM at a time (memory: one BAM's molecules)
    bad_files, bad_sel, overlap = [], [], []
    rows_of = collections.defaultdict(list)
    for r in rows:
        rows_of[os.path.realpath(r["bam"])].append(r)
    for bam, rs in rows_of.items():
        mv_ok, seen_ = set(bams[bam]["movies"]), set()
        for r in rs:
            for f in (r["bam"], r["ref"], r["ref"] + ".fai", r["reads"]):
                if not os.path.isfile(f):
                    bad_files.append(f)
            if not os.path.isfile(r["reads"]):
                continue
            lines = gzip.open(r["reads"], "rt").read().split("\n")[:-1]
            names_ = set(lines)
            parsed = [x.split("/") for x in lines]
            prob = ("empty" if not lines else "not unique" if len(names_) != len(lines) else
                    "a name is not movie/zmw" if any(len(p) != 2 or not p[1].isdigit() for p in parsed) else
                    "a movie is not in the BAM" if any(p[0] not in mv_ok for p in parsed) else
                    f"{len(lines)} names but the catalog selects {r['_n']}" if len(lines) != r["_n"] else "")
            if prob:
                bad_sel.append(f"{r['reads']}: {prob}")
            if seen_ & names_:
                overlap.append(r["reads"])
            seen_ |= names_
    check("every bam / ref / ref.fai / reads file exists", not bad_files, ", ".join(bad_files[:5]))
    check("every selection: non-empty, unique, movie/zmw names, movies in its BAM, count = catalog", not bad_sel,
          "; ".join(bad_sel[:5]))
    check("no molecule in two rows of one BAM (mot's rule)", not overlap, ", ".join(overlap[:3]))
    # the catalog agrees with the selection files on disk (every selected molecule is in exactly one file)
    import pyarrow.parquet as pq
    diff = []
    for s in stats:
        t = pq.read_table(os.path.join(lib, "catalog", "molecules", f"{safe(s)}.parquet"), columns=["selected"])
        n_cat = int(sum(t.column("selected").to_pylist()))
        d = os.path.join(lib, "assignments", safe(s))
        n_files = sum(sum(1 for _ in gzip.open(os.path.join(r, f), "rt")) for r, _, fs in os.walk(d) for f in fs
                      if f.endswith(".reads.txt.gz"))
        if n_cat != n_files:
            diff.append(f"{s}: catalog {n_cat}, files {n_files}")
    check("selected molecules in the catalog = molecules in the selection files", not diff, "; ".join(diff[:5]))
    cram_counts = {s: (st["molecules_on_mags"], sum(st["molecules_" + x] for x in STATUS)) for s, st in stats.items()}
    check("every molecule on a MAG has one status (counts add up)", all(x == y for x, y in cram_counts.values()))
    ok_mags = sum(i["status"] == "ok" for i in mags.values())
    check("at least one MAG exported", ok_mags > 0)
    check("at least one mot row", bool(rows))

    # ---- write the table, catalog, qc, manifest
    tsv(table, MOT_COLUMNS, [tuple(r[k] for k in MOT_COLUMNS) for r in rows])
    cat = os.path.join(lib, "catalog")
    tsv(os.path.join(cat, "mags.tsv"),
        ["MAG", "file", "MGS", "is_MGS_rep", "status", "sample", "assembly", "n_contigs", "bp", "n_circular",
         "completeness", "contamination", "fasta", "md5", "checked_vs_genecat", "notes"],
        [(m, i["file"], i["MGS"], i["is_rep"], i["status"], i["sample"], i["asm"], i["genome"].get("n_contigs", len(i["contigs"])),
          i["genome"].get("bp", ""), i["genome"].get("n_circular", ""), i["completeness"], i["contamination"], i["fasta"],
          i["genome"].get("md5", ""), int(bool(i["genome"].get("genecat_fasta"))),
          "; ".join(x for x in (f"missing {i['genome']['n_missing']} contig(s): {','.join(i['genome']['missing'])}"
                                if i["genome"].get("n_missing") else "",
                                f"contig(s) also in another MAG: {','.join(i['shared'][:3])}" if i.get("shared") else "",
                                f".circ sequence differs: {','.join(i['genome']['circ_differs'])}"
                                if i["genome"].get("circ_differs") else "") if x))
         for m, i in mags.items()])
    with gzip.open(os.path.join(cat, "contigs.tsv.gz.tmp"), "wt", compresslevel=6) as o:
        o.write("MAG\tcontig\tassembly\tlength\tmd5\tcircular\n")
        for m, c, L, h, circ in sorted(contig_rows):
            o.write(f"{m}\t{c}\t{mags[m]['asm']}\t{L}\t{h}\t{circ}\n")
    os.replace(os.path.join(cat, "contigs.tsv.gz.tmp"), os.path.join(cat, "contigs.tsv.gz"))
    tsv(os.path.join(cat, "samples.tsv"),
        ["sample", "status", "assembly_group", "assembly", "cram", "n_bams", "bam_provenance", "co_assembled_with",
         "records", "non_pacbio_names", "qcfail_records", "dup_records", "retained_alignments", "molecules",
         "molecules_on_mags", "molecules_selected", "notes"],
        [(s, v["status"], v["group"], v["asm"], v["cram"], len(v["bams"]), v["provenance"],
          ",".join(x for x, w in samples.items() if w["asm"] == v["asm"] and x != s and v["asm"]),
          stats.get(s, {}).get("records", ""), stats.get(s, {}).get("non_pacbio_names", ""),
          stats.get(s, {}).get("qcfail_records", ""), stats.get(s, {}).get("dup_records", ""),
          retained(v, stats.get(s)), stats.get(s, {}).get("molecules", ""), stats.get(s, {}).get("molecules_on_mags", ""),
          stats.get(s, {}).get("molecules_selected", ""),
          "; ".join(v["notes"] + ([f"other inputs: {len(v['other_inputs'])} (e.g. {v['other_inputs'][0]})"]
                                  if v["other_inputs"] else [])))
         for s, v in samples.items()])
    tsv(os.path.join(cat, "bams.tsv"),
        ["bam", "samples", "opens", "movies", "chemistry", "read_type", "kinetics", "mot_ready", "why"],
        [(b["path"], ",".join(b["samples"]), int(b["ok"]), ",".join(b["movies"]), b["chemistry"], b["read_types"],
          int(b["kinetics"]), int(b["mot_ready"]), b["why"]) for b in bams.values()])
    q = os.path.join(lib, "qc")
    tsv(os.path.join(q, "assignment_summary.tsv"),
        ["sample", "MAG", "MGS", "bam", "molecules_on_MAG", "unique", "ambiguous_low_mapq", "split_across_mags",
         "split_mag_unbinned", "contig_in_two_mags", "reference_mismatch", "no_primary", "qc_fail", "duplicate",
         "selected", "mot_row", "why_no_row"], amag)
    tsv(os.path.join(q, "unresolved.tsv"), ["sample", "MAG", "bam", "reason"], sorted(set(unresolved)))
    # ---- totals + status
    T = collections.Counter()
    for st in stats.values():
        for k in ("records", "primary", "supplementary", "secondary", "qcfail_records", "dup_records", "molecules",
                  "molecules_on_mags", "molecules_off_mag", "molecules_selected") + tuple("molecules_" + x for x in STATUS):
            T[k] += st.get(k, 0)
        for c in st["per_mag"].values():
            T["qc_fail_on_mags"] += c.get("qc_fail", 0); T["duplicate_on_mags"] += c.get("duplicate", 0)
            T["selected_qc_fail"] += c.get("selected_qc_fail", 0); T["selected_duplicate"] += c.get("selected_duplicate", 0)
            T["bam_unresolved"] += c.get("bam_movie_not_in_bams", 0) + c.get("bam_movie_in_two_bams", 0)
    hifi = [s for s, v in samples.items() if any(bams[os.path.realpath(p)]["mot_ready"] for p in v["bams"])]
    no_ret = [s for s, v in samples.items() if s in stats and retained(v, stats[s]) == "no evidence"]
    failed = [c for c in checks if not c["ok"]]
    status = "FAILED" if failed else "PARTIAL" if unresolved else "COMPLETE"
    combos = len(amag)
    summary = [
        ("MF4 run", D["run"]), ("GeneCat", D["gc"]), ("binner", D["binner"]), ("maps", ", ".join(D["maps"])),
        ("samples discovered", len(samples)), ("samples with mot-ready HiFi BAMs", len(hifi)),
        ("samples whose CRAM was read", len(stats)),
        ("MAGs discovered", len(mags)), ("MAGs exported", ok_mags),
        ("MAG FASTAs checked vs GeneCat", sum(bool(i["genome"].get("genecat_fasta")) and i["status"] == "ok" for i in mags.values())),
        ("raw BAMs", len(bams)), ("raw BAMs mot-ready", sum(b["mot_ready"] for b in bams.values())),
        ("CRAMs processed", len(stats)),
        ("CRAM records (mapped)", T["records"]), ("molecules (mapped)", T["molecules"]),
        ("molecules on final MAGs", T["molecules_on_mags"]), ("molecules off MAGs", T["molecules_off_mag"]),
        ("  of those on MAGs: QC-passing (no 0x200)", T["molecules_on_mags"] - T["qc_fail_on_mags"]),
        ("  QC-failed (0x200)", T["qc_fail_on_mags"]), ("  duplicate-flagged (0x400)", T["duplicate_on_mags"]),
        ("unique MAG assignments", T["molecules_unique"]),
        ("ambiguous: low MAPQ", T["molecules_ambiguous_low_mapq"]),
        ("ambiguous: split across MAGs", T["molecules_split_across_mags"]),
        ("ambiguous: MAG + unbinned", T["molecules_split_mag_unbinned"]),
        ("ambiguous: contig in two MAGs", T["molecules_contig_in_two_mags"]),
        ("reference mismatch", T["molecules_reference_mismatch"]), ("no primary", T["molecules_no_primary"]),
        ("unresolved BAM (movie)", T["bam_unresolved"]),
        ("selected for mot", T["molecules_selected"]), ("  of which QC-failed / duplicate", f"{T['selected_qc_fail']} / {T['selected_duplicate']}"),
        ("sample/MAG/BAM combinations", combos), ("mot rows", len(rows)),
        ("mot samples", len({r["sample"] for r in rows})),
        ("mot units (sample x genome)", len({(r["sample"], r["genome"]) for r in rows})),
        ("  with an MGS / with GeneCat abundance", f"{len({(r['sample'], r['genome']) for r in rows if r['mgs'] != '.'})} / "
                                                   f"{len({(r['sample'], r['genome']) for r in rows if r['abundance'] != '.'})}"),
        ("GeneCat abundance table", D.get("abundance_path") or "none"),
        ("min MAPQ", a.min_mapq), ("QC-failed molecules selected", "no" if a.exclude_qcfail else "yes"),
        ("CRAMs without retained-alignment evidence", len(no_ret)),
        ("unresolved items", len(set(unresolved))), ("status", status),
    ]
    tsv(os.path.join(q, "export_summary.tsv"), ["item", "value"], summary)
    write_atomic(os.path.join(q, "validation.json"), json.dumps(dict(status=status, checks=checks,
                                                                     unresolved=len(set(unresolved))), indent=1))
    manifest = dict(
        tool="mmt export-mot", mmt_version=a.version, method=METHOD, created=datetime.datetime.now().isoformat(timespec="seconds"),
        table=table, library=lib, status=status,
        inputs=dict(mf4_run=D["run"], genecat=D["gc"], binner=D["binner"], maps=D["maps"], maps_from=D["map_src"],
                    magvsgc=os.path.join(D["gc"], f"Bin_{D['binner']}", "MAGvsGC.txt.gz"),
                    magvsgc_stamp=stamp(os.path.join(D["gc"], f"Bin_{D['binner']}", "MAGvsGC.txt.gz")),
                    abundance=D.get("abundance_path", ""), abundance_stamp=stamp(D.get("abundance_path", ""))),
        table_format=dict(columns=MOT_COLUMNS, mot="0.7 input table (mot help input)",
                          abundance="GeneCat's MGS abundance of the MAG's MGS in the sample; '.' = no MGS, or no value "
                                    "> 0 for the pair in GeneCat's table"),
        settings=dict(min_mapq=a.min_mapq, exclude_qcfail=a.exclude_qcfail, threads=a.threads),
        assignment=dict(scope="competitive within the (co-)assembly each sample's reads were mapped to by MATAFILER4 "
                              "(<run>/<S>/mapping/<S>-smd.cram); not across the cohort",
                        statuses={k: v for k, v in STATUS_DOC.items()}, selected=SELECTED_DOC),
        counts={k: v for k, v in summary},
        crams={s: dict(cram=st["cram"], stamp=stamp(st["cram"]), retained_alignments=retained(samples[s], st),
                       has_m5=st["cram_has_m5"]) for s, st in stats.items()},
        seconds=round(time.time() - t_start, 1))
    write_atomic(os.path.join(lib, "manifest.json"), json.dumps(manifest, indent=1))

    # ---- terminal report
    w = max(len(k) for k, _ in summary) + 2
    print("\nMMT - export MATAFILER4 cohort to Methyl-o-Type\n")
    for k, v in summary:
        print(f"  {k + ':':<{w}}{v}")
    print(f"\n  output:   {table}\n  library:  {lib}/")
    if no_ret:
        print(f"\n  NOTE: {len(no_ret)} CRAM(s) show no retained (0x200 / 0x400) alignments and no .retained.size file: "
              f"made without -retainAllCRAM, or nothing failed MATAFILER4's filter (e.g. {no_ret[0]}). "
              "Rejected alignments that MATAFILER4 discarded cannot be recovered.")
    print("  NOTE: assignment is competitive within each sample's own (co-)assembly only - MATAFILER4 mapped every sample "
          "to its own assembly, so the CRAMs cannot tell between MAGs of different assemblies.")
    if failed:
        print("\n  FAILED checks:\n" + "\n".join(f"    - {c['check']}: {c['detail']}" for c in failed))
    if unresolved:
        u = sorted(set(unresolved))
        print(f"\n  not exported ({len(u)}; all in qc/unresolved.tsv):")
        for x in u[:12]:
            print("    - " + " | ".join(y for y in x if y))
        if len(u) > 12:
            print(f"    ... {len(u) - 12} more")
    print(f"\nStatus: {status}")
    sys.exit({"COMPLETE": 0, "PARTIAL": 3, "FAILED": 1}[status])


def retained(v, st):
    if st is None:
        return ""
    if os.path.isfile(v["cram"] + ".retained.size"):
        return "yes (.retained.size)"
    if st.get("qcfail_records") or st.get("dup_records"):
        return "yes (flagged records)"
    return "no evidence"


STATUS_DOC = {
    "unique": "primary + supplementary records all on contigs of one final MAG, primary MAPQ >= min MAPQ",
    "ambiguous_low_mapq": "one MAG, but primary MAPQ < min MAPQ: a near-equal placement elsewhere in the assembly",
    "split_across_mags": "primary / supplementary pieces on two or more final MAGs",
    "split_mag_unbinned": "pieces on one final MAG and on contigs of no final MAG",
    "contig_in_two_mags": "on a contig that two final MAGs of the assembly claim",
    "reference_mismatch": "on a contig whose CRAM header length / M5 differs from the exported sequence",
    "no_primary": "only supplementary records of the molecule are in the CRAM",
}
SELECTED_DOC = ("status unique, its movie in exactly one of the sample's raw BAMs, that BAM mot-ready (opens, kinetics, "
                "READTYPE CCS), the MAG exported; QC-failed (0x200) and duplicate-flagged (0x400) molecules included "
                "unless --exclude-qcfail")


if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:                     # mf4.stop("ERROR: ...") = the export cannot be made at all
        if isinstance(e.code, str):
            print(e.code, file=sys.stderr)
            print("Status: FAILED")
            sys.exit(1)
        raise
