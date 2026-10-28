#!/usr/bin/env python3
"""MMT derep-ctg: each MGS as a short list of collinear blocks ("dereplicated contigs").

  prep     once per MMT index: catalogue members + prodigal genes of every MGS-MAG contig    -> PREP/
  plan     which MGS need (re)computing                                                      -> stdout
  mgs      one MGS: families -> anchors -> blocks -> instances -> representatives -> links -> lift  -> OUT/<m>/
  summary  OUT/summary.tsv from every MGS folder

A block is an ordered run of gene families that sits in the same order and orientation in the MAGs of an MGS.
Method and reasons: docs/derep-ctg.md. The rules follow Panaroo (context-aware merging of gene variants,
low-support and contig-end filtering, paralog care) and PPanGGOLiN (multi-copy families bridged, persistent backbone
vs variable regions, incomplete-aware comparison at contig ends), on MATAFILER4's own cohort-wide gene catalogue.
"""
import argparse, bisect, collections, csv, gzip, hashlib, itertools, json, os, re, shutil, statistics, subprocess, sys

csv.field_size_limit(1 << 30)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
RE_ID, RE_PART = re.compile(r"(?:^|;)ID=\d+_(\d+)"), re.compile(r"partial=(\d\d)")
RE_CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CODE = dict(zip((a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"), "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def die(msg): sys.exit(f"ERROR: {msg}")
def log(msg): print(msg, file=sys.stderr, flush=True)
def opener(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)
def rc(s): return s.translate(COMP)[::-1]
def md5(s, n=16): return hashlib.md5(s.encode()).hexdigest()[:n]
def first_file(*ps): return next((p for p in ps if p and os.path.isfile(p) and os.path.getsize(p) > 0), None)
def tsv(path):
    with opener(path) as fh: return list(csv.DictReader(fh, delimiter="\t"))
def num(x, d=float("nan")):   # "NA", "", "nan" -> d
    try: v = float(x)
    except (TypeError, ValueError): return d
    return d if v != v else v


# column types by name, so that every table has the same schema in every MGS (also when empty): string, float, bool,
# int32, int8 ; any other column is int64 (a `types` argument overrides)
PQ_TYPES = dict(str="MGS MAG contig block block_a block_b tier samples strand start_status end_status orientation catalogue rep_gene "
                    "rep_MAG rep_contig rep_strand role gene cls path MAGs seq nb5_key nb3_key",
                float="frac identity aligned_frac prevalence prevalence_adj median_identity shared_frac completeness contamination "
                      "frac_core frac_blocks gc gc_skew",
                bool="complete is_rep circular rep_complete rep_is_MGS_rep partial terminal is_MGS_rep mixed dup_in_mag",
                int32="slot slot_a slot_b slot_after instance host_instance host_slot nb5_instance nb3_instance c_start c_end b_start "
                      "w_start w_end first_slot last_slot",
                int16="a c g t", int8="end_a end_b ori")
DELTA = {"c_start", "c_end", "b_start", "w_start", "w_end", "start", "end"}   # positions: delta encoding, no dictionary


def pq_type(c, types=None):
    import pyarrow as pa
    by = {x: getattr(pa, {"str": "string", "float": "float64", "bool": "bool_"}.get(k, k))() for k, v in PQ_TYPES.items() for x in v.split()}
    return (types or {}).get(c) or by.get(c, pa.int64())


def pq_writer(path, schema):
    """parquet writer with MMT's encodings: positions delta-encoded, sequences plain, the rest dictionary-encoded"""
    import pyarrow as pa, pyarrow.parquet as pq
    d = [f.name for f in schema if f.name in DELTA and pa.types.is_integer(f.type)]
    try: return pq.ParquetWriter(path, schema, compression="zstd", use_dictionary=[f.name for f in schema if f.name not in d and f.name != "seq"],
                                 column_encoding={c: "DELTA_BINARY_PACKED" for c in d})
    except (TypeError, ValueError, NotImplementedError, pa.ArrowNotImplementedError): return pq.ParquetWriter(path, schema, compression="zstd")


def write_pq(cols, path, types=None):
    """dict of equal-length lists -> zstd parquet (temp name, then rename)"""
    import pyarrow as pa
    t = pa.table({k: pa.array(v, type=pq_type(k, types)) for k, v in cols.items()})
    tmp = f"{path}.{os.getpid()}.tmp"
    w = pq_writer(tmp, t.schema)
    w.write_table(t); w.close()
    os.replace(tmp, path)


def rel(p, root):   # paths inside the prep folder are stored relative to it (the folder is built elsewhere, then renamed)
    return os.path.relpath(p, root) if os.path.abspath(p).startswith(os.path.abspath(root) + os.sep) else p


def write_text(path, text):
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "w") as o: o.write(text)
    os.replace(tmp, path)


# ============================================================================================== prep
def catalogue_files(gcd, gc_id):
    idx = next((p for t in ("clstr", "ctsv") for p in [first_file(f"{gcd}/compl.incompl.{gc_id}.fna.{t}.idx",
                                                                  f"{gcd}/compl.incompl.{gc_id}.fna.{t}.idx.gz")] if p), None)
    prot = first_file(f"{gcd}/compl.incompl.{gc_id}.prot.faa", f"{gcd}/compl.incompl.{gc_id}.prot.faa.gz")
    if not idx: die(f"no gene catalogue cluster index compl.incompl.{gc_id}.fna.clstr.idx in {gcd} (DEREP_GC_ID?)")
    if not prot: die(f"no catalogue proteins compl.incompl.{gc_id}.prot.faa in {gcd}")
    return idx, prot


def fasta_access(fa, own, samtools, need=()):
    """(fasta, fai, gzi) for pysam, never writing next to MATAFILER4's / GeneCat's files: a current .fai holding every
    name in `need` is used as is ; else MMT's own index own.fai (+ own.gzi for bgzip) ; plain gzip -> own copy own.fa
    (of the `need` records only, when given)"""
    fai = fa + ".fai"
    if (not fa.endswith(".gz") and os.path.isfile(fai) and os.path.getmtime(fai) >= os.path.getmtime(fa)
            and (not need or set(need) <= {l.split("\t", 1)[0] for l in open(fai)})):
        return fa, fai, ""
    if fa.endswith(".gz"):
        if subprocess.run([samtools, "faidx", "--fai-idx", own + ".fai", "--gzi-idx", own + ".gzi", fa], capture_output=True).returncode == 0:
            return fa, own + ".fai", own + ".gzi"
        if need:                                    # only the MAG contigs: never a full copy of every assembly
            with gzip.open(fa, "rt") as i, open(own + ".fa", "w") as o:
                keep = False
                for l in i:
                    if l[0] == ">": keep = l[1:].split(None, 1)[0] in need
                    if keep: o.write(l)
        else:
            with gzip.open(fa, "rb") as i, open(own + ".fa", "wb") as o: shutil.copyfileobj(i, o, 1 << 24)
        fa = own + ".fa"
    r = subprocess.run([samtools, "faidx", "--fai-idx", own + ".fai", fa], capture_output=True, text=True)
    if r.returncode: die(f"samtools faidx {fa}: {r.stderr.strip()[:300]}")
    return fa, own + ".fai", ""


def fasta_open(fa, fai, gzi=""):
    import pysam
    return pysam.FastaFile(fa, filepath_index=fai, filepath_index_compressed=gzi or None)


def protein_file(P, cats, cat_f, tmp):
    """the catalogue proteins of `cats` through an index of these names only (htslib loads a .fai whole: ~4 GB and
    ~45 s for 15 M catalogue proteins in every process otherwise)"""
    import pyarrow.parquet as pq
    T = pq.read_table(os.path.join(P, cat_f["proteins_index"]), filters=[("name", "in", sorted({str(c) for c in cats}))]).to_pydict()
    fai = f"{tmp}.prot.fai"
    with open(fai, "w") as o: o.writelines("\t".join(map(str, r)) + "\n" for r in zip(*T.values()))
    try: return fasta_open(os.path.join(P, cat_f["proteins"]), fai, os.path.join(P, cat_f["proteins_gzi"]) if cat_f["proteins_gzi"] else "")
    finally: os.remove(fai)


def translate(q): return "".join(CODE.get(q[i:i + 3], "X") for i in range(0, len(q) - 2, 3)).rstrip("*")
def kmers(p, k=5): return {p[i:i + k] for i in range(len(p) - k + 1)}


def prep_assembly(job):
    """one assembly: genes on its MGS-MAG contigs (prodigal GFF + catalogue member), checked -> genes/<key>.parquet"""
    import pyarrow as pa
    key, asm, host, P, check_n, samtools = job
    gff = first_file(*(f"{asm}/genePred/genes{b}.gff{z}" for b in ("", ".bac") for z in ("", ".gz")))
    fna = first_file(*(f"{asm}/genePred/genes{b}.shrtHD.fna{z}" for b in ("", ".bac") for z in ("", ".gz")))
    fa = first_file(f"{asm}/scaffolds.fasta.filt", f"{asm}/scaffolds.fasta.filt.gz")
    if not (gff and fa): return key, f"no genePred/genes.gff or scaffolds.fasta.filt in {asm}", None
    fasta, fai, gzi = fasta_access(fa, f"{P}/fasta/{key}", samtools, host)
    lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
    miss = [c for c in host if c not in lens]
    if miss: return key, f"{len(miss)} MAG contigs not in {fa} (e.g. {miss[0]})", None
    if len(lens) > len(host):           # MGS jobs load a .fai whole: an index of the MAG contigs only
        with open(fai) as i, open(f"{P}/fasta/{key}.mag.fai", "w") as o: o.writelines(l for l in i if l.split("\t", 1)[0] in host)
        if fai == f"{P}/fasta/{key}.fai": os.remove(fai)
        fai = f"{P}/fasta/{key}.mag.fai"
    mem = {}
    for l in open(f"{P}/members/{key}.tsv"):
        c, n, i = l.rstrip("\n").split("\t")
        if mem.setdefault((c, int(n)), int(i)) != int(i): return key, f"gene {c}_{n} in two catalogue clusters", None
    rows, first, seen = [], [], set()
    with opener(gff) as fh:
        for l in fh:
            if l[0] == "#": continue
            x = l.rstrip("\n").split("\t")
            if len(x) < 9 or x[2] != "CDS": continue
            m = RE_ID.search(x[8])
            if not m: return key, f"{gff}: CDS without prodigal ID=<k>_<n>: {l[:120]}", None
            c, n, s, e = x[0], int(m.group(1)), int(x[3]), int(x[4])
            if c in host:
                if len(first) < check_n: first.append((c, n, s, e, x[6]))
                if (c, n) in seen: return key, f"{gff}: gene {c}_{n} twice", None
                seen.add((c, n)); p = RE_PART.search(x[8])
                rows.append((c, n, s, e, 1 if x[6] == "+" else -1, bool(p and p.group(1) != "00"), mem.pop((c, n), -1)))
    if mem:
        c, n = next(iter(mem))
        return key, f"{len(mem)} catalogue member genes not in {gff} (e.g. {c}_{n}) - catalogue from another gene prediction?", None
    checked, pcheck = 0, []
    if first and fna:   # gene name = <contig>_<n> and the GFF coordinates must reproduce MATAFILER4's gene FASTA
        want = {f"{c}_{n}": (c, s, e, st) for c, n, s, e, st in first}
        got, name, todo = {}, None, set(want)
        with opener(fna) as fh:
            for l in fh:
                if l[0] == ">":
                    if not todo: break
                    name = l[1:].split()[0]
                    if name in todo: todo.discard(name); got[name] = []
                    else: name = None
                elif name: got[name].append(l.strip())
        if not got: return key, f"none of the first {len(want)} GFF genes (<contig>_<n>, e.g. {next(iter(want))}) is in {fna}", None
        with fasta_open(fasta, fai, gzi) as F:
            bad = sum("".join(got[g]).upper() != (F.fetch(c, s - 1, e).upper() if st == "+" else rc(F.fetch(c, s - 1, e).upper()))
                      for g, (c, s, e, st) in want.items() if g in got)
        if bad: return key, f"{bad} of {len(got)} genes differ between {fna} and their GFF coordinates on {fa}", None
        checked = len(got)
        info = {(r[0], r[1]): (r[5], r[6]) for r in rows}   # complete genes >= 300 nt: protein vs their catalogue gene's (in prep)
        pcheck = [(info[(c, n)][1], translate("".join(got[f"{c}_{n}"]).upper())) for c, n, *_ in first
                  if f"{c}_{n}" in got and not info[(c, n)][0] and info[(c, n)][1] >= 0 and len("".join(got[f"{c}_{n}"])) >= 300]
    rows.sort(key=lambda r: (host[r[0]][1], host[r[0]][0], r[0], r[2], r[3]))
    col = lambda i: [r[i] for r in rows]
    write_pq(dict(MGS=[host[r[0]][1] for r in rows], MAG=[host[r[0]][0] for r in rows], contig=col(0), n=col(1), start=col(2),
                  end=col(3), strand=col(4), partial=col(5), cat=col(6)), f"{P}/genes/{key}.parquet",
             dict(n=pa.int32(), start=pa.int64(), end=pa.int64(), strand=pa.int8(), cat=pa.int64()))
    write_pq(dict(contig=list(host), length=[lens[c] for c in host], MAG=[v[0] for v in host.values()], MGS=[v[1] for v in host.values()]),
             f"{P}/contigs/{key}.parquet")
    return key, None, dict(key=key, assembly=asm, fasta=rel(fasta, P), fai=rel(fai, P), gzi=rel(gzi, P) if gzi else "", genes=len(rows),
                           no_catalogue=sum(r[6] < 0 and r[3] - r[2] >= 99 for r in rows), checked=checked, pcheck=pcheck)


def cmd_prep(a):
    from concurrent.futures import ProcessPoolExecutor
    import pyarrow.csv as pcsv, pyarrow.parquet as pq
    P = a.out
    for d in ("members", "genes", "contigs", "fasta"): os.makedirs(f"{P}/{d}", exist_ok=True)
    host, asm_ctg = {}, collections.defaultdict(dict)
    with gzip.open(f"{a.idx}/contigs.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = (r["assembly"], r["MAG"], r["MGS"])
            if host.setdefault(r["contig"], v) != v: die(f"contig {r['contig']} in two MAGs - contig names not unique across assemblies")
            asm_ctg[r["assembly"]][r["contig"]] = (r["MAG"], r["MGS"])
    if not host: die("the MMT index lists no contig of a usable MGS MAG")
    key = {asm: "a" + md5(asm, 12) for asm in asm_ctg}
    cidx, prot = catalogue_files(a.genecat, a.gc_id)
    log(f"prep: {len(host)} MAG contigs in {len(asm_ctg)} assemblies ; catalogue {cidx}")
    # 1) one pass over the catalogue: members on MGS-MAG contigs (one file per assembly) ; genes found in >= 2 MGS
    import resource
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    if soft != resource.RLIM_INFINITY and soft < len(key) + 64:
        resource.setrlimit(resource.RLIMIT_NOFILE, (hard if hard == resource.RLIM_INFINITY else min(hard, len(key) + 1024),) * 2)
    outs = {k: open(f"{P}/members/{k}.tsv", "w") for k in key.values()}
    n_mem = n_cat = n_other = 0
    with opener(cidx) as fh, open(f"{P}/.cat_mgs.tsv", "w") as cm:
        cm.write("cat\tn_mgs\n")
        for line in fh:
            i, _, mem = line.rstrip("\n").partition("\t")
            if not i.isdigit(): continue
            ms = set(); n_cat += 1
            for m in mem.split(","):
                ctg, _, n = (m[1:] if m[:1] == ">" else m).rpartition("_")
                h = host.get(ctg)
                if h is None: n_other += 1; continue
                outs[key[h[0]]].write(f"{ctg}\t{n}\t{i}\n"); ms.add(h[2]); n_mem += 1
            if len(ms) > 1: cm.write(f"{i}\t{len(ms)}\n")
    for o in outs.values(): o.close()
    if n_mem == 0: die(f"no member of {cidx} lies on an MGS MAG contig - member names <contig>_<n> do not match the bin files?")
    log(f"prep: {n_cat} catalogue genes ; {n_mem} members on MGS MAG contigs, {n_other} elsewhere")
    # 2) per assembly, in parallel: genes + checks
    meta, errs = [], []
    with ProcessPoolExecutor(max_workers=max(1, a.threads)) as ex:
        for k, err, m in ex.map(prep_assembly, [(key[s], s, c, P, a.check, a.samtools) for s, c in sorted(asm_ctg.items())]):
            (errs.append(f"{k}: {err}") if err else meta.append(m))
    if errs: die(f"prep failed for {len(errs)} assemblies: " + " | ".join(errs[:5]))
    tot, noc = sum(m["genes"] for m in meta), sum(m["no_catalogue"] for m in meta)
    if noc > 0.02 * tot: die(f"{noc} of {tot} genes (>= 100 nt) on MGS MAG contigs have no catalogue entry - catalogue made without these assemblies?")
    pq.write_table(pcsv.read_csv(f"{P}/.cat_mgs.tsv", parse_options=pcsv.ParseOptions(delimiter="\t")), f"{P}/cat_mgs.parquet", compression="zstd")
    os.remove(f"{P}/.cat_mgs.tsv")
    # catalogue proteins: their .fai as a table (each MGS job indexes only its own catalogue genes)
    import pyarrow as pa
    pfa, pfai, pgzi = fasta_access(prot, f"{P}/fasta/catalogue_proteins", a.samtools)
    pq.write_table(pcsv.read_csv(pfai, read_options=pcsv.ReadOptions(column_names=["name", "length", "offset", "line_bases", "line_width"]),
                                 parse_options=pcsv.ParseOptions(delimiter="\t"), convert_options=pcsv.ConvertOptions(column_types={"name": pa.string()})),
                   f"{P}/catalogue_proteins_fai.parquet", compression="zstd")
    cat_f = dict(proteins=rel(pfa, P), proteins_gzi=rel(pgzi, P) if pgzi else "", proteins_index="catalogue_proteins_fai.parquet")
    # the proteins of genes checked against MATAFILER4's gene FASTA must resemble their catalogue gene's protein: a
    # catalogue made from another gene prediction of the same contigs passes every name-based check
    want = {c for m in meta for c, _ in m["pcheck"]}
    if want:
        with protein_file(P, want, cat_f, f"{P}/fasta/check") as F:
            names, bad = set(F.references), []
            for m in meta:
                t = [(c, q) for c, q in m["pcheck"] if str(c) in names]
                unlike = sum(len(kmers(q) & kmers(F.fetch(str(c)).upper().rstrip("*"))) < 0.2 * max(1, len(kmers(q))) for c, q in t)
                if len(m["pcheck"]) >= 10 and len(t) < len(m["pcheck"]) / 2:
                    bad.append(f"{m['assembly']}: {len(m['pcheck']) - len(t)} catalogue genes have no protein in {prot} (its headers must be the catalogue ids, a .fai next to it complete)")
                elif len(t) >= 10 and unlike > 0.05 * len(t):
                    bad.append(f"{m['assembly']}: {unlike} of {len(t)} checked genes do not resemble the protein of their catalogue gene "
                               "(catalogue made from another gene prediction - genePred redone after GeneCat?)")
        if bad: die(f"catalogue check failed for {len(bad)} assemblies: " + " | ".join(bad[:3]))
    smp = collections.defaultdict(list)
    for r in tsv(f"{a.idx}/samples.tsv"): smp[r["assembly"]].append(r["sample"])
    write_text(f"{P}/assemblies.tsv", "key\tassembly\tsamples\tfasta\tfai\tgzi\tgenes\tno_catalogue\tchecked\tprotein_checked\n" + "".join(
        f"{m['key']}\t{m['assembly']}\t{','.join(smp.get(m['assembly'], []))}\t{m['fasta']}\t{m['fai']}\t{m['gzi']}\t{m['genes']}\t"
        f"{m['no_catalogue']}\t{m['checked']}\t{len(m['pcheck'])}\n" for m in sorted(meta, key=lambda m: m["key"])))
    write_text(f"{P}/catalogue.tsv", f"cluster_index\t{cidx}\n" + "".join(f"{k}\t{v}\n" for k, v in cat_f.items()))
    shutil.rmtree(f"{P}/members")
    log(f"prep: done - {tot} genes ({noc} of >= 100 nt without catalogue entry), {sum(m['checked'] for m in meta)} checked against "
        f"genes.shrtHD.fna, {sum(len(m['pcheck']) for m in meta)} against their catalogue protein")


# ============================================================================================== plan
def cmd_plan(a):
    fp, nm = collections.defaultdict(list), collections.Counter()     # base_fp covers prep (genes, catalogue) + settings
    for r in tsv(f"{a.idx}/mags.tsv"):
        if r["status"] == "ok": nm[r["MGS"]] += 1; fp[r["MGS"]].append("|".join(r[k] for k in ("MAG", "is_rep", "assembly", "completeness", "contamination")))
    known = [r["MGS"] for r in tsv(f"{a.idx}/mgs.tsv")]
    sel = [m for m in known if m in fp]
    if a.mgs:
        want = [l.split()[0] for l in open(a.mgs) if l.strip() and not l.startswith("#")]
        unk = [m for m in want if m not in set(known)]
        if unk: die(f"--mgs: not in MAGvsGC: {' '.join(unk[:10])}")
        if any(m not in fp for m in want): log(f"WARNING: no usable MAG, skipped: {' '.join(m for m in want if m not in fp)}")
        sel = [m for m in dict.fromkeys(want) if m in fp]
    if a.part:
        i, n = (int(x) for x in a.part.split("/"))
        if not 1 <= i <= n: die(f"--part {a.part}: needs I/N with 1 <= I <= N")
        load, part = [0] * n, {}                    # largest first, each to the least loaded part (cost ~ MAGs + 1)
        for m in sorted(sel, key=lambda m: (-nm[m], m)):
            k = min(range(n), key=lambda j: (load[j], j)); part[m] = k; load[k] += nm[m] + 1
        sel = [m for m in sel if part[m] == i - 1]
    if not (a.mgs or a.part):          # result folders (they hold a .fp) of MGS that are gone
        for d in sorted(os.listdir(a.out)) if os.path.isdir(a.out) else []:
            if not d.startswith((".", "_")) and d not in fp and os.path.isfile(f"{a.out}/{d}/.fp"): print(f"{d}\tstale\t-\t0")
    for m in sel:
        f = md5(f"{a.base_fp}|{m}|" + ",".join(sorted(fp[m])))
        d = f"{a.out}/{m}"
        cur = not a.force and os.path.isfile(f"{d}/summary.tsv") and os.path.isfile(f"{d}/.fp") and open(f"{d}/.fp").read().strip() == f
        print(f"{m}\t{'current' if cur else 'todo'}\t{f}\t{nm[m]}")


# ============================================================================================== mgs
def side_r(st): return 1 if st > 0 else 0          # gene end facing the NEXT gene along the contig (1 = 3' end)
def side_l(st): return 0 if st > 0 else 1          # gene end facing the PREVIOUS gene


class MGS:
    """one MGS: every gene of every member MAG, in contig order"""
    def __init__(s, a):
        import pyarrow as pa, pyarrow.parquet as pq
        s.a, s.m, s.S = a, a.mgs, json.loads(a.settings)
        at = lambda x: os.path.join(a.prep, x)                       # relative = inside the prep folder
        s.asm = {r["assembly"]: dict(r, fasta=at(r["fasta"]), fai=at(r["fai"]), gzi=r["gzi"] and at(r["gzi"])) for r in tsv(f"{a.prep}/assemblies.tsv")}
        s.MAG = {r["MAG"]: r for r in tsv(f"{a.idx}/mags.tsv") if r["MGS"] == s.m and r["status"] == "ok"}
        if not s.MAG: die(f"{s.m}: no usable MAG in the index")
        s.mags = sorted(s.MAG); mi = {g: i for i, g in enumerate(s.mags)}
        keys = sorted({s.asm[s.MAG[g]["assembly"]]["key"] for g in s.mags})
        rd = lambda d: pa.concat_tables([pq.read_table(f"{a.prep}/{d}/{k}.parquet", filters=[("MGS", "=", s.m)]) for k in keys]).to_pydict()
        T, C = rd("genes"), rd("contigs")
        s.clen, s.cmag = dict(zip(C["contig"], C["length"])), {c: mi[g] for c, g in zip(C["contig"], C["MAG"])}
        cats = sorted({c for c in T["cat"] if c >= 0})
        s.shared = dict(zip(*pq.read_table(f"{a.prep}/cat_mgs.parquet", filters=[("cat", "in", cats)]).to_pydict().values())) if cats else {}
        s.cat_f = {k: v and at(v) for k, v in (l.rstrip("\n").split("\t") for l in open(f"{a.prep}/catalogue.tsv"))}
        s.N = len(T["contig"])
        s.gm, s.ctg, s.n = [mi[x] for x in T["MAG"]], T["contig"], T["n"]
        s.st, s.en, s.sd, s.cat = T["start"], T["end"], T["strand"], T["cat"]
        s.masked = [p or c < 0 for p, c in zip(T["partial"], T["cat"])]
        s.partial = T["partial"]
        s.contigs = []                                  # (mag index, contig, [gene indices along the contig])
        for g in sorted(range(s.N), key=lambda g: (s.gm[g], s.ctg[g], s.st[g], s.en[g])):
            if not s.contigs or s.contigs[-1][1] != s.ctg[g]: s.contigs.append((s.gm[g], s.ctg[g], []))
            s.contigs[-1][2].append(g)
        s.comp = [min(1.0, max(0.05, num(s.MAG[g]["completeness"], 100) / 100)) for g in s.mags]
        s.qual = [num(s.MAG[g]["completeness"], 0) - 5 * num(s.MAG[g]["contamination"], 0) for g in s.mags]
        s.smp = [s.asm[s.MAG[g]["assembly"]]["samples"] for g in s.mags]
        s.nm = len(s.mags)
        s.edge_min, s.fam_min = min(s.S["min_edge_mags"], s.nm), min(s.S["min_family_mags"], s.nm)
        # mixed MAG = co-binned strains: > mixed_mag of its catalogue genes in several copies. Its copies are strains, not
        # paralogs: they neither make linkers nor stop the merging of alleles that co-occur in it
        cc = [collections.Counter() for _ in range(s.nm)]
        for g in range(s.N):
            if not s.masked[g]: cc[s.gm[g]][s.cat[g]] += 1
        s.mixed = {i for i, c in enumerate(cc) if c and sum(k > 1 for k in c.values()) > s.S["mixed_mag"] * len(c)}

    # ------------------------------------------------------------------ 1. families
    def families(s):
        """catalogue genes merged into one family when they fill the same slot (identical flanking gene ends) in
        different MAGs, never co-occur in a MAG, and their proteins are >= fam_min_id identical over the shorter one
        (length ratio >= fam_len_ratio) - Panaroo collapse_families / MG-STK locus groups"""
        S, parent = s.S, {}
        def find(x):
            while parent.get(x, x) != x: parent[x] = parent.get(parent[x], parent[x]); x = parent[x]
            return x
        cm = collections.defaultdict(int)                          # catalogue gene -> MAGs (bitset) with a complete copy
        for g in range(s.N):
            if not s.masked[g]: cm[s.cat[g]] |= 0 if s.gm[g] in s.mixed else 1 << s.gm[g]
        fmags, frep, s.merges = dict(cm), {c: c for c in cm}, []
        if S["fam_min_id"] > 0:
            import edlib
            prot, cache = [], {}
            def pseq(c):
                if not prot:                                       # only an MGS with a merge candidate reads proteins
                    os.makedirs(s.a.tmp, exist_ok=True)
                    prot.append(protein_file(s.a.prep, {x for x in s.cat if x >= 0}, s.cat_f, f"{s.a.tmp}/derep_{s.m}_{os.getpid()}"))
                try: return prot[0].fetch(str(c)).rstrip("*").upper()
                except KeyError: return ""
            def similar(x, y):
                k = (min(x, y), max(x, y))
                if k not in cache:
                    p, q = sorted((pseq(x), pseq(y)), key=len); cache[k] = 0.0
                    if p and len(p) >= S["fam_len_ratio"] * len(q):
                        d = edlib.align(p, q, mode="HW", task="distance", k=int(len(p) * (1 - S["fam_min_id"])) + 1)["editDistance"]
                        cache[k] = 1 - d / len(p) if d >= 0 else 0.0
                return cache[k]
            for _ in range(3):                                     # merges can create new identical slots
                slots = collections.defaultdict(set)
                for mi, c, gl in s.contigs:
                    t = [(find(s.cat[g]), s.sd[g]) for g in gl if not s.masked[g]]
                    for i in range(1, len(t) - 1):
                        L, R, o = (t[i - 1][0], side_r(t[i - 1][1])), (t[i + 1][0], side_l(t[i + 1][1])), t[i][1]
                        if (R, L) < (L, R): L, R, o = R, L, -o
                        slots[(L, R, o)].add(t[i][0])
                cand = sorted({(x, y) for v in slots.values() for x in v for y in v if x < y})
                n_new = 0
                for x, y in cand:
                    rx, ry = find(x), find(y)
                    if rx == ry or fmags[rx] & fmags[ry]: continue
                    ident = similar(frep[rx], frep[ry])
                    if ident < S["fam_min_id"]: continue
                    s.merges.append((frep[rx], frep[ry], round(ident, 4))); n_new += 1
                    r, o = min(rx, ry), max(rx, ry)
                    parent[o] = r; fmags[r] |= fmags.pop(o)
                    frep[r] = max(frep[r], frep.pop(o), key=lambda c: (cm[c].bit_count(), -c))
                if not n_new: break
            for x in prot: x.close()
        s.fam = [find(c) if c >= 0 else -1 for c in s.cat]
        s.fam_cats = collections.defaultdict(set)
        for g in range(s.N):
            if s.cat[g] >= 0: s.fam_cats[s.fam[g]].add(s.cat[g])
        s.fam_shared = {f: max((s.shared.get(c, 1) for c in cs), default=1) for f, cs in s.fam_cats.items()}

    # ------------------------------------------------------------------ 2. family classes
    def classes(s):
        """linker = multi-copy (>= max_copy copies in a MAG, or > 1 copy in >= dup_margin of its MAGs and in >= 2 MAGs;
        mixed MAGs not counted: PPanGGOLiN multigenic) or hub (>= hub_partners different neighbours at BOTH ends, none in half of them: insertion
        elements) ; rare = in < min_family_mags MAGs ; the rest are anchors. Linkers and rare genes are bridged."""
        S, fam = s.S, s.fam
        s.nmag, s.multi, s.maxc, s.pres = (collections.Counter() for _ in range(4))
        part, endm = collections.defaultdict(collections.Counter), collections.Counter()
        for mi, cl in itertools.groupby(s.contigs, key=lambda x: x[0]):
            cl = list(cl); c = collections.Counter(fam[g] for _, _, gl in cl for g in gl if not s.masked[g])
            for f, k in c.items():
                s.nmag[f] += 1
                if mi not in s.mixed: s.multi[f] += k > 1; s.maxc[f] = max(s.maxc[f], k)
            s.pres.update({fam[g] for _, _, gl in cl for g in gl if fam[g] >= 0})   # partial copies at contig ends count
        s.linker = {f for f in s.nmag if s.maxc[f] >= S["max_copy"] or s.multi[f] >= max(2, S["dup_margin"] * s.nmag[f])}
        for mi, cl in itertools.groupby(s.contigs, key=lambda x: x[0]):
            pairs = set()
            for _, _, gl in cl:
                t = [g for g in gl if not s.masked[g] and fam[g] not in s.linker]
                for g, h in zip(t, t[1:]):
                    pairs.add(((fam[g], side_r(s.sd[g])), fam[h])); pairs.add(((fam[h], side_l(s.sd[h])), fam[g]))
            for e, f in pairs: part[e][f] += 1
            endm.update({e for e, _ in pairs})
        dispersed = lambda e: len(part[e]) >= S["hub_partners"] and max(part[e].values()) < 0.5 * endm[e]
        s.hub = {f for f in s.nmag if f not in s.linker and S["hub_partners"] > 0 and dispersed((f, 0)) and dispersed((f, 1))}
        s.tot_comp = sum(s.comp)
        s.padj = {f: s.pres[f] / s.tot_comp for f in s.pres}
        s.anchor = {f for f in s.nmag if f not in s.linker and f not in s.hub and s.nmag[f] >= s.fam_min}
        s.core_anchor = {f for f in s.anchor if s.padj[f] >= S["core_prev"]}
        L = collections.defaultdict(list)
        for g in range(s.N):
            if not s.masked[g]: L[fam[g]].append(s.en[g] - s.st[g] + 1)
        s.mlen = {f: statistics.median_low(v) for f, v in L.items()}

    # ------------------------------------------------------------------ 3. blocks
    def build(s, ok, sep):
        """unitigs of the adjacency graph of anchors `ok` (other genes bridged; genes of `sep` cut the path). Two gene
        ends x, y are joined (solid) when each is the other's most frequent neighbour (ties: no join), the pair is seen
        in >= min_edge_mags MAGs, at BOTH ends the MAGs agreeing with the order make up >= synteny of the MAGs with any
        neighbour there (contig ends do not count), and an insert of more than max_insert_bp beyond the usual spacing
        occurs in at most (1 - synteny) of them. Agreeing = the same neighbour, or another one that only skips (a
        deletion of) or adds (an insertion of) <= max_skip genes: it faces y's far end, or y, within max_skip anchors in
        some MAG - so presence/absence of a few genes never splits a block, rearrangements do."""
        S, fam, t = s.S, s.fam, max(1, s.S["max_skip"])                                  # t = 1 for the k == 1 pairs (E)
        T, E = collections.defaultdict(int), collections.defaultdict(lambda: [0, []])     # MAG bitsets ; [bitset, gaps]
        nb, obs = collections.defaultdict(set), collections.defaultdict(int)
        key = lambda x, y: (x, y) if x < y else (y, x)
        def walk(mi, path):
            for i, g in enumerate(path):
                x = (fam[g], side_r(s.sd[g]))
                for k in range(1, t + 1):
                    if i + k >= len(path): break
                    h = path[i + k]
                    if fam[h] == fam[g]: continue
                    y, b = (fam[h], side_l(s.sd[h])), 1 << mi
                    if k <= S["max_skip"]: T[key(x, y)] |= b
                    if k == 1:
                        d = E[key(x, y)]; d[0] |= b; d[1].append(s.st[h] - s.en[g] - 1)
                        nb[x].add(y); nb[y].add(x); obs[x] |= b; obs[y] |= b
        for mi, c, gl in s.contigs:
            path = []
            for g in gl:
                if s.masked[g]: continue
                if fam[g] in sep: walk(mi, path); path = []
                elif fam[g] in ok: path.append(g)
            walk(mi, path)
        far = lambda e: (e[0], 1 - e[1])
        sup = lambda x, y: E[key(x, y)][0]
        def best(x):
            r = sorted(nb[x], key=lambda y: -sup(x, y).bit_count())
            return r[0] if r and (len(r) == 1 or sup(x, r[0]).bit_count() > sup(x, r[1]).bit_count()) else None
        seen_t = lambda x, y: T.get(key(x, y), 0).bit_count() >= s.edge_min
        def share(x, y):
            agree = sup(x, y)
            for z in nb[x]:
                if z != y and (seen_t(far(y), z) or seen_t(far(z), y)): agree |= sup(x, z)
            return agree.bit_count() / obs[x].bit_count()
        solid = {}
        for (x, y), (b, gaps) in list(E.items()):
            n = b.bit_count()
            if n < s.edge_min or best(x) != y or best(y) != x: continue
            if share(x, y) < S["synteny"] or share(y, x) < S["synteny"]: continue
            gaps = sorted(gaps); med = gaps[(len(gaps) - 1) // 2]
            if sum(v - med > S["max_insert_bp"] for v in gaps) > (1 - S["synteny"]) * len(gaps): continue
            solid[x], solid[y] = (y, med, n), (x, med, n)
        blocks, seen = [], set()
        for f in sorted(ok & set(s.nmag)):
            if f in seen: continue
            left, cur, circ = [], (f, 0), False
            while cur in solid:                         # walk left: g entered at end sg from its right side
                (g, sg), med, n = solid[cur]
                if g == f: circ = True; break
                left.append((g, 1 if sg == 1 else -1, med)); cur = (g, 1 - sg)
            right, cur, close = [], (f, 1), None
            while cur in solid:                         # walk right: g entered at end sg from its left side
                (g, sg), med, n = solid[cur]
                if g == f: close = (med, n); break
                right.append((g, 1 if sg == 0 else -1, med, n)); cur = (g, 1 - sg)
            if circ:                                    # a cycle: cut at its weakest junction (fewest MAGs)
                ring = [(f, 1)] + [(g, o) for g, o, _, _ in right]
                link = [(m_, n_) for _, _, m_, n_ in right] + [close]      # link[i] joins ring[i] -> ring[i+1]
                c = min(range(len(link)), key=lambda i: (link[i][1], i)) + 1
                path, gaps = ring[c:] + ring[:c], [m_ for m_, _ in link[c:] + link[:c]][:-1]
            else:
                path = [(g, o) for g, o, _ in reversed(left)] + [(f, 1)] + [(g, o) for g, o, _, _ in right]
                gaps = [x[2] for x in reversed(left)] + [x[2] for x in right]
            seen.update(g for g, _ in path)
            if path[0][0] > path[-1][0]: path, gaps = [(g, -o) for g, o in reversed(path)], gaps[::-1]
            bp = sum(s.mlen[g] for g, _ in path) + sum(max(0, x) for x in gaps)
            if len(path) >= S["min_block_genes"] and bp >= S["min_block_bp"]: blocks.append((path, gaps, circ))
        return blocks

    def blocks(s):
        """tier 1 (core): anchors with completeness-adjusted prevalence >= core_prev - the persistent backbone, all
        variable content bridged ; tier 2 (accessory): the other anchors, between core-block genes"""
        bl = [(p, g, c, "core") for p, g, c in s.build(s.core_anchor, set())]
        in_core = {f for p, _, _, _ in bl for f, _ in p}
        bl += [(p, g, c, "accessory") for p, g, c in s.build(s.anchor - in_core, in_core)]
        bl.sort(key=lambda b: (b[3] != "core", -len(b[0]), b[0][0][0]))
        s.B, s.bid = bl, [f"B{i + 1:04d}" for i in range(len(bl))]
        s.bpos, s.cs, s.ce = {}, [], []                 # family -> (block, pos, orientation) ; consensus gene coords
        for bi, (path, gaps, _, _) in enumerate(bl):
            x, a_, b_ = 0, [], []
            for j, (f, o) in enumerate(path):
                s.bpos[f] = (bi, j, o); a_.append(x); x += s.mlen[f]; b_.append(x)
                if j < len(gaps): x += max(0, gaps[j])
            s.cs.append(a_); s.ce.append(b_)

    # ------------------------------------------------------------------ 4. instances
    def instances(s):
        """a block on a contig: a run of its genes in one direction, <= max_skip of its genes missing in a row and
        <= max_insert_bp more sequence than the consensus spacing between two of them. Each end of an instance is
        complete, contig_end (the contig ends before the missing part could fit: incomplete, not absent) or open."""
        S, I = s.S, []
        for mi, c, gl in s.contigs:
            occ = collections.defaultdict(list)
            for g in gl:
                x = None if s.masked[g] else s.bpos.get(s.fam[g])
                if x: occ[x[0]].append((g, x[1], s.sd[g] * x[2]))
            for bi, L in occ.items():
                cs, ce, nb = s.cs[bi], s.ce[bi], len(s.B[bi][0])
                def emit(run, d, dup):
                    strand, ps = d or run[0][2], [p for _, p, _ in run]
                    lo, hi = min(ps), max(ps)
                    s0, e0 = min(s.st[g] for g, _, _ in run), max(s.en[g] for g, _, _ in run)
                    ins = sum(max(0, (s.st[h] - s.en[g] - 1) - (cs[max(p, q)] - ce[min(p, q)]))
                              for (g, p, _), (h, q, _) in zip(run, run[1:]))
                    # the block goes on along this contig beyond an end (L is in contig order): its missing genes are not cut off
                    room_l = float("inf") if s.st[L[0][0]] < s0 else s0 - 1
                    room_r = float("inf") if s.en[L[-1][0]] > e0 else s.clen[c] - e0
                    room_lo, room_hi = (room_l, room_r) if strand > 0 else (room_r, room_l)
                    stat = lambda miss, room, need: "complete" if miss == 0 else ("contig_end" if room < need else "open")
                    I.append(dict(block=bi, mag=mi, contig=c, start=s0, end=e0, strand=strand, first=lo, last=hi,
                                  present=len(set(ps)), flipped=sum(r != strand for _, _, r in run), dup=dup, insert_bp=ins,
                                  start_status=stat(lo, room_lo, cs[lo]), end_status=stat(nb - 1 - hi, room_hi, ce[-1] - ce[hi]),
                                  genes=[(g, p) for g, p, _ in run]))
                run, d, dup = [L[0]], 0, 0
                for e in L[1:]:
                    p, dp = run[-1][1], e[1] - run[-1][1]
                    near = (s.st[e[0]] - s.en[run[-1][0]] - 1) - (cs[max(p, e[1])] - ce[min(p, e[1])]) <= S["max_insert_bp"]
                    if dp == 0 and near: dup += 1; continue  # a second copy inside the run: not part of it
                    # a jump over missing positions must read in the run's direction (else: first gene of an inverted piece)
                    ok = dp != 0 and near and (d == 0 or (dp > 0) == (d > 0)) and abs(dp) <= S["max_skip"] + 1 and \
                        (abs(dp) == 1 or e[2] == (1 if dp > 0 else -1))
                    if ok: d = d or (1 if dp > 0 else -1); run.append(e)
                    else: emit(run, d, dup); run, d, dup = [e], 0, 0
                emit(run, d, dup)
        I.sort(key=lambda r: (r["block"], r["mag"], r["contig"], r["start"]))
        for i, r in enumerate(I): r["id"], r["host"], r["host_pos"] = i, -1, -1
        s.I = I
        s.by_block = collections.defaultdict(list)
        for r in I: s.by_block[r["block"]].append(r)
        # every gene: its instance (block genes: own run ; other genes: the innermost instance whose span holds them)
        s.g_inst, s.g_slot, s.g_pos = {}, {}, {}
        for r in I:
            for g, p in r["genes"]: s.g_inst[g], s.g_pos[g] = r["id"], p
        by_ctg = collections.defaultdict(list)
        for r in I: by_ctg[r["contig"]].append(r)
        s.by_ctg, s.links = by_ctg, collections.defaultdict(lambda: [set(), [], 0])
        cg = {c: gl for _, c, gl in s.contigs}
        cst = {c: [s.st[g] for g in gl] for c, gl in cg.items()}        # gene starts along each contig (sorted)
        for c, rs in by_ctg.items():
            rs.sort(key=lambda r: (r["start"], -r["end"]))
            top = []
            for r in rs:                                # nesting: an instance inside another one's span = insert
                h = next((h for h in reversed(top) if h["start"] <= r["start"] and r["end"] <= h["end"]), None)
                if h is None: top.append(r); continue
                before = [p for g, p in h["genes"] if (s.en[g] <= r["start"] if h["strand"] > 0 else s.st[g] >= r["end"])]
                r["host"], r["host_pos"] = h["id"], max(before, default=-1)
            for r in sorted(rs, key=lambda r: r["end"] - r["start"]):
                gp = {g: p for g, p in r["genes"]}
                last, inside = -1, cg[c][bisect.bisect_left(cst[c], r["start"]):bisect.bisect_right(cst[c], r["end"])]
                for g in (inside if r["strand"] > 0 else reversed(inside)):
                    if not (r["start"] <= s.st[g] and s.en[g] <= r["end"]): continue
                    if g in gp: last = gp[g]
                    elif g not in s.g_inst: s.g_inst[g], s.g_slot[g] = r["id"], last
            for x, y in zip(top, top[1:]):              # links: block ends meeting along the contig
                ea = (x["block"], x["last"], 3) if x["strand"] > 0 else (x["block"], x["first"], 5)
                eb = (y["block"], y["first"], 5) if y["strand"] > 0 else (y["block"], y["last"], 3)
                k = (ea, eb) if ea <= eb else (eb, ea)
                term = all((p == 0 and sd == 5) or (p == len(s.B[b][0]) - 1 and sd == 3) for b, p, sd in k)
                L = s.links[k]; L[0].add(x["mag"]); L[1].append(y["start"] - x["end"] - 1); L[2] = term
            for k_, x in enumerate(top):                # neighbours along the contig (top-level instances)
                x["nbL"], x["nbR"] = (top[k_ - 1] if k_ else None), (top[k_ + 1] if k_ + 1 < len(top) else None)

    # ------------------------------------------------------------------ 5. representatives + sequences
    def representatives(s):
        """per block the instance that is complete (all genes, both ends), then fewest flipped genes and inserts, then
        best MAG quality (completeness - 5 x contamination), then the MGS representative MAG; its sequence on the
        block's orientation is the block's sequence"""
        fa = {}
        def seq_of(r):
            k = s.asm[s.MAG[s.mags[r["mag"]]]["assembly"]]
            if k["key"] not in fa: fa[k["key"]] = fasta_open(k["fasta"], k["fai"], k["gzi"])
            q = fa[k["key"]].fetch(r["contig"], r["start"] - 1, r["end"]).upper()
            return q if r["strand"] > 0 else rc(q)
        s.seq_of, s.fa = seq_of, fa
        s.rep = {}
        for bi, rs in s.by_block.items():
            nb = len(s.B[bi][0])
            s.rep[bi] = max(rs, key=lambda r: (r["present"] == nb and r["start_status"] == r["end_status"] == "complete", r["present"],
                                               -r["flipped"], -r["insert_bp"], round(s.qual[r["mag"]], 3),
                                               s.MAG[s.mags[r["mag"]]]["is_rep"] == "1", s.mags[r["mag"]], -r["start"]))
        s.seqs = {bi: seq_of(r) for bi, r in s.rep.items()}

    # ------------------------------------------------------------------ 6. base-level lift, polished block sequences
    def lift(s):
        """every instance aligned to its block's representative sequence (minimap2 asm20 --eqx ; per instance the hits
        by score, each column used once on both axes). The block sequence is then POLISHED: a substitution wherever more
        than half of the aligned instances (representative included) carry another base - so no one MAG's private
        differences shape every comparison. lift = runs of identical bases between each instance and its polished block
        (exact-match segments, 1-based): block_pos = b_start + (pos - c_start) (ori +1) or b_start + (c_end - pos) (-1).
        Segments of different instances never overlap on a contig (the smaller instance keeps the bases: nested inserts,
        overlapping genes at block ends). identity = identical / aligned bases ; aligned_frac = aligned / instance bases."""
        import numpy as np
        s.ident, s.seg, s.npol, s.pseq = {}, [], {}, dict(s.seqs)
        if not s.B: return
        hits = collections.defaultdict(list)
        if s.S["identity"] and s.a.minimap2:
            tmp = f"{s.a.tmp}/derep_{s.m}_{os.getpid()}"; os.makedirs(tmp, exist_ok=True)
            try: hits = s._align(tmp)
            finally:
                if not s.a.keep_tmp: shutil.rmtree(tmp, ignore_errors=True)
        CODE_ = np.full(256, 4, np.uint8)
        for i, ch in enumerate(b"ACGT"): CODE_[ch] = i
        COMPC, BASES = np.array([3, 2, 1, 0, 4], np.uint8), np.frombuffer(b"ACGTN", np.uint8)
        enc = lambda q: CODE_[np.frombuffer(q.encode(), np.uint8)]

        def chunks(r, L):
            """aligned columns of instance r on its block: [(qa, ta, sgn)] ; qa = 0-based on the instance's sequence in
            block orientation, ta = 0-based on the block, sgn = +1 / -1 (hit on the reverse strand: qa runs backwards)"""
            if r is s.rep[r["block"]]:
                return [(np.arange(L, dtype=np.int64), np.arange(L, dtype=np.int64), 1)]
            out, qlen = [], r["end"] - r["start"] + 1
            qc, tc = np.zeros(qlen, bool), np.zeros(L, bool)
            for _, qs, qe, sg, ts, cg in sorted(hits.get(r["id"], []), key=lambda h: -h[0]):
                qa_, ta_, t, q = [], [], ts, (qs if sg > 0 else qe - 1)
                run = 0
                for n, op in RE_CIG.findall(cg + "0I"):     # sentinel flushes the last aligned run
                    n = int(n)
                    if op in "M=X": run += n; continue
                    if run:
                        ta_.append(np.arange(t, t + run)); qa_.append(q + sg * np.arange(run)); t += run; q += sg * run; run = 0
                    if op == "D": t += n
                    elif op == "I": q += sg * n
                if not ta_: continue
                qa, ta = np.concatenate(qa_), np.concatenate(ta_)
                k = ~qc[qa] & ~tc[ta]
                qa, ta = qa[k], ta[k]; qc[qa] = True; tc[ta] = True
                if len(qa): out.append((qa, ta, sg))
            return out

        for bi in range(len(s.B)):
            rs, R = s.by_block[bi], s.seqs[bi]
            L, Rc = len(R), enc(R)
            cnt, buf, nbuf = np.zeros(5 * L, np.int64), [], 0  # pass 1: votes per block column (base * L + column)
            for r in rs:
                qcode = enc(s.seq_of(r))
                for qa, ta, sg in chunks(r, L):
                    b = qcode[qa] if sg > 0 else COMPC[qcode[qa]]
                    buf.append(b.astype(np.int64) * L + ta); nbuf += len(ta)
                    if nbuf > 5e7: cnt += np.bincount(np.concatenate(buf), minlength=5 * L); buf, nbuf = [], 0
            if buf: cnt += np.bincount(np.concatenate(buf), minlength=5 * L)
            cnt = cnt.reshape(5, L)
            tot, best = cnt[:4].sum(0), cnt[:4].argmax(0)
            sub = (cnt[:4].max(0) * 2 > tot) & (best != Rc)
            P = Rc.copy(); P[sub] = best[sub]
            s.npol[bi] = int(sub.sum()); s.pseq[bi] = BASES[P].tobytes().decode()
            for r in rs:                                       # pass 2: exact-match runs against the polished block
                qcode, n_al, n_eq = enc(s.seq_of(r)), 0, 0
                for qa, ta, sg in chunks(r, L):
                    b = qcode[qa] if sg > 0 else COMPC[qcode[qa]]
                    eq = (b == P[ta]) & (b < 4)
                    n_al += len(qa); n_eq += int(eq.sum())
                    cont = (np.diff(ta) == 1) & (np.diff(qa) == sg)
                    st = np.flatnonzero(eq & ~np.r_[False, eq[:-1] & cont])
                    en = np.flatnonzero(eq & ~np.r_[eq[1:] & cont, False])
                    ori = sg * r["strand"]
                    for i, j in zip(st, en):
                        c0 = r["start"] + qa[i] if r["strand"] > 0 else r["end"] - qa[i]
                        c1 = r["start"] + qa[j] if r["strand"] > 0 else r["end"] - qa[j]
                        s.seg.append((r["id"], int(min(c0, c1)), int(max(c0, c1)), int(ta[i]) + 1, ori))   # b_start = lowest block pos
                if n_al: s.ident[r["id"]] = (n_eq / n_al, n_al / (r["end"] - r["start"] + 1))
        s._resolve()

    def _align(s, tmp):
        """minimap2 asm20 --eqx of every non-representative instance (block orientation) against the representatives ->
        {instance: [(AS, q_start, q_end, +1/-1, t_start, cigar)]} on its own block only"""
        with open(f"{tmp}/ref.fa", "w") as o:
            for bi in range(len(s.B)): o.write(f">{bi}\n{s.seqs[bi]}\n")
        n, hits = 0, collections.defaultdict(list)
        with open(f"{tmp}/q.fa", "w") as o:
            for r in s.I:
                if r is not s.rep[r["block"]]: o.write(f">{r['id']}\n{s.seq_of(r)}\n"); n += 1
        if not n: return hits
        p = subprocess.run([s.a.minimap2, "-x", "asm20", "-c", "--eqx", "--secondary=no", "-t", str(s.a.threads), f"{tmp}/ref.fa",
                            f"{tmp}/q.fa"], capture_output=True, text=True)
        if p.returncode: die(f"{s.m}: minimap2 failed: {p.stderr.strip()[-300:]}")
        for l in p.stdout.splitlines():
            x = l.split("\t")
            if int(x[5]) != s.I[int(x[0])]["block"]: continue
            tag = {f[:2]: f[5:] for f in x[12:]}
            hits[int(x[0])].append((int(tag.get("AS", x[9])), int(x[2]), int(x[3]), 1 if x[4] == "+" else -1, int(x[7]), tag.get("cg", "")))
        return hits

    def _resolve(s):
        """segments of different instances never overlap on a contig: instances in order of span (smallest first) keep
        their bases, larger ones lose the overlap ; then every segment is checked against the contig sequence"""
        inst = {r["id"]: r for r in s.I}
        by = collections.defaultdict(list)
        for g in s.seg: by[g[0]].append(g)
        claimed, out = collections.defaultdict(list), []      # contig -> sorted disjoint [start, end]
        for i in sorted(by, key=lambda i: (inst[i]["end"] - inst[i]["start"], i)):
            c, mine = inst[i]["contig"], []
            cl = claimed[c]
            for _, a_, b_, bs, ori in by[i]:
                pieces, x = [], a_
                k = bisect.bisect_right(cl, [a_, float("inf")]) - 1
                k = max(k, 0)
                while k < len(cl) and cl[k][0] <= b_:
                    u, v = cl[k]
                    if v >= x and u <= b_:
                        if u > x: pieces.append((x, u - 1))
                        x = max(x, v + 1)
                    k += 1
                if x <= b_: pieces.append((x, b_))
                for u, v in pieces:
                    mine.append((i, c, u, v, bs + (u - a_) if ori > 0 else bs + (b_ - v), ori))
            out += mine
            for _, _, u, v, _, _ in mine:                         # claim (kept sorted and merged)
                k = bisect.bisect_left(cl, [u, v])
                cl.insert(k, [u, v])
            merged = []
            for u, v in cl:
                if merged and u <= merged[-1][1] + 1: merged[-1][1] = max(merged[-1][1], v)
                else: merged.append([u, v])
            claimed[c] = merged
        out.sort(key=lambda g: (inst[g[0]]["mag"], g[1], g[2]))
        bad = 0
        for i, c, u, v, bs, ori in out:                           # every base of every segment = its block base
            r = inst[i]
            k = s.asm[s.MAG[s.mags[r["mag"]]]["assembly"]]
            q = s.fa[k["key"]].fetch(c, u - 1, v).upper()
            if (q if ori > 0 else rc(q)) != s.pseq[r["block"]][bs - 1:bs - 1 + v - u + 1]: bad += 1
        if bad: die(f"{s.m}: {bad} of {len(out)} lift segments do not match their block sequence (bug)")
        s.seg = out

    # ------------------------------------------------------------------ 7. outputs
    def write(s, out):
        import pyarrow as pa
        m, I, B, bid, nb = s.m, s.I, s.B, s.bid, range(len(s.B))
        tier = lambda bi: B[bi][3]
        with gzip.open(f"{out}/blocks.fna.gz", "wt", compresslevel=6) as o:
            for bi in nb:
                r = s.rep[bi]
                o.write(f">{bid[bi]} MGS={m} tier={tier(bi)} genes={len(B[bi][0])} len={len(s.pseq[bi])} "
                        f"rep={s.mags[r['mag']]}:{r['contig']}:{r['start']}-{r['end']}:{'+' if r['strand'] > 0 else '-'} "
                        f"polished={s.npol.get(bi, 0)}\n")
                q = s.pseq[bi]; o.write("".join(q[i:i + 80] + "\n" for i in range(0, len(q), 80)))
        F = lambda f: [f(r) for r in I]
        isrep = lambda r: r is s.rep[r["block"]]
        complete = lambda r: r["present"] == len(B[r["block"]][0]) and r["start_status"] == r["end_status"] == "complete"
        grp = collections.defaultdict(list)                   # copies of a block in one MAG ; dup = they share block genes
        for r in I: grp[(r["block"], r["mag"])].append(r)
        copies = {k: len(v) for k, v in grp.items()}
        dupin = {r["id"] for v in grp.values() if len(v) > 1 for r in v for y in v
                 if y is not r and y["first"] <= r["last"] and r["first"] <= y["last"]}

        def nbr(r, five):
            """neighbour at the block's 5' (five) or 3' end of instance r: (instance, key, gap bp, MAGs with that link)"""
            if r["host"] >= 0: return r["host"], f"host:{bid[s.I[r['host']]['block']]}", None, None
            left = (r["strand"] > 0) == five
            y = r.get("nbL") if left else r.get("nbR")
            if y is None: return None, "contig_end", (r["start"] - 1 if left else s.clen[r["contig"]] - r["end"]), None
            if left:
                ey = (y["block"], y["last"], 3) if y["strand"] > 0 else (y["block"], y["first"], 5)
                er = (r["block"], r["first"], 5) if r["strand"] > 0 else (r["block"], r["last"], 3)
                gap = r["start"] - y["end"] - 1
            else:
                ey = (y["block"], y["first"], 5) if y["strand"] > 0 else (y["block"], y["last"], 3)
                er = (r["block"], r["last"], 3) if r["strand"] > 0 else (r["block"], r["first"], 5)
                gap = y["start"] - r["end"] - 1
            k = (ey, er) if ey <= er else (er, ey)
            return y["id"], f"{bid[ey[0]]}:{ey[1]}:{ey[2]}", gap, (len(s.links[k][0]) if k in s.links else None)
        N5, N3 = {r["id"]: nbr(r, True) for r in I}, {r["id"]: nbr(r, False) for r in I}
        write_pq(dict(MGS=[m] * len(I), instance=F(lambda r: r["id"]), block=F(lambda r: bid[r["block"]]), tier=F(lambda r: tier(r["block"])),
                      MAG=F(lambda r: s.mags[r["mag"]]), samples=F(lambda r: s.smp[r["mag"]]), contig=F(lambda r: r["contig"]),
                      contig_length=F(lambda r: s.clen[r["contig"]]),
                      c_start=F(lambda r: r["start"]), c_end=F(lambda r: r["end"]), strand=F(lambda r: "+" if r["strand"] > 0 else "-"),
                      first_slot=F(lambda r: r["first"]), last_slot=F(lambda r: r["last"]), n_present=F(lambda r: r["present"]),
                      frac=F(lambda r: round(r["present"] / len(B[r["block"]][0]), 4)), complete=F(complete),
                      start_status=F(lambda r: r["start_status"]), end_status=F(lambda r: r["end_status"]),
                      n_flipped=F(lambda r: r["flipped"]), n_dup=F(lambda r: r["dup"]), insert_bp=F(lambda r: r["insert_bp"]),
                      host_instance=F(lambda r: r["host"]), host_slot=F(lambda r: r["host_pos"]),
                      nb5_instance=F(lambda r: N5[r["id"]][0]), nb5_key=F(lambda r: N5[r["id"]][1]), nb5_gap=F(lambda r: N5[r["id"]][2]),
                      nb5_mags=F(lambda r: N5[r["id"]][3]),
                      nb3_instance=F(lambda r: N3[r["id"]][0]), nb3_key=F(lambda r: N3[r["id"]][1]), nb3_gap=F(lambda r: N3[r["id"]][2]),
                      nb3_mags=F(lambda r: N3[r["id"]][3]),
                      n_mag_copies=F(lambda r: copies[(r["block"], r["mag"])]), dup_in_mag=F(lambda r: r["id"] in dupin),
                      identity=F(lambda r: round(s.ident[r["id"]][0], 5) if r["id"] in s.ident else None),
                      aligned_frac=F(lambda r: round(s.ident[r["id"]][1], 4) if r["id"] in s.ident else None), is_rep=F(isrep)),
                 f"{out}/instances.parquet")
        # lift: exact-match segments instance <-> polished block sequence
        G_ = s.seg
        write_pq(dict(MGS=[m] * len(G_), MAG=[s.mags[s.I[g[0]]["mag"]] for g in G_], instance=[g[0] for g in G_],
                      block=[bid[s.I[g[0]]["block"]] for g in G_], contig=[g[1] for g in G_], c_start=[g[2] for g in G_],
                      c_end=[g[3] for g in G_], b_start=[g[4] for g in G_], ori=[g[5] for g in G_]), f"{out}/lift.parquet")
        bs = []
        for bi in nb:
            rs = s.by_block[bi]; ms = {r["mag"] for r in rs}; ids = [s.ident[r["id"]][0] for r in rs if r["id"] in s.ident and r is not s.rep[bi]]
            r0, path, q = s.rep[bi], B[bi][0], s.pseq[bi]
            g_, c_ = q.count("G"), q.count("C")
            bs.append(dict(MGS=m, block=bid[bi], tier=tier(bi), n_genes=len(path), seq_genes=r0["present"], length=len(q),
                           consensus_length=s.ce[bi][-1], gc=round((g_ + c_) / max(1, len(q)), 4),
                           gc_skew=round((g_ - c_) / (g_ + c_), 4) if g_ + c_ else None, n_polished=s.npol.get(bi, 0), seq=q,
                           circular=B[bi][2], n_instances=len(rs), n_mags=len(ms), n_mags_complete=len({r["mag"] for r in rs if complete(r)}),
                           prevalence=round(len(ms) / s.nm, 4), prevalence_adj=round(len(ms) / s.tot_comp, 4),
                           n_samples=len({x for i in ms for x in s.smp[i].split(",") if x}),
                           median_identity=round(statistics.median(ids), 5) if ids else None,
                           shared_frac=round(sum(s.fam_shared.get(f, 1) > 1 for f, _ in path) / len(path), 4),
                           max_mgs=max(s.fam_shared.get(f, 1) for f, _ in path),
                           rep_MAG=s.mags[r0["mag"]], rep_contig=r0["contig"], rep_c_start=r0["start"], rep_c_end=r0["end"],
                           rep_strand="+" if r0["strand"] > 0 else "-", rep_complete=complete(r0),
                           rep_is_MGS_rep=s.MAG[s.mags[r0["mag"]]]["is_rep"] == "1"))
        cols = lambda rows, keys: {k: [r[k] for r in rows] for k in keys}
        BK = ["MGS", "block", "tier", "n_genes", "seq_genes", "length", "consensus_length", "gc", "gc_skew", "n_polished", "circular",
              "n_instances", "n_mags", "n_mags_complete", "prevalence", "prevalence_adj", "n_samples", "median_identity", "shared_frac",
              "max_mgs", "rep_MAG", "rep_contig", "rep_c_start", "rep_c_end", "rep_strand", "rep_complete", "rep_is_MGS_rep", "seq"]
        write_pq(cols(bs, BK), f"{out}/blocks.parquet")
        # block genes: consensus layout coordinates (cons_*) and block-sequence coordinates (seq_* = the frame of block_pos)
        npos = collections.Counter((r["block"], p) for r in I for p in {p for _, p in r["genes"]})
        bg = []
        for bi in nb:
            r0 = s.rep[bi]; rg = {p: g for g, p in r0["genes"]}
            for j, (f, o) in enumerate(B[bi][0]):
                g = rg.get(j)
                if g is None: rs_ = re_ = rgn = None
                elif r0["strand"] > 0: rs_, re_, rgn = s.st[g] - r0["start"] + 1, s.en[g] - r0["start"] + 1, f"{s.ctg[g]}_{s.n[g]}"
                else: rs_, re_, rgn = r0["end"] - s.en[g] + 1, r0["end"] - s.st[g] + 1, f"{s.ctg[g]}_{s.n[g]}"
                bg.append(dict(MGS=m, block=bid[bi], slot=j, family=f, orientation="+" if o > 0 else "-",
                               catalogue=sorted(s.fam_cats[f]), n_instances=npos[(bi, j)], median_len=s.mlen[f],
                               cons_start=s.cs[bi][j] + 1, cons_end=s.ce[bi][j], rep_gene=rgn, seq_start=rs_, seq_end=re_,
                               shared_mgs=s.fam_shared.get(f, 1)))
        GK = ["MGS", "block", "slot", "family", "orientation", "catalogue", "n_instances", "median_len", "cons_start", "cons_end",
              "rep_gene", "seq_start", "seq_end", "shared_mgs"]
        write_pq(cols(bg, GK), f"{out}/block_genes.parquet", dict(catalogue=pa.list_(pa.int64())))
        # every gene of every MAG: role, family, block position, instance -> fully traceable
        fam = s.fam
        def role(g):
            f = fam[g]
            if s.masked[g]: return "masked"
            if f in s.linker: return "linker"
            if f in s.hub: return "hub"
            if f not in s.anchor: return "rare"
            return tier(s.bpos[f][0]) if f in s.bpos else "unblocked"
        G = range(s.N)
        write_pq(dict(MGS=[m] * s.N, MAG=[s.mags[s.gm[g]] for g in G], contig=s.ctg, gene=[f"{c}_{n}" for c, n in zip(s.ctg, s.n)],
                      c_start=s.st, c_end=s.en, strand=["+" if x > 0 else "-" for x in s.sd], partial=s.partial, catalogue=s.cat, family=fam,
                      role=[role(g) for g in G], block=[bid[s.bpos[fam[g]][0]] if g in s.g_pos else None for g in G],
                      slot=[s.g_pos.get(g) for g in G],
                      instance=[s.g_inst.get(g) for g in G], slot_after=[s.g_slot.get(g) for g in G]),
                 f"{out}/genes.parquet", dict(catalogue=pa.int64()))
        Fm = sorted(s.pres)
        cls = lambda f: ("linker" if f in s.linker else "hub" if f in s.hub else "partial_only" if not s.nmag[f] else
                         "rare" if f not in s.anchor else "core_anchor" if f in s.core_anchor else "anchor")
        write_pq(dict(MGS=[m] * len(Fm), family=Fm, catalogue=[sorted(s.fam_cats[f]) for f in Fm],
                      n_catalogue=[len(s.fam_cats[f]) for f in Fm], n_mags=[s.nmag[f] for f in Fm],
                      prevalence_adj=[round(s.padj[f], 4) for f in Fm], max_copies=[s.maxc[f] for f in Fm],
                      multi_copy_mags=[s.multi[f] for f in Fm], median_len=[s.mlen.get(f) for f in Fm],
                      cls=[cls(f) for f in Fm], block=[bid[s.bpos[f][0]] if f in s.bpos else None for f in Fm],
                      slot=[s.bpos[f][1] if f in s.bpos else None for f in Fm], shared_mgs=[s.fam_shared.get(f, 1) for f in Fm]),
                 f"{out}/families.parquet", dict(catalogue=pa.list_(pa.int64())))
        write_pq(dict(cat_a=[x for x, _, _ in s.merges], cat_b=[y for _, y, _ in s.merges], identity=[z for _, _, z in s.merges]),
                 f"{out}/merged_catalogue_genes.parquet")
        Lk = sorted(s.links.items(), key=lambda kv: (-len(kv[1][0]), kv[0]))
        write_pq(dict(MGS=[m] * len(Lk), block_a=[bid[k[0][0]] for k, _ in Lk], slot_a=[k[0][1] for k, _ in Lk], end_a=[k[0][2] for k, _ in Lk],
                      block_b=[bid[k[1][0]] for k, _ in Lk], slot_b=[k[1][1] for k, _ in Lk], end_b=[k[1][2] for k, _ in Lk],
                      terminal=[v[2] for _, v in Lk], n_mags=[len(v[0]) for _, v in Lk], median_gap=[statistics.median_low(v[1]) for _, v in Lk],
                      MAGs=[",".join(sorted(s.mags[i] for i in v[0])) for _, v in Lk]),
                 f"{out}/links.parquet")
        def union_bp(iv):
            tot, a_, b_ = 0, None, None
            for x, y in sorted(iv):
                if b_ is None or x > b_ + 1:
                    tot += (b_ - a_ + 1) if b_ is not None else 0; a_, b_ = x, y
                else: b_ = max(b_, y)
            return tot + ((b_ - a_ + 1) if b_ is not None else 0)
        rows, s.mcov = [], collections.defaultdict(lambda: [0, 0])
        cg = {c: gl for _, c, gl in s.contigs}
        for mi, c in sorted((s.cmag[c], c) for c in s.clen):
            gl = cg.get(c, [])
            rs = sorted(s.by_ctg.get(c, []), key=lambda r: (r["start"], -r["end"]))
            path = " ".join(("[" if r["host"] >= 0 else "") + bid[r["block"]] + ("+" if r["strand"] > 0 else "-") + ("]" if r["host"] >= 0 else "") for r in rs)
            bc, ba = union_bp([(r["start"], r["end"]) for r in rs if tier(r["block"]) == "core"]), union_bp([(r["start"], r["end"]) for r in rs])
            s.mcov[mi][0] += bc; s.mcov[mi][1] += ba
            rows.append(dict(MGS=m, MAG=s.mags[mi], contig=c, length=s.clen[c], n_genes=len(gl), n_masked=sum(s.masked[g] for g in gl),
                             bp_core=bc, bp_blocks=ba, path=path))
        write_pq(cols(rows, ["MGS", "MAG", "contig", "length", "n_genes", "n_masked", "bp_core", "bp_blocks", "path"]), f"{out}/contigs.parquet")
        s.mbp = collections.Counter()
        for c, L in s.clen.items(): s.mbp[s.cmag[c]] += L
        nbk = collections.Counter(mg for mg, _ in {(r["mag"], r["block"]) for r in I})
        R = range(s.nm)
        write_pq(dict(MGS=[m] * s.nm, MAG=s.mags, samples=s.smp, is_MGS_rep=[s.MAG[g]["is_rep"] == "1" for g in s.mags],
                      completeness=[num(s.MAG[g]["completeness"]) for g in s.mags], contamination=[num(s.MAG[g]["contamination"]) for g in s.mags],
                      bp=[s.mbp[i] for i in R], bp_core=[s.mcov[i][0] for i in R], bp_blocks=[s.mcov[i][1] for i in R],
                      frac_core=[round(s.mcov[i][0] / max(1, s.mbp[i]), 4) for i in R], frac_blocks=[round(s.mcov[i][1] / max(1, s.mbp[i]), 4) for i in R],
                      n_blocks=[nbk[i] for i in R], mixed=[i in s.mixed for i in R]), f"{out}/mags.parquet")
        content = sorted((sum(r["end"] - r["start"] + 1 for r in s.by_block[bi]) for bi in nb), reverse=True)
        def n_for(fr):
            acc = 0
            for i, x in enumerate(content):
                acc += x
                if acc >= fr * sum(content): return i + 1
            return 0
        ids = [v[0] for i, v in s.ident.items() if s.I[i] is not s.rep[s.I[i]["block"]]]
        tbp = max(1, sum(s.mbp.values()))
        s.composition(out)                                   # + contig_seqs.parquet
        summ = dict(MGS=m, status="done" if B else "no_blocks", n_mags=s.nm, n_contigs=len(s.contigs), bp=sum(s.mbp.values()), n_genes=s.N,
                    n_masked=sum(s.masked), n_mixed=len(s.mixed), n_families=len(s.pres), n_merged=len(s.merges), n_linker=len(s.linker), n_hub=len(s.hub),
                    n_rare=sum(1 for f in s.pres if cls(f) == "rare"), n_core_anchor=len(s.core_anchor),
                    n_blocks_core=sum(b[3] == "core" for b in B), n_blocks_accessory=sum(b[3] != "core" for b in B),
                    bp_blocks=sum(map(len, s.seqs.values())), n_polished=sum(s.npol.values()),
                    frac_genes_in_seq=round(sum(s.rep[bi]["present"] for bi in nb) / max(1, sum(len(B[bi][0]) for bi in nb)), 4),
                    n_lift=len(s.seg), frac_bp_core=round(sum(v[0] for v in s.mcov.values()) / tbp, 4),
                    frac_bp_blocks=round(sum(v[1] for v in s.mcov.values()) / tbp, 4), blocks_for_50pct=n_for(0.5), blocks_for_90pct=n_for(0.9),
                    n_instances=len(I), median_identity=round(statistics.median(ids), 5) if ids else "", n_links=len(s.links),
                    edge_min=s.edge_min, family_min=s.fam_min, settings=s.a.settings_txt, mmt=s.a.mmt_version)
        write_text(f"{out}/summary.tsv", "\t".join(summ) + "\n" + "\t".join(map(str, summ.values())) + "\n")
        return summ

    def composition(s, out):
        """base counts of every MAG contig (contig + strand) in 1 kb windows: GC skew, cumulative skew (replichores,
        ori / ter) and composition context of any position can be derived from it ; contig_seqs.parquet: the sequence of
        every MAG contig (upper case, as in the assembly), written in batches (an MGS may hold Gb of contigs)"""
        import numpy as np, pyarrow as pa
        rows = collections.defaultdict(list)
        sch = pa.schema([(c, pq_type(c)) for c in ("MGS", "MAG", "contig", "length", "seq")])
        tmp = f"{out}/contig_seqs.parquet.{os.getpid()}.tmp"
        w, buf, nb = pq_writer(tmp, sch), collections.defaultdict(list), 0
        def flush():
            if buf["contig"]: w.write_table(pa.table({f.name: pa.array(buf[f.name], f.type) for f in sch}))
            buf.clear()
        for c in sorted(s.clen, key=lambda c: (s.cmag[c], c)):
            k = s.asm[s.MAG[s.mags[s.cmag[c]]]["assembly"]]
            if k["key"] not in s.fa: s.fa[k["key"]] = fasta_open(k["fasta"], k["fai"], k["gzi"])
            sq = s.fa[k["key"]].fetch(c).upper()
            for f, v in (("MGS", s.m), ("MAG", s.mags[s.cmag[c]]), ("contig", c), ("length", len(sq)), ("seq", sq)): buf[f].append(v)
            nb += len(sq)
            if nb >= 32 << 20: flush(); nb = 0
            q = np.frombuffer(sq.encode(), np.uint8)
            L, W = len(q), 1000
            n = (L + W - 1) // W
            x = np.zeros(n * W, np.uint8); x[:L] = q; x = x.reshape(n, W)
            st = np.arange(n) * W + 1
            rows["MGS"] += [s.m] * n; rows["MAG"] += [s.mags[s.cmag[c]]] * n; rows["contig"] += [c] * n
            rows["w_start"] += st.tolist(); rows["w_end"] += np.minimum(st + W - 1, L).tolist()
            for b, col in ((65, "a"), (67, "c"), (71, "g"), (84, "t")): rows[col] += (x == b).sum(1).tolist()
        flush(); w.close(); os.replace(tmp, f"{out}/contig_seqs.parquet")
        write_pq({k: rows[k] for k in ("MGS", "MAG", "contig", "w_start", "w_end", "a", "c", "g", "t")}, f"{out}/composition.parquet")


def same_prep(a):
    f = f"{a.prep}/prep.fp"
    if a.prep_fp and (open(f).read().strip() if os.path.isfile(f) else "") != a.prep_fp:
        die(f"{a.mgs}: {a.prep} was rebuilt by another run (other inputs or settings) meanwhile - nothing written, rerun")


def cmd_mgs(a):
    same_prep(a)
    x = MGS(a)
    out = f"{a.out}/.{a.mgs}.{os.getpid()}.tmp"; shutil.rmtree(out, ignore_errors=True); os.makedirs(out)
    try:
        x.families(); x.classes(); x.blocks(); x.instances(); x.representatives(); x.lift()
        summ = x.write(out)
        for f in x.fa.values(): f.close()
        write_text(f"{out}/.fp", a.fp + "\n")
        same_prep(a)
    except BaseException:
        shutil.rmtree(out, ignore_errors=True); raise
    final, old = f"{a.out}/{a.mgs}", f"{a.out}/.{a.mgs}.{os.getpid()}.old"
    if os.path.isdir(final): os.replace(final, old)     # swap by renames: readers never see a half-deleted folder
    os.replace(out, final)
    shutil.rmtree(old, ignore_errors=True)
    log(f"{a.mgs}: {summ['n_mags']} MAGs, {summ['n_families']} families ({summ['n_merged']} merged, {summ['n_linker']} linker, "
        f"{summ['n_hub']} hub) -> {summ['n_blocks_core']} core + {summ['n_blocks_accessory']} accessory blocks, {summ['n_instances']} "
        f"instances ; blocks cover {summ['frac_bp_blocks']:.3f} of the MAG bp (core {summ['frac_bp_core']:.3f})")


def cmd_summary(a):
    rows, head = [], None
    for d in sorted(os.listdir(a.out)):
        f = f"{a.out}/{d}/summary.tsv"
        if not d.startswith((".", "_")) and os.path.isfile(f"{a.out}/{d}/.fp") and os.path.isfile(f):
            h, v = open(f).read().rstrip("\n").split("\n")[:2]
            rows.append(dict(zip(h.split("\t"), v.split("\t"))))
    head = list(dict.fromkeys(k for r in rows for k in r))      # MGS done by different MMT versions: union of columns
    write_text(f"{a.out}/summary.tsv", "\t".join(head) + "\n" + "".join("\t".join(r.get(k, "") for k in head) + "\n" for r in rows)
               if rows else "")
    print(len(rows))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("prep"); p.add_argument("--idx", required=True); p.add_argument("--genecat", required=True)
    p.add_argument("--gc-id", default="95"); p.add_argument("--out", required=True); p.add_argument("--threads", type=int, default=4)
    p.add_argument("--check", type=int, default=200); p.add_argument("--samtools", default="samtools")
    p = sp.add_parser("plan"); p.add_argument("--idx", required=True); p.add_argument("--out", required=True)
    p.add_argument("--base-fp", required=True); p.add_argument("--mgs"); p.add_argument("--part"); p.add_argument("--force", action="store_true")
    p = sp.add_parser("mgs"); p.add_argument("--idx", required=True); p.add_argument("--prep", required=True)
    p.add_argument("--mgs", required=True); p.add_argument("--out", required=True); p.add_argument("--tmp", default="/tmp")
    p.add_argument("--threads", type=int, default=1); p.add_argument("--minimap2", default="")
    p.add_argument("--settings", required=True); p.add_argument("--settings-txt", default=""); p.add_argument("--fp", default="")
    p.add_argument("--mmt-version", default=""); p.add_argument("--prep-fp", default="")
    p.add_argument("--keep-tmp", action="store_true", help="keep the instance-alignment work folder (mmt --debug)")
    p = sp.add_parser("summary"); p.add_argument("--out", required=True)
    a = ap.parse_args()
    {"prep": cmd_prep, "plan": cmd_plan, "mgs": cmd_mgs, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    main()
