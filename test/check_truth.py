#!/usr/bin/env python3
"""check_truth.py MOCKDIR OUT INDEXDIR MODE... : compare MMT MGScall results on the mock run (mock_mf4.py) with the truth
(INDEXDIR: the index cache in the work space, <tmp>/mmt/<OUT>_<key>/index).
Prints one line per check (PASS / FAIL) and exits 1 if any check failed. Needs duckdb on PATH or $DUCKDB_PY (python duckdb
not required: reads parquet with pyarrow if present, else with the duckdb CLI given in $DUCKDB_BIN)."""
import csv, glob, io, json, os, subprocess, sys

mock, out, idx, modes = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
T = json.load(open(os.path.join(mock, "truth.json")))
DUCK = os.environ.get("DUCKDB_BIN", "duckdb")
fails = 0


def check(ok, msg):
    global fails
    print(("PASS  " if ok else "FAIL  ") + msg)
    fails += not ok


def q(sql):   # duckdb CLI -> list of dict rows
    r = subprocess.run([DUCK, "-csv", "-c", sql], capture_output=True, text=True, check=True).stdout
    return list(csv.DictReader(io.StringIO(r)))


COMP = str.maketrans("ACGT", "TGCA")
def motif_sites(seqs, names, motif):   # (contig, pos 1-based, strand) of the A at motif position 2, both strands
    s = set()
    rcm = motif.translate(COMP)[::-1]
    for name, seq in zip(names, seqs):
        i = seq.find(motif)
        while i >= 0:
            s.add((name, i + 2, "+")); i = seq.find(motif, i + 1)
        i = seq.find(rcm)
        while i >= 0:
            s.add((name, i + len(motif) - 1, "-")); i = seq.find(rcm, i + 1)
    return s


SPECIES = {"MGS.0": {"A", "A2", "A3"}, "MGS.1": {"B"}}
REP = {"MGS.0": (T["seq"]["A"], ["S1__C1", "S1__C2"], "GAAG"),
       "MGS.1": (T["seq"]["B"], ["S2.run2M2__C3", "S2.run2M2__C4"], "GATC")}
G2A, G2B = "S2.run2M2__SemiBin_3", "S2.run2M2__SemiBin_5.fa.gz"
EXPECT = {"mgs-rep": {("MGS.0", s, "") for s in ("S1", "S2", "S2.run2", "S3")} | {("MGS.1", s, "") for s in ("S1", "S2", "S2.run2")},
          "mag-rep": {("MGS.0", "S1", "S1__SemiBin_1"), ("MGS.0", "S1", "S1__SemiBin_6"), ("MGS.0", "S2", G2A),
                      ("MGS.0", "S2.run2", G2A), ("MGS.0", "S3", "S3__SemiBin_2"), ("MGS.1", "S1", "S1__SemiBin_4"),
                      ("MGS.1", "S2", G2B), ("MGS.1", "S2.run2", G2B)}}
EXPECT["mag-mag"] = EXPECT["mag-rep"]
# mag-mag: every MAG's own contigs (index) and their sequences (the mock assemblies)
MAGCTG = {}
for r in q(f"SELECT MAG, contig FROM read_csv('{idx}/contigs.tsv.gz', delim='\\t', header=true, all_varchar=true) ORDER BY contig"):
    MAGCTG.setdefault(r["MAG"], []).append(r["contig"])
SEQ = {}
for fa in glob.glob(os.path.join(mock, "**", "scaffolds.fasta.filt"), recursive=True):
    for blk in open(fa).read().split(">")[1:]:
        h, _, b = blk.partition("\n"); SEQ[h.split()[0]] = b.replace("\n", "").upper()
MAGS = {("MGS.0", "S1"): "S1__SemiBin_1,S1__SemiBin_6", ("MGS.0", "S2"): G2A, ("MGS.0", "S2.run2"): G2A,
        ("MGS.0", "S3"): "S3__SemiBin_2", ("MGS.1", "S1"): "S1__SemiBin_4", ("MGS.1", "S2"): G2B, ("MGS.1", "S2.run2"): G2B}


def parts_of(mode):   # result files -> {(MGS, sample, MAG or ""): path}
    D = os.path.join(out, mode)
    if mode != "mgs-rep":
        return {(p.split("/")[-3], p.split("/")[-2], os.path.basename(p)[:-len(".parquet")]): p
                for p in glob.glob(os.path.join(D, "*", "*", "*.parquet"))}
    return {(p.split("/")[-2], os.path.basename(p)[:-len(".parquet")], ""): p for p in glob.glob(os.path.join(D, "MGS*", "*.parquet"))}

# ---- read -> MAG index
for f in sorted(glob.glob(os.path.join(idx, "reads", "*.parquet"))):
    s = os.path.basename(f)[:-len(".parquet")]
    rows = q(f"SELECT read, MGS, MAG, is_rep, dup FROM read_parquet('{f}')")
    wrong = [r for r in rows if T["reads"].get(r["read"]) not in SPECIES.get(r["MGS"], ())]
    flagged = [r for r in rows if T["flag"].get(r["read"], 0) & 0x200]
    dupok = all((r["dup"] == "true") == bool(T["flag"].get(r["read"], 0) & 0x400) for r in rows)
    check(not wrong and not flagged and dupok and rows, f"reads {s}: {len(rows)} assigned, {len(wrong)} to a wrong species, "
          f"{len(flagged)} QC-fail kept, duplicate flags {'right' if dupok else 'WRONG'}")
    for m, sp in SPECIES.items():   # every unflagged, primary-mapped truth read of the species is found
        want = {r for r, n in T["reads"].items() if r.startswith(f"m84_{s}/") and n in sp
                and r in T["flag"] and not T["flag"][r] & 0x200}
        got = {r["read"] for r in rows if r["MGS"] == m}
        if want:
            check(len(got & want) >= 0.97 * len(want), f"reads {s} {m}: {len(got & want)}/{len(want)} truth reads assigned")

for mode in modes:
    D = os.path.join(out, mode)
    P = parts_of(mode)
    parts = set(P)
    check(parts == EXPECT[mode], f"{mode}: results for {sorted(parts)}")
    glob_ = f"{D}/*/*/*.parquet" if mode != "mgs-rep" else f"{D}/MGS*/*.parquet"
    cols = q(f"DESCRIBE SELECT * FROM read_parquet('{glob_}')")
    want_cols = ["MGS", "MAG", "sample", "contig", "pos", "strand", "base", "coverage", "score", "ipdRatio", "tMean", "tErr",
                 "modelPrediction", "frac", "fracLow", "fracUp", "call", "identificationQv", "pile_depth", "n_ref", "top_alt",
                 "n_alt", "ctx_alt_pct", "ctx_alt_n", "clip_ctx_n"]
    check([c["column_name"] for c in cols] == want_cols, f"{mode}: columns {[c['column_name'] for c in cols]}")
    fpat = "'%/' || MGS || '/' || sample || '/' || MAG || '.parquet'" if mode != "mgs-rep" else "'%/' || MGS || '/' || sample || '.parquet'"
    bad = q(f"""SELECT count(*) n FROM read_parquet('{glob_}', filename = true)
               WHERE (base NOT IN ('A','C') AND call = 'none') OR base NOT IN ('A','C','G','T') OR coverage < 3
                  OR strand NOT IN ('+','-') OR pile_depth < coverage OR filename NOT LIKE {fpat}""")[0]["n"]
    check(bad == "0", f"{mode}: rows of a G/T without a kinetic call, coverage < 3, bad strand, pile_depth < coverage or "
          f"MGS/sample/MAG not matching the file: {bad}")
    gt = q(f"SELECT count(*) n FROM read_parquet('{glob_}') WHERE base IN ('G','T')")[0]["n"]
    check(int(gt) > 0, f"{mode}: G/T rows with a kinetic signal (modified_base) kept: {gt}")
    if mode == "mgs-rep":
        mg = {(r["MGS"], r["sample"]): r["MAG"] for r in q(f"SELECT DISTINCT MGS, sample, MAG FROM read_parquet('{glob_}')")}
        check(mg == MAGS, f"mgs-rep: MAG column = every MAG of the MGS in the sample")
    uns = q(f"""SELECT count(*) n FROM (SELECT filename, pos, contig, lag(contig) OVER w pc, lag(pos) OVER w pp
               FROM read_parquet('{glob_}', filename = true, file_row_number = true)
               WINDOW w AS (PARTITION BY filename ORDER BY file_row_number)) WHERE contig = pc AND pos < pp""")[0]["n"]
    check(uns == "0", f"{mode}: rows out of contig / position order inside files: {uns}")
    meta = q(f"SELECT DISTINCT key, value FROM parquet_kv_metadata('{glob_}')")
    check(any(r["key"] == "mmt" and "1-based" in r["value"] for r in meta), f"{mode}: column conventions in the file metadata (pos 1-based)")
    # state files agree with the results
    for s in sorted({s for _, s, _ in EXPECT[mode]}):
        st = list(csv.DictReader(open(os.path.join(D, "stats", f"{s}.tsv")), delimiter="\t"))
        done = {(r["MGS"], s, r["MAG"] if mode != "mgs-rep" else "") for r in st if r["status"] in ("done", "linked")}
        check(done == {p for p in parts if p[1] == s}, f"{mode} {s}: state done/linked = {sorted(done)}")
    # methylation at the motif = truth. ipdSummary gives frac only where it identified m6A, so the fraction is also
    # read from the mean ipdRatio: ipdRatio - 1 grows linearly with the methylated fraction (scaled to a 0.9 sample)
    stat = {}
    for (m, s, g), path in sorted(P.items()):
        seq, names, motif = REP[m]
        if mode == "mag-mag":   # the MAG's own contigs
            names = MAGCTG[g]; seq = [SEQ[c] for c in names]
        sites = motif_sites(seq, names, motif)
        rows = q(f"""SELECT contig, pos, strand, frac, call, ipdRatio FROM read_parquet('{path}')
                    WHERE base = 'A' AND coverage >= 11 AND ctx_alt_pct = 0""")
        on = [r for r in rows if (r["contig"], int(r["pos"]), r["strand"]) in sites]
        off = [r for r in rows if (r["contig"], int(r["pos"]), r["strand"]) not in sites]
        fr = 0.9 if m == "MGS.1" else T["FR"][s]
        fs = [float(r["frac"]) for r in on if r["frac"]]
        stat[(m, s, g)] = dict(n=len(on), fr=fr, ipd=sum(float(r["ipdRatio"]) for r in on) / max(1, len(on)),
                               frac=sum(fs) / max(1, len(fs)), det=sum(r["call"] != "none" for r in on) / max(1, len(on)),
                               false=sum(r["call"] == "m6A" for r in off) / max(1, len(off)), noff=len(off))
    for (m, s, g), v in sorted(stat.items()):
        if v["n"] <= 50 and g == "S1__SemiBin_6":   # the A3 strain (one contig, 10x): too few clean motif sites at >= 11 reads
            check(True, f"{mode} {m} {s} {g}: {v['n']} motif A (low-depth strain, not scored)"); continue
        ref = [w["ipd"] for (mm, _, _), w in stat.items() if mm == m and w["fr"] == 0.9 and w["n"] > 50]
        est = 0.9 * (v["ipd"] - 1) / (sum(ref) / len(ref) - 1)
        ok = v["n"] > 50 and abs(est - v["fr"]) < 0.1 and v["false"] < 0.01 and v["det"] > (0.85 if v["fr"] == 0.9 else 0.6)
        ok = ok and (v["fr"] != 0.9 or v["frac"] > 0.85)
        check(ok, f"{mode} {m} {s}{' ' + g if g else ''}: {v['n']} motif A, fraction from ipdRatio {est:.2f} (truth {v['fr']}), with a call "
                  f"{v['det']:.2f}, mean frac where given {v['frac']:.2f}, m6A at {v['false']:.4f} of {v['noff']} other A")

# ---- strain differences: S3's reads (strain A2) on the MGS.0 representative
if "mgs-rep" in modes:
    A, A2 = T["seq"]["A"], T["seq"]["A2"]
    ia, il = T["island"]
    snp = {("S1__C1", i + 1) for i in range(len(A[0])) if A2[0][i if i < ia else i + il] != A[0][i]}
    snp |= {("S1__C2", i + 1) for i in range(len(A[1])) if A2[1][i] != A[1][i]}
    f = os.path.join(out, "mgs-rep", "MGS.0", "S3.parquet")
    rows = q(f"SELECT contig, pos, strand, pile_depth, n_ref, ctx_alt_pct, clip_ctx_n FROM read_parquet('{f}') WHERE coverage >= 5")
    at = [r for r in rows if (r["contig"], int(r["pos"])) in snp]
    flagged = sum(int(r["n_ref"]) <= 0.2 * int(r["pile_depth"]) for r in at) / max(1, len(at))
    near = {(c, p + d) for c, p in snp for d in range(-10, 11)}
    far = [r for r in rows if (r["contig"], int(r["pos"])) not in near and abs(int(r["pos"]) - ia) > 20]
    clean = sum(r["ctx_alt_pct"] == "0" for r in far) / max(1, len(far))
    ctx = sum(int(r["ctx_alt_pct"]) >= 50 for r in rows if (r["contig"], int(r["pos"])) in near) / \
        max(1, sum(1 for r in rows if (r["contig"], int(r["pos"])) in near))
    junction = [r for r in rows if r["contig"] == "S1__C1" and abs(int(r["pos"]) - ia) <= 5]
    check(len(at) > 100 and flagged > 0.95, f"S3 strain SNPs: {len(at)} rows at SNPs, {flagged:.3f} show the other base")
    check(ctx > 0.95 and clean > 0.99, f"S3 context flag: {ctx:.3f} of rows within 10 bp of a SNP flagged, "
          f"{clean:.4f} of {len(far)} rows far from SNPs clean")
    check(any(int(r["clip_ctx_n"]) > 0 for r in junction), f"S3 island junction: reads clipped at representative position ~{ia}")
    # --mgs-rep S1 pools the A strain (representative bin) and the A3 strain (own bin) ; --mag-rep keeps them apart
    if "mag-rep" in modes:
        mol = {}
        for mode in ("mgs-rep", "mag-rep"):
            for r in csv.DictReader(open(os.path.join(out, mode, "stats", "S1.tsv")), delimiter="\t"):
                if r["MGS"] == "MGS.0":
                    mol[(mode, r["MAG"])] = int(r["molecules"])
        nA3 = sum(1 for r, n in T["reads"].items() if r.startswith("m84_S1/") and n == "A3" and not T["flag"].get(r, 0) & 0x600)
        a3, a1 = mol.get(("mag-rep", "S1__SemiBin_6"), 0), mol.get(("mag-rep", "S1__SemiBin_1"), 0)
        check(mol.get(("mgs-rep", "S1__SemiBin_1,S1__SemiBin_6")) == a1 + a3 and abs(a3 - nA3) <= 0.03 * nA3,
              f"S1 MGS.0 molecules: mgs-rep {mol.get(('mgs-rep', 'S1__SemiBin_1,S1__SemiBin_6'))} = mag-rep SemiBin_1 {a1} + "
              f"SemiBin_6 {a3} (A3 strain reads, not QC-fail / duplicate: {nA3})")
        # the A3 MAG's table: its reads cover contig 2 only and carry A3's SNPs
        f3 = os.path.join(out, "mag-rep", "MGS.0", "S1", "S1__SemiBin_6.parquet")
        r3 = q(f"SELECT contig, count(*) n FROM read_parquet('{f3}') WHERE coverage >= 5 GROUP BY contig")
        n3 = {r["contig"]: int(r["n"]) for r in r3}
        A3 = T["seq"]["A3"][0]
        snp3 = {("S1__C2", i + 1) for i in range(len(A[1])) if A3[i] != A[1][i]}
        rows3 = q(f"SELECT contig, pos, pile_depth, n_ref FROM read_parquet('{f3}') WHERE coverage >= 5")
        at3 = [r for r in rows3 if (r["contig"], int(r["pos"])) in snp3]
        fl3 = sum(int(r["n_ref"]) <= 0.2 * int(r["pile_depth"]) for r in at3) / max(1, len(at3))
        check(n3.get("S1__C2", 0) > 1000 and n3.get("S1__C1", 0) < 0.02 * n3.get("S1__C2", 1) and len(at3) > 50 and fl3 > 0.95,
              f"mag-rep S1 SemiBin_6 (A3 strain): rows on contig 2 {n3.get('S1__C2', 0)}, contig 1 {n3.get('S1__C1', 0)} ; "
              f"{len(at3)} rows at A3 SNPs, {fl3:.3f} show the other base")
# ---- mag-mag: each MAG's reads on its own contigs
if "mag-mag" in modes:
    G = f"{out}/mag-mag/*/*/*.parquet"
    off = q(f"SELECT DISTINCT MAG, contig FROM read_parquet('{G}')")
    check(off and all(r["contig"] in MAGCTG.get(r["MAG"], ()) for r in off),
          f"mag-mag: every table on its own MAG's contigs ({len(off)} MAG x contig pairs)")
    rows = q(f"SELECT contig, pos, strand, base FROM read_parquet('{G}')")
    bad = sum((SEQ[r["contig"]][int(r["pos"]) - 1] if r["strand"] == "+" else SEQ[r["contig"]][int(r["pos"]) - 1].translate(COMP)) != r["base"]
              for r in rows)
    check(rows and bad == 0, f"mag-mag: base = the MAG contig's base (- strand: complement) on {len(rows)} rows, {bad} differ")
    a = f"{out}/mag-rep/MGS.0/S1/S1__SemiBin_1.parquet"; b = f"{out}/mag-mag/MGS.0/S1/S1__SemiBin_1.parquet"
    # fracLow / fracUp: ipdSummary's bootstrap interval - random, differs between any two runs
    x = "SELECT * EXCLUDE (fracLow, fracUp) FROM read_parquet"
    d = q(f"SELECT (SELECT count(*) FROM ({x}('{a}') EXCEPT ALL {x}('{b}'))) + (SELECT count(*) FROM ({x}('{b}') EXCEPT ALL {x}('{a}'))) n")[0]["n"]
    check(d == "0", f"mag-mag = mag-rep for the representative MAG in its own sample (S1 SemiBin_1, same reads and sequence; "
                    f"fracLow/fracUp left out): {d} rows differ")
    for g, f in (("S1__SemiBin_6", "MGS.0/S1/S1__SemiBin_6"), ("S3__SemiBin_2", "MGS.0/S3/S3__SemiBin_2")):
        r = q(f"SELECT count(*) n, count(*) FILTER (WHERE n_ref <= 0.2 * pile_depth) alt, count(*) FILTER (WHERE ctx_alt_pct >= 50) ctx "
              f"FROM read_parquet('{out}/mag-mag/{f}.parquet') WHERE coverage >= 5")[0]
        check(int(r["n"]) > 1000 and int(r["alt"]) <= 0.001 * int(r["n"]) and int(r["ctx"]) <= 0.002 * int(r["n"]),
              f"mag-mag {g} (a strain other than the representative) on itself: {r['n']} rows, {r['alt']} show another base, "
              f"{r['ctx']} context-flagged")
print(f"{'ALL PASS' if not fails else f'{fails} FAILED'}")
sys.exit(1 if fails else 0)
