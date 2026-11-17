#!/usr/bin/env python3
"""check_motexport.py MOCK OUT.tsv MODE - checks of `mmt export-mot` against the mock's truth (motexport_mock.py).
MODE full   = the whole mock (PARTIAL expected) ; clean = map_clean.txt (COMPLETE expected) ;
     noqc   = full with --exclude-qcfail. Prints PASS/FAIL lines ; exit 1 if any check fails."""
import collections, gzip, hashlib, json, os, sys
import pyarrow.parquet as pq
import pysam

mock, table, mode = sys.argv[1:4]
T = json.load(open(f"{mock}/truth.json"))
lib = table[:-4] + ".motlib"
NP = NF = 0


def t(name, ok, detail=""):
    global NP, NF
    print(("PASS  " if ok else "FAIL  ") + name + ("" if ok else f"   [{detail}]"))
    NP += bool(ok); NF += not ok


rows = [l.rstrip("\n").split("\t") for l in open(table)]
head, rows = rows[0], [dict(zip(rows[0], r)) for r in rows[1:]]
t("table header = mot columns", head == ["sample", "bam", "ref", "role", "reads", "group"], head)
reads = {r["reads"]: gzip.open(r["reads"], "rt").read().split("\n")[:-1] for r in rows}
man = json.load(open(f"{lib}/manifest.json"))
want_status = {"full": "PARTIAL", "noqc": "PARTIAL", "clean": "COMPLETE"}[mode]
t(f"status {want_status}", man["status"] == want_status, man["status"])
val = json.load(open(f"{lib}/qc/validation.json"))
t("every validation check passed", all(c["ok"] for c in val["checks"]), [c for c in val["checks"] if not c["ok"]])
mags = {l.split("\t")[0]: dict(zip(open(f"{lib}/catalog/mags.tsv").readline().rstrip("\n").split("\t"),
                                   l.rstrip("\n").split("\t")))
        for l in list(open(f"{lib}/catalog/mags.tsv"))[1:]}
cat = {}
for f in os.listdir(f"{lib}/catalog/molecules"):
    if f.endswith(".parquet"):
        for r in pq.read_table(f"{lib}/catalog/molecules/{f}").to_pylist():
            cat[(r["sample"], f"{r['movie']}/{r['zmw']}")] = r
truth_sel = {k: v for k, v in T["selected"].items() if mode != "clean" or k.split("|")[0] in ("S1", "S2", "S3", "S4", "S5")}
if mode == "noqc":                                                  # S1's 0x200 molecules (ZMW 3, 4) leave the selection
    k = [x for x in truth_sel if x.startswith("S1|")][0]
    truth_sel[k] = [m for m in truth_sel[k] if m.split("/")[1] not in ("3", "4")]


def fa(path):
    return {r.name: r.sequence for r in pysam.FastxFile(path)}


def got(s, mag, lab):
    r = [x for x in rows if x["group"] == s and x["reads"].endswith(f"/{lab}.reads.txt.gz")
         and os.path.basename(os.path.dirname(x["reads"])) == mags[mag]["file"]]
    return reads[r[0]["reads"]] if r else []


# ---- the selections are exactly the truth
for k, v in sorted(truth_sel.items()):
    s, mag, lab = k.split("|")
    t(f"selection {s} {mag} {lab}: exactly the true molecules ({len(v)})", got(s, mag, lab) == sorted(
        v, key=lambda x: (x.split("/")[0], int(x.split("/")[1]))), (len(got(s, mag, lab)), len(v)))
t("no row beyond the truth", len(rows) == len(truth_sel), (len(rows), len(truth_sel)))

# A simple sample
r1 = [r for r in rows if r["group"] == "S1"]
t("A: S1 has one row (one BAM, one MAG) named S1__S1__SemiBin_1", [r["sample"] for r in r1] == ["S1__S1__SemiBin_1"])
t("A: S1's BAM resolved through the map (input_raw.txt names a moved file)",
  r1 and r1[0]["bam"].endswith("/raw/c1/S1.hifi_reads.bam"))
# B multi-contig MAG ; C circular
f2 = fa(mags["S2__SemiBin_1"]["fasta"])
src = fa(f"{mock}/mf4/run/S2/assemblies/metag/scaffolds.fasta.filt")
t("B: S2__SemiBin_1 FASTA = its 3 contigs, sequences identical to the assembly",
  sorted(f2) == sorted(T["mags"]["S2__SemiBin_1"]) and all(f2[c] == src[c] for c in f2))
t("B: one FASTA per MAG, 3 records, not concatenated", len(f2) == 3)
t("B: checked against GeneCat's own FASTA (representative)", mags["S2__SemiBin_1"]["checked_vs_genecat"] == "1")
f4 = fa(mags["S2__SemiBin_2"]["fasta"])
circ = [l.split("\t") for l in gzip.open(f"{lib}/catalog/contigs.tsv.gz", "rt")]
t("C: circular contig exported as in the assembly (no added overlap), flagged circular",
  f4 == {"S2__C4": src["S2__C4"]} and ["S2__SemiBin_2", "S2__C4"] == [c[0] for c in circ if c[1] == "S2__C4"] + ["S2__C4"]
  and [c[5].strip() for c in circ if c[1] == "S2__C4"] == ["1"] and mags["S2__SemiBin_2"]["n_circular"] == "1")
t("C: MAG without MGS ('?') is exported too", mags["S2__SemiBin_2"]["status"] == "ok" and mags["S2__SemiBin_2"]["MGS"] == "")
t("C: contig md5 in the catalog = md5 of the sequence", all(
    c[4] == hashlib.md5(src[c[1]].upper().encode()).hexdigest() for c in circ[1:] if c[1].startswith("S2__")))
# D three BAMs, one sample x MAG
r3 = [r for r in rows if r["group"] == "S3" and "SemiBin_1" in r["sample"]]
t("D: S3 x SemiBin_1 = 3 rows, one per BAM, same sample name and same ref",
  len(r3) == 3 and len({r["sample"] for r in r3}) == 1 and len({r["ref"] for r in r3}) == 1
  and len({r["bam"] for r in r3}) == 3)
t("D: each row's molecules are of its own BAM's movie", all(
    {x.split("/")[0] for x in reads[r["reads"]]} == {pysam.AlignmentFile(r["bam"], check_sq=False).header.to_dict()["RG"][0]["PU"]}
    for r in r3))
# E co-assembly
r45 = [r for r in rows if r["group"] in ("S4", "S5")]
t("E: co-assembled S4 and S5 -> two mot samples on one MAG FASTA, each with its own BAM",
  sorted(r["sample"] for r in r45) == ["S4__S5M2__SemiBin_1", "S5__S5M2__SemiBin_1"]
  and len({r["ref"] for r in r45}) == 1 and len({r["bam"] for r in r45}) == 2)
if mode != "clean":
    # F QC-flagged alignments kept with their status
    qc = {k[1]: v for k, v in cat.items() if k[0] == "S1"}
    zq = [v for m, v in qc.items() if m.split("/")[1] in ("3", "4")]
    zd = [v for m, v in qc.items() if m.split("/")[1] == "5"]
    t("F: 0x200 molecules in the catalog, flagged qc_fail, status unique",
      len(zq) == 2 and all(v["qc_fail"] and v["status"] == "unique" for v in zq))
    t("F: 0x400 molecule flagged duplicate", len(zd) == 1 and zd[0]["duplicate"])
    t(f"F: QC-flagged molecules {'not ' if mode == 'noqc' else ''}selected",
      all(v["selected"] == (mode != "noqc") for v in zq) and zd[0]["selected"])
    t("F: off-MAG molecules (low-quality bin, unbinned) not in the catalog, counted",
      not any(m.split("/")[1] in ("31", "36") for m in qc) and
      [l for l in open(f"{lib}/qc/export_summary.tsv") if l.startswith("molecules off MAGs")][0].split("\t")[1].strip() != "0")
    # G / H / statuses
    for k, st in T["status"].items():
        s, m = k.split("|")
        t(f"status {s} {m} = {st}", cat.get((s, m), {}).get("status") == st, cat.get((s, m), {}).get("status"))
    sp = [v for k, v in cat.items() if k[0] == "S3" and v["status"] == "split_across_mags"][0]
    t("G: split molecule names the other MAG and is in no selection",
      sp["other_MAGs"] in ("S3__SemiBin_3", "S3__SemiBin_2") and not sp["selected"]
      and f"{sp['movie']}/{sp['zmw']}" not in {x for v in reads.values() for x in v})
    h = [k[1] for k, v in cat.items() if k[0] == "S3" and v["n_supplementary"] == 1 and v["status"] == "unique"]
    t("H: primary + supplementary molecule listed once", len(h) == 1 and sum(x == h[0] for v in reads.values() for x in v) == 1)
    bs = [v for k, v in cat.items() if k[0] == "S3" and v["n_primary"] == 2]
    t("H: by-strand fwd + rev records = one molecule", len(bs) == 1 and bs[0]["selected"])
    allsel = [x for v in reads.values() for x in v]
    t("H: no molecule twice in the selections of one BAM", len(allsel) == len(set(allsel)))
    un = open(f"{lib}/qc/unresolved.tsv").read()
    # I missing BAM ; J other technologies ; K changed reference
    t("I: S6 (no raw BAM) reported, no row", "S6\t" in un and not [r for r in rows if r["group"] == "S6"])
    t("I: S6 molecules in the catalog with bam_status movie_not_in_bams",
      all(v["bam_status"] == "movie_not_in_bams" for k, v in cat.items() if k[0] == "S6") and any(k[0] == "S6" for k in cat))
    t("J: Illumina (S7) and ONT (S8) reads: no row, not reported as a problem, no molecule in the catalog",
      not [r for r in rows if r["group"] in ("S7", "S8")] and "S7\t" not in un and "S8\t" not in un
      and not any(k[0] in ("S7", "S8") for k in cat))
    t("J: PacBio BAM without kinetics (S9) reported, no row", "S9\t" in un and not [r for r in rows if r["group"] == "S9"])
    t("K: S10 contig changed after mapping -> reported, its molecules reference_mismatch, never selected",
      "S10\t" in un and all(not v["selected"] for k, v in cat.items() if k[0] == "S10" and v["status"] == "reference_mismatch"))
    smp = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{lib}/catalog/samples.tsv")}
    t("retained alignments evidence recorded per CRAM", smp["S1"][12].startswith("yes"))
# references usable by mot: bgzip + fai
for r in rows[:1]:
    t("ref is bgzip with .fai and .gzi", open(r["ref"], "rb").read(4)[:2] == b"\x1f\x8b"
      and os.path.isfile(r["ref"] + ".fai") and os.path.isfile(r["ref"] + ".gzi"))
print(f"{NP} passed, {NF} failed")
sys.exit(1 if NF else 0)
