#!/usr/bin/env python3
"""Self-check of mmt_landscape.py on a mock MMT output written with MMT's own schema (lib/mmtpq.py).

  python3 analysis/test_landscape.py [--duckdb BIN] [--keep DIR]

Planted truth: 3 MGS, 10 samples (2 chemistries S01-S06 / S07-S10), 8 patients (S06 + S07 and S09 + S10 = one person
each, S08 missing from the map), --run pairs S09~S10 in-batch and S06~S07 across-batch,
coverage 2x to 90x. MGS.1: the same 8 % of A sites (m6A) and 3 % of C sites (m4C) methylated in every sample, a 5-kb
A island (40 %), sample S06 half-methylated. MGS.2: methylated sites drawn anew in every patient (null for agreement and
constitutive sites ; S07 / S10 repeat S06's / S09's draw = technical replicates, so between-patient change >> technical
bands in MGS.2, not in MGS.1: fig7, --between 3). MGS.3: one sample. Detection 1 - exp(-coverage * frac / 6), masked rows, G/T modified_base rows
(must never be counted). Every summary number is re-derived here with pandas from the parquet files and must match
exactly; the tests must find what was planted and nothing else; chunking must not change any result.
"""
import argparse, glob, os, shutil, stat, subprocess, sys, tempfile
import numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq
from scipy.stats import binom

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import mmtpq  # noqa: E402

COMP = str.maketrans("ACGT", "TGCA")
LEVELS = np.float32([0.05, 0.25, 0.75, 0.95])
FAILS = []
MC = 11


def check(ok, what):
    print(("PASS " if ok else "FAIL ") + what)
    if not ok:
        FAILS.append(what)


def mock(out, rng):
    samples = [f"S{i:02d}" for i in range(1, 11)]
    lam = dict(zip(samples, [2, 4, 8, 12, 18, 25, 35, 50, 70, 90]))          # mean coverage per strand
    chem = {s: ("102-739-100/102-118-800" if i < 6 else "103-496-900/103-496-700") for i, s in enumerate(samples)}
    genomes = {"MGS.1": {"c1": 120000, "c2": 20000}, "MGS.2": {"k1": 30000}, "MGS.3": {"z1": 8000}}
    present = {"MGS.1": samples, "MGS.2": samples[2:], "MGS.3": ["S09"]}
    stats = {s: [] for s in samples}
    truth, drawn = {}, {}
    for mgs, ctgs in genomes.items():
        seqs = {c: "".join(rng.choice(list("ACGT"), L)) for c, L in ctgs.items()}
        sites = []
        for c, sq in seqs.items():
            for i, b in enumerate(sq):
                sites.append((c, i + 1, "+", b))
                sites.append((c, i + 1, "-", b.translate(COMP)))
        S = pd.DataFrame(sites, columns=["contig", "pos", "strand", "base"])
        u = rng.random(len(S))
        S["meth"] = ((S["base"] == "A") & (u < 0.08)) | ((S["base"] == "C") & (u < 0.03))
        if mgs == "MGS.1":   # a 5-kb methylation island on c1 (spatial clustering)
            isl = (S["contig"] == "c1") & (S["pos"] <= 5000) & (S["base"] == "A")
            S.loc[isl, "meth"] = u[isl.values] < 0.40
        S["site_eff"] = rng.gamma(4, 0.25, len(S))                          # per-site coverage variation
        truth[mgs] = S[["contig", "pos", "strand", "base", "meth"]]
        os.makedirs(os.path.join(out, "mgs-rep", mgs), exist_ok=True)
        for s in samples:
            if s not in present[mgs]:
                continue
            d = S.copy()
            if mgs == "MGS.2":   # null: methylated sites drawn anew in every patient (replicates repeat the draw)
                v = drawn.setdefault({"S10": "S09", "S07": "S06"}.get(s, s), rng.random(len(d)))
                d["meth"] = ((d["base"] == "A") & (v < 0.08)) | ((d["base"] == "C") & (v < 0.03))
            d["coverage"] = rng.poisson(lam[s] * d["site_eff"])
            d = d[d["coverage"] >= 3].copy()
            tf = np.where(d["meth"], 0.5 if (s == "S06" and mgs == "MGS.1") else 1.0, 0.0)
            pcall = np.where(d["meth"], 1 - np.exp(-d["coverage"] * tf / 6), 0.002)
            called = rng.random(len(d)) < pcall
            possible = (~called) & (rng.random(len(d)) < np.where(d["meth"], 0.3, 0.01))
            d["call"] = np.where(called, np.where(d["base"] == "A", "m6A", np.where(d["base"] == "C", "m4C", "none")),
                                 np.where(possible, "modified_base", "none"))
            gt = d["base"].isin(["G", "T"])
            d.loc[gt, "call"] = np.where(rng.random(gt.sum()) < 0.005, "modified_base", "none")
            d = d[~gt | (d["call"] != "none")]                                # MMT keeps G/T only with a call
            isc = d["call"].isin(["m6A", "m4C"])
            tf = np.where(d["meth"], 0.5 if (s == "S06" and mgs == "MGS.1") else 1.0, 0.0)
            d["frac"] = np.where(isc & (d["coverage"] > 10), np.round(np.clip(tf + rng.normal(0, 0.03, len(d)), 0, 1), 3), np.nan)
            d["score"] = np.where(isc, 20 + 2 * d["coverage"] * tf, rng.integers(0, 15, len(d))).astype(int)
            d["identificationQv"] = np.where(isc, 30.0, np.nan)
            d["pile_depth"] = d["coverage"]
            d["n_ref"] = d["coverage"]
            m = rng.random(len(d))
            d["ctx_alt_pct"] = np.where(m < 0.03, 40, 0)
            d.loc[(m >= 0.03) & (m < 0.04), "n_ref"] = (d["coverage"] * 0.5).astype(int)
            d["clip_ctx_n"] = np.where((m >= 0.04) & (m < 0.05), 2, 0)
            d["top_alt"], d["n_alt"], d["ctx_alt_n"] = "", d["pile_depth"] - d["n_ref"], 0
            d["ipdRatio"] = 1.0; d["tMean"] = 1.0; d["tErr"] = 0.1; d["modelPrediction"] = 1.0
            d["fracLow"] = d["frac"]; d["fracUp"] = d["frac"]
            d["MGS"], d["MAG"], d["sample"] = mgs, f"{s}__bin1", s
            d = d.sort_values(["contig", "pos", "strand"])
            sch = mmtpq.schema("mgscall")
            t = pa.Table.from_pandas(d[[f.name for f in sch]], preserve_index=False)
            mmtpq.write(os.path.join(out, "mgs-rep", mgs, s + ".parquet"), mmtpq.cast(t, sch), sch)
            stats[s].append([mgs, f"{s}__bin1", "done", "fp", "1", "100", "100", "100", str(len(d)), f"S01__{mgs}", "now", "0.6.0"])
    stats["S01"].append(["MGS.9", "S01__bin9", "no_reads", "fp", "1", "0", "0", "0", "0", "S01__bin9", "now", "0.6.0"])
    os.makedirs(os.path.join(out, "mgs-rep", "stats"), exist_ok=True)
    for s in samples:
        with open(os.path.join(out, "mgs-rep", "stats", s + ".tsv"), "w") as o:
            o.write("MGS\tMAG\tstatus\tfingerprint\tabundance\tmolecules\tpseudo_reads\taligned\trows\trep_MAG\tupdated\tmmt\n")
            o.writelines("\t".join(r) + "\n" for r in stats[s])
        k1, k2 = chem[s].split("/")
        with open(os.path.join(out, "mgs-rep", "stats", s + ".run.tsv"), "w") as o:
            o.write(f"sample\t{s}\nmode\tmgs-rep\nmovie\tm84000_{s}\nbinding_kit\t{k1}\nsequencing_kit\t{k2}\n"
                    f"kinetics_model\tSP2-C2\n")
    return samples, truth


def write_map(path, pat):
    with open(path, "w") as o:
        o.write("#OutPath\t/x/\n#RunID\tT\n#SmplID\tPath\tSmplPrefix\tAssmblGrps\n")
        o.writelines(f"{s}\t/raw/\t{s}\t{p}\n" for s, p in pat.items())


def rows_of(out, mgs):
    """every A / C row (>= 3 reads) of an MGS with clean, category, class - re-derived from the parquet files"""
    d = pd.concat([pq.read_table(f).to_pandas() for f in sorted(glob.glob(os.path.join(out, "mgs-rep", mgs, "*.parquet")))])
    d = d[d["base"].isin(["A", "C"]) & (d["coverage"] >= 3)].copy()
    d["clean"] = ((d["pile_depth"] > 0) & (d["n_ref"] >= 0.95 * d["pile_depth"]) & (d["ctx_alt_pct"] < 10) &
                  (d["clip_ctx_n"] == 0))
    want = np.where(d["base"] == "A", "m6A", "m4C")
    d["cat"] = np.where(d["call"] == want, "called", np.where(d["call"] == "none", "none",
                        np.where(d["call"] == "modified_base", "possible", "other")))
    lv = np.searchsorted(LEVELS, d["frac"].astype(np.float32).values, side="right").astype(float)
    g = d["clean"] & (d["coverage"] >= MC)
    d["cls"] = np.where(g & (d["cat"] == "none"), 0, np.where(g & (d["cat"] == "called") & d["frac"].notna(), lv, np.nan))
    d["graded"] = g
    return d


def units_ref(d):
    rows = []
    for (mgs, s, b), x in d.groupby(["MGS", "sample", "base"]):
        c = x[x["clean"]]
        gr = c[c["coverage"] >= MC]
        r = dict(MGS=mgs, sample=s, base=b, sites=len(x), masked=int((~x["clean"]).sum()), clean=len(c),
                 tier_low=int((c["coverage"] < MC).sum()), tier_mid=int(((c["coverage"] >= MC) & (c["coverage"] < 20)).sum()),
                 tier_high=int((c["coverage"] >= 20).sum()), not_graded=int((c["coverage"] < MC).sum()),
                 possible=int(gr["cat"].isin(["possible", "other"]).sum()),
                 called_no_frac=int(((gr["cat"] == "called") & gr["frac"].isna()).sum()))
        for k in range(5):
            r[f"cls{k}"] = int((gr["cls"] == k).sum())
        cf = c[(c["cat"] == "called") & c["frac"].notna()]
        n, p, z = cf["coverage"].astype(float), cf["frac"].astype(np.float32).astype(float), 1.96
        cen = (p + z * z / (2 * n)) / (1 + z * z / n)
        h = z * np.sqrt(np.clip(p * (1 - p), 0, None) / n + z * z / (4 * n * n)) / (1 + z * z / n)
        lvl = lambda v: np.searchsorted(LEVELS.astype(float), v.values, side="right")
        span = np.clip(lvl(cen + h) - lvl(cen - h) + 1, None, 3)
        for k in (1, 2, 3):
            r[f"span{k}"] = int((span == k).sum())
        rows.append(r)
    return pd.DataFrame(rows)


def sites_ref(d, rank, patient):
    """per site: the best-covered graded sample of every patient (ties: sample name order) -> n, k0..k4, e4, e0"""
    g = d[d["graded"]].copy()
    g["code"] = g["cls"].fillna(9).astype(int)
    g["patient"] = g["sample"].map(patient)
    g["rank"] = g["sample"].map(rank)
    g = g.sort_values(["coverage", "rank"], ascending=[False, True])
    pr = g.drop_duplicates(["patient", "contig", "pos", "strand", "base"]).copy()
    pr["st"] = np.select([pr["coverage"] >= 50, pr["coverage"] >= 30, pr["coverage"] >= 20, pr["coverage"] >= 15], [50, 30, 20, 15], MC)
    v = pr[pr["code"] != 9]
    pi = v.groupby(["patient", "base", "st"]).agg(p4=("code", lambda c: (c == 4).mean()), p0=("code", lambda c: (c == 0).mean()))
    v = v.join(pi, on=["patient", "base", "st"])
    key = ["base", "contig", "pos", "strand"]
    sc = v.groupby(key).agg(n=("code", "size"), e4=("p4", "sum"), e0=("p0", "sum"), max_cov=("coverage", "max"))
    for k in range(5):
        sc[f"k{k}"] = v.assign(x=v["code"] == k).groupby(key)["x"].sum()
    sc["n_possible"] = pr[pr["code"] == 9].groupby(key).size()
    sc = sc.fillna({"n_possible": 0}).reset_index()
    return sc[sc["n"] >= 2], pr


def alpha_ref(sc):
    """Krippendorff's alpha from its definition: per unit sum over ordered value pairs d2 / (m - 1), over expected"""
    K = sc[[f"k{k}" for k in range(5)]].values.astype(float)
    m = K.sum(1)
    nc = K.sum(0)
    n = nc.sum()
    D = np.zeros((5, 5))
    for c in range(5):
        for e in range(5):
            lo, hi = min(c, e), max(c, e)
            D[c, e] = 0 if c == e else (nc[lo:hi + 1].sum() - (nc[lo] + nc[hi]) / 2) ** 2
    Do = (np.einsum("uc,cd,ud->u", K, D, K) / (m - 1)).sum() / n
    De = (nc @ D @ nc) / (n * (n - 1))
    return 1 - Do / De


def bh(p):
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    r = np.empty(len(p))
    r[o] = np.minimum(q, 1)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duckdb", default=os.environ.get("MMT_DUCKDB", "duckdb"))
    ap.add_argument("--keep", help="run in this folder and keep it")
    a = ap.parse_args()
    top = a.keep or tempfile.mkdtemp(prefix="landscape_test_")
    shutil.rmtree(top, ignore_errors=True)
    os.makedirs(top)
    out, work, rep = (os.path.join(top, x) for x in ("OUT", "WORK", "REPORT"))
    rng = np.random.default_rng(7)
    samples, truth = mock(out, rng)
    pmap = {s: f"P{s[1:]}" for s in samples if s != "S08"}
    pmap["S10"] = "P09"                                                      # S09 + S10: one person (same chemistry)
    pmap["S07"] = "P06"                                                      # S06 + S07: one person (2 chemistries)
    mp = os.path.join(top, "map.txt")
    write_map(mp, pmap)
    tool = [sys.executable, os.path.join(HERE, "mmt_landscape.py")]
    run = lambda *x: subprocess.run(tool + list(x), capture_output=True, text=True)
    summ = lambda w, *x: run("summarise", "-o", out, "-w", w, "--map", mp, "--duckdb", a.duckdb, "--jobs", "2",
                             "--threads", "2", "--mem-gb", "2", *x)
    p = summ(work)
    print(p.stdout[-1500:], p.stderr[-800:])
    check(p.returncode == 0, "summarise exits 0")
    p2 = summ(os.path.join(top, "WORK_chunked"), "--chunk-rows", "3000")
    check(p2.returncode == 0 and "MGS.1: done" in p2.stdout and "4 chunks" in p2.stdout,
          "summarise with --chunk-rows 3000: MGS.1 in 4 chunks (c1 split in 3 blocks)")
    p = run("plot", "-w", work, "-r", rep, "--min-sites", "2000")
    print(p.stdout[-800:], p.stderr[-1500:])
    check(p.returncode == 0, "plot exits 0")
    if FAILS:
        sys.exit(1)
    part = lambda mgs, x, w=work: pq.read_table(os.path.join(w, "parts", mgs, x + ".parquet")).to_pandas()

    # ---- patients
    pt = pd.read_csv(os.path.join(work, "patients.tsv"), sep="\t")
    patient, rank = dict(zip(pt["sample"], pt["patient"])), dict(zip(pt["sample"], pt["rank"]))
    check(patient["S09"] == patient["S10"] == "P09" and patient["S06"] == patient["S07"] == "P06" and
          patient["S08"] == "S08" and len(set(patient.values())) == 8,
          "patients: S06 + S07 and S09 + S10 one patient each, S08 (not in map) its own, 8 patients")
    check("1 samples not in" in summ(work).stdout, "warning names samples missing from the map")

    # ---- every part holds A / C only ; G/T rows never counted
    for x in ("cube", "frac", "joint", "score", "sites", "power", "npat", "spatial", "alpha", "site_summary", "track",
              "site_classes"):
        bases = set()
        for f in glob.glob(os.path.join(work, "parts", "*", x + ".parquet")):
            bases |= set(pq.read_table(f, columns=["base"]).column("base").to_pylist())
        check(bases <= {"A", "C"} and bases, f"{x}: bases {sorted(bases)} only A / C")

    # ---- exact per MGS x sample x base numbers against the independent re-derivation
    D = {m: rows_of(out, m) for m in ("MGS.1", "MGS.2", "MGS.3")}
    ref = pd.concat([units_ref(d) for d in D.values()], ignore_index=True)
    u = pd.read_csv(os.path.join(rep, "tables", "per_MGS_x_sample.tsv"), sep="\t")
    m = ref.merge(u, on=["MGS", "sample", "base"], suffixes=("_ref", ""))
    check(len(m) == len(ref) == len(u), f"per MGS x sample x base rows {len(u)} = reference {len(ref)}")
    for col in ("sites", "masked", "clean", "tier_low", "tier_mid", "tier_high", "not_graded", "possible", "called_no_frac",
                "cls0", "cls1", "cls2", "cls3", "cls4"):
        bad = (m[col + "_ref"] != m[col]).sum()
        check(bad == 0, f"{col}: exact in all {len(m)} rows" + (f" ({bad} differ)" if bad else ""))
    fr = pd.concat([pq.read_table(f).to_pandas() for f in glob.glob(os.path.join(work, "parts", "*", "frac.parquet"))])
    fr = fr[fr["fbin"] >= 0].assign(span=lambda x: x["ci_span"].clip(upper=3))
    ci = fr.groupby(["MGS", "sample", "base", "span"])["n"].sum().unstack(fill_value=0).add_prefix("got").reset_index()
    m2 = ref.merge(ci, on=["MGS", "sample", "base"], how="left").fillna(0)
    check(all((m2[f"got{k}"] == m2[f"span{k}"]).all() for k in (1, 2, 3)),
          "frac 95 % CI class span exact (independent Wilson formula)")
    ed = D["MGS.1"]
    e95 = ed[(ed["cat"] == "called") & (ed["frac"].astype(np.float32) == np.float32(0.95)) & ed["graded"]]
    check(len(e95) > 0 and (e95["cls"] == 4).all(), f"frac stored as 0.950 -> near-complete ({len(e95)} rows)")

    # ---- per-site patient counts, chance rates, constitutive labels, alpha: exact
    for mgs in ("MGS.1", "MGS.2"):
        sr, _ = sites_ref(D[mgs], rank, patient)
        got = part(mgs, "site_classes")
        key = ["base", "contig", "pos", "strand"]
        mm = sr.merge(got, on=key, suffixes=("_ref", ""), how="outer", indicator=True)
        check((mm["_merge"] == "both").all(), f"{mgs}: same {len(got)} sites with >= 2 patients ({len(sr)} re-derived)")
        cols = ["n", "n_possible", "max_cov"] + [f"k{k}" for k in range(5)]
        check(all((mm[c + "_ref"] == mm[c]).all() for c in cols), f"{mgs}: n, k0..k4, n_possible, max_cov exact per site")
        check(np.allclose(mm["e4_ref"], mm["e4"]) and np.allclose(mm["e0_ref"], mm["e0"]),
              f"{mgs}: chance expectations e4 / e0 exact (per patient x stratum class shares)")
        for b in "AC":
            g = mm[mm["base"] == b]
            t = g[g["n"] >= 3]
            p = binom.sf(t["k4"] - 1, t["n"], np.clip(t["e4_ref"] / t["n"], 0, 1))
            check(np.allclose(bh(p), t["q_meth"], equal_nan=True), f"{mgs} {b}: binomial p + BH q (always methylated) exact")
            al = part(mgs, "alpha")
            ar = alpha_ref(g)
            check(abs(al[al["base"] == b]["alpha"].iloc[0] - ar) < 1e-9, f"{mgs} {b}: alpha {ar:.4f} = definition")

    # ---- chunking changes nothing
    w2 = os.path.join(top, "WORK_chunked")
    same = True
    for mgs in ("MGS.1", "MGS.2", "MGS.3"):
        for x in ("cube", "power", "npat", "spatial", "alpha", "site_classes", "track"):
            a1, a2 = part(mgs, x), part(mgs, x, w2)
            fl = [c for c in a1.columns if a1[c].dtype.kind == "f" and  # nullable integers (e.g. cls) read as float: keys
                  not np.array_equal(a1[c].dropna(), a1[c].dropna().round())]
            k = [c for c in a1.columns if c not in fl]
            a1, a2 = (z.sort_values(k).reset_index(drop=True) for z in (a1, a2))
            same &= (a1.shape == a2.shape and a1[k].equals(a2[k]) and   # floats: summation order differs between chunkings
                     all(np.allclose(a1[c], a2[c], rtol=1e-9, atol=0, equal_nan=True) for c in fl))
    check(same, "chunks of 3000 rows (each 50-kb block alone) = one chunk: every per-site / test table identical "
                "(counts + labels exact, floats to 1e-9)")

    # ---- planted truth found, nothing else
    st = pd.read_csv(os.path.join(rep, "stats_per_MGS.tsv"), sep="\t")
    g = lambda m_, b_, c_: st[(st["MGS"] == m_) & (st["base"] == b_)][c_].iloc[0]
    for b in "AC":
        check(g("MGS.1", b, "alpha_sub") > 0.5 and g("MGS.1", b, "q_alpha") < 0.05,
              f"fixed sites (MGS.1 {b}): alpha {g('MGS.1', b, 'alpha_sub'):.2f}, null {g('MGS.1', b, 'alpha_null_mean'):+.3f}, "
              f"q {g('MGS.1', b, 'q_alpha'):.3f}")
        check(abs(g("MGS.2", b, "alpha_sub") - g("MGS.2", b, "alpha_null_mean")) < 0.03 and g("MGS.2", b, "q_alpha") >= 0.05,
              f"random sites (MGS.2 {b}): alpha {g('MGS.2', b, 'alpha_sub'):+.3f}, null {g('MGS.2', b, 'alpha_null_mean'):+.3f}, "
              f"q {g('MGS.2', b, 'q_alpha'):.2f}")
    check(st[st["MGS"] == "MGS.3"]["alpha_sub"].isna().all(), "one-patient MGS.3: agreement not tested")
    for mgs in ("MGS.1", "MGS.2"):
        sc = part(mgs, "site_classes").merge(truth[mgs], on=["contig", "pos", "strand", "base"])
        am, au = sc[sc["label"] == "always_methylated"], sc[sc["label"] == "always_unmethylated"]
        if mgs == "MGS.1":
            tm = sc[sc["meth"] & (sc["n"] >= 3) & (sc["base"] == "A")]
            check(len(am) > 100 and am["meth"].all(), f"MGS.1: {len(am)} always-methylated sites, all planted methylated")
            check((tm["label"] == "always_methylated").mean() > 0.5,
                  f"MGS.1: {100 * (tm['label'] == 'always_methylated').mean():.0f} % of planted A sites with >= 3 patients found")
            check((~au["meth"]).all(), f"MGS.1: {len(au)} always-unmethylated sites, none planted methylated")
        else:
            check(len(am) == 0, f"MGS.2 (random): {len(am)} always-methylated sites")
    check(g("MGS.1", "A", "log2_dispersion") > 1 and g("MGS.1", "A", "q_spatial") < 0.05,
          f"island (MGS.1 A): log2 dispersion {g('MGS.1', 'A', 'log2_dispersion'):.2f}, q {g('MGS.1', 'A', 'q_spatial'):.1e}")
    for m_ in ("MGS.2", "MGS.3"):
        check(g(m_, "A", "q_spatial") >= 0.05, f"no island ({m_} A): log2 dispersion {g(m_, 'A', 'log2_dispersion'):+.2f}, "
                                               f"q {g(m_, 'A', 'q_spatial'):.2f}")
    # spatial statistic re-derived for MGS.1 A
    _, pr = sites_ref(D["MGS.1"], rank, patient)
    v = pr[(pr["code"] != 9) & (pr["base"] == "A")]
    s_ = v.groupby(["contig", "pos", "strand"]).agg(f=("code", lambda c: (c >= 1).mean()), mx=("coverage", "max")).reset_index()
    s_["k"] = (s_["f"] >= 0.5).astype(int)
    s_["st"] = np.select([s_["mx"] >= 50, s_["mx"] >= 30, s_["mx"] >= 20, s_["mx"] >= 15], [50, 30, 20, 15], MC)
    s_["p"] = s_.groupby("st")["k"].transform("mean")
    w = s_.assign(win=(s_["pos"] - 1) // 1000, v=s_["p"] * (1 - s_["p"])).groupby(["contig", "win"]).agg(
        n=("k", "size"), k=("k", "sum"), e=("p", "sum"), v=("v", "sum"))
    w = w[w["n"] >= 10]
    dd = ((w["k"] - w["e"]) ** 2 / w["v"]).sum()
    check(abs(dd - g("MGS.1", "A", "spatial_d")) < 1e-6 * dd and len(w) == g("MGS.1", "A", "spatial_windows"),
          f"spatial statistic MGS.1 A re-derived ({dd:.1f} on {len(w)} windows)")
    ac = pd.read_csv(os.path.join(rep, "stats_across_MGS.tsv"), sep="\t")
    check(len(ac) == 6 and set(ac["test"]) == {"agreement", "spatial", "constitutive"}, "stats_across_MGS.tsv: 3 tests x A / C")

    # ---- per-MGS values
    mg = pd.read_csv(os.path.join(rep, "tables", "per_MGS.tsv"), sep="\t")
    uu = u[u["usable"]]
    exp = (uu.assign(f=uu["methylated"] / uu["graded"]).groupby(["MGS", "base"])["f"].mean() * 100).rename("e").reset_index()
    m3 = mg.merge(exp, on=["MGS", "base"])
    check(len(m3) > 0 and np.allclose(m3["methylated_graded_pct"], m3["e"]), "per-MGS methylated % = mean over usable samples")
    v2 = mg[(mg["MGS"] == "MGS.2") & (mg["base"] == "A")]["methylated_graded_pct"].iloc[0]
    check(6.5 <= v2 <= 8.5, f"MGS.2 A methylated {v2:.2f} % of graded sites (planted 8 % x detection 0.84-1)")
    lv = mg[(mg["MGS"] == "MGS.1") & (mg["base"] == "A")]["level_near-complete_pct"].iloc[0]
    check(lv > 80, f"MGS.1 A: {lv:.0f} % of methylated sites near-complete (planted frac 1, sd 0.03)")
    s6 = u[(u["MGS"] == "MGS.1") & (u["sample"] == "S06") & (u["base"] == "A")].iloc[0]
    oth = u[(u["MGS"] == "MGS.1") & (u["sample"].isin(["S07", "S08", "S09", "S10"])) & (u["base"] == "A")]
    check(s6["lvl2"] > 0.8 and (oth["lvl4"] > 0.8).all(),
          f"planted half-methylated S06: medium {s6['lvl2']:.2f} of its methylated sites ; others near-complete "
          f">= {oth['lvl4'].min():.2f}")
    s1 = u[(u["sample"] == "S01") & (u["MGS"] == "MGS.1")]
    check((s1["f_not_graded"] > 0.95).all(), f"2x sample S01: {s1['f_not_graded'].min():.2f} not graded")
    n95t = pd.read_csv(os.path.join(rep, "tables", "n95.tsv"), sep="\t")
    n95 = n95t[(n95t["base"] == "A") & (n95t["chemistry"] == "all")]["n95"].iloc[0]
    check(14 <= n95 <= 24, f"n95 diagnostic A = {n95} (planted 1-exp(-c/6) -> 18)")

    # ---- variance (--run): validation, only the variance pass runs, exact pair tables, planted difference
    runf = os.path.join(top, "run.txt")
    write_run = lambda rows: open(runf, "w").write("# sample\tset\n" + "".join(f"{s}\t{x}\n" for s, x in rows))
    for rows, why, msg in (([("S09", "in-batch"), ("S10", "in batch")], "label typo", "want: sample TAB in-batch|across-batch"),
                           ([("S09", "in-batch"), ("S10", "across-batch")], "one sample per set", "is the only"),
                           ([("S09", "in-batch"), ("S99", "in-batch")], "unknown sample", "S99 is not in the map"),
                           ([("S09", "in-batch"), ("S09", "in-batch")], "listed twice", "S09 listed twice")):
        write_run(rows)
        p = summ(work, "--run", runf)
        check(p.returncode != 0 and msg in p.stderr, f"--run {why}: stops before any work ({msg})")
    p = run("summarise", "-o", out, "-w", work, "--run", runf, "--duckdb", a.duckdb)
    check(p.returncode != 0 and "--run needs --map" in p.stderr, "--run without --map: stops")
    write_run([("S09", "in-batch"), ("S10", "in-batch"), ("S06", "across-batch"), ("S07", "across-batch")])
    stamp0 = {f: os.stat(f).st_mtime_ns for f in glob.glob(os.path.join(work, "parts", "*", "cube.parquet"))}
    p = summ(work, "--run", runf)
    print(p.stdout[-600:], p.stderr[-600:])
    check(p.returncode == 0 and "0 to summarise" in p.stdout and "2 variance only" in p.stdout and
          all(os.stat(f).st_mtime_ns == t for f, t in stamp0.items()),
          "--run added: only the variance pass runs (summaries untouched)")
    check(os.path.isfile(os.path.join(work, "parts", "MGS.1", "variance.parquet")) and
          not os.path.exists(os.path.join(work, "parts", "MGS.3", "variance.parquet")),
          "variance tables for MGS.1 / MGS.2 ; none for MGS.3 (one sample)")
    var = pd.concat([part(m_, "variance") for m_ in ("MGS.1", "MGS.2")], ignore_index=True)
    edges = [50, 30, 20, 15]
    rows, bad = [], 0
    for mgs in ("MGS.1", "MGS.2"):
        g_ = D[mgs][D[mgs]["graded"]].copy()
        g_["code"] = g_["cls"].fillna(9).astype(int)
        g_["val"] = np.where(g_["cls"].isna(), np.nan, np.where(g_["cat"] == "none", 0.0, g_["frac"].astype(np.float32).astype(float)))
        for x, s1, s2 in (("in-batch", "S09", "S10"), ("across-batch", "S06", "S07")):
            m_ = g_[g_["sample"] == s1].merge(g_[g_["sample"] == s2], on=["contig", "pos", "strand", "base"])
            lo = np.minimum(m_["coverage_x"], m_["coverage_y"])
            m_["st"] = np.select([lo >= e for e in edges], edges, MC)
            m_["d"] = m_["val_x"] - m_["val_y"]
            r_ = m_.groupby(["base", "st", "code_x", "code_y"]).agg(n=("d", "size"), d_sum=("d", "sum"),
                                                                     d2_sum=("d", lambda z: (z * z).sum())).reset_index()
            got = var[(var["MGS"] == mgs) & (var["set"] == x)].merge(
                r_, left_on=["base", "st", "code_a", "code_b"], right_on=["base", "st", "code_x", "code_y"], how="outer")
            bad += int(got["n_x"].ne(got["n_y"]).sum()) + int((~np.isclose(got["d2_sum_x"].fillna(0), got["d2_sum_y"])).sum())
            rows.append(m_.assign(set=x, MGS=mgs))
    check(bad == 0, "pair tables exact: every site x class x class x coverage stratum count and frac sum (pandas re-derivation)")
    rp = os.path.join(top, "REPORT_var")
    p = run("plot", "-w", work, "-r", rp, "--min-sites", "2000", "--cards", "0")
    check(p.returncode == 0 and all(os.path.isfile(os.path.join(rp, f"variance_{x}.tsv")) for x in ("summary", "by_coverage", "pairs")),
          "plot: variance_summary / _by_coverage / _pairs tables")
    vs = pd.read_csv(os.path.join(rp, "variance_summary.tsv"), sep="\t")
    check(list(zip(vs["set"], vs["base"])) == [("in-batch", "A"), ("in-batch", "C"), ("across-batch", "A"), ("across-batch", "C")],
          "variance_summary: the two result sets x A / C")
    allm = pd.concat(rows)
    for x in ("in-batch", "across-batch"):
        z = allm[(allm["set"] == x) & (allm["base"] == "A") & (allm["code_x"] < 9) & (allm["code_y"] < 9)]
        K = np.zeros((len(z), 5))
        K[np.arange(len(z)), z["code_x"]] += 1
        K[np.arange(len(z)), z["code_y"]] += 1
        r = vs[(vs["set"] == x) & (vs["base"] == "A")].iloc[0]
        check(abs(alpha_ref(pd.DataFrame(K, columns=[f"k{k}" for k in range(5)])) - r["alpha"]) < 1e-9 and r["sites"] == len(z)
              and abs(r["same_class_pct"] - 100 * (z["code_x"] == z["code_y"]).mean()) < 1e-9
              and abs(r["frac_sd_per_sample"] - np.sqrt((z["d"] ** 2).mean() / 2)) < 1e-9,
              f"{x} A pooled: alpha {r['alpha']:.4f}, same class {r['same_class_pct']:.2f} %, frac sd {r['frac_sd_per_sample']:.4f} "
              f"on {r['sites']:,} sites = re-derived")
    va = vs.set_index(["set", "base"])
    check(va.loc[("across-batch", "A"), "frac_sd_per_sample"] > 1.2 * va.loc[("in-batch", "A"), "frac_sd_per_sample"] and
          va.loc[("across-batch", "A"), "alpha"] < va.loc[("in-batch", "A"), "alpha"],
          f"planted: across-batch pair (S06 half-methylated) noisier than in-batch: frac sd "
          f"{va.loc[('across-batch', 'A'), 'frac_sd_per_sample']:.3f} vs {va.loc[('in-batch', 'A'), 'frac_sd_per_sample']:.3f}")
    # ---- --between (fig7): pair choice, exact switched_pct, planted biological change in MGS.2 only
    p = run("summarise", "-o", out, "-w", work, "--between", "3", "--duckdb", a.duckdb)
    check(p.returncode != 0 and "--between needs --map" in p.stderr, "--between without --map: stops")
    p = summ(work, "--run", runf, "--between", "3")
    check(p.returncode == 0 and "0 to summarise" in p.stdout and "2 variance only" in p.stdout and
          all(os.stat(f).st_mtime_ns == t for f, t in stamp0.items()), "--between added: only the variance pass runs")
    ut = pd.read_csv(os.path.join(work, "units.tsv"), sep="\t")
    for mgs in ("MGS.1", "MGS.2"):
        bp = part(mgs, "variance").query("set == 'between-patient'")[["sample_a", "sample_b"]].drop_duplicates()
        u_ = ut[ut["MGS"] == mgs].assign(patient=lambda d: d["sample"].map(patient)).sort_values(["rows", "sample"],
                                                                                                ascending=[False, True])
        want = set(u_.drop_duplicates("patient")["sample"].head(6))
        got = list(bp["sample_a"]) + list(bp["sample_b"])
        check(len(bp) == 3 and set(got) == want and len(set(map(patient.get, got))) == 6,
              f"{mgs} --between 3: 3 disjoint pairs of 6 different patients, best-covered sample each ({sorted(got)})")
    rb = os.path.join(top, "REPORT_between")
    p = run("plot", "-w", work, "-r", rb, "--min-sites", "2000", "--cards", "0")
    print(p.stdout[-800:], p.stderr[-800:])
    check(p.returncode == 0 and all(os.path.isfile(os.path.join(rb, f)) for f in
          ("fig7_variation.pdf", "fig7_variation.png", "variance_evidence.tsv", "variance_bands.tsv")),
          "plot: fig7_variation + variance_evidence / variance_bands")
    vp = pd.read_csv(os.path.join(rb, "variance_pairs.tsv"), sep="\t")
    bad = 0
    for mgs in ("MGS.1", "MGS.2"):
        g_ = D[mgs][D[mgs]["graded"]].copy()
        g_["code"] = g_["cls"].fillna(9).astype(int)
        for _, r in vp[vp["MGS"] == mgs].iterrows():
            m_ = g_[(g_["sample"] == r["sample_a"]) & (g_["base"] == r["base"])].merge(
                g_[g_["sample"] == r["sample_b"]], on=["contig", "pos", "strand", "base"])
            m_ = m_[(np.minimum(m_["coverage_x"], m_["coverage_y"]) >= 30) & (m_["code_x"] < 9) & (m_["code_y"] < 9)]
            un = int(((m_["code_x"] > 0) | (m_["code_y"] > 0)).sum())
            sw = int((((m_["code_x"] == 0) & (m_["code_y"] > 0)) | ((m_["code_y"] == 0) & (m_["code_x"] > 0))).sum())
            bad += (un != r["methylated_either_hicov"]) or (un > 0 and abs(100 * sw / un - r["switched_pct_hicov"]) > 1e-9)
    check(bad == 0, f"switched_pct at >= 30 reads exact for all {len(vp)} pair x MGS x base rows (pandas re-derivation)")
    ev = pd.read_csv(os.path.join(rb, "variance_evidence.tsv"), sep="\t").set_index(["MGS", "base"])
    e2, e1 = ev.loc[("MGS.2", "A")], ev.loc[("MGS.1", "A")]
    tech2 = np.nanmax([e2["in_batch_median"], e2["across_batch_median"]])
    check(e2["between_patient_median"] > 3 * tech2 and bool(e2["above_band"]) and not bool(e1["above_band"]),
          f"planted: MGS.2 between patients {e2['between_patient_median']:.1f} % switched vs replicates {tech2:.1f} % "
          f"(above band) ; MGS.1 {e1['between_patient_median']:.1f} % (within band)")
    bd = pd.read_csv(os.path.join(rb, "variance_bands.tsv"), sep="\t").set_index(["set", "base"])
    nb = vp[vp["methylated_either_hicov"] >= 30].groupby(["set", "base"]).size()
    check(all(bd.loc[k, "pairs"] == nb.get(k, 0) for k in bd.index),
          "variance_bands: pair counts = pairs x MGS with >= 30 methylated sites at >= 30 reads")
    p = run("plot", "-w", work, "-r", rb, "--cards", "0", "--pair-min-cov", "25")
    check(p.returncode != 0 and "stratum edge" in p.stderr, "--pair-min-cov not a stratum edge: stops")
    write_run([("S09", "in-batch"), ("S10", "in-batch")])
    vst =os.stat(os.path.join(work, "parts", "MGS.1", "variance.parquet")).st_mtime_ns
    p = summ(work, "--run", runf)
    check(p.returncode == 0 and "0 to summarise" in p.stdout and vst != os.stat(os.path.join(work, "parts", "MGS.1",
          "variance.parquet")).st_mtime_ns and all(os.stat(f).st_mtime_ns == t for f, t in stamp0.items()),
          "run file edited: only the variance pass is redone")
    p = summ(work, "--run", runf)
    check(p.returncode == 0 and "0 to summarise (0.0 GB), 0 variance only" in p.stdout, "same run file again: nothing redone")
    p = summ(work)
    check(p.returncode == 0 and not glob.glob(os.path.join(work, "parts", "*", "variance.*")) and
          not os.path.exists(os.path.join(work, "variance_pairs.tsv")), "summarise without --run: variance tables removed")

    # ---- outputs
    for f in ("fig1_data_quality", "fig2_coverage_effects", "fig3_landscape", "fig4_levels", "fig5_confidence",
              "fig6_non_random"):
        check(os.path.getsize(os.path.join(rep, f + ".pdf")) > 5000 and os.path.isfile(os.path.join(rep, f + ".png")), f"{f} pdf + png")
    import re
    pages = len(re.findall(rb"/Type\s*/Page[^s]", open(os.path.join(rep, "mgs_cards.pdf"), "rb").read()))
    want = mg.groupby("MGS")["samples_usable"].max().ge(1).sum()
    check(pages == want, f"cards: {pages} pages = MGS with a usable sample ({want})")
    cs = pq.read_table(os.path.join(rep, "constitutive_sites.parquet")).to_pandas()
    check(set(cs["label"]) <= {"always_methylated", "always_unmethylated"} and len(cs) > 100, f"constitutive_sites: {len(cs)} rows")

    # ---- resume, change, map change, drop, corrupt, out of memory
    stamp = {f: os.stat(f).st_mtime_ns for f in glob.glob(os.path.join(work, "parts", "*", "cube.parquet"))}
    p = summ(work)
    check(p.returncode == 0 and "0 to summarise" in p.stdout and
          all(os.stat(f).st_mtime_ns == t for f, t in stamp.items()), "rerun: nothing redone")
    f2 = os.path.join(out, "mgs-rep", "MGS.2", "S05.parquet")
    t = pq.read_table(f2)
    pq.write_table(t.slice(0, t.num_rows // 2), f2)
    p = summ(work)
    redo = [f for f, t in stamp.items() if os.stat(f).st_mtime_ns != t]
    check(p.returncode == 0 and len(redo) == 1 and "/MGS.2/" in redo[0], f"changed file: only MGS.2 redone ({len(redo)})")
    pmap["S08"] = "P08"
    write_map(mp, pmap)
    p = summ(work)
    check(p.returncode == 0 and "3 to summarise" in p.stdout, "map changed (patients): every MGS redone")
    shutil.rmtree(os.path.join(out, "mgs-rep", "MGS.3"))
    p = summ(work)
    check(p.returncode == 0 and not os.path.exists(os.path.join(work, "parts", "MGS.3")), "MGS gone from OUT -> dropped")
    fake = os.path.join(top, "duckdb_oom")
    with open(fake, "w") as o:
        o.write("#!/bin/sh\nkill -9 $$\n")
    os.chmod(fake, os.stat(fake).st_mode | stat.S_IEXEC)
    p = run("summarise", "-o", out, "-w", os.path.join(top, "WORK_oom"), "--map", mp, "--duckdb", fake)
    check(p.returncode != 0 and "killed by signal 9" in p.stdout and "out of memory" in p.stdout,
          "DuckDB killed (out of memory): reported as signal 9 with advice")
    bad = os.path.join(out, "mgs-rep", "MGS.2", "S06.parquet")
    open(bad, "wb").write(b"not parquet")
    p = summ(work)
    check(p.returncode != 0 and "MGS.2" in p.stderr and os.path.isfile(os.path.join(work, "parts", "MGS.1", "fp")),
          "corrupt file: exits non-zero naming MGS.2, other MGS kept")

    print(f"\n{'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILED'}  ({top})")
    if not a.keep and not FAILS:
        shutil.rmtree(top, ignore_errors=True)
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
