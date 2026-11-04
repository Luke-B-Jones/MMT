#!/usr/bin/env python3
"""Self-check of mmt_landscape.py on a mock MMT output written with MMT's own schema (lib/mmtpq.py).

  python3 analysis/test_landscape.py [--duckdb BIN] [--keep DIR]

Planted truth: 3 MGS, 10 samples (2 chemistries), coverage from 2x to 90x, 8 % of A sites methylated (m6A) and 3 % of
C sites (m4C) in every sample, sample S03 half-methylated on MGS.1, detection 1 - exp(-coverage * frac / 6) (-> n95
near 18), masked rows, G/T modified_base rows (must never be counted). Every summary number is re-derived here with
pandas from the parquet files; the run must match exactly. Then resume, change, drop and the figures are checked.
"""
import argparse, glob, os, shutil, subprocess, sys, tempfile
import numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import mmtpq  # noqa: E402

COMP = str.maketrans("ACGT", "TGCA")
FAILS = []


def check(ok, what):
    print(("PASS " if ok else "FAIL ") + what)
    if not ok:
        FAILS.append(what)


def mock(out, rng):
    samples = [f"S{i:02d}" for i in range(1, 11)]
    lam = dict(zip(samples, [2, 4, 8, 12, 18, 25, 35, 50, 70, 90]))          # mean coverage per strand
    chem = {s: ("102-739-100/102-118-800" if i < 6 else "103-496-900/103-496-700") for i, s in enumerate(samples)}
    genomes = {"MGS.1": {"c1": 30000, "c2": 20000}, "MGS.2": {"k1": 30000}, "MGS.3": {"z1": 8000}}
    present = {"MGS.1": samples, "MGS.2": samples[2:], "MGS.3": ["S09"]}
    stats = {s: [] for s in samples}
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
        os.makedirs(os.path.join(out, "mgs-rep", mgs), exist_ok=True)
        for s in samples:
            if s not in present[mgs]:
                continue
            d = S.copy()
            if mgs == "MGS.2":   # null for site consistency: methylated sites drawn anew in every sample
                v = rng.random(len(d))
                d["meth"] = ((d["base"] == "A") & (v < 0.08)) | ((d["base"] == "C") & (v < 0.03))
            d["coverage"] = rng.poisson(lam[s] * d["site_eff"])
            d = d[d["coverage"] >= 3].copy()
            tf = np.where(d["meth"], 0.5 if (s == "S03" and mgs == "MGS.1") else 1.0, 0.0)
            pcall = np.where(d["meth"], 1 - np.exp(-d["coverage"] * tf / 6), 0.002)
            r = rng.random(len(d))
            called = r < pcall
            possible = (~called) & (rng.random(len(d)) < np.where(d["meth"], 0.3, 0.01))
            d["call"] = np.where(called, np.where(d["base"] == "A", "m6A", np.where(d["base"] == "C", "m4C", "none")),
                                 np.where(possible, "modified_base", "none"))
            gt = d["base"].isin(["G", "T"])
            d.loc[gt, "call"] = np.where(rng.random(gt.sum()) < 0.005, "modified_base", "none")
            d = d[~gt | (d["call"] != "none")]                                # MMT keeps G/T only with a call
            isc = d["call"].isin(["m6A", "m4C"])
            truth = np.where(d["meth"], 0.5 if (s == "S03" and mgs == "MGS.1") else 1.0, 0.0)
            d["frac"] = np.where(isc & (d["coverage"] > 10), np.round(np.clip(truth + rng.normal(0, 0.06, len(d)), 0, 1), 3), np.nan)
            d["score"] = np.where(isc, 20 + 2 * d["coverage"] * truth, rng.integers(0, 15, len(d))).astype(int)
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
    return samples


def reference(out, n95):
    """independent pandas re-derivation of the per sample x MGS x base numbers"""
    rows = []
    for f in glob.glob(os.path.join(out, "mgs-rep", "MGS.*", "*.parquet")):
        d = pq.read_table(f).to_pandas()
        d = d[d["base"].isin(["A", "C"]) & (d["coverage"] >= 3)]
        clean = (d["pile_depth"] > 0) & (d["n_ref"] >= 0.95 * d["pile_depth"]) & (d["ctx_alt_pct"] < 10) & (d["clip_ctx_n"] == 0)
        for b in "AC":
            x = d[(d["base"] == b)]
            c = x[clean[x.index]]
            want = "m6A" if b == "A" else "m4C"
            called = c["call"] == want
            hi = c["coverage"] >= n95[b]
            fr = c[called & c["frac"].notna() & hi]
            z = 1.96
            n, p = c["coverage"], c["frac"]
            cen = (p + z * z / (2 * n)) / (1 + z * z / n)
            h = z * np.sqrt(np.clip(p * (1 - p), 0, None) / n + z * z / (4 * n * n)) / (1 + z * z / n)
            cic = c[called & c["frac"].notna()]
            up, lo = (cen + h)[cic.index], (cen - h)[cic.index]
            rows.append(dict(MGS=x["MGS"].iloc[0] if len(x) else "", sample=os.path.basename(f)[:-8], base=b,
                             sites=len(x), masked=int((~clean[x.index]).sum()), clean=len(c), called=int(called.sum()),
                             possible=int((c["call"].isin(["modified_base"]) | (~c["call"].isin(["none", want, "modified_base"]))).sum()),
                             unmethylated=int(((c["call"] == "none") & hi).sum()),
                             undetermined=int(((c["call"] == "none") & ~hi).sum()),
                             full=int((fr["frac"] >= 0.9).sum()), high_partial=int(((fr["frac"] >= 0.5) & (fr["frac"] < 0.9)).sum()),
                             low_partial=int((fr["frac"] < 0.5).sum()),
                             ci_partial=int((up < 0.9).sum()), ci_full=int(((up >= 0.9) & (lo >= 0.5)).sum())))
    return pd.DataFrame(rows)


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
    mock(out, rng)
    tool = [sys.executable, os.path.join(HERE, "mmt_landscape.py")]
    run = lambda *x: subprocess.run(tool + list(x), capture_output=True, text=True)
    p = run("summarise", "-o", out, "-w", work, "--duckdb", a.duckdb, "--jobs", "2", "--threads", "2", "--mem-gb", "2")
    print(p.stdout[-800:], p.stderr[-800:])
    check(p.returncode == 0, "summarise exits 0")
    p = run("plot", "-w", work, "-r", rep, "--min-sites", "2000")
    print(p.stdout[-800:], p.stderr[-1500:])
    check(p.returncode == 0, "plot exits 0")
    if FAILS:
        sys.exit(1)

    # ---- every part holds A / C only ; G/T rows never counted
    for x in ("cube", "frac", "joint", "score", "sites", "power", "consistency", "spatial"):
        bases = set()
        for f in glob.glob(os.path.join(work, "parts", "*", x + ".parquet")):
            bases |= set(pq.read_table(f, columns=["base"]).column("base").to_pylist())
        check(bases <= {"A", "C"}, f"{x}: bases {sorted(bases)} only A / C")

    n95t = pd.read_csv(os.path.join(rep, "tables", "n95.tsv"), sep="\t")
    n95 = {b: int(n95t[(n95t["base"] == b) & (n95t["chemistry"] == "all")]["n95_used"].iloc[0]) for b in "AC"}
    check(14 <= n95["A"] <= 24, f"n95 A = {n95['A']} (planted 1-exp(-c/6) -> 18)")
    check(14 <= n95["C"] <= 32, f"n95 C = {n95['C']} (planted 1-exp(-c/6) -> 18; few C references in the mock)")

    # ---- exact counts against the independent re-derivation
    ref = reference(out, n95)
    u = pd.read_csv(os.path.join(rep, "tables", "per_MGS_x_sample.tsv"), sep="\t")
    m = ref.merge(u, on=["MGS", "sample", "base"], suffixes=("_ref", ""))
    check(len(m) == len(ref) == len(u), f"per MGS x sample x base rows {len(u)} = reference {len(ref)}")
    for col in ("sites", "masked", "clean", "called", "possible", "unmethylated", "undetermined",
                "full", "high_partial", "low_partial"):
        bad = (m[col + "_ref"] != m[col]).sum()
        check(bad == 0, f"{col}: exact in all {len(m)} rows" + (f" ({bad} differ)" if bad else ""))
    fr = pd.concat([pq.read_table(f).to_pandas() for f in glob.glob(os.path.join(work, "parts", "*", "frac.parquet"))])
    ci = fr.groupby(["MGS", "sample", "base", "ci"])["n"].sum().unstack(fill_value=0).add_prefix("got_").reset_index()
    m2 = ref.merge(ci, on=["MGS", "sample", "base"])
    check(len(m2) == len(ref) and (m2["got_partial"] == m2["ci_partial"]).all() and (m2["got_full"] == m2["ci_full"]).all(),
          "Wilson CI classes exact (independent formula)")

    # ---- per-MGS means = mean over usable samples
    mg = pd.read_csv(os.path.join(rep, "tables", "per_MGS.tsv"), sep="\t")
    uu = u[u["usable"]]
    exp = (uu.assign(f=uu["called"] / uu["clean"]).groupby(["MGS", "base"])["f"].mean() * 100).rename("e").reset_index()
    m3 = mg.merge(exp, on=["MGS", "base"])
    check(np.allclose(m3["called_pct"], m3["e"]), "per-MGS called % = mean over usable samples")
    s3 = u[(u["MGS"] == "MGS.1") & (u["sample"] == "S03") & (u["base"] == "A")].iloc[0]
    others = u[(u["MGS"] == "MGS.1") & (u["sample"] != "S03") & (u["base"] == "A") & (u["frac_n"] > 50)]
    check(s3["f_full"] < 0.2 < others["f_full"].min(), f"planted half-methylated S03 seen: full {s3['f_full']:.2f} vs others "
                                                         f">= {others['f_full'].min():.2f}")
    a_called = mg[mg["base"] == "A"]["called_or_possible_pct"]
    check(a_called.notna().all(), "per-MGS A values present")
    s1 = u[(u["sample"] == "S01") & (u["MGS"] == "MGS.1")]
    check((s1["f_undetermined"] > 0.8).all(), f"2x sample S01: mostly undetermined ({s1['f_undetermined'].min():.2f})")
    hi = mg[(mg["MGS"] == "MGS.2") & (mg["base"] == "A")]["called_hi_pct"].iloc[0]
    check(7.0 <= hi <= 8.5, f"coverage-controlled called % MGS.2 A = {hi:.2f} (planted 8 % x detection >= 0.95)")
    sl = pd.read_csv(os.path.join(rep, "tables", "fig2_coverage_effects_A_skew.tsv"), sep="\t")
    check(len(sl) > 0 and sl["slope"].median() > 0, f"planted coverage effect seen as positive skew (median {sl['slope'].median():+.3f})")
    pw = pq.read_table(os.path.join(work, "parts", "MGS.1", "power.parquet")).to_pandas()   # fixed sites: true negatives
    un = pw[(pw["kind"] == "unmethylated") & (pw["base"] == "A")]
    check(un["n"].sum() > 1000 and un["n_called"].sum() / un["n"].sum() < 0.01,
          f"false-call rate at unmethylated references {un['n_called'].sum() / max(un['n'].sum(), 1):.4f} (planted 0.002)")

    # ---- non-random methylation tests
    st = pd.read_csv(os.path.join(rep, "stats_per_MGS.tsv"), sep="\t")
    g = lambda m_, b_, c_: st[(st["MGS"] == m_) & (st["base"] == b_)][c_].iloc[0]
    check(g("MGS.1", "A", "rho") > 0.5 and g("MGS.1", "A", "q_consistency") < 0.05,
          f"fixed sites (MGS.1 A): rho {g('MGS.1', 'A', 'rho'):.2f}, q {g('MGS.1', 'A', 'q_consistency'):.1e}")
    check(abs(g("MGS.2", "A", "rho")) < 0.05 and g("MGS.2", "A", "q_consistency") >= 0.05,
          f"random sites (MGS.2 A): rho {g('MGS.2', 'A', 'rho'):.3f}, q {g('MGS.2', 'A', 'q_consistency'):.2f}")
    check(g("MGS.1", "A", "log2_dispersion") > 1 and g("MGS.1", "A", "q_spatial") < 0.05,
          f"island (MGS.1 A): log2 dispersion {g('MGS.1', 'A', 'log2_dispersion'):.2f}, q {g('MGS.1', 'A', 'q_spatial'):.1e}")
    check(g("MGS.3", "A", "log2_dispersion") < 0.3 and g("MGS.3", "A", "q_spatial") >= 0.05,
          f"fixed sites, no island (MGS.3 A): log2 dispersion {g('MGS.3', 'A', 'log2_dispersion'):.2f}, "
          f"q {g('MGS.3', 'A', 'q_spatial'):.2f}")
    check(st[st["MGS"] == "MGS.3"]["rho"].isna().all(), "single-sample MGS.3: consistency not tested")
    # rho re-derived independently for MGS.1 A
    rows = []
    for f in glob.glob(os.path.join(out, "mgs-rep", "MGS.1", "*.parquet")):
        d = pq.read_table(f).to_pandas()
        d = d[(d["base"] == "A") & (d["coverage"] >= 30) & (d["pile_depth"] > 0) & (d["n_ref"] >= 0.95 * d["pile_depth"])
              & (d["ctx_alt_pct"] < 10) & (d["clip_ctx_n"] == 0)]
        rows.append(d.assign(c=(d["call"] == "m6A").astype(int))[["sample", "contig", "pos", "strand", "c"]])
    d = pd.concat(rows)
    d["p"] = d.groupby("sample")["c"].transform("mean")
    d["v"], d["a"] = d["p"] * (1 - d["p"]), np.sqrt(d["p"] * (1 - d["p"]))
    s_ = d.groupby(["contig", "pos", "strand"]).agg(n=("c", "size"), k=("c", "sum"), e=("p", "sum"), v=("v", "sum"), a=("a", "sum"))
    s_ = s_[s_["n"] >= 2]
    rho = ((s_["k"] - s_["e"]) ** 2 - s_["v"]).sum() / (s_["a"] ** 2 - s_["v"]).sum()
    check(abs(rho - g("MGS.1", "A", "rho")) < 1e-9 and len(s_) == g("MGS.1", "A", "consistency_sites"),
          f"rho MGS.1 A re-derived {rho:.6f} on {len(s_)} sites")
    ac = pd.read_csv(os.path.join(rep, "stats_across_MGS.tsv"), sep="\t")
    check(len(ac) == 4 and set(ac["test"]) == {"consistency", "spatial"}, "stats_across_MGS.tsv: 2 tests x A / C")

    # ---- outputs
    for f in ("fig1_data_quality", "fig2_coverage_effects", "fig3_landscape", "fig4_completeness", "fig5_confidence",
              "fig6_non_random"):
        check(os.path.getsize(os.path.join(rep, f + ".pdf")) > 5000 and os.path.isfile(os.path.join(rep, f + ".png")), f"{f} pdf + png")
    import re
    pages = len(re.findall(rb"/Type\s*/Page[^s]", open(os.path.join(rep, "mgs_cards.pdf"), "rb").read()))
    want = mg.groupby("MGS")["samples_usable"].max().ge(1).sum()
    check(pages == want, f"cards: {pages} pages = MGS with a usable sample ({want})")

    # ---- resume, change, drop
    stamp = {f: os.stat(f).st_mtime_ns for f in glob.glob(os.path.join(work, "parts", "*", "cube.parquet"))}
    p = run("summarise", "-o", out, "-w", work, "--duckdb", a.duckdb)
    check(p.returncode == 0 and "0 to summarise" in p.stdout and
          all(os.stat(f).st_mtime_ns == t for f, t in stamp.items()), "rerun: nothing redone")
    f2 = os.path.join(out, "mgs-rep", "MGS.2", "S05.parquet")
    t = pq.read_table(f2)
    pq.write_table(t.slice(0, t.num_rows // 2), f2)
    p = run("summarise", "-o", out, "-w", work, "--duckdb", a.duckdb)
    redo = [f for f, t in stamp.items() if os.stat(f).st_mtime_ns != t]
    check(p.returncode == 0 and len(redo) == 1 and "/MGS.2/" in redo[0], f"changed file: only MGS.2 redone ({len(redo)})")
    shutil.rmtree(os.path.join(out, "mgs-rep", "MGS.3"))
    p = run("summarise", "-o", out, "-w", work, "--duckdb", a.duckdb)
    check(p.returncode == 0 and not os.path.exists(os.path.join(work, "parts", "MGS.3")), "MGS gone from OUT -> dropped")
    bad = os.path.join(out, "mgs-rep", "MGS.2", "S06.parquet")
    open(bad, "wb").write(b"not parquet")
    p = run("summarise", "-o", out, "-w", work, "--duckdb", a.duckdb)
    check(p.returncode != 0 and "MGS.2" in p.stderr and os.path.isfile(os.path.join(work, "parts", "MGS.1", "fp")),
          "corrupt file: exits non-zero naming MGS.2, other MGS kept")

    print(f"\n{'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILED'}  ({top})")
    if not a.keep and not FAILS:
        shutil.rmtree(top, ignore_errors=True)
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
