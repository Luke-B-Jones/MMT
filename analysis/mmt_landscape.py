#!/usr/bin/env python3
"""MMT landscape: cohort-wide methylation statistics + publication figures from MMT MGScall --mgs-rep results.

  summarise  one pass per MGS over <OUT>/mgs-rep/<MGS>/<sample>.parquet with the DuckDB CLI (the binary MMT uses)
             -> small tables in WORK/parts/<MGS>/ + WORK/units.tsv, samples.tsv (from <OUT>/mgs-rep/stats/) and, with
             --genecat, genome.tsv (A / C sites of each representative). Resumes: an MGS whose files did not change is
             skipped ; MGS no longer in OUT are dropped from WORK.
  plot       WORK -> REPORT: fig1..fig5 (pdf + png), mgs_cards.pdf (one page per MGS), tables/*.tsv (every plotted number)

  mmt_landscape.py summarise -o OUT -w WORK [--duckdb BIN] [--jobs 4] [--threads 4] [--mem-gb 16] [--genecat GC --binner SC]
  mmt_landscape.py plot -w WORK -r REPORT [--taxonomy TSV] [--n95 A=20,C=25] [--min-sites 20000] [--min-breadth 0.1]

Every statistic is computed for A sites (-> m6A) and C sites (-> m4C) separately. G/T rows are never read.
  site          one base on one strand of the representative with an MMT row (>= 3 kinetic reads, MMT's minimum)
  clean         >= --min-ref-pct % of the reads carry the reference base, no position within MMT's context window where
                >= --max-ctx-pct % of the reads differ (ctx_alt_pct), no read clip in that window (clip_ctx_n = 0)
  state         clean sites only: called (m6A on A, m4C on C) | possible (modified_base, or the other type) |
                unmethylated (no call, coverage >= n95) | undetermined (no call, coverage < n95)
  n95           lowest coverage from which >= 95 % of reference-methylated sites are called. Reference sites: called in
                >= 90 % of >= 3 samples with >= 30 reads, chosen in one half of the samples (sorted names, odd / even)
                and measured in the other half. Reference-unmethylated sites (never any call there; a 5 % subset) give
                the false-call rate.
  usable        an MGS x sample with >= --min-sites clean A+C sites covering >= --min-breadth of the representative's
                A+C sites (--genecat), else of the union of sites seen in any sample
  frac classes  called sites with coverage >= n95: full (frac >= 0.9) | high partial (0.5-0.9) | low partial (< 0.5)
  consistency   are the same bases methylated across samples? Sites with >= --sig-cov reads in >= 2 samples; k = samples
                calling the site, expected sum of p_s (p_s = the sample's call rate at such sites). rho = sum((k - E)^2 - V)
                / sum(cross terms): 0 = calls placed at random, 1 = the same sites in every sample. p: chi-square
                overdispersion sum((k - E)^2 / V) on the tested sites ; q: BH across MGS. Not tested (NaN) with < 100
                such sites, or when only one sample makes calls there.
  spatial       are methylated sites spread at random along the genome, or clustered in regions? One map per MGS: a
                site with >= --sig-cov reads is methylated when >= half of the samples covering it call it (samples share
                the genome, so per-sample tests are not independent and are not pooled). --window bp windows with >= 10
                such sites: dispersion sum((k - n p)^2 / n p (1 - p)) / (windows - 1) ; log2 0 = random, > 0 = clustered
                (methylation islands, or regions unusually rich / poor in the target motifs), < 0 = evenly spaced.
                p: chi-square ; q: BH across MGS.
  frac CI       Wilson 95 % interval of frac with n = coverage: confident partial (upper < 0.9) | consistent with full
                (upper >= 0.9, lower >= 0.5) | uncertain (the rest). ipdSummary's fracLow / fracUp are not used: they
                come from an unseeded bootstrap and change between runs of the same data.
"""
import argparse, concurrent.futures as cf, glob, gzip, hashlib, json, math, os, shutil, subprocess, sys, time

VERSION = "1"
COV_EDGES = list(range(3, 31)) + [31, 35, 40, 50, 60, 80, 100, 150, 200, 300, 500, 1000]   # fine coverage bins (lower edges)
COARSE = [3, 11, 15, 20, 30, 50, 100]                                                       # coarse bins (lower edges)
PARTS = ("cube", "frac", "joint", "score", "sites", "power", "consistency", "spatial")
MOD = {"A": "m6A", "C": "m4C"}


def die(msg):
    sys.exit("ERROR: " + msg)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def case_bins(col, edges):
    """SQL: lower edge of the bin holding col (edges ascending ; below the first edge -> 0)"""
    return "CASE " + " ".join(f"WHEN {col} >= {e} THEN {e}" for e in reversed(edges)) + " ELSE 0 END"


def read_kv(path):
    """MMT run.tsv (key TAB value lines; first value of a repeated key wins)"""
    kv = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t", 1)
            if len(f) == 2 and f[0] not in kv:
                kv[f[0]] = f[1]
    return kv


# ============================================================================ summarise
def mgs_sql(mgs, src, dst, a):
    """one MGS: every part as a COPY from a table of clean / state / bins per A or C row"""
    q = lambda s: s.replace("'", "''")
    clean = (f"pile_depth > 0 AND n_ref >= {a.min_ref_pct / 100.0} * pile_depth AND ctx_alt_pct < {a.max_ctx_pct} "
             f"AND clip_ctx_n = 0")
    out = lambda name: f"'{q(os.path.join(dst, name + '.parquet'))}' (FORMAT parquet, COMPRESSION zstd)"
    m = q(mgs)
    # Wilson 95 % bound of frac with n = coverage: (p + z2/2n +- z sqrt(p(1-p)/n + z2/4n2)) / (1 + z2/n)
    wilson = lambda sg: (f"((frac + 1.9208 / coverage {sg} 1.96 * sqrt(greatest(frac * (1 - frac), 0) / coverage "
                         f"+ 0.9604 / (coverage * coverage))) / (1 + 3.8416 / coverage))")
    return f"""
SET threads = {a.threads}; SET memory_limit = '{a.mem_gb}GB'; SET preserve_insertion_order = false;
SET temp_directory = '{q(os.path.join(dst, 'duck.tmp'))}';
CREATE TABLE t AS
SELECT sample, contig, pos, strand, base, coverage, score, frac,
       (CASE WHEN {clean} THEN 1 ELSE 0 END)::UTINYINT AS clean,
       CASE WHEN call = (CASE base WHEN 'A' THEN 'm6A' ELSE 'm4C' END) THEN 'called'
            WHEN call = 'none' THEN 'none' WHEN call = 'modified_base' THEN 'possible' ELSE 'other' END AS cat,
       ({case_bins('coverage', COV_EDGES)})::USMALLINT AS cb,
       ({case_bins('coverage', COARSE)})::USMALLINT AS cc
FROM read_parquet('{q(src)}') WHERE base IN ('A', 'C') AND coverage >= 3;
COPY (SELECT '{m}' AS MGS, sample, base, cb, clean, cat, count(*)::BIGINT AS n, sum(coverage)::BIGINT AS cov_sum
      FROM t GROUP BY ALL) TO {out('cube')};
CREATE TABLE fr AS
SELECT sample, base, cb, cc, coverage, score,
       CASE WHEN frac IS NULL THEN -1 ELSE least(floor(frac * 20), 19)::INTEGER END AS fbin,
       CASE WHEN frac IS NULL THEN 'no_frac'
            WHEN {wilson('+')} < 0.9 THEN 'partial' WHEN {wilson('-')} >= 0.5 THEN 'full' ELSE 'uncertain' END AS ci
FROM t WHERE clean = 1 AND cat = 'called';
COPY (SELECT '{m}' AS MGS, sample, base, cb, fbin, ci, count(*)::BIGINT AS n FROM fr GROUP BY ALL) TO {out('frac')};
COPY (SELECT '{m}' AS MGS, base, cc, fbin, least(floor(score / 5), 50)::INTEGER * 5 AS sbin, count(*)::BIGINT AS n
      FROM fr WHERE fbin >= 0 GROUP BY ALL) TO {out('joint')};
COPY (SELECT '{m}' AS MGS, base, cat, cc, least(floor(score / 5), 50)::INTEGER * 5 AS sbin, count(*)::BIGINT AS n
      FROM t WHERE clean = 1 GROUP BY ALL) TO {out('score')};
COPY (SELECT '{m}' AS MGS, base, count(*)::BIGINT AS sites_any, count(*) FILTER (WHERE c = 1)::BIGINT AS sites_clean_any
      FROM (SELECT contig, pos, strand, base, max(clean) AS c FROM t GROUP BY ALL) GROUP BY ALL) TO {out('sites')};
CREATE TABLE folds AS SELECT sample, (row_number() OVER (ORDER BY sample) % 2)::UTINYINT AS fold
      FROM (SELECT DISTINCT sample FROM t);
CREATE TABLE ref AS
SELECT contig, pos, strand, base, fold,
       CASE WHEN avg((cat = 'called')::INTEGER) >= 0.9 THEN 'methylated' ELSE 'unmethylated' END AS kind
FROM t JOIN folds USING (sample) WHERE clean = 1 AND coverage >= {a.ref_cov}
GROUP BY contig, pos, strand, base, fold
HAVING count(*) >= {a.ref_n} AND (avg((cat = 'called')::INTEGER) >= 0.9
       OR (sum((cat <> 'none')::INTEGER) = 0 AND hash(contig, pos, strand) % 20 = 0));
COPY (SELECT '{m}' AS MGS, t.sample, t.base, r.kind, t.cb, count(*)::BIGINT AS n,
             sum((t.cat = 'called')::INTEGER)::BIGINT AS n_called, sum((t.cat IN ('possible', 'other'))::INTEGER)::BIGINT AS n_possible
      FROM t JOIN folds f USING (sample)
             JOIN ref r ON r.contig = t.contig AND r.pos = t.pos AND r.strand = t.strand AND r.base = t.base AND r.fold <> f.fold
      WHERE t.clean = 1 GROUP BY ALL) TO {out('power')};
CREATE TABLE ps AS SELECT sample, base, avg((cat = 'called')::INTEGER) AS p
      FROM t WHERE clean = 1 AND coverage >= {a.sig_cov} GROUP BY ALL;
COPY (SELECT '{m}' AS MGS, base, count(*)::BIGINT AS sites, sum(n_s)::BIGINT AS obs, sum(k)::BIGINT AS called,
             sum(CASE WHEN v > 0 THEN (k - e) * (k - e) / v END) AS d, count(*) FILTER (WHERE v > 0)::BIGINT AS n_test,
             sum((k - e) * (k - e) - v) AS s_excess, sum(a * a - v) AS s_cross,
             count(*) FILTER (WHERE k = n_s)::BIGINT AS sites_all, count(*) FILTER (WHERE k = 0)::BIGINT AS sites_none,
             max(n_s)::BIGINT AS max_samples
      FROM (SELECT t.contig, t.pos, t.strand, t.base, count(*) AS n_s, sum((t.cat = 'called')::INTEGER) AS k,
                   sum(ps.p) AS e, sum(ps.p * (1 - ps.p)) AS v, sum(sqrt(ps.p * (1 - ps.p))) AS a
            FROM t JOIN ps USING (sample, base) WHERE t.clean = 1 AND t.coverage >= {a.sig_cov}
            GROUP BY ALL HAVING count(*) >= 2) GROUP BY ALL) TO {out('consistency')};
COPY (WITH s AS (SELECT contig, pos, strand, base, avg((cat = 'called')::INTEGER) AS f
                 FROM t WHERE clean = 1 AND coverage >= {a.sig_cov} GROUP BY ALL),
           w AS (SELECT base, contig, floor((pos - 1) / {a.window})::INTEGER AS win, count(*) AS n,
                        sum((f >= 0.5)::INTEGER) AS k FROM s GROUP BY ALL HAVING count(*) >= 10),
           r AS (SELECT base, sum(k) / sum(n) AS p FROM w GROUP BY ALL)
      SELECT '{m}' AS MGS, base, count(*)::BIGINT AS windows, sum(n)::BIGINT AS sites, sum(k)::BIGINT AS methylated,
             sum(CASE WHEN p > 0 AND p < 1 THEN (k - n * p) * (k - n * p) / (n * p * (1 - p)) END) AS d
      FROM w JOIN r USING (base) GROUP BY ALL) TO {out('spatial')};
"""


def mgs_files(d):
    fs = sorted(glob.glob(os.path.join(d, "*.parquet")))
    return fs, hashlib.md5(json.dumps([[os.path.basename(f), os.path.getsize(f), os.stat(f).st_mtime_ns] for f in fs])
                           .encode()).hexdigest()


def run_mgs(mgs, d, work, a, fp):
    dst = os.path.join(work, "parts", mgs)
    tmp = dst + ".tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    t0 = time.time()
    p = subprocess.run([a.duckdb, "-bail", os.path.join(tmp, "work.duckdb")], input=mgs_sql(mgs, os.path.join(d, "*.parquet"), tmp, a),
                       text=True, capture_output=True)
    if p.returncode != 0 or not all(os.path.isfile(os.path.join(tmp, x + ".parquet")) for x in PARTS):
        return mgs, f"duckdb failed: {(p.stderr or p.stdout).strip()[-600:]}"
    for x in ("work.duckdb", "work.duckdb.wal", "duck.tmp"):
        x = os.path.join(tmp, x)
        shutil.rmtree(x, ignore_errors=True) if os.path.isdir(x) else (os.path.exists(x) and os.remove(x))
    with open(os.path.join(tmp, "fp"), "w") as fh:
        fh.write(fp + "\n")
    shutil.rmtree(dst, ignore_errors=True)
    os.rename(tmp, dst)
    return mgs, f"done in {time.time() - t0:.0f} s"


def genome_counts(path):
    """A / C sites on both strands of a FASTA: A sites = A + T, C sites = C + G"""
    n = dict.fromkeys("ACGT", 0)
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if not line.startswith(">"):
                s = line.strip().upper()
                for b in "ACGT":
                    n[b] += s.count(b)
    return n["A"] + n["T"], n["C"] + n["G"], sum(n.values())


def summarise(a):
    root = os.path.join(a.out, "mgs-rep")
    if not os.path.isdir(root):
        die(f"no mgs-rep results in {a.out} (expected {root}/<MGS>/<sample>.parquet)")
    if shutil.which(a.duckdb) is None and not os.access(a.duckdb, os.X_OK):
        die(f"DuckDB CLI not found: {a.duckdb} (--duckdb, or MMT_DUCKDB)")
    work = a.work
    os.makedirs(os.path.join(work, "parts"), exist_ok=True)
    settings = (f"v{VERSION}|ref{a.min_ref_pct}|ctx{a.max_ctx_pct}|refcov{a.ref_cov}|refn{a.ref_n}|sig{a.sig_cov}"
                f"|win{a.window}")

    # ---- per-sample tables MMT keeps next to the results (presence, read counts, chemistry, kinetics model)
    units, samples = [], []
    for f in sorted(glob.glob(os.path.join(root, "stats", "*.run.tsv"))):
        kv = read_kv(f)
        s = kv.get("sample", os.path.basename(f)[:-8])
        samples.append([s, f"{kv.get('binding_kit', 'NA')}/{kv.get('sequencing_kit', 'NA')}",
                        kv.get("kinetics_model", "NA"), kv.get("movie", "NA")])
    for f in sorted(glob.glob(os.path.join(root, "stats", "*.tsv"))):
        if f.endswith(".run.tsv"):
            continue
        s = os.path.basename(f)[:-4]
        with open(f) as fh:
            hd = fh.readline().rstrip("\n").split("\t")
            ix = {k: hd.index(k) for k in ("MGS", "status", "molecules", "aligned", "rows", "rep_MAG") if k in hd}
            for line in fh:
                x = line.rstrip("\n").split("\t")
                units.append([s] + [x[ix[k]] if k in ix and ix[k] < len(x) else "NA"
                                    for k in ("MGS", "status", "molecules", "aligned", "rows", "rep_MAG")])
    with open(os.path.join(work, "units.tsv"), "w") as o:
        o.write("sample\tMGS\tstatus\tmolecules\taligned\trows\trep_MAG\n")
        o.writelines("\t".join(r) + "\n" for r in units)
    with open(os.path.join(work, "samples.tsv"), "w") as o:
        o.write("sample\tchemistry\tkinetics_model\tmovie\n")
        o.writelines("\t".join(r) + "\n" for r in samples)
    log(f"{len(samples)} samples, {len(units)} sample x MGS rows from {root}/stats")

    # ---- representative A / C sites (optional): GeneCat's own FASTA of the MGS representative (MMT's naming rule)
    if a.genecat:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
        from mf4 import rep_fasta_name
        rep = {}
        for r in units:
            if r[6] not in ("", "NA"):
                rep.setdefault(r[1], r[6])
        gdir = os.path.join(a.genecat, f"Bin_{a.binner}", "Genomes", "MGS_ctg")
        old = {}
        gf = os.path.join(work, "genome.tsv")
        if os.path.isfile(gf):
            with open(gf) as fh:
                next(fh)
                for line in fh:
                    x = line.rstrip("\n").split("\t")
                    old[(x[0], x[1], x[6])] = x
        rows, miss = [], 0
        for mgs, mag in sorted(rep.items()):
            fa = os.path.join(gdir, rep_fasta_name(mgs, mag))
            if not os.path.isfile(fa):
                miss += 1
                continue
            st = f"{os.path.getsize(fa)}:{os.stat(fa).st_mtime_ns}"
            if (mgs, mag, st) in old:
                rows.append(old[(mgs, mag, st)]); continue
            na, nc, L = genome_counts(fa)
            rows.append([mgs, mag, str(na), str(nc), str(L), fa, st])
        with open(gf, "w") as o:
            o.write("MGS\trep_MAG\tA_sites\tC_sites\tlength\tfasta\tstamp\n")
            o.writelines("\t".join(r) + "\n" for r in rows)
        log(f"representative genomes: {len(rows)} counted, {miss} FASTA missing in {gdir}")

    # ---- per MGS: skip unchanged, drop gone, run the rest (largest first: shortest tail)
    todo, keep = [], set()
    for d in sorted(glob.glob(os.path.join(root, "*", ""))):
        mgs = os.path.basename(os.path.dirname(d))
        if mgs == "stats" or (a.select and mgs not in a.select):
            continue
        fs, fp = mgs_files(d)
        if not fs:
            continue
        keep.add(mgs)
        fp = settings + "|" + fp
        pf = os.path.join(work, "parts", mgs, "fp")
        if os.path.isfile(pf) and open(pf).read().strip() == fp:
            continue
        todo.append((sum(os.path.getsize(f) for f in fs), mgs, d, fp))
    if not a.select:
        for d in glob.glob(os.path.join(work, "parts", "*")):
            if os.path.basename(d) not in keep:
                shutil.rmtree(d, ignore_errors=True)
                log(f"dropped {os.path.basename(d)} (no longer in {root})")
    todo.sort(reverse=True)
    log(f"{len(keep)} MGS with results ; {len(todo)} to summarise ({sum(t[0] for t in todo) / 2**30:.1f} GB), "
        f"{len(keep) - len(todo)} unchanged")
    fails = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_mgs, mgs, d, work, a, fp) for _, mgs, d, fp in todo]
        for i, fu in enumerate(cf.as_completed(futs), 1):
            mgs, msg = fu.result()
            log(f"{i}/{len(todo)} {mgs}: {msg}")
            if not msg.startswith("done"):
                fails.append(mgs)
    if fails:
        die(f"{len(fails)} MGS failed: {' '.join(fails[:20])} - fix and rerun (finished MGS are kept)")
    log(f"summaries complete in {work}")


# ============================================================================ plot: data
def load(work):
    import pandas as pd, pyarrow.parquet as pq
    parts = {}
    for x in PARTS:
        fs = sorted(glob.glob(os.path.join(work, "parts", "*", x + ".parquet")))
        if not fs:
            die(f"no {x}.parquet in {work}/parts - run summarise first")
        parts[x] = pd.concat([pq.read_table(f).to_pandas() for f in fs], ignore_index=True)
    units = pd.read_csv(os.path.join(work, "units.tsv"), sep="\t", dtype=str)
    samples = pd.read_csv(os.path.join(work, "samples.tsv"), sep="\t", dtype=str)
    gf = os.path.join(work, "genome.tsv")
    genome = pd.read_csv(gf, sep="\t", dtype={"MGS": str}) if os.path.isfile(gf) else None
    return parts, units, samples, genome


def bh(p):
    """Benjamini-Hochberg q-values (NaN stays NaN)"""
    import numpy as np
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        v = p[ok]
        o = np.argsort(v)
        q = np.minimum.accumulate((v[o] * len(v) / np.arange(1, len(v) + 1))[::-1])[::-1]
        r = np.empty(len(v))
        r[o] = np.minimum(q, 1)
        out[ok] = r
    return out


def bin_mid(cb):
    """representative coverage of a fine bin (exact up to 30, geometric middle above)"""
    import numpy as np
    cb = np.asarray(cb, float)
    up = {e: (COV_EDGES[i + 1] if i + 1 < len(COV_EDGES) else e * 2) for i, e in enumerate(COV_EDGES)}
    hi = np.array([up.get(int(c), c + 1) for c in cb], float)
    return np.where(cb <= 30, cb, np.sqrt(cb * (hi - 1)))


def isotonic(y, w):
    """weighted pool-adjacent-violators: the closest non-decreasing fit (detection cannot fall with coverage)"""
    blocks = []
    for yi, wi in zip(y, w):
        blocks.append([yi, wi, 1])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            v2, w2, c2 = blocks.pop()
            v1, w1, c1 = blocks.pop()
            blocks.append([(v1 * w1 + v2 * w2) / (w1 + w2), w1 + w2, c1 + c2])
    return [v for v, _, c in blocks for _ in range(c)]


def n95_from_power(pw, min_n):
    """lowest fine bin where the isotonic fit of the call rate (bins with >= min_n sites) reaches 95 % ; None = never"""
    g = pw.groupby("cb")[["n", "n_called"]].sum().sort_index()
    g = g[g["n"] >= min_n]
    if g.empty:
        return None
    hit = [e for e, f in zip(g.index, isotonic((g["n_called"] / g["n"]).values, g["n"].values)) if f >= 0.95]
    return int(hit[0]) if hit else None


def analyse(parts, units, samples, genome, a):
    """-> dict of tables every figure and card is drawn from"""
    import numpy as np, pandas as pd
    cube, frac, power = parts["cube"], parts["frac"], parts["power"]
    chem = dict(zip(samples["sample"], samples["chemistry"]))
    power = power.assign(chemistry=power["sample"].map(chem).fillna("NA"))

    # ---- n95 per base (pooled over MGS and samples), also per chemistry for the figure
    n95, rows = {}, []
    for b in "AC":
        pm = power[(power["base"] == b) & (power["kind"] == "methylated")]
        auto = n95_from_power(pm, a.power_min_n)
        n95[b] = a.n95.get(b, auto)
        rows.append([b, "all", auto, n95[b], int(pm["n"].sum())])
        for c, g in pm.groupby("chemistry"):
            rows.append([b, c, n95_from_power(g, a.power_min_n), None, int(g["n"].sum())])
    n95_tab = pd.DataFrame(rows, columns=["base", "chemistry", "n95_estimated", "n95_used", "reference_rows"])
    for b in "AC":
        if n95[b] is None:
            log(f"WARNING: {b}: reference-methylated sites never reach 95 % calls - no site counts as unmethylated")
            n95[b] = 10 ** 9

    # ---- per MGS x sample x base: tiers, states, coverage
    c = cube.assign(cbv=cube["cb"].astype(int), possible=cube["cat"].isin(["possible", "other"]))
    c["n95"] = c["base"].map(n95)
    key = ["MGS", "sample", "base"]
    cl = c[c["clean"] == 1]
    agg = lambda m: cl[m].groupby(key)["n"].sum()
    u = pd.DataFrame({
        "sites": c.groupby(key)["n"].sum(),
        "masked": c[c["clean"] == 0].groupby(key)["n"].sum(),
        "clean": cl.groupby(key)["n"].sum(),
        "cov_sum": cl.groupby(key)["cov_sum"].sum(),
        "tier_low": agg(cl["cbv"] < 11),
        "tier_mid": agg((cl["cbv"] >= 11) & (cl["cbv"] < cl["n95"])),
        "tier_high": agg(cl["cbv"] >= cl["n95"]),
        "called": agg(cl["cat"] == "called"),
        "possible": agg(cl["possible"]),
        "unmethylated": agg((cl["cat"] == "none") & (cl["cbv"] >= cl["n95"])),
        "undetermined": agg((cl["cat"] == "none") & (cl["cbv"] < cl["n95"])),
        "other_type": agg(cl["cat"] == "other"),
        "called_low": agg((cl["cat"] == "called") & (cl["cbv"] < 11)),
        "called_hi": agg((cl["cat"] == "called") & (cl["cbv"] >= cl["n95"])),
        "possible_hi": agg(cl["possible"] & (cl["cbv"] >= cl["n95"])),
    }).fillna(0).reset_index()
    # median coverage of clean sites (bin lower edges ; exact up to 30)
    s = cl.groupby(key + ["cbv"])["n"].sum().reset_index().sort_values(key + ["cbv"])
    s["cum"] = s.groupby(key)["n"].cumsum()
    s["tot"] = s.groupby(key)["n"].transform("sum")
    med = s[s["cum"] >= s["tot"] / 2].groupby(key)["cbv"].first().rename("cov_median")
    u = u.merge(med, on=key, how="left")
    u["cov_mean"] = u["cov_sum"] / u["clean"].where(u["clean"] > 0)
    for st in ("called", "possible", "unmethylated", "undetermined"):
        u["f_" + st] = u[st] / u["clean"].where(u["clean"] > 0)
    for t in ("low", "mid", "high"):
        u["f_tier_" + t] = u["tier_" + t] / u["clean"].where(u["clean"] > 0)
    hi_ok = u["tier_high"].where(u["tier_high"] >= a.min_hi_sites)     # coverage-controlled: sites with >= n95 reads
    u["f_called_hi"] = u["called_hi"] / hi_ok
    u["f_possible_hi"] = u["possible_hi"] / hi_ok
    u["f_unmethylated_hi"] = u["unmethylated"] / hi_ok
    u["f_masked"] = u["masked"] / u["sites"].where(u["sites"] > 0)

    # ---- denominators and usability (A+C together)
    sites = parts["sites"]
    if genome is not None and len(genome):
        den = pd.concat([genome[["MGS"]].assign(base="A", den=genome["A_sites"]),
                         genome[["MGS"]].assign(base="C", den=genome["C_sites"])])
        den_src = "representative genome"
    else:
        den = sites.rename(columns={"sites_clean_any": "den"})[["MGS", "base", "den"]]
        den_src = "union of clean sites over samples"
    u = u.merge(den, on=["MGS", "base"], how="left")
    u["breadth"] = u["clean"] / u["den"]
    tot = u.groupby(["MGS", "sample"]).agg(clean_AC=("clean", "sum"), den_AC=("den", "sum")).reset_index()
    tot["usable"] = (tot["clean_AC"] >= a.min_sites) & (tot["clean_AC"] / tot["den_AC"] >= a.min_breadth)
    u = u.merge(tot[["MGS", "sample", "usable"]], on=["MGS", "sample"])
    u["chemistry"] = u["sample"].map(chem).fillna("NA")

    # ---- coverage skew per usable unit: logit(call rate) ~ log2(coverage), bins >= n95 with >= 50 sites, WLS
    hi = cl[cl["cbv"] >= cl["n95"]]
    cc = pd.DataFrame({"n": hi.groupby(key + ["cbv"])["n"].sum(),
                       "k": hi[hi["cat"] == "called"].groupby(key + ["cbv"])["n"].sum()}).fillna(0).reset_index()
    cc = cc[cc["n"] >= 50]
    cc["p"] = (cc["k"] + 0.5) / (cc["n"] + 1)
    cc["y"] = np.log(cc["p"] / (1 - cc["p"]))
    cc["x"] = np.log2(bin_mid(cc["cbv"]))
    cc["w"] = cc["n"] * cc["p"] * (1 - cc["p"])
    sl = []
    for k_, g in cc.groupby(key):
        if len(g) < 3 or g["x"].max() - g["x"].min() < 1:
            continue
        w, x, y = g["w"].values, g["x"].values, g["y"].values
        xm, ym = np.average(x, weights=w), np.average(y, weights=w)
        sxx = np.sum(w * (x - xm) ** 2)
        b1 = np.sum(w * (x - xm) * (y - ym)) / sxx
        res = y - ym - b1 * (x - xm)
        disp = max(1.0, np.sum(w * res ** 2) / max(len(g) - 2, 1))   # overdispersion (never below binomial)
        se = math.sqrt(disp / sxx)
        sl.append(list(k_) + [b1, se, len(g)])
    slope = pd.DataFrame(sl, columns=key + ["slope", "se", "bins"])
    if len(slope):
        from scipy.stats import norm
        slope["p"] = 2 * norm.sf(np.abs(slope["slope"] / slope["se"]))
        slope["q"] = bh(slope["p"].values)
        slope = slope.merge(u[key + ["usable"]], on=key)
        slope = slope[slope["usable"]].drop(columns="usable")

    # ---- frac classes (coverage >= n95) and CI classes per unit
    f = frac.assign(cbv=frac["cb"].astype(int))
    f["n95"] = f["base"].map(n95)
    fh = f[(f["fbin"] >= 0) & (f["cbv"] >= f["n95"])]
    fc = pd.DataFrame({
        "frac_n": fh.groupby(key)["n"].sum(),
        "full": fh[fh["fbin"] >= 18].groupby(key)["n"].sum(),
        "high_partial": fh[(fh["fbin"] >= 10) & (fh["fbin"] < 18)].groupby(key)["n"].sum(),
        "low_partial": fh[fh["fbin"] < 10].groupby(key)["n"].sum(),
        "called_no_frac": f[f["fbin"] < 0].groupby(key)["n"].sum(),
    }).fillna(0).reset_index()
    u = u.merge(fc, on=key, how="left").fillna({x: 0 for x in ("frac_n", "full", "high_partial", "low_partial", "called_no_frac")})
    for x in ("full", "high_partial", "low_partial"):
        u["f_" + x] = u[x] / u["frac_n"].where(u["frac_n"] > 0)

    # ---- per MGS (means over usable samples ; coverage over all samples with rows)
    un = units.copy()
    un["has_rows"] = un["status"].isin(["done", "linked"])
    pres = un.groupby("MGS").agg(samples_listed=("sample", "nunique"), samples_with_rows=("has_rows", "sum")).reset_index()
    rows = []
    for (mgs, b), g in u.groupby(["MGS", "base"]):
        gu = g[g["usable"]]
        r = {"MGS": mgs, "base": b, "mod": MOD[b], "samples_with_table": g["sample"].nunique(), "samples_usable": len(gu),
             "cov_median_of_samples": g["cov_median"].median(), "cov_q25": g["cov_median"].quantile(.25),
             "cov_q75": g["cov_median"].quantile(.75), "cov_min": g["cov_median"].min(), "cov_max": g["cov_median"].max(),
             "breadth_median": g["breadth"].median(), "masked_pct": 100 * g["masked"].sum() / max(g["sites"].sum(), 1),
             "low_cov_pct": 100 * g["tier_low"].sum() / max(g["clean"].sum(), 1),
             "samples_low_cov_only": int(((g["tier_mid"] + g["tier_high"]) == 0).sum())}
        for st in ("called", "possible", "unmethylated", "undetermined", "full", "high_partial", "low_partial"):
            v = gu["f_" + st]
            r[st + "_pct"] = 100 * v.mean() if len(gu) else np.nan
            r[st + "_pct_sd"] = 100 * v.std() if len(gu) > 1 else np.nan
        r["called_or_possible_pct"] = r["called_pct"] + r["possible_pct"] if len(gu) else np.nan
        gh = gu[gu["f_called_hi"].notna()]
        r["samples_hi"] = len(gh)
        for st in ("called", "possible", "unmethylated"):
            r[st + "_hi_pct"] = 100 * gh["f_" + st + "_hi"].mean() if len(gh) else np.nan
        r["called_or_possible_hi_pct"] = r["called_hi_pct"] + r["possible_hi_pct"] if len(gh) else np.nan
        sg = slope[(slope["MGS"] == mgs) & (slope["base"] == b)] if len(slope) else slope
        r["skew_slope_median"] = sg["slope"].median() if len(sg) else np.nan
        r["skew_samples_tested"] = len(sg)
        r["skew_samples_q05"] = int((sg["q"] < 0.05).sum()) if len(sg) else 0
        rows.append(r)
    mg = pd.DataFrame(rows).merge(pres, on="MGS", how="left")
    if a.taxonomy:
        tx = pd.read_csv(a.taxonomy, sep="\t", header=None, names=["MGS", "taxonomy"], dtype=str, comment="#")
        mg = mg.merge(tx, on="MGS", how="left")
    # ---- non-random methylation: site consistency across samples, spatial clustering within samples
    from scipy.stats import chi2, wilcoxon
    cs = parts["consistency"].copy()
    # testable: >= 100 sites seen in >= 2 samples with some variance, and enough between-sample information (cross term
    # >= 0.01 per site; ~0 when only one sample makes calls at these sites -> rho would be noise / 0)
    ok = (cs["n_test"] >= 100) & (cs["s_cross"] >= 0.01 * cs["n_test"])
    cs["rho"] = (cs["s_excess"] / cs["s_cross"]).where(ok)
    cs["p_consistency"] = [chi2.sf(d, n) if t else np.nan for d, n, t in zip(cs["d"], cs["n_test"], ok)]
    # one test per MGS on its consensus map (samples share the genome: their window statistics are not independent)
    spm = parts["spatial"].rename(columns={"windows": "spatial_windows", "sites": "spatial_sites",
                                           "methylated": "spatial_methylated", "d": "spatial_d"})
    spm = spm[(spm["spatial_windows"] >= 3) & spm["spatial_d"].notna()].copy()
    spm["spatial_df"] = spm["spatial_windows"] - 1
    spm["log2_dispersion"] = np.log2(spm["spatial_d"] / spm["spatial_df"])
    spm["p_spatial"] = chi2.sf(spm["spatial_d"], spm["spatial_df"])
    st = cs.rename(columns={"sites": "consistency_sites", "max_samples": "consistency_max_samples",
                            "sites_all": "sites_called_in_all", "sites_none": "sites_called_in_none"})[
        ["MGS", "base", "consistency_sites", "consistency_max_samples", "sites_called_in_all", "sites_called_in_none", "rho",
         "p_consistency"]].merge(spm, on=["MGS", "base"], how="outer")
    st.insert(2, "mod", st["base"].map(MOD))
    for t in ("consistency", "spatial"):
        st["q_" + t] = np.nan
        for b in "AC":
            k_ = (st["base"] == b) & st["p_" + t].notna()
            st.loc[k_, "q_" + t] = bh(st.loc[k_, "p_" + t].values)
    rows = []
    for t, eff in (("consistency", "rho"), ("spatial", "log2_dispersion")):
        for b in "AC":
            d = st[(st["base"] == b) & st[eff].notna()]
            e = d[eff].values
            pw_ = wilcoxon(e, alternative="greater").pvalue if len(e) >= 5 and np.any(e != 0) else np.nan
            rows.append(dict(test=t, base=b, mod=MOD[b], effect=eff, MGS_tested=len(d), MGS_significant=int((d["q_" + t] < 0.05).sum()),
                             pct_significant=100 * (d["q_" + t] < 0.05).mean() if len(d) else np.nan,
                             effect_median=np.median(e) if len(e) else np.nan,
                             effect_q25=np.percentile(e, 25) if len(e) else np.nan,
                             effect_q75=np.percentile(e, 75) if len(e) else np.nan, p_across_MGS=pw_,
                             null="rho = 0 (random placement)" if t == "consistency" else "log2 dispersion = 0 (random along the genome)",
                             test_across_MGS="one-sided Wilcoxon signed-rank, MGS as units"))
    stats_across = pd.DataFrame(rows)
    mg = mg.merge(st[["MGS", "base", "rho", "q_consistency", "log2_dispersion", "q_spatial"]], on=["MGS", "base"], how="left")
    order = (mg.groupby("MGS")["samples_usable"].max().sort_values(ascending=False))
    return dict(n95=n95, n95_tab=n95_tab, units=u, mgs=mg, slope=slope, power=power, den_src=den_src, order=list(order.index),
                frac=f, joint=parts["joint"], score=parts["score"], cube=c, stats=st, stats_across=stats_across,
                min_sites=a.min_sites, min_breadth=a.min_breadth)


# ============================================================================ plot: style
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
BASE_COL = {"A": "#2a78d6", "C": "#d95926"}
BASE_LIGHT = {"A": "#86b6ef", "C": "#f0a070"}
BASE_DARK = {"A": "#104281", "C": "#8a3412"}
GREY_UNMETH, GREY_UNDET = "#8f8e88", "#e2e1dc"
TIER = ["#adaca6", "#7d7c77", "#4a4945"]
LS = ["-", "--", ":", "-."]
COARSE_LAB = ["3-10", "11-14", "15-19", "20-29", "30-49", "50-99", "100+"]


def style():
    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
        "axes.titlesize": 7.5, "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.labelsize": 7,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "axes.linewidth": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "axes.axisbelow": True,
        "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "legend.fontsize": 6.5, "legend.frameon": False, "lines.linewidth": 1.2, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.dpi": 300, "figure.dpi": 150, "text.color": INK})


def title(ax, letter, text):
    ax.set_title(f"{letter}   {text}" if letter else text, loc="left")


def logx(ax, lo=3, hi=1000):
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([t for t in (1, 3, 10, 30, 100, 300, 1000) if lo <= t <= hi]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())


def logy(ax):
    from matplotlib.ticker import FuncFormatter
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def stack_bars(ax, x, ys, colors, labels):
    import numpy as np
    bot = np.zeros(len(x))
    for y, c, lab in zip(ys, colors, labels):
        y = np.nan_to_num(np.asarray(y, float))
        ax.bar(x, y, bottom=bot, width=1.0, color=c, label=lab, linewidth=0)
        bot += y


def rank_axis(ax, n, what, label=True):
    ax.set_xlim(-0.5, max(n - 0.5, 0.5))
    ax.set_xlabel(f"MGS, ranked by {what} (n = {n})" if label else "")
    ax.grid(axis="x", visible=False)


def violins(ax, data, colors, labels):
    """violin + median dot + IQR bar per group (groups with < 2 values: dot only)"""
    import numpy as np
    pos = list(range(1, len(data) + 1))
    for p, d, c in zip(pos, data, colors):
        d = np.asarray(d, float)
        d = d[~np.isnan(d)]
        if len(d) > 1 and np.ptp(d) > 0:
            vp = ax.violinplot([d], positions=[p], widths=0.8, showextrema=False)
            vp["bodies"][0].set_facecolor(c); vp["bodies"][0].set_alpha(0.8); vp["bodies"][0].set_edgecolor("none")
        if len(d):
            q = np.percentile(d, [25, 50, 75])
            ax.plot([p, p], [q[0], q[2]], color=INK, lw=1.2, solid_capstyle="butt")
            ax.plot(p, q[1], "o", ms=3.5, color="white", mec=INK, mew=0.8, zorder=3)
    ax.set_xticks(pos)
    ax.set_xticklabels(labels)
    ax.set_xlim(0.4, len(data) + 0.6)
    ax.grid(axis="x", visible=False)


def coarse_of(cbv):
    import numpy as np
    e = np.asarray(cbv)
    return np.searchsorted(np.array(COARSE), e, side="right") - 1


def save(fig, report, name, tables=None):
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(report, "tables"), exist_ok=True)
    fig.savefig(os.path.join(report, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(report, name + ".png"), bbox_inches="tight")
    for k, df in (tables or {}).items():
        df.to_csv(os.path.join(report, "tables", f"{name}_{k}.tsv"), sep="\t", index=False)
    plt.close(fig)


def rate_curves(cube, b, cat_mask):
    """per MGS: rate of a state among clean sites per fine coverage bin"""
    import pandas as pd
    c = cube[(cube["base"] == b) & (cube["clean"] == 1)]
    return pd.DataFrame({"n": c.groupby(["MGS", "cbv"])["n"].sum(),
                         "k": c[cat_mask(c)].groupby(["MGS", "cbv"])["n"].sum()}).fillna(0).reset_index()


# ============================================================================ plot: figures
def fig1(R, report):
    """data quality: samples per MGS, coverage, coverage tiers, masked sites (x = MGS ranked by usable samples)"""
    import numpy as np, matplotlib.pyplot as plt
    mg = R["mgs"]
    fig = plt.figure(figsize=(7.1, 7.6))
    gs = fig.add_gridspec(4, 2, height_ratios=[0.8, 1, 1, 0.8])
    m0 = mg[mg["base"] == "A"].set_index("MGS").reindex(R["order"]).reset_index()
    m0 = m0[m0["samples_with_table"] > 0]
    x0 = np.arange(len(m0))
    ax = fig.add_subplot(gs[0, :])
    ax.bar(x0, m0["samples_with_table"], width=1.0, color=GREY_UNDET, label="samples with a table", linewidth=0)
    ax.bar(x0, m0["samples_usable"], width=1.0, color=INK2, label=f"usable (>= {R['min_sites']:,} clean A+C sites, "
                                                                  f">= {R['min_breadth']:.0%} breadth)", linewidth=0)
    ax.set_ylabel("samples")
    title(ax, "a", "Samples per MGS")
    ax.legend(loc="upper right")
    rank_axis(ax, len(m0), "usable samples")
    tabs = {"per_MGS": mg}
    for j, b in enumerate("AC"):
        m = mg[mg["base"] == b].set_index("MGS").reindex(m0["MGS"]).reset_index()
        x = np.arange(len(m))
        ax = fig.add_subplot(gs[1, j])
        ax.fill_between(x, m["cov_min"], m["cov_max"], step="mid", color=BASE_LIGHT[b], alpha=0.35, linewidth=0, label="min-max")
        ax.fill_between(x, m["cov_q25"], m["cov_q75"], step="mid", color=BASE_LIGHT[b], alpha=0.9, linewidth=0, label="IQR")
        ax.step(x, m["cov_median_of_samples"], where="mid", color=BASE_DARK[b], lw=1.0, label="median")
        logy(ax)
        ax.axhline(R["n95"][b], color=INK, lw=0.7, ls="--", label="n95")
        ax.axhline(11, color=INK2, lw=0.6, ls=":", label="11 (frac from here)")
        ax.set_ylabel("median coverage per sample\n(clean sites, reads / strand)")
        title(ax, "bc"[j], f"{b} sites ({MOD[b]}): coverage (n95 = {R['n95'][b]})")
        if j == 1:
            ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), title="per MGS, over samples", title_fontsize=6)
        rank_axis(ax, len(m), "usable samples", label=False)
        ax = fig.add_subplot(gs[2, j])
        t = R["units"][R["units"]["base"] == b].groupby("MGS")[["tier_low", "tier_mid", "tier_high", "clean"]].sum()
        t = t.reindex(m["MGS"]).fillna(0)
        tt = t[["tier_low", "tier_mid", "tier_high"]].div(t["clean"].where(t["clean"] > 0), axis=0).fillna(0) * 100
        stack_bars(ax, x, [tt["tier_high"], tt["tier_mid"], tt["tier_low"]], TIER[::-1],
                   [">= n95 reads", "11 to n95 - 1", "3-10 reads (no frac possible)"])
        ax.set_ylim(0, 100)
        ax.set_ylabel("% of clean sites (all samples)")
        title(ax, "de"[j], f"{b}: coverage tiers of clean sites")
        rank_axis(ax, len(m), "usable samples", label=False)
        if j == 1:
            ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
        tabs[f"{b}_tiers"] = tt.reset_index()
        ax = fig.add_subplot(gs[3, j])
        ax.fill_between(x, 0, m["masked_pct"], step="mid", color=GREY_UNMETH, alpha=0.5, linewidth=0)
        ax.step(x, m["masked_pct"], where="mid", color=INK2, lw=0.8)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("% of sites masked")
        title(ax, "fg"[j], f"{b}: masked sites (sample differs from rep.)")
        rank_axis(ax, len(m), "usable samples")
    fig.align_ylabels()
    fig.tight_layout(h_pad=0.8, w_pad=1.2)
    save(fig, report, "fig1_data_quality", tabs)


def fig2(R, report, a):
    """coverage effects: call rates vs coverage per MGS, calls at reference sites (power, false calls), skew slopes"""
    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    usable = R["units"][R["units"]["usable"]][["MGS", "sample", "base"]]
    cube = R["cube"].merge(usable, on=["MGS", "sample", "base"])
    fig, axs = plt.subplots(2, 4, figsize=(7.1, 3.9))
    tabs = {"n95": R["n95_tab"]}
    for i, b in enumerate("AC"):
        for j, (what, mask, key) in enumerate([(f"called {MOD[b]}", lambda c: c["cat"] == "called", "called"),
                                               ("modified_base", lambda c: c["possible"], "possible")]):
            ax = axs[i, j]
            d = rate_curves(cube, b, mask)
            d = d[d["n"] >= 200].copy()
            d["rate"] = 100 * d["k"] / d["n"]
            for _, g in d.groupby("MGS"):
                ax.plot(g["cbv"], g["rate"], color=BASE_LIGHT[b], lw=0.3, alpha=0.45)
            nm = d.groupby("cbv")["MGS"].nunique()
            q = d[d["cbv"].isin(nm[nm >= 3].index)].groupby("cbv")["rate"].quantile([.25, .5, .75]).unstack()
            if len(q):
                ax.fill_between(q.index, q[0.25], q[0.75], color=BASE_COL[b], alpha=0.25, linewidth=0, label="IQR over MGS")
                ax.plot(q.index, q[0.5], color=BASE_DARK[b], lw=1.4, label="median over MGS")
            ax.axvline(R["n95"][b], color=INK, lw=0.6, ls="--")
            ax.axvline(11, color=INK2, lw=0.5, ls=":")
            logx(ax)
            ax.set_ylim(bottom=0)
            ax.set_xlabel("coverage (reads / strand)")
            ax.set_ylabel(f"% of clean {b} sites")
            title(ax, "abcdefgh"[i * 4 + j], f"{b}: {what}")
            if len(q):
                ax.legend(loc="lower right" if j == 0 else "upper right", fontsize=5.5)
            tabs[f"{b}_{key}_rate"] = d
        ax = axs[i, 2]
        pw = R["power"][R["power"]["base"] == b]
        chems = sorted(pw["chemistry"].unique())
        for ci, chm in enumerate(chems):
            for kind, col in (("methylated", BASE_COL[b]), ("unmethylated", GREY_UNMETH)):
                g = pw[(pw["kind"] == kind) & (pw["chemistry"] == chm)].groupby("cb")[["n", "n_called"]].sum()
                g = g[g["n"] >= a.power_min_n]
                if len(g):
                    ax.plot(g.index.astype(int), 100 * g["n_called"] / g["n"], color=col, ls=LS[ci % 4], lw=1.0)
        ax.axhline(95, color=INK2, lw=0.5, ls=":")
        ax.axvline(R["n95"][b], color=INK, lw=0.6, ls="--")
        ax.text(R["n95"][b] * 1.08, 12, f"n95 = {R['n95'][b]}", fontsize=6, color=INK, va="center")
        logx(ax)
        ax.set_ylim(-2, 102)
        ax.set_xlabel("coverage (reads / strand)")
        ax.set_ylabel(f"% called {MOD[b]}")
        title(ax, "abcdefgh"[i * 4 + 2], f"{b}: calls at reference sites")
        h = [Line2D([], [], color=BASE_COL[b], lw=1.2, label="always methyl."),
             Line2D([], [], color=GREY_UNMETH, lw=1.2, label="never methyl.")]
        if len(chems) > 1:
            h += [Line2D([], [], color=INK2, ls=LS[k % 4], lw=1.0, label=c.split("/")[-1]) for k, c in enumerate(chems)]
        ax.legend(handles=h, fontsize=5, loc="lower right", bbox_to_anchor=(1.0, 0.18),
                  title="reference sites" + (" / seq. kit" if len(chems) > 1 else ""), title_fontsize=5)
        tabs[f"{b}_power"] = pw.groupby(["kind", "chemistry", "cb"])[["n", "n_called", "n_possible"]].sum().reset_index()
        ax = axs[i, 3]
        s = R["slope"][R["slope"]["base"] == b] if len(R["slope"]) else R["slope"]
        if len(s):
            lim = max(np.nanpercentile(np.abs(s["slope"]), 99), 0.1)
            bins = np.linspace(-lim, lim, 31)
            ax.hist(s["slope"].clip(-lim, lim), bins=bins, color=BASE_LIGHT[b], edgecolor="white", linewidth=0.4)
            ax.hist(s.loc[s["q"] < 0.05, "slope"].clip(-lim, lim), bins=bins, color=BASE_COL[b], edgecolor="white", linewidth=0.4)
            ax.axvline(0, color=INK, lw=0.6)
            ax.set_ylim(top=ax.get_ylim()[1] * 1.45)
            ax.text(0.03, 0.97, f"{len(s)} samples x MGS\n{(s['q'] < 0.05).mean() * 100:.0f} % q < 0.05 (dark)\n"
                                f"median {s['slope'].median():+.3f}", transform=ax.transAxes, va="top", fontsize=5.5, color=INK2)
        ax.set_xlabel("slope: logit(call rate) per doubling")
        ax.set_ylabel("samples x MGS")
        title(ax, "abcdefgh"[i * 4 + 3], f"{b}: coverage skew above n95")
        tabs[f"{b}_skew"] = s
    fig.tight_layout(w_pad=0.6, h_pad=1.0)
    save(fig, report, "fig2_coverage_effects", tabs)


def fig3(R, report):
    """what fraction of clean sites is methylated: state composition per MGS + coverage-controlled values across MGS"""
    import numpy as np, matplotlib.pyplot as plt
    mg = R["mgs"]
    fig = plt.figure(figsize=(7.1, 6.6))
    gs = fig.add_gridspec(4, 2, width_ratios=[3.0, 1.25], height_ratios=[1.5, 0.75, 1.5, 0.75])
    tabs = {}
    for i, b in enumerate("AC"):
        m = mg[(mg["base"] == b) & (mg["samples_usable"] > 0)].sort_values(["called_pct", "possible_pct"], ascending=False)
        x = np.arange(len(m))
        ax = fig.add_subplot(gs[2 * i, 0])
        stack_bars(ax, x, [m["called_pct"], m["possible_pct"], m["unmethylated_pct"], m["undetermined_pct"]],
                   [BASE_COL[b], BASE_LIGHT[b], GREY_UNMETH, GREY_UNDET],
                   [f"called {MOD[b]}", "possible (modified_base)", "unmethylated (no call, >= n95 reads)",
                    "undetermined (no call, < n95 reads)"])
        ax.set_ylim(0, 100)
        ax.set_ylabel(f"% of clean {b} sites\n(mean over usable samples)")
        title(ax, "ac"[i], f"{b} sites: state of every clean site, per MGS (n95 = {R['n95'][b]})")
        rank_axis(ax, len(m), "% called", label=False)
        hd, lb = ax.get_legend_handles_labels()
        ax2 = fig.add_subplot(gs[2 * i + 1, 0], sharex=ax)
        stack_bars(ax2, x, [m["called_pct"], m["possible_pct"]], [BASE_COL[b], BASE_LIGHT[b]], ["", ""])
        top = np.nanmax((m["called_pct"] + m["possible_pct"]).values) if len(m) else 1
        ax2.set_ylim(0, max(top * 1.1, 0.5))
        ax2.set_ylabel("zoom: % called\n+ possible")
        rank_axis(ax2, len(m), "% called")
        ax2.legend(hd, lb, loc="upper center", bbox_to_anchor=(0.5, -0.42), ncol=2, fontsize=5.5, handlelength=1.2,
                   columnspacing=1.2)
        ax = fig.add_subplot(gs[2 * i:2 * i + 2, 1])
        violins(ax, [m["called_hi_pct"], m["possible_hi_pct"], m["called_or_possible_hi_pct"]],
                [BASE_COL[b], BASE_LIGHT[b], BASE_DARK[b]], ["called", "possible", "either"])
        ax.set_ylim(bottom=0)
        ax.set_ylabel(f"% of clean {b} sites with >= n95 reads")
        title(ax, "bd"[i], f"{b}: well-covered sites")
        ax.text(0.98, 0.98, f"n = {m['called_hi_pct'].notna().sum()} MGS\nunmethylated = 100 - either",
                transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=INK2)
        tabs[f"{b}_per_MGS"] = m[["MGS", "samples_usable", "samples_hi", "called_pct", "possible_pct", "unmethylated_pct",
                                   "undetermined_pct", "called_hi_pct", "possible_hi_pct", "unmethylated_hi_pct",
                                   "called_or_possible_hi_pct"]]
    fig.tight_layout(h_pad=0.4, w_pad=1.2)
    save(fig, report, "fig3_landscape", tabs)


def fig4(R, report):
    """completeness: frac classes per MGS (coverage >= n95), frac by coverage (pooled), partial share across MGS"""
    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    mg, f = R["mgs"], R["frac"]
    fig = plt.figure(figsize=(7.1, 4.9))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.1, 1.55, 0.9])
    tabs = {}
    for i, b in enumerate("AC"):
        m = mg[(mg["base"] == b) & (mg["samples_usable"] > 0) & mg["full_pct"].notna()].sort_values("full_pct", ascending=False)
        x = np.arange(len(m))
        ax = fig.add_subplot(gs[i, 0])
        stack_bars(ax, x, [m["full_pct"], m["high_partial_pct"], m["low_partial_pct"]],
                   [BASE_DARK[b], BASE_COL[b], BASE_LIGHT[b]], ["full (frac >= 0.9)", "partial 0.5-0.9", "partial < 0.5"])
        ax.set_ylim(0, 100)
        ax.set_ylabel(f"% of called {MOD[b]} sites\n(>= n95 = {R['n95'][b]} reads)")
        title(ax, "ad"[i], f"{b}: how complete is {MOD[b]}, per MGS")
        rank_axis(ax, len(m), "% full")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=5.5, handlelength=1.2)
        ax = fig.add_subplot(gs[i, 1])
        g = f[(f["base"] == b) & (f["fbin"] >= 0)].assign(cc=lambda d: coarse_of(d["cbv"]))
        g = g.groupby(["cc", "fbin"])["n"].sum().unstack(fill_value=0).reindex(columns=range(20), fill_value=0)
        g = g[g.sum(axis=1) >= 100]
        if len(g):
            z = (g.T / g.sum(axis=1)).values * 100
            cmap = LinearSegmentedColormap.from_list("b", ["#ffffff", BASE_LIGHT[b], BASE_COL[b], BASE_DARK[b]])
            im = ax.pcolormesh(np.arange(len(g) + 1), np.linspace(0, 1, 21), z, cmap=cmap, shading="flat", rasterized=True)
            ax.set_xticks(np.arange(len(g)) + 0.5)
            ax.set_xticklabels([COARSE_LAB[k] for k in g.index], rotation=45, ha="right")
            cb = fig.colorbar(im, ax=ax, pad=0.03, fraction=0.07)
            cb.set_label("% of the column", fontsize=6)
            cb.ax.tick_params(labelsize=5.5)
            tabs[f"{b}_frac_by_coverage"] = g.rename(index=lambda k: COARSE_LAB[k]).reset_index()
        ax.grid(False)
        ax.set_xlabel("coverage (reads / strand)")
        ax.set_ylabel(f"frac ({MOD[b]})")
        title(ax, "be"[i], f"{b}: frac by coverage, all MGS")
        ax = fig.add_subplot(gs[i, 2])
        violins(ax, [100 - m["full_pct"]], [BASE_COL[b]], ["partial\n(frac < 0.9)"])
        ax.set_ylabel(f"% of called {MOD[b]} sites per MGS")
        title(ax, "cf"[i], "across MGS")
        if len(m):
            ax.text(0.97, 0.97, f"median {np.median(100 - m['full_pct']):.0f} %\nn = {len(m)} MGS", transform=ax.transAxes,
                    ha="right", va="top", fontsize=5.5, color=INK2)
        tabs[f"{b}_per_MGS"] = m[["MGS", "samples_usable", "full_pct", "high_partial_pct", "low_partial_pct"]]
    fig.tight_layout(w_pad=1.0, h_pad=1.0)
    save(fig, report, "fig4_completeness", tabs)


def fig5(R, report):
    """confidence: frac x score (pooled), score by state, is frac resolvable (Wilson CI) by coverage"""
    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, LinearSegmentedColormap
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.5), gridspec_kw={"width_ratios": [1.15, 1, 1.15]})
    tabs = {}
    usable = R["units"][R["units"]["usable"]][["MGS", "sample", "base"]]
    for i, b in enumerate("AC"):
        ax = axs[i, 0]
        j = R["joint"][R["joint"]["base"] == b].groupby(["fbin", "sbin"])["n"].sum().unstack(fill_value=0)
        if len(j):
            j = j.reindex(index=range(20), columns=range(0, 255, 5), fill_value=0)
            cmap = LinearSegmentedColormap.from_list("b", [BASE_LIGHT[b], BASE_COL[b], BASE_DARK[b], "#000000"])
            z = np.ma.masked_equal(j.values.T.astype(float), 0)
            im = ax.pcolormesh(np.linspace(0, 1, 21), np.arange(0, 260, 5), z, cmap=cmap,
                               norm=LogNorm(vmin=1, vmax=max(z.max(), 2)), shading="flat", rasterized=True)
            cb = fig.colorbar(im, ax=ax, pad=0.03, fraction=0.07)
            cb.set_label("called sites", fontsize=6)
            cb.ax.tick_params(labelsize=5.5)
            tabs[f"{b}_frac_x_score"] = j.reset_index()
        ax.grid(False)
        ax.set_xlabel(f"frac ({MOD[b]})")
        ax.set_ylabel("score (-10 log10 p, capped 250)")
        title(ax, "ad"[i], f"{b}: frac vs score")
        ax = axs[i, 1]
        s = R["score"][R["score"]["base"] == b].copy()
        s["state"] = s["cat"].replace({"other": "possible"})
        for st, col, lab in (("none", GREY_UNMETH, "no call"), ("possible", BASE_LIGHT[b], "modified_base"),
                             ("called", BASE_COL[b], f"called {MOD[b]}")):
            g = s[s["state"] == st].groupby("sbin")["n"].sum()
            if g.sum() == 0:
                continue
            g = g.reindex(range(0, 255, 5), fill_value=0)
            ax.step(g.index, 100 * g / g.sum(), where="post", color=col, lw=1.1, label=f"{lab} (n = {g.sum():,})")
        logy(ax)
        ax.set_xlabel("score")
        ax.set_ylabel("% of the state's clean sites")
        title(ax, "be"[i], f"{b}: score by state")
        ax.set_ylim(top=ax.get_ylim()[1] * 30)
        ax.legend(fontsize=5.2, loc="upper right")
        tabs[f"{b}_score_by_state"] = s.groupby(["state", "cc", "sbin"])["n"].sum().reset_index()
        ax = axs[i, 2]
        f = R["frac"].merge(usable, on=["MGS", "sample", "base"])
        f = f[(f["base"] == b) & (f["ci"] != "no_frac")].assign(cc=lambda d: coarse_of(d["cbv"]))
        g = f.groupby(["cc", "ci"])["n"].sum().unstack(fill_value=0).reindex(columns=["partial", "uncertain", "full"], fill_value=0)
        g = g[g.sum(axis=1) >= 100]
        if len(g):
            p = g.div(g.sum(axis=1), axis=0) * 100
            xs = np.arange(len(p))
            bot = np.zeros(len(p))
            for col, c_, lab in (("partial", BASE_DARK[b], "confident partial (upper CI < 0.9)"),
                                 ("uncertain", GREY_UNDET, "uncertain"),
                                 ("full", BASE_LIGHT[b], "consistent with full")):
                ax.bar(xs, p[col], bottom=bot, width=0.85, color=c_, label=lab, edgecolor="white", linewidth=0.5)
                bot += p[col].values
            ax.set_xticks(xs)
            ax.set_xticklabels([COARSE_LAB[k] for k in p.index], rotation=45, ha="right")
            tabs[f"{b}_ci_by_coverage"] = g.rename(index=lambda k: COARSE_LAB[k]).reset_index()
            ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=5.2)
        ax.set_ylim(0, 100)
        ax.grid(axis="x", visible=False)
        ax.set_xlabel("coverage (reads / strand)")
        ax.set_ylabel(f"% of called {MOD[b]} sites with frac")
        title(ax, "cf"[i], f"{b}: is frac resolvable?")
    fig.tight_layout(w_pad=0.8, h_pad=1.2)
    save(fig, report, "fig5_confidence", tabs)


def fig6(R, report):
    """non-random methylation: site consistency across samples + spatial clustering, per MGS and across MGS"""
    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    st, ac = R["stats"], R["stats_across"]
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.4), gridspec_kw={"width_ratios": [1.5, 1.5, 1.1]})
    for i, b in enumerate("AC"):
        s = st[st["base"] == b]
        for j, (eff, q, lab, null, ylab) in enumerate([
                ("rho", "q_consistency", "same sites in every sample?", 0.0, "site consistency rho\n(0 = random, 1 = same sites)"),
                ("log2_dispersion", "q_spatial", "clustered along the genome?", 0.0, "spatial log2 dispersion\n(0 = random)")]):
            ax = axs[i, j]
            d = s[s[eff].notna()].sort_values(eff, ascending=False)
            x = np.arange(len(d))
            sig = (d[q] < 0.05).values
            ax.axhline(null, color=INK2, lw=0.6)
            ax.scatter(x[~sig], d[eff].values[~sig], s=9, color=GREY_UNMETH, linewidths=0, label="q >= 0.05")
            ax.scatter(x[sig], d[eff].values[sig], s=9, color=BASE_COL[b], linewidths=0, label="q < 0.05")
            rank_axis(ax, len(d), eff.replace("_", " "))
            ax.set_ylabel(ylab)
            title(ax, "abcdef"[i * 3 + j], f"{b} ({MOD[b]}): {lab}")
            r = ac[(ac["base"] == b) & (ac["test"] == ("consistency" if j == 0 else "spatial"))]
            if len(r):
                r = r.iloc[0]
                ax.text(0.98, 0.97, f"{int(r['MGS_significant'])} / {int(r['MGS_tested'])} MGS q < 0.05\n"
                                    f"median {r['effect_median']:+.2f}\nacross MGS p = {r['p_across_MGS']:.1e}",
                        transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=INK2)
            ax.legend(loc="lower left", fontsize=5.5, markerscale=1.2)
        ax = axs[i, 2]
        d = s[s["rho"].notna() & s["log2_dispersion"].notna()]
        both = ((d["q_consistency"] < 0.05) & (d["q_spatial"] < 0.05)).values
        one = ((d["q_consistency"] < 0.05) ^ (d["q_spatial"] < 0.05)).values
        ax.axhline(0, color=INK2, lw=0.6)
        ax.axvline(0, color=INK2, lw=0.6)
        for msk, col, sz in ((~both & ~one, GREY_UNDET, 7), (one, BASE_LIGHT[b], 9), (both, BASE_COL[b], 11)):
            ax.scatter(d["rho"].values[msk], d["log2_dispersion"].values[msk], s=sz, color=col, linewidths=0)
        ax.set_xlabel("site consistency (rho)")
        ax.set_ylabel("spatial log2 dispersion")
        title(ax, "cf"[i], f"{b}: both tests")
        ax.legend(handles=[Line2D([], [], ls="", marker="o", ms=3, color=c, label=t) for c, t in
                           ((BASE_COL[b], "both q < 0.05"), (BASE_LIGHT[b], "one"), (GREY_UNDET, "neither"))],
                  fontsize=5.2, loc="upper left")
    fig.tight_layout(w_pad=0.8, h_pad=1.2)
    save(fig, report, "fig6_non_random", {"per_MGS": st, "across_MGS": ac})


def cards(R, report, a):
    """one page per MGS: rows A / C ; coverage per sample, state per sample, calls vs coverage, frac distribution"""
    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    u, mg, cube, f, st = R["units"], R["mgs"], R["cube"], R["frac"], R["stats"]
    order = [m for m in R["order"] if mg.loc[mg["MGS"] == m, "samples_usable"].max() >= a.cards_min_usable]
    if a.cards > 0:
        order = order[:a.cards]
    usable = u[u["usable"]][["MGS", "sample", "base"]]
    cu = cube.merge(usable, on=["MGS", "sample", "base"])
    ref = {}
    for b in "AC":   # cohort reference curve: median over MGS of the called rate per coverage bin (usable samples)
        d = rate_curves(cu, b, lambda c: c["cat"] == "called")
        d = d[d["n"] >= 200]
        ref[b] = (100 * d["k"] / d["n"]).groupby(d["cbv"]).median()
    frac_by, cu_by = dict(tuple(f.groupby("MGS"))), dict(tuple(cu.groupby("MGS")))
    fmt = lambda v, s="{:.2f}": ("NA" if v is None or (isinstance(v, float) and np.isnan(v))
                                 else "< 1e-300" if s.endswith("e}") and v == 0 else s.format(v))
    path = os.path.join(report, "mgs_cards.pdf")
    with PdfPages(path) as pdf:
        for mgs in order:
            fig, axs = plt.subplots(2, 4, figsize=(8.27, 5.4), gridspec_kw={"width_ratios": [1.2, 1.5, 1.3, 1.1]})
            info = mg[mg["MGS"] == mgs].set_index("base")
            sinfo = st[st["MGS"] == mgs].set_index("base")
            h = info.iloc[0]
            tax = h["taxonomy"] if "taxonomy" in info.columns and isinstance(h.get("taxonomy"), str) else ""
            lines = [f"{mgs}" + (f"   {tax}" if tax else ""),
                     f"samples with a table {int(h['samples_with_table'])}, usable {int(h['samples_usable'])}   |   median "
                     f"coverage {fmt(h['cov_median_of_samples'], '{:.0f}')}   |   n95 A {R['n95']['A']}, C {R['n95']['C']}"]
            for b in "AC":
                if b in info.index:
                    r, sr = info.loc[b], (sinfo.loc[b] if b in sinfo.index else {})
                    lines.append(f"{b}: >= n95 reads: called {fmt(r['called_hi_pct'])} % {MOD[b]}, possible {fmt(r['possible_hi_pct'])} %"
                                 f"  |  full {fmt(r['full_pct'], '{:.0f}')} % of called  |  site consistency rho "
                                 f"{fmt(sr.get('rho', np.nan))} (q {fmt(sr.get('q_consistency', np.nan), '{:.1e}')})  |  spatial "
                                 f"log2 dispersion {fmt(sr.get('log2_dispersion', np.nan))} (q {fmt(sr.get('q_spatial', np.nan), '{:.1e}')})")
            fig.text(0.01, 0.995, lines[0], fontsize=8, fontweight="bold", va="top")
            fig.text(0.01, 0.962, "\n".join(lines[1:]), fontsize=6, va="top", color=INK2, linespacing=1.5)
            fm, cm = frac_by.get(mgs), cu_by.get(mgs)
            for i, b in enumerate("AC"):
                g = u[(u["MGS"] == mgs) & (u["base"] == b)].sort_values("cov_median")
                y = np.arange(len(g))
                ax = axs[i, 0]
                ax.barh(y, g["cov_median"], color=[BASE_COL[b] if v else GREY_UNDET for v in g["usable"]], height=0.8, linewidth=0)
                ax.axvline(R["n95"][b], color=INK, lw=0.6, ls="--")
                logx(ax, 1, 1000)
                ax.set_yticks(y if len(g) <= 40 else [])
                ax.set_yticklabels(g["sample"], fontsize=4.5)
                ax.set_xlabel("median coverage (grey: not usable)")
                ax.grid(axis="y", visible=False)
                title(ax, "", f"{b}: coverage per sample")
                ax = axs[i, 1]
                left = np.zeros(len(g))
                for stt, col in (("called", BASE_COL[b]), ("possible", BASE_LIGHT[b]), ("unmethylated", GREY_UNMETH),
                                 ("undetermined", GREY_UNDET)):
                    v = 100 * g["f_" + stt].fillna(0).values
                    ax.barh(y, v, left=left, color=col, height=0.8, linewidth=0, label=stt)
                    left += v
                ax.set_xlim(0, 100)
                ax.set_yticks([])
                ax.set_xlabel(f"% of clean {b} sites")
                ax.grid(axis="y", visible=False)
                title(ax, "", f"{b}: state per sample (same order)")
                if i == 1:
                    ax.legend(fontsize=5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=4, handlelength=1, columnspacing=0.6)
                ax = axs[i, 2]
                if len(ref[b]):
                    ax.plot(ref[b].index, ref[b].values, color=GREY_UNMETH, lw=0.9, label="cohort median (called)")
                if cm is not None:
                    d = rate_curves(cm, b, lambda x_: x_["cat"] == "called")
                    d = d[d["n"] >= 50]
                    dp = rate_curves(cm, b, lambda x_: x_["possible"])
                    dp = dp[dp["n"] >= 50]
                    if len(d):
                        ax.plot(d["cbv"], 100 * d["k"] / d["n"], "o-", color=BASE_COL[b], ms=2, lw=0.9, label=f"called {MOD[b]}")
                    if len(dp):
                        ax.plot(dp["cbv"], 100 * dp["k"] / dp["n"], "s--", color=BASE_LIGHT[b], ms=1.8, lw=0.8, label="possible")
                ax.axvline(R["n95"][b], color=INK, lw=0.6, ls="--")
                logx(ax)
                ax.set_ylim(bottom=0)
                ax.set_xlabel("coverage (reads / strand)")
                ax.set_ylabel(f"% of clean {b} sites")
                title(ax, "", f"{b}: calls vs coverage")
                ax.legend(fontsize=5, loc="best")
                ax = axs[i, 3]
                if fm is not None:
                    fb = fm.merge(usable, on=["MGS", "sample", "base"])
                    fb = fb[(fb["base"] == b) & (fb["fbin"] >= 0)]
                    for sel, col, lab in ((fb["cbv"] >= R["n95"][b], BASE_COL[b], ">= n95"), (fb["cbv"] < R["n95"][b], BASE_LIGHT[b], "< n95")):
                        hh = fb[sel].groupby("fbin")["n"].sum().reindex(range(20), fill_value=0)
                        if hh.sum():
                            ax.step(np.arange(21) / 20, np.r_[100 * hh / hh.sum(), 100 * hh.iloc[-1] / hh.sum()], where="post",
                                    color=col, lw=1.0, label=f"{lab} (n = {hh.sum():,})")
                ax.set_xlim(0, 1)
                ax.set_xlabel(f"frac ({MOD[b]})")
                ax.set_ylabel("% of called sites")
                title(ax, "", f"{b}: frac distribution")
                ax.legend(fontsize=5, loc="upper left")
            fig.tight_layout(rect=(0, 0, 1, 0.86), w_pad=0.6, h_pad=0.9)
            pdf.savefig(fig)
            plt.close(fig)
    log(f"cards: {len(order)} MGS -> {path}")


def plot(a):
    style()
    parts, units, samples, genome = load(a.work)
    R = analyse(parts, units, samples, genome, a)
    os.makedirs(os.path.join(a.report, "tables"), exist_ok=True)
    R["mgs"].to_csv(os.path.join(a.report, "tables", "per_MGS.tsv"), sep="\t", index=False)
    R["units"].to_csv(os.path.join(a.report, "tables", "per_MGS_x_sample.tsv"), sep="\t", index=False)
    R["n95_tab"].to_csv(os.path.join(a.report, "tables", "n95.tsv"), sep="\t", index=False)
    log(f"n95: A {R['n95']['A']}, C {R['n95']['C']} ; usable MGS x sample {int(R['units']['usable'].sum() / 2)} ; "
        f"breadth denominator: {R['den_src']}")
    R["stats"].to_csv(os.path.join(a.report, "stats_per_MGS.tsv"), sep="\t", index=False)
    R["stats_across"].to_csv(os.path.join(a.report, "stats_across_MGS.tsv"), sep="\t", index=False)
    for fn in (fig1, fig3, fig4, fig5, fig6):
        fn(R, a.report)
    fig2(R, a.report, a)
    if a.cards != 0:
        cards(R, a.report, a)
    with open(os.path.join(a.report, "README.txt"), "w") as o:
        o.write(f"mmt_landscape.py v{VERSION} plot {time.strftime('%Y-%m-%d %H:%M')}\nWORK {os.path.abspath(a.work)}\n"
                f"n95 A {R['n95']['A']} C {R['n95']['C']} (tables/n95.tsv) ; breadth denominator: {R['den_src']}\n"
                f"usable: >= {a.min_sites} clean A+C sites and breadth >= {a.min_breadth}\n\n" + __doc__)
    log(f"report in {a.report}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("summarise", help="MMT results -> small summary tables (resumable)")
    s.add_argument("-o", "--out", required=True, help="MMT output folder (MMT -o): holds mgs-rep/")
    s.add_argument("-w", "--work", required=True, help="summary folder (created)")
    s.add_argument("--duckdb", default=os.environ.get("MMT_DUCKDB", "duckdb"), help="DuckDB CLI (default $MMT_DUCKDB or duckdb)")
    s.add_argument("--jobs", type=int, default=4, help="MGS at once")
    s.add_argument("--threads", type=int, default=4, help="DuckDB threads per MGS")
    s.add_argument("--mem-gb", type=int, default=16, help="DuckDB memory per MGS (spills to disk beyond)")
    s.add_argument("--genecat", help="GeneCat folder: representative genome sizes for breadth (optional)")
    s.add_argument("--binner", default="SB", help="GeneCat Bin_<binner> (default SB)")
    s.add_argument("--select", type=lambda x: set(x.split(",")), help="only these MGS (comma-separated)")
    s.add_argument("--min-ref-pct", type=float, default=95, help="clean: %% of reads with the reference base (95)")
    s.add_argument("--max-ctx-pct", type=int, default=10, help="clean: ctx_alt_pct below this (10)")
    s.add_argument("--ref-cov", type=int, default=30, help="reference sites: coverage of the defining samples (30)")
    s.add_argument("--ref-n", type=int, default=3, help="reference sites: samples per half (3)")
    s.add_argument("--sig-cov", type=int, default=30, help="non-randomness tests: sites with >= this many reads (30)")
    s.add_argument("--window", type=int, default=1000, help="spatial test: window (bp) along each contig (1000)")
    q = sp.add_parser("plot", help="summary tables -> figures, cards, tables")
    q.add_argument("-w", "--work", required=True)
    q.add_argument("-r", "--report", required=True, help="figure / table folder (created)")
    q.add_argument("--taxonomy", help="TSV MGS<TAB>label (no header) for card titles and per_MGS.tsv")
    q.add_argument("--n95", default="", type=lambda x: {k: int(v) for k, v in (y.split("=") for y in x.split(",") if y)},
                   help="override n95, e.g. A=20,C=25")
    q.add_argument("--min-sites", type=int, default=20000, help="usable: clean A+C sites (20000)")
    q.add_argument("--min-breadth", type=float, default=0.10, help="usable: fraction of the denominator (0.10)")
    q.add_argument("--min-hi-sites", type=int, default=1000,
                   help="coverage-controlled %% need this many clean sites with >= n95 reads in the sample (1000)")
    q.add_argument("--power-min-n", type=int, default=30, help="reference sites per coverage bin for power / n95 (30)")
    q.add_argument("--cards", type=int, default=-1, help="cards: -1 all, 0 none, N the first N MGS")
    q.add_argument("--cards-min-usable", type=int, default=1, help="cards only for MGS with >= this many usable samples")
    a = p.parse_args()
    summarise(a) if a.cmd == "summarise" else plot(a)


if __name__ == "__main__":
    main()
