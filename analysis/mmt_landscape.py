#!/usr/bin/env python3
"""MMT landscape: cohort-wide methylation statistics + publication figures from MMT MGScall --mgs-rep results.

  summarise  per MGS, two DuckDB CLI passes (the binary MMT uses) over <OUT>/mgs-rep/<MGS>/<sample>.parquet + per-site
             tests in python -> small tables in WORK/parts/<MGS>/, WORK/units.tsv, samples.tsv, patients.tsv (--map),
             genome.tsv (--genecat), settings.json. Pass 1 streams per-sample summaries from the files (nothing stored);
             pass 2 builds the patient x site table one genome chunk at a time (--chunk-rows, integer keys), so memory
             stays bounded whatever the MGS size (120 M rows / 60 samples: DuckDB peak 1.1 GB at --mem-gb 3). The log
             gives every MGS's DuckDB peak. Resumes: an MGS whose files and settings did not change is skipped ; MGS
             no longer in OUT are dropped from WORK ; a failed MGS keeps its SQL + DuckDB log in WORK/parts/<MGS>.tmp.
  plot       WORK -> REPORT: fig1..fig6 (+ fig7 with --run / --between) (pdf + png), mgs_cards.pdf (one page per MGS), stats_per_MGS.tsv,
             stats_across_MGS.tsv, constitutive_sites.parquet, tables/*.tsv (every plotted number)

  mmt_landscape.py summarise -o OUT -w WORK --map MAP [--run FILE] [--between N] [--duckdb BIN] [--jobs 2] [--threads 8]
                             [--mem-gb 24] [--select MGS,..] [--genecat GC --binner SC] [--min-cov 11]
  mmt_landscape.py plot -w WORK -r REPORT [--taxonomy TSV] [--min-sites 20000] [--min-breadth 0.1]

Every statistic is for A sites (-> m6A) and C sites (-> m4C) separately. G/T rows are never read.
  site          one base on one strand of the representative with an MMT row (>= 3 kinetic reads, MMT's minimum)
  clean         >= --min-ref-pct % of the reads carry the reference base, no position within MMT's context window where
                >= --max-ctx-pct % of the reads differ (ctx_alt_pct), no read clip in that window (clip_ctx_n = 0)
  graded        clean and >= --min-cov reads (default 11: ipdSummary reports frac only above 10 reads). One threshold for
                every analysis ; 3 to min-cov - 1 reads = not graded (counted, never classed).
  class         graded sites: near-complete unmethylated (no call, or frac < 5 %) | low 5-25 % | medium 25-75 % |
                high 75-95 % | near-complete >= 95 % (frac of an m6A call on A / m4C call on C) ; possible =
                modified_base (or the other type): a kinetic signal of unknown level, kept apart.
  patient       --map AssmblGrps (one person: time points and replicate cells of a person share it) ; a sample not in
                the map is its own patient. Per site and patient the best-covered graded sample counts (ties: sample
                name order), so time points / replicates never count twice.
  stratum       coverage of that sample at the site: min-cov-14, 15-19, 20-29, 30-49, >= 50 reads
  agreement     Krippendorff's alpha (ordinal) of the 5 classes across patients, sites with >= 2 patients: 0 = what
                random placement gives, 1 = every patient has the same class at every site. p: permutation (shuffle each
                patient's classes among its sites of the same stratum, on <= --alpha-sites sampled sites, --perm
                rounds), so coverage effects and each patient's class mix are part of the null. q: BH across MGS.
  constitutive  per site with >= --min-patients patients: k patients near-complete (or near-complete unmethylated) vs
                chance = sum over its patients of that patient's class share in the same stratum ; one-sided binomial
                p (conservative for this Poisson-binomial), BH across the MGS's sites -> always_methylated /
                always_unmethylated (q < 0.05 and >= 80 % of patients in that class) | other | few_patients
  spatial       are methylated sites (class >= low in >= half of a site's patients) spread at random along the genome?
                --window bp windows with >= 10 sites: sum((k - E)^2 / V) with E, V from the stratum-specific rates
                (coverage-adjusted), chi-square on windows - 1 ; log2 dispersion 0 = random, > 0 = clustered
                (methylation islands, or regions rich / poor in the target motifs), < 0 = evenly spaced. q: BH.
  n95           diagnostic only: coverage from which >= 95 % of reference-methylated sites are called (sites called in
                >= 90 % of >= --ref-n samples with >= --ref-cov reads in one half of the samples, measured in the other)
  usable        an MGS x sample with >= --min-sites graded A+C sites covering >= --min-breadth of the representative's A+C
                sites (--genecat), else of the union of clean sites seen in any sample
  frac CI       Wilson 95 % interval of frac with n = coverage: how many classes it spans (1, 2, >= 3). ipdSummary's
                fracLow / fracUp are not used: they come from an unseeded bootstrap and change between runs.
  variance      summarise --run FILE (lines 'sample TAB in-batch|across-batch'): within each AssmblGrps of --map, every
                two listed samples with the same label are compared on the sites graded in both (clean, >= min-cov reads
                in each) -> WORK/parts/<MGS>/variance.parquet (own fingerprint: adding or editing --run never redoes the
                summaries). plot -> variance_summary.tsv (per set x base: pooled over every pair and MGS, + medians over
                pairs x MGS with >= --var-min-sites sites), variance_by_coverage.tsv (by the lower coverage of the
                two), variance_pairs.tsv (per MGS x pair x base). Measures: same_class_pct (5 classes),
                methylated_discordant_pct (one sample < 5 %, the other >= 5 %), alpha (Krippendorff, ordinal, the pair as
                two raters), frac_rms_diff (no call = 0), frac_sd_per_sample = frac_rms_diff / sqrt 2 (the noise of one
                sample) ; sites_possible_either = sites where one sample has a 'possible' call or a call without frac.
  between       summarise --between N: per MGS, N disjoint pairs of DIFFERENT patients (the best-covered sample of each
                patient, the 2N best-covered patients, paired at random with a seed fixed by the MGS name) -> the same
                pair tables, set 'between-patient' (variance pass only).
  fig7          change per MGS: switched_pct = % of the sites methylated (>= 5 %) in either sample of a pair that are
                unmethylated (< 5 %) in the other, on sites with >= --pair-min-cov reads in both samples, pairs with
                >= --pair-min-meth such sites. Technical bands = mean +- 1 SD over replicate pairs x MGS (in-batch,
                across-batch) ; threshold = across-batch mean + 2 SD (two random patients are mostly from different
                runs, so across-batch is their technical floor) ; above_band = the MGS's median between-patient value
                exceeds it. Paired one-sided Wilcoxon across MGS: between-patient median > across-batch median.
                -> variance_evidence.tsv (per MGS x base), variance_bands.tsv (per set x base + the test).
"""
import argparse, collections, concurrent.futures as cf, glob, gzip, hashlib, itertools, json, math, os, random, shutil, \
    subprocess, sys, time

VERSION = "2"
COV_EDGES = list(range(3, 31)) + [31, 35, 40, 50, 60, 80, 100, 150, 200, 300, 500, 1000]   # fine coverage bins (lower edges)
COARSE = [3, 11, 15, 20, 30, 50, 100]                                                       # coarse bins (lower edges)
LEVELS = [0.05, 0.25, 0.75, 0.95]                                                           # class edges of frac
CLASSES = ["unmethylated", "low", "medium", "high", "near-complete"]
CLASS_LAB = ["near-complete unmethylated (< 5 %)", "low (5-25 %)", "medium (25-75 %)", "high (75-95 %)",
             "near-complete (>= 95 %)"]
PARTS = ("cube", "frac", "joint", "score", "sites", "contigs", "power", "npat", "spatial", "alpha", "site_summary",
         "track", "constitutive")
MOD = {"A": "m6A", "C": "m4C"}


def die(msg):
    sys.exit("ERROR: " + msg)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def case_bins(col, edges):
    """SQL: lower edge of the bin holding col (edges ascending ; below the first edge -> 0)"""
    return "CASE " + " ".join(f"WHEN {col} >= {e} THEN {e}" for e in reversed(edges)) + " ELSE 0 END"


def level_sql(x):
    """SQL: class 0-4 of a fraction"""
    # FLOAT edges: frac is stored as float32, so a frac printed 0.950 by ipdSummary is >= 0.95 (near-complete)
    return "CASE " + " ".join(f"WHEN {x} < {e}::FLOAT THEN {i}" for i, e in enumerate(LEVELS)) + " ELSE 4 END"


def strata(mc):
    return [mc] + [e for e in (15, 20, 30, 50) if e > mc]


def read_kv(path):
    """MMT run.tsv (key TAB value lines; first value of a repeated key wins)"""
    kv = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t", 1)
            if len(f) == 2 and f[0] not in kv:
                kv[f[0]] = f[1]
    return kv


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


def kalpha(K):
    """Krippendorff's alpha, ordinal metric, from a sites x classes count matrix (sites with < 2 values ignored)"""
    import numpy as np
    K = np.asarray(K, float)
    m = K.sum(1)
    K, m = K[m >= 2], m[m >= 2]
    if len(K) == 0:
        return np.nan
    W = K / (m - 1)[:, None]
    return alpha_coinc(W.T @ K - np.diag(W.sum(0)))       # coincidence matrix


def alpha_coinc(O):
    """Krippendorff's alpha (ordinal) from a coincidence matrix ; two samples compared site by site: O = C + C.T of their
    class cross-table C"""
    import numpy as np
    O = np.asarray(O, float)
    nc = O.sum(1)
    n = nc.sum()
    cum = np.concatenate([[0], np.cumsum(nc)])
    k = len(nc)
    d2 = np.zeros((k, k))
    for c in range(k):
        for e in range(c + 1, k):
            d2[c, e] = d2[e, c] = (cum[e + 1] - cum[c] - (nc[c] + nc[e]) / 2) ** 2
    de = (np.outer(nc, nc) * d2).sum()
    return 1 - (n - 1) * (O * d2).sum() / de if de > 0 else np.nan


# ============================================================================ summarise
BLOCK = 50000   # bp: unit of the genome chunks of pass 2


def sql_head(dst, a):
    mem = max(1, int(a.mem_gb * 0.7))     # DuckDB's own limit; the rest of the budget covers what it allocates outside it
    return (f"SET threads = {a.threads}; SET memory_limit = '{mem}GB'; SET preserve_insertion_order = false;\n"
            f"SET temp_directory = '{sq(os.path.join(dst, 'duck.tmp'))}';\n")


def sq(s):
    return s.replace("'", "''")


def view_sql(src, a, name="v"):
    """every A / C row of the MGS (>= 3 reads) with clean flag, call category, class and coverage bins - a view, never stored"""
    clean = (f"pile_depth > 0 AND n_ref >= {a.min_ref_pct / 100.0} * pile_depth AND ctx_alt_pct < {a.max_ctx_pct} "
             f"AND clip_ctx_n = 0")
    return f"""
CREATE OR REPLACE VIEW {name} AS
SELECT *, (CASE WHEN clean = 1 AND coverage >= {a.min_cov} THEN
             CASE WHEN cat = 'none' THEN 0 WHEN cat = 'called' AND frac IS NOT NULL THEN {level_sql('frac')} END
          END)::UTINYINT AS cls
FROM (SELECT sample, contig, pos, strand, base, coverage, score, frac,
             (CASE WHEN {clean} THEN 1 ELSE 0 END)::UTINYINT AS clean,
             CASE WHEN call = (CASE base WHEN 'A' THEN 'm6A' ELSE 'm4C' END) THEN 'called'
                  WHEN call = 'none' THEN 'none' WHEN call = 'modified_base' THEN 'possible' ELSE 'other' END AS cat,
             ({case_bins('coverage', COV_EDGES)})::USMALLINT AS cb,
             ({case_bins('coverage', COARSE)})::USMALLINT AS cc
      FROM read_parquet('{sq(src)}') WHERE base IN ('A', 'C') AND coverage >= 3);
"""


def pass1_sql(mgs, src, dst, a):
    """pass 1: per-sample summaries (small GROUP BYs streamed from the parquet files), the site list (one row per site),
    contig ids, graded rows per genome block, n95 diagnostic"""
    out = lambda name: f"'{sq(os.path.join(dst, name + '.parquet'))}' (FORMAT parquet, COMPRESSION zstd)"
    m = sq(mgs)
    # Wilson 95 % bound of frac with n = coverage: (p + z2/2n +- z sqrt(p(1-p)/n + z2/4n2)) / (1 + z2/n)
    wilson = lambda sg: (f"((frac::DOUBLE + 1.9208 / coverage {sg} 1.96 * sqrt(greatest(frac::DOUBLE * (1 - frac::DOUBLE), 0) "
                         f"/ coverage + 0.9604 / (coverage * coverage))) / (1 + 3.8416 / coverage))")
    return sql_head(dst, a) + view_sql(src, a) + f"""
COPY (SELECT '{m}' AS MGS, sample, base, cb, clean, cat, cls, count(*)::BIGINT AS n, sum(coverage)::BIGINT AS cov_sum
      FROM v GROUP BY ALL) TO {out('cube')};
COPY (SELECT '{m}' AS MGS, sample, base, cb,
             CASE WHEN frac IS NULL THEN -1 ELSE least(floor(frac * 20), 19)::INTEGER END AS fbin,
             CASE WHEN frac IS NULL THEN 0 ELSE ({level_sql(wilson('+'))}) - ({level_sql(wilson('-'))}) + 1 END AS ci_span,
             count(*)::BIGINT AS n
      FROM v WHERE clean = 1 AND cat = 'called' GROUP BY ALL) TO {out('frac')};
COPY (SELECT '{m}' AS MGS, base, cc, least(floor(frac * 20), 19)::INTEGER AS fbin,
             least(floor(score / 5), 50)::INTEGER * 5 AS sbin, count(*)::BIGINT AS n
      FROM v WHERE clean = 1 AND cat = 'called' AND frac IS NOT NULL GROUP BY ALL) TO {out('joint')};
COPY (SELECT '{m}' AS MGS, base, cat, cc, least(floor(score / 5), 50)::INTEGER * 5 AS sbin, count(*)::BIGINT AS n
      FROM v WHERE clean = 1 GROUP BY ALL) TO {out('score')};
CREATE TABLE site AS SELECT contig, pos, strand, base, max(clean) AS c,
       sum((clean = 1 AND coverage >= {a.min_cov})::INTEGER) AS g FROM v GROUP BY ALL;
COPY (SELECT '{m}' AS MGS, base, count(*)::BIGINT AS sites_any, count(*) FILTER (WHERE c = 1)::BIGINT AS sites_clean_any
      FROM site GROUP BY ALL) TO {out('sites')};
CREATE TABLE ctg AS SELECT contig, max(pos)::BIGINT AS length, (row_number() OVER (ORDER BY contig))::INTEGER AS cid
      FROM site GROUP BY contig;
COPY (SELECT '{m}' AS MGS, contig, length FROM ctg) TO {out('contigs')};
COPY (SELECT contig, cid FROM ctg) TO {out('ctg')};
COPY (SELECT contig, ((pos - 1) // {BLOCK})::INTEGER AS blk, sum(g)::BIGINT AS rows FROM site GROUP BY ALL
      HAVING sum(g) > 0) TO {out('blocks')};
DROP TABLE site;
CREATE TABLE folds AS SELECT sample, (row_number() OVER (ORDER BY sample) % 2)::UTINYINT AS fold
      FROM (SELECT DISTINCT sample FROM read_parquet('{sq(os.path.join(dst, 'cube.parquet'))}'));
CREATE TABLE ref AS
SELECT contig, pos, strand, base, fold,
       CASE WHEN avg((cat = 'called')::INTEGER) >= 0.9 THEN 'methylated' ELSE 'unmethylated' END AS kind
FROM v JOIN folds USING (sample) WHERE clean = 1 AND coverage >= {a.ref_cov}
GROUP BY contig, pos, strand, base, fold
HAVING count(*) >= {a.ref_n} AND (avg((cat = 'called')::INTEGER) >= 0.9
       OR (sum((cat <> 'none')::INTEGER) = 0 AND hash(contig, pos, strand) % 20 = 0));
COPY (SELECT '{m}' AS MGS, t.sample, t.base, r.kind, t.cb, count(*)::BIGINT AS n,
             sum((t.cat = 'called')::INTEGER)::BIGINT AS n_called, sum((t.cat IN ('possible', 'other'))::INTEGER)::BIGINT AS n_possible
      FROM v t JOIN folds f USING (sample)
             JOIN ref r ON r.contig = t.contig AND r.pos = t.pos AND r.strand = t.strand AND r.base = t.base AND r.fold <> f.fold
      WHERE t.clean = 1 GROUP BY ALL) TO {out('power')};
"""


def chunks_of(dst, chunk_rows):
    """genome blocks (contig order) -> chunks of <= chunk_rows graded rows (a block is never split) ; writes chunks.parquet,
    returns [(chunk, scan filter)] - the filter only lets DuckDB skip row groups (MMT files are sorted by contig, pos)"""
    import pyarrow as pa, pyarrow.parquet as pq
    b = pq.read_table(os.path.join(dst, "blocks.parquet")).to_pandas().sort_values(["contig", "blk"])
    k, cur, ch = 0, 0, []
    for r in b["rows"].values:
        if cur and cur + r > chunk_rows:
            k, cur = k + 1, 0
        ch.append(k)
        cur += r
    b["chunk"] = ch
    pq.write_table(pa.Table.from_pandas(b[["contig", "blk", "chunk"]], preserve_index=False),
                   os.path.join(dst, "chunks.parquet"))
    out = []
    for k, g in b.groupby("chunk"):
        lo, hi = g["contig"].iloc[0], g["contig"].iloc[-1]
        f = (f"v.contig = '{sq(lo)}' AND v.pos BETWEEN {g['blk'].min() * BLOCK + 1} AND {(g['blk'].max() + 1) * BLOCK}"
             if lo == hi else f"v.contig >= '{sq(lo)}' AND v.contig <= '{sq(hi)}'")
        out.append((int(k), f))
    return out or [(0, "false")]


def pass2_sql(mgs, src, dst, a, patients, chunks):
    """pass 2: per patient x site the best-covered graded sample, one genome chunk at a time (bounded memory) ->
    pr/<chunk>.parquet (integer keys) ; then per-site counts, chance rates, the agreement subset, spatial windows"""
    out = lambda name: f"'{sq(os.path.join(dst, name + '.parquet'))}' (FORMAT parquet, COMPRESSION zstd)"
    m, mc = sq(mgs), a.min_cov
    base = "CASE bc WHEN 1 THEN 'C' ELSE 'A' END"
    s = [sql_head(dst, a), view_sql(src, a), f"""
CREATE TABLE smp AS SELECT sample, rank, (dense_rank() OVER (ORDER BY patient))::INTEGER AS pid
      FROM read_csv('{sq(patients)}', delim = '\t', header = true,
                    columns = {{'sample': 'VARCHAR', 'patient': 'VARCHAR', 'rank': 'INTEGER'}});
CREATE TABLE ctg AS SELECT * FROM read_parquet('{sq(os.path.join(dst, 'ctg.parquet'))}');
CREATE TABLE bm AS SELECT * FROM read_parquet('{sq(os.path.join(dst, 'chunks.parquet'))}');
"""]
    for k, f in chunks:
        s.append(f"""
COPY (SELECT s.pid, c.cid, v.pos::INTEGER AS pos, (v.strand = '-')::UTINYINT AS sd, (v.base = 'C')::UTINYINT AS bc,
             arg_max(coalesce(v.cls, 9), v.coverage::BIGINT * 8192 + (8191 - s.rank))::UTINYINT AS code,
             max(v.coverage)::INTEGER AS cov
      FROM v JOIN smp s USING (sample) JOIN ctg c ON c.contig = v.contig
             JOIN bm b ON b.contig = v.contig AND b.blk = (v.pos - 1) // {BLOCK}
      WHERE b.chunk = {k} AND v.clean = 1 AND v.coverage >= {mc} AND {f}
      GROUP BY ALL) TO {out(f'pr/{k}')};""")
    s.append(f"""
CREATE VIEW pr AS SELECT *, ({case_bins('cov', strata(mc))})::USMALLINT AS st
      FROM read_parquet('{sq(os.path.join(dst, 'pr', '*.parquet'))}');
COPY (SELECT '{m}' AS MGS, {base} AS base, count(DISTINCT pid)::BIGINT AS patients FROM pr WHERE code <> 9 GROUP BY ALL)
      TO {out('npat')};
CREATE TABLE pi AS SELECT pid, bc, st, avg((code = 4)::INTEGER) AS p4, avg((code = 0)::INTEGER) AS p0
      FROM pr WHERE code <> 9 GROUP BY ALL;
CREATE TABLE sc AS
SELECT r.bc, r.cid, r.pos, r.sd,
       count(*) FILTER (WHERE r.code <> 9)::INTEGER AS n, count(*) FILTER (WHERE r.code = 9)::INTEGER AS n_possible,
       {", ".join(f"count(*) FILTER (WHERE r.code = {k})::INTEGER AS k{k}" for k in range(5))},
       sum(pi.p4) FILTER (WHERE r.code <> 9) AS e4, sum(pi.p0) FILTER (WHERE r.code <> 9) AS e0,
       max(r.cov) FILTER (WHERE r.code <> 9) AS max_cov
FROM pr r LEFT JOIN pi USING (pid, bc, st)
GROUP BY r.bc, r.cid, r.pos, r.sd HAVING count(*) FILTER (WHERE r.code <> 9) >= 2;
COPY (SELECT {base} AS base, c.contig, s.pos, CASE s.sd WHEN 1 THEN '-' ELSE '+' END AS strand,
             s.* EXCLUDE (bc, cid, pos, sd)
      FROM sc s JOIN ctg c USING (cid)) TO {out('site_counts')};
CREATE TABLE sub AS
SELECT s.bc, s.cid, s.pos, s.sd FROM sc s JOIN (SELECT bc, count(*) AS nb FROM sc GROUP BY bc) c USING (bc)
WHERE hash(s.cid, s.pos, s.sd) % 1000000 < 1000000.0 * {a.alpha_sites} / c.nb;
DROP TABLE sc;
COPY (SELECT {base} AS base, r.cid AS contig, r.pos, r.sd AS strand, r.pid AS patient, r.code, r.st
      FROM pr r JOIN sub USING (bc, cid, pos, sd) WHERE r.code <> 9) TO {out('alpha_sub')};
COPY (WITH s AS (SELECT bc, cid, pos, avg((code >= 1)::INTEGER) AS f, max(cov) AS mx
                 FROM pr WHERE code <> 9 GROUP BY bc, cid, pos, sd),
           s2 AS (SELECT bc, cid, pos, (f >= 0.5)::INTEGER AS k, {case_bins('mx', strata(mc))} AS st FROM s),
           ps AS (SELECT bc, st, avg(k) AS p FROM s2 GROUP BY ALL),
           w AS (SELECT bc, cid, floor((pos - 1) / {a.window})::INTEGER AS win, count(*) AS n, sum(k) AS k,
                        sum(p) AS e, sum(p * (1 - p)) AS v
                 FROM s2 JOIN ps USING (bc, st) GROUP BY ALL HAVING count(*) >= 10)
      SELECT '{m}' AS MGS, {base} AS base, count(*)::BIGINT AS windows, sum(n)::BIGINT AS sites, sum(k)::BIGINT AS methylated,
             sum(CASE WHEN v > 0 THEN (k - e) * (k - e) / v END) AS d
      FROM w GROUP BY ALL) TO {out('spatial')};
""")
    return "".join(s)


def post_mgs(mgs, dst, a):
    """per-site constitutive tests, alpha + its permutation null, summaries for the plots"""
    import numpy as np, pandas as pd, pyarrow.parquet as pq
    from scipy.stats import binom
    sc = pq.read_table(os.path.join(dst, "site_counts.parquet")).to_pandas()
    rows_a, summ, tracks, keep = [], [], [], []
    sub = pq.read_table(os.path.join(dst, "alpha_sub.parquet")).to_pandas()
    sub = sub.sort_values(["base", "contig", "pos", "strand", "patient"], ignore_index=True)   # DuckDB row order varies: same input -> same p
    for b in "AC":
        g = sc[sc["base"] == b].copy()
        tested = g["n"] >= a.min_patients
        for k, e, lab in ((4, "e4", "meth"), (0, "e0", "unmeth")):
            p = np.full(len(g), np.nan)
            t = tested.values
            p[t] = binom.sf(g[f"k{k}"].values[t] - 1, g["n"].values[t], np.clip(g[e].values[t] / g["n"].values[t], 0, 1))
            g["p_" + lab] = p
            g["q_" + lab] = bh(p)
        g["label"] = np.where(~tested, "few_patients",
                     np.where((g["q_meth"] < 0.05) & (g["k4"] >= 0.8 * g["n"]), "always_methylated",
                     np.where((g["q_unmeth"] < 0.05) & (g["k0"] >= 0.8 * g["n"]), "always_unmethylated", "other")))
        keep.append(g)
        cnt = g["label"].value_counts()
        summ.append(dict(MGS=mgs, base=b, sites_2plus=len(g), tested=int(tested.sum()),
                         **{x: int(cnt.get(x, 0)) for x in ("always_methylated", "always_unmethylated", "other", "few_patients")}))
        if len(g):
            tr = g[tested].assign(win=(g["pos"] - 1) // a.window)
            tr = tr.groupby(["contig", "win", "label"]).size().unstack(fill_value=0).reset_index()
            tracks.append(tr.assign(MGS=mgs, base=b))
        # agreement: alpha on every site, permutation null on the sampled sites
        r = dict(MGS=mgs, base=b, alpha_sites=int((g["n"] >= 2).sum()),
                 alpha=kalpha(g[[f"k{k}" for k in range(5)]].values) if len(g) else np.nan)
        s = sub[sub["base"] == b]
        if s["contig"].nunique() and len(s) > 1:
            site = s.groupby(["contig", "pos", "strand"], sort=False).ngroup().values
            grp = s.groupby(["patient", "st"], sort=False).ngroup().values
            code = s["code"].values.astype(int)
            S = site.max() + 1
            alpha_of = lambda c: kalpha(np.bincount(site * 5 + c, minlength=S * 5).reshape(S, 5))
            obs = alpha_of(code)
            rng = np.random.default_rng(int(hashlib.md5(f"{mgs}|{b}".encode()).hexdigest()[:8], 16))
            o0 = np.argsort(grp, kind="stable")
            null = []
            for _ in range(a.perm):
                cp = np.empty_like(code)
                cp[o0] = code[np.lexsort((rng.random(len(code)), grp))]
                null.append(alpha_of(cp))
            null = np.array(null)
            r.update(alpha_sub_sites=S, alpha_sub=obs, alpha_null_mean=np.nanmean(null), alpha_null_sd=np.nanstd(null),
                     p_alpha=(1 + np.sum(null >= obs)) / (1 + len(null)) if not np.isnan(obs) else np.nan)
        rows_a.append(r)
    allg = pd.concat(keep, ignore_index=True).assign(MGS=mgs)
    allg.to_parquet(os.path.join(dst, "site_classes.parquet"), index=False)
    allg[allg["label"].isin(["always_methylated", "always_unmethylated"])][
        ["MGS", "base", "contig", "pos", "strand", "n", "k4", "k0", "q_meth", "q_unmeth", "label"]].to_parquet(
        os.path.join(dst, "constitutive.parquet"), index=False)
    pd.DataFrame(rows_a).to_parquet(os.path.join(dst, "alpha.parquet"), index=False)
    pd.DataFrame(summ).to_parquet(os.path.join(dst, "site_summary.parquet"), index=False)
    tr = pd.concat(tracks, ignore_index=True) if tracks else pd.DataFrame(columns=["MGS", "base", "contig", "win"])
    for x in ("always_methylated", "always_unmethylated", "other"):
        if x not in tr:
            tr[x] = 0
    tr.to_parquet(os.path.join(dst, "track.parquet"), index=False)
    for x in ("site_counts", "alpha_sub"):
        os.remove(os.path.join(dst, x + ".parquet"))


VAR_VERSION = "1"
SETS = ("in-batch", "across-batch")
ALLSETS = SETS + ("between-patient",)


def read_run(path, groups):
    """--run FILE: lines 'sample<TAB>in-batch|across-batch' -> pairs [(set, group, sample_a, sample_b)]: within each
    AssmblGrps of the map, every two listed samples with the same label (sample_a < sample_b)"""
    lab, bad = {}, []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            f = line.rstrip("\r\n").split("\t")
            if not f[0].strip() or f[0].startswith("#"):
                continue
            if len(f) < 2 or f[1].strip() not in SETS:
                bad.append(f"line {n}: '{line.strip()}' - want: sample TAB in-batch|across-batch")
                continue
            s = f[0].strip()
            if s in lab:
                bad.append(f"line {n}: {s} listed twice")
            elif s not in groups:
                bad.append(f"line {n}: {s} is not in the map (or has no AssmblGrps)")
            else:
                lab[s] = f[1].strip()
    by = collections.defaultdict(list)
    for s, x in lab.items():
        by[(groups[s], x)].append(s)
    bad += [f"{v[0]} is the only {x} sample of AssmblGrps {g} - nothing to compare it with"
            for (g, x), v in sorted(by.items()) if len(v) < 2]
    if bad:
        die(f"--run {path}: {len(bad)} problem(s):\n  " + "\n  ".join(bad[:30]))
    return sorted((x, g, s1, s2) for (g, x), v in by.items() for s1, s2 in itertools.combinations(sorted(v), 2))


def between_pairs(mgs, have, pm, nrows, n):
    """--between N: N disjoint pairs of different patients of this MGS -> [('between-patient', 'pa|pb', s1, s2)]"""
    best = {}
    for s in sorted(have):   # best-covered sample per patient (ties: name order)
        p = pm.get(s, s)
        if p not in best or nrows.get((mgs, s), 0) > nrows.get((mgs, best[p]), 0):
            best[p] = s
    top = sorted(best, key=lambda p: (-nrows.get((mgs, best[p]), 0), p))[:2 * n]
    random.Random(hashlib.md5(mgs.encode()).hexdigest()).shuffle(top)
    return sorted(("between-patient", "|".join(sorted((pa, pb))), *sorted((best[pa], best[pb])))
                  for pa, pb in zip(top[0::2], top[1::2]))


def variance_sql(mgs, d, dst, a, pairs):
    """per pair of samples: sites graded (clean, >= min-cov reads) in both -> class of each x lower coverage stratum,
    with the sum and sum of squares of frac(a) - frac(b) (no call = 0 ; 'possible' / called without frac = class 9, no
    frac)"""
    m, mc = sq(mgs), a.min_cov
    g = lambda name: (f"(SELECT contig, pos, strand, base, coverage, coalesce(cls, 9) AS code, CASE WHEN cls IS NULL THEN NULL "
                      f"WHEN cat = 'none' THEN 0 ELSE frac::DOUBLE END AS val FROM {name} WHERE clean = 1 AND coverage >= {mc})")
    s = [sql_head(dst, a), """
CREATE TABLE res (set VARCHAR, grp VARCHAR, sample_a VARCHAR, sample_b VARCHAR, base VARCHAR, st USMALLINT, code_a UTINYINT,
                  code_b UTINYINT, n BIGINT, d_sum DOUBLE, d2_sum DOUBLE);"""]
    for x, grp, s1, s2 in pairs:
        s.append(view_sql(os.path.join(d, s1 + ".parquet"), a, "xa") + view_sql(os.path.join(d, s2 + ".parquet"), a, "xb"))
        s.append(f"""
INSERT INTO res SELECT '{x}', '{sq(grp)}', '{sq(s1)}', '{sq(s2)}', x.base, ({case_bins('least(x.coverage, y.coverage)', strata(mc))}),
       x.code, y.code, count(*), sum(x.val - y.val), sum((x.val - y.val) * (x.val - y.val))
FROM {g('xa')} x JOIN {g('xb')} y USING (contig, pos, strand, base) GROUP BY ALL;""")
    s.append(f"\nCOPY (SELECT '{m}' AS MGS, * FROM res) TO '{sq(os.path.join(dst, 'variance.part'))}' "
             f"(FORMAT parquet, COMPRESSION zstd);\n")
    return "".join(s)


def run_variance(mgs, d, dst, a, pairs, vfp):
    """the pairs of this MGS -> dst/variance.parquet + variance.fp ; -> error text or None"""
    vt = os.path.join(dst, "var.tmp")
    shutil.rmtree(vt, ignore_errors=True)
    os.makedirs(vt)
    err, _ = duck(a, os.path.join(vt, "v.duckdb"), variance_sql(mgs, d, vt, a, pairs), vt, "variance")
    if err:
        return err + f" (files kept in {vt})"
    os.replace(os.path.join(vt, "variance.part"), os.path.join(dst, "variance.parquet"))
    with open(os.path.join(dst, "variance.fp"), "w") as fh:
        fh.write(vfp + "\n")
    shutil.rmtree(vt)
    return None


def drop_variance(dst):
    for x in ("variance.parquet", "variance.fp"):
        if os.path.exists(os.path.join(dst, x)):
            os.remove(os.path.join(dst, x))


def job(mgs, d, work, a, fp, main, pairs, vfp):
    """main summaries (when changed), then the variance pairs of this MGS (when changed or after the main redo)"""
    msg = ""
    if main:
        mgs, msg = run_mgs(mgs, d, work, a, fp)
        if not msg.startswith("done"):
            return mgs, msg
    dst = os.path.join(work, "parts", mgs)
    if not pairs:
        drop_variance(dst)
        return mgs, msg
    vf = os.path.join(dst, "variance.fp")
    if main or not (os.path.isfile(vf) and open(vf).read().strip() == vfp):
        err = run_variance(mgs, d, dst, a, pairs, vfp)
        if err:
            return mgs, f"variance failed: {err}"
        msg = (msg + " ; " if msg else "done: ") + f"variance {len(pairs)} pair{'s' * (len(pairs) > 1)}"
    return mgs, msg


def mgs_files(d):
    fs = sorted(glob.glob(os.path.join(d, "*.parquet")))
    return fs, hashlib.md5(json.dumps([[os.path.basename(f), os.path.getsize(f), os.stat(f).st_mtime_ns] for f in fs])
                           .encode()).hexdigest()


def duck(a, db, sql, tmp, tag):
    """one DuckDB CLI script -> (error text or None, peak RSS in GB of that DuckDB process)"""
    f = os.path.join(tmp, tag + ".sql")
    with open(f, "w") as o:
        o.write(sql)
    peak = 0   # DuckDB's own high-water mark (/proc VmHWM ; rusage would include python's memory at the spawn)
    with open(f) as i, open(f + ".log", "w") as e:
        p = subprocess.Popen([a.duckdb, "-bail", db], stdin=i, stdout=e, stderr=subprocess.STDOUT)
        while p.poll() is None:
            try:
                with open(f"/proc/{p.pid}/status") as s:
                    peak = max(peak, next(int(x.split()[1]) for x in s if x.startswith("VmHWM")))
            except (OSError, StopIteration, ValueError):
                pass
            time.sleep(0.2)
    rc = p.returncode
    peak /= 2**20
    if rc == 0:
        return None, peak
    why = open(f + ".log").read().strip()[-600:]
    if rc < 0 or rc == 137:
        sig = -rc if rc < 0 else 9
        why = (f"killed by signal {sig}" + (" (9 = usually out of memory: fewer --jobs, or a larger session)" if sig == 9
               else "") + (f": {why}" if why else ""))
    return f"duckdb {tag} failed (exit {rc}, peak {peak:.1f} GB): {why or 'no message'}", peak


def run_mgs(mgs, d, work, a, fp):
    dst = os.path.join(work, "parts", mgs)
    tmp = dst + ".tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(os.path.join(tmp, "pr"))
    t0 = time.time()
    src = os.path.join(d, "*.parquet")
    err, pk1 = duck(a, os.path.join(tmp, "pass1.duckdb"), pass1_sql(mgs, src, tmp, a), tmp, "pass1")
    t1 = time.time() - t0
    if err:
        return mgs, err + f" (files kept in {tmp})"
    chunks = chunks_of(tmp, a.chunk_rows or a.mem_gb * 2_000_000)
    err, pk2 = duck(a, os.path.join(tmp, "pass2.duckdb"),
                    pass2_sql(mgs, src, tmp, a, os.path.join(work, "patients.tsv"), chunks), tmp, "pass2")
    if err:
        return mgs, err + f" (files kept in {tmp})"
    for x in os.listdir(tmp):
        if x.startswith(("pass1", "pass2", "duck.tmp", "blocks", "chunks", "ctg.")) or x == "pr":
            x = os.path.join(tmp, x)
            shutil.rmtree(x) if os.path.isdir(x) else os.remove(x)
    try:
        post_mgs(mgs, tmp, a)
    except Exception as e:   # one MGS failing must not stop the others
        return mgs, f"site tests failed: {type(e).__name__}: {e} (files kept in {tmp})"
    with open(os.path.join(tmp, "fp"), "w") as fh:
        fh.write(fp + "\n")
    shutil.rmtree(dst, ignore_errors=True)
    os.rename(tmp, dst)
    return mgs, (f"done in {time.time() - t0:.0f} s (pass 1 {t1:.0f} s), {len(chunks)} chunk{'s' * (len(chunks) > 1)}, "
                 f"DuckDB peak {max(pk1, pk2):.2f} GB")


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


def map_patients(path):
    """MATAFILER map -> SmplID -> AssmblGrps (the column named in the #SmplID header line)"""
    pat, col = {}, None
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").rstrip("\r").split("\t")
            if f[0] == "#SmplID":
                col = f.index("AssmblGrps") if "AssmblGrps" in f else None
            elif f[0] and not f[0].startswith("#") and col is not None and col < len(f) and f[col]:
                pat[f[0]] = f[col]
    if col is None:
        die(f"no '#SmplID ... AssmblGrps' header line in {path}")
    return pat


def summarise(a):
    root = os.path.join(a.out, "mgs-rep")
    if not os.path.isdir(root):
        die(f"no mgs-rep results in {a.out} (expected {root}/<MGS>/<sample>.parquet)")
    if shutil.which(a.duckdb) is None and not os.access(a.duckdb, os.X_OK):
        die(f"DuckDB CLI not found: {a.duckdb} (--duckdb, or MMT_DUCKDB)")
    work = a.work
    os.makedirs(os.path.join(work, "parts"), exist_ok=True)

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

    # ---- patients: AssmblGrps of the map ; samples not in the map are their own patient
    names = sorted({r[0] for r in samples} | {r[0] for r in units} |
                   {os.path.basename(f)[:-8] for f in glob.glob(os.path.join(root, "*", "*.parquet"))})
    pm = map_patients(a.map) if a.map else {}
    if not a.map:
        log("WARNING: no --map: every sample is its own patient (time points / replicates of a person then count twice)")
    miss = [s for s in names if s not in pm]
    if a.map and miss:
        log(f"WARNING: {len(miss)} samples not in {a.map} are their own patient: {' '.join(miss[:10])}")
    with open(os.path.join(work, "patients.tsv"), "w") as o:
        o.write("sample\tpatient\trank\n")
        o.writelines(f"{s}\t{pm.get(s, s)}\t{i}\n" for i, s in enumerate(names))
    log(f"{len(set(pm.get(s, s) for s in names))} patients for {len(names)} samples")
    vpairs, vfile = [], os.path.join(work, "variance_pairs.tsv")
    if a.run:
        if not a.map:
            die("--run needs --map (the assembly groups)")
        vpairs = read_run(a.run, pm)
        with open(vfile, "w") as o:
            o.write("set\tgroup\tsample_a\tsample_b\n")
            o.writelines("\t".join(p) + "\n" for p in vpairs)
        log(f"--run: {sum(p[0] == SETS[0] for p in vpairs)} in-batch and {sum(p[0] == SETS[1] for p in vpairs)} "
            f"across-batch pairs in {len({p[1] for p in vpairs})} assembly groups -> {vfile}")
    elif os.path.exists(vfile):
        os.remove(vfile)
    if a.between and not a.map:
        die("--between needs --map (one sample per patient)")
    nrows = {(r[1], r[0]): int(r[5]) for r in units if r[5].isdigit()}
    pmd5 = hashlib.md5(open(os.path.join(work, "patients.tsv"), "rb").read()).hexdigest()
    settings = dict(version=VERSION, min_cov=a.min_cov, min_ref_pct=a.min_ref_pct, max_ctx_pct=a.max_ctx_pct,
                    ref_cov=a.ref_cov, ref_n=a.ref_n, window=a.window, alpha_sites=a.alpha_sites, perm=a.perm,
                    min_patients=a.min_patients, strata=strata(a.min_cov), levels=LEVELS, patients_md5=pmd5)
    with open(os.path.join(work, "settings.json"), "w") as o:
        json.dump(settings, o, indent=1)
    sfp = hashlib.md5(json.dumps(settings, sort_keys=True).encode()).hexdigest()

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
        rows, nmiss = [], 0
        for mgs, mag in sorted(rep.items()):
            fa = os.path.join(gdir, rep_fasta_name(mgs, mag))
            if not os.path.isfile(fa):
                nmiss += 1
                continue
            st = f"{os.path.getsize(fa)}:{os.stat(fa).st_mtime_ns}"
            if (mgs, mag, st) in old:
                rows.append(old[(mgs, mag, st)]); continue
            na, nc, L = genome_counts(fa)
            rows.append([mgs, mag, str(na), str(nc), str(L), fa, st])
        with open(gf, "w") as o:
            o.write("MGS\trep_MAG\tA_sites\tC_sites\tlength\tfasta\tstamp\n")
            o.writelines("\t".join(r) + "\n" for r in rows)
        log(f"representative genomes: {len(rows)} counted, {nmiss} FASTA missing in {gdir}")

    # ---- per MGS: skip unchanged, drop gone, run the rest (largest first: shortest tail)
    todo, keepm, nbetween = [], set(), 0
    for d in sorted(glob.glob(os.path.join(root, "*", ""))):
        mgs = os.path.basename(os.path.dirname(d))
        if mgs == "stats" or (a.select and mgs not in a.select):
            continue
        fs, fp = mgs_files(d)
        if not fs:
            continue
        keepm.add(mgs)
        have = {os.path.basename(f)[:-8] for f in fs}
        pairs = [p for p in vpairs if p[2] in have and p[3] in have]
        if a.between:
            pairs += between_pairs(mgs, have, pm, nrows, a.between)
            nbetween += sum(p[0] == ALLSETS[2] for p in pairs)
        vfp = hashlib.md5(json.dumps([VAR_VERSION, fp, pairs, a.min_cov, a.min_ref_pct, a.max_ctx_pct, LEVELS])
                          .encode()).hexdigest()
        fp = f"v{VERSION}|{sfp}|{fp}"
        pf, vf = os.path.join(work, "parts", mgs, "fp"), os.path.join(work, "parts", mgs, "variance.fp")
        main = not (os.path.isfile(pf) and open(pf).read().strip() == fp)
        var = bool(pairs) and not (os.path.isfile(vf) and open(vf).read().strip() == vfp)
        if not main and not var:
            if not pairs:
                drop_variance(os.path.join(work, "parts", mgs))
            continue
        todo.append((sum(os.path.getsize(f) for f in fs) if main else 0, mgs, d, fp, main, pairs, vfp))
    if not a.select:
        for d in glob.glob(os.path.join(work, "parts", "*")):
            if os.path.basename(d) not in keepm:
                shutil.rmtree(d, ignore_errors=True)
                log(f"dropped {os.path.basename(d)} (no longer in {root})")
    if a.between:
        log(f"--between {a.between}: {nbetween} between-patient pairs over {len(keepm)} MGS")
    todo.sort(key=lambda t: (t[0], t[1]), reverse=True)
    nm = sum(t[4] for t in todo)
    log(f"{len(keepm)} MGS with results ; {nm} to summarise ({sum(t[0] for t in todo) / 2**30:.1f} GB), "
        f"{len(todo) - nm} variance only, {len(keepm) - len(todo)} unchanged ; {a.jobs} at once x {a.threads} threads "
        f"x {a.mem_gb} GB")
    fails = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(job, mgs, d, work, a, fp, main, pairs, vfp) for _, mgs, d, fp, main, pairs, vfp in todo]
        for i, fu in enumerate(cf.as_completed(futs), 1):
            mgs, msg = fu.result()
            log(f"{i}/{len(todo)} {mgs}: {msg}")
            if not msg.startswith("done"):
                fails.append(mgs)
    if fails:
        die(f"{len(fails)} MGS failed: {' '.join(fails[:20])} - fix and rerun (finished MGS are kept)")
    import resource
    log(f"summaries complete in {work} (python peak {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.1f} GB)")


# ============================================================================ plot: data
def load(work):
    import pandas as pd, pyarrow.parquet as pq
    parts = {}
    for x in PARTS:
        fs = sorted(glob.glob(os.path.join(work, "parts", "*", x + ".parquet")))
        if not fs:
            die(f"no {x}.parquet in {work}/parts - run summarise first (with this version)")
        parts[x] = pd.concat([pq.read_table(f).to_pandas() for f in fs], ignore_index=True)
    units = pd.read_csv(os.path.join(work, "units.tsv"), sep="\t", dtype=str)
    samples = pd.read_csv(os.path.join(work, "samples.tsv"), sep="\t", dtype=str)
    gf = os.path.join(work, "genome.tsv")
    genome = pd.read_csv(gf, sep="\t", dtype={"MGS": str}) if os.path.isfile(gf) else None
    settings = json.load(open(os.path.join(work, "settings.json")))
    return parts, units, samples, genome, settings


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


def analyse(parts, units, samples, genome, st_, a):
    """-> dict of tables every figure and card is drawn from"""
    import numpy as np, pandas as pd
    from scipy.stats import chi2, wilcoxon
    mc = st_["min_cov"]
    cube, frac, power = parts["cube"], parts["frac"], parts["power"]
    chem = dict(zip(samples["sample"], samples["chemistry"]))
    power = power.assign(chemistry=power["sample"].map(chem).fillna("NA"))

    # ---- n95 (diagnostic): per base, pooled over MGS and samples, also per chemistry
    n95, rows = {}, []
    for b in "AC":
        pm = power[(power["base"] == b) & (power["kind"] == "methylated")]
        n95[b] = n95_from_power(pm, a.power_min_n)
        rows.append([b, "all", n95[b], int(pm["n"].sum())])
        for c_, g in pm.groupby("chemistry"):
            rows.append([b, c_, n95_from_power(g, a.power_min_n), int(g["n"].sum())])
    n95_tab = pd.DataFrame(rows, columns=["base", "chemistry", "n95", "reference_rows"])

    # ---- per MGS x sample x base: tiers, classes, coverage
    c = cube.assign(cbv=cube["cb"].astype(int), possible=cube["cat"].isin(["possible", "other"]))
    key = ["MGS", "sample", "base"]
    cl = c[c["clean"] == 1]
    gr = cl[cl["cbv"] >= mc]
    agg = lambda m: cl[m].groupby(key)["n"].sum()
    d = {"sites": c.groupby(key)["n"].sum(), "masked": c[c["clean"] == 0].groupby(key)["n"].sum(),
         "clean": cl.groupby(key)["n"].sum(), "cov_sum": cl.groupby(key)["cov_sum"].sum(),
         "tier_low": agg(cl["cbv"] < mc), "tier_mid": agg((cl["cbv"] >= mc) & (cl["cbv"] < 20)),
         "tier_high": agg(cl["cbv"] >= 20), "not_graded": agg(cl["cbv"] < mc),
         "possible": gr[gr["possible"]].groupby(key)["n"].sum(),
         "called_no_frac": gr[(gr["cat"] == "called") & gr["cls"].isna()].groupby(key)["n"].sum()}
    for k in range(5):
        d[f"cls{k}"] = gr[gr["cls"] == k].groupby(key)["n"].sum()
    u = pd.DataFrame(d).fillna(0).reset_index()
    s = cl.groupby(key + ["cbv"])["n"].sum().reset_index().sort_values(key + ["cbv"])
    s["cum"] = s.groupby(key)["n"].cumsum()
    s["tot"] = s.groupby(key)["n"].transform("sum")
    u = u.merge(s[s["cum"] >= s["tot"] / 2].groupby(key)["cbv"].first().rename("cov_median"), on=key, how="left")
    nz = lambda x: x.where(x > 0)
    u["cov_mean"] = u["cov_sum"] / nz(u["clean"])
    u["graded"] = sum(u[f"cls{k}"] for k in range(5)) + u["possible"] + u["called_no_frac"]
    u["methylated"] = sum(u[f"cls{k}"] for k in range(1, 5))
    for k in range(5):                                            # composition of every clean site
        u[f"f_cls{k}"] = u[f"cls{k}"] / nz(u["clean"])
    u["f_possible"] = (u["possible"] + u["called_no_frac"]) / nz(u["clean"])
    u["f_not_graded"] = u["not_graded"] / nz(u["clean"])
    u["g_methylated"] = u["methylated"] / nz(u["graded"])          # among graded sites (>= min_cov reads)
    u["g_possible"] = (u["possible"] + u["called_no_frac"]) / nz(u["graded"])
    u["g_unmethylated"] = u["cls0"] / nz(u["graded"])
    for k in range(1, 5):                                         # level among methylated graded sites
        u[f"lvl{k}"] = u[f"cls{k}"] / nz(u["methylated"])
    for t in ("low", "mid", "high"):
        u["f_tier_" + t] = u["tier_" + t] / nz(u["clean"])
    u["f_masked"] = u["masked"] / nz(u["sites"])

    # ---- denominators and usability (A+C together)
    if genome is not None and len(genome):
        den = pd.concat([genome[["MGS"]].assign(base="A", den=genome["A_sites"]),
                         genome[["MGS"]].assign(base="C", den=genome["C_sites"])])
        den_src = "representative genome"
    else:
        den = parts["sites"].rename(columns={"sites_clean_any": "den"})[["MGS", "base", "den"]]
        den_src = "union of clean sites over samples"
    u = u.merge(den, on=["MGS", "base"], how="left")
    u["breadth"] = u["graded"] / u["den"]
    tot = u.groupby(["MGS", "sample"]).agg(g_AC=("graded", "sum"), den_AC=("den", "sum")).reset_index()
    tot["usable"] = (tot["g_AC"] >= a.min_sites) & (tot["g_AC"] / tot["den_AC"] >= a.min_breadth)
    u = u.merge(tot[["MGS", "sample", "usable"]], on=["MGS", "sample"])
    u["chemistry"] = u["sample"].map(chem).fillna("NA")

    # ---- coverage skew per usable unit: logit(call rate) ~ log2(coverage) over graded bins with >= 50 sites, WLS
    cc = pd.DataFrame({"n": gr.groupby(key + ["cbv"])["n"].sum(),
                       "k": gr[gr["cat"] == "called"].groupby(key + ["cbv"])["n"].sum()}).fillna(0).reset_index()
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
        sl.append(list(k_) + [b1, math.sqrt(disp / sxx), len(g)])
    slope = pd.DataFrame(sl, columns=key + ["slope", "se", "bins"])
    if len(slope):
        from scipy.stats import norm
        slope["p"] = 2 * norm.sf(np.abs(slope["slope"] / slope["se"]))
        slope["q"] = bh(slope["p"].values)
        slope = slope.merge(u[key + ["usable"]], on=key)
        slope = slope[slope["usable"]].drop(columns="usable")

    # ---- per MGS (means over usable samples ; coverage over all samples with a table)
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
             "not_graded_pct": 100 * g["not_graded"].sum() / max(g["clean"].sum(), 1),
             "samples_not_graded_only": int((g["graded"] == 0).sum())}
        mean = lambda col: 100 * gu[col].mean() if len(gu) else np.nan
        for k in range(5):
            r[f"{CLASSES[k]}_pct"] = mean(f"f_cls{k}")
        r["possible_pct"], r["not_graded_pct_usable"] = mean("f_possible"), mean("f_not_graded")
        r["methylated_graded_pct"], r["possible_graded_pct"] = mean("g_methylated"), mean("g_possible")
        r["unmethylated_graded_pct"] = mean("g_unmethylated")
        r["either_graded_pct"] = r["methylated_graded_pct"] + r["possible_graded_pct"]
        gl = gu[gu["methylated"] > 0]
        for k in range(1, 5):
            r[f"level_{CLASSES[k]}_pct"] = 100 * gl[f"lvl{k}"].mean() if len(gl) else np.nan
        sg = slope[(slope["MGS"] == mgs) & (slope["base"] == b)] if len(slope) else slope
        r["skew_slope_median"] = sg["slope"].median() if len(sg) else np.nan
        r["skew_samples_tested"] = len(sg)
        r["skew_samples_q05"] = int((sg["q"] < 0.05).sum()) if len(sg) else 0
        rows.append(r)
    mg = pd.DataFrame(rows).merge(pres, on="MGS", how="left")
    if a.taxonomy:
        tx = pd.read_csv(a.taxonomy, sep="\t", header=None, names=["MGS", "taxonomy"], dtype=str, comment="#")
        mg = mg.merge(tx, on="MGS", how="left")

    # ---- non-random: agreement (alpha), constitutive sites, spatial clustering ; q across MGS per base
    st = parts["alpha"].merge(parts["npat"], on=["MGS", "base"], how="outer")
    st = st.merge(parts["site_summary"], on=["MGS", "base"], how="outer")
    sp = parts["spatial"].rename(columns={"windows": "spatial_windows", "sites": "spatial_sites",
                                          "methylated": "spatial_methylated", "d": "spatial_d"})
    sp = sp[(sp["spatial_windows"] >= 3) & sp["spatial_d"].notna()].copy()
    sp["spatial_df"] = sp["spatial_windows"] - 1
    sp["log2_dispersion"] = np.log2(sp["spatial_d"] / sp["spatial_df"])
    sp["p_spatial"] = chi2.sf(sp["spatial_d"], sp["spatial_df"])
    st = st.merge(sp, on=["MGS", "base"], how="outer")
    st.insert(2, "mod", st["base"].map(MOD))
    for x in ("always_methylated", "always_unmethylated", "other"):
        st[x + "_pct"] = 100 * st[x] / st["tested"].where(st["tested"] > 0)
    for t_, p_ in (("alpha", "p_alpha"), ("spatial", "p_spatial")):
        st["q_" + t_] = np.nan
        for b in "AC":
            k_ = (st["base"] == b) & st[p_].notna()
            st.loc[k_, "q_" + t_] = bh(st.loc[k_, p_].values)
    st["alpha_excess"] = st["alpha_sub"] - st["alpha_null_mean"]
    rows = []
    for t_, eff, qq, null in (("agreement", "alpha_excess", "q_alpha", "alpha - its permutation null = 0"),
                              ("spatial", "log2_dispersion", "q_spatial", "log2 dispersion = 0 (random along the genome)"),
                              ("constitutive", "always_methylated_pct", None, "share of tested sites always methylated")):
        for b in "AC":
            dd = st[(st["base"] == b) & st[eff].notna()]
            e = dd[eff].values
            pw_ = wilcoxon(e, alternative="greater").pvalue if t_ != "constitutive" and len(e) >= 5 and np.any(e != 0) else np.nan
            rows.append(dict(test=t_, base=b, mod=MOD[b], effect=eff, MGS_tested=len(dd),
                             MGS_significant=int((dd[qq] < 0.05).sum()) if qq else np.nan,
                             pct_significant=100 * (dd[qq] < 0.05).mean() if qq and len(dd) else np.nan,
                             effect_median=np.median(e) if len(e) else np.nan,
                             effect_q25=np.percentile(e, 25) if len(e) else np.nan,
                             effect_q75=np.percentile(e, 75) if len(e) else np.nan, p_across_MGS=pw_, null=null,
                             test_across_MGS="one-sided Wilcoxon signed-rank, MGS as units" if t_ != "constitutive" else ""))
    stats_across = pd.DataFrame(rows)
    mg = mg.merge(st[["MGS", "base", "alpha", "q_alpha", "always_methylated_pct", "always_unmethylated_pct",
                      "log2_dispersion", "q_spatial"]], on=["MGS", "base"], how="left")
    order = mg.groupby("MGS")["samples_usable"].max().sort_values(ascending=False)
    return dict(n95=n95, n95_tab=n95_tab, units=u, mgs=mg, slope=slope, power=power, den_src=den_src, order=list(order.index),
                frac=frac.assign(cbv=frac["cb"].astype(int)), joint=parts["joint"], score=parts["score"], cube=c,
                stats=st, stats_across=stats_across, track=parts["track"], contigs=parts["contigs"],
                constitutive=parts["constitutive"], min_cov=mc, min_sites=a.min_sites, min_breadth=a.min_breadth)


# ============================================================================ plot: style
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
BASE_COL = {"A": "#2a78d6", "C": "#d95926"}
BASE_LIGHT = {"A": "#86b6ef", "C": "#f0a070"}
BASE_DARK = {"A": "#104281", "C": "#8a3412"}
RAMP = {"A": ["#86b6ef", "#3987e5", "#184f95", "#0b2c57"], "C": ["#f0a070", "#d65a1f", "#963a13", "#4a1c08"]}  # low..near-complete
POSSIBLE, GREY_UNMETH, GREY_UNDET = "#6250d6", "#8f8e88", "#e2e1dc"
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


def headroom(ax, f=0.3):
    """free band above the data for annotations"""
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + (hi - lo) * f)


def coarse_of(cbv):
    import numpy as np
    return np.searchsorted(np.array(COARSE), np.asarray(cbv), side="right") - 1


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
    mg, mc = R["mgs"], R["min_cov"]
    fig = plt.figure(figsize=(7.1, 7.6))
    gs = fig.add_gridspec(4, 2, height_ratios=[0.8, 1, 1, 0.8])
    m0 = mg[mg["base"] == "A"].set_index("MGS").reindex(R["order"]).reset_index()
    m0 = m0[m0["samples_with_table"] > 0]
    x0 = np.arange(len(m0))
    ax = fig.add_subplot(gs[0, :])
    ax.bar(x0, m0["samples_with_table"], width=1.0, color=GREY_UNDET, label="samples with a table", linewidth=0)
    ax.bar(x0, m0["samples_usable"], width=1.0, color=INK2, label=f"usable (>= {R['min_sites']:,} graded A+C sites, "
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
        ax.axhline(mc, color=INK, lw=0.7, ls="--", label=f"{mc} (graded from here)")
        ax.set_ylabel("median coverage per sample\n(clean sites, reads / strand)")
        title(ax, "bc"[j], f"{b} sites ({MOD[b]}): coverage")
        if j == 1:
            ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), title="per MGS, over samples", title_fontsize=6)
        rank_axis(ax, len(m), "usable samples", label=False)
        ax = fig.add_subplot(gs[2, j])
        t = R["units"][R["units"]["base"] == b].groupby("MGS")[["tier_low", "tier_mid", "tier_high", "clean"]].sum()
        t = t.reindex(m["MGS"]).fillna(0)
        tt = t[["tier_low", "tier_mid", "tier_high"]].div(t["clean"].where(t["clean"] > 0), axis=0).fillna(0) * 100
        stack_bars(ax, x, [tt["tier_high"], tt["tier_mid"], tt["tier_low"]], TIER[::-1],
                   [">= 20 reads", f"{mc}-19 reads", f"3-{mc - 1} reads (not graded)"])
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
    mc = R["min_cov"]
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
            ax.axvline(mc, color=INK, lw=0.6, ls="--")
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
        ax.axvline(mc, color=INK, lw=0.6, ls="--")
        if R["n95"][b]:
            ax.text(0.97, 0.06, f"95 % detection from {R['n95'][b]} reads", transform=ax.transAxes, ha="right",
                    fontsize=5.5, color=INK)
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
        title(ax, "abcdefgh"[i * 4 + 3], f"{b}: coverage skew, >= {mc} reads")
        tabs[f"{b}_skew"] = s
    fig.tight_layout(w_pad=0.6, h_pad=1.0)
    save(fig, report, "fig2_coverage_effects", tabs)


def fig3(R, report):
    """what fraction of clean sites is methylated, at which level: class composition per MGS + values across MGS"""
    import numpy as np, matplotlib.pyplot as plt
    mg, mc = R["mgs"], R["min_cov"]
    fig = plt.figure(figsize=(7.1, 6.8))
    gs = fig.add_gridspec(4, 2, width_ratios=[3.0, 1.25], height_ratios=[1.5, 0.75, 1.5, 0.75])
    tabs = {}
    for i, b in enumerate("AC"):
        m = mg[(mg["base"] == b) & (mg["samples_usable"] > 0)].copy()
        m["meth_pct"] = sum(m[f"{CLASSES[k]}_pct"] for k in range(1, 5))
        m = m.sort_values(["meth_pct", "possible_pct"], ascending=False)
        x = np.arange(len(m))
        cols = RAMP[b][::-1] + [POSSIBLE, GREY_UNMETH, GREY_UNDET]
        labs = [CLASS_LAB[k] for k in (4, 3, 2, 1)] + ["possible (modified_base)", CLASS_LAB[0], f"not graded (< {mc} reads)"]
        ax = fig.add_subplot(gs[2 * i, 0])
        stack_bars(ax, x, [m[f"{CLASSES[k]}_pct"] for k in (4, 3, 2, 1)] + [m["possible_pct"], m["unmethylated_pct"],
                                                                             m["not_graded_pct_usable"]], cols, labs)
        ax.set_ylim(0, 100)
        ax.set_ylabel(f"% of clean {b} sites\n(mean over usable samples)")
        title(ax, "ac"[i], f"{b} sites: methylation class of every clean site, per MGS")
        rank_axis(ax, len(m), "% methylated", label=False)
        hd, lb = ax.get_legend_handles_labels()
        ax2 = fig.add_subplot(gs[2 * i + 1, 0], sharex=ax)
        stack_bars(ax2, x, [m[f"{CLASSES[k]}_pct"] for k in (4, 3, 2, 1)] + [m["possible_pct"]], cols[:5], [""] * 5)
        top = np.nanmax((m["meth_pct"] + m["possible_pct"]).values) if len(m) else 1
        ax2.set_ylim(0, max(top * 1.1, 0.5))
        ax2.set_ylabel("zoom: methylated\n+ possible")
        rank_axis(ax2, len(m), "% methylated")
        ax2.legend(hd, lb, loc="upper center", bbox_to_anchor=(0.5, 0), borderaxespad=4.6, ncol=3, fontsize=5.5, handlelength=1.2,
                   columnspacing=1.2)
        ax = fig.add_subplot(gs[2 * i:2 * i + 2, 1])
        violins(ax, [m["methylated_graded_pct"], m["possible_graded_pct"], m["either_graded_pct"]],
                [BASE_COL[b], POSSIBLE, BASE_DARK[b]], ["methylated\n(>= 5 %)", "possible", "either"])
        ax.set_ylim(bottom=0)
        headroom(ax, 0.2)
        ax.set_ylabel(f"% of graded {b} sites (>= {mc} reads)")
        title(ax, "bd"[i], f"{b}: graded sites")
        ax.text(0.98, 0.98, f"n = {m['methylated_graded_pct'].notna().sum()} MGS\nunmethylated = 100 - either",
                transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=INK2)
        tabs[f"{b}_per_MGS"] = m[["MGS", "samples_usable"] + [f"{c}_pct" for c in CLASSES] +
                                 ["possible_pct", "not_graded_pct_usable", "methylated_graded_pct", "possible_graded_pct",
                                  "unmethylated_graded_pct", "either_graded_pct"]]
    fig.tight_layout(h_pad=0.4, w_pad=1.2)
    save(fig, report, "fig3_landscape", tabs)


def fig4(R, report):
    """methylation level of methylated sites per MGS, frac by coverage (pooled), near-complete share across MGS"""
    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    mg, f, mc = R["mgs"], R["frac"], R["min_cov"]
    fig = plt.figure(figsize=(7.1, 4.9))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.1, 1.55, 0.9])
    tabs = {}
    for i, b in enumerate("AC"):
        m = mg[(mg["base"] == b) & mg["level_near-complete_pct"].notna()].sort_values("level_near-complete_pct",
                                                                                         ascending=False)
        x = np.arange(len(m))
        ax = fig.add_subplot(gs[i, 0])
        stack_bars(ax, x, [m[f"level_{CLASSES[k]}_pct"] for k in (4, 3, 2, 1)], RAMP[b][::-1],
                   [CLASS_LAB[k] for k in (4, 3, 2, 1)])
        ax.set_ylim(0, 100)
        ax.set_ylabel(f"% of methylated {b} sites\n(>= {mc} reads, frac >= 5 %)")
        title(ax, "ad"[i], f"{b}: methylation level of {MOD[b]} sites, per MGS")
        rank_axis(ax, len(m), "% near-complete")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=5.5, handlelength=1.2)
        ax = fig.add_subplot(gs[i, 1])
        g = f[(f["base"] == b) & (f["fbin"] >= 0)].assign(cc=lambda d: coarse_of(d["cbv"]))
        g = g.groupby(["cc", "fbin"])["n"].sum().unstack(fill_value=0).reindex(columns=range(20), fill_value=0)
        g = g[g.sum(axis=1) >= 100]
        if len(g):
            z = (g.T / g.sum(axis=1)).values * 100
            cmap = LinearSegmentedColormap.from_list("b", ["#ffffff", BASE_LIGHT[b], BASE_COL[b], BASE_DARK[b]])
            im = ax.pcolormesh(np.arange(len(g) + 1), np.linspace(0, 1, 21), z, cmap=cmap, shading="flat", rasterized=True)
            for e in LEVELS:
                ax.axhline(e, color=INK2, lw=0.4, ls=":")
            ax.set_xticks(np.arange(len(g)) + 0.5)
            ax.set_xticklabels([COARSE_LAB[k] for k in g.index], rotation=45, ha="right")
            cb = fig.colorbar(im, ax=ax, pad=0.03, fraction=0.07)
            cb.set_label("% of the column", fontsize=6)
            cb.ax.tick_params(labelsize=5.5)
            tabs[f"{b}_frac_by_coverage"] = g.rename(index=lambda k: COARSE_LAB[k]).reset_index()
        ax.grid(False)
        ax.set_xlabel("coverage (reads / strand)")
        ax.set_ylabel(f"frac ({MOD[b]}); dotted: class edges")
        title(ax, "be"[i], f"{b}: frac by coverage, all MGS")
        ax = fig.add_subplot(gs[i, 2])
        violins(ax, [m["level_near-complete_pct"]], [BASE_DARK[b]], ["near-complete\n(>= 95 %)"])
        headroom(ax, 0.25)
        ax.set_ylabel(f"% of methylated {b} sites per MGS")
        title(ax, "cf"[i], "across MGS")
        if len(m):
            ax.text(0.97, 0.97, f"median {np.median(m['level_near-complete_pct']):.0f} %\nn = {len(m)} MGS",
                    transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=INK2)
        tabs[f"{b}_per_MGS"] = m[["MGS", "samples_usable"] + [f"level_{CLASSES[k]}_pct" for k in range(1, 5)]]
    fig.tight_layout(w_pad=1.0, h_pad=1.0)
    save(fig, report, "fig4_levels", tabs)


def fig5(R, report):
    """confidence: frac x score (pooled), score by state, how many classes the frac 95 % interval spans by coverage"""
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
        for st_, col, lab in (("none", GREY_UNMETH, "no call"), ("possible", POSSIBLE, "modified_base"),
                              ("called", BASE_COL[b], f"called {MOD[b]}")):
            g = s[s["state"] == st_].groupby("sbin")["n"].sum()
            if g.sum() == 0:
                continue
            g = g.reindex(range(0, 255, 5), fill_value=0)
            n_ = g.sum()
            ax.step(g.index, 100 * g / n_, where="post", color=col, lw=1.1,
                    label=f"{lab} ({n_ / 1e6:.1f} M)" if n_ >= 1e6 else f"{lab} ({n_ / 1e3:.0f} k)" if n_ >= 1e3 else f"{lab} ({n_})")
        logy(ax)
        ax.set_xlabel("score")
        ax.set_ylabel("% of the state's clean sites")
        title(ax, "be"[i], f"{b}: score by state")
        ax.set_ylim(top=ax.get_ylim()[1] * 30)
        ax.legend(fontsize=5.2, loc="upper right")
        tabs[f"{b}_score_by_state"] = s.groupby(["state", "cc", "sbin"])["n"].sum().reset_index()
        ax = axs[i, 2]
        f = R["frac"].merge(usable, on=["MGS", "sample", "base"])
        f = f[(f["base"] == b) & (f["fbin"] >= 0)].assign(cc=lambda d: coarse_of(d["cbv"]),
                                                         span=lambda d: d["ci_span"].clip(upper=3))
        g = f.groupby(["cc", "span"])["n"].sum().unstack(fill_value=0).reindex(columns=[1, 2, 3], fill_value=0)
        g = g[g.sum(axis=1) >= 100]
        if len(g):
            p = g.div(g.sum(axis=1), axis=0) * 100
            xs = np.arange(len(p))
            bot = np.zeros(len(p))
            for col, c_, lab in ((1, BASE_DARK[b], "within one class"), (2, BASE_COL[b], "spans 2 classes"),
                                 (3, GREY_UNDET, "spans >= 3 classes")):
                ax.bar(xs, p[col], bottom=bot, width=0.85, color=c_, label=lab, edgecolor="white", linewidth=0.5)
                bot += p[col].values
            ax.set_xticks(xs)
            ax.set_xticklabels([COARSE_LAB[k] for k in p.index], rotation=45, ha="right")
            tabs[f"{b}_ci_span_by_coverage"] = g.rename(index=lambda k: COARSE_LAB[k]).reset_index()
            ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=5.2, title="frac 95 % CI", title_fontsize=5.5)
        ax.set_ylim(0, 100)
        ax.grid(axis="x", visible=False)
        ax.set_xlabel("coverage (reads / strand)")
        ax.set_ylabel(f"% of called {MOD[b]} sites with frac")
        title(ax, "cf"[i], f"{b}: is the class resolvable?")
    fig.tight_layout(w_pad=0.8, h_pad=1.2)
    save(fig, report, "fig5_confidence", tabs)


def fig6(R, report):
    """non-random methylation: agreement across patients, constitutive sites, spatial clustering"""
    import numpy as np, matplotlib.pyplot as plt
    st, ac = R["stats"], R["stats_across"]
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.6), gridspec_kw={"width_ratios": [1.3, 1.5, 1.3]})
    for i, b in enumerate("AC"):
        s = st[st["base"] == b]
        get = lambda t_: ac[(ac["base"] == b) & (ac["test"] == t_)].iloc[0] if len(ac[(ac["base"] == b) & (ac["test"] == t_)]) else None
        ax = axs[i, 0]
        d = s[s["alpha_sub"].notna()].sort_values("alpha_sub", ascending=False)
        x = np.arange(len(d))
        ax.scatter(x, d["alpha_null_mean"], s=7, color=GREY_UNDET, linewidths=0, label="shuffled null")
        sig = (d["q_alpha"] < 0.05).values
        ax.scatter(x[~sig], d["alpha_sub"].values[~sig], s=9, color=GREY_UNMETH, linewidths=0, label="q >= 0.05")
        ax.scatter(x[sig], d["alpha_sub"].values[sig], s=9, color=BASE_COL[b], linewidths=0, label="q < 0.05")
        ax.axhline(0, color=INK2, lw=0.6)
        rank_axis(ax, len(d), "alpha")
        lo = min(-0.1, np.nanmin(d["alpha_sub"]) - 0.05) if len(d) else -0.1
        ax.set_ylim(lo, 1 + (1 - lo) * 0.42)                   # band above 1 (alpha <= 1) for the notes
        ax.set_yticks([t for t in np.arange(-1, 1.01, 0.2) if t >= lo])
        ax.set_ylabel(f"Krippendorff's alpha ({MOD[b]},\n5 classes, 1 = same class everywhere)")
        title(ax, "ad"[i], f"{b}: agreement across patients")
        r = get("agreement")
        if r is not None:
            ax.text(0.98, 0.97, f"{int(r['MGS_significant'])} / {int(r['MGS_tested'])} MGS q < 0.05\n"
                                f"median excess {r['effect_median']:+.2f}\nacross MGS p = {r['p_across_MGS']:.1e}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=INK2)
        ax.legend(loc="upper left", fontsize=5.2, markerscale=1.2, handletextpad=0.2, borderaxespad=0.2)
        ax = axs[i, 1]
        d = s[s["tested"] > 0].sort_values(["always_methylated_pct", "always_unmethylated_pct"], ascending=False)
        x = np.arange(len(d))
        stack_bars(ax, x, [d["always_methylated_pct"], d["always_unmethylated_pct"]], [RAMP[b][3], GREY_UNMETH],
                   ["always near-complete methylated", "always near-complete unmethylated"])
        ax.set_ylim(bottom=0)
        headroom(ax, 0.3)
        ax.set_ylabel(f"% of tested sites (>= {R['min_patients']} patients;\nthe rest: not constitutive)")
        rank_axis(ax, len(d), "% always methylated")
        title(ax, "be"[i], f"{b}: constitutive sites")
        ax.legend(loc="upper right", fontsize=5.2)
        ax = axs[i, 2]
        d = s[s["log2_dispersion"].notna()].sort_values("log2_dispersion", ascending=False)
        x = np.arange(len(d))
        sig = (d["q_spatial"] < 0.05).values
        ax.axhline(0, color=INK2, lw=0.6)
        ax.scatter(x[~sig], d["log2_dispersion"].values[~sig], s=9, color=GREY_UNMETH, linewidths=0, label="q >= 0.05")
        ax.scatter(x[sig], d["log2_dispersion"].values[sig], s=9, color=BASE_COL[b], linewidths=0, label="q < 0.05")
        rank_axis(ax, len(d), "log2 dispersion")
        ax.set_ylabel("spatial log2 dispersion\n(0 = random, coverage-adjusted)")
        title(ax, "cf"[i], f"{b}: clustered along the genome?")
        r = get("spatial")
        if r is not None:
            ax.text(0.98, 0.97, f"{int(r['MGS_significant'])} / {int(r['MGS_tested'])} MGS q < 0.05\n"
                                f"median {r['effect_median']:+.2f}", transform=ax.transAxes, ha="right", va="top",
                    fontsize=5.5, color=INK2)
        ax.legend(loc="best", fontsize=5.2, markerscale=1.2)
    fig.tight_layout(w_pad=0.8, h_pad=1.2)
    save(fig, report, "fig6_non_random", {"per_MGS": st, "across_MGS": ac})


def cards(R, report, a):
    """one page per MGS: rows A / C (coverage, class per sample, calls vs coverage, frac) + constitutive-site track"""
    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    u, mg, cube, f, st, mc = R["units"], R["mgs"], R["cube"], R["frac"], R["stats"], R["min_cov"]
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
    tr_by, ctg_by = dict(tuple(R["track"].groupby("MGS"))), dict(tuple(R["contigs"].groupby("MGS")))
    fmt = lambda v, s="{:.2f}": ("NA" if v is None or (isinstance(v, float) and np.isnan(v))
                                 else "< 1e-300" if s.endswith("e}") and v == 0 else s.format(v))
    path = os.path.join(report, "mgs_cards.pdf")
    with PdfPages(path) as pdf:
        for mgs in order:
            fig = plt.figure(figsize=(8.27, 7.0))
            gs = fig.add_gridspec(3, 4, width_ratios=[1.2, 1.5, 1.3, 1.1], height_ratios=[1, 1, 0.55])
            axs = np.array([[fig.add_subplot(gs[r, c]) for c in range(4)] for r in range(2)])
            info = mg[mg["MGS"] == mgs].set_index("base")
            sinfo = st[st["MGS"] == mgs].set_index("base")
            h = info.iloc[0]
            tax = h["taxonomy"] if "taxonomy" in info.columns and isinstance(h.get("taxonomy"), str) else ""
            lines = [f"{mgs}" + (f"   {tax}" if tax else ""),
                     f"samples with a table {int(h['samples_with_table'])}, usable {int(h['samples_usable'])}   |   median "
                     f"coverage {fmt(h['cov_median_of_samples'], '{:.0f}')}   |   graded from {mc} reads"]
            for b in "AC":
                if b in info.index:
                    r, sr = info.loc[b], (sinfo.loc[b] if b in sinfo.index else {})
                    lines.append(f"{b}: methylated {fmt(r['methylated_graded_pct'])} % of graded sites (possible "
                                 f"{fmt(r['possible_graded_pct'])} %)  |  near-complete {fmt(r['level_near-complete_pct'], '{:.0f}')} "
                                 f"% of methylated  |  alpha {fmt(sr.get('alpha_sub', np.nan))} (q {fmt(sr.get('q_alpha', np.nan), '{:.1e}')})"
                                 f"  |  always methylated {fmt(sr.get('always_methylated_pct', np.nan), '{:.1f}')} % of "
                                 f"{fmt(sr.get('tested', np.nan), '{:.0f}')} sites  |  spatial log2 "
                                 f"{fmt(sr.get('log2_dispersion', np.nan))} (q {fmt(sr.get('q_spatial', np.nan), '{:.1e}')})")
            fig.text(0.01, 0.995, lines[0], fontsize=8, fontweight="bold", va="top")
            fig.text(0.01, 0.972, "\n".join(lines[1:]), fontsize=5.6, va="top", color=INK2, linespacing=1.5)
            fm, cm = frac_by.get(mgs), cu_by.get(mgs)
            for i, b in enumerate("AC"):
                g = u[(u["MGS"] == mgs) & (u["base"] == b)].sort_values("cov_median")
                y = np.arange(len(g))
                ax = axs[i, 0]
                ax.barh(y, g["cov_median"], color=[BASE_COL[b] if v else GREY_UNDET for v in g["usable"]], height=0.8, linewidth=0)
                ax.axvline(mc, color=INK, lw=0.6, ls="--")
                logx(ax, 1, 1000)
                if len(g) <= 40:   # sample names only while they stay legible
                    ax.set_yticks(y)
                    ax.set_yticklabels(g["sample"], fontsize=4.5)
                else:
                    ax.set_yticks([])
                    ax.set_ylabel(f"{len(g)} samples")
                ax.set_xlabel("median coverage (grey: not usable)")
                ax.grid(axis="y", visible=False)
                title(ax, "", f"{b}: coverage per sample")
                ax = axs[i, 1]
                left = np.zeros(len(g))
                cols = RAMP[b][::-1] + [POSSIBLE, GREY_UNMETH, GREY_UNDET]
                parts_ = [g[f"f_cls{k}"] for k in (4, 3, 2, 1)] + [g["f_possible"], g["f_cls0"], g["f_not_graded"]]
                labs = ["near-complete", "high", "medium", "low", "possible", "unmethylated", "not graded"]
                for v, col, lab in zip(parts_, cols, labs):
                    v = 100 * v.fillna(0).values
                    ax.barh(y, v, left=left, color=col, height=0.8, linewidth=0, label=lab)
                    left += v
                ax.set_xlim(0, 100)
                ax.set_yticks([])
                ax.set_xlabel(f"% of clean {b} sites")
                ax.grid(axis="y", visible=False)
                title(ax, "", f"{b}: class per sample (same order)")
                if i == 1:
                    ax.legend(fontsize=4.8, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=4, handlelength=1,
                              columnspacing=0.6)
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
                        ax.plot(dp["cbv"], 100 * dp["k"] / dp["n"], "s--", color=POSSIBLE, ms=1.8, lw=0.8, label="possible")
                ax.axvline(mc, color=INK, lw=0.6, ls="--")
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
                    hh = fb.groupby("fbin")["n"].sum().reindex(range(20), fill_value=0)
                    if hh.sum():
                        ax.step(np.arange(21) / 20, np.r_[100 * hh / hh.sum(), 100 * hh.iloc[-1] / hh.sum()], where="post",
                                color=BASE_COL[b], lw=1.0, label=f"n = {hh.sum():,}")
                        ax.legend(fontsize=5, loc="upper left")
                for e in LEVELS:
                    ax.axvline(e, color=INK2, lw=0.4, ls=":")
                ax.set_xlim(0, 1)
                ax.set_xlabel(f"frac ({MOD[b]}); dotted: class edges")
                ax.set_ylabel("% of called sites")
                title(ax, "", f"{b}: frac distribution")
            ax = fig.add_subplot(gs[2, :])     # constitutive sites along the genome (contigs end to end, name order)
            tr, ctg = tr_by.get(mgs), ctg_by.get(mgs)
            if tr is not None and ctg is not None and len(tr):
                ctg = ctg.sort_values("contig")
                off = dict(zip(ctg["contig"], np.r_[0, np.cumsum(ctg["length"].values)[:-1]]))
                sc_, un = (1e6, "Mb") if ctg["length"].sum() >= 2e6 else (1e3, "kb")
                for b, sign in (("A", 1), ("C", -1)):
                    t = tr[tr["base"] == b]
                    tot = (t["always_methylated"] + t["always_unmethylated"] + t["other"]).where(lambda v: v > 0)
                    xk = (t["contig"].map(off) + t["win"] * a.window + a.window / 2) / sc_
                    am, au = 100 * t["always_methylated"] / tot, 100 * t["always_unmethylated"] / tot
                    ax.bar(xk, sign * am, width=a.window / sc_, color=RAMP[b][3], linewidth=0, label=f"{b}: always methylated")
                    ax.bar(xk, sign * au, bottom=sign * am, width=a.window / sc_, color=GREY_UNMETH, linewidth=0,
                           label="always unmethylated" if b == "A" else None)
                for o_ in list(off.values())[1:]:
                    ax.axvline(o_ / sc_, color=INK2, lw=0.3)
                ax.axhline(0, color=INK, lw=0.5)
                ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):g}"))
                headroom(ax, 0.25)
                ax.legend(fontsize=5, loc="upper right", ncol=3)
                ax.set_xlabel(f"representative genome ({un}; contigs end to end, thin lines = contig ends; "
                              f"{a.window} bp windows)")
            ax.set_ylabel("% of tested sites\nA up / C down")
            ax.grid(axis="x", visible=False)
            title(ax, "", f"Constitutive sites along the genome (>= {R['min_patients']} patients)")
            fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=0.6, h_pad=0.9)
            pdf.savefig(fig)
            plt.close(fig)
    log(f"cards: {len(order)} MGS -> {path}")


def variance_tables(work, mc, a):
    """WORK/parts/*/variance.parquet -> dict of tables (summary per set x base, per coverage stratum, per pair, fig7:
    pairs at >= --pair-min-cov, per stratum, technical bands + test, evidence per MGS x base) or None"""
    import numpy as np, pandas as pd, pyarrow.parquet as pq
    from scipy.stats import wilcoxon
    fs = sorted(glob.glob(os.path.join(work, "parts", "*", "variance.parquet")))
    if not fs:
        return None
    v = pd.concat([pq.read_table(f).to_pandas() for f in fs], ignore_index=True)
    both = (v["code_a"] < 9) & (v["code_b"] < 9)
    disc = ((v["code_a"] == 0) & (v["code_b"].between(1, 4))) | ((v["code_b"] == 0) & (v["code_a"].between(1, 4)))
    v = v.assign(both=both, sites=v["n"].where(both, 0), same=v["n"].where(both & (v["code_a"] == v["code_b"]), 0),
                 mdisc=v["n"].where(both & disc, 0), union=v["n"].where(both & ((v["code_a"] > 0) | (v["code_b"] > 0)), 0),
                 possible=v["n"].where(~both, 0), d_sum=v["d_sum"].fillna(0), d2_sum=v["d2_sum"].fillna(0),
                 cell=v["code_a"] * 5 + v["code_b"])
    edges = strata(mc)
    lab = {e: (f"{e}-{edges[i + 1] - 1}" if i + 1 < len(edges) else f"{e}+") for i, e in enumerate(edges)}
    order = {s: i for i, s in enumerate(ALLSETS)}
    by_set = lambda df: df.sort_values(["set", "base"], key=lambda c: c.map(order) if c.name == "set" else c)

    def summ(df, keys):
        g = df.groupby(keys)[["sites", "same", "mdisc", "union", "possible", "d_sum", "d2_sum"]].sum()
        cells = df[df["both"]].groupby(keys + ["cell"])["n"].sum().unstack(fill_value=0).reindex(columns=range(25), fill_value=0)
        cells = cells.reindex(g.index, fill_value=0)
        alpha = [alpha_coinc(C + C.T) for C in cells.values.reshape(-1, 5, 5)]
        nz = g["sites"].where(g["sites"] > 0)
        out = pd.DataFrame({"sites": g["sites"], "sites_possible_either": g["possible"],
                            "same_class_pct": 100 * g["same"] / nz, "methylated_discordant_pct": 100 * g["mdisc"] / nz,
                            "methylated_either": g["union"], "switched_pct": 100 * g["mdisc"] / g["union"].where(g["union"] > 0),
                            "alpha": alpha, "frac_rms_diff": np.sqrt(g["d2_sum"] / nz)}, index=g.index)
        out["frac_sd_per_sample"] = out["frac_rms_diff"] / np.sqrt(2)
        return out, g

    pair_keys = ["MGS", "set", "grp", "sample_a", "sample_b", "base"]
    pairs, g = summ(v, pair_keys)
    pairs["frac_mean_diff_a_minus_b"] = g["d_sum"] / g["sites"].where(g["sites"] > 0)
    hi, _ = summ(v[v["st"] >= a.pair_min_cov], pair_keys)   # fig7: both samples >= --pair-min-cov reads
    hi = hi[["sites", "methylated_either", "switched_pct", "alpha"]].add_suffix("_hicov")
    pairs = pairs.join(hi).reset_index().rename(columns={"grp": "group"})
    pairs.insert(6, "mod", pairs["base"].map(MOD))
    cov, _ = summ(v, ["set", "base", "st"])
    cov = cov.reset_index()
    cov.insert(2, "mod", cov["base"].map(MOD))
    cov["coverage"] = cov["st"].map(lab)
    cov = cov.drop(columns="st")
    sm, _ = summ(v, ["set", "base"])
    sm = sm.reset_index()
    sm.insert(2, "mod", sm["base"].map(MOD))
    cnt = (pairs.assign(pair=pairs["sample_a"] + "~" + pairs["sample_b"]).groupby(["set", "base"])
           .agg(pairs=("pair", "nunique"), groups=("group", "nunique"), MGS=("MGS", "nunique")))
    ok = pairs[pairs["sites"] >= a.var_min_sites]
    med = ok.groupby(["set", "base"]).agg(pair_x_MGS_with_min_sites=("sites", "size"),
                                          median_same_class_pct=("same_class_pct", "median"),
                                          median_alpha=("alpha", "median"),
                                          median_frac_sd_per_sample=("frac_sd_per_sample", "median"))
    sm = by_set(sm.merge(cnt.reset_index(), on=["set", "base"]).merge(med.reset_index(), on=["set", "base"], how="left"))

    # ---- fig7: pair x MGS values (enough methylated sites at high coverage), per stratum, bands, evidence
    pv = pairs[pairs["methylated_either_hicov"] >= a.pair_min_meth]
    pst = v.groupby(pair_keys + ["st"])[["mdisc", "union"]].sum().reset_index()
    pst = pst[pst["union"] >= a.pair_min_meth].assign(switched_pct=lambda d: 100 * d["mdisc"] / d["union"])
    q = lambda p_: (lambda x: np.percentile(x, p_))
    ev = (pv.groupby(["MGS", "base", "set"])["switched_pct_hicov"].agg(["size", "median", q(25), q(75)])
          .set_axis(["pairs", "median", "q25", "q75"], axis=1).unstack("set"))
    ev.columns = [f"{x.replace('-', '_')}_{k}" for k, x in ev.columns]
    for x in ALLSETS:   # every set's columns, also when a set has no pairs
        for k in ("pairs", "median", "q25", "q75"):
            c = f"{x.replace('-', '_')}_{k}"
            if c not in ev.columns:
                ev[c] = 0 if k == "pairs" else np.nan
    ev = ev.reset_index()
    bands = []
    for b in "AC":
        row = {}
        for x in ALLSETS:
            z = pv[(pv["base"] == b) & (pv["set"] == x)]["switched_pct_hicov"]
            bands.append(dict(set=x, base=b, mod=MOD[b], pairs=len(z), MGS=pv[(pv["base"] == b) & (pv["set"] == x)]["MGS"].nunique(),
                              mean=z.mean(), sd=z.std(ddof=1), median=z.median()))
            row[x] = bands[-1]
        t = row["across-batch"]["mean"] + 2 * row["across-batch"]["sd"] if row["across-batch"]["pairs"] >= 2 else np.nan
        e = ev[ev["base"] == b]
        sel = e["between_patient_median"].notna() & e["across_batch_median"].notna()
        dif = (e.loc[sel, "between_patient_median"] - e.loc[sel, "across_batch_median"]).values
        p_ = wilcoxon(dif, alternative="greater").pvalue if len(dif) >= 5 and np.any(dif != 0) else np.nan
        tested = e["between_patient_median"].notna()
        above = tested & (e["between_patient_median"] > t)
        ev.loc[e.index, "above_band"] = above.where(tested)
        for r in bands[-3:]:
            r.update(threshold_mean_2sd=t, MGS_between_tested=int(tested.sum()), MGS_above_threshold=int(above.sum()),
                     MGS_paired=int(sel.sum()), MGS_between_gt_across=int((dif > 0).sum()), wilcoxon_p=p_)
    bands = by_set(pd.DataFrame(bands))
    ev.insert(2, "mod", ev["base"].map(MOD))
    ev["switched_pct_cov"] = a.pair_min_cov
    return dict(variance_summary=sm, variance_by_coverage=cov, variance_pairs=pairs, variance_evidence=ev,
                variance_bands=bands, _pv=pv, _pst=pst, _lab=lab)


SET_COL = {"in-batch": "#2e9b72", "across-batch": "#7a5bc4"}


def fig7(V, report, a):
    """methylome change per MGS between patients, against the technical noise of replicate pairs"""
    import numpy as np, matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ev, pv, pst, bd, lab = V["variance_evidence"], V["_pv"], V["_pst"], V["variance_bands"], V["_lab"]
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 5.2), gridspec_kw={"width_ratios": [2.2, 1.05, 1.25]})
    ylab = f"methylated sites switched (%)\n(>= 5 % in one, < 5 % in the other ;\nboth >= {a.pair_min_cov} reads)"
    for i, b in enumerate("AC"):
        bb = bd[bd["base"] == b].set_index("set")
        thr = bb["threshold_mean_2sd"].iloc[0]
        ax = axs[i, 0]
        e = ev[(ev["base"] == b) & ev["between_patient_median"].notna()].sort_values(["between_patient_median", "MGS"])
        x = np.arange(len(e))
        pos = dict(zip(e["MGS"], x))
        for s_ in SETS:   # technical band: mean +- 1 SD over replicate pairs x MGS
            if bb.loc[s_, "pairs"] >= 2:
                m_, sd = bb.loc[s_, "mean"], bb.loc[s_, "sd"]
                ax.axhspan(max(m_ - sd, 0), m_ + sd, color=SET_COL[s_], alpha=0.16, lw=0, zorder=0)
                ax.axhline(m_, color=SET_COL[s_], lw=0.7, zorder=1)
        if np.isfinite(thr):
            ax.axhline(thr, color=SET_COL["across-batch"], lw=0.8, ls="--", zorder=1)
        up = e["above_band"].fillna(False).astype(bool).values
        ax.vlines(x, e["between_patient_q25"], e["between_patient_q75"], color=GREY_UNDET, lw=0.6, zorder=2)
        ax.scatter(x, e["between_patient_median"], s=6, c=[BASE_COL[b] if u else GREY_UNMETH for u in up], lw=0, zorder=3)
        for s_, mk in zip(SETS, ("o", "D")):   # each replicate pair at its MGS
            r = pv[(pv["base"] == b) & (pv["set"] == s_) & pv["MGS"].isin(pos)]
            ax.scatter(r["MGS"].map(pos), r["switched_pct_hicov"], s=7, marker=mk, facecolors="none",
                       edgecolors=SET_COL[s_], lw=0.6, zorder=4)
        rank_axis(ax, len(e), "median between-patient change")
        ax.set_ylim(bottom=0)
        headroom(ax, 0.32)
        ax.set_ylabel(ylab)
        title(ax, "ad"[i], f"{b} ({MOD[b]}): between-patient change per MGS")
        r0 = bb.loc["between-patient"]
        if np.isfinite(thr):
            top3 = ", ".join(e[up]["MGS"].iloc[::-1].head(3))
            ax.text(0.02, 0.97, f"{int(r0['MGS_above_threshold'])} / {int(r0['MGS_between_tested'])} MGS: median above "
                                f"across-batch mean + 2 SD ({thr:.1f} %)" + (f"\nlargest: {top3}" if top3 else ""),
                    transform=ax.transAxes, ha="left", va="top", fontsize=5.5, color=INK2)
        ax = axs[i, 1]
        ax.sharey(axs[i, 0])
        cols = [SET_COL["in-batch"], SET_COL["across-batch"], BASE_COL[b]]
        violins(ax, [pv[(pv["base"] == b) & (pv["set"] == s_)]["switched_pct_hicov"] for s_ in ALLSETS], cols,
                ["in-\nbatch", "across-\nbatch", "between\npatients"])
        ax.tick_params(axis="x", labelsize=5.5)
        title(ax, "be"[i], f"{b}: all pairs x MGS")
        note = [f"{s_}: n = {int(bb.loc[s_, 'pairs'])}" for s_ in ALLSETS]
        if np.isfinite(r0["wilcoxon_p"]):
            note += [f"between > across-batch in {int(r0['MGS_between_gt_across'])} / {int(r0['MGS_paired'])} MGS",
                     f"paired Wilcoxon p = {r0['wilcoxon_p']:.1e}"]
        ax.text(0.98, 0.97, "\n".join(note), transform=ax.transAxes, ha="right", va="top", fontsize=5, color=INK2)
        ax = axs[i, 2]
        sts = sorted(lab)
        for k, (s_, c) in enumerate(zip(ALLSETS, cols)):
            z = pst[(pst["base"] == b) & (pst["set"] == s_)]
            m_ = [z[z["st"] == t_]["switched_pct"] for t_ in sts]
            xs = np.arange(len(sts)) + (k - 1) * 0.22
            md = [np.median(w) if len(w) else np.nan for w in m_]
            lo_ = [np.percentile(w, 25) if len(w) else np.nan for w in m_]
            hi_ = [np.percentile(w, 75) if len(w) else np.nan for w in m_]
            ax.vlines(xs, lo_, hi_, color=c, lw=0.9)
            ax.plot(xs, md, "o-", ms=2.8, color=c, lw=0.8, label=s_)
        ax.set_xticks(range(len(sts)))
        ax.set_xticklabels([lab[t_] for t_ in sts], rotation=30)
        ax.set_xlabel("reads (lower of the two samples)")
        ax.set_ylim(bottom=0)
        headroom(ax, 0.25)
        ax.set_ylabel("methylated sites switched (%)")
        ax.grid(axis="x", visible=False)
        title(ax, "cf"[i], f"{b}: by coverage (median, IQR)")
        ax.legend(loc="upper right", fontsize=5)
    fig.tight_layout(w_pad=0.8, h_pad=1.2, rect=(0, 0.07, 1, 1))
    fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=3, color=BASE_COL["A"], label="between patients: MGS median above threshold"),
                        Line2D([], [], marker="o", ls="", ms=3, color=GREY_UNMETH, label="MGS median within threshold (bar = IQR of its pairs)"),
                        Line2D([], [], color=SET_COL["across-batch"], ls="--", lw=0.8, label="threshold: across-batch mean + 2 SD"),
                        Patch(color=SET_COL["in-batch"], alpha=0.35, label="in-batch replicate pairs: mean +- 1 SD"),
                        Patch(color=SET_COL["across-batch"], alpha=0.35, label="across-batch replicate pairs: mean +- 1 SD"),
                        Line2D([], [], marker="o", ls="", ms=3, mfc="none", color=SET_COL["in-batch"], label="one in-batch pair"),
                        Line2D([], [], marker="D", ls="", ms=3, mfc="none", color=SET_COL["across-batch"], label="one across-batch pair")],
               loc="lower center", ncol=3, fontsize=5.5, handletextpad=0.3, columnspacing=1.2, bbox_to_anchor=(0.5, 0))
    save(fig, report, "fig7_variation")


def plot(a):
    style()
    parts, units, samples, genome, st_ = load(a.work)
    a.window, a.min_patients = st_["window"], st_["min_patients"]
    R = analyse(parts, units, samples, genome, st_, a)
    R["min_patients"] = st_["min_patients"]
    os.makedirs(os.path.join(a.report, "tables"), exist_ok=True)
    R["mgs"].to_csv(os.path.join(a.report, "tables", "per_MGS.tsv"), sep="\t", index=False)
    R["units"].to_csv(os.path.join(a.report, "tables", "per_MGS_x_sample.tsv"), sep="\t", index=False)
    R["n95_tab"].to_csv(os.path.join(a.report, "tables", "n95.tsv"), sep="\t", index=False)
    R["stats"].to_csv(os.path.join(a.report, "stats_per_MGS.tsv"), sep="\t", index=False)
    R["stats_across"].to_csv(os.path.join(a.report, "stats_across_MGS.tsv"), sep="\t", index=False)
    R["constitutive"].to_parquet(os.path.join(a.report, "constitutive_sites.parquet"), index=False)
    if a.pair_min_cov not in strata(st_["min_cov"]):
        die(f"--pair-min-cov {a.pair_min_cov}: must be a coverage stratum edge {strata(st_['min_cov'])}")
    vt = variance_tables(a.work, st_["min_cov"], a)
    if vt is not None:
        names = [k for k in vt if not k.startswith("_")]
        for name in names:
            vt[name].to_csv(os.path.join(a.report, name + ".tsv"), sep="\t", index=False)
        log(f"variance: {len(vt['variance_pairs'])} pair x MGS x base rows -> {', '.join(n + '.tsv' for n in names)}")
    log(f"graded from {R['min_cov']} reads ; n95 (diagnostic) A {R['n95']['A']}, C {R['n95']['C']} ; usable MGS x sample "
        f"{int(R['units']['usable'].sum() / 2)} ; breadth denominator: {R['den_src']}")
    for fn in (fig1, fig3, fig4, fig5, fig6):
        fn(R, a.report)
    fig2(R, a.report, a)
    if vt is None:
        log("fig7 skipped: no pair tables (summarise --run FILE --between N)")
    else:
        n = vt["variance_bands"].set_index(["set", "base"])["pairs"]
        if n.xs("between-patient").sum() and max(n.xs("in-batch").max(), n.xs("across-batch").max()) >= 2:
            fig7(vt, a.report, a)
        else:
            log("fig7 skipped: needs between-patient pairs (--between) and >= 2 replicate pairs (--run) with >= "
                f"{a.pair_min_meth} methylated sites at >= {a.pair_min_cov} reads")
    if a.cards != 0:
        cards(R, a.report, a)
    with open(os.path.join(a.report, "README.txt"), "w") as o:
        o.write(f"mmt_landscape.py v{VERSION} plot {time.strftime('%Y-%m-%d %H:%M')}\nWORK {os.path.abspath(a.work)}\n"
                f"settings {json.dumps(st_)}\nusable: >= {a.min_sites} graded A+C sites and breadth >= {a.min_breadth} "
                f"({R['den_src']})\n\n" + __doc__)
    log(f"report in {a.report}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("summarise", help="MMT results -> small summary tables (resumable)")
    s.add_argument("-o", "--out", required=True, help="MMT output folder (MMT -o): holds mgs-rep/")
    s.add_argument("-w", "--work", required=True, help="summary folder (created)")
    s.add_argument("--map", help="MATAFILER map: AssmblGrps = patient (strongly recommended)")
    s.add_argument("--run", help="variance pairs: lines 'sample TAB in-batch|across-batch' ; within each AssmblGrps, every "
                                 "two samples with the same label are compared (needs --map)")
    s.add_argument("--between", type=int, default=0, help="per MGS, this many pairs of different patients for the "
                                                          "biological comparison in fig7 (0 = none ; needs --map)")
    s.add_argument("--duckdb", default=os.environ.get("MMT_DUCKDB", "duckdb"), help="DuckDB CLI (default $MMT_DUCKDB or duckdb)")
    s.add_argument("--jobs", type=int, default=2, help="MGS at once (2)")
    s.add_argument("--threads", type=int, default=8, help="DuckDB threads per MGS (8)")
    s.add_argument("--mem-gb", type=int, default=24, help="memory budget per MGS; DuckDB gets 70 %%, spills to disk beyond "
                                                          "(24 ; jobs x mem-gb must fit the session)")
    s.add_argument("--chunk-rows", type=int, default=0, help="pass 2: graded rows per genome chunk (0 = 2 M x --mem-gb ; "
                                                             "results do not depend on it)")
    s.add_argument("--genecat", help="GeneCat folder: representative genome sizes for breadth (optional)")
    s.add_argument("--binner", default="SB", help="GeneCat Bin_<binner> (default SB)")
    s.add_argument("--select", type=lambda x: set(x.split(",")), help="only these MGS (comma-separated)")
    s.add_argument("--min-cov", type=int, default=11, help="graded from this many reads (11: frac exists above 10)")
    s.add_argument("--min-ref-pct", type=float, default=95, help="clean: %% of reads with the reference base (95)")
    s.add_argument("--max-ctx-pct", type=int, default=10, help="clean: ctx_alt_pct below this (10)")
    s.add_argument("--ref-cov", type=int, default=30, help="n95 diagnostic: coverage defining reference sites (30)")
    s.add_argument("--ref-n", type=int, default=3, help="n95 diagnostic: samples per half (3)")
    s.add_argument("--window", type=int, default=1000, help="spatial test + genome track: window (bp) (1000)")
    s.add_argument("--min-patients", type=int, default=3, help="constitutive-site test: patients per site (3)")
    s.add_argument("--alpha-sites", type=int, default=20000, help="agreement permutation: sampled sites per base (20000)")
    s.add_argument("--perm", type=int, default=199, help="agreement permutations (199: smallest p = 0.005)")
    q = sp.add_parser("plot", help="summary tables -> figures, cards, tables")
    q.add_argument("-w", "--work", required=True)
    q.add_argument("-r", "--report", required=True, help="figure / table folder (created)")
    q.add_argument("--taxonomy", help="TSV MGS<TAB>label (no header) for card titles and per_MGS.tsv")
    q.add_argument("--min-sites", type=int, default=20000, help="usable: graded A+C sites (20000)")
    q.add_argument("--min-breadth", type=float, default=0.10, help="usable: fraction of the denominator (0.10)")
    q.add_argument("--power-min-n", type=int, default=30, help="reference sites per coverage bin for the n95 curve (30)")
    q.add_argument("--cards", type=int, default=-1, help="cards: -1 all, 0 none, N the first N MGS")
    q.add_argument("--cards-min-usable", type=int, default=1, help="cards only for MGS with >= this many usable samples")
    q.add_argument("--var-min-sites", type=int, default=1000, help="variance medians: pairs x MGS with >= this many "
                                                                  "sites graded in both samples (1000)")
    q.add_argument("--pair-min-cov", type=int, default=30, help="fig7: sites with >= this many reads in both samples "
                                                               "(a stratum edge: 11, 15, 20, 30, 50 ; default 30)")
    q.add_argument("--pair-min-meth", type=int, default=30, help="fig7: pairs with >= this many methylated sites (either "
                                                                "sample) at that coverage (30)")
    a = p.parse_args()
    summarise(a) if a.cmd == "summarise" else plot(a)


if __name__ == "__main__":
    main()
