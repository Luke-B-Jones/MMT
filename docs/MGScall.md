# MGScall - methylation of each MGS per sample (on its representative MAG, or on each MAG itself)

`mmt MGScall (--mgs-rep | --mag-rep | --mag-mag) -o OUT [--select MGS,...] [--samples S,...] [-l LOGDIR] [--tmp DIR]`

For every sample and every selected MGS (GeneCat species bin), MGScall takes the sample's PacBio HiFi reads that
MATAFILER4 placed on that MGS, aligns them to the MGS's representative MAG (mag-mag: each MAG's reads to that MAG), and
calls m6A and m4C with ipdSummary from the read kinetics. Result: one table per MGS x sample (mag-rep, mag-mag: per
MAG), every A and C of the reference on both strands with >= 3 kinetic reads, as ipdSummary reports it - nothing is filtered or masked; the columns say what
each value can be trusted for.

## Modes (one per run; each has its own folder OUT/<mode>/)

| Mode | Reads used for an MGS in a sample | Tables |
|---|---|---|
| `--mgs-rep` | reads on ANY MAG of the MGS in that sample's assembly, pooled (the species as present in the sample) | `mgs-rep/<MGS>/<S>.parquet` |
| `--mag-rep` | reads on EACH MAG of the MGS in that sample's assembly, apart (strains apart) | `mag-rep/<MGS>/<S>/<MAG>.parquet` |
| `--mag-mag` | the same reads as mag-rep, each MAG's on that MAG's own contigs (the strain's own sequence; coordinates differ per MAG) | `mag-mag/<MGS>/<S>/<MAG>.parquet` |

All modes cover the same samples: every sample whose assembly has a MAG of the MGS (co-assembly: every sample of the
group). Where a sample has one MAG of the MGS, mgs-rep and mag-rep use the same reads and give the same table: the mode run
second hard-links the first one's file (status `linked`; a copy across file systems) instead of computing it again -
run mgs-rep first, then mag-rep only computes samples with several MAGs of an MGS.

## Method

1. **Layout** (index: a cache in the work space `<tmp>/mmt/<OUT name>_<key>/index/`, rebuilt when missing or when the
   map, `MAGvsGC.txt.gz` or GeneCat's abundance table changes). From the MATAFILER4 map, run folder and GeneCat
   `Bin_<binner>/MAGvsGC.txt.gz`, with GeneCat's own rules: representative = starred MAG of the MGS, else the first
   non-canopy MAG; MAG `<sample or alias>__<bin>` -> sample (co-assembly alias `<SmplID>M<n>`) -> assembly
   (`assemblies/metag/assembly.txt`) -> contigs (`Binning/<binner>/<SmplID>`). The representative's sequence is
   GeneCat's own `Genomes/MGS_ctg/<MGS>.ctgs.<MAG>[.fna].gz` - the MAG's contigs exactly as assembled and binned
   (MMT stops if they differ from the bin file). Contig names stay MATAFILER4's, so every row traces back to MF4.
2. **Presence.** A sample is processed for an MGS when its assembly has a MAG of the MGS (any reads; rows need >= 3).
   GeneCat's abundance table (`Bin_<binner>/Annotation/Abundance/MGS.matL0.txt`, MGS x samples) is read too: its
   value is in every stats row (`abundance`), and where it has the MGS but no MAG of it was binned in the sample, the
   mmt log (PLAN) and `mmt status` list the sample as "present but not binned" (no reads are assigned to it there).
3. **Read -> MAG** (cached per sample next to the index). Each primary alignment in MATAFILER4's
   `<S>/mapping/<S>-smd.cram` (minimap2 map-hifi to the sample's own assembly) whose contig belongs to an MGS MAG.
   Alignments MATAFILER4 flagged as failed (0x200: > 5 % edits, < 50 % of the read aligned or MAPQ < 30) are never
   used; duplicates (0x400) only with `MGSCALL_KEEP_DUP=1`. No new competitive mapping: MATAFILER4 already decided.
4. **Kinetics reads.** The selected reads are taken from the sample's raw kinetics BAM (map `#DirPath` +
   `SmplPrefix`.bam; one pass), split into one pseudo-read per strand (`ccs-kinetics-bystrandify`) and grouped by MGS
   (mag-rep, mag-mag: by MAG).
5. **Alignment** of each group to its MGS's representative only (mag-mag: to its MAG only - the MAG's contigs taken from
   the sample's assembly `scaffolds.fasta.filt[.gz]`, the same contigs GeneCat binned) (`pbmm2 --preset CCS -y MGSCALL_MIN_IDENTITY`).
   Kept: primary, MAPQ >= `MGSCALL_MIN_MAPQ` (drops repeats within the representative), no hard clip, aligned part
   >= `MGSCALL_MIN_PCT_ALIGNED` % of the read and >= `MGSCALL_MIN_ALN_BP` bp. Reads of a strain that differs from
   the representative are kept as long as they align at that identity; the sequence-agreement columns show where
   they differ.
6. **Kinetics** - ipdSummary (`--identify m6A,m4C --methylFraction --mapQvThreshold MGSCALL_MIN_MAPQ --minCoverage 1
   --maxCoverage 100000`, plus `MMT_IPD_ARGS` / `--ipd-args` as they are) on all MGS of the sample at once (mag-rep,
   mag-mag: one run per "layer" - the n-th MAG of every MGS in run n, so two MAGs never share a reference) and, from the same
   aligned reads, `samtools mpileup` for the sequence-agreement columns. Without extra arguments ipdSummary picks the
   kinetics model matching each BAM's chemistry; the model it loaded is recorded per sample (`stats/<S>.run.tsv`,
   `kinetics_model`) with the chemistry of the raw BAM (binding kit, sequencing kit, basecaller).
7. **Tables.** Every row ipdSummary writes for an A or C (on the row's strand), and for a G or T where it flags a
   kinetic signal (`call` = `modified_base`: other or unknown modifications, e.g. phosphorothioate), with kinetic
   coverage >= `MGSCALL_MIN_ROW_COV` (3 = every ipdSummary row: it writes none below 3 reads).

## Output

```
OUT/mmt_<date>.log                     the folder's log: every run, its plan per MGS, one line per sample with each
                                       MGS's (mag-rep: each MAG's) outcome, failures with their reason
OUT/mgs-rep/<MGS>/<S>.parquet          one table per MGS x sample (zstd), sorted by contig, pos, strand
OUT/mag-rep/<MGS>/<S>/<MAG>.parquet    one table per MGS x sample x MAG
OUT/mag-mag/<MGS>/<S>/<MAG>.parquet    the same on the MAG's own contigs
OUT/<mode>/stats/<S>.tsv               every MGS (mag-rep, mag-mag: every MAG) of the sample: status, abundance, counts, fingerprint
OUT/<mode>/stats/<S>.run.tsv           the sample's last run: raw BAM, CRAM, chemistry, ipdSummary arguments and the
                                       kinetics model it used, read accounting
```

Run logs (every step, tool versions, timings, ipdSummary's model line and warnings) are kept only with `-l DIR`:
`DIR/<mode>/<sample>.log` for a one-sample run (appended), else `DIR/<mode>/<date>_<job>.log`.

Every file carries its own `MGS`, `MAG` and `sample` columns (constant columns cost < 0.1 bit per row in parquet), so
files can be read singly or all together:

```sql
-- duckdb
SELECT MGS, sample, count(*) FILTER (WHERE call = 'm6A') AS m6A_sites
FROM read_parquet('OUT/mgs-rep/*/*.parquet') GROUP BY ALL;
SELECT MGS, sample, MAG, count(*) FROM read_parquet('OUT/mag-rep/*/*/*.parquet') GROUP BY ALL;
SELECT * FROM parquet_kv_metadata('OUT/mgs-rep/MGS.0/S1.parquet');   -- the column conventions (key "mmt")
```

```python
import pandas as pd
d = pd.read_parquet("OUT/mgs-rep/MGS.12")            # every sample of one MGS
```

### Columns of a table

ipdSummary's per-site statistics (all of them; it reports no others for native DNA), then MMT's sequence agreement:

| Column | Meaning |
|---|---|
| `MGS`, `MAG`, `sample` | the MGS, the MAG(s) whose reads made the table (mgs-rep: every MAG of the MGS in the sample, comma-separated) and the sample |
| `contig` | representative contig (MATAFILER4 name) |
| `pos` | **1-based** position on the contig (ipdSummary `tpl`) |
| `strand` | `+` / `-`: the strand of the base whose kinetics are reported (`-` = the reverse complement) |
| `base` | that base, on that strand (on `-` the complement of the contig base): A or C, or G/T with a kinetic call |
| `coverage` | kinetic reads (molecules, one IPD each) ipdSummary used on this strand here |
| `score` | -10 log10 p of the observed vs expected IPD (significance: grows with coverage) |
| `tMean`, `tErr` | mean IPD of the reads at this site (outliers capped) and its standard error |
| `modelPrediction` | IPD expected for an unmodified base in this sequence context (the kinetics model) |
| `ipdRatio` | tMean / modelPrediction (effect size: 1 = as expected for an unmodified base) |
| `frac`, `fracLow`, `fracUp` | estimated methylated fraction and its interval (ipdSummary). Only where ipdSummary identified m6A / m4C with > 10 reads; NULL elsewhere, which does NOT mean 0. A partly methylated site is often `modified_base` or `none` with no `frac` - use `ipdRatio` there |
| `call` | `m6A`, `m4C`, `modified_base` (significant kinetic signal, type not identified) or `none` |
| `identificationQv` | confidence that the signal is the identified type (m6A vs m4C vs other) |
| `pile_depth` | reads at this position on this strand (same reads; deletions included) - >= `coverage` |
| `n_ref` | of those, reads with the representative's base |
| `top_alt`, `n_alt` | most common other observation (A/C/G/T/N on this strand, `del`, `ins`) and its read count; NULL when all agree |
| `ctx_alt_pct`, `ctx_alt_n` | within +- `MGSCALL_CTX` bp (10) on this strand: highest % (and count) of reads differing from the representative (substitution, indel, or a soft clip >= `MGSCALL_CLIP_MIN` bp); 0 = clean context |
| `clip_ctx_n` | reads with a clip point within +- `MGSCALL_CTX` bp |

Why the agreement columns: ipdSummary uses every aligned base, mismatches included, and predicts the expected kinetics
from the representative's sequence (kineticsTools' model uses 10 bases on one side of the site and 4 on the other). A sample whose strain differs at
or near a site gives kinetics that are not comparable there; `n_ref / pile_depth` (the site) and `ctx_alt_pct` (the
context, +- 10 bp covers the model's window either way) show where. Values are never changed: downstream decides.

Notes: `modified_base` rises with coverage (significance, not effect size: small systematic deviations become
significant with many reads) - compare samples on `ipdRatio` / `frac` with their coverage, not on call counts. 5mC
gives no reliable kinetic signal on native DNA and is not called; its sites keep their row (C) with whatever signal
there is, so nothing is lost if a 5mC-aware method is applied later.

### `stats/<S>.tsv`

`MGS  MAG  status  fingerprint  abundance  molecules  pseudo_reads  aligned  rows  rep_MAG  updated  mmt`

One row per MGS (mag-rep: per MAG). `abundance` = GeneCat's `MGS.matL0.txt` value for the sample (0 = absent there).

| status | meaning |
|---|---|
| `done` | table written |
| `linked` | the other mode's identical table, hard-linked (one MAG of the MGS in the sample) |
| `no_reads` | no read of this sample on the MAG(s) (after the 0x200 / 0x400 filters) |
| `no_kinetics` | reads found, but no kinetics left after bystrandify |
| `no_aligned` | no read passed the alignment filters on the representative |
| `no_rows` | aligned, but no site with coverage >= `MGSCALL_MIN_ROW_COV` |

`molecules` = HiFi reads selected, `pseudo_reads` = per-strand reads, `aligned` = pseudo-reads kept after alignment.

## Re-running into the same OUT

Every result carries a fingerprint of what it depends on: MMT method version, mode, the `MGSCALL_*` settings and
ipdSummary arguments, the sample's CRAM and raw BAM, the representative's FASTA, and the contigs of the sample's MAGs
of the MGS. A run computes only what is missing or whose fingerprint changed; everything else is kept. `--select` /
`--samples` runs add to the same folders. A result is dropped when its MGS disappears from GeneCat or no longer has a
MAG in the sample. `--force` redoes the selection (computed, never linked). A failed sample does not stop the others; rerun the same command
to redo it (`grep FAIL OUT/mmt_*.log` says why it failed). `mmt status -o OUT` lists, per MGS, the samples done /
linked / not yet and the ones present but not binned.

## HPC notes

- The unit of work is a sample: its raw kinetics BAM is read once per run, whatever the number of MGS.
- `--tmp`: use a shared scratch folder. It holds the index and the reference indexes (caches, shared by all jobs)
  and the work folders (`<tmp>/mmt/<OUT name>_<key>/<mode>/<S>`, removed when the sample ends - done, failed, or
  stopped by scancel / the time limit; `--debug` keeps them, and a rerun resumes from the last finished step).
- Threads and memory come from the SLURM allocation (`--cpus-per-task`, `--mem`) unless `-t` / `--mem` are given.
- Jobs may run at the same time on the same OUT: a sample is locked while it is processed (`stats/<S>.lock`); a lock
  of a job that no longer exists in `squeue` is taken over.
- Terminal / `sbatch` output: a start line, the plan and an end line (after a failure without `-l`, the run log's last
  lines). A pipeline can set `MMT_LOG_DATE` once so every job of one run writes the same `mmt_<date>.log`.
