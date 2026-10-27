# MGScall - methylation of each MGS on its representative MAG

`mmt MGScall (--mgs-rep | --mag-rep) -o OUT [--select MGS,...] [--samples S,...]`

For every sample and every selected MGS (GeneCat species bin), MGScall takes the sample's PacBio HiFi reads that
MATAFILER4 placed on that MGS, aligns them to the MGS's representative MAG, and calls m6A and m4C with ipdSummary
from the read kinetics. Result: one table per MGS x sample, every A and C of the representative on both strands.

## Modes (one per run; each has its own results folder)

| Mode | Reads used for an MGS in a sample |
|---|---|
| `--mgs-rep` | reads on ANY MAG of the MGS in that sample's assembly (the species as present in the sample) |
| `--mag-rep` | reads on the REPRESENTATIVE MAG only - so only the sample(s) whose assembly produced it (co-assembly: every sample of the group) |

## Method

1. **Layout** (`OUT/index/`, rebuilt when the map or `MAGvsGC.txt.gz` changes). From the MATAFILER4 map, run folder
   and GeneCat `Bin_<binner>/MAGvsGC.txt.gz`, with GeneCat's own rules: representative = starred MAG of the MGS,
   else the first non-canopy MAG; MAG `<sample or alias>__<bin>` -> sample (co-assembly alias `<SmplID>M<n>`) ->
   assembly (`assemblies/metag/assembly.txt`) -> contigs (`Binning/<binner>/<SmplID>`). The representative's
   sequence is GeneCat's own `Genomes/MGS_ctg/<MGS>.ctgs.<MAG>[.fna].gz`; its contigs must equal the bin-file contigs,
   otherwise MMT stops (its naming rules would not match this GeneCat version).
2. **Read -> MAG** (`OUT/index/reads/sample=<S>/`, shared by both modes). Each primary alignment in MATAFILER4's
   `<S>/mapping/<S>-smd.cram` (minimap2 map-hifi to the sample's own assembly) whose contig belongs to an MGS MAG.
   Alignments MATAFILER4 flagged as failed (0x200: > 5 % edits, < 50 % of the read aligned or MAPQ < 30) are never
   used; duplicates (0x400) only with `MGSCALL_KEEP_DUP=1`. No new competitive mapping: MATAFILER4 already decided.
3. **Kinetics reads.** The selected reads are taken from the sample's raw kinetics BAM (map `#DirPath` +
   `SmplPrefix`.bam; one pass), split into one pseudo-read per strand (`ccs-kinetics-bystrandify`) and grouped by MGS.
4. **Alignment** of each MGS's reads to its own representative only (`pbmm2 --preset CCS -y MGSCALL_MIN_IDENTITY`).
   Kept: primary, MAPQ >= `MGSCALL_MIN_MAPQ` (drops repeats within the representative), no hard clip, aligned part
   >= `MGSCALL_MIN_PCT_ALIGNED` % of the read and >= `MGSCALL_MIN_ALN_BP` bp. Reads of a strain that differs from
   the representative are kept as long as they align at that identity; the sequence-agreement columns show where
   they differ.
5. **Kinetics** - one ipdSummary run per sample on all its MGS (`--identify m6A,m4C --methylFraction
   --mapQvThreshold MGSCALL_MIN_MAPQ --minCoverage 1 --maxCoverage 100000`) and, from the same aligned reads,
   `samtools mpileup` for the sequence-agreement columns.
6. **Tables.** Rows for A and C (on the row's strand) with kinetic coverage >= `MGSCALL_MIN_ROW_COV`, written per MGS.

## Output

```
OUT/MGScall/<mode>/MGS=<MGS>/sample=<S>/data.parquet   one table per MGS x sample
OUT/MGScall/<mode>/_state/<S>.tsv                      every MGS of the sample: status, read counts, fingerprint
OUT/MGScall/<mode>/_state/<S>.stats.tsv                read accounting of the sample's last run
```

Read a whole mode (the folder names become the columns `MGS` and `sample`):

```sql
-- duckdb
SELECT MGS, sample, count(*) FILTER (WHERE call = 'm6A') AS m6A_sites
FROM read_parquet('OUT/MGScall/mgs-rep/*/*/*.parquet', hive_partitioning = true)
GROUP BY ALL;
```

```python
import pandas as pd
d = pd.read_parquet("OUT/MGScall/mgs-rep", filters=[("MGS", "=", "MGS.12")])   # pyarrow reads the partitions
```

### Columns of `data.parquet`

| Column | Meaning |
|---|---|
| `contig` | representative contig (original MATAFILER4 name) |
| `pos` | 1-based position |
| `strand` | `+` / `-` (strand of the base whose kinetics are reported) |
| `base` | A or C on that strand |
| `coverage` | reads whose kinetics ipdSummary used here |
| `score` | ipdSummary modification score (-10 log10 p) |
| `ipdRatio`, `tMean`, `tErr`, `modelPrediction` | observed / expected inter-pulse duration (ipdSummary) |
| `frac`, `fracLow`, `fracUp` | estimated methylated fraction and 95 % interval (ipdSummary). Given only where ipdSummary identified m6A / m4C with > 10 reads; NULL elsewhere, which does NOT mean 0. A partly methylated site is often called `modified_base` or `none` with no `frac` - use `ipdRatio` there |
| `call` | `m6A`, `m4C`, `modified_base` (signal, type not identified) or `none` |
| `identificationQv` | confidence of the m6A / m4C identification |
| `pile_depth` | reads at this position on this strand in the pileup (deletions included) |
| `n_ref` | of those, reads with the representative's base |
| `top_alt`, `n_alt` | most common other observation (A/C/G/T/N on this strand, `del`, `ins`) and its read count; NULL when all agree |
| `ctx_alt_pct`, `ctx_alt_n` | within +- `MGSCALL_CTX` bp on this strand: highest % (and count) of reads differing from the representative (substitution, indel, or a soft clip >= `MGSCALL_CLIP_MIN` bp); 0 = clean context |
| `clip_ctx_n` | reads with a clip point within +- `MGSCALL_CTX` bp |

Why the agreement columns: ipdSummary uses every aligned base, mismatches included, and predicts the expected kinetics
from the representative's sequence. A sample whose strain differs at or near a site gives kinetics that cannot be
compared with other samples there; `n_ref / pile_depth` and `ctx_alt_pct` show where.

### `_state/<S>.tsv`

`MGS  status  fingerprint  molecules  pseudo_reads  aligned  rows  updated  mmt`

| status | meaning |
|---|---|
| `done` | table written |
| `no_reads` | no read of this sample lies on the MGS's MAG(s) (after the 0x200 / 0x400 filters) |
| `no_kinetics` | reads found, but no kinetics left after bystrandify |
| `no_aligned` | no read passed the alignment filters on the representative |
| `no_rows` | aligned, but no A/C with coverage >= `MGSCALL_MIN_ROW_COV` |

`molecules` = HiFi reads selected, `pseudo_reads` = per-strand reads, `aligned` = pseudo-reads kept after alignment.

## Re-running into the same OUT

Every result carries a fingerprint of what it depends on: MMT method version, mode, the `MGSCALL_*` settings, the
index (map + MAGvsGC), the sample's CRAM and raw BAM, and the representative's FASTA. A run computes only what is
missing or whose fingerprint changed; everything else is kept. `--select` / `--samples` runs add to the same folders.
A result is dropped when its MGS disappears from GeneCat or no longer has a MAG in the sample. `--force` redoes the
selection. A failed sample does not stop the others; rerun the same command to resume.

## HPC notes

- The unit of work is a sample: its raw kinetics BAM is read once per run, whatever the number of MGS. Selecting
  more MGS per run costs little extra; one MGS per run over all samples costs a full pass over every raw BAM.
- `--tmp "$TMPDIR"` (node-local) is fastest. A shared scratch folder instead keeps the reference indexes
  (`<tmp>/mmt/refs/`) between jobs and lets a killed job resume its sample from the last finished step.
- Threads and memory come from the SLURM allocation (`--cpus-per-task`, `--mem`) unless `-t` / `--mem` are given.
- Jobs may run at the same time on the same OUT: a sample is locked while it is processed (`_state/<S>.lock`); a lock
  of a job that no longer exists in `squeue` is taken over.
- Terminal / `sbatch` output: a start line, the plan and an end line. Everything else is in
  `OUT/logs/MGScall/<mode>/<date>_<job>.log` (+ `ipdSummary/<S>.log`).
