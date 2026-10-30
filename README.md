# MMT - MF4 Methylation Toolkit

Methylation analyses (PacBio HiFi kinetics: m6A, m4C) that run on top of a finished
[MATAFILER4](https://github.com/hildebra/MATAFILER) run and its GeneCat MGS stage. MMT reads MATAFILER4's own outputs
(assemblies, bins, read alignments) and GeneCat's species bins (MGS) and representative MAGs; it never changes them.

Version: see `VERSION` and `CHANGELOG.md`. `mmt version` also prints the git commit.

## Commands

| Command | What |
|---|---|
| `mmt MGScall --mgs-rep` | per MGS x sample: reads of every MAG of the MGS in the sample, pooled, aligned to the MGS's representative MAG -> m6A / m4C per site |
| `mmt MGScall --mag-rep` | the same per MAG: each of the sample's MAGs of the MGS apart (strains apart) ; one-MAG samples link the mgs-rep table |
| `mmt MGScall --mag-mag` | each MAG's reads on that MAG's own contigs (not the representative) |
| `mmt derep-ctg` | per MGS: a short list of collinear blocks ("dereplicated contigs") from all its MAGs, with every occurrence in every MAG down to the base (lift table), block sequences, block order along the contigs, genes, gene families and the MAG contig sequences (old name `dereplicate-contigs` still works) |
| `mmt ctg-meth` | per MGS x sample: every MAG's own contigs -> m6A / m4C per site, placed on the derep-ctg blocks; plus read sets per junction (other conformations of a contig, and the as-assembled one at the same breakpoint) with their own calls on the junction flanks |
| `mmt check` | settings in effect, tool and input checks; with `-o OUT` also builds the index and reports usable samples / MGS |
| `mmt status -o OUT` | per MGS and command: the samples it is planned for and their state (done / linked / not yet ...), with names, and where GeneCat has it without a MAG; also written to the mmt log |
| `mmt help [command]` | help |

Method, outputs and every column: [docs/MGScall.md](docs/MGScall.md), [docs/derep-ctg.md](docs/derep-ctg.md),
[docs/ctg-meth.md](docs/ctg-meth.md).

## Install (HPC)

```bash
cd /hpc-home/$USER/tools
git clone https://github.com/Luke-B-Jones/MMT.git
chmod +x MMT/mmt                                   # only if git did not keep the executable bit
cp MMT/config/site.example.conf MMT/config/site.conf   # then edit the paths (site.conf is not tracked by git)

# put `mmt` on your PATH (once): a link in ~/bin, and ~/bin on PATH in ~/.bashrc unless it is there already
mkdir -p ~/bin && ln -sf "$PWD/MMT/mmt" ~/bin/mmt
[[ :$PATH: == *":$HOME/bin:"* ]] || echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
mmt version                                        # works from any folder now

mmt check -o /path/to/results                      # tools, inputs, index report
```

Use PATH, not an alias: SLURM jobs (`sbatch --wrap 'mmt ...'`) inherit PATH from the shell they are submitted from,
but aliases exist only in interactive shells. Without the PATH step, call MMT by its full path
(`/hpc-home/$USER/tools/MMT/mmt`).

Update: `git -C /hpc-home/$USER/tools/MMT pull` (or, to match GitHub exactly: `git fetch origin && git reset --hard
origin/main` - `config/site.conf` is not tracked and stays). Results already in an output folder stay valid unless a
change says otherwise in `CHANGELOG.md` (then the next run redoes only what changed).

### Requirements

| Tool | From |
|---|---|
| python 3 + pysam, pyarrow, numpy | micromamba env (`MMT_ENV_PREFIX`, e.g. the micromamba root = base env) |
| derep-ctg, ctg-meth: python >= 3.10 with edlib (+ pandas for the tests); derep-ctg: minimap2 | same env; `MMT_MINIMAP2` (e.g. the MF4 env's) |
| samtools >= 1.15 | micromamba env (`MMT_SAMTOOLS`) |
| DuckDB command-line binary (tested 1.1.3) | `MMT_DUCKDB` |
| SMRT Link: `ipdSummary`, `ccs-kinetics-bystrandify`, `pbmm2`, `pbindex` | installed separately, `MMT_SMRT_BIN` = its `smrtcmds/bin` |
| MATAFILER4 run + GeneCat MGS stage | `MMT_MAP` (the run's map), `MMT_GENECAT`; the run folder comes from the map's `#OutPath`/`#RunID` (or `-i`) |

The raw kinetics BAMs (with `fi/fp/ri/rp`) are found through the map: `#DirPath` + `SmplPrefix` + `.bam`.

## Settings

`config/mmt.conf` documents every setting with its default. Later sources win:
`config/mmt.conf` -> `config/site.conf` -> `$MMT_CONFIG` -> `--config FILE` -> command-line options.
Threads and memory default to the SLURM allocation.

## Running

```bash
# one MGS, all samples that have it, one job ; then each MAG apart
mmt MGScall --mgs-rep -o /path/results --select MGS.12
mmt MGScall --mag-rep -o /path/results --select MGS.12
mmt MGScall --mag-mag -o /path/results --select MGS.12

# see what would run
mmt MGScall --mgs-rep -o /path/results --dry-run

# dereplicated contigs (blocks) of two MGS, then their per-MAG contig methylation (all samples, one job)
mmt derep-ctg -o /path/results --mgs MGS.12,MGS.40
mmt ctg-meth -o /path/results --select MGS.12,MGS.40
```

SLURM array, one sample per task (the efficient unit: each raw BAM is read once per run); `--tmp` on shared scratch
(it holds the index and reference caches the jobs share), run logs kept with `-l`:

```bash
OUT=/path/results; SEL=MGS.12,MGS.40; TMP=/scratch/me/mmt_tmp; LOGS=/scratch/me/mmt_logs
mmt MGScall --mgs-rep -o $OUT --select $SEL --tmp $TMP --list-samples > todo.txt
sbatch --array=1-$(wc -l < todo.txt)%20 -c 16 --mem=64G -t 1-00:00:00 --wrap \
  'S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" todo.txt); mmt MGScall --mgs-rep -o '$OUT' --select '$SEL' --samples "$S" --tmp '$TMP' -l '$LOGS
```

`ctg-meth` runs the same way as an array (`mmt ctg-meth -o $OUT --select $SEL --list-samples`, then one sample per
task); it needs a current `derep-ctg` result for each MGS.

The output folder is meant to be reused: every run adds what is missing, keeps what is current and redoes what
changed (inputs, settings, MMT method). Failed samples do not stop the others; rerun the same command to redo them.

- **One log per output folder and date: `OUT/mmt_<date>.log`.** Every run of every command appends: its start (job,
  version, selection), its plan per MGS (in how many samples the MGS has a MAG, to do / up to date, where GeneCat's
  abundance table has it without a MAG), one line per sample with each MGS's outcome (`done` with rows and molecules,
  `linked`, `no_reads`, ...), failures with their reason, and its end. `grep MGS.12 OUT/mmt_*.log` is the history of
  an MGS; `mmt status -o OUT` adds a summary per MGS. A pipeline sets `MMT_LOG_DATE` once: one file per pipeline run.
- **Run logs only with `-l DIR`** (`DIR/<analysis>/<sample>.log` for a one-sample run, appended; ipdSummary's model
  line and warnings included). Without `-l` no run log is kept; after a failure its last lines go to the terminal.
- **Presence.** A sample is processed for an MGS when its assembly has a MAG of it. GeneCat's abundance table
  (`MGS.matL0.txt`) is recorded per result and lists where the MGS is present but was not binned.
- **Work space (`--tmp`; on a cluster a shared scratch folder).** Each sample's work folder is removed when the sample
  ends - done, failed or stopped. `--debug` keeps it (and a rerun resumes from it). Caches stay for the next jobs: the
  index of the MATAFILER4 / GeneCat layout and the per-sample read -> MAG tables (`<tmp>/mmt/<OUT>_<key>/index`) and
  the reference indexes (`<tmp>/mmt/refs*/`); they are rebuilt when missing - delete `<tmp>/mmt` when done.
- **Kinetics model.** ipdSummary picks the model matching each BAM's chemistry; `--ipd-args "--ipdModel FILE"` (or
  `MMT_IPD_ARGS`) forces one. The model used is recorded per sample (`<analysis>/stats/<S>.run.tsv`).

## Output folder

```
OUT/
  mmt_<date>.log          the folder's log: every run, plan per MGS, each sample x MGS outcome, failures, status
  mgs-rep/                <MGS>/<S>.parquet ; stats/<S>.tsv (per MGS) + stats/<S>.run.tsv (the sample's last run)
  mag-rep/                <MGS>/<S>/<MAG>.parquet ; stats/
  mag-mag/                <MGS>/<S>/<MAG>.parquet ; stats/
  derep-ctg/              <MGS>/ (blocks.fna.gz + parquet tables incl. lift, genes, contig_seqs), summary.tsv, _prep/
  ctg-meth/               <MGS>/<S>.parquet (sites), <MGS>/junctions|edges|reads/<S>.parquet, stats/
```

Nothing else: no copies of MATAFILER4 / GeneCat tables (they live in the work space as a cache), no logs (`-l`),
no documentation files (this repository's `docs/`). Folder names are plain (`MGS.12`); every table carries its own
`MGS` / `MAG` / `sample` columns, so `read_parquet('OUT/mgs-rep/*/*.parquet')` or
`read_parquet('OUT/ctg-meth/MGS.12/*.parquet')` reads a whole set. Positions are 1-based; each table's column
conventions are in its parquet metadata (key `mmt`).

## Repository

```
mmt                 wrapper (commands, help, version)
VERSION, CHANGELOG.md
config/             mmt.conf (defaults, documented), site.example.conf (NBI HPC paths)
lib/common.sh       logging, settings, tools, locks, fingerprints, index
lib/mf4.py          MATAFILER4 / GeneCat layout (index), read assignment, read splitting
lib/reads.sh        the shared read index (read -> MAG per sample, a cache next to the index)
lib/mmtpq.py        parquet schemas and writers of the per-sample tables
lib/derep.py        derep-ctg engine (gene families, blocks, instances, lift, block sequences)
lib/ctgmeth.py      ctg-meth engine (MAG references, junctions, windows, lift onto the blocks)
lib/kinetics.sh     ipdSummary + pileup -> one kinetics table
modules/            one file per command
docs/               method and output description per command
test/               synthetic MATAFILER4 run with known methylation + end-to-end test
```

Adding a command: `modules/<Name>.sh` with `module_help` and `module_main`, one line in `module_file` in `mmt`,
`docs/<Name>.md`, a `CHANGELOG.md` entry. Results go to `OUT/<name>/` (+ `stats/`), run logs to `<-l DIR>/<name>/`.

### Test

`test/run_test.sh WORKDIR` builds a small synthetic MATAFILER4 run (two species, strains, a co-assembly, QC-fail and
duplicate flags, known m6A fractions) and runs every command against it, checking results against the truth, reruns,
selection, settings changes, failure + resume and locks. It needs samtools, minimap2, a python with pysam, DuckDB and a
kineticsTools `ipdSummary`; SMRT Link's other tools are replaced by small stand-ins (`test/standins.sh`).
`test/run_derep_test.sh WORKDIR` does the same for `derep-ctg` (synthetic MGS with known blocks; needs
samtools, minimap2 and a python >= 3.10 with pysam, pyarrow, pandas, edlib; ~1 min).
`test/run_ctgmeth_test.sh WORKDIR` for `ctg-meth` (synthetic MATAFILER4 run with planted rearrangements, strain
differences and conformation-specific methylation; same needs plus a kineticsTools `ipdSummary`).
