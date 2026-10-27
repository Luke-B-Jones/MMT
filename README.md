# MMT - MF4 Methylation Toolkit

Methylation analyses (PacBio HiFi kinetics: m6A, m4C) that run on top of a finished
[MATAFILER4](https://github.com/hildebra/MATAFILER) run and its GeneCat MGS stage. MMT reads MATAFILER4's own outputs
(assemblies, bins, read alignments) and GeneCat's species bins (MGS) and representative MAGs; it never changes them.

Version: see `VERSION` and `CHANGELOG.md`. `mmt version` also prints the git commit.

## Commands

| Command | What |
|---|---|
| `mmt MGScall --mgs-rep` | per MGS x sample: reads of every MAG of the MGS, aligned to the MGS's representative MAG -> m6A / m4C per site |
| `mmt MGScall --mag-rep` | the same with the representative MAG's own reads only |
| `mmt check` | settings in effect, tool and input checks; with `-o OUT` also builds the index and reports usable samples / MGS |
| `mmt help [command]` | help |

Method, outputs and every column: [docs/MGScall.md](docs/MGScall.md) (also copied into each results folder).

## Install (HPC)

```bash
cd /hpc-home/$USER/tools
git clone https://github.com/Luke-B-Jones/MMT.git
chmod +x MMT/mmt                                   # only if git did not keep the executable bit
mkdir -p ~/bin && ln -sf "$PWD/MMT/mmt" ~/bin/mmt  # ~/bin on PATH
cp MMT/config/site.example.conf MMT/config/site.conf   # then edit the paths (site.conf is not tracked by git)
mmt check -o /path/to/results                      # tools, inputs, index report
```

Update: `git -C /hpc-home/$USER/tools/MMT pull`. Results already in an output folder stay valid unless a change
says otherwise in `CHANGELOG.md` (then the next run redoes only what changed).

### Requirements

| Tool | From |
|---|---|
| python 3 + pysam | micromamba env (`MMT_ENV_PREFIX`, e.g. the micromamba root = base env) |
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
# one MGS, all samples that have it, one job
mmt MGScall --mgs-rep -o /path/results --select MGS.12

# see what would run
mmt MGScall --mgs-rep -o /path/results --dry-run
```

SLURM array, one sample per task (the efficient unit: each raw BAM is read once per run):

```bash
OUT=/path/results; SEL=MGS.12,MGS.40
mmt MGScall --mgs-rep -o $OUT --select $SEL --list-samples > todo.txt
sbatch --array=1-$(wc -l < todo.txt)%20 -c 16 --mem=64G -t 1-00:00 --wrap \
  'S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" todo.txt); mmt MGScall --mgs-rep -o '$OUT' --select '$SEL' --samples "$S" --tmp "$TMPDIR"'
```

The output folder is meant to be reused: every run adds what is missing, keeps what is current and redoes what
changed (inputs, settings, MMT method). Failed samples do not stop the others; rerun the same command to resume.

## Output folder

```
OUT/
  README.md               what each folder holds (written by MMT)
  mmt_runs.tsv            one line per run: time, version, job, command, samples done / failed, log
  index/                  MATAFILER4 / GeneCat layout as MMT reads it (+ reads/sample=<S>/ : read -> MAG)
  MGScall/<mode>/         MGS=<MGS>/sample=<S>/data.parquet, _state/
  logs/                   one log per run and command
```

## Repository

```
mmt                 wrapper (commands, help, version)
VERSION, CHANGELOG.md
config/             mmt.conf (defaults, documented), site.example.conf (NBI HPC paths)
lib/common.sh       logging, settings, tools, locks, fingerprints, index
lib/mf4.py          MATAFILER4 / GeneCat layout (index), read assignment, read splitting
lib/kinetics.sh     ipdSummary + pileup -> one kinetics table
modules/            one file per command
docs/               method and output description per command
test/               synthetic MATAFILER4 run with known methylation + end-to-end test
```

Adding a command: `modules/<Name>.sh` with `module_help` and `module_main`, one line in `module_file` in `mmt`,
`docs/<Name>.md`, a `CHANGELOG.md` entry. Results go to `OUT/<Name>/`, logs to `OUT/logs/<Name>/`.

### Test

`test/run_test.sh WORKDIR` builds a small synthetic MATAFILER4 run (two species, strains, a co-assembly, QC-fail and
duplicate flags, known m6A fractions) and runs every command against it, checking results against the truth, reruns,
selection, settings changes, failure + resume and locks. It needs samtools, minimap2, a python with pysam, DuckDB and a
kineticsTools `ipdSummary`; SMRT Link's other tools are replaced by small stand-ins (`test/standins.sh`).
