# Changelog

Versions follow `MAJOR.MINOR.PATCH` (file `VERSION`). Results record the version that made them (`_state/*.tsv`,
`mmt_runs.tsv`). A change that alters MGScall results also bumps `MGSCALL_METHOD` in `modules/MGScall.sh` (for
dereplicate-contigs: `DEREP_METHOD` in `modules/dereplicate-contigs.sh`), so the next run redoes them.

## 0.2.0 - 2026-10-26

- `mmt dereplicate-contigs`: each MGS as a short list of collinear gene blocks ("dereplicated contigs") built from
  all its MAGs - gene families (GeneCat catalogue + context-constrained merging of divergent alleles), core and
  accessory blocks, every instance in every MAG (incomplete-aware ends, nesting), one representative sequence per block
  (`blocks.fna.gz`), instance identity and alignments (minimap2), block links along contigs, and parquet tables that
  trace every gene, contig and MAG. Shared preparation once per index (`--prep-only`); `--mgs`, `--part I/N` (SLURM
  arrays, balanced by MAG count), `--list-mgs`, `--dry-run`, `--force`. Settings `DEREP_*` in `config/mmt.conf`;
  method and every column in `docs/dereplicate-contigs.md`.
- Requirements for it: minimap2 (`MMT_MINIMAP2`), python >= 3.10 with pyarrow and edlib. `mmt check` reports them.
- Test: `test/run_derep_test.sh` (synthetic run with known blocks, rearrangements, islands, alleles, a MAG of co-binned
  strains; reruns, parts, settings and input changes, locks, failure + resume); `test/derep_scale.py` (timing).
- Index (all commands; `INDEX_METHOD` 2, so the next run of any command rebuilds it once - results are redone only
  where their inputs changed): also rebuilt when a bin file or `assembly.txt` changes (`index/inputs.tsv`); an MGS
  whose rows are not grouped in MAGvsGC gets GeneCat's starred MAG as representative (was: the first MAG of its last
  group of rows); maps with Windows line ends; a warning on the terminal when member MAGs cannot be used.
- Locks: a stale lock is taken over by one job only (several jobs starting after a killed one could all take it); a
  busy SLURM controller no longer makes a live job's lock look stale.
- Usage errors reach the terminal also after the run log has started; ID lists (`--mgs`, `--select`, `--samples`)
  may also be separated by spaces.

## 0.1.0 - 2026-10-26

First version.

- `mmt` wrapper: commands, help, version (with git commit), settings from `config/mmt.conf` -> `config/site.conf` ->
  `$MMT_CONFIG` -> `--config` -> options; one log per run, a few lines on the terminal.
- `mmt check`: settings in effect, tool checks, input checks, index report.
- `mmt MGScall --mgs-rep | --mag-rep`: per-MGS x sample methylation (m6A, m4C) on GeneCat's representative MAG,
  reads as placed by MATAFILER4 (CRAM), kinetics from the raw BAM; `--select`, `--samples`, `--list-samples`,
  `--dry-run`, `--force`, `--keep-tmp`. Results updated in place; fingerprints decide what is redone.
- Index of the MATAFILER4 / GeneCat layout (`OUT/index/`), built with GeneCat's own naming rules and checked against
  GeneCat's representative FASTAs.
- Test: synthetic MATAFILER4 run with known methylation (`test/`).
