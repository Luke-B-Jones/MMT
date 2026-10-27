# Changelog

Versions follow `MAJOR.MINOR.PATCH` (file `VERSION`). Results record the version that made them (`_state/*.tsv`,
`mmt_runs.tsv`). A change that alters MGScall results also bumps `MGSCALL_METHOD` in `modules/MGScall.sh`, so the next
run redoes them.

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
