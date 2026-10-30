# Changelog

Versions follow `MAJOR.MINOR.PATCH` (file `VERSION`). Results record the version that made them (`stats/*.tsv`,
the mmt log). A change that alters results also bumps the command's method version (`MGSCALL_METHOD` in
`modules/MGScall.sh`, `DEREP_METHOD` in `modules/derep-ctg.sh`, `CTGMETH_METHOD` / `LIFT_METHOD` in
`modules/ctg-meth.sh`, `ASSIGN_METHOD` in `lib/reads.sh`), so the next run redoes them.

## 0.6.0 - 2026-10-29

- **`MGScall --mag-mag`** - each MAG's reads (as in mag-rep) called on that MAG's own contigs instead of the
  representative: `mag-mag/<MGS>/<S>/<MAG>.parquet`, `contig` = the MAG's contig (MATAFILER4 name). The MAG references
  are cut from the sample's assembly (`scaffolds.fasta.filt[.gz]`) and cached in the work space (shared with ctg-meth).
  Same samples, filters, columns, stats and layers as mag-rep; nothing is linked (other reference). No change to the
  results of the other modes (no re-run).

## 0.5.0 - 2026-10-29

Output folder: only results. `-o OUT` holds one folder per analysis and the dated log - nothing else.
- **Layout** - `OUT/mgs-rep/`, `OUT/mag-rep/` (no `MGScall/` level), `OUT/derep-ctg/`, `OUT/ctg-meth/`; `_state/` is
  now `stats/` (`<S>.tsv` per MGS, `<S>.run.tsv` = the sample's last run, was `<S>.stats.tsv`). No README files in OUT
  (documentation lives in `docs/`). No `OUT/index/`: the MATAFILER4 / GeneCat index and the per-sample read -> MAG
  tables are a cache in the work space (`<tmp>/mmt/<OUT name>_<key>/index/`, rebuilt when missing; use shared scratch
  for `--tmp` on a cluster). The 1 kb read-depth tables are no longer written (nothing used them).
- **Logs** - the folder's log is `OUT/mmt_<date>.log` (today's date, or `MMT_LOG_DATE` - a pipeline sets it once so a
  whole run writes one file). Run logs only with **`-l/--log DIR`** (all commands): `DIR/<analysis>/<sample>.log` for a
  one-sample run (appended), else `DIR/<analysis>/<date>_<job>.log`; ipdSummary's log is no longer kept apart - the
  kinetics model it loaded, its warnings and, after a failure, its last lines go into the run log. Without `-l` no run
  log is kept; after a failure its last 25 lines are printed to the terminal (sbatch .out).
- **mag-rep = every MAG apart** - for every sample with a MAG of the MGS (was: the representative's own sample only),
  each MAG's reads on the representative: `mag-rep/<MGS>/<S>/<MAG>.parquet` (strains apart, cross-sample). Where a
  sample has one MAG of the MGS the table equals mgs-rep's: the mode run second hard-links it (status `linked`) instead
  of computing it (not with `--force`: that always computes). Several MAGs of one MGS are called in separate
  ipdSummary runs ("layers").
- **`MAG` column** in every MGScall table (mgs-rep: the MAGs pooled, comma-separated). Column conventions (1-based
  positions, strand / base, rows kept) are stored in each file's parquet metadata (key `mmt`).
- **G/T rows with a kinetic signal** - besides every A and C, G/T rows that ipdSummary flags (`modified_base`) are kept
  (other / unknown modifications); before, they were dropped. (ctg-meth: A/C only, as before.)
- **Presence** - the 10-molecule rule is gone (MGScall: any reads; ctg-meth `CTGMETH_PRESENT_MOL` default 1).
  GeneCat's MGS abundance table (`Bin_<binner>/Annotation/Abundance/MGS.matL0.txt`) is read into the index: its value
  is in every stats row (`abundance`); the PLAN lines and `mmt status` list samples where the MGS is present but was
  not binned. MGScall status `absent` is now `no_reads`.
- **`--ipd-args "ARGS"` / `MMT_IPD_ARGS`** - extra ipdSummary arguments, as they are (e.g. `--ipdModel FILE`); part of
  every result's fingerprint. The kinetics model ipdSummary used and the raw BAM's chemistry are recorded per sample
  (`stats/<S>.run.tsv`); `mmt status` lists the models.
- `mmt status` reads the index (rebuilds it if missing) and reports `linked`, MAG tables, not-binned samples, models.
- Removed: migration of MMT <= 0.2 layouts, `modules/dereplicate-contigs.sh` (old copy; `dereplicate-contigs` stays an
  alias of `derep-ctg`).
- Upgrading: `MGSCALL_METHOD` 2 - every MGScall result is redone (new columns and rows); the results of 0.4 in
  `OUT/MGScall/`, `OUT/index/`, `OUT/logs/`, `OUT/README.md` and `OUT/mmt.log` are not used any more: delete them.

## 0.4.0 - 2026-10-28

- **`OUT/mmt.log`** - one human-readable log per output folder, appended by every run of every command (one line per
  event, each line stands alone: `<date> <time>  <command>  <EVENT>  <text>`): START (job, host, version, selection,
  run log), PLAN per MGS (in how many samples it has a MAG / the representative, to do, up to date, and where it is not),
  one SAMPLE line per sample with every MGS's outcome (`done` with rows and molecules, `absent`, `no_kinetics`, ...,
  dropped results), FAIL with the reason and the run log, END (processed / failed / up to date, minutes), STOP for runs
  ended early (error, scancel, time limit), INDEX when the index is (re)built, WARN. Written under a lock: array tasks
  on many nodes append safely. Replaces `mmt_runs.tsv` (no longer written). Dry runs are not logged.
- **`mmt status -o OUT [--select IDS]`** (new) - per MGS and command (MGScall per mode, ctg-meth): the samples it is
  planned for and the state of each (done / absent / ... / not yet) with the sample names; derep-ctg: which MGS are
  done. Printed and appended to `OUT/mmt.log`. Reads only OUT (no tools, no inputs).
- **Presence check** (MGScall, ctg-meth) - before a sample's raw BAM is read, the molecules MATAFILER4 placed on each
  MGS's MAG(s) are counted from the read index; fewer than `MGSCALL_PRESENT_MOL` / `CTGMETH_PRESENT_MOL` (10) = state
  `absent` (replaces `no_reads`): recorded, not processed, no table; when every MGS of a sample is absent its raw BAM
  is not read. (A sample is planned for an MGS only when its assembly has a MAG of it, as before.)
- **Work folders removed whatever happens** - a sample's work folder is removed when the sample ends: done, failed,
  or stopped by scancel / the time limit (SIGTERM). `--debug` (all commands; `--keep-tmp` is the same) keeps them and
  a rerun resumes from the last finished step; a work folder left by a job killed outright (SIGKILL) is resumed too.
  Work folders are named `<tmp>/mmt/<OUT name>_<key>/` (was `out_<key>`). derep-ctg: `--debug` keeps the instance
  alignments.
- Upgrading an output folder: the presence setting is part of every MGScall / ctg-meth result's fingerprint, so the
  first 0.4 run redoes the results of 0.3 (MGScall states `no_reads` become `absent`). Work folders of 0.3
  (`<tmp>/mmt/out_<key>/`) can be deleted.
- Tests: MGScall, ctg-meth and derep-ctg suites check mmt.log (START/PLAN/SAMPLE/FAIL/END, START = END + STOP),
  work folders removed after failures without `--debug` and kept with it, resume from a kept folder, absent MGS (state,
  no file, raw BAM not read, rerun current), `mmt status`.

## 0.3.0 - 2026-10-28

- `mmt ctg-meth` (new): per MGS x sample, m6A / m4C of every MAG's own contigs from the reads MATAFILER4 placed on the
  MAG, each row placed on the derep-ctg blocks (block, block_pos, block_strand, instance, ctx_dist) and the MAG's genes
  (gene, up_gene) ; plus read sets per breakpoint: large rearrangements (split alignments, CIGAR indels >=
  `CTGMETH_MIN_SV`, short inversions the aligner passed through, contig ends joined elsewhere) grouped by breakpoint and
  partner - each conformation with enough molecules (ALT) and the as-assembled one at the same breakpoint (REF) called
  from its own reads on the flanks all of them cover. Per MGS: `<S>.parquet` (sites), `junctions/`, `edges/`, `reads/`
  (`<S>.parquet` each) ; `_state/`. Settings `CTGMETH_*`; `docs/ctg-meth.md` (method, every column, comparison rules,
  DuckDB recipes for block patterns, conformations, a chromosome frame, GC skew + depth, sequence context, genes).
  Needs python >= 3.10 with pysam, pyarrow, numpy, edlib.
- `mmt derep-ctg` = `dereplicate-contigs` renamed (the old name still works); results in `OUT/derep-ctg/<MGS>/`
  (`OUT/dereplicate-contigs/MGS=<MGS>/` of MMT <= 0.2 is moved by the first run). New tables: `lift.parquet` (every base
  of every instance on its block sequence), `composition.parquet` (base counts per 1 kb window), `contig_seqs.parquet`
  (MAG contig sequences). Block sequences are polished to the majority base of their instances. Gene-order indices are
  named `slot` (`block_genes`, `genes`, `families`: `pos` -> `slot` ; `instances`: `first_pos` / `last_pos` / `host_pos`
  -> `first_slot` / `last_slot` / `host_slot` ; `links`: `pos_a` / `pos_b` -> `slot_a` / `slot_b`) ; `catalogue` in
  `block_genes` / `families` is a list of catalogue gene IDs. `DEREP_METHOD` 3: every MGS is redone.
- MGScall: `OUT/MGScall/<mode>/<MGS>/<S>.parquet` with columns `MGS` and `sample` (plain folder names, read with
  `read_parquet('OUT/MGScall/mgs-rep/*/*.parquet')`); results of MMT <= 0.2 are moved into this layout, not recomputed.
- Read index shared by MGScall and ctg-meth: `OUT/index/reads/<S>.parquet` (+ `sample`; duplicates kept, flagged
  `dup`, each command filters them by its setting) and new `OUT/index/depth/<S>.parquet` (read depth of every MGS-MAG
  contig in 1 kb bins, supplementary alignments included) ; rebuilt once. `OUT/index/movies.tsv`: read groups and
  chemistry of every raw BAM (batch key) ; the index is rebuilt once (`INDEX_METHOD` 3).
- Positions are signed 32-bit integers in every table (differences never overflow).
- Work folders are kept per output folder (`<tmp>/mmt/out_<key>/`): output folders can share one `--tmp`; reference
  caches are keyed by their source (MGS + representative FASTA, MAG + assembly).
- Upgrading an output folder: run one MGScall command on it before an array (e.g. the `--list-samples` of the array
  recipe) - it moves the 0.2 results in one go (a damaged old file is renamed `*.unreadable` and recomputed). Work
  folders of MMT <= 0.2 (`<tmp>/mmt/MGScall_<mode>/`) are no longer used and can be deleted.
- Tests: `test/run_ctgmeth_test.sh` (+ `ctgmeth_mock.py`, `check_ctgmeth.py`: planted inversions, deletion,
  insertion, misassembly, contig break, circular plasmid, strain differences, conformation-specific methylation);
  MGScall and derep-ctg tests cover the new layouts and their migration.

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
