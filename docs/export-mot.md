# mmt export-mot

A finished MATAFILER4 run and its GeneCat MGS stage become **one table with exactly one row per final MAG x sample**,
including MGS identity, normalized abundance, coverage evidence and paired raw-BAM/read-selection lists. One FASTA
per final MAG and the molecule selections remain in the accompanying library. Nothing in the MF4 run is changed;
no BAM or CRAM is copied.

**Schema 2 (MMT 0.9.0) requires an updated MOT reader.** Older MOT rejects the new columns and JSON lists. There is
one exported input/abundance table, rather than a separate metadata table.

```bash
mmt export-mot -i /path/METH_PearlAge -gc /path/GC_METH_PearlAge -o /path/MOT/mot_input.tsv --binner SB -t 8
# Feed this table to MOT after its reader supports schema 2 (contract below).
```

Run it as a job (one node, `-t` = CRAMs read at once; memory about 0.5 GB per million CRAM records per CRAM read at
once, e.g. `-t 8 --mem 32G`). A rerun redoes only samples / assemblies whose inputs or settings changed. Exit status:
0 COMPLETE, 3 PARTIAL, 1 FAILED.

## What it reads (all found from `-i` and `-gc`)

| What | Where | Why that file |
|---|---|---|
| maps | `--map`, else the maps GeneCat was run with (`<GC>/LOGandSUB/GCmaps.ori`), else GeneCat's copies (`LOGandSUB/inmap.txt` -> `map.N.txt`) | samples, co-assembly groups (aliases), BAM search rule |
| samples | map rows with a folder `<run>/<SmplID>` | the run's samples (rows of other runs in a combined GeneCat are skipped) |
| assembly | first line of `<run>/<S>/assemblies/metag/assembly.txt` | shared by every member of a co-assembly group |
| final MAGs | every non-canopy row of `<GC>/Bin_<binner>/MAGvsGC.txt.gz` (or `.txt`; MGS `?` included) | MAG and MGS found by header name; conflicting identities fail |
| MAG contigs | `<assembly>/Binning/<binner>/<SmplID>` (or `.gz`; MAG `<sample or alias>__<bin>`, alias `<SmplID>M<n>` -> SmplID) | GeneCat's own rules (`lib/mf4.py bin_contigs`, shared with the MMT index) |
| sequence | `<assembly>/scaffolds.fasta.filt(.gz)` | the binning input and MATAFILER4's mapping reference |
| circular | contigs listed in `<assembly>/scaffolds.fasta.circ` | MATAFILER4's record; sequences are never re-circularised or trimmed |
| check | `<GC>/Bin_<binner>/Genomes/MGS_ctg/<MGS>.ctgs.<MAG>` where it exists (representatives) | contig set and every base must equal |
| raw BAMs | `<run>/<S>/input_raw.txt` (the files MATAFILER4 read), else the map (`#DirPath` + `<SmplPrefix>*.bam`) | a sample may have several; a BAM that does not open is not used |
| reads | `<run>/<S>/mapping/<S>-smd.cram` (+ `.retained.size`, `-smd.reference.stat`) | the sample's reads aligned to its own (co-)assembly |
| abundance | sample `mapping/<S>-smd.bam.coverage(.gz)`, then `-smd.cram.coverage(.gz)`; fallback sample `assemblies/metag/ContigStats/Coverage.percontig(.gz)` | never co-assembly pooled depth or MGS abundance |

The binner is automatic only when exactly one `Bin_*` contains a MAGvsGC table. With several, **`--binner TAG` is
required**: a configured SB default cannot silently select the wrong result. Every MAG identity and bin membership
uses that binner. Assembly paths are canonicalized; relative pointers resolve against `assembly.txt`.

Relocated or future MF4 coverage layouts can be supplied with `--coverage-map FILE`: a TSV with `sample` and
`coverage` columns, optionally `format` (`bedgraph` or `percontig`; otherwise detected). Paths are relative to the TSV,
or absolute. Unknown/repeated sample IDs and missing paths fail before processing. Overrides are authoritative;
an invalid override is not silently replaced. Discovery stays within the named sample, never recursively choosing
another sample's coverage or a support/co-assembly track.

## Outputs

`OUT.tsv` has row key **`group` + `MAG`**. It contains every discovered final MAG x run sample combination, including
those with no selected molecules. `sample` = `<SmplID>__<MAG>` (MAG without `.fa.gz`; characters outside
`A-Za-z0-9._-` -> `_`; over 100 characters: cut + checksum). `group` = original SmplID; `role` = native.

| Columns | Meaning |
|---|---|
| `sample`, `group`, `MAG`, `MGS`, `binner` | Calling ID, original sample, exact MAG, GeneCat MGS (`?` if none), selected binner |
| `bam`, `reads` | **JSON arrays in every row**, including singleton arrays; index `i` pairs BAM `i` with selection file `i`. Empty arrays `[]` when no inputs are callable. |
| `ref`, `role` | Exported MAG FASTA (empty if unavailable), `native` |
| `mot_status`, `mot_reason` | `ready` rows have validated BAM/selection pairs; otherwise `not_profiled`, `sample_unavailable`, `MAG_unavailable`, or `no_selected_molecules` |
| `normalized_abundance`, `abundance_unit` | Mean MAG depth per Gbp contributing to this sample's assembly depth; unit `fold_coverage_per_Gbp_assembly_mapped` |
| `abundance_status`, `abundance_reason` | `ok`, `not_profiled`, `unavailable`, `invalid_coverage`, `reference_unavailable`, `reference_mismatch`, `MAG_unavailable`, `shared_contigs`, or `zero_sample_depth` |
| `mag_mean_depth`, `mag_breadth_1x`, `mag_depth_sd` | Length-weighted mean depth, fraction of MAG bases at depth >=1, base-level depth SD; zero-depth positions included |
| `mag_mapped_bases`, `sample_mapped_bases` | Sum of interval length x depth within MAG and across the full sample assembly |
| `mag_reference_bp`, `sample_reference_bp` | MAG length and full assembly length |
| `mapped_base_fraction` | MAG depth integral / sample depth integral: a mapped DNA-base fraction, not a cell fraction |
| `coverage_source`, `coverage_format` | Exact evidence file and `bedgraph`/`percontig` format |

Unavailable numeric values are **`NA`, never zero**. An observed MAG with zero depth and a nonzero sample denominator
has numeric zero. A sample mapped only against its own assembly has **no estimate for another assembly's MAG**:
these pairs are `not_profiled`, even when another MAG belongs to the same MGS. MGS abundance is never copied to
individual member MAGs.

### Abundance definition

```text
sample_mapped_bases    = sum(interval_length * depth) across ALL assembly contigs
mag_mapped_bases       = the same sum restricted to this MAG
mag_mean_depth         = mag_mapped_bases / mag_reference_bp
normalized_abundance  = mag_mean_depth * 1e9 / sample_mapped_bases
```

Unbinned contigs and lower-quality bins contribute to the denominator. Fragmented MAGs are length-weighted.
Coverage inherits MF4's existing depth filters, independently of export's `--min-mapq`, `--exclude-qcfail`, BAM
byte sizes, BAM count and selected-molecule counts. Support-library coverage is not mixed into primary coverage.

This is **coverage per Gbp of assembly-mapped sequencing effort**, not per Gbp of all sequenced DNA, a genome
percentage, or physical copies per gram. Unmapped/unassembled sequence is excluded. Assembly recovery and mapping
ambiguity can affect comparisons; existing own-assembly mappings are not cohort-wide competitive quantification.

Full zero-inclusive bedGraph is validated against every assembly contig and length. Gaps, overlaps, unknown
contigs, truncated references, non-finite/negative depths and incomplete tracks yield `NA` and PARTIAL. Contig-mean
fallback requires every assembly contig too; rounding is inherited and breadth/SD remain `NA`. Shared MAG contigs
are not used as independent abundance evidence. Abundance caches include source path/size/nanosecond mtime,
assembly identity, method and MAG membership; changed coverage refreshes abundance without redoing unchanged
molecule assignments.

### Contract for the MOT reader update

1. Read columns by name. Preserve `group` + `MAG` as row key and `MGS` as metadata.
2. Call only `mot_status=ready` rows. Decode `bam` and `reads` as equal-length, nonempty JSON arrays. Process each
   matched pair under the row's single `sample` and `ref`, combining BAM evidence within that sample.
3. Preserve abundance metadata. `abundance_status=ok` identifies usable estimates; `NA`/unprofiled is not absence.
   Molecule inputs can be callable even when abundance is unavailable.
4. Abundance describes the MAG x sample once; do not sum it once per BAM. It is not a methylation probability
   or a reason to alter the kinetic null.

`OUT.motlib/`

| Path | Content |
|---|---|
| `genomes/<MAG>.fa.gz` (+ `.fai`, `.gzi`) | one bgzip FASTA per final MAG: its contigs as MATAFILER4 named them, in assembly order, one record each |
| `assignments/<SmplID>/<MAG>/<BAM>.reads.txt.gz` | the molecules (`movie/zmw`, sorted, unique) of that MAG in that BAM - mot's read selection |
| `catalog/mags.tsv` | MAG, file, MGS, is_MGS_rep, status, sample, assembly, contigs, bp, circular contigs, CheckM, FASTA, md5, checked vs GeneCat |
| `catalog/contigs.tsv.gz` | MAG, contig, assembly, length, md5 (upper case, as CRAM M5), circular (1 / 0 / NA = no `.circ` file) |
| `catalog/samples.tsv` | per sample: status, group, assembly, CRAM, BAMs and their provenance, co-assembled samples, record counts, retained-alignment evidence |
| `catalog/bams.tsv` | per raw BAM: samples, movies, chemistry, read type, kinetics, mot-ready, why not |
| `catalog/molecules/<SmplID>.parquet` | every molecule with a record on a final MAG (columns below) |
| `qc/export_summary.tsv`, `qc/assignment_summary.tsv` (sample x MAG x BAM counts), `qc/unresolved.tsv`, `qc/validation.json` | QC |
| `manifest.json` | inputs, settings, method, status definitions, counts, CRAM stamps |

Molecule catalog columns: sample, MAG, movie, zmw, source_bam, status, bam_status, selected, mapq (primary), primary_mapped,
qc_fail (0x200), duplicate (0x400), n_primary, n_supplementary, n_secondary, contig, pos (primary), other_MAGs,
secondary_MAGs.

## Assignment

MATAFILER4 aligns each sample's reads to its own (co-)assembly. So the CRAM is a competitive mapping **within that
assembly** (every contig of it, binned or not): a molecule's best place among those contigs. It is **not** a cohort-wide
competitive mapping: a sample's CRAM never saw the MAGs of other assemblies, and callable inputs exist only where the CRAM
has evidence (a sample x MAG of its own assembly). Co-assembled samples each have their own CRAM against the shared
assembly: each becomes its own mot sample on the shared MAG FASTA.

Per molecule (`movie/zmw`; by-strand `/fwd` `/rev` records and supplementary records are the same molecule), from its
primary and supplementary records:

| status | meaning | selected |
|---|---|---|
| unique | all on contigs of one final MAG, primary MAPQ >= `--min-mapq` (10) | yes |
| ambiguous_low_mapq | one MAG, primary MAPQ below: a near-equal placement elsewhere in the assembly | no |
| split_across_mags | pieces on two or more MAGs | no |
| split_mag_unbinned | pieces on one MAG and on contigs of no final MAG | no |
| contig_in_two_mags | on a contig two final MAGs claim | no |
| reference_mismatch | on a contig whose CRAM header length / M5 differs from the exported sequence | no |
| no_primary | only supplementary records | no |

Molecules only on contigs of no final MAG (lower-quality bins, unbinned) are counted, not listed. Secondary records
never decide (MAPQ already reflects them) but are listed (`secondary_MAGs`).

**QC-flagged alignments**: with `-retainAllCRAM 1` MATAFILER4 keeps alignments its filter rejected (0x200) and
duplicates (0x400). They are kept in the catalog with their flags and selected when unique - mot aligns and filters
the molecules itself (`--exclude-qcfail` leaves 0x200 out). A CRAM without such records and without
`.retained.size` is reported: rejected alignments that were never written cannot be recovered.

**Source BAM**: a molecule belongs to the sample's raw BAM whose read groups / first records carry its movie. A movie in
none of the sample's BAMs (`movie_not_in_bams`) or in two (`movie_in_two_bams`) is reported and not selected. Names that
are not PacBio (`movie/zmw` with a PacBio movie name, or a movie of the sample's BAMs) - Illumina, ONT - are counted and
never exported. A BAM is mot-ready when it opens, carries kinetics (`fi`/`ri`/`ip`) and READTYPE CCS.

## Status

- **COMPLETE**: every check passed and nothing was left out.
- **PARTIAL**: the table holds every row that validated; `qc/unresolved.tsv` lists what was left out and why: a sample
  with PacBio reads but no raw BAM, a BAM without kinetics, a CRAM made against another version of a contig, a MAG
  that could not be exported (contigs missing, sequence differs from GeneCat), molecules whose BAM cannot be named,
  several chemistries in one mot sample, one BAM in two samples with the same molecules. Samples without PacBio data
  (Illumina, ONT) and assemblies without a final MAG are not problems; they are listed in `catalog/samples.tsv`.
- **FAILED**: a check of the written files failed, or nothing can be exported (no GeneCat table, MMT's naming rules
  disagree with GeneCat's own FASTA, ...).

Checks before COMPLETE / PARTIAL (`qc/validation.json`): mot sample names valid; one selection per calling sample x BAM; one
reference per mot sample; every bam / ref / .fai / reads file exists; every selection non-empty, unique, `movie/zmw`,
its movies in its BAM, its count = the catalog's; selected molecules in the catalog = lines in the selection files; no
molecule in two rows of one BAM; every molecule has one status. The MAG FASTAs were checked base by base against
GeneCat's where GeneCat has one.

## Limits

- Assignment is per (co-)assembly, not cohort-wide (above). A cohort-wide competitive assignment needs a new alignment
  of every sample's reads to all MAGs (or one representative per MGS) - not done here.
- The source BAM of a molecule is named by its movie, not by scanning the BAM for every ZMW; mot reports how many
  listed molecules it found.
- A polished contig of the same length is caught only through the CRAM's M5 (CRAM headers normally carry it).
- The schema-2 table adds MAG / MGS / abundance columns and JSON BAM/read arrays. MOT must support this contract
  before using the table. MMT does not modify MOT itself.
- The full MAG x sample Cartesian table can be large. Unprofiled pairs are explicit `NA` rows, not extra mappings.
