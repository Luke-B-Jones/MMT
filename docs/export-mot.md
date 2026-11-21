# mmt export-mot

A finished MATAFILER4 run and its GeneCat MGS stage become one input table for Methyl-o-Type (`mot`): one FASTA per
final MAG, and per sample x MAG x raw BAM the PacBio molecules MATAFILER4 aligned to that MAG. Nothing in the run is
changed; no BAM or CRAM is copied.

```bash
mmt export-mot -i /path/METH_PearlAge -gc /path/GC_METH_PearlAge -o /path/MOT/mot_input.tsv -t 8
mot call -i /path/MOT/mot_input.tsv -b sprq.motbin -o /path/MOT/calls --check    # mot (0.7+) validates it
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
| final MAGs | every non-canopy row of `<GC>/Bin_<binner>/MAGvsGC.txt.gz` (MGS `?` included) | the bins that passed GeneCat's quality filter; lower bins are not MAGs |
| MAG contigs | `<assembly>/Binning/<binner>/<SmplID>` (MAG `<sample or alias>__<bin>`, alias `<SmplID>M<n>` -> SmplID) | GeneCat's own rules (`lib/mf4.py bin_contigs`, shared with the MMT index) |
| sequence | `<assembly>/scaffolds.fasta.filt(.gz)` | the binning input and MATAFILER4's mapping reference |
| circular | contigs listed in `<assembly>/scaffolds.fasta.circ` | MATAFILER4's record; sequences are never re-circularised or trimmed |
| check | `<GC>/Bin_<binner>/Genomes/MGS_ctg/<MGS>.ctgs.<MAG>` where it exists (representatives) | contig set and every base must equal |
| raw BAMs | `<run>/<S>/input_raw.txt` (the files MATAFILER4 read), else the map (`#DirPath` + `<SmplPrefix>*.bam`) | a sample may have several; a BAM that does not open is not used |
| reads | `<run>/<S>/mapping/<S>-smd.cram` (+ `.retained.size`, `-smd.reference.stat`) | the sample's reads aligned to its own (co-)assembly |
| abundance | `--abundance`, else `<GC>/Bin_<binner>/Annotation/Abundance/MGS.matL0.txt` (MGS x samples, or the transpose) | the table's `abundance` column: GeneCat's MGS abundance in the sample |

The binner is the only `Bin_*` with a MAGvsGC table (several: `--binner`, else `MMT_BINNER` if it is one of them).

## Outputs

`OUT.tsv` - mot 0.7's input table (`mot help input`): `sample genome bam reads ref mgs abundance role`, one row per
sample x MAG x raw BAM that has selected molecules. mot writes `OUT/<sample>/<genome>/`.

| Column | Value |
|---|---|
| `sample` | the SmplID: every MAG of a sample shares one mot sample (mot fits one adapter per sample, on all its genomes) |
| `genome` | the MAG (its file name, without `.fa.gz`); the rows of one sample x MAG from several BAMs are one mot unit |
| `bam` | the raw HiFi BAM with kinetics |
| `reads` | that BAM's molecules of that MAG (`assignments/...reads.txt.gz`); never empty |
| `ref` | the MAG FASTA (`genomes/<MAG>.fa.gz`) |
| `mgs` | the MAG's MGS (GeneCat); `.` for a MAG without one (MGS `?`) |
| `abundance` | GeneCat's abundance of that MGS in that sample (`MGS.matL0.txt`, as written there); `.` without an MGS, or where the table has no value above 0 for the pair (or no table) |
| `role` | `native` |

Names follow mot's rule (characters outside `A-Za-z0-9._-` -> `_`; over 100 characters: cut + checksum).

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
competitive mapping: a sample's CRAM never saw the MAGs of other assemblies, and mot rows are made only where the CRAM
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
  several chemistries in one sample (mot needs one per sample), one BAM in two samples with the same molecules. Samples without PacBio data
  (Illumina, ONT) and assemblies without a final MAG are not problems; they are listed in `catalog/samples.tsv`.
- **FAILED**: a check of the written files failed, or nothing can be exported (no GeneCat table, MMT's naming rules
  disagree with GeneCat's own FASTA, ...).

Checks before COMPLETE / PARTIAL (`qc/validation.json`): mot sample and genome names valid; one row per sample x
genome x BAM; one reference, mgs and abundance per sample x genome; one SmplID per sample name and one MAG per genome
name; one chemistry per sample; abundance a number or `.`; every bam / ref / .fai / reads file exists; every selection non-empty, unique, `movie/zmw`,
its movies in its BAM, its count = the catalog's; selected molecules in the catalog = lines in the selection files; no
molecule in two rows of one BAM; every molecule has one status. The MAG FASTAs were checked base by base against
GeneCat's where GeneCat has one.

## Limits

- Assignment is per (co-)assembly, not cohort-wide (above). A cohort-wide competitive assignment needs a new alignment
  of every sample's reads to all MAGs (or one representative per MGS) - not done here.
- The source BAM of a molecule is named by its movie, not by scanning the BAM for every ZMW; mot reports how many
  listed molecules it found.
- A polished contig of the same length is caught only through the CRAM's M5 (CRAM headers normally carry it).
- `abundance` is GeneCat's per-MGS value, not a per-MAG measurement: two MAGs of one MGS in a sample carry the same
  value. The per-MAG molecule counts are in `qc/assignment_summary.tsv`.
- The table is mot 0.7's (`sample` + `genome`); mot 0.6 does not read it. For mot 0.6, export with MMT 0.8.0.
