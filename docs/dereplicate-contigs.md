# mmt dereplicate-contigs

Each MGS (GeneCat species bin) as a short list of **blocks**: runs of genes that sit in the same order and orientation
in the MAGs of the MGS. A block is a "dereplicated contig": one representative sequence stands for every place the
block occurs in every MAG (its **instances**), and the order in which blocks follow each other along the contigs
(**links**) shows where MAGs differ by rearrangements, insertions and deletions. The block list is cohort-wide per
MGS; the per-sample (per-MAG) view is the instance table.

Purpose: stable, reusable units between "a gene" and "a MAG" - for methylation of genes, intergenic regions and gene
neighbourhoods across all MAGs of a species, for block presence/order across strains, and as a compact reference.
No reads are used; everything comes from MATAFILER4's assemblies, gene predictions and bins and GeneCat's gene
catalogue and MGS.

## Quick start

```bash
mmt dereplicate-contigs -o OUT --mgs MGS.12,MGS.40     # two MGS
mmt dereplicate-contigs -o OUT                          # every MGS with a usable MAG, one job
mmt dereplicate-contigs -o OUT --dry-run                # what would run
```

SLURM: build the shared preparation once, then run parts of the MGS list as an array (MGS are spread over the parts
by MAG count, so parts take similar time):

```bash
sbatch -c 16 --mem=16G -t 2:00:00 --wrap 'mmt dereplicate-contigs -o OUT --prep-only -t 16'
sbatch --array=1-10 -c 4 --mem=8G -t 4:00:00 --wrap \
  'mmt dereplicate-contigs -o OUT --part ${SLURM_ARRAY_TASK_ID}/10 --tmp $TMPDIR'
mmt dereplicate-contigs -o OUT        # afterwards: nothing left to do, refreshes OUT/dereplicate-contigs/summary.tsv
```

`--list-mgs` prints the MGS that still have work (one per line) for your own job layout. Every run skips MGS whose
result is current, redoes those whose inputs, settings or MMT method changed, and a failed MGS does not stop the
others (rerun the same command to resume; 5 failures in a row stop the run: a problem common to all MGS). Runs that
overlap (a refresh started while the array still runs) wait for each other's MGS instead of doing them twice.

## Inputs

| What | Where (found by MMT's index of the MATAFILER4 run, `OUT/index/`) |
|---|---|
| MGS members, completeness, contamination, representative MAG | GeneCat `Bin_<binner>/MAGvsGC.txt.gz` (Canopy MAGs and MGS `?` are ignored) |
| MAG contigs | `<assembly>/Binning/<binner>/<SmplID>` (co-assemblies: the alias in the map) |
| contig sequences | `<assembly>/scaffolds.fasta.filt[.gz]` |
| genes | prodigal `<assembly>/genePred/genes[.bac].gff[.gz]`, checked against `genes[.bac].shrtHD.fna[.gz]` |
| gene identity across MAGs | GeneCat catalogue `compl.incompl.<DEREP_GC_ID>.fna.clstr.idx` (or `.ctsv.idx`, `.gz`) + `compl.incompl.<id>.prot.faa` |

A MAG is used when the index marks it `ok` (in the map, assembly and bin file found). An MGS is used when it has at
least one such MAG (GeneCat's representative FASTA is not needed). Nothing is written next to MATAFILER4's files.

The preparation (once per index; `OUT/dereplicate-contigs/_prep/`) reads the catalogue index once and keeps, for every
gene on a contig of an MGS MAG: contig, prodigal gene number, coordinates, strand, partial flag and catalogue gene.
It stops when: a catalogue member is missing from the GFF (catalogue from another gene prediction), more than 2 % of
the genes of >= 100 nt (GeneCat's minimum gene length) have no catalogue entry, a MAG contig is missing from the
assembly, a gene is listed twice, the first `DEREP_PREP_CHECK` genes of an assembly do not reproduce MATAFILER4's
gene FASTA from their GFF coordinates (a GFF that belongs to other contigs), or more than 5 % of those (complete,
>= 300 nt) do not resemble the protein of their catalogue gene (a catalogue made from another gene prediction of the
same contigs, e.g. genePred redone after GeneCat - the member names alone cannot show it). Plain, gzip and bgzip
FASTA are read; MMT's own indexes live in `_prep/` (a plain-gzip assembly is copied there, MAG contigs only).
The preparation is rebuilt only when what it reads changes: the MAG contigs and their assemblies, the assemblies'
`scaffolds.fasta.filt` / `genes.gff`, the catalogue files (a touched map or MAGvsGC, a new map column or a `.fai` next
to an assembly redo nothing).

## Method

Per MGS, on the genes of all its MAGs, in contig order:

1. **Gene families.** GeneCat's catalogue already groups genes at 95 % nucleotide identity (cohort-wide). Divergent
   alleles of one gene can still fall into two catalogue genes, which would look like a gene replacement. Two catalogue
   genes are merged into one family when they fill the same slot in different MAGs (identical flanking gene ends and
   relative orientation), never occur together in one MAG, and their proteins are >= `DEREP_FAM_MIN_ID` identical over
   the shorter one (length ratio >= `DEREP_FAM_LEN_RATIO`); up to three rounds (Panaroo's context-constrained family
   collapse, on MATAFILER4's catalogue). Real replacements (different genes in the same slot) stay apart.
2. **Family classes.** Families that cannot anchor a gene order are *bridged* (their neighbours are joined across them):
   *linker* = multi-copy (>= `DEREP_MAX_COPY` copies in a MAG, or > 1 copy in >= `DEREP_DUP_MARGIN` of its MAGs and in
   >= 2 MAGs: transposases, paralogs - PPanGGOLiN's multigenic families); *hub* = >= `DEREP_HUB_PARTNERS` different
   neighbours at both ends with none dominant (mobile elements); *rare* = complete copies in fewer than
   `DEREP_MIN_FAMILY_MAGS` MAGs; *partial_only* = only partial copies. The rest are *anchors*. Partial genes (contig
   ends) and genes without catalogue entry are masked (not used for order). Prevalence is completeness-adjusted:
   MAGs with the family (partial copies count) / sum of the MAGs' completeness, so a gene missing from incomplete MAGs
   is not mistaken for accessory (unknown completeness counts as 100 %).
   A *mixed* MAG (co-binned strains: more than `DEREP_MIXED_MAG` of its catalogue genes in several copies) keeps all
   its genes and instances, but its copies are strains, not paralogs: they make no linker and do not keep alleles
   apart in step 1. Without this one MAG of three strains turns every gene of its MGS into a linker.
3. **Blocks.** Gene-end adjacency graph of the anchors (bridged genes skipped). Two gene ends are joined when each is
   the other's most frequent neighbour (a tie joins nothing), the pair is seen in >= `DEREP_MIN_EDGE_MAGS` MAGs, at
   BOTH ends >= `DEREP_SYNTENY` of the MAGs that have a neighbour there agree (contig ends do not count: a fragmented
   assembly never splits a block), and at most (1 - `DEREP_SYNTENY`) of them carry an insert of more than
   `DEREP_MAX_INSERT_BP` beyond the usual spacing. *Agree* includes a neighbour that is explained by the deletion or
   insertion of up to `DEREP_MAX_SKIP` genes (seen in >= `DEREP_MIN_EDGE_MAGS` MAGs), so gene presence/absence does
   not split blocks, rearrangements do. Blocks are the maximal unbranched paths (unitigs) of joined genes; a ring is
   cut at its weakest junction (`circular` = true). Orientation: the block is read so that its first family id is
   smaller than its last.
   - **Core blocks** are built from anchors with prevalence >= `DEREP_CORE_PREV`, with all other genes bridged: the
     persistent backbone, cut only where MAGs disagree on order or orientation (PPanGGOLiN's persistent genome).
   - **Accessory blocks** are built from the remaining anchors (accessory genes, and core genes left out of a core
     block because their run was too short), cut at core-block genes: recurrent variable regions (islands, prophages,
     plasmids; PPanGGOLiN's regions of genome plasticity / Panaroo's accessory paths).
   - Blocks need >= `DEREP_MIN_BLOCK_GENES` genes and >= `DEREP_MIN_BLOCK_BP` consensus bp. Genes in no block are kept
     in `genes.parquet` (role `unblocked`, `rare`, `linker`, ...).
4. **Instances.** Every run of a block's genes along a contig in one direction, missing at most `DEREP_MAX_SKIP` block
   genes in a row and with no insert of more than `DEREP_MAX_INSERT_BP` beyond the consensus spacing (a larger insert,
   a rearrangement or a contig end starts a new instance; a gene read against the run's direction after missing
   positions starts a new, inverted instance, so a small inversion gives its own instance and two links; a single gene
   flipped in place stays, counted in `n_flipped`). Each end is `complete` (first / last block gene present),
   `contig_end` (genes missing, but the contig ends before they could fit: incomplete, not absent) or `open` (genes
   missing although there was room, or the block goes on along the same contig: a real difference). Another copy of a
   block gene within `DEREP_MAX_INSERT_BP` of the run is counted in `n_dup`; a farther one starts its own instance.
   An instance inside another one's span is nested: an insertion
   (`host_instance`, inserted after block position `host_pos` of the host - PPanGGOLiN's spot). Consecutive top-level
   instances on a contig give the links between block ends.
5. **Representatives.** Per block the instance that is complete, then has the most genes, the fewest genes in the
   wrong orientation and the fewest inserted bp, then the best MAG (completeness - 5 x contamination; unknown = worst),
   then the MGS's representative MAG. Its sequence (first gene start to last gene end, on the block's orientation) is
   the block's sequence in `blocks.fna.gz`.
6. **Identity.** Every other instance is aligned to its block's sequence (minimap2 `-x asm20 -c`): identity = matches /
   aligned bases (substitutions only; inserts and deletions show in `insert_bp` and `aligned_frac`), aligned fraction
   and the alignments themselves (CIGAR: coordinates on an instance can be lifted to the block).

Results are deterministic for given inputs and settings.

## Settings (`DEREP_*` in config/mmt.conf; any change redoes every MGS on the next run)

| Setting | Default | Meaning |
|---|---|---|
| `DEREP_GC_ID` | 95 | GeneCat catalogue identity (file names `compl.incompl.<id>.*`) |
| `DEREP_PREP_CHECK` | 200 | genes per assembly checked against `genes.shrtHD.fna` (0 = no check) |
| `DEREP_FAM_MIN_ID` | 0.7 | protein identity to merge two catalogue genes into one family (0 = never merge; then edlib is not needed) |
| `DEREP_FAM_LEN_RATIO` | 0.8 | shorter / longer protein length for a merge |
| `DEREP_MAX_COPY` | 3 | linker: this many copies in one MAG |
| `DEREP_DUP_MARGIN` | 0.05 | linker: > 1 copy in this share of its MAGs (and in >= 2 MAGs) |
| `DEREP_MIXED_MAG` | 0.1 | a MAG with more than this share of its catalogue genes in several copies is mixed (co-binned strains) |
| `DEREP_HUB_PARTNERS` | 5 | hub: this many different neighbours at both ends (0 = no hubs) |
| `DEREP_MIN_FAMILY_MAGS` | 2 | rare: complete copies in fewer MAGs (capped at the number of MAGs) |
| `DEREP_CORE_PREV` | 0.7 | completeness-adjusted prevalence of core-block families |
| `DEREP_SYNTENY` | 0.9 | share of MAGs that must agree at both gene ends of a join (> 0.5) |
| `DEREP_MIN_EDGE_MAGS` | 2 | MAGs that must show a join (capped at the number of MAGs) |
| `DEREP_MAX_INSERT_BP` | 10000 | extra bp between two block genes that still counts as the same block / instance |
| `DEREP_MAX_SKIP` | 3 | block genes that may be missing in a row (deletions / insertions tolerated by joins and instances) |
| `DEREP_MIN_BLOCK_GENES` | 3 | smallest block, genes |
| `DEREP_MIN_BLOCK_BP` | 2000 | smallest block, consensus bp |
| `DEREP_IDENTITY` | 1 | 1 = instance identity + alignments (minimap2, `MMT_MINIMAP2`); 0 = skip |

Fewer, longer blocks: lower `DEREP_SYNTENY` (e.g. 0.8) or raise `DEREP_MAX_SKIP` / `DEREP_MAX_INSERT_BP`. Small
accessory regions as blocks too: lower `DEREP_MIN_BLOCK_GENES` / `DEREP_MIN_BLOCK_BP`. Stricter core: raise
`DEREP_CORE_PREV`. Keep `DEREP_MIN_FAMILY_MAGS` and `DEREP_MIN_EDGE_MAGS` at 2: at 1, genes and joins seen in a
single MAG count (MAG-private blocks, and private genes are no longer bridged).

## Output (`OUT/dereplicate-contigs/`)

```
README.md                 this file
summary.tsv               one row per MGS (rebuilt at the end of every run)
_prep/                    shared preparation (genes of all MGS MAGs, catalogue links) - rebuilt when inputs change
MGS=<MGS>/                one folder per MGS (replaced as a whole when redone)
  blocks.fna.gz           one sequence per block: >B0001 MGS=.. tier=.. genes=.. len=.. rep=MAG:contig:start-end:strand
  blocks.parquet          one row per block
  block_genes.parquet     one row per block position (gene family)
  instances.parquet       one row per block occurrence in a MAG contig
  alignments.parquet      alignments of the instances to their block's sequence
  links.parquet           block ends that meet along contigs
  genes.parquet           every gene of every MAG of the MGS
  families.parquet        every gene family
  merged_catalogue_genes.parquet   catalogue genes merged into one family
  contigs.parquet         every MAG contig with its block path
  mags.parquet            every MAG with its block coverage
  summary.tsv             the MGS's row of summary.tsv (+ the settings used)
logs/dereplicate-contigs/ one log per run (OUT/logs)
```

Coordinates: contig coordinates are 1-based inclusive. "On the block's orientation": `strand` `+` = the block reads
along the contig; `-` = along the reverse complement, and the instance's sequence is the reverse complement of
`contig[start..end]`. Block IDs (`B0001`, ...) are numbered core first, then by number of genes; they are stable for unchanged
inputs and settings, not across changes.

### blocks.parquet

| Column | |
|---|---|
| MGS, block | IDs |
| tier | `core` / `accessory` |
| n_genes | block positions (gene families) |
| length | bp of the block's sequence (the representative instance) |
| consensus_length | median gene lengths + median spacing |
| circular | the block closed into a ring (cut at its weakest junction) |
| n_instances, n_mags, n_mags_complete | instances; MAGs with an instance; MAGs with a complete instance |
| prevalence, prevalence_adj | MAGs with an instance / MAGs; / sum of MAG completeness |
| n_samples | samples behind those MAGs (a co-assembly counts all its samples) |
| median_identity | median instance identity to the block's sequence |
| shared_frac, max_mgs | share of the block's families whose catalogue genes also occur in MAGs of other MGS; the most MGS sharing one of its genes |
| rep_MAG, rep_contig, rep_start, rep_end, rep_strand | the representative instance |
| rep_complete, rep_is_MGS_rep | it is complete; its MAG is the MGS's representative MAG (GeneCat's starred MAG; when that is a Canopy MAG, the first other MAG) |

### block_genes.parquet

MGS, block, pos (0-based position), family, orientation (of the gene on the block), catalogue (catalogue gene IDs of the
family), n_instances (instances holding this position), median_len, cons_start / cons_end (consensus coordinates,
1-based), rep_gene / rep_start / rep_end (the gene in the representative instance and its coordinates on the block's
sequence, 1-based; NULL when the representative lacks it), shared_mgs.

### instances.parquet

| Column | |
|---|---|
| MGS, instance, block, tier | instance = integer ID within the MGS |
| MAG, samples, contig, start, end, strand | where (first gene start .. last gene end) |
| first_pos, last_pos, n_present, frac | block positions covered; genes present; n_present / n_genes |
| complete | all genes present and both ends complete |
| start_status, end_status | `complete` / `contig_end` / `open` for the block's start and end (see Method 4) |
| n_flipped, n_dup | genes in the other orientation; extra copies of a block gene inside the run |
| insert_bp | bp beyond the consensus spacing, summed over the instance |
| host_instance, host_pos | nested in that instance, after its block position host_pos (-1 = not nested) |
| identity, aligned_frac | to the block's sequence: matches / aligned bases; share of the instance aligned (NULL for the representative itself, or without alignment) |
| is_rep | the block's representative instance |

### alignments.parquet

instance, block, q_start / q_end (0-based, end-exclusive, on the instance's block-oriented sequence), strand, b_start /
b_end (on the block's sequence), matches, aln_len, cigar (minimap2 `cg:Z`, PAF conventions). Several rows per instance
when its alignment is split (inserts, rearranged pieces). Lifting a position x of contig c into block coordinates:
take the instance of the gene holding x from genes.parquet (between genes: the innermost instance whose span holds x),
q = x - start (`+`) or end - x (`-`), then the alignment row with q_start <= q < q_end and walk its CIGAR from
b_start (a `-` row - a piece flipped in that MAG - runs backwards along the instance: from q_end - 1 down); for the
representative instance (is_rep) the block position is q. Positions in a few instance ends beyond large inserts have
no alignment row (aligned_frac < 1).

### links.parquet

block_a, pos_a, end_a and block_b, pos_b, end_b: the two gene ends that face each other between consecutive (not
nested) instances on a contig (end 5 / 3 = the 5' / 3' end of that block position on the block's orientation);
terminal (both are block ends: a junction between whole blocks; otherwise an instance is incomplete there), n_mags,
median_gap (bp between the two instances; negative = overlap), MAGs. An inversion shows as one block joined to another
at two different ends; a translocation as a block end with several partners.

### genes.parquet

Every gene on a contig of an MGS MAG: MAG, contig, gene (`<contig>_<n>`, as in MATAFILER4's gene files), start, end,
strand, partial, catalogue (-1 = none), family, role (`core` / `accessory` = its family is in a block of that tier,
`unblocked` = anchor in no block, `rare`, `linker`, `hub`, `masked` = partial or no catalogue entry), block + pos (when it
is a member of an instance), instance (its own instance, or the innermost instance whose span holds it), slot_after
(a gene inside an instance that is not one of its block genes: the block position just before it on the block's
orientation, -1 = before the first).

### families.parquet

family (= its smallest catalogue gene ID), catalogue, n_catalogue, n_mags (MAGs with a complete copy), prevalence_adj,
max_copies, multi_copy_mags, median_len, cls (`core_anchor`, `anchor`, `rare`, `partial_only`, `linker`, `hub`),
block + pos, shared_mgs (MGS whose MAGs carry one of its catalogue genes).

### merged_catalogue_genes.parquet

cat_a, cat_b, identity: each merge of two families (Method 1): the two catalogue genes whose proteins were compared
(each its family's most common catalogue gene at that point) and their protein identity.

### contigs.parquet and mags.parquet

contigs: MGS, MAG, contig, length, n_genes, n_masked, bp_core / bp_blocks (contig bp inside core / any instances), path
(the blocks along the contig, e.g. `B0003+ [B0012-] B0001-`; brackets = nested). mags: MGS, MAG, samples, is_MGS_rep,
completeness, contamination, bp, bp_core, bp_blocks, frac_core, frac_blocks, n_blocks, mixed (co-binned strains,
Method 2).

### summary.tsv

MGS, status (`done` / `no_blocks`), n_mags, n_contigs, bp (all MAG contigs), n_genes, n_masked, n_mixed, n_families, n_merged,
n_linker, n_hub, n_rare, n_core_anchor, n_blocks_core, n_blocks_accessory, bp_blocks (sum of block sequence lengths =
size of the dereplicated set), frac_bp_core, frac_bp_blocks (share of MAG bp inside core / any instances),
blocks_for_50pct, blocks_for_90pct (blocks, largest first, holding 50 / 90 % of all instance bp), n_instances,
median_identity, n_links, edge_min / family_min (the capped supports used), settings, mmt (MMT version that made it).

### Reading the results

```python
import pandas as pd
b = pd.read_parquet("OUT/dereplicate-contigs/MGS=MGS.12/blocks.parquet")
i = pd.read_parquet("OUT/dereplicate-contigs/MGS=MGS.12/instances.parquet")
i.pivot_table(index="MAG", columns="block", values="frac", aggfunc="max")       # block presence per MAG
```

```sql
-- duckdb: every block of every MGS with its prevalence
SELECT MGS, block, tier, n_genes, length, prevalence_adj
FROM read_parquet('OUT/dereplicate-contigs/MGS=*/blocks.parquet', union_by_name = true);
-- which block holds position 152340 of a contig
SELECT block, strand, start, end, host_instance FROM 'OUT/dereplicate-contigs/MGS=MGS.12/instances.parquet'
WHERE contig = 'PDB31__C12_L=200345=' AND start <= 152340 AND 152340 <= "end" ORDER BY "end" - start;
```

## Behaviour worth knowing

- Only binned contigs of MGS MAGs are used; unbinned contigs and MAGs outside MGS are not.
- Small MGS: supports are capped at the number of MAGs. With one MAG, every contig (of >= `DEREP_MIN_BLOCK_GENES`
  genes) is a block. With two or three MAGs every join must be seen in two of them: blocks stay conservative.
- Genes in several copies in >= `DEREP_DUP_MARGIN` of their MAGs (and in >= 2 MAGs) are bridged: a duplicated
  region is then part of the surrounding block's sequence, but its genes have no block position.
- Two alleles are merged only when both neighbours agree: an allele always next to an extra gene (e.g. an insertion
  element) stays a family of its own.
- Genes of two blocks interleaved along a contig make one instance nested in the other although it holds genes of
  the host; genes.parquet still gives every gene its own instance.
- A circular block (all junctions of a ring agree) is cut at its weakest junction: an instance that crosses the cut
  is reported as two instances.
- The block sequence is one real instance, not a consensus; `rep_complete` = false when no MAG holds the block complete
  (then rep coordinates of the missing positions are NULL in block_genes).
- Very divergent instances (below ~80 % identity) may not align (identity NULL).
- Sequence between blocks (junctions, bridged genes at block ends) belongs to no block: links describe it.

## Resources

Preparation: one job, 15-30 min and ~4 GB for a cohort catalogue of ~15 M genes, parallel over assemblies (`-t`).
Per MGS: one python process (~1 s for most MGS, ~1 GB and a minute or two for an MGS with 150+ MAGs) plus minimap2
with `-t` threads, which takes ~90 % of the time of large MGS (`DEREP_IDENTITY=0` skips it). Each MGS process reads
only its own catalogue proteins and MAG contigs. Temp space (`--tmp`): up to ~1 GB per job (the instances of the
largest MGS). A plain-gzip assembly adds a copy of its MAG contigs to `_prep/` (bgzip or plain: none).
Results are written to a temp folder and renamed into place; a lock per MGS keeps two jobs apart; a killed job leaves
nothing half-written; a job whose preparation was rebuilt by another run meanwhile writes nothing.

## Test

`test/run_derep_test.sh WORKDIR` builds a synthetic MATAFILER4 + GeneCat run with known blocks (`test/derep_mock.py`:
three MGS, strains with an inversion, a translocation, islands, a prophage shared with another MGS, a plasmid,
divergent alleles, a gene replacement, deletions, a transposase, contamination, contig breaks, a co-assembly) and
checks the results against the truth (`test/check_derep.py`), then reruns, selection, parts, settings changes, input
changes (gzip / bgzip files, other catalogue naming, a catalogue that does not fit the genes, re-binning, what redoes
what), locks, failure + resume and a MAG of three co-binned strains. `test/derep_scale.py` makes a large run for timing.
