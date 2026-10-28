# ctg-meth - methylation of every MAG's own contigs, per sample, on the dereplicated contigs

`mmt ctg-meth -o OUT [--select MGS,...] [--samples S,...] [--derep DIR]`

For every sample and every selected MGS with a MAG in the sample's assembly, ctg-meth takes the sample's PacBio HiFi
reads that MATAFILER4 placed on that MAG, aligns them to the **MAG's own contigs** and calls m6A and m4C with ipdSummary.
It writes:

1. **read set `all`** - every read: the whole-contig methylation of the MAG in that sample;
2. **junction read sets** - reads that show a large rearrangement (another conformation of the genome than the assembly,
   or a contig end continuing elsewhere) are grouped by breakpoint and partner; each conformation with enough molecules
   (ALT), and the as-assembled conformation at the same breakpoint (REF), gets its own methylation calls, from those
   reads only, on the contig flanks that every one of them covers.

Every row is placed on the **derep-ctg blocks** (`mmt derep-ctg`: block, base position on the block, block strand,
instance) and on the MAG's genes, so the same base can be compared across samples, MAGs and conformations.
No SNVs are called: only large rearrangements make read sets.

## Method

1. **Reads.** The sample's read index (`OUT/index/reads/<S>.parquet`, shared with MGScall): each read whose primary
   alignment in MATAFILER4's CRAM lies on a contig of an MGS MAG. QC-failed alignments (0x200) are never used;
   duplicates (0x400) only with `CTGMETH_KEEP_DUP=1`. The reads of the selected MGS are taken from the raw kinetics BAM
   (one pass), split into one pseudo-read per strand (`ccs-kinetics-bystrandify`), and grouped by MAG.
2. **Alignment** of each MAG's pseudo-reads to that MAG's own contigs (`pbmm2 --preset CCS -y CTGMETH_MIN_IDENTITY`,
   gap-compressed identity; supplementary alignments kept). Read set `all` = primary alignments with MAPQ >=
   `CTGMETH_MIN_MAPQ`, no hard clip, >= `CTGMETH_MIN_PCT_ALIGNED` % and >= `CTGMETH_MIN_ALN_BP` bp of the read aligned.
3. **Breakpoints.** For one pseudo-read per molecule (the forward one; the reverse one if the forward one did not align):
   its alignments are cut into **pieces** at every large event - a split alignment, a CIGAR deletion / insertion >=
   `CTGMETH_MIN_SV` bp, or a short inversion that the aligner passed straight through (a stretch with >= 15 differences
   in 50 bp whose read sequence is the reverse complement of the contig, >= 90 % identity). Pieces need MAPQ >=
   `CTGMETH_MIN_MAPQ` and >= 100 aligned bp. Consecutive pieces along the read join two **breakpoint ends**
   (contig, position, side: `R` = the sequence ending at that base continues elsewhere, `L` = the sequence starting
   there comes from elsewhere); a read end with >= 500 unaligned bp joins an end to "unplaced" sequence.
   Ends within `CTGMETH_TOL` bp are one; a **breakpoint** is the gap after a base of a contig.
4. **Read sets.** ALT = one pair of ends (or an end and unplaced sequence, or two ends with an insert of another
   length class) with its molecules; listed from 2 molecules (single molecules are only counted). REF = the molecules
   that span a breakpoint contiguously with >= `CTGMETH_ANCHOR` aligned bp on both sides, at breakpoints where an
   ALT or unplaced set has >= `CTGMETH_MIN_MOL` molecules and >= `CTGMETH_MIN_FRAC` of the molecules there.
   A molecule that shows an ALT conformation at a breakpoint is never in its REF set. **Called** (methylation rows written):
   ALT sets placed in the MAG ("claimable placement": both ends on its contigs) and REF sets, each with
   >= `CTGMETH_MIN_MOL` molecules and >= `CTGMETH_MIN_FRAC` of the breakpoint's molecules. Unplaced sets are listed,
   never called. `CTGMETH_JUNCTIONS=0` lists every set and calls none.
   Event classes: `deletion`, `insertion`, `replacement`, `duplication`, `inversion`, `translocation` (two
   contigs), `contig_join` (two contig ends: the assembly broke there), `circular` (a contig's end continues at its
   start), `unplaced`, `ref`. Flags: `linked:<set>` (>= 80 % of the molecules shared, e.g. the two ends of an
   inversion), `side<k>_not_aligned` (that side is a short inversion the aligner passed through: no rows there),
   `assembly_suspect` (ALT >= 90 % and no REF set: the reads deny the assembled conformation), `fold_back` (an
   inversion of < 300 bp without a partner: possibly a read folding back on itself).
5. **Callable distance.** Per called set and side, n(d) = the molecules whose piece reaches at least d + M bp into
   that side from the breakpoint (M = `CTGMETH_CTX` + 5). D = the d <= `CTGMETH_WINDOW` with n(d) >= `CTGMETH_MIN_MOL`
   that maximises n(d) x d. Only those n(D) molecules are used for that side, and rows are written for 0 <= dist <= D:
   **every read of the side covers every reported base** (within `CTGMETH_TOL` bp of the breakpoint a read's own
   breakpoint may lie a few bp off the set's).
6. **Windows.** Each called side becomes a small reference: the flank (D + M bp) plus, at the breakpoint, the
   `CTGMETH_CTX` + 5 bp that follow it in that conformation (the partner's sequence; an insert as N), so ipdSummary
   predicts the kinetics next to the junction from the right context. Both pseudo-reads of every molecule used get
   exactly one record there: the one covering the base next to the breakpoint, clipped to the window (parts outside ->
   soft clip). Hard-clipped (supplementary) records are rebuilt from the unaligned pseudo-read, kinetics in native
   orientation (PacBio convention, checked on every run).
7. **Kinetics** - one ipdSummary run per sample on every MAG contig (read set `all`) and every window
   (`--identify m6A,m4C --methylFraction --mapQvThreshold CTGMETH_MIN_MAPQ --minCoverage 1 --maxCoverage 100000`),
   plus, from the same reads, `samtools mpileup` for the sequence-agreement columns (as in MGScall).
8. **Lift** onto the derep-ctg blocks (`OUT/derep-ctg/<MGS>/lift.parquet`: exact-match segments between every
   instance and its block sequence) and the MAG's genes (`genes.parquet`). Every lifted row is checked: its base must be
   the block's base at `block_pos` on `block_strand`, else the run stops. A newer derep-ctg result only re-annotates
   the files (re-lift: no kinetics recomputed). The derep-ctg result of an MGS must hold every usable MAG of the MGS in
   the current index; otherwise the MGS is skipped with a warning (rerun `mmt derep-ctg` after adding samples).

## Output

```
OUT/ctg-meth/<MGS>/<S>.parquet            sites: read set 'all', then the junction read sets (long format)
OUT/ctg-meth/<MGS>/junctions/<S>.parquet  one row per read set x side
OUT/ctg-meth/<MGS>/edges/<S>.parquet      every instance end and contig end of the sample's MAGs
OUT/ctg-meth/<MGS>/reads/<S>.parquet      every aligned piece of every molecule used
OUT/ctg-meth/_state/<S>.tsv               per MGS: status, fingerprints, counts ; _state/<S>.stats.tsv read accounting
```

The four files exist for every MGS x sample processed (possibly empty, always with the same columns). Each has the
columns `MGS` and `sample`: `read_parquet('OUT/ctg-meth/MGS.345/*.parquet')` is every sample of an MGS (it does not
read the sub-folders), `read_parquet('OUT/ctg-meth/*/junctions/*.parquet')` every junction table. Positions are signed
32-bit integers everywhere (differences never overflow). Shared inputs of the analyses: the derep-ctg tables of the MGS
(`OUT/derep-ctg/<MGS>/`: blocks, instances, lift, genes, block_genes, families, contigs, contig_seqs, composition) and
the sample's read depth (`OUT/index/depth/<S>.parquet`, 1 kb bins of every MGS-MAG contig).

### Sites `<MGS>/<S>.parquet`

One row per base x strand with an A or C on that strand and >= `CTGMETH_MIN_ROW_COV` kinetic reads (no row = fewer
reads). Rows of read set `all` come first, sorted by block, instance, block_pos, block strand (rows off the blocks:
by MAG, contig, pos, strand), then the junction read sets by read set, side, dist.

| Column | Meaning |
|---|---|
| MGS, sample | the file's MGS and sample |
| read_set | `all`, or a junction read set (`R0001` = REF, `J0001` = ALT; see the junction table) |
| side, dist | read sets: side of the breakpoint (1 / 2) and bp from the breakpoint into that side (0 = the base next to it). NULL for `all` |
| MAG, contig, pos, strand | MAG contig position (1-based) and strand (`+` / `-`) of the base |
| base | A or C on that strand |
| block, instance | derep-ctg block and instance this base lies on (the innermost one: an island inside its host gives the island). NULL: not on a block, or a base that differs from the block |
| block_pos, block_strand | the base's position on the block sequence (`blocks.parquet` seq, 1-based) and strand relative to the block |
| ctx_dist | the largest k such that the k bases on each side equal the block's (0 = a neighbouring base differs, or an indel or the instance's end is next to it), capped at 255. Comparable across instances when ctx_dist >= 10 (`CTGMETH_CTX`) |
| gene | the MAG gene (`genes.parquet` gene, `<contig>_<n>`) holding the base (overlapping genes: the one starting last) |
| up_gene, up_dist | the gene whose start lies 1..300 bp downstream of the base in that gene's orientation (promoter zone), and that distance. One gene: the nearest start; a tie goes to the `+` gene. Can also be set inside an upstream gene's body (`gene` not NULL) |
| coverage, score, ipdRatio, tMean, tErr, modelPrediction, frac, fracLow, fracUp, call, identificationQv | ipdSummary, as in MGScall (`coverage` = kinetic reads of this strand, capped by ipdSummary's subsampling - not depth: see `OUT/index/depth/`). `fracLow` / `fracUp` are bootstrap bounds: they vary slightly between runs; everything else is deterministic |
| pile_depth, n_ref, top_alt, n_alt, ctx_alt_pct, ctx_alt_n, clip_ctx_n | sequence agreement of the same reads with the MAG, as in MGScall |
| kin_fp, lift_fp | fingerprints of the kinetics (one value per sample) and of the derep-ctg result the lift came from (one value per MGS folder when every sample is current: check before comparing samples, see below) |

### Junction table `<MGS>/junctions/<S>.parquet` (one row per read set x side)

| Column | Meaning |
|---|---|
| MGS, sample, MAG, read_set, side | the read set and side (REF: 1 = left of the breakpoint, 2 = right; ALT: one row per placed end). Read-set IDs are numbered per sample (all MGS of the sample) and change when a result is recomputed: join sites on (MGS, sample, read_set, side); key stored results on (MAG, contig, locus, bp_side) for REF sets, and with partner_kind, partner_contig, partner_pos, partner_side, insert_bp for ALT sets (several ALT sets can share a breakpoint end) |
| kind, event | `ref` / `alt` ; event class (method 4) |
| contig, pos, bp_side | this side's breakpoint end: `R` = the anchor is the sequence ending at pos ; `L` = starting at pos |
| locus | the breakpoint: the gap after base `locus` of the contig (ends within `CTGMETH_TOL` bp are one breakpoint). ALT and REF sets of one breakpoint share (MAG, contig, locus): the pairing key |
| partner_kind, partner_contig, partner_pos, partner_side | where the read goes on: `contig` (REF: the contig continues), `placed` (another end in the MAG), `unplaced` |
| insert_bp, event_bp | unaligned bp between the two ends (insertions; N in the window) ; event size (deletion, insertion, duplication, inversion) |
| n_molecules | molecules of the read set |
| locus_molecules, locus_ref | molecules at this side's breakpoint (REF + every ALT set there) ; of those, REF (as assembled; NULL where REF was not counted: below the thresholds of method 4, or a contig end). 1 - locus_ref / locus_molecules = share of the molecules not as assembled, whatever the conformation |
| frac | n_molecules / locus_molecules (largest of the set's two breakpoints) = this set's share. The two junctions of one inversion are two sets sharing breakpoints, so each holds part of the inversion's molecules: use locus_ref for the event's prevalence |
| called, call_bp, n_used, n_rows | rows written for this side ; callable distance D ; molecules used ; rows |
| div_kb | mismatches + small indels per kb of the side's reads inside the window (a high value = another strain rather than another conformation) |
| flags | method 4 |
| molecules | the molecules used for this side (all molecules of an uncalled set, on side 1) |
| anchor_instance, anchor_block, anchor_block_pos, anchor_end, anchor_end_dist | the instance this side runs into: the smallest one holding the breakpoint end (an island included), else the nearest one on the anchor side ; its block ; the breakpoint end's position on the block (NULL when that base is not on the instance's lift) ; the end (5 / 3) of the instance facing the breakpoint and the bp from the breakpoint to that end (<= 0: inside the instance, > 0: outside). An incomplete instance's end is not its block's end (`instances.parquet` start_status / end_status) |
| partner_instance, partner_block, partner_block_pos, partner_end, partner_end_dist | the same for the partner end |
| kin_fp, lift_fp | as in the sites |

### Edges `<MGS>/edges/<S>.parquet`

| Column | Meaning |
|---|---|
| kind | `instance_end` (an end of a derep-ctg instance, islands included) or `contig_end` |
| MAG, contig, pos, bp_side | the end: last (`R`) or first (`L`) base |
| instance, block, block_end | the instance and its block end (5 / 3) |
| nb_key, gap | the MAG's own neighbour at that end (`instances.parquet` nb5 / nb3 key; an island: `host:<block>`) and the bp to it |
| n_span | molecules spanning the end contiguously (>= `CTGMETH_ANCHOR` bp on both sides) |
| n_event | molecules (each counted once) of listed ALT sets with an end facing the same way in the gap after this end |
| n_runoff, n_end | contig ends: molecules running off the contig (>= 500 unaligned bp, or continuing on another piece) / ending there |

`n_event = 0` with a good `n_span` means the sample is homogeneous at that end: its read set `all` rows there belong to
the MAG's own conformation.

### Reads `<MGS>/reads/<S>.parquet`

One row per aligned piece of every molecule (molecule orientation = the raw read): MAG, molecule (`movie/zmw/ccs`),
piece (order along the molecule), contig, r_start, r_end, strand, q_start, q_end, q_len, mapq, supplementary, synthetic
(a short inversion found in a passed-through alignment), ev_start / ev_end (what bounds the piece at r_start / r_end:
`junction`, `overhang`, `end`), div (mismatches + small indels), in_all (the record of the listed pseudo-read - the
forward one, or the reverse one when the forward one did not align - is in read set `all`), read_sets.
For read-level work later (which molecules cover which blocks; per-molecule paths through the blocks).

### State `_state/<S>.tsv`

Per MGS: status (`done`, `absent` = fewer than `CTGMETH_PRESENT_MOL` molecules on the MGS's MAG(s) in the sample:
not processed, no files ; `no_kinetics`, `no_aligned`, `no_rows`), kin_fp, lift_fp, molecules,
pseudo_reads, aligned (pseudo-reads in `all`), rows, read_sets (read sets listed), called (read sets with rows),
updated, mmt, ipdsummary (version used).

## Comparing across samples, MAGs and conformations

Rules that keep comparisons fair (use the same rules for every group compared):

- **One derep-ctg result**: `SELECT count(DISTINCT lift_fp) FROM read_parquet('OUT/ctg-meth/MGS.345/*.parquet')` must
  be 1 (block and instance IDs change when derep-ctg is rerun; `mmt ctg-meth` over all samples re-lifts in minutes and
  warns when an MGS is mixed).
- **Same bases**: compare on (block, block_pos, block_strand) with `ctx_dist >= 10`; that also excludes the
  representative's private differences (block sequences are polished to the majority base of their instances).
- **Coverage**: `frac` exists only for called sites with > 10 reads; compare uncalled sites with `ipdRatio` or `score`,
  stratified by motif (context: recipe 6).
- **Strain mixtures**: `ctx_alt_pct < 10` (no other base in the +-10 bp context among the sample's reads).
- **Junction context**: `dist >= 10` (the kinetics next to a junction depend on what follows it); `clip_ctx_n = 0`
  for read set `all`.
- **Read set `all` near a rearrangement** holds each molecule's primary alignment, so where a sample has another
  conformation it mixes both: for whole-block comparisons use instance ends with `n_event = 0` (edges) or leave out
  the bases within `call_bp` of listed ALT breakpoints. At a circular contig's ends `all` covers about half the
  molecules; the `circular` read set has the same molecules again (do not add the two).
- **Batch**: kinetics depend on chemistry - `OUT/index/movies.tsv` gives each sample's movies, binding / sequencing kit
  and basecaller.
- **Unit**: a sample x MAG; co-assembly MAGs appear in several samples, and repeated time points of one person are not
  independent.
- **IDs**: `instance`, `family` and read-set IDs are per MGS (read sets: per sample); across MGS key on (MGS, ...) and,
  for genes, on `catalogue` (GeneCat catalogue gene number = the key of GeneCat's annotation tables).
- **Strand fractions cannot show hemimethylation**: after replication each daughter duplex has its new strand on
  opposite sides, so population fractions per strand stay symmetric; that needs per-molecule calls (not made yet).

### DuckDB recipes

```sql
-- (1) the whole-contig (block) pattern per sample: one column per sample:instance
PIVOT (SELECT sample || ':' || instance AS col, block_pos, block_strand, ipdRatio
       FROM read_parquet('OUT/ctg-meth/MGS.345/*.parquet')
       WHERE read_set = 'all' AND block = 'B0003' AND ctx_dist >= 10 AND coverage >= 11)
ON col USING first(ipdRatio) GROUP BY block_pos, block_strand ORDER BY block_pos;

-- (1b) does the block's pattern differ between instances (samples / MAGs)?
WITH s AS (SELECT sample, instance, block, block_pos, block_strand, ipdRatio FROM read_parquet('OUT/ctg-meth/MGS.345/*.parquet')
           WHERE read_set = 'all' AND ctx_dist >= 10 AND coverage >= 11)
SELECT a.block, a.sample s1, a.instance i1, b.sample s2, b.instance i2, count(*) n, corr(a.ipdRatio, b.ipdRatio) r
FROM s a JOIN s b USING (block, block_pos, block_strand)
WHERE (a.sample, a.instance) < (b.sample, b.instance) GROUP BY ALL;

-- (1c) position of the instance in its MAG (neighbours, contig, strand) next to its methylation
WITH s AS (SELECT sample, MAG, instance, block, avg(ipdRatio) ipd FROM read_parquet('OUT/ctg-meth/MGS.345/*.parquet')
           WHERE read_set = 'all' AND ctx_dist >= 10 AND coverage >= 11 GROUP BY ALL)
SELECT s.*, i.contig, i.c_start, i.c_end, i.strand, i.contig_length, i.nb5_key, i.nb3_key, i.nb5_mags, i.nb3_mags
FROM s JOIN read_parquet('OUT/derep-ctg/MGS.345/instances.parquet') i USING (instance);

-- (1d) genes of a block across instances (block-sequence frame: block_genes seq_start..seq_end) - a join example;
--      an analysis would stratify by motif
SELECT s.sample, s.instance, g.slot, g.family, count(*) sites, avg(s.ipdRatio) ipd
FROM read_parquet('OUT/ctg-meth/MGS.345/*.parquet') s
JOIN read_parquet('OUT/derep-ctg/MGS.345/block_genes.parquet') g
  ON g.block = s.block AND s.block_pos BETWEEN g.seq_start AND g.seq_end
WHERE s.read_set = 'all' AND s.ctx_dist >= 10 GROUP BY ALL;

-- (2) the 3' end of block B0003 in every conformation, sample and MAG (breakpoints within 1 kb of the instance's 3' end;
--     anchor_end_dist <= 0 = inside the instance)
SELECT j.sample, j.MAG, j.read_set, j.kind, j.event, j.locus, j.partner_block, j.partner_end, j.n_molecules, j.frac,
       s.dist, s.block_pos, s.block_strand, s.base, s.coverage, s.ipdRatio, s.frac
FROM read_parquet('OUT/ctg-meth/MGS.345/junctions/*.parquet') j
JOIN read_parquet('OUT/ctg-meth/MGS.345/*.parquet') s USING (sample, read_set, side)
WHERE j.anchor_block = 'B0003' AND j.anchor_end = 3 AND j.anchor_end_dist BETWEEN -1000 AND 1000 AND j.called AND s.dist >= 10;

-- (3) effect of the partner (what B0003's 3' end is joined to), by distance from the junction ...
WITH x AS (SELECT j.sample, j.MAG, j.kind, coalesce(j.partner_block || ':' || j.partner_end, j.partner_kind) AS partner,
                  s.block_pos, s.block_strand, s.dist, s.base, s.ipdRatio
           FROM read_parquet('OUT/ctg-meth/MGS.345/junctions/*.parquet') j
           JOIN read_parquet('OUT/ctg-meth/MGS.345/*.parquet') s USING (sample, read_set, side)
           WHERE j.anchor_block = 'B0003' AND j.anchor_end = 3 AND j.anchor_end_dist BETWEEN -1000 AND 1000 AND s.dist >= 10)
SELECT partner, kind, dist // 500 * 500 AS dist_bin, count(DISTINCT sample) samples, avg(ipdRatio) FILTER (base = 'A') ipd_A
FROM x GROUP BY ALL ORDER BY partner, dist_bin;

-- (3b) ... and paired within a sample: ALT minus REF at the same breakpoint and flank, on the same bases
WITH j AS (SELECT * FROM read_parquet('OUT/ctg-meth/MGS.345/junctions/*.parquet') WHERE called),
     s AS (SELECT * FROM read_parquet('OUT/ctg-meth/MGS.345/*.parquet') WHERE read_set <> 'all' AND dist >= 10)
SELECT a.sample, a.MAG, a.contig, a.locus, a.bp_side, a.read_set AS alt_set, a.event, count(*) n, avg(sa.ipdRatio - sr.ipdRatio) d_ipd
FROM j a JOIN j r ON r.sample = a.sample AND r.MAG = a.MAG AND r.contig = a.contig AND r.locus = a.locus AND r.bp_side = a.bp_side
                 AND a.kind = 'alt' AND r.kind = 'ref'
JOIN s sa ON sa.sample = a.sample AND sa.read_set = a.read_set AND sa.side = a.side
JOIN s sr ON sr.sample = r.sample AND sr.read_set = r.read_set AND sr.side = r.side
         AND sr.contig = sa.contig AND sr.pos = sa.pos AND sr.strand = sa.strand
GROUP BY ALL;

-- (4) one chromosome frame per MGS: every sample's block rows on the contig coordinates of a frame MAG (e.g. a closed
--     one: a called 'circular' junction, contigs.parquet) through the frame MAG's lift ; blocks the frame MAG holds in
--     several copies (instances.dup_in_mag) would give several positions: left out
WITH f AS (SELECT l.* FROM read_parquet('OUT/derep-ctg/MGS.345/lift.parquet') l
           JOIN read_parquet('OUT/derep-ctg/MGS.345/instances.parquet') i USING (instance)
           WHERE l.MAG = 'S1__SemiBin_1' AND NOT i.dup_in_mag)
SELECT s.sample, f.contig AS frame_contig,
       CASE WHEN f.ori > 0 THEN f.c_start + (s.block_pos - f.b_start) ELSE f.c_end - (s.block_pos - f.b_start) END AS frame_pos,
       CASE WHEN (f.ori > 0) = (s.block_strand = '+') THEN '+' ELSE '-' END AS frame_strand, s.base, s.ipdRatio, s.frac
FROM read_parquet('OUT/ctg-meth/MGS.345/*.parquet') s
JOIN f ON f.block = s.block AND s.block_pos BETWEEN f.b_start AND f.b_start + (f.c_end - f.c_start)
WHERE s.read_set = 'all' AND s.ctx_dist >= 10;

-- (5) cumulative GC skew (ori near its minimum, ter near its maximum on a closed chromosome) and each sample's read
--     depth (growth-rate gradients), in the same 1 kb bins of the MAG contigs ; leave out bins within one read length
--     of a linear contig end (depth falls there)
SELECT c.MAG, c.contig, c.w_start,
       sum((c.g - c.c) / greatest(c.g + c.c, 1)) OVER (PARTITION BY c.contig ORDER BY c.w_start) AS cum_skew,
       d.sample, d.bp / (d.c_end - d.c_start + 1) AS depth
FROM read_parquet('OUT/derep-ctg/MGS.345/composition.parquet') c
LEFT JOIN read_parquet('OUT/index/depth/*.parquet') d ON d.contig = c.contig AND d.c_start = c.w_start;

-- (6) sequence context of any site (+-5 bp, on the site's strand), also off the blocks
SELECT s.sample, s.MAG, s.contig, s.pos, s.strand, s.ipdRatio,
       CASE s.strand WHEN '+' THEN substr(q.seq, s.pos - 5, 11)
            ELSE translate(reverse(substr(q.seq, s.pos - 5, 11)), 'ACGT', 'TGCA') END AS ctx11
FROM read_parquet('OUT/ctg-meth/MGS.345/S1.parquet') s
JOIN read_parquet('OUT/derep-ctg/MGS.345/contig_seqs.parquet') q ON q.contig = s.contig
WHERE s.read_set = 'all' AND s.pos > 5;

-- (7) gene level: position in the gene (1 = first base of the start codon) and coding strand, for profiles of one gene
--     family across MAGs and samples (copies of identical sequence: group by md5 of the gene's contig_seqs substring)
SELECT s.sample, s.MAG, g.gene, g.family, g.catalogue,
       CASE g.strand WHEN '+' THEN s.pos - g.c_start + 1 ELSE g.c_end - s.pos + 1 END AS gene_pos,
       s.strand = g.strand AS coding_strand, s.base, s.ipdRatio
FROM read_parquet('OUT/ctg-meth/MGS.345/*.parquet') s
JOIN read_parquet('OUT/derep-ctg/MGS.345/genes.parquet') g ON g.gene = s.gene          -- gene IDs are unique in an MGS
WHERE s.read_set = 'all';

-- (7b) promoter windows of every gene (both genes of a divergent pair get the shared region; up_gene keeps one)
SELECT g.gene, s.sample, s.pos, s.strand, s.ipdRatio,
       CASE g.strand WHEN '+' THEN g.c_start - s.pos ELSE s.pos - g.c_end END AS upstream_bp
FROM read_parquet('OUT/ctg-meth/MGS.345/*.parquet') s
JOIN read_parquet('OUT/derep-ctg/MGS.345/genes.parquet') g
  ON g.contig = s.contig AND ((g.strand = '+' AND g.c_start - s.pos BETWEEN 1 AND 300)
                           OR (g.strand = '-' AND s.pos - g.c_end BETWEEN 1 AND 300))
WHERE s.read_set = 'all';

-- breadth: share of an instance's bases with a row (each base has one A/C strand)
SELECT s.sample, s.instance, count(*) / any_value(l.bp) AS breadth
FROM read_parquet('OUT/ctg-meth/MGS.345/*.parquet') s
JOIN (SELECT instance, sum(c_end - c_start + 1) bp FROM 'OUT/derep-ctg/MGS.345/lift.parquet' GROUP BY 1) l USING (instance)
WHERE s.read_set = 'all' GROUP BY ALL;
```

Intergenic promoter rows: `gene IS NULL AND up_gene IS NOT NULL`. Gene-relative positions shift with indels between
copies: compare them only between copies of identical sequence, and leave out sites within 10 bp of a gene end (the
flanks differ).

## Settings (`config/mmt.conf`)

`CTGMETH_KEEP_DUP`, `CTGMETH_MIN_IDENTITY` (98), `CTGMETH_MIN_MAPQ` (10), `CTGMETH_MIN_PCT_ALIGNED` (50),
`CTGMETH_MIN_ALN_BP` (1000), `CTGMETH_MIN_ROW_COV` (3), `CTGMETH_CTX` (10), `CTGMETH_CLIP_MIN` (10),
`CTGMETH_JUNCTIONS` (1), `CTGMETH_MIN_SV` (300), `CTGMETH_ANCHOR` (1000), `CTGMETH_TOL` (20), `CTGMETH_MIN_MOL` (5),
`CTGMETH_MIN_FRAC` (0.02), `CTGMETH_WINDOW` (5000). A changed setting recomputes the results it affects; `kin_fp`
records them. Fixed parts of the method (covered by its version): pieces >= 100 bp, contig-end zone 1000 bp, overhang
500 bp, short-inversion scan 15 differences in 50 bp at >= 90 % reverse-complement identity, promoter zone 300 bp.

## Behaviour

- Unit of work = one sample (its raw BAM is read once per run); SLURM array: `--list-samples`, then one sample per task.
- Presence check: before a sample's raw BAM is read, the molecules MATAFILER4 placed on each MGS's MAG(s) there are
  counted (read index); an MGS with fewer than `CTGMETH_PRESENT_MOL` is `absent` (state only, no files) - when all are
  absent the raw BAM is not read.
- Work folders (`<tmp>/mmt/<OUT name>_<key>/ctg-meth/<S>`, one per output folder, so several output folders can share
  one `--tmp`) are removed when the sample ends, done or failed. `--debug` keeps them, and a rerun then continues from
  the last finished step (so does a run after a job was killed outright). Finished steps delete the read copies they
  no longer need. Work-space disk per sample: about three copies of the
  sample's selected kinetics reads at the peak, plus ipdSummary's CSV while it is loaded (~80 bytes per base x strand
  with >= 3 reads: ~30 GB for 200 Mb of covered MAG contigs) - put `--tmp` where that fits.
- Results are redone when their reads, assembly, MAG contigs, settings or MMT method change (`kin_fp`); only
  re-annotated when the derep-ctg result changes (`lift_fp`); dropped when the MGS has no usable MAG any more or no MAG
  in the sample's assembly. An MGS without a current derep-ctg result (missing, or lacking a MAG of the index) is
  skipped (warning) and its results are kept. A run warns when an MGS's samples lie on different derep-ctg results.
- Needs: SMRT Link (ipdSummary, ccs-kinetics-bystrandify, pbmm2, pbindex), samtools >= 1.15, DuckDB, python >= 3.10
  with pysam, pyarrow, numpy, edlib.
- Not done (yet): placement of unplaced overhangs outside the MAG (other contigs of the assembly), per-read methylation
  calls, per-molecule phasing of orientation mixtures.
