# shellcheck shell=bash
# kinetics_table BAM REF.fa CONTIGS.txt DIR LOG
#   ipdSummary (m6A + m4C, methylated fraction) and, from the SAME aligned reads, samtools mpileup -> sequence agreement.
#   Taken unchanged in substance from the validated Ref-Based pipeline (2026_PA_RefMeth.sh, 20 Oct 2026 fixes included).
#   Needs: SAMTOOLS IPDSUMMARY DUCKDB THREADS DUCK_MEM MIN_MAPQ CTX CLIP_MIN ; writes DIR/k.duckdb (table k) + DIR/stats.tsv:
#     contig pos(1-based) strand(0 = +, 1 = -) base(on that strand) score tMean tErr modelPrediction ipdRatio coverage
#     frac fracLow fracUp call(m6A/m4C/modified_base/none) identificationQv
#     pile_depth  reads whose kinetics are on this strand at this position (= reads aligned to the OTHER strand:
#                 kineticsTools reports a read's IPDs for its template strand), deletions included, not capped
#     n_ref       of those, reads carrying the reference base
#     top_alt     most common other observation A/C/G/T/N (same strand as base) / del / ins ; NULL when all agree
#     n_alt       reads showing top_alt
#     ctx_alt_pct / ctx_alt_n   within +-CTX bp on the same strand, the position with the highest % of reads that
#                 disagree with the reference (substitution, deletion, insertion, or a soft clip >= CLIP_MIN bp = junction
#                 to other sequence) ; positions with < 2 such reads ignored ; 0/0 = clean context
#     clip_ctx_n  reads of this strand with a clip point within +-CTX bp
#   Why the pileup: ipdSummary uses every aligned base (mismatches included) and predicts the expected IPD from the
#   REFERENCE context, so strain SNPs / indels at or near a site are invisible in its own output.
#   ipdSummary conventions: strand = 1 - read orientation ; a CSV row needs >= 3 reads ; frac only above 10 reads ;
#   --maxCoverage default -1 silently drops one read per site -> set to 100000 (20 Oct fix).
kinetics_table() {
    local bam=$1 ref=$2 ctgs=$3 d=$4 lg=$5 pp n_csv chk n_lt mdiff
    awk -F'\t' 'NR==FNR{k[$1]; next} ($1 in k){print $1 "\t0\t" $2}' "$ctgs" "$ref.fai" > "$d/regions.bed"
    #   pileup strand = ipdSummary's: reverse-aligned reads (, a #) -> strand 0 (+), base as printed ;
    #   forward-aligned reads (. A *) -> strand 1 (-), base complemented
    ( "$SAMTOOLS" mpileup -B -Q 0 -d 1000000 --reverse-del -l "$d/regions.bed" -f "$ref" "$bam" 2> "$d/mpileup.log" \
      | awk 'BEGIN { FS = "\t"; OFS = "\t" }
        $4 == 0 { next }
        {
            s = $5
            gsub(/\^./, "", s)
            gsub(/\$/, "", s)
            if (s !~ /[^.,]/) {
                nf = gsub(/\./, ".", s); nr = gsub(/,/, ",", s)
                if (nr) print $1, $2, 0, nr, nr, "", 0
                if (nf) print $1, $2, 1, nf, nf, "", 0
                next
            }
            split("", F); split("", R); df = dr = mf = mr = 0
            n = length(s); i = 1
            while (i <= n) {
                c = substr(s, i, 1)
                if (c == ".")      { mf++; df++ }
                else if (c == ",") { mr++; dr++ }
                else if (c == "*") { df++; F["del"]++ }
                else if (c == "#") { dr++; R["del"]++ }
                else if (c == "+" || c == "-") {
                    j = i + 1; while (substr(s, j, 1) ~ /[0-9]/) j++
                    L = substr(s, i + 1, j - i - 1) + 0
                    if (c == "+") { if (substr(s, j, 1) ~ /[ACGTN*]/) F["ins"]++; else R["ins"]++ }
                    i = j + L; continue
                }
                else if (c ~ /[ACGTN]/) { df++; F[c]++ }
                else if (c ~ /[acgtn]/) { dr++; R[toupper(c)]++ }
                i++
            }
            if (dr) { top(R); print $1, $2, 0, dr, mr, tb, tn }
            if (df) { top(F); print $1, $2, 1, df, mf, cm(tb), tn }
        }
        function cm(x) { return x == "A" ? "T" : x == "C" ? "G" : x == "G" ? "C" : x == "T" ? "A" : x }
        function top(A,   k, o) {
            tb = ""; tn = 0; split("A C G T N del ins", o, " ")
            for (k = 1; k <= 7; k++) if ((o[k] in A) && A[o[k]] > tn) { tn = A[o[k]]; tb = o[k] }
        }' > "$d/pile.tsv" ) &
    pp=$!
    #   clip points (soft clip >= CLIP_MIN bp) of the same reads, on their kinetics strand
    "$SAMTOOLS" view -F 2308 "$bam" | awk -v m="$CLIP_MIN" 'BEGIN { FS = OFS = "\t" } {
        c = $6; st = (int($2 / 16) % 2) ? 0 : 1; L = R = rl = 0
        if (match(c, /^[0-9]+S/)) L = substr(c, 1, RLENGTH - 1) + 0
        if (match(c, /[0-9]+S$/)) R = substr(c, RSTART, RLENGTH - 1) + 0
        while (match(c, /[0-9]+[MIDNSHP=X]/)) { n = substr(c, RSTART, RLENGTH - 1) + 0; o = substr(c, RSTART + RLENGTH - 1, 1)
            if (o ~ /[MDN=X]/) rl += n; c = substr(c, RSTART + RLENGTH) }
        if (L >= m) print $3, $4, st
        if (R >= m) print $3, $4 + rl - 1, st }' > "$d/clips.tsv"
    "$IPDSUMMARY" "$bam" --reference "$ref" -W "$ctgs" \
        --identify m6A,m4C --methylFraction --mapQvThreshold "$MIN_MAPQ" --minCoverage 1 --maxCoverage 100000 -j "$THREADS" \
        --gff "$d/calls.gff" --csv "$d/bases.csv" --log-file "$lg"
    wait "$pp" || die "pileup (samtools mpileup | awk) failed: $(tail -n 3 "$d/mpileup.log")"
    grep -v '^#' "$d/calls.gff" > "$d/calls.tsv" || true
    "$DUCKDB" -bail "$d/k.duckdb" -c "
        SET threads=$THREADS; SET memory_limit='$DUCK_MEM'; SET temp_directory='$d/duck_tmp'; SET preserve_insertion_order=false;
        CREATE TABLE calls AS
            SELECT seqid, pos, CASE gstrand WHEN '+' THEN 0 ELSE 1 END AS strand, any_value(type) AS call,
                   any_value(TRY_CAST(NULLIF(regexp_extract(attr, 'identificationQv=([0-9.]+)', 1), '') AS FLOAT)) AS identificationQv
            FROM read_csv('$d/calls.tsv', delim='\t', header=false, auto_detect=false,
                 columns={'seqid':'VARCHAR','src':'VARCHAR','type':'VARCHAR','pos':'BIGINT','pend':'BIGINT',
                          'gscore':'VARCHAR','gstrand':'VARCHAR','phase':'VARCHAR','attr':'VARCHAR'})
            GROUP BY seqid, pos, gstrand;
        CREATE TABLE pile AS SELECT * FROM read_csv('$d/pile.tsv', delim='\t', header=false, auto_detect=false,
                 columns={'ctg':'VARCHAR','pos':'BIGINT','strand':'INTEGER','pile_depth':'INTEGER','n_ref':'INTEGER',
                          'top_alt':'VARCHAR','n_alt':'INTEGER'});
        CREATE TABLE clips AS SELECT ctg, pos, strand, count(*)::INTEGER AS cn FROM read_csv('$d/clips.tsv', delim='\t', header=false,
                 auto_detect=false, columns={'ctg':'VARCHAR','pos':'BIGINT','strand':'INTEGER'}) GROUP BY ALL;
        CREATE TABLE ctx AS
            WITH ev AS (SELECT coalesce(p.ctg, c.ctg) AS ctg, coalesce(p.strand, c.strand) AS strand, coalesce(p.pos, c.pos) AS vpos,
                               coalesce(greatest(p.pile_depth - p.n_ref, p.n_alt), 0) + coalesce(c.cn, 0) AS an0,
                               greatest(coalesce(p.pile_depth, 0), coalesce(c.cn, 0)) AS depth, coalesce(c.cn, 0) AS cn
                        FROM pile p FULL JOIN clips c ON p.ctg = c.ctg AND p.pos = c.pos AND p.strand = c.strand),
                 e2 AS (SELECT ctg, strand, vpos, least(an0, depth) AS an, depth, cn FROM ev WHERE an0 >= 2 OR cn >= 1)
            SELECT ctg, strand, vpos + o AS pos,
                   coalesce(max(floor(100.0 * an / depth)) FILTER (WHERE an >= 2), 0)::UTINYINT AS ctx_alt_pct,
                   coalesce(arg_max(an, floor(100.0 * an / depth)) FILTER (WHERE an >= 2), 0)::INTEGER AS ctx_alt_n,
                   sum(cn)::INTEGER AS clip_ctx_n
            FROM e2, range(-$CTX, $CTX + 1) t(o) GROUP BY ALL;
        CREATE TABLE k AS
            SELECT c.refName AS contig, c.tpl AS pos, c.strand, c.base, c.score, c.tMean::FLOAT AS tMean, c.tErr::FLOAT AS tErr,
                   c.modelPrediction::FLOAT AS modelPrediction, c.ipdRatio::FLOAT AS ipdRatio, c.coverage::INTEGER AS coverage,
                   c.frac::FLOAT AS frac, c.fracLow::FLOAT AS fracLow, c.fracUp::FLOAT AS fracUp,
                   coalesce(g.call, 'none') AS call, g.identificationQv,
                   p.pile_depth, p.n_ref, p.top_alt, p.n_alt,
                   coalesce(x.ctx_alt_n, 0)::INTEGER AS ctx_alt_n, coalesce(x.ctx_alt_pct, 0)::UTINYINT AS ctx_alt_pct,
                   coalesce(x.clip_ctx_n, 0)::INTEGER AS clip_ctx_n
            FROM read_csv('$d/bases.csv', header=true,
                 types={'refName':'VARCHAR','tpl':'BIGINT','strand':'INTEGER','base':'VARCHAR','score':'INTEGER',
                        'tMean':'DOUBLE','tErr':'DOUBLE','modelPrediction':'DOUBLE','ipdRatio':'DOUBLE','coverage':'BIGINT',
                        'frac':'DOUBLE','fracLow':'DOUBLE','fracUp':'DOUBLE'}) c
            LEFT JOIN calls g ON g.seqid = c.refName AND g.pos = c.tpl AND g.strand = c.strand
            LEFT JOIN pile p ON p.ctg = c.refName AND p.pos = c.tpl AND p.strand = c.strand
            LEFT JOIN ctx x ON x.ctg = c.refName AND x.pos = c.tpl AND x.strand = c.strand;
        DROP TABLE calls; DROP TABLE pile; DROP TABLE ctx; DROP TABLE clips;"
    # integrity: every CSV row exactly once and with pileup counts (same reads -> always)
    n_csv=$(( $(wc -l < "$d/bases.csv") - 1 ))
    chk=$("$DUCKDB" -readonly "$d/k.duckdb" -noheader -list -c "SELECT count(*) || ' ' || count(*) FILTER (WHERE pile_depth IS NULL) FROM k")
    [[ $chk == "$n_csv 0" ]] || die "kinetics table rows / rows without pileup = $chk, ipdSummary CSV rows = $n_csv"
    #   the pileup counts the reads whose IPDs ipdSummary used -> pile_depth >= coverage (where ipdSummary did not
    #   subsample, < 1500 reads) ; anything else = strand / read-set mismatch -> stop
    chk=$("$DUCKDB" -readonly "$d/k.duckdb" -noheader -list -c "SELECT count(*) FILTER (WHERE pile_depth < coverage AND coverage < 1490)
        || ' ' || coalesce(round(avg(pile_depth - coverage) FILTER (WHERE coverage < 1490), 4), 0) FROM k")
    read -r n_lt mdiff <<< "$chk"
    printf 'rows_pile_lt_cov\t%s\nmean_pile_minus_cov\t%s\nclip_points\t%s\n' "$n_lt" "$mdiff" "$(wc -l < "$d/clips.tsv")" >> "$d/stats.tsv"
    (( n_lt * 1000 <= n_csv )) || die "pile_depth < ipdSummary coverage in $n_lt of $n_csv rows (read sets differ)"
    rm -rf "$d/bases.csv" "$d/calls.gff" "$d/calls.tsv" "$d/pile.tsv" "$d/mpileup.log" "$d/regions.bed" "$d/clips.tsv" "$d/duck_tmp"
}
