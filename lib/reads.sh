# shellcheck shell=bash
# ensure_reads S W : the sample's read index, shared by every read-based command (MGScall, ctg-meth)
#   OUT/index/reads/<S>.parquet  sample, read, MGS, MAG, is_rep, contig, pos, mapq, dup - one row per read whose primary
#                                alignment in MATAFILER4's CRAM lies on an MGS-MAG contig (QC-fail 0x200 never ;
#                                duplicates 0x400 kept with dup = true: each command filters them by its own setting)
#   OUT/index/depth/<S>.parquet  sample, MGS, MAG, contig, c_start, c_end, bp - read depth of every MGS-MAG contig in 1 kb
#                                bins: reference bp of the primary and supplementary alignments that are neither QC-fail
#                                nor duplicate (supplementary: a read crossing a circular contig's start counts on both
#                                ends) ; depth = bp / (c_end - c_start + 1)
#   Needs CRAM[S] ASM[S] ASM_MD5[assembly] IDX SAMTOOLS PYTHON DUCKDB THREADS DUCK_MEM ; sets READS_PQ.
#   Several commands / jobs may build it at once: same content, private temp names, last rename wins.
ASSIGN_METHOD=3   # bump when the read index changes -> rebuilt on the next run of any read-based command

ensure_reads() {
    local S=$1 W=$2 dir fp fpf tp td tf dp
    dir="$MMT_OUT/index/reads"; READS_PQ="$dir/$S.parquet"; dp="$MMT_OUT/index/depth/$S.parquet"
    fpf="$dir/_fp/$S"
    fp=$(md5str "$ASSIGN_METHOD|${CRAM[$S]}:$(stamp "${CRAM[$S]}")|${ASM[$S]}:${ASM_MD5[${ASM[$S]}]:-none}")
    [[ -s $READS_PQ && -s $dp && $(cat "$fpf" 2> /dev/null) == "$fp" ]] && return 0
    mkdir -p "$dir/_fp" "$dir/_stats" "$MMT_OUT/index/depth" "$W"
    log "$S: assigning reads from ${CRAM[$S]} (assembly ${ASM[$S]})"
    "$SAMTOOLS" view -h -@ 2 -F 0x104 --input-fmt-option required_fields=0x3f "${CRAM[$S]}" \
        | "$PYTHON" "$MMT_HOME/lib/mf4.py" assign --contigs "$IDX/contigs.tsv.gz" --assembly "${ASM[$S]}" \
            --out "$W/assign.tsv" --stats "$W/assign.stats" --depth "$W/depth.tsv"
    tp=$(tmpname "$READS_PQ"); td=$(tmpname "$dp")
    "$DUCKDB" -bail -c "SET threads=$THREADS; SET memory_limit='$DUCK_MEM';
        COPY (SELECT '$S' AS sample, read, MGS, MAG, is_rep::BOOLEAN AS is_rep, contig, pos::INTEGER AS pos, mapq::UTINYINT AS mapq,
                     dup::BOOLEAN AS dup
              FROM read_csv('$W/assign.tsv', delim='\t', header=true, all_varchar=true) ORDER BY MGS, read)
        TO '$tp' (FORMAT parquet, COMPRESSION zstd);
        COPY (SELECT '$S' AS sample, MGS, MAG, contig, c_start::INTEGER AS c_start, c_end::INTEGER AS c_end, bp::BIGINT AS bp
              FROM read_csv('$W/depth.tsv', delim='\t', header=true, all_varchar=true) ORDER BY MGS, MAG, contig, c_start)
        TO '$td' (FORMAT parquet, COMPRESSION zstd)"
    publish "$tp" "$READS_PQ"; publish "$td" "$dp"
    tf=$(tmpname "$dir/_stats/$S.tsv"); cp "$W/assign.stats" "$tf"; publish "$tf" "$dir/_stats/$S.tsv"
    tf=$(tmpname "$fpf"); echo "$fp" > "$tf"; publish "$tf" "$fpf"
    rm -rf "$W/assign.tsv" "$W/depth.tsv" "$dir/sample=$S"          # (MMT <= 0.2: index/reads/sample=<S>/data.parquet)
    log "$S: $(awk -F'\t' '$1 == "assigned_reads" {print $2}' "$W/assign.stats") reads on MGS MAGs ($(tr '\n\t' '; ' < "$W/assign.stats"))"
}
