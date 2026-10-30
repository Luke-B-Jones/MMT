# shellcheck shell=bash
# ensure_reads S W : the sample's read index (a cache in the work space, next to the MF4 index), shared by every
#   read-based command (MGScall, ctg-meth):
#   <IDX>/reads/<S>.parquet  sample, read, MGS, MAG, is_rep, contig, pos, mapq, dup - one row per read whose primary
#                            alignment in MATAFILER4's CRAM lies on an MGS-MAG contig (QC-fail 0x200 never ;
#                            duplicates 0x400 kept with dup = true: each command filters them by its own setting)
#   <IDX>/reads/_stats/<S>.tsv  read accounting of the CRAM
#   Needs CRAM[S] ASM[S] ASM_MD5[assembly] IDX SAMTOOLS PYTHON DUCKDB THREADS DUCK_MEM ; sets READS_PQ.
#   Several commands / jobs may build it at once: same content, private temp names, last rename wins.
ASSIGN_METHOD=3   # bump when the read index changes -> rebuilt on the next run of any read-based command

ensure_reads() {
    local S=$1 W=$2 dir fp fpf tp tf
    dir="$IDX/reads"; READS_PQ="$dir/$S.parquet"
    fpf="$dir/_fp/$S"
    fp=$(md5str "$ASSIGN_METHOD|${CRAM[$S]}:$(stamp "${CRAM[$S]}")|${ASM[$S]}:${ASM_MD5[${ASM[$S]}]:-none}")
    [[ -s $READS_PQ && $(cat "$fpf" 2> /dev/null) == "$fp" ]] && return 0
    mkdir -p "$dir/_fp" "$dir/_stats" "$W"
    log "$S: assigning reads from ${CRAM[$S]} (assembly ${ASM[$S]})"
    "$SAMTOOLS" view -h -@ 2 -F 0x104 --input-fmt-option required_fields=0x3f "${CRAM[$S]}" \
        | "$PYTHON" "$MMT_HOME/lib/mf4.py" assign --contigs "$IDX/contigs.tsv.gz" --assembly "${ASM[$S]}" \
            --out "$W/assign.tsv" --stats "$W/assign.stats" --depth "$W/depth.tsv"
    tp=$(tmpname "$READS_PQ")
    "$DUCKDB" -bail -c "SET threads=$THREADS; SET memory_limit='$DUCK_MEM';
        COPY (SELECT '$S' AS sample, read, MGS, MAG, is_rep::BOOLEAN AS is_rep, contig, pos::INTEGER AS pos, mapq::UTINYINT AS mapq,
                     dup::BOOLEAN AS dup
              FROM read_csv('$W/assign.tsv', delim='\t', header=true, all_varchar=true) ORDER BY MGS, read)
        TO '$tp' (FORMAT parquet, COMPRESSION zstd)"
    publish "$tp" "$READS_PQ"
    tf=$(tmpname "$dir/_stats/$S.tsv"); cp "$W/assign.stats" "$tf"; publish "$tf" "$dir/_stats/$S.tsv"
    tf=$(tmpname "$fpf"); echo "$fp" > "$tf"; publish "$tf" "$fpf"
    rm -f "$W/assign.tsv" "$W/depth.tsv"
    log "$S: $(awk -F'\t' '$1 == "assigned_reads" {print $2}' "$W/assign.stats") reads on MGS MAGs ($(tr '\n\t' '; ' < "$W/assign.stats"))"
}

# ---- MAG references (ctg-meth, MGScall --mag-mag): a MAG = its contigs in the sample's assembly
asm_fasta() { local f; for f in "$1/scaffolds.fasta.filt" "$1/scaffolds.fasta.filt.gz"; do [[ -s $f ]] && { echo "$f"; return 0; }; done; return 1; }
mag_ref() { echo "$TMPROOT/refs_mag/$(md5str "$1|${MAGREF_ASM:?}")"; }   # a MAG's reference (MAG names may hold any character)

ensure_mag_refs() {   # ensure_mag_refs S ASSEMBLY_FASTA MAG... : TMPROOT/refs_mag/<md5>/{mag.fa,.fai,mag.mmi,stamp} for each MAG
    local S=$1 asmf=$2 mag r st tmp; shift 2
    local -a todo=()
    declare -A want=()
    for mag in "$@"; do
        r=$(mag_ref "$mag")
        st="$asmf:$(stamp "$asmf")|$(awk -F'\t' -v g="$mag" '$2 == g {print $1}' "$W/contigs.tsv" | sort | md5sum | cut -c1-16)"
        want[$mag]=$st
        [[ $(cat "$r/stamp" 2> /dev/null) == "$st" && -s $r/mag.mmi ]] || todo+=("$mag")
    done
    (( ${#todo[@]} )) || return 0
    tmp=$(tmpname "$TMPROOT/refs_mag/build"); rm -rf "$tmp"; mkdir -p "$tmp"
    for mag in "${todo[@]}"; do awk -F'\t' -v g="$mag" -v n="$(md5str "$mag")" '$2 == g {print $1 "\t" g "\t" n}' "$W/contigs.tsv"; done > "$tmp/contigs.tsv"
    "$PYTHON" "$MMT_HOME/lib/ctgmeth.py" refs --fasta "$asmf" --contigs "$tmp/contigs.tsv" --out "$tmp"
    for mag in "${todo[@]}"; do
        r=$(mag_ref "$mag"); local b; b="$tmp/$(md5str "$mag")"
        "$SAMTOOLS" faidx "$b/mag.fa"
        "$PBMM2" index --preset CCS -j "$THREADS" "$b/mag.fa" "$b/mag.mmi"
        echo "${want[$mag]}" > "$b/stamp"
        [[ -d $r && $(cat "$r/stamp" 2> /dev/null) != "${want[$mag]}" ]] && rm -rf "$r"     # an outdated cache entry
        mv -T "$b" "$r" 2> /dev/null || rm -rf "$b"                                          # another job was faster
        [[ -s $r/mag.mmi ]] || die "reference of $mag missing after build ($r)"
    done
    rm -rf "$tmp"
}
