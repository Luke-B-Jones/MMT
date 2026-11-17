# shellcheck shell=bash
# MMT export-mot - a finished MATAFILER4 run + GeneCat -> Methyl-o-Type input table + library. Sourced by the mmt wrapper.

module_help() {
    cat << 'EOF'
mmt export-mot - turn a finished MATAFILER4 run + its GeneCat MGS stage into a Methyl-o-Type (mot) input table

Usage: mmt export-mot -i MF4_RUN -gc GENECAT -o OUT.tsv [options]

  -i, --mf4-run DIR    the MATAFILER4 run folder (#OutPath/#RunID: one folder per sample)
  -gc, --genecat DIR   its GeneCat folder (holds Bin_<binner>/MAGvsGC.txt.gz and LOGandSUB/)
  -o, --out FILE       the mot input table to write; the library goes next to it: OUT.motlib/
Options
  --map FILE           the run's map (default: the maps GeneCat was run with, LOGandSUB/GCmaps.ori, else its copies)
  --binner TAG         GeneCat Bin_<TAG> (default: the only one there; several -> MMT_BINNER if it is one of them)
  --min-mapq N         primary MAPQ for a unique assignment [EXPORTMOT_MIN_MAPQ, 10]
  --exclude-qcfail     do not select molecules MATAFILER4 flagged QC-fail (0x200) (default: selected, flagged in
                       the catalog - mot filters alignments itself)
  -t, --threads N      CRAMs / assemblies read at once [the SLURM allocation, else 4]
  --config FILE        extra settings file

What it writes
  OUT.tsv              sample  bam  ref  role  reads  group - one row per sample x MAG x raw BAM with selected
                       molecules; sample = <SmplID>__<MAG>, group = <SmplID>, reads = that BAM's molecules of that MAG
  OUT.motlib/          genomes/<MAG>.fa.gz (+ .fai .gzi): one FASTA per final MAG (MATAFILER4's own contigs)
                       assignments/<SmplID>/<MAG>/<BAM>.reads.txt.gz: movie/zmw, one per line (mot's read selection)
                       catalog/ mags.tsv contigs.tsv.gz samples.tsv bams.tsv molecules/<SmplID>.parquet (every
                       molecule on a final MAG: status, MAPQ, QC-fail / duplicate flags, source BAM, selected)
                       qc/ export_summary.tsv assignment_summary.tsv unresolved.tsv validation.json ; manifest.json
Assignment is competitive within each sample's own (co-)assembly (MATAFILER4's CRAM), not across the cohort.
A rerun reuses every sample and assembly whose inputs and settings did not change.
Exit status: 0 COMPLETE, 3 PARTIAL (what validated is written; the rest is in qc/unresolved.tsv), 1 FAILED.
Needs python with pysam, pyarrow, numpy and samtools (MMT_ENV_PREFIX / MMT_PYTHON / MMT_SAMTOOLS). Details:
docs/export-mot.md
EOF
}

module_main() {
    local TABLE="" BINNER="" MAPQ=${EXPORTMOT_MIN_MAPQ:-10} QCF="" lib rc=0
    while (( $# )); do
        case $1 in
            -o|--out|-gc|--binner|--min-mapq) [[ $# -ge 2 && -n $2 ]] || usage_die "$1 needs a value" ;;
        esac
        case $1 in
            -o|--out)         TABLE=$2; shift 2; continue ;;      # a file here (common_opt's -o is a folder)
            -gc)              MMT_GENECAT=$2; shift 2; continue ;;
            --binner)         BINNER=$2; shift 2; continue ;;
            --min-mapq)       MAPQ=$2; shift 2; continue ;;
            --exclude-qcfail) QCF=--exclude-qcfail; shift; continue ;;
            -h|--help)        module_help; exit 0 ;;
        esac
        if common_opt "$@"; then shift "$COMMON_SHIFT"; continue; fi
        usage_die "unknown option '$1'"
    done
    [[ -n ${MMT_MF4_RUN:-} ]] || usage_die "-i MF4_RUN is required (the MATAFILER4 run folder)"
    [[ -d $MMT_MF4_RUN ]] || die "MATAFILER4 run folder not found: $MMT_MF4_RUN"
    [[ -n ${MMT_GENECAT:-} ]] || usage_die "-gc GENECAT is required (or MMT_GENECAT in config/site.conf)"
    [[ -d $MMT_GENECAT ]] || die "GeneCat folder not found: $MMT_GENECAT"
    [[ -n $TABLE ]] || usage_die "-o OUT.tsv is required (the mot input table to write)"
    [[ $MAPQ =~ ^[0-9]+$ ]] || usage_die "--min-mapq must be a number (got '$MAPQ')"
    [[ -z ${MMT_MAP:-} || -f $MMT_MAP ]] || die "map not found: $MMT_MAP"
    TABLE=$(abspath "$TABLE"); lib="${TABLE%.tsv}.motlib"
    mkdir -p "$lib" || die "cannot create $lib"
    if [[ -n ${MMT_ENV_PREFIX:-} ]]; then
        [[ -d $MMT_ENV_PREFIX/bin ]] || die "MMT_ENV_PREFIX=$MMT_ENV_PREFIX has no bin/ folder"
        PATH="$MMT_ENV_PREFIX/bin:$PATH"
    fi
    need_tool PYTHON "$MMT_PYTHON" python
    need_tool SAMTOOLS "$MMT_SAMTOOLS" samtools
    "$PYTHON" -c 'import pysam, pyarrow, numpy' 2> /dev/null || die "python $PYTHON cannot import pysam, pyarrow and numpy"
    export REF_PATH=/nonexistent/mmt REF_CACHE=/nonexistent/mmt   # never let htslib fetch a CRAM reference
    resolve_resources
    lock_take "$lib/.lock" || die "another export-mot writes $lib (lock $lib/.lock) - wait for it, or remove the lock if none runs"
    "$PYTHON" "$MMT_HOME/lib/motexport.py" --run "$(abspath "$MMT_MF4_RUN")" --genecat "$(abspath "$MMT_GENECAT")" \
        --table "$TABLE" ${MMT_MAP:+--map "$(abspath "$MMT_MAP")"} ${BINNER:+--binner "$BINNER"} \
        --binner-default "$MMT_BINNER" --min-mapq "$MAPQ" $QCF --threads "$THREADS" --samtools "$SAMTOOLS" \
        --version "$(mmt_version)" || rc=$?
    lock_release "$lib/.lock"
    exit "$rc"
}
