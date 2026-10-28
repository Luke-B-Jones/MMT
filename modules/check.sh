# shellcheck shell=bash
# MMT check - settings in effect, tools, inputs. Sourced by the mmt wrapper. Prints to the terminal (no log).

module_help() {
    cat << 'EOF2'
mmt check - show the settings in effect and check tools and inputs

Usage: mmt check [-o OUT] [options]

  Prints the settings after config/mmt.conf -> config/site.conf -> $MMT_CONFIG -> --config -> options, then checks
  every tool (path, version, required features) and the MATAFILER4 / GeneCat inputs.
  With -o OUT it also builds (or refreshes) OUT/index/ and prints its report: samples, MGS and MAGs MMT can use.
  Runs no analysis. Takes the same input options as MGScall (-i, --genecat, --map, --config, -t, --mem).
EOF2
}

module_main() {
    local v
    while (( $# )); do
        if common_opt "$@"; then shift "$COMMON_SHIFT"; continue; fi
        case $1 in
            -h|--help) module_help; exit 0 ;;
            *) usage_die "unknown option '$1'" ;;
        esac
    done
    echo "MMT $(mmt_version) in $MMT_HOME"
    echo "config files: $MMT_HOME/config/mmt.conf$([[ -f $MMT_HOME/config/site.conf ]] && echo " + config/site.conf" || echo " (no config/site.conf)")${MMT_CONFIG:+ + $MMT_CONFIG}"
    echo
    echo "== settings"
    for v in MMT_MF4_RUN MMT_GENECAT MMT_MAP MMT_BINNER MMT_OUT MMT_TMP MMT_SMRT_BIN MMT_ENV_PREFIX MMT_PYTHON MMT_SAMTOOLS \
             MMT_DUCKDB MMT_MINIMAP2 MMT_THREADS MMT_MEM_GB; do
        printf '  %-16s %s\n' "$v" "${!v:-(empty)}"
    done
    while read -r v; do printf '  %-24s %s\n' "$v" "${!v}"; done < <(compgen -v MGSCALL_; compgen -v DEREP_; compgen -v CTGMETH_)
    resolve_resources
    echo "  -> resources: $THREADS threads, $MEM_GB GB${SLURM_JOB_ID:+ (SLURM job $SLURM_JOB_ID)}"
    echo
    echo "== tools"
    setup_tools
    tool_versions | sed 's/^/  /'
    echo "  all tools found, samtools has mpileup --reverse-del, python has pysam + pyarrow + numpy"
    "$PYTHON" -c 'import sys, pyarrow, edlib; assert sys.version_info >= (3, 10)' 2> /dev/null \
        && echo "  derep-ctg / ctg-meth: python >= 3.10 with pyarrow + edlib" \
        || echo "  WARNING: derep-ctg and ctg-meth need python >= 3.10 with pyarrow + edlib ($PYTHON)"
    command -v "${MMT_MINIMAP2:-minimap2}" > /dev/null && echo "  minimap2   $(command -v "${MMT_MINIMAP2:-minimap2}")" \
        || echo "  WARNING: minimap2 '${MMT_MINIMAP2:-minimap2}' not found (derep-ctg with DEREP_IDENTITY=1)"
    echo
    echo "== inputs"
    resolve_inputs
    echo "  map          $MMT_MAP ($(grep -vc '^#' "$MMT_MAP" || true) lines without #)"
    echo "  MF4 run      $MMT_MF4_RUN ($(find "$MMT_MF4_RUN" -mindepth 1 -maxdepth 1 -type d | wc -l) folders)"
    echo "  GeneCat MGS  $MVG ($(zcat "$MVG" | tail -n +2 | wc -l) MAG rows)"
    if [[ -n ${MMT_OUT:-} ]]; then
        resolve_out
        MLOG_CMD="check"; MLOG_ON=1          # an index (re)built here is noted in OUT/mmt.log
        echo
        echo "== index ($MMT_OUT/index)"
        ensure_index > /dev/null
        sed 's/^/  /' "$IDX/report.txt"
    else
        echo
        echo "(add -o OUT to build the index and see which samples / MGS are usable)"
    fi
}
