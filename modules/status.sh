# shellcheck shell=bash
# MMT status - what is done per MGS in an output folder. Sourced by the mmt wrapper. Reads OUT only (no inputs, no
# tools); prints the summary and appends it to OUT/mmt.log.

module_help() {
    cat << 'EOF'
mmt status - what is done per MGS in an output folder (printed, and added to OUT/mmt.log)

Usage: mmt status -o OUT [--select IDS]

For MGScall (each mode) and ctg-meth: per MGS, the samples it is planned for (a MAG of the MGS in the sample's
assembly; mag-rep: the representative's own assembly) and the state of each - done, absent (fewer molecules than
*_PRESENT_MOL), no_kinetics, no_aligned, no_rows, or 'not yet' (no result: not run yet, or failed - see mmt.log) -
with the sample names. derep-ctg: per MGS whether its result exists.
  --select IDS   only these MGS ("MGS.12,MGS.40" or a file) ; default every MGS with any result
Reads only OUT/index and the _state tables: quick, needs no SLURM job.
EOF
}

module_main() {
    local SELECT="" a lab kind f m
    local -a lines=() sf=() sel=()
    while (( $# )); do
        if common_opt "$@"; then shift "$COMMON_SHIFT"; continue; fi
        case $1 in
            --select|--mgs) [[ $# -ge 2 && -n $2 ]] || usage_die "$1 needs MGS IDs or a file"; SELECT=$2; shift ;;
            -h|--help) module_help; exit 0 ;;
            *) usage_die "unknown option '$1'" ;;
        esac
        shift
    done
    resolve_out
    IDX="$MMT_OUT/index"
    [[ -s $IDX/samples.tsv && -s $IDX/units.tsv ]] || die "no index in $MMT_OUT/index - run 'mmt check -o $MMT_OUT' (or any command) first"
    [[ -n $SELECT ]] && mapfile -t sel < <(read_ids "$SELECT")
    for a in "MGScall/mgs-rep|MGScall --mgs-rep|all" "MGScall/mag-rep|MGScall --mag-rep|rep" "ctg-meth|ctg-meth|all"; do
        IFS='|' read -r f lab kind <<< "$a"
        [[ -d $MMT_OUT/$f/_state ]] || continue
        sf=()
        for m in "$MMT_OUT/$f"/_state/*.tsv; do [[ $m == *.stats.tsv || ! -f $m ]] || sf+=("$m"); done
        mapfile -t -O "${#lines[@]}" lines < <(state_summary "$lab" "$kind" "${sel[*]:-}" "${sf[@]}")
    done
    if [[ -d $MMT_OUT/derep-ctg ]]; then
        local -a dm=() dn=()
        if (( ${#sel[@]} )); then dm=("${sel[@]}")
        else mapfile -t dm < <(find "$MMT_OUT/derep-ctg" -mindepth 1 -maxdepth 1 -type d -name 'MGS*' -printf '%f\n' | sort -V); fi
        for m in "${dm[@]}"; do [[ -s $MMT_OUT/derep-ctg/$m/summary.tsv ]] && dn+=("$m"); done
        (( ${#dm[@]} == 0 )) || lines+=("derep-ctg | done $(id_list "${dn[@]}") ; of ${#dm[@]} MGS")
    fi
    (( ${#lines[@]} )) || lines=("no results in $MMT_OUT yet")
    printf '%s\n' "${lines[@]}"
    MLOG_CMD="status"; MLOG_ON=1
    mlog STATUS "${lines[@]}"
}

state_summary() {   # state_summary LABEL KIND "MGS..." STATEFILE... : per MGS, the planned samples and their state
    local lab=$1 kind=$2 sel=$3; shift 3
    awk -F'\t' -v lab="$lab" -v kind="$kind" -v sel="$sel" '
        FILENAME == ARGV[1] { if (FNR > 1 && $2 == "ok") { ns++; ord[ns] = $1; asm[$1] = $5 }; next }
        FILENAME == ARGV[2] { if (FNR > 1 && (kind != "rep" || $4 == 1)) U[$1 "|" $2] = 1; next }
        FNR == 1 { s = FILENAME; sub(/.*\//, "", s); sub(/\.tsv$/, "", s); next }
        { st[$1 "|" s] = $2; if (sel == "") M[$1] = 1 }
        END {
            if (sel != "") { n = split(sel, x, " "); for (i = 1; i <= n; i++) M[x[i]] = 1 }
            no = split("done absent no_kinetics no_aligned no_rows no_reads not_yet", O, " ")
            for (m in M) {
                split("", C); split("", L); split("", seen); n = 0
                for (i = 1; i <= ns; i++) {
                    k = ord[i]
                    if (!((asm[k] "|" m) in U)) continue
                    n++; c = ((m "|" k) in st) ? st[m "|" k] : "not_yet"
                    C[c]++; L[c] = L[c] (L[c] == "" ? "" : ", ") k
                }
                line = lab " | " m ": planned for " n " sample(s) (" (kind == "rep" ? "its representative" : "a MAG of it") " in the sample'"'"'s assembly)"
                sep = " -> "
                for (j = 1; j <= no; j++) if (O[j] in C) { line = line sep O[j] " " C[O[j]]; sep = ", "; seen[O[j]] = 1 }
                for (c in C) if (!(c in seen)) { line = line sep c " " C[c]; sep = ", " }
                if (n == 0) line = line " -> no sample"
                print m "\t0\t" line
                r = 1
                for (j = 1; j <= no; j++) if (O[j] in C) print m "\t" r++ "\t" lab " | " m " " O[j] " " C[O[j]] ": " L[O[j]]
                for (c in C) if (!(c in seen)) print m "\t" r++ "\t" lab " | " m " " c " " C[c] ": " L[c]
            }
        }' "$IDX/samples.tsv" "$IDX/units.tsv" "$@" | sort -t$'\t' -k1,1V -k2,2n | cut -f3- | sed 's/not_yet/not yet/'
}
