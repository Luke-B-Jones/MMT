# shellcheck shell=bash
# MMT status - what is done per MGS in an output folder. Sourced by the mmt wrapper. Reads OUT and the index (built if
# missing, in the work space); prints the summary and appends it to the mmt log.

module_help() {
    cat << 'EOF'
mmt status - what is done per MGS in an output folder (printed, and added to the mmt log)

Usage: mmt status -o OUT [--select IDS] [--tmp DIR]

For mgs-rep, mag-rep, mag-mag and ctg-meth: per MGS, the samples it is planned for (a MAG of the MGS in the sample's assembly),
the state of each - done, linked (mag-rep / mgs-rep: the other mode's identical table), no_reads, no_kinetics,
no_aligned, no_rows (ctg-meth: absent), or 'not yet' (no result: not run yet, or failed - see the mmt log) - with the
sample names, and the samples where GeneCat's abundance table has the MGS but no MAG of it was binned. mag-rep and
mag-mag also count their MAG tables. The kinetics models ipdSummary used are listed. derep-ctg: per MGS whether its result exists.
  --select IDS   only these MGS ("MGS.12,MGS.40" or a file) ; default every MGS with any result
Needs the index (rebuilt in the work space --tmp if missing: same inputs as the other commands) ; no SLURM job.
EOF
}

module_main() {
    local SELECT="" a lab f sc m
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
    resolve_inputs
    need_tool PYTHON "$MMT_PYTHON" python
    ensure_index > /dev/null
    [[ -n $SELECT ]] && mapfile -t sel < <(read_ids "$SELECT")
    for a in "mgs-rep|3" "mag-rep|3" "mag-mag|3" "ctg-meth|2"; do
        IFS='|' read -r f sc <<< "$a"
        [[ -d $MMT_OUT/$f/stats ]] || continue
        sf=()
        for m in "$MMT_OUT/$f"/stats/*.tsv; do [[ $m == *.run.tsv || ! -f $m ]] || sf+=("$m"); done
        mapfile -t -O "${#lines[@]}" lines < <(state_summary "$f" "$sc" "${sel[*]:-}" "${sf[@]}")
        lab=$(cat "$MMT_OUT/$f"/stats/*.run.tsv 2> /dev/null | awk -F'\t' '$1 ~ /^kinetics_model/ {n[$2]++} END {for (k in n) printf "%s%s (%d)", (c++ ? ", " : ""), k, n[k]}')
        [[ -z $lab ]] || lines+=("$f | kinetics model(s): $lab")
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

state_summary() {   # state_summary LABEL STATUS_COLUMN "MGS..." STATEFILE... : per MGS, the planned samples and their state
    local lab=$1 sc=$2 sel=$3; shift 3
    awk -F'\t' -v lab="$lab" -v sc="$sc" -v sel="$sel" '
        FILENAME == ARGV[1] { if (FNR > 1 && $2 == "ok") { ns++; ord[ns] = $1; asm[$1] = $5 }; next }
        FILENAME == ARGV[2] { if (FNR > 1) U[$1 "|" $2] = 1; next }
        FILENAME == ARGV[3] { if (FNR > 1) A[$1 "|" $2] = 1; next }
        FNR == 1 { s = FILENAME; sub(/.*\//, "", s); sub(/\.tsv$/, "", s); next }
        {   # mag-rep has one row per MAG: the sample is done when any MAG table exists (counted), else its first state
            k = $1 "|" s; x = $sc
            if (x == "done" || x == "linked") { nt[$1]++; if (!(k in st) || (st[k] != "done" && st[k] != "linked")) st[k] = x }
            else if (!(k in st)) st[k] = x
            if (sel == "") M[$1] = 1
        }
        END {
            if (sel != "") { n = split(sel, xs, " "); for (i = 1; i <= n; i++) M[xs[i]] = 1 }
            no = split("done linked absent no_reads no_kinetics no_aligned no_rows not_yet", O, " ")
            for (m in M) {
                split("", C); split("", L); split("", seen); n = 0; nb = 0; NB = ""
                for (i = 1; i <= ns; i++) {
                    k = ord[i]
                    if (!((asm[k] "|" m) in U)) { if ((k "|" m) in A) { nb++; NB = NB (NB == "" ? "" : ", ") k }; continue }
                    n++; c = ((m "|" k) in st) ? st[m "|" k] : "not_yet"
                    C[c]++; L[c] = L[c] (L[c] == "" ? "" : ", ") k
                }
                line = lab " | " m ": planned for " n " sample(s) (a MAG of it in the sample'"'"'s assembly)"
                sep = " -> "
                for (j = 1; j <= no; j++) if (O[j] in C) { line = line sep O[j] " " C[O[j]]; sep = ", "; seen[O[j]] = 1 }
                for (c in C) if (!(c in seen)) { line = line sep c " " C[c]; sep = ", " }
                if (n == 0) line = line " -> no sample"
                if ((lab == "mag-rep" || lab == "mag-mag") && nt[m]) line = line " ; " nt[m] " MAG table(s)"
                line = line " ; present (GeneCat abundance) but not binned: " nb
                print m "\t0\t" line
                r = 1
                for (j = 1; j <= no; j++) if (O[j] in C) print m "\t" r++ "\t" lab " | " m " " O[j] " " C[O[j]] ": " L[O[j]]
                for (c in C) if (!(c in seen)) print m "\t" r++ "\t" lab " | " m " " c " " C[c] ": " L[c]
                if (nb) print m "\t" r++ "\t" lab " | " m " not binned " nb ": " NB
            }
        }' "$IDX/samples.tsv" "$IDX/units.tsv" "$IDX/abundance.tsv" "$@" | sort -t$'\t' -k1,1V -k2,2n | cut -f3- | sed 's/not_yet/not yet/'
}
