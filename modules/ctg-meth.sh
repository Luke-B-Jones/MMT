# shellcheck shell=bash
# MMT ctg-meth - methylation of every MAG's own contigs, per sample, lifted onto the derep-ctg blocks, with read sets
# per junction (conformation). Sourced by the mmt wrapper. Method, outputs, columns: docs/ctg-meth.md.
# Engine: lib/ctgmeth.py ; kinetics: lib/kinetics.sh ; read index: lib/reads.sh.
# shellcheck source=../lib/kinetics.sh
source "$MMT_HOME/lib/kinetics.sh"
# shellcheck source=../lib/reads.sh
source "$MMT_HOME/lib/reads.sh"

CTGMETH_METHOD=1   # bump when a code change alters the kinetics / read sets -> the next run recomputes every result
LIFT_METHOD=1      # bump when the block annotation (lift) changes -> the next run re-lifts every result (no kinetics)

module_help() {
    cat << 'EOF'
mmt ctg-meth - methylation of every MAG's own contigs per sample, lifted onto the derep-ctg blocks

Usage: mmt ctg-meth -o OUT [options]          (needs mmt derep-ctg results in OUT/derep-ctg, or --derep DIR)

For every sample and every selected MGS with a MAG in the sample's assembly: the sample's PacBio HiFi reads that
MATAFILER4 placed on that MAG (its CRAM), aligned to the MAG's own contigs, ipdSummary m6A + m4C from the kinetics:
  read set 'all'    every read -> the whole-contig methylation of the MAG in that sample
  junction sets     reads that show a large rearrangement (split alignment, CIGAR indel >= CTGMETH_MIN_SV, short
                    inversion, contig end joined elsewhere) are clustered by breakpoint and partner; for each
                    conformation (ALT) with enough molecules, and the as-assembled conformation (REF) at the same
                    breakpoint, methylation is called from those reads only, on the contig flanks every one of them
                    covers ("callable distance")
Every row is placed on the derep-ctg blocks (block, block_pos, block strand, instance, distance to the nearest sequence
difference) and on the MAG's genes. No SNV calling: only large rearrangements make read sets.

Inputs (defaults from config/site.conf)
  -i, --mf4-run DIR / --genecat DIR / --map FILE   as for every command
  --derep DIR         derep-ctg results (default OUT/derep-ctg) ; an MGS without a current result there is skipped
Selection
  --select IDS        only these MGS ("MGS.12,MGS.40" or a file) ; --mgs is the same
  --samples IDS       only these samples ; HPC: one sample per array task (each raw BAM is read once per run)
Other
  --list-samples      print the samples that have work (one per line, nothing else) and stop - for SLURM arrays
  --dry-run           show what would run and stop
  --force             redo the selected results even if up to date
  --debug             keep each sample's work folder (success or failure; a rerun resumes from it) - without it the
                      work folder is removed when the sample ends, whatever happened (--keep-tmp: the same)
  -t, --threads N / --mem GB / --tmp DIR / -o, --out DIR / --config FILE   as for every command
Settings: CTGMETH_* in config/mmt.conf. A changed setting redoes the results it affects ; a newer derep-ctg result
only re-annotates them (no kinetics recomputed).
Presence check: before a sample's raw BAM is read, the molecules MATAFILER4 placed on each MGS's MAG(s) are counted;
an MGS with fewer than CTGMETH_PRESENT_MOL (default 1) is recorded as 'absent' (no files) and not processed.
Every run is summarised in the mmt log OUT/mmt_<date>.log (plan per MGS, each sample x MGS outcome, failures);
'mmt status -o OUT'. Run logs only with -l DIR (DIR/ctg-meth/<sample>.log).

Output (in OUT/ctg-meth)
  <MGS>/<S>.parquet            sites: read set 'all' then the junction read sets (long format, column read_set)
  <MGS>/junctions/<S>.parquet  one row per read set x side: breakpoint, partner, molecules, callable distance, blocks
  <MGS>/edges/<S>.parquet      every instance end and contig end: molecules spanning it / with a breakpoint in its gap
  <MGS>/reads/<S>.parquet      every aligned piece of every molecule used (placement, read sets)
  stats/<S>.tsv                per MGS: status, fingerprints, counts ; stats/<S>.run.tsv: the sample's last run
Read all samples of an MGS: SELECT * FROM read_parquet('OUT/ctg-meth/MGS.12/*.parquet')

Examples
  mmt ctg-meth -o OUT --select MGS.12 --samples S1          # one MGS, one sample
  mmt ctg-meth -o OUT --select MGS.12,MGS.40 --list-samples > todo.txt
  sbatch --array=1-$(wc -l < todo.txt)%20 -c 16 --mem=64G -t 1-00:00 --wrap \
    'S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" todo.txt); mmt ctg-meth -o OUT --select MGS.12,MGS.40 --samples "$S" --tmp "$TMPDIR"'
EOF
}

ctgmeth_settings() {   # SETTINGS (fingerprints, logs), CM_JSON (lib/ctgmeth.py) ; checks every value
    local k v txt="" json=""
    local -a keys=(keep_dup min_identity min_mapq min_pct_aligned min_aln_bp min_row_cov ctx clip_min junctions min_sv anchor tol
                   min_mol min_frac window)
    for k in "${keys[@]}"; do
        v=CTGMETH_${k^^}; v=${!v:-}
        [[ $v =~ ^(0|[1-9][0-9]*)(\.[0-9]+)?$ ]] || usage_die "CTGMETH_${k^^}='$v' is not a number"
        [[ $k == min_frac || $v =~ ^[0-9]+$ ]] || usage_die "CTGMETH_${k^^}='$v' must be a whole number"
        txt+="${txt:+ }$k=$v"; json+="${json:+,}\"$k\":$v"
    done
    (( CTGMETH_MIN_IDENTITY <= 100 && CTGMETH_MIN_PCT_ALIGNED <= 100 )) || usage_die "CTGMETH_MIN_IDENTITY and _MIN_PCT_ALIGNED are percentages"
    (( CTGMETH_MIN_ROW_COV >= 3 )) || usage_die "CTGMETH_MIN_ROW_COV must be >= 3 (ipdSummary writes rows from 3 reads)"
    (( CTGMETH_MIN_MOL >= 3 && CTGMETH_ANCHOR >= 100 && CTGMETH_WINDOW >= 100 && CTGMETH_MIN_SV >= 50 )) \
        || usage_die "CTGMETH_MIN_MOL >= 3, CTGMETH_ANCHOR >= 100, CTGMETH_WINDOW >= 100 and CTGMETH_MIN_SV >= 50 are required"
    awk -v f="$CTGMETH_MIN_FRAC" 'BEGIN {exit !(f <= 1)}' || usage_die "CTGMETH_MIN_FRAC must be within 0..1"
    SETTINGS=$txt
    # fixed parts of the method (CTGMETH_METHOD covers them): smallest piece, contig-end zone, overhang, dense-region scan
    CM_JSON="{$json,\"min_piece\":100,\"end_tol\":1000,\"overhang\":500,\"dense_n\":15,\"dense_w\":50,\"inv_id\":0.9}"
}

module_main() {
    local SELECT="" SAMPLES="" FORCE=0 DRY=0 LIST=0 t0 S msg n_ok=0 n_skip=0 n_units=0 n_cur=0 err desc
    local -a failed=() todo_samples=()
    MMT_DEBUG=${MMT_DEBUG:-0}
    DEREP=""
    while (( $# )); do
        if common_opt "$@"; then shift "$COMMON_SHIFT"; continue; fi
        case $1 in
            --select|--mgs) [[ $# -ge 2 && -n $2 ]] || usage_die "$1 needs MGS IDs or a file"; SELECT=$2; shift ;;
            --samples) [[ $# -ge 2 && -n $2 ]] || usage_die "--samples needs sample IDs or a file"; SAMPLES=$2; shift ;;
            --derep)   [[ $# -ge 2 && -n $2 ]] || usage_die "--derep needs a folder"; DEREP=$2; shift ;;
            --force)    FORCE=1 ;;
            --dry-run)  DRY=1 ;;
            --list-samples) LIST=1 ;;
            -h|--help)  module_help; exit 0 ;;
            *) usage_die "unknown option '$1'" ;;
        esac
        shift
    done
    (( LIST )) && MMT_QUIET=1
    [[ $CTGMETH_PRESENT_MOL =~ ^[0-9]+$ ]] || usage_die "CTGMETH_PRESENT_MOL='$CTGMETH_PRESENT_MOL' must be a whole number"
    t0=$(date +%s)
    ctgmeth_settings
    resolve_out
    job_tag
    local lname
    if [[ -n $SAMPLES && ! -f $SAMPLES && $SAMPLES != *[,\ ]* ]]; then lname=$SAMPLES
    else lname="$(date +%Y%m%d-%H%M%S)_$JOBTAG$( (( LIST )) && echo _list || true)"; fi
    start_log ctg-meth "$lname"
    say "MMT $(mmt_version) ctg-meth -> $MMT_OUT/ctg-meth${MMT_LOG_DIR:+ (log: $MMT_LOG)}"
    log "command: mmt ctg-meth ${SELECT:+--select $SELECT }${SAMPLES:+--samples $SAMPLES }(force=$FORCE dry=$DRY list=$LIST debug=$MMT_DEBUG) on $(hostname), job $JOBTAG"
    desc="${SELECT:+--select $SELECT }${SAMPLES:+--samples $SAMPLES }"
    (( FORCE )) && desc+="--force "
    (( LIST )) && desc+="(listing the samples with work) "
    (( MMT_DEBUG )) && desc+="--debug "
    (( DRY )) || mlog_start "ctg-meth" "${desc% }"
    trap 'run_exit $?' EXIT
    trap 'exit 143' TERM INT
    resolve_inputs
    resolve_resources
    setup_tools
    "$PYTHON" -c 'import sys, edlib, numpy, pyarrow; assert sys.version_info >= (3, 10)' 2> /dev/null \
        || die "python $PYTHON: ctg-meth needs python >= 3.10 with edlib, numpy and pyarrow"
    DUCK_MEM="$(( MEM_GB * 6 / 10 > 1 ? MEM_GB * 6 / 10 : 1 ))GB"
    SORT_T=$(( THREADS < 4 ? THREADS : 4 ))
    SORT_MEM="$(( MEM_GB * 1024 * 3 / 10 / SORT_T > 256 ? MEM_GB * 1024 * 3 / 10 / SORT_T : 256 ))M"
    MIN_MAPQ=$CTGMETH_MIN_MAPQ CTX=$CTGMETH_CTX CLIP_MIN=$CTGMETH_CLIP_MIN
    log "inputs: MF4 run $MMT_MF4_RUN | GeneCat $MMT_GENECAT (Bin_$MMT_BINNER) | map $MMT_MAP"
    log "resources: $THREADS threads, $MEM_GB GB (DuckDB $DUCK_MEM, sort $SORT_T x $SORT_MEM) | work space $WROOT"
    log "settings: $SETTINGS present=$CTGMETH_PRESENT_MOL | method $CTGMETH_METHOD / lift $LIFT_METHOD"
    tool_versions | sed 's/^/tool: /'
    IPD_VERSION=$("$IPDSUMMARY" --version 2>&1 | head -1 | tr '\t' ' ' || true)

    ensure_index
    D="$MMT_OUT/ctg-meth"; STATE="$D/stats"
    DEREP=$(abspath "${DEREP:-$MMT_OUT/derep-ctg}")
    mkdir -p "$STATE"
    pick_mgs
    pick_samples
    log "selection: ${#MGS_SEL[@]} MGS, ${#SAMPLE_SEL[@]} samples ; derep-ctg results: $DEREP"
    declare -A P_TODO=() P_CUR=()
    local m
    for S in "${SAMPLE_SEL[@]}"; do
        plan_sample "$S"
        n_units=$(( n_units + ${#TODO_KIN[@]} + ${#TODO_LIFT[@]} )); n_cur=$(( n_cur + N_CURRENT ))
        for m in "${TODO_KIN[@]}" "${TODO_LIFT[@]}"; do P_TODO[$m]+="$S "; done
        for m in "${CURRENT[@]}"; do P_CUR[$m]+="$S "; done
        if (( ${#TODO_KIN[@]} + ${#TODO_LIFT[@]} + ${#STALE[@]} )); then todo_samples+=("$S"); else n_skip=$(( n_skip + 1 )); fi
        msg="plan $S: ${#TODO_KIN[@]} to compute, ${#TODO_LIFT[@]} to re-lift, $N_CURRENT up to date${STALE[*]:+, to drop: ${STALE[*]}}"
        if (( DRY )); then say "$msg"; else log "$msg"; fi
    done
    say "plan: ${#todo_samples[@]} sample(s) with work ($n_units MGS x sample results), $n_skip sample(s) up to date, $n_cur result(s) current"
    if [[ -z $SAMPLES ]]; then plan_log; fi
    if (( LIST )); then
        (( ${#todo_samples[@]} == 0 )) || printf '%s\n' "${todo_samples[@]}" >&3
        mlog END "listed ${#todo_samples[@]} sample(s) with work ($n_units MGS x sample results) - nothing processed by a listing run"
        RUN_ENDED=1; return 0
    fi
    if (( DRY )); then say "dry run - nothing done"; RUN_ENDED=1; return 0; fi
    for S in "${todo_samples[@]}"; do
        isolated process_sample "$S"
        if (( RC == 0 )); then n_ok=$(( n_ok + 1 )); else
            failed+=("$S"); err=$(last_error)
            log "== $S FAILED (exit $RC)"
            mlog FAIL "$S: ${err:-exit $RC} - $(log_ref)$( (( MMT_DEBUG )) && echo " ; work folder kept: $WROOT/ctg-meth/$S")"
        fi
    done
    mixed_lift
    mlog END "$n_ok sample(s) processed, ${#failed[@]} failed$( (( ${#failed[@]} )) && echo " ($(id_list "${failed[@]}"))"), $n_skip already up to date ; $(( ($(date +%s) - t0) / 60 )) min"
    RUN_ENDED=1
    if (( ${#failed[@]} )); then
        say "finished with errors after $(( ($(date +%s) - t0) / 60 )) min: ${#failed[@]} sample(s) failed: ${failed[*]} - rerun the same command to resume"
        exit 1
    fi
    say "done in $(( ($(date +%s) - t0) / 60 )) min: $n_ok sample(s) processed, $n_skip already up to date"
}

# ============================================================================ selection + plan
pick_mgs() {   # MGS_SEL = selected MGS with a usable MAG and a current derep-ctg result ; USABLE[MGS] ; IS_SEL[MGS] ; DFP[MGS]
    declare -gA USABLE=() IS_SEL=() DFP=()
    declare -ga MGS_SEL=()
    local m want=() unk=() nod=()
    declare -A known=()
    while IFS=$'\t' read -r m _; do known[$m]=1; done < <(tail -n +2 "$IDX/mgs.tsv")
    while IFS=$'\t' read -r _ m _ st _; do [[ $st == ok ]] && USABLE[$m]=1; done < <(tail -n +2 "$IDX/mags.tsv")
    if [[ -n $SELECT ]]; then
        mapfile -t want < <(read_ids "$SELECT")
        (( ${#want[@]} )) || usage_die "--select: no MGS IDs in '$SELECT'"
        for m in "${want[@]}"; do [[ -n ${known[$m]:-} ]] || unk+=("$m"); done
        (( ${#unk[@]} == 0 )) || die "--select: not in MAGvsGC: ${unk[*]}"
    else
        mapfile -t want < <(printf '%s\n' "${!USABLE[@]}" | sort -V)
    fi
    local cand=() old=() lst
    for m in "${want[@]}"; do
        [[ -n ${USABLE[$m]:-} ]] || { say "WARNING: $m has no usable MAG (see $IDX/report.txt) - skipped"; continue; }
        if [[ -s $DEREP/$m/.fp && -s $DEREP/$m/lift.parquet && -s $DEREP/$m/mags.parquet && -s $DEREP/$m/summary.tsv ]]; then cand+=("$m")
        else nod+=("$m"); fi
    done
    # the derep-ctg result must hold every usable MAG of the MGS in the index (else rows of a new MAG stay off the blocks)
    declare -A dmag=()
    if (( ${#cand[@]} )); then
        lst=$(tmpname "$STATE/derep_mgs"); printf '%s\n' "${cand[@]}" > "$lst"
        "$PYTHON" "$MMT_HOME/lib/ctgmeth.py" derep_mags --derep "$DEREP" --mgs "$lst" > "$lst.mags" \
            || { rm -f "$lst" "$lst.mags"; die "reading the MAG lists of the derep-ctg results in $DEREP failed"; }
        while IFS=$'\t' read -r m g; do dmag[$m|$g]=1; done < "$lst.mags"
        rm -f "$lst" "$lst.mags"
        declare -A miss=()
        while IFS=$'\t' read -r g m _ st _; do
            [[ $st == ok && -n ${USABLE[$m]:-} && -z ${dmag[$m|$g]:-} ]] && miss[$m]=1
        done < <(tail -n +2 "$IDX/mags.tsv")
        for m in "${cand[@]}"; do
            if [[ -n ${miss[$m]:-} ]]; then old+=("$m"); else DFP[$m]=$(< "$DEREP/$m/.fp"); MGS_SEL+=("$m"); IS_SEL[$m]=1; fi
        done
    fi
    (( ${#nod[@]} == 0 )) || say "WARNING: no derep-ctg result for ${#nod[@]} MGS (${nod[*]:0:10}) in $DEREP - run 'mmt derep-ctg' for them; skipped"
    (( ${#old[@]} == 0 )) || say "WARNING: the derep-ctg result of ${#old[@]} MGS (${old[*]:0:10}) lacks MAGs of the current index (made before samples were added?) - rerun 'mmt derep-ctg' for them; skipped, their results kept"
    (( ${#MGS_SEL[@]} )) || die "no MGS selected with a usable MAG and a current derep-ctg result"
}

pick_samples() {   # SAMPLE_SEL ; RAW[S] CRAM[S] ASM[S] ; UNITS[S] = selected MGS with a MAG in its assembly ;
                   # UMD5[asm|MGS] (MAG contigs) ; UMAGS[asm|MGS] ; ASM_MD5[asm]
    declare -gA RAW=() CRAM=() ASM=() UNITS=() UMD5=() UMAGS=() ASM_MD5=()
    declare -ga SAMPLE_SEL=() OKS=()
    local s st raw cram asm m mags mall want=() bad=() unk=()
    declare -A seen=() au=()
    while IFS=$'\t' read -r asm m mags _ mall _; do
        [[ -n ${IS_SEL[$m]:-} ]] || continue
        au[$asm]+="$m "; UMD5[$asm|$m]=$mall; UMAGS[$asm|$m]=$mags
    done < <(tail -n +2 "$IDX/units.tsv")
    while IFS=$'\t' read -r asm mall; do ASM_MD5[$asm]=$mall; done < <(tail -n +2 "$IDX/assemblies.tsv")
    while IFS=$'\t' read -r s st raw cram asm; do
        seen[$s]=$st
        [[ $st == ok ]] || continue
        OKS+=("$s"); RAW[$s]=$raw; CRAM[$s]=$cram; ASM[$s]=$asm; UNITS[$s]=${au[$asm]:-}
    done < <(tail -n +2 "$IDX/samples.tsv")
    if [[ -n $SAMPLES ]]; then
        mapfile -t want < <(read_ids "$SAMPLES")
        (( ${#want[@]} )) || usage_die "--samples: no sample IDs in '$SAMPLES'"
        for s in "${want[@]}"; do
            if [[ -z ${seen[$s]:-} ]]; then unk+=("$s")
            elif [[ ${seen[$s]} != ok ]]; then bad+=("$s(${seen[$s]})")
            else SAMPLE_SEL+=("$s"); fi
        done
        (( ${#unk[@]} == 0 )) || die "--samples: not in the map: ${unk[*]}"
        (( ${#bad[@]} == 0 )) || say "WARNING: skipped, sample not usable: ${bad[*]} (see $IDX/report.txt)"
    else
        while IFS=$'\t' read -r s st _; do
            [[ $st == ok ]] && [[ -n ${UNITS[$s]:-} || -f $STATE/$s.tsv ]] && SAMPLE_SEL+=("$s")
        done < <(tail -n +2 "$IDX/samples.tsv")
    fi
    return 0
}


plan_sample() {   # plan_sample S -> TODO_KIN (compute), TODO_LIFT (re-lift only), CURRENT, STALE (drop), N_CURRENT ;
                  # KFP / LFP[MGS]
    local S=$1 base m st kf lf
    declare -gA KFP=() LFP=()
    declare -ga TODO_KIN=() TODO_LIFT=() STALE=() CURRENT=()
    declare -A ok=() okl=() old_st=()
    N_CURRENT=0
    base="$CTGMETH_METHOD|$ASSIGN_METHOD|$SETTINGS|present=$CTGMETH_PRESENT_MOL|ipd=${MMT_IPD_ARGS:-}|${CRAM[$S]}:$(stamp "${CRAM[$S]}")|${RAW[$S]}:$(stamp "${RAW[$S]}")"
    base+="|$(asm_fasta "${ASM[$S]}" || echo none):$(stamp "$(asm_fasta "${ASM[$S]}" || echo none)")"
    if [[ -f $STATE/$S.tsv ]]; then
        while IFS=$'\t' read -r m st kf lf _; do ok[$m]=$kf; okl[$m]=$lf; old_st[$m]=$st; done < <(tail -n +2 "$STATE/$S.tsv")
    fi
    for m in ${UNITS[$S]:-}; do
        KFP[$m]=$(md5str "$base|${UMD5[${ASM[$S]}|$m]}")
        LFP[$m]=$(md5str "$LIFT_METHOD|${DFP[$m]}")
        if (( ! FORCE )) && [[ ${ok[$m]:-} == "${KFP[$m]}" && ${old_st[$m]:-} == absent ]]; then   # absent: no files to lift
            N_CURRENT=$(( N_CURRENT + 1 )); CURRENT+=("$m")
        elif (( FORCE )) || [[ ${ok[$m]:-} != "${KFP[$m]}" ]] || ! files_ok "$m" "$S"; then TODO_KIN+=("$m")
        elif [[ ${okl[$m]:-} != "${LFP[$m]}" ]]; then TODO_LIFT+=("$m")
        else N_CURRENT=$(( N_CURRENT + 1 )); CURRENT+=("$m"); fi
    done
    for m in "${!ok[@]}"; do      # results of MGS no longer usable, or selected but without a MAG in this assembly now
        [[ " ${UNITS[$S]:-} " == *" $m "* ]] && continue
        if [[ -z ${USABLE[$m]:-} || -n ${IS_SEL[$m]:-} ]]; then STALE+=("$m"); fi
    done
    return 0
}

mixed_lift() {   # warn for selected MGS whose results (all samples) were lifted from different derep-ctg results
    local f m mixed=() sf=()
    for f in "$STATE"/*.tsv; do [[ $f == *.run.tsv || ! -f $f ]] || sf+=("$f"); done
    (( ${#sf[@]} )) || return 0
    while read -r m; do [[ -n ${IS_SEL[$m]:-} ]] && mixed+=("$m"); done < <(awk -F'\t' 'FNR > 1 && $2 == "done" && !s[$1 FS $4]++ {n[$1]++}
        END {for (m in n) if (n[m] > 1) print m}' "${sf[@]}" | sort -V)
    (( ${#mixed[@]} == 0 )) || say "WARNING: ${#mixed[@]} MGS (${mixed[*]:0:10}) have samples lifted from different derep-ctg results: run 'mmt ctg-meth' over all samples (re-lift only, no kinetics) before comparing samples"
}

files_ok() { [[ -f $D/$1/$2.parquet && -f $D/$1/junctions/$2.parquet && -f $D/$1/edges/$2.parquet && -f $D/$1/reads/$2.parquet ]]; }

# ============================================================================ one sample
process_sample() {
    local S=$1 LOCK m ts; LOCK="$STATE/$S.lock"
    ts=$(date +%s)
    lock_take "$LOCK" || die "$S: another job is working on this sample ($(cat "$LOCK/owner" 2> /dev/null)) - remove $LOCK if not"
    W="$WROOT/ctg-meth/$S"                 # this sample's work folder: ours while we hold the lock
    # whatever ends this sample (done, error, scancel / time limit): work folder removed (kept with --debug), lock released
    trap 'rc=$?; work_done "'"$W"'" $(( rc == 0 )); lock_release "'"$LOCK"'"' EXIT
    trap 'exit 143' TERM INT
    declare -gA MOL=() PSE=() ALN=() ROWS=() NSETS=() NCALL=()
    plan_sample "$S"                       # again, under the lock
    if (( ${#TODO_KIN[@]} + ${#TODO_LIFT[@]} == 0 )); then
        if (( ${#STALE[@]} )); then W="" write_state "$S" "$ts"; else log "== $S: up to date (another job finished it)"; fi
        return 0
    fi
    log "== $S: ${#TODO_KIN[@]} MGS to compute, ${#TODO_LIFT[@]} to re-lift | raw ${RAW[$S]} | CRAM ${CRAM[$S]}"
    (( ${#TODO_KIN[@]} )) && compute_sample "$S"
    for m in "${TODO_LIFT[@]}"; do
        "$PYTHON" "$MMT_HOME/lib/ctgmeth.py" relift --out "$D" --sample "$S" --mgs "$m" --derep "$DEREP" --kin-fp "${KFP[$m]}" \
            --lift-fp "${LFP[$m]}" --derep-fp "${DFP[$m]}" --settings "$CM_JSON" --mags "${UMAGS[${ASM[$S]}|$m]}" > /dev/null \
            || die "$S $m: re-lift failed"
        log "$S: $m re-lifted onto the current derep-ctg result"
    done
    write_state "$S" "$ts"
    log "== $S: done"                      # (work folder: removed by the EXIT trap, kept with --debug)
}

compute_sample() {
    local S=$1 t m mag n wfp mof dupf asmf want found
    wfp=$(md5str "${TODO_KIN[*]}|$(for m in "${TODO_KIN[@]}"; do printf '%s ' "${KFP[$m]}"; done)")
    [[ $(cat "$W/work.fp" 2> /dev/null) == "$wfp" ]] && (( ! FORCE )) || rm -rf "$W"   # other work than last time, or --force: start clean
    [[ -d $W ]] && log "$S: resuming from the work folder of an earlier run: $W"
    mkdir -p "$W"; echo "$wfp" > "$W/work.fp"
    rm -f "$W/written.ok"                 # writing is cheap: always redone (the files may have been removed since)
    [[ -f $W/run.stats ]] || : > "$W/run.stats"
    asmf=$(asm_fasta "${ASM[$S]}") || die "$S: no scaffolds.fasta.filt(.gz) in ${ASM[$S]}"
    MAGREF_ASM=$asmf                      # MAG reference cache key: MAG + its assembly (OUTs may share --tmp)

    # 1. read -> MAG (MATAFILER4's CRAM) ; 2. presence + the reads of the MGS to compute, and the contigs of their MAGs:
    #    an MGS with fewer than CTGMETH_PRESENT_MOL molecules on its MAG(s) in this sample is absent, not processed
    t=$(date +%s); ensure_reads "$S" "$W"; log "$S: reads index ok ($(( $(date +%s) - t )) s)"
    if [[ ! -f $W/select.ok ]]; then
        printf '%s\n' "${TODO_KIN[@]}" > "$W/todo.txt"
        dupf=""; (( CTGMETH_KEEP_DUP )) || dupf="AND NOT dup"
        "$DUCKDB" -bail -c "SET threads=$THREADS; SET memory_limit='$DUCK_MEM';
            CREATE TEMP TABLE s AS SELECT read, MAG, MGS FROM read_parquet('$READS_PQ')
                WHERE MGS IN (SELECT m FROM read_csv('$W/todo.txt', header=false, columns={'m':'VARCHAR'})) $dupf;
            COPY (SELECT t.m AS MGS, count(s.read) AS molecules FROM read_csv('$W/todo.txt', header=false, columns={'m':'VARCHAR'}) t
                  LEFT JOIN s ON s.MGS = t.m GROUP BY t.m ORDER BY t.m) TO '$W/molecules.tsv' (HEADER, DELIMITER '\t');
            COPY (SELECT read, MAG, MGS FROM s WHERE MGS IN (SELECT MGS FROM s GROUP BY MGS HAVING count(*) >= $CTGMETH_PRESENT_MOL)
                  ORDER BY read) TO '$W/sel.tsv' (HEADER, DELIMITER '\t')"
        tail -n +2 "$W/sel.tsv" | cut -f1 > "$W/names.txt"
        printf 'selected_molecules\t%s\n' "$(wc -l < "$W/names.txt")" >> "$W/run.stats"
        awk -F'\t' -v k="$CTGMETH_PRESENT_MOL" 'NR > 1 && $2 >= k {print $1}' "$W/molecules.tsv" > "$W/present.txt"
        gzip -dc "$IDX/contigs.tsv.gz" | awk -F'\t' -v a="${ASM[$S]}" 'NR == FNR {t[$1]; next} $1 == a && ($4 in t) {print $2 "\t" $3 "\t" $4}' \
            "$W/present.txt" - > "$W/contigs.tsv"
        log "$S: presence (molecules on the MGS's MAGs in this sample, >= $CTGMETH_PRESENT_MOL = present): $(awk -F'\t' -v k="$CTGMETH_PRESENT_MOL" \
            'NR > 1 {printf "%s %d%s  ", $1, $2, ($2 >= k ? "" : " ABSENT")}' "$W/molecules.tsv")"
        [[ -s $W/names.txt ]] || log "$S: every MGS to compute is absent - raw BAM not read"
        touch "$W/select.ok"
    fi
    local -a present=()
    mapfile -t present < "$W/present.txt"

    # 3. those reads from the raw kinetics BAM -> one pseudo-read per strand -> one BAM per MAG
    if [[ -s $W/names.txt && ! -f $W/extract.ok ]]; then
        t=$(date +%s)
        "$SAMTOOLS" view -@ "$THREADS" -b -N "$W/names.txt" -o "$W/sel.bam" "${RAW[$S]}"
        want=$(wc -l < "$W/names.txt"); found=$("$SAMTOOLS" view -@ "$THREADS" -c "$W/sel.bam")
        printf 'raw_reads_found\t%s\n' "$found" >> "$W/run.stats"
        (( found * 100 >= want * 99 )) || die "$S: only $found of $want selected reads are in ${RAW[$S]} - read names differ between CRAM and raw BAM?"
        "$BYSTRANDIFY" "$W/sel.bam" "$W/bys.bam"
        rm -f "$W/sel.bam"; rm -rf "$W/byMAG"
        tail -n +2 "$W/sel.tsv" | cut -f1,2 | sed '1i read\tMAG' > "$W/sel_mag.tsv"
        mof=$(( $(ulimit -n) - 100 )); (( mof > 2000 )) && mof=2000; (( mof < 50 )) && mof=50
        "$PYTHON" "$MMT_HOME/lib/mf4.py" split --bam "$W/bys.bam" --reads "$W/sel_mag.tsv" --outdir "$W/byMAG" --threads "$SORT_T" --max-open "$mof"
        rm -f "$W/bys.bam"
        log "$S: $found molecules extracted and split into $(( $(wc -l < "$W/byMAG/counts.tsv") - 1 )) MAGs ($(( $(date +%s) - t )) s)"
        touch "$W/extract.ok"
    fi

    # 4. each MAG's reads aligned to the MAG's own contigs ; 'all' = primary, MAPQ, no hard clip, % / bp aligned
    if [[ -f $W/extract.ok && ! -f $W/align.ok ]]; then
        t=$(date +%s); mkdir -p "$W/aln"; : > "$W/aln/counts.tsv"
        ensure_mag_refs "$S" "$asmf" $(tail -n +2 "$W/byMAG/counts.tsv" | cut -f1)
        while IFS=$'\t' read -r mag mol n; do
            local r; r=$(mag_ref "$mag")
            "$PBMM2" align --preset CCS --sample "$S" -y "$CTGMETH_MIN_IDENTITY" -j "$THREADS" "$r/mag.mmi" "$W/byMAG/$mag.bam" "$W/aln/$mag.raw.bam"
            "$SAMTOOLS" sort -n -@ "$SORT_T" -m "$SORT_MEM" -T "$W/aln/$mag.ns" -o "$W/aln/$mag.n.bam" "$W/aln/$mag.raw.bam"
            "$SAMTOOLS" view -u -F 2308 -q "$CTGMETH_MIN_MAPQ" "$W/aln/$mag.raw.bam" \
                | "$SAMTOOLS" view -u -e "hclen==0 && (qlen-sclen)*100 >= $CTGMETH_MIN_PCT_ALIGNED*qlen && (qlen-sclen) >= $CTGMETH_MIN_ALN_BP" - \
                | "$SAMTOOLS" sort -@ "$SORT_T" -m "$SORT_MEM" -T "$W/aln/$mag.srt" -o "$W/aln/$mag.all.bam" -
            printf '%s\t%s\t%s\t%s\t%s\n' "$mag" "$mol" "$n" "$("$SAMTOOLS" view -c -F 2308 "$W/aln/$mag.raw.bam")" \
                "$("$SAMTOOLS" view -c "$W/aln/$mag.all.bam")" >> "$W/aln/counts.tsv"
            rm -f "$W/aln/$mag.raw.bam"
        done < <(tail -n +2 "$W/byMAG/counts.tsv")
        log "$S: aligned to $(wc -l < "$W/aln/counts.tsv") MAGs: $(awk -F'\t' '{s += $5} END {print s + 0}' "$W/aln/counts.tsv") pseudo-reads in 'all' ($(( $(date +%s) - t )) s)"
        touch "$W/align.ok"
    fi

    # 5. junctions: breakpoints, read sets, windows (lib/ctgmeth.py)
    if [[ -f $W/align.ok && ! -f $W/junctions.ok ]]; then
        t=$(date +%s); rm -rf "$W/jx"
        : > "$W/mags.tsv"
        while IFS=$'\t' read -r mag _; do
            m=$(awk -F'\t' -v g="$mag" '$2 == g {print $3; exit}' "$W/contigs.tsv")
            printf '%s\t%s\t%s\t%s\t%s\n' "$mag" "$m" "$(mag_ref "$mag")/mag.fa" "$W/aln/$mag.n.bam" "$W/byMAG/$mag.bam" >> "$W/mags.tsv"
        done < "$W/aln/counts.tsv"
        "$PYTHON" "$MMT_HOME/lib/ctgmeth.py" junctions --mags "$W/mags.tsv" --settings "$CM_JSON" --out "$W/jx"
        "$SAMTOOLS" sort -@ "$SORT_T" -m "$SORT_MEM" -T "$W/jx/w.srt" -o "$W/jx/win.s.bam" "$W/jx/win.bam"; rm -f "$W/jx/win.bam"
        cat "$W/jx/stats.tsv" >> "$W/run.stats"
        log "$S: junctions: $(awk -F'\t' '$1 == "junction_sets_alt" {a = $2} $1 == "junction_sets_ref" {r = $2} $1 == "junction_sets_called" {c = $2} $1 == "junction_windows" {w = $2} END {printf "%d ALT + %d REF read sets, %d called, %d windows", a, r, c, w}' "$W/jx/stats.tsv") ($(( $(date +%s) - t )) s)"
        touch "$W/junctions.ok"
    fi
    [[ -f $W/junctions.ok ]] && rm -rf "$W/byMAG" "$W"/aln/*.n.bam      # read copies no longer needed (work-space disk)

    # 6. one BAM (every MAG's 'all' reads + the window records) + one reference -> ipdSummary + pileup once
    local nwin=0 nall=0
    if [[ -f $W/junctions.ok ]]; then
        nwin=$("$SAMTOOLS" view -c "$W/jx/win.s.bam"); nall=$(awk -F'\t' '{s += $5} END {print s + 0}' "$W/aln/counts.tsv")
    fi
    if [[ -f $W/junctions.ok && ! -f $W/kinetics.ok ]] && (( nall + nwin > 0 )); then
        t=$(date +%s)
        if [[ ! -f $W/merge.ok ]]; then
            awk -F'\t' -v d="$W/aln" '$5 > 0 {print d "/" $1 ".all.bam"}' "$W/aln/counts.tsv" > "$W/bams.txt"
            (( nwin > 0 )) && echo "$W/jx/win.s.bam" >> "$W/bams.txt"
            "$SAMTOOLS" merge -@ "$THREADS" -c -p -f -o "$W/all.bam" -b "$W/bams.txt"
            "$SAMTOOLS" index -@ "$THREADS" "$W/all.bam"
            "$PBINDEX" "$W/all.bam"
            : > "$W/ref.fa"
            while IFS=$'\t' read -r mag _; do cat "$(mag_ref "$mag")/mag.fa" >> "$W/ref.fa"; done < "$W/aln/counts.tsv"
            cat "$W/jx/windows.fa" >> "$W/ref.fa"
            "$SAMTOOLS" faidx "$W/ref.fa"
            touch "$W/merge.ok"
        fi
        rm -f "$W"/aln/*.all.bam                                           # merged into all.bam
        rm -rf "$W/k"; mkdir -p "$W/k"
        "$SAMTOOLS" idxstats "$W/all.bam" | awk '$3 > 0 {print $1}' > "$W/k/contigs.txt"
        kinetics_table "$W/all.bam" "$W/ref.fa" "$W/k/contigs.txt" "$W/k"
        cat "$W/k/stats.tsv" >> "$W/run.stats"
        log "$S: kinetics table ready ($(( $(date +%s) - t )) s)"
        touch "$W/kinetics.ok"
    fi

    # 7. per MGS present: sites (lifted onto the blocks), junction table, edges, reads
    if [[ ! -f $W/written.ok ]] && (( ${#present[@]} )); then
        t=$(date +%s)
        : > "$W/mgs_list.tsv"
        for m in "${present[@]}"; do printf '%s\t%s\t%s\t%s\t%s\n' "$m" "${KFP[$m]}" "${LFP[$m]}" "${UMAGS[${ASM[$S]}|$m]}" "${DFP[$m]}" >> "$W/mgs_list.tsv"; done
        if [[ -f $W/aln/counts.tsv ]]; then
            while IFS=$'\t' read -r mag _; do cat "$(mag_ref "$mag")/mag.fa.fai"; done < "$W/aln/counts.tsv"
        fi | awk -F'\t' 'NR == FNR {n[$1] = $2; next} ($1 in n) {print $1 "\t" $2 "\t" $3 "\t" n[$1]}' - "$W/contigs.tsv" > "$W/wcontigs.tsv"
        mkdir -p "$W/jx"
        local ex=""
        if [[ -f $W/kinetics.ok ]]; then export_rows "$S"; ex="$W/k/export.parquet"; fi
        "$PYTHON" "$MMT_HOME/lib/ctgmeth.py" write ${ex:+--export "$ex"} --jx "$W/jx" --derep "$DEREP" --out "$D" --sample "$S" \
            --mgs-list "$W/mgs_list.tsv" --settings "$CM_JSON" --contigs "$W/wcontigs.tsv" > "$W/written.tsv"
        [[ $(cut -f1 "$W/written.tsv" | sort | tr '\n' ' ') == "$(printf '%s\n' "${present[@]}" | sort | tr '\n' ' ')" ]] \
            || die "$S: not every MGS was written (see $W/written.tsv)"
        if [[ -f $W/k/rows.tsv ]]; then
            awk -F'\t' 'NR == FNR {n[$1] = $2; next} n[$1] != $2 {bad = 1; print "rows " $1 ": export " n[$1] ", written " $2 > "/dev/stderr"} END {exit bad}' \
                "$W/k/rows.tsv" "$W/written.tsv" || die "$S: rows written differ from the kinetics export"
        fi
        log "$S: written ($(( $(date +%s) - t )) s)"
        touch "$W/written.ok"
    fi
    # counts for the state (molecules: the presence count, absent MGS included)
    while IFS=$'\t' read -r m n; do MOL[$m]=$n; done < <(tail -n +2 "$W/molecules.tsv")
    if [[ -f $W/aln/counts.tsv ]]; then
        while IFS=$'\t' read -r m n; do PSE[$m]=$n; done < <(awk -F'\t' 'NR == FNR {g[$2] = $3; next} {p[g[$1]] += $3} END {for (m in p) print m "\t" p[m]}' "$W/contigs.tsv" "$W/aln/counts.tsv")
        while IFS=$'\t' read -r m n; do ALN[$m]=$n; done < <(awk -F'\t' 'NR == FNR {g[$2] = $3; next} {p[g[$1]] += $5} END {for (m in p) print m "\t" p[m]}' "$W/contigs.tsv" "$W/aln/counts.tsv")
    fi
    if [[ -f $W/written.tsv ]]; then
        while IFS=$'\t' read -r m n a b; do ROWS[$m]=$n; NSETS[$m]=$a; NCALL[$m]=$b; done < "$W/written.tsv"
    fi
}

export_rows() {   # export_rows S : kinetics table -> W/k/export.parquet (MGS, sample, read_set, side, dist, MAG, contig, pos, ...)
    local S=$1 K="$W/k/k.duckdb" chk
    "$DUCKDB" -bail "$K" -c "SET threads=$THREADS; SET memory_limit='$DUCK_MEM'; SET temp_directory='$W/k/duck_tmp';
        DROP TABLE IF EXISTS ctg; DROP TABLE IF EXISTS win; DROP TABLE IF EXISTS ex;
        CREATE TABLE ctg AS SELECT * FROM read_csv('$W/wcontigs.tsv', delim='\t', header=false,
            columns={'contig':'VARCHAR','MAG':'VARCHAR','MGS':'VARCHAR','len':'BIGINT'});
        CREATE TABLE win AS SELECT * FROM read_csv('$W/jx/windows.tsv', delim='\t', header=true,
            types={'wname':'VARCHAR','read_set':'VARCHAR','side':'INTEGER','MGS':'VARCHAR','MAG':'VARCHAR','contig':'VARCHAR',
                   'w_c_start':'BIGINT','w_c_end':'BIGINT','prefix':'BIGINT','bp':'BIGINT','bp_side':'VARCHAR','D':'BIGINT'});
        CREATE TABLE ex AS
        SELECT c.MGS, '$S' AS sample, 'all' AS read_set, NULL::TINYINT AS side, NULL::INTEGER AS dist, c.MAG, k.contig,
               k.pos::INTEGER AS pos, CASE k.strand WHEN 0 THEN '+' ELSE '-' END AS strand, k.base, k.coverage, k.score, k.ipdRatio,
               k.tMean, k.tErr, k.modelPrediction, k.frac, k.fracLow, k.fracUp, k.call, k.identificationQv, k.pile_depth, k.n_ref,
               k.top_alt, k.n_alt, k.ctx_alt_pct, k.ctx_alt_n, k.clip_ctx_n::USMALLINT AS clip_ctx_n
        FROM k JOIN ctg c ON c.contig = k.contig
        WHERE k.base IN ('A', 'C') AND k.coverage >= $CTGMETH_MIN_ROW_COV
        UNION ALL
        SELECT MGS, '$S', read_set, side, dist, MAG, wcontig, p::INTEGER, CASE strand WHEN 0 THEN '+' ELSE '-' END, base, coverage,
               score, ipdRatio, tMean, tErr, modelPrediction, frac, fracLow, fracUp, call, identificationQv, pile_depth, n_ref, top_alt,
               n_alt, ctx_alt_pct, ctx_alt_n, clip_ctx_n
        FROM (SELECT w.MGS, w.read_set, w.side::TINYINT AS side, w.MAG, w.contig AS wcontig, w.w_c_start + k.pos - 1 - w.prefix AS p,
                     (CASE WHEN w.bp_side = 'R' THEN w.bp - (w.w_c_start + k.pos - 1 - w.prefix)
                           ELSE (w.w_c_start + k.pos - 1 - w.prefix) - w.bp END)::INTEGER AS dist, w.D, w.prefix, w.w_c_end,
                     k.pos AS wpos, k.strand, k.base, k.coverage, k.score, k.ipdRatio, k.tMean, k.tErr, k.modelPrediction, k.frac, k.fracLow,
                     k.fracUp, k.call, k.identificationQv, k.pile_depth, k.n_ref, k.top_alt, k.n_alt, k.ctx_alt_pct, k.ctx_alt_n,
                     k.clip_ctx_n::USMALLINT AS clip_ctx_n
              FROM k JOIN win w ON w.wname = k.contig) x
        WHERE base IN ('A', 'C') AND coverage >= $CTGMETH_MIN_ROW_COV AND wpos > prefix AND p <= w_c_end AND dist BETWEEN 0 AND D;
        COPY (SELECT * FROM ex ORDER BY MGS, read_set, side, contig, pos, strand) TO '$W/k/export.parquet' (FORMAT parquet, COMPRESSION zstd);"
    # every kinetics row lies on a MAG contig or a window ; 'all' rows = every A/C row of the MAG contigs
    chk=$("$DUCKDB" -readonly "$K" -noheader -list -c "SELECT
        (SELECT count(*) FROM k WHERE contig NOT IN (SELECT contig FROM ctg) AND contig NOT IN (SELECT wname FROM win)) || ' ' ||
        (SELECT count(*) FROM k WHERE contig IN (SELECT contig FROM ctg) AND base IN ('A','C') AND coverage >= $CTGMETH_MIN_ROW_COV) || ' ' ||
        (SELECT count(*) FROM ex WHERE read_set = 'all')")
    read -r a b c <<< "$chk"
    [[ $a == 0 && $b == "$c" ]] || die "$S: kinetics rows on unknown contigs $a ; A/C rows on MAG contigs $b, exported $c"
    "$DUCKDB" -readonly "$K" -noheader -csv -separator $'\t' -c "SELECT MGS, count(*) FROM ex GROUP BY MGS ORDER BY MGS" > "$W/k/rows.tsv"
}

write_state() {   # write_state S [T0] : stats/<S>.tsv = old lines minus TODO_KIN / TODO_LIFT (replaced) and STALE (dropped
                  # with their files) ; one SAMPLE line in OUT/mmt.log
    local S=$1 ts=${2:-} m st now t kind line="" x
    local -a out=()
    now=$(date '+%F %T')
    for m in "${STALE[@]}"; do
        for kind in "" junctions/ edges/ reads/; do rm -f "$D/$m/$kind$S.parquet"; done
        log "$S: $m dropped (no longer usable, or no MAG of it in this sample's assembly)"
    done
    t=$(tmpname "$STATE/$S.tsv")
    {
        printf 'MGS\tstatus\tkin_fp\tlift_fp\tmolecules\tpseudo_reads\taligned\trows\tread_sets\tcalled\tupdated\tmmt\tipdsummary\n'
        if [[ -f $STATE/$S.tsv ]]; then
            awk -F'\t' -v OFS='\t' -v drop=" ${TODO_KIN[*]} ${STALE[*]} " -v lift=" ${TODO_LIFT[*]} " -v now="$now" \
                -v lfp="$(for m in "${TODO_LIFT[@]}"; do printf '%s=%s ' "$m" "${LFP[$m]}"; done)" '
                BEGIN {n = split(lfp, x, " "); for (i = 1; i <= n; i++) {split(x[i], y, "="); L[y[1]] = y[2]}}
                FNR > 1 && index(drop, " " $1 " ") == 0 { if (index(lift, " " $1 " ")) {$4 = L[$1]; $11 = now}; print }' "$STATE/$S.tsv"
        fi
        for m in "${TODO_KIN[@]}"; do
            if (( ${ROWS[$m]:-0} > 0 )); then st=done
            elif (( ${MOL[$m]:-0} < CTGMETH_PRESENT_MOL || ${MOL[$m]:-0} == 0 )); then st=absent
                for kind in "" junctions/ edges/ reads/; do rm -f "$D/$m/$kind$S.parquet"; done   # an older result no longer holds
            elif (( ${PSE[$m]:-0} == 0 )); then st=no_kinetics
            elif (( ${ALN[$m]:-0} == 0 )); then st=no_aligned
            else st=no_rows; fi
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$m" "$st" "${KFP[$m]}" "${LFP[$m]}" "${MOL[$m]:-0}" \
                "${PSE[$m]:-0}" "${ALN[$m]:-0}" "${ROWS[$m]:-0}" "${NSETS[$m]:-0}" "${NCALL[$m]:-0}" "$now" "$MMT_VERSION" "${IPD_VERSION:-}"
            if [[ $st == done ]]; then out+=("$m done (${ROWS[$m]} rows, ${NSETS[$m]:-0} junction read sets, ${MOL[$m]:-0} molecules)")
            elif [[ $st == absent ]]; then out+=("$m absent (${MOL[$m]:-0} molecules < $CTGMETH_PRESENT_MOL)")
            else out+=("$m $st (${MOL[$m]:-0} molecules)"); fi
        done
    } > "$t"
    publish "$t" "$STATE/$S.tsv"
    if [[ -f ${W:-}/run.stats ]]; then
        t=$(tmpname "$STATE/$S.run.tsv")
        { printf 'sample\t%s\nupdated\t%s\nipdsummary_args\t%s\n' "$S" "$now" "${MMT_IPD_ARGS:-}"
          cat "$IDX/reads/_stats/$S.tsv" "$W/run.stats" 2> /dev/null || true; } > "$t"
        publish "$t" "$STATE/$S.run.tsv"
    fi
    log "$S: state: $(awk -F'\t' 'FNR > 1 {n[$2]++} END {for (s in n) printf "%s %d  ", s, n[s]}' "$STATE/$S.tsv")${STALE[*]:+ dropped: ${STALE[*]}}"
    for m in "${TODO_LIFT[@]}"; do out+=("$m re-lifted onto the current derep-ctg result"); done
    for m in "${STALE[@]}"; do out+=("$m dropped (no MAG of it in this sample any more, or no longer usable)"); done
    for x in "${out[@]}"; do line+="${line:+ ; }$x"; done
    mlog SAMPLE "$S: $line${ts:+ ; $(( ($(date +%s) - ts + 59) / 60 )) min}"
}
