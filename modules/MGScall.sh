# shellcheck shell=bash
# MMT MGScall - methylation of each MGS per sample: the species' reads pooled on its representative MAG (mgs-rep), each of
# its MAGs in the sample apart on the representative (mag-rep) or on the MAG itself (mag-mag). Sourced by the mmt wrapper.
# Method, outputs and columns: docs/MGScall.md.
# shellcheck source=../lib/kinetics.sh
source "$MMT_HOME/lib/kinetics.sh"
# shellcheck source=../lib/reads.sh
source "$MMT_HOME/lib/reads.sh"

MGSCALL_METHOD=2   # bump when a code change alters MGScall results -> the next run redoes every result
                   # (2: MAG column, G/T rows with a kinetic signal, mag-rep = every MAG apart)

module_help() {
    cat << 'EOF'
mmt MGScall - methylation of each MGS per sample (ipdSummary m6A + m4C)

Usage: mmt MGScall (--mgs-rep | --mag-rep | --mag-mag) -o OUT [options]

Mode (one per run; each writes its own folder OUT/<mode>/)
  --mgs-rep           per sample, the reads of EVERY MAG of the MGS in it, pooled, on the representative
  --mag-rep           per sample, the reads of EACH MAG of the MGS in it, apart, on the representative (strains apart)
  --mag-mag           per sample, the reads of EACH MAG of the MGS in it, apart, on that MAG's own contigs
All modes: every sample whose assembly has a MAG of the MGS. Which reads belong to which MAG comes from MATAFILER4's
own alignment (the sample's CRAM): no re-assignment. Where a sample has one MAG of the MGS, mgs-rep and mag-rep give
the same table: the second of them run links (hard link) the first one's file instead of computing it again.

Inputs (defaults from config/site.conf)
  -i, --mf4-run DIR   MATAFILER4 run folder (default: the map's #OutPath/#RunID)
  --genecat DIR       GeneCat folder (with Bin_<binner>/MAGvsGC.txt.gz)
  --map FILE          MATAFILER4 map (raw kinetics BAMs: #DirPath + SmplPrefix + .bam)
Output, logs and work space
  -o, --out DIR       results folder - reuse it: finished results are kept, new ones added, changed inputs redone
  -l, --log DIR       keep the run logs: DIR/<mode>/<sample>.log (one sample per run) or DIR/<mode>/<date>_<job>.log,
                      ipdSummary's key lines included. Without it no log is kept (a failure prints the log's last lines)
  --tmp DIR           work space: work folders (removed when a sample ends), the MF4 index and reference indexes
                      (caches). On a cluster use a shared scratch folder, so jobs share the caches
Selection
  --select IDS        only these MGS: "MGS.12,MGS.40" or a file (one ID per line); default all usable MGS
  --samples IDS       only these samples (same format); default every usable sample with work
                      HPC: one sample per array task is the efficient unit (each raw BAM is read once per run)
Kinetics
  --ipd-args "ARGS"   extra ipdSummary arguments, as they are (e.g. "--ipdModel /path/model.npz.gz" to force a kinetics
                      model) ; default MMT_IPD_ARGS. Without: ipdSummary picks the model matching each BAM's chemistry. The
                      model used is recorded per sample (stats/<S>.run.tsv) ; changing the arguments redoes the results
Resources (default: the SLURM allocation)
  -t, --threads N     CPUs to use
  --mem GB            memory to plan for (DuckDB, sorting)
Other
  --dry-run           show what would run (per sample: MGS to do / up to date) and stop
  --list-samples      print the samples that have work (one per line, nothing else) and stop - for SLURM arrays
  --force             redo the selected results even if up to date (computed, never linked)
  --debug             keep each sample's work folder and run log (a rerun resumes from the work folder)
  --config FILE       extra settings file (after config/site.conf)

Output (in OUT)
  mmt_<date>.log                         the folder's log: every run, its plan per MGS (incl. where GeneCat's abundance
                                         table has the MGS but no MAG was binned), each sample x MGS outcome, failures
  mgs-rep/<MGS>/<S>.parquet              one row per A/C (and flagged G/T) of the representative x strand with >= 3
                                         kinetic reads ; columns MGS, MAG, sample, contig, pos (1-based), strand, ...
  mag-rep/<MGS>/<S>/<MAG>.parquet        the same, from one MAG's reads
  mag-mag/<MGS>/<S>/<MAG>.parquet        one MAG's reads on its own contigs (contig = the MAG's contig)
  <mode>/stats/<S>.tsv                   per MGS (mag-rep, mag-mag: per MAG): status, abundance, molecules, rows, fingerprint
  <mode>/stats/<S>.run.tsv               the sample's last run: raw BAM, chemistry, kinetics model, read accounting
Read all of a mode with DuckDB: SELECT * FROM read_parquet('OUT/mgs-rep/*/*.parquet') ;
mag-rep, mag-mag: read_parquet('OUT/mag-rep/*/*/*.parquet') (or OUT/mag-mag/...)

Examples
  one job, all samples one after another:
    mmt MGScall --mgs-rep -o OUT --select MGS.12,MGS.40
  SLURM array, one sample per task (list once, before submitting):
    mmt MGScall --mgs-rep -o OUT --select MGS.12,MGS.40 --list-samples > todo.txt
    sbatch --array=1-$(wc -l < todo.txt)%20 job.sh
    # job.sh:  S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" todo.txt)
    #          mmt MGScall --mgs-rep -o OUT --select MGS.12,MGS.40 --samples "$S" --tmp /scratch/me/mmt_tmp -l LOGS
  what is done per MGS (printed and added to the mmt log):
    mmt status -o OUT --select MGS.12,MGS.40
EOF
}

module_main() {
    MODE=""
    local SELECT="" SAMPLES="" FORCE=0 DRY=0 LIST=0 t0 S msg n_ok=0 n_skip=0 err desc lname
    local -a failed=() todo_samples=()
    MMT_DEBUG=${MMT_DEBUG:-0}
    while (( $# )); do
        if common_opt "$@"; then shift "$COMMON_SHIFT"; continue; fi
        case $1 in
            --mgs-rep|--mag-rep|--mag-mag)
                [[ -z $MODE || $MODE == "${1#--}" ]] || usage_die "one mode per run: --mgs-rep, --mag-rep or --mag-mag"
                MODE=${1#--} ;;
            --select)  [[ $# -ge 2 && -n $2 ]] || usage_die "--select needs MGS IDs or a file"; SELECT=$2; shift ;;
            --samples) [[ $# -ge 2 && -n $2 ]] || usage_die "--samples needs sample IDs or a file"; SAMPLES=$2; shift ;;
            --force)    FORCE=1 ;;
            --dry-run)  DRY=1 ;;
            --list-samples) LIST=1 ;;
            -h|--help)  module_help; exit 0 ;;
            *) usage_die "unknown option '$1'" ;;
        esac
        shift
    done
    [[ -n $MODE ]] || usage_die "choose a mode: --mgs-rep, --mag-rep or --mag-mag"
    case $MODE in mgs-rep) OTHER=mag-rep ;; mag-rep) OTHER=mgs-rep ;; *) OTHER="" ;; esac   # the mode whose tables can be linked
    (( LIST )) && MMT_QUIET=1
    t0=$(date +%s)
    resolve_out
    job_tag
    if [[ -n $SAMPLES && ! -f $SAMPLES && $SAMPLES != *[,\ ]* ]]; then lname=$SAMPLES      # one sample: its own log
    else lname="$(date +%Y%m%d-%H%M%S)_$JOBTAG$( (( LIST )) && echo _list || true)"; fi
    start_log "$MODE" "$lname"
    say "MMT $(mmt_version) MGScall --$MODE -> $MMT_OUT/$MODE${MMT_LOG_DIR:+ (log: $MMT_LOG)}"
    log "command: mmt MGScall --$MODE ${SELECT:+--select $SELECT }${SAMPLES:+--samples $SAMPLES }(force=$FORCE dry=$DRY list=$LIST debug=$MMT_DEBUG) on $(hostname), job $JOBTAG"
    desc="${SELECT:+--select $SELECT }${SAMPLES:+--samples $SAMPLES }"
    (( FORCE )) && desc+="--force "
    (( LIST )) && desc+="(listing the samples with work) "
    (( MMT_DEBUG )) && desc+="--debug "
    (( DRY )) || mlog_start "$MODE" "${desc% }"
    # a run that stops early (an error, scancel, the time limit) says so in the mmt log ; the run log is removed (no -l)
    trap 'run_exit $?' EXIT
    trap 'exit 143' TERM INT
    resolve_inputs
    resolve_resources
    setup_tools
    DUCK_MEM="$(( MEM_GB * 6 / 10 > 1 ? MEM_GB * 6 / 10 : 1 ))GB"
    SORT_T=$(( THREADS < 4 ? THREADS : 4 ))
    SORT_MEM="$(( MEM_GB * 1024 * 3 / 10 / SORT_T > 256 ? MEM_GB * 1024 * 3 / 10 / SORT_T : 256 ))M"
    MIN_MAPQ=$MGSCALL_MIN_MAPQ CTX=$MGSCALL_CTX CLIP_MIN=$MGSCALL_CLIP_MIN
    SETTINGS="keepdup=$MGSCALL_KEEP_DUP id=$MGSCALL_MIN_IDENTITY mapq=$MGSCALL_MIN_MAPQ pct=$MGSCALL_MIN_PCT_ALIGNED bp=$MGSCALL_MIN_ALN_BP rowcov=$MGSCALL_MIN_ROW_COV ctx=$MGSCALL_CTX clip=$MGSCALL_CLIP_MIN ipd=${MMT_IPD_ARGS:-}"
    log "inputs: MF4 run $MMT_MF4_RUN | GeneCat $MMT_GENECAT (Bin_$MMT_BINNER) | map $MMT_MAP"
    log "resources: $THREADS threads, $MEM_GB GB (DuckDB $DUCK_MEM, sort $SORT_T x $SORT_MEM) | work space $WROOT"
    log "settings: $SETTINGS | method $MGSCALL_METHOD"
    tool_versions | sed 's/^/tool: /'

    ensure_index
    D="$MMT_OUT/$MODE"; STATE="$D/stats"
    mkdir -p "$STATE"
    pick_mgs
    pick_samples
    log "selection: ${#MGS_SEL[@]} MGS, ${#SAMPLE_SEL[@]} samples"

    # plan: which MGS of which sample need (re)computing
    local n_units=0 n_cur=0 m
    declare -A P_TODO=() P_CUR=()
    for S in "${SAMPLE_SEL[@]}"; do
        plan_sample "$S"
        n_units=$(( n_units + ${#TODO[@]} )); n_cur=$(( n_cur + N_CURRENT ))
        for m in "${TODO[@]}"; do P_TODO[$m]+="$S "; done
        for m in "${CURRENT[@]}"; do P_CUR[$m]+="$S "; done
        if (( ${#TODO[@]} + ${#STALE[@]} )); then todo_samples+=("$S"); else n_skip=$(( n_skip + 1 )); fi
        msg="plan $S: ${#TODO[@]} to do, $N_CURRENT up to date${STALE[*]:+, to drop: ${STALE[*]}}"
        if (( DRY )); then say "$msg"; else log "$msg"; fi
    done
    say "plan: ${#todo_samples[@]} sample(s) with work ($n_units MGS x sample results to make), $n_skip sample(s) up to date, $n_cur result(s) already current"
    if [[ -z $SAMPLES ]]; then plan_log; fi       # per MGS: where it has a MAG / is present without one, to do / up to date
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
            mlog FAIL "$S: ${err:-exit $RC} - $(log_ref)$( (( MMT_DEBUG )) && echo " ; work folder kept: $WROOT/$MODE/$S")"
        fi
    done
    mlog END "$n_ok sample(s) processed, ${#failed[@]} failed$( (( ${#failed[@]} )) && echo " ($(id_list "${failed[@]}"))"), $n_skip already up to date ; $(( ($(date +%s) - t0) / 60 )) min"
    RUN_ENDED=1
    if (( ${#failed[@]} )); then
        say "finished with errors after $(( ($(date +%s) - t0) / 60 )) min: ${#failed[@]} sample(s) failed: ${failed[*]} - rerun the same command to resume"
        exit 1
    fi
    say "done in $(( ($(date +%s) - t0) / 60 )) min: $n_ok sample(s) processed, $n_skip already up to date"
}

# ============================================================================ selection + plan
pick_mgs() {   # MGS_SEL (usable, selected) ; REP_FASTA[MGS] ; REP_MAG[MGS] ; REP_STAMP[MGS] ; OK_MGS[MGS] ; IS_SEL[MGS]
    declare -gA REP_FASTA=() REP_MAG=() REP_STAMP=() OK_MGS=() IS_SEL=()
    declare -ga MGS_SEL=()
    local m st r fa want=() bad=() unk=()
    declare -A seen=()
    while IFS=$'\t' read -r m st r _ _ fa _; do
        seen[$m]=$st
        [[ $st == ok ]] || continue
        OK_MGS[$m]=1; REP_FASTA[$m]=$fa; REP_MAG[$m]=$r
    done < <(tail -n +2 "$IDX/mgs.tsv")
    if [[ -n $SELECT ]]; then
        mapfile -t want < <(read_ids "$SELECT")
        (( ${#want[@]} )) || usage_die "--select: no MGS IDs in '$SELECT'"
        for m in "${want[@]}"; do
            if [[ -z ${seen[$m]:-} ]]; then unk+=("$m")
            elif [[ -z ${OK_MGS[$m]:-} ]]; then bad+=("$m(${seen[$m]})")
            else MGS_SEL+=("$m"); fi
        done
        (( ${#unk[@]} == 0 )) || die "--select: not in MAGvsGC: ${unk[*]}"
        (( ${#bad[@]} == 0 )) || { say "WARNING: skipped, representative not usable: ${bad[*]} (see $IDX/report.txt)"
                                   mlog WARN "skipped, representative not usable: ${bad[*]} (see $IDX/report.txt)"; }
    else
        mapfile -t MGS_SEL < <(awk -F'\t' 'NR > 1 && $2 == "ok" {print $1}' "$IDX/mgs.tsv")
    fi
    (( ${#MGS_SEL[@]} )) || die "no usable MGS selected"
    for m in "${MGS_SEL[@]}"; do REP_STAMP[$m]=$(stamp "${REP_FASTA[$m]}"); IS_SEL[$m]=1; done
}

pick_samples() {   # SAMPLE_SEL ; OKS (every usable sample) ; RAW[S] CRAM[S] ASM[S] ; UNITS[S] = the selected MGS with a MAG
                   # in that sample's assembly ; UMAGS[assembly|MGS] = those MAGs (comma-separated) ; UMD5[assembly|MGS] =
                   # checksum of their contigs ; ASM_MD5 ; ABUND[S|MGS] = GeneCat abundance (> 0 only)
    declare -gA RAW=() CRAM=() ASM=() UNITS=() UMAGS=() UMD5=() ASM_MD5=() ABUND=()
    declare -ga SAMPLE_SEL=() OKS=()
    local s st raw cram asm m mags mall v want=() bad=() unk=()
    declare -A seen=() asm_units=()
    while IFS=$'\t' read -r asm m mags _ mall _; do
        [[ -n ${IS_SEL[$m]:-} ]] || continue
        asm_units[$asm]+="$m "; UMAGS[$asm|$m]=$mags; UMD5[$asm|$m]=$mall
    done < <(tail -n +2 "$IDX/units.tsv")
    while IFS=$'\t' read -r asm mall; do ASM_MD5[$asm]=$mall; done < <(tail -n +2 "$IDX/assemblies.tsv")
    while IFS=$'\t' read -r s m v; do [[ -n ${IS_SEL[$m]:-} ]] && ABUND[$s|$m]=$v; done < <(tail -n +2 "$IDX/abundance.tsv")
    while IFS=$'\t' read -r s st raw cram asm; do
        seen[$s]=$st
        [[ $st == ok ]] || continue
        OKS+=("$s"); RAW[$s]=$raw; CRAM[$s]=$cram; ASM[$s]=$asm; UNITS[$s]=${asm_units[$asm]:-}
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
        (( ${#bad[@]} == 0 )) || { say "WARNING: skipped, sample not usable: ${bad[*]} (see $IDX/report.txt)"
                                   mlog WARN "skipped, sample not usable: ${bad[*]} (see $IDX/report.txt)"; }
    else
        while IFS=$'\t' read -r s st _; do   # samples with work, or with older results that may need dropping
            [[ $st == ok ]] && [[ -n ${UNITS[$s]:-} || -f $STATE/$s.tsv ]] && SAMPLE_SEL+=("$s")
        done < <(tail -n +2 "$IDX/samples.tsv")
    fi
    return 0
}

declare -A SFP=()
unit_fp() {   # unit_fp MODE S MGS -> UFP = what a result depends on (method, mode, settings, the sample's CRAM + raw BAM,
              # the representative's FASTA (mag-mag: the assembly FASTA), the contigs of the sample's MAGs of the MGS) ;
              # the sample part is computed once
    local k="$1|$2" a=""
    if [[ -z ${SFP[$k]:-} ]]; then
        [[ $1 != mag-mag ]] || { a=$(asm_fasta "${ASM[$2]}" || echo none); a="|$a:$(stamp "$a")"; }
        SFP[$k]=$(md5str "$MGSCALL_METHOD|$1|$SETTINGS|${CRAM[$2]}:$(stamp "${CRAM[$2]}")|${RAW[$2]}:$(stamp "${RAW[$2]}")$a")
    fi
    UFP="${SFP[$k]}:${REP_STAMP[$3]}:${UMD5[${ASM[$2]}|$3]}"
}
res_file() {   # res_file MODE MGS S MAG : the result file of one unit
    if [[ $1 == mgs-rep ]]; then echo "$MMT_OUT/mgs-rep/$2/$3.parquet"; else echo "$MMT_OUT/$1/$2/$3/$4.parquet"; fi
}

plan_sample() {   # plan_sample S -> TODO (MGS to (re)do), CURRENT (up to date), STALE (old results to drop), N_CURRENT ;
                  # FP[MGS] = new fingerprint. A result is current when every state row of its MGS has the new fingerprint
                  # and every done / linked row still has its file.
    local S=$1 m g st fp
    declare -gA FP=()
    declare -ga TODO=() CURRENT=() STALE=()
    declare -A old_fp=() bad=()
    N_CURRENT=0
    if [[ -f $STATE/$S.tsv ]]; then
        while IFS=$'\t' read -r m g st fp _; do
            old_fp[$m]=$fp
            [[ $st == done || $st == linked ]] && [[ ! -s $(res_file "$MODE" "$m" "$S" "$g") ]] && bad[$m]=1
        done < <(tail -n +2 "$STATE/$S.tsv")
    fi
    for m in ${UNITS[$S]:-}; do
        unit_fp "$MODE" "$S" "$m"; FP[$m]=$UFP
        if (( ! FORCE )) && [[ ${old_fp[$m]:-} == "${FP[$m]}" && -z ${bad[$m]:-} ]]; then
            N_CURRENT=$(( N_CURRENT + 1 )); CURRENT+=("$m")
        else
            TODO+=("$m")
        fi
    done
    #   an old result is dropped when its MGS is no longer usable, or is selected but has no MAG in this assembly any more
    for m in "${!old_fp[@]}"; do
        [[ " ${UNITS[$S]:-} " == *" $m "* ]] && continue
        if [[ -z ${OK_MGS[$m]:-} || -n ${IS_SEL[$m]:-} ]]; then STALE+=("$m"); fi
    done
    return 0
}

# ============================================================================ one sample
process_sample() {
    local S=$1
    local W LOCK="$STATE/$S.lock" t m g k r mol pseudo n f mof keyf="MGS" ts L
    ts=$(date +%s)
    lock_take "$LOCK" || die "$S: another job is working on this sample ($(cat "$LOCK/owner" 2> /dev/null)) - remove $LOCK if not"
    W="$WROOT/$MODE/$S"                       # this sample's work folder: ours while we hold the lock
    # whatever ends this sample (done, error, scancel / time limit): work folder removed (kept with --debug), lock released
    trap 'rc=$?; work_done "'"$W"'" $(( rc == 0 )); lock_release "'"$LOCK"'"' EXIT
    trap 'exit 143' TERM INT
    declare -gA MOL=() PSE=() ALN=() ROWS=() LINKED=()
    plan_sample "$S"                          # again, under the lock (another job may have finished some)
    if (( ${#TODO[@]} == 0 )); then
        if (( ${#STALE[@]} )); then write_state "$S" "" "$ts"; else log "== $S: up to date (another job finished it)"; fi
        return 0
    fi
    log "== $S: ${#TODO[@]} MGS to do ($MODE) | raw ${RAW[$S]} | CRAM ${CRAM[$S]}"
    if [[ $MODE == mag-mag ]]; then           # the MAG references come from the sample's assembly (cache key: MAG + assembly)
        MAGREF_ASM=$(asm_fasta "${ASM[$S]}") || die "$S: no scaffolds.fasta.filt(.gz) in ${ASM[$S]}"
    fi

    # 0. one MAG of the MGS in this sample: mgs-rep and mag-rep give the same table -> link the other one's current file
    local -a COMP=()
    for m in "${TODO[@]}"; do
        if link_other "$S" "$m"; then LINKED[$m]=1; else COMP+=("$m"); fi
    done
    (( ${#LINKED[@]} == 0 )) || log "$S: linked from $OTHER (one MAG, same reads): ${!LINKED[*]}"
    if (( ${#COMP[@]} == 0 )); then write_state "$S" "" "$ts"; log "== $S: done"; return 0; fi

    local wfp
    wfp=$(md5str "${COMP[*]}|$(for m in "${COMP[@]}"; do printf '%s ' "${FP[$m]}"; done)")
    if [[ $(cat "$W/work.fp" 2> /dev/null) != "$wfp" ]] || (( FORCE )); then rm -rf "$W"; fi   # other work than last time, or --force: start clean
    [[ -d $W ]] && log "$S: resuming from the work folder of an earlier run: $W"
    mkdir -p "$W"; echo "$wfp" > "$W/work.fp"
    [[ -f $W/run.stats ]] || : > "$W/run.stats"

    # 1. read -> MAG -> MGS for this sample (shared by both modes and ctg-meth; a cache next to the index)
    t=$(date +%s); ensure_reads "$S" "$W"; log "$S: reads index ok ($(( $(date +%s) - t )) s)"

    # 2. the reads to use, per unit: the MGS (mgs-rep) or each of its MAGs in this sample (mag-rep, mag-mag). keys.tsv =
    #    unit, MGS, layer (kinetics run: one unit per MGS in each ; mag-*: the n-th MAG of every MGS in run n)
    [[ $MODE == mgs-rep ]] || keyf="MAG"
    if [[ ! -f $W/select.ok ]]; then
        for m in "${COMP[@]}"; do
            if [[ $MODE != mgs-rep ]]; then
                L=0; for g in $(tr ',' '\n' <<< "${UMAGS[${ASM[$S]}|$m]}" | sort); do L=$(( L + 1 )); printf '%s\t%s\t%s\n' "$g" "$m" "$L"; done
            else printf '%s\t%s\t1\n' "$m" "$m"; fi
        done > "$W/keys.tsv"
        #   mag-mag: contig, MAG, MGS of these MAGs (their references)
        [[ $MODE != mag-mag ]] || gzip -dc "$IDX/contigs.tsv.gz" | awk -F'\t' -v a="${ASM[$S]}" \
            'NR == FNR {u[$1]; next} $1 == a && ($3 in u) {print $2 "\t" $3 "\t" $4}' "$W/keys.tsv" - > "$W/contigs.tsv"
        local dupf=""; (( MGSCALL_KEEP_DUP )) || dupf="AND NOT dup"
        "$DUCKDB" -bail -c "SET threads=$THREADS; SET memory_limit='$DUCK_MEM';
            CREATE TEMP TABLE k AS SELECT * FROM read_csv('$W/keys.tsv', header=false, delim='\t', columns={'ukey':'VARCHAR','m':'VARCHAR','layer':'INTEGER'});
            CREATE TEMP TABLE s AS SELECT read, $keyf AS ukey FROM read_parquet('$READS_PQ')
                WHERE MGS IN (SELECT m FROM k) AND $keyf IN (SELECT ukey FROM k) $dupf;
            COPY (SELECT k.ukey AS unit, count(s.read) AS molecules FROM k LEFT JOIN s ON s.ukey = k.ukey GROUP BY k.ukey ORDER BY k.ukey)
                TO '$W/molecules.tsv' (HEADER, DELIMITER '\t');
            COPY (SELECT read, ukey AS unit FROM s ORDER BY read) TO '$W/sel.tsv' (HEADER, DELIMITER '\t')"
        tail -n +2 "$W/sel.tsv" | cut -f1 > "$W/names.txt"
        printf 'selected_molecules\t%s\n' "$(wc -l < "$W/names.txt")" >> "$W/run.stats"
        log "$S: molecules per ${keyf} (MATAFILER4's CRAM$( (( MGSCALL_KEEP_DUP )) || echo ", duplicates left out")): $(awk -F'\t' 'NR > 1 {printf "%s %d  ", $1, $2}' "$W/molecules.tsv")"
        [[ -s $W/names.txt ]] || log "$S: no reads on these MAGs - raw BAM not read"
        touch "$W/select.ok"
    fi

    # 3. those reads from the raw kinetics BAM (one pass) -> one pseudo-read per strand -> one BAM per unit
    if [[ -s $W/names.txt && ! -f $W/extract.ok ]]; then
        t=$(date +%s)
        "$SAMTOOLS" view -@ "$THREADS" -b -N "$W/names.txt" -o "$W/sel.bam" "${RAW[$S]}"
        local want found
        want=$(wc -l < "$W/names.txt"); found=$("$SAMTOOLS" view -@ "$THREADS" -c "$W/sel.bam")
        printf 'raw_reads_found\t%s\n' "$found" >> "$W/run.stats"
        (( found * 100 >= want * 99 )) \
            || die "$S: only $found of $want selected reads are in ${RAW[$S]} - read names differ between CRAM and raw BAM?"
        "$BYSTRANDIFY" "$W/sel.bam" "$W/bys.bam"
        rm -f "$W/sel.bam"
        rm -rf "$W/byKey"
        mof=$(( $(ulimit -n) - 100 )); (( mof > 2000 )) && mof=2000; (( mof < 50 )) && mof=50   # BAMs open at once
        "$PYTHON" "$MMT_HOME/lib/mf4.py" split --bam "$W/bys.bam" --reads "$W/sel.tsv" --outdir "$W/byKey" \
            --threads "$SORT_T" --max-open "$mof"
        rm -f "$W/bys.bam"
        log "$S: $found molecules extracted and split into $(( $(wc -l < "$W/byKey/counts.tsv") - 1 )) ${keyf} ($(( $(date +%s) - t )) s)"
        touch "$W/extract.ok"
    fi

    # 4. each unit's reads aligned to its MGS's representative only (mag-mag: to its MAG only ; MATAFILER4 already placed them)
    #    kept: primary, MAPQ >= MIN_MAPQ, identity >= MIN_IDENTITY (pbmm2 -y), no hard clip, >= MIN_PCT_ALIGNED % of the
    #    read and >= MIN_ALN_BP aligned
    if [[ -f $W/extract.ok && ! -f $W/align.ok ]]; then
        t=$(date +%s); mkdir -p "$W/aln"; : > "$W/aln/counts.tsv"
        [[ $MODE != mag-mag ]] || ensure_mag_refs "$S" "$MAGREF_ASM" $(tail -n +2 "$W/byKey/counts.tsv" | cut -f1)
        while IFS=$'\t' read -r k mol pseudo; do
            m=$(awk -F'\t' -v k="$k" '$1 == k {print $2; exit}' "$W/keys.tsv")
            if [[ $MODE == mag-mag ]]; then r="$(mag_ref "$k")/mag.mmi"; else ensure_ref "$m"; r="$(ref_dir "$m")/rep.mmi"; fi
            "$PBMM2" align --preset CCS --sample "$S" -y "$MGSCALL_MIN_IDENTITY" -j "$THREADS" \
                "$r" "$W/byKey/$k.bam" "$W/aln/$k.raw.bam"
            "$SAMTOOLS" view -u -F 2308 -q "$MGSCALL_MIN_MAPQ" "$W/aln/$k.raw.bam" \
                | "$SAMTOOLS" view -u -e "hclen==0 && (qlen-sclen)*100 >= $MGSCALL_MIN_PCT_ALIGNED*qlen && (qlen-sclen) >= $MGSCALL_MIN_ALN_BP" - \
                | "$SAMTOOLS" sort -@ "$SORT_T" -m "$SORT_MEM" -T "$W/aln/$k.srt" -o "$W/aln/$k.bam" -
            n=$("$SAMTOOLS" view -c "$W/aln/$k.bam")
            rm -f "$W/aln/$k.raw.bam"
            (( n > 0 )) || rm -f "$W/aln/$k.bam"
            printf '%s\t%s\t%s\t%s\n' "$k" "$mol" "$pseudo" "$n" >> "$W/aln/counts.tsv"
        done < <(tail -n +2 "$W/byKey/counts.tsv")
        rm -rf "$W/byKey"
        log "$S: aligned $(awk -F'\t' '{s += $4} END {print s + 0}' "$W/aln/counts.tsv") pseudo-reads of $(awk -F'\t' '$4 > 0' "$W/aln/counts.tsv" | wc -l) ${keyf} to their references ($(( $(date +%s) - t )) s)"
        touch "$W/align.ok"
    fi

    # 5. per layer one BAM + one reference -> ipdSummary + pileup once (lib/kinetics.sh) -> one table per unit
    if [[ -f $W/align.ok ]]; then
        for L in $(cut -f3 "$W/keys.tsv" | sort -nu); do
            awk -F'\t' -v L="$L" 'NR == FNR {if ($3 == L) u[$1]; next} ($1 in u) && $4 > 0 {print $1}' "$W/keys.tsv" "$W/aln/counts.tsv" > "$W/layer$L.txt"
            [[ -s $W/layer$L.txt ]] || continue
            if [[ ! -f $W/k$L.ok ]]; then
                t=$(date +%s)
                if [[ ! -f $W/merge$L.ok ]]; then
                    sed "s|^|$W/aln/|; s|\$|.bam|" "$W/layer$L.txt" > "$W/bams$L.txt"
                    "$SAMTOOLS" merge -@ "$THREADS" -c -p -f -o "$W/all$L.bam" -b "$W/bams$L.txt"
                    "$SAMTOOLS" index -@ "$THREADS" "$W/all$L.bam"
                    "$PBINDEX" "$W/all$L.bam"
                    : > "$W/ref$L.fa"
                    [[ $MODE != mag-mag ]] || ensure_mag_refs "$S" "$MAGREF_ASM" $(cat "$W/layer$L.txt")
                    while read -r k; do
                        if [[ $MODE == mag-mag ]]; then cat "$(mag_ref "$k")/mag.fa" >> "$W/ref$L.fa"; continue; fi
                        m=$(awk -F'\t' -v k="$k" '$1 == k {print $2; exit}' "$W/keys.tsv"); ensure_ref "$m"; cat "$(ref_dir "$m")/rep.fa" >> "$W/ref$L.fa"
                    done < "$W/layer$L.txt"
                    "$SAMTOOLS" faidx "$W/ref$L.fa"
                    xargs rm -f < "$W/bams$L.txt"     # merged: the per-unit BAMs are no longer needed (disk)
                    touch "$W/merge$L.ok"
                fi
                rm -rf "$W/k$L"; mkdir -p "$W/k$L"
                "$SAMTOOLS" idxstats "$W/all$L.bam" | awk '$3 > 0 {print $1}' > "$W/k$L/contigs.txt"
                kinetics_table "$W/all$L.bam" "$W/ref$L.fa" "$W/k$L/contigs.txt" "$W/k$L"
                log "$S: kinetics table ready$( [[ $MODE != mgs-rep ]] && echo " (MAG layer $L: $(wc -l < "$W/layer$L.txt") MAG(s))") ($(( $(date +%s) - t )) s)"
                touch "$W/k$L.ok"
            fi
            write_layer "$S" "$L"
        done
    fi

    # 6. state of every MGS done in this run (work folder: removed by the EXIT trap)
    write_state "$S" "$W" "$ts"
    log "== $S: done"
}

link_other() {   # link_other S MGS : the sample has one MAG of the MGS and the other mode holds its current table -> hard link
                 # it here (a copy across file systems) ; 0 = linked
    # runs as an `if` condition (no errexit): every step that can fail returns 1 = compute instead
    local S=$1 m=$2 g o of ofp ost orows omol opse oaln dst t
    [[ -n $OTHER ]] && (( ! FORCE )) || return 1      # no mode to link from (mag-mag), or --force: compute
    g=${UMAGS[${ASM[$S]}|$m]}
    [[ $g != *,* ]] || return 1
    of="$MMT_OUT/$OTHER/stats/$S.tsv"; [[ -f $of ]] || return 1
    o=$(res_file "$OTHER" "$m" "$S" "$g"); [[ -s $o ]] || return 1
    IFS=$'\t' read -r ost ofp omol opse oaln orows < <(awk -F'\t' -v m="$m" -v g="$g" '$1 == m && $2 == g {print $3 "\t" $4 "\t" $6 "\t" $7 "\t" $8 "\t" $9; exit}' "$of")
    unit_fp "$OTHER" "$S" "$m"
    [[ ${ost:-} == done || ${ost:-} == linked ]] && [[ $ofp == "$UFP" ]] || return 1
    dst=$(res_file "$MODE" "$m" "$S" "$g"); mkdir -p "$(dirname "$dst")" || return 1
    if [[ ! $o -ef $dst ]]; then                       # already the same file (linked before): nothing to do
        t=$(tmpname "$dst")
        { ln -f "$o" "$t" 2> /dev/null || cp "$o" "$t"; } && publish "$t" "$dst" || { rm -f "$t"; return 1; }
    fi
    [[ $MODE == mag-rep ]] || g=$m                    # the unit key: the MAG (mag-rep) or the MGS (mgs-rep)
    MOL[$g]=$omol PSE[$g]=$opse ALN[$g]=$oaln ROWS[$g]=$orows
    return 0
}

ref_dir() { echo "$TMPROOT/refs/$1.$(md5str "${REP_FASTA[$1]}" | cut -c1-8)"; }   # per MGS + FASTA (shared --tmp)

ensure_ref() {   # ensure_ref MGS : ref_dir/{rep.fa, rep.fa.fai, rep.mmi} from GeneCat's representative FASTA
    local m=$1 r tmp
    r=$(ref_dir "$1")
    [[ $(cat "$r/stamp" 2> /dev/null) == "${REP_STAMP[$m]}" && -s $r/rep.mmi ]] && return 0
    tmp=$(tmpname "$TMPROOT/refs/$m"); rm -rf "$tmp"; mkdir -p "$tmp"
    gzip -dcf "${REP_FASTA[$m]}" | awk '/^>/ {print $1; next} {print toupper($0)}' > "$tmp/rep.fa"
    "$SAMTOOLS" faidx "$tmp/rep.fa"
    "$PBMM2" index --preset CCS -j "$THREADS" "$tmp/rep.fa" "$tmp/rep.mmi"
    echo "${REP_STAMP[$m]}" > "$tmp/stamp"
    [[ -d $r && $(cat "$r/stamp" 2> /dev/null) != "${REP_STAMP[$m]}" ]] && rm -rf "$r"   # an outdated cache entry
    mv -T "$tmp" "$r" 2> /dev/null || rm -rf "$tmp"                                     # another job was faster
    [[ -s $r/rep.mmi ]] || die "reference index for $m missing after build ($r)"
}

write_layer() {   # write_layer S L : the kinetics table of layer L -> one parquet per unit (temp name, checked, renamed)
    local S=$1 L=$2 K n_ac n_out n_nomgs k m n; K="$W/k$L/k.duckdb"
    if [[ ! -f $W/k$L/written.ok ]]; then
        # unit of each MGS in this layer: its MAG(s) (MAG column) and its file
        awk -F'\t' -v L="$L" 'NR == FNR {u[$1]; next} $3 == L && ($1 in u) {print $1 "\t" $2}' "$W/layer$L.txt" "$W/keys.tsv" | while IFS=$'\t' read -r k m; do
            if [[ $MODE != mgs-rep ]]; then printf '%s\t%s\t%s\n' "$m" "$k" "$m/$S/$k.parquet"
            else printf '%s\t%s\t%s\n' "$m" "${UMAGS[${ASM[$S]}|$m]}" "$m/$S.parquet"; fi
        done > "$W/k$L/units.tsv"
        cut -f1,3 "$W/k$L/units.tsv" > "$W/k$L/paths.tsv"
        local rows="(k.base IN ('A', 'C') OR k.call <> 'none') AND k.coverage >= $MGSCALL_MIN_ROW_COV"
        #   contig -> MGS: the representatives' contigs, or (mag-mag) the contigs of this layer's MAGs (one MAG per MGS)
        local ctg="SELECT MGS, contig FROM read_csv('$IDX/rep_contigs.tsv.gz', delim='\t', header=true, all_varchar=true) WHERE MGS IN (SELECT MGS FROM u)"
        [[ $MODE != mag-mag ]] || ctg="SELECT MGS, contig FROM read_csv('$W/contigs.tsv', delim='\t', header=false,
            columns={'contig':'VARCHAR','MAG':'VARCHAR','MGS':'VARCHAR'}) WHERE MAG IN (SELECT MAG FROM u)"
        "$DUCKDB" -bail "$K" -c "SET threads=$THREADS; SET memory_limit='$DUCK_MEM'; SET temp_directory='$W/k$L/duck_tmp';
            DROP TABLE IF EXISTS ctg; DROP TABLE IF EXISTS u; DROP TABLE IF EXISTS out;
            CREATE TABLE u AS SELECT * FROM read_csv('$W/k$L/units.tsv', header=false, delim='\t', columns={'MGS':'VARCHAR','MAG':'VARCHAR','path':'VARCHAR'});
            CREATE TABLE ctg AS $ctg;
            CREATE TABLE out AS
                SELECT c.MGS, u.MAG, k.contig, k.pos::INTEGER AS pos, CASE k.strand WHEN 0 THEN '+' ELSE '-' END AS strand, k.base,
                       k.coverage, k.score, k.ipdRatio, k.tMean, k.tErr, k.modelPrediction, k.frac, k.fracLow, k.fracUp,
                       k.call, k.identificationQv, k.pile_depth, k.n_ref, k.top_alt, k.n_alt,
                       k.ctx_alt_pct, k.ctx_alt_n, k.clip_ctx_n::USMALLINT AS clip_ctx_n
                FROM k LEFT JOIN ctg c ON c.contig = k.contig LEFT JOIN u ON u.MGS = c.MGS
                WHERE $rows;"
        read -r n_ac n_out n_nomgs < <("$DUCKDB" -readonly "$K" -noheader -list -c "SELECT
            (SELECT count(*) FROM k WHERE $rows) || ' ' || (SELECT count(*) FROM out) || ' ' ||
            (SELECT count(*) FROM out WHERE MGS IS NULL OR MAG IS NULL)")
        [[ $n_ac == "$n_out" && $n_nomgs == 0 ]] || die "$S: kept rows $n_ac, table rows $n_out, rows without MGS / MAG $n_nomgs"
        "$DUCKDB" -readonly "$K" -noheader -csv -separator $'\t' -c "SELECT MGS, count(*) FROM out GROUP BY MGS ORDER BY MGS" > "$W/k$L/rows.tsv"
        # one sorted export, split per MGS by lib/mmtpq.py (DuckDB's partitioned write neither keeps the order inside a
        # partition nor one file per partition)
        "$DUCKDB" -bail -readonly "$K" -c "SET threads=$THREADS; SET memory_limit='$DUCK_MEM'; SET temp_directory='$W/k$L/duck_tmp';
            COPY (SELECT MGS, MAG, '$S' AS sample, * EXCLUDE (MGS, MAG) FROM out ORDER BY MGS, contig, pos, strand)
            TO '$W/k$L/export.parquet' (FORMAT parquet, COMPRESSION zstd)"
        "$PYTHON" "$MMT_HOME/lib/mmtpq.py" split --input "$W/k$L/export.parquet" --out "$MMT_OUT/$MODE" --sample "$S" \
            --paths "$W/k$L/paths.tsv" | sort > "$W/k$L/written.tsv"
        cmp -s "$W/k$L/written.tsv" <(sort "$W/k$L/rows.tsv") || die "$S: rows written per MGS differ from the kinetics table: $(diff "$W/k$L/written.tsv" <(sort "$W/k$L/rows.tsv") | head -3 | tr '\n' ' ')"
        rm -f "$W/k$L/export.parquet"
        touch "$W/k$L/written.ok"
        awk -F'\t' -v OFS='\t' -v L="$L" '{if (L > 1) $1 = $1 "_layer" L; print}' "$W/k$L/stats.tsv" >> "$W/run.stats"
    fi
    while IFS=$'\t' read -r m n; do
        k=$m; [[ $MODE != mgs-rep ]] && k=$(awk -F'\t' -v m="$m" '$1 == m {print $2; exit}' "$W/k$L/units.tsv")
        ROWS[$k]=$n
    done < "$W/k$L/rows.tsv"
}

write_state() {   # write_state S [W] [T0] : stats/<S>.tsv = old lines, minus TODO (replaced) and STALE (dropped with their
                  # result) ; stats/<S>.run.tsv ; one SAMPLE line in the mmt log. Uses TODO STALE FP LINKED and MOL PSE ALN
                  # ROWS (per unit: the MGS in mgs-rep, the MAG in mag-rep / mag-mag)
    local S=$1 W=${2:-} ts=${3:-} m g u f x n st mol pseudo aln rows now t ab
    local -a out=() parts=()
    now=$(date '+%F %T')
    if [[ -n $W && -f $W/molecules.tsv ]]; then
        while IFS=$'\t' read -r u n; do MOL[$u]=$n; done < <(tail -n +2 "$W/molecules.tsv")
    fi
    if [[ -n $W && -f $W/aln/counts.tsv ]]; then
        while IFS=$'\t' read -r u mol pseudo aln; do PSE[$u]=$pseudo; ALN[$u]=$aln; done < "$W/aln/counts.tsv"
    fi
    for m in "${STALE[@]}"; do
        if [[ $MODE != mgs-rep ]]; then rm -rf "${D:?}/$m/$S"; else rm -f "$D/$m/$S.parquet"; fi
        rmdir "$D/$m" 2> /dev/null || true
        log "$S: $m dropped (no longer usable, or no MAG of it in this sample's assembly)"
    done
    t=$(tmpname "$STATE/$S.tsv")
    {
        printf 'MGS\tMAG\tstatus\tfingerprint\tabundance\tmolecules\tpseudo_reads\taligned\trows\trep_MAG\tupdated\tmmt\n'
        if [[ -f $STATE/$S.tsv ]]; then
            awk -F'\t' -v drop=" ${TODO[*]} ${STALE[*]} " 'FNR > 1 && index(drop, " " $1 " ") == 0' "$STATE/$S.tsv"
        fi
        for m in "${TODO[@]}"; do
            ab=${ABUND[$S|$m]:-0}
            parts=()
            if [[ $MODE != mgs-rep ]]; then
                for f in "$D/$m/$S"/*.parquet; do      # files of MAGs no longer in this sample's assembly
                    [[ -f $f && ,${UMAGS[${ASM[$S]}|$m]}, != *",$(basename "$f" .parquet),"* ]] && rm -f "$f"
                done
                for g in $(tr ',' '\n' <<< "${UMAGS[${ASM[$S]}|$m]}" | sort); do
                    state_row "$m" "$g" "$g"
                done
                rmdir "$D/$m/$S" 2> /dev/null || true
            else
                state_row "$m" "${UMAGS[${ASM[$S]}|$m]}" "$m"
            fi
            x=""; for g in "${parts[@]}"; do x+="${x:+, }$g"; done
            out+=("$m $x")
        done
    } > "$t"
    publish "$t" "$STATE/$S.tsv"
    if (( ${#TODO[@]} )); then                         # the sample's last run (also when every table was linked)
        t=$(tmpname "$STATE/$S.run.tsv")
        { printf 'sample\t%s\nmode\t%s\nupdated\t%s\nmmt\t%s\nraw_bam\t%s\ncram\t%s\n' "$S" "$MODE" "$now" "$MMT_VERSION" "${RAW[$S]}" "${CRAM[$S]}"
          awk -F'\t' -v s="$S" 'NR == 1 {for (i = 2; i <= NF; i++) h[i] = $i; next} $1 == s {for (i = 2; i <= NF; i++) print h[i] "\t" $i}' "$IDX/movies.tsv"
          printf 'ipdsummary_args\t%s\n' "${MMT_IPD_ARGS:-}"
          [[ -z $OTHER ]] || printf 'linked_from_%s\t%s\n' "$OTHER" "${!LINKED[*]}"
          cat "$IDX/reads/_stats/$S.tsv" ${W:+"$W/run.stats"} 2> /dev/null || true; } > "$t"
        publish "$t" "$STATE/$S.run.tsv"
    fi
    for m in "${STALE[@]}"; do out+=("$m dropped (no MAG of it in this sample any more, or no longer usable)"); done
    log "$S: state: $(awk -F'\t' -v todo=" ${TODO[*]} " 'FNR > 1 && index(todo, " " $1 " ") {n[$3]++} END {for (s in n) printf "%s %d  ", s, n[s]}' "$STATE/$S.tsv")${STALE[*]:+ dropped: ${STALE[*]}}"
    local line="" x
    for x in "${out[@]}"; do line+="${line:+ ; }$x"; done
    mlog SAMPLE "$S: $line${ts:+ ; $(( ($(date +%s) - ts + 59) / 60 )) min}"
}

state_row() {   # state_row MGS MAG UNIT : one stats line (stdout) + its text in parts[] (write_state's locals)
    local m=$1 g=$2 u=$3 f
    mol=${MOL[$u]:-0} pseudo=${PSE[$u]:-0} aln=${ALN[$u]:-0} rows=${ROWS[$u]:-0}
    f=$(res_file "$MODE" "$m" "$S" "$g")
    if [[ -n ${LINKED[$m]:-} ]]; then st=linked
    elif (( rows > 0 )); then st=done
    else
        if (( mol == 0 )); then st=no_reads; elif (( pseudo == 0 )); then st=no_kinetics
        elif (( aln == 0 )); then st=no_aligned; else st=no_rows; fi
        rm -f "$f"                                    # an older result for this unit no longer holds
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$m" "$g" "$st" "${FP[$m]}" "$ab" "$mol" "$pseudo" "$aln" "$rows" "${REP_MAG[$m]}" "$now" "$MMT_VERSION"
    local pre=""; [[ $MODE == mgs-rep ]] || pre="$g "
    if [[ $st == done || $st == linked ]]; then parts+=("$pre$st ($rows sites, $mol molecules)")
    else parts+=("$pre$st ($mol molecules, $pseudo strand reads, $aln aligned)"); fi
}
