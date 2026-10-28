# shellcheck shell=bash
# MMT MGScall - methylation of each MGS on its representative MAG, one table per MGS x sample.
# Sourced by the mmt wrapper. Method, outputs and columns: docs/MGScall.md (copied into <OUT>/MGScall/README.md).
# shellcheck source=../lib/kinetics.sh
source "$MMT_HOME/lib/kinetics.sh"
# shellcheck source=../lib/reads.sh
source "$MMT_HOME/lib/reads.sh"

MGSCALL_METHOD=1   # bump when a code change alters MGScall results -> the next run redoes every result

module_help() {
    cat << 'EOF'
mmt MGScall - methylation of each MGS on its representative MAG, per sample (ipdSummary m6A + m4C)

Usage: mmt MGScall (--mgs-rep | --mag-rep) -o OUT [options]

Mode (one per run)
  --mgs-rep           reads of EVERY MAG of the MGS (each sample that has a MAG of it), aligned to the representative
  --mag-rep           reads of the REPRESENTATIVE MAG only (its own sample, or its co-assembly samples)
Which reads belong to which MAG comes from MATAFILER4's own alignment (the sample's CRAM): no re-assignment.

Inputs (defaults from config/site.conf)
  -i, --mf4-run DIR   MATAFILER4 run folder (default: the map's #OutPath/#RunID)
  --genecat DIR       GeneCat folder (with Bin_<binner>/MAGvsGC.txt.gz)
  --map FILE          MATAFILER4 map (raw kinetics BAMs: #DirPath + SmplPrefix + .bam)
Output and work space
  -o, --out DIR       results folder - reuse it: finished results are kept, new ones added, changed inputs redone
  --tmp DIR           work space (default $TMPDIR, else /tmp); a shared scratch folder also caches the reference indexes
Selection
  --select IDS        only these MGS: "MGS.12,MGS.40" or a file (one ID per line); default all usable MGS
  --samples IDS       only these samples (same format); default every usable sample with work
                      HPC: one sample per array task is the efficient unit (each raw BAM is read once per run)
Resources (default: the SLURM allocation)
  -t, --threads N     CPUs to use
  --mem GB            memory to plan for (DuckDB, sorting)
Other
  --dry-run           show what would run (per sample: MGS to do / up to date) and stop
  --list-samples      print the samples that have work (one per line, nothing else) and stop - for SLURM arrays
  --force             redo the selected results even if up to date
  --debug             keep each sample's work folder (success or failure; a rerun resumes from it). Without it the
                      work folder is removed when the sample ends, whatever happened (--keep-tmp: the same)
  --config FILE       extra settings file (after config/site.conf)

Presence check: a sample is planned for an MGS only when its assembly has a MAG of the MGS (mag-rep: the
representative). Before its raw BAM is read, the molecules MATAFILER4 placed on those MAGs are counted: an MGS with
fewer than MGSCALL_PRESENT_MOL is recorded as 'absent' and not processed; when every selected MGS of a sample is
absent its raw BAM is not read at all.

Output (in OUT)
  mmt.log                            the folder's log: every run, its plan per MGS, each sample x MGS outcome, failures
  MGScall/<mode>/<MGS>/<S>.parquet   one row per A/C site x strand of the representative (columns MGS, sample, ...)
  MGScall/<mode>/_state/<S>.tsv      what is done for each MGS of that sample (+ read counts)
  index/                             samples, MGS, MAGs, contigs; index/reads/<S>.parquet = read -> MAG
  logs/MGScall/<mode>/               one detailed log per run (the terminal gets a start, the plan and an end line)
Read all of a mode with DuckDB: SELECT * FROM read_parquet('OUT/MGScall/mgs-rep/*/*.parquet')
Results of MMT <= 0.2 (MGS=<MGS>/sample=<S>/data.parquet) are moved into this layout by the first run (no recompute).

Examples
  one job, all samples one after another:
    mmt MGScall --mgs-rep -o OUT --select MGS.12,MGS.40
  SLURM array, one sample per task (list once, before submitting):
    mmt MGScall --mgs-rep -o OUT --select MGS.12,MGS.40 --list-samples > todo.txt
    sbatch --array=1-$(wc -l < todo.txt)%20 job.sh
    # job.sh:  S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" todo.txt)
    #          mmt MGScall --mgs-rep -o OUT --select MGS.12,MGS.40 --samples "$S" --tmp "$TMPDIR"
  what is done per MGS (printed and added to OUT/mmt.log):
    mmt status -o OUT --select MGS.12,MGS.40
EOF
}

module_main() {
    MODE=""
    local SELECT="" SAMPLES="" FORCE=0 DRY=0 LIST=0 t0 S msg n_ok=0 n_skip=0 err desc
    local -a failed=() todo_samples=()
    MMT_DEBUG=${MMT_DEBUG:-0}
    while (( $# )); do
        if common_opt "$@"; then shift "$COMMON_SHIFT"; continue; fi
        case $1 in
            --mgs-rep|--mag-rep)
                [[ -z $MODE || $MODE == "${1#--}" ]] || usage_die "one mode per run: --mgs-rep or --mag-rep"
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
    [[ -n $MODE ]] || usage_die "choose a mode: --mgs-rep or --mag-rep"
    [[ $MGSCALL_PRESENT_MOL =~ ^[0-9]+$ ]] || usage_die "MGSCALL_PRESENT_MOL='$MGSCALL_PRESENT_MOL' must be a whole number"
    (( LIST )) && MMT_QUIET=1
    t0=$(date +%s)
    resolve_out
    job_tag
    start_log "$MMT_OUT/logs/MGScall/$MODE/$(date +%Y%m%d-%H%M%S)_$JOBTAG.log"
    say "MMT $(mmt_version) MGScall --$MODE -> $MMT_OUT (log: $MMT_LOG)"
    log "command: mmt MGScall --$MODE ${SELECT:+--select $SELECT }${SAMPLES:+--samples $SAMPLES }(force=$FORCE dry=$DRY list=$LIST debug=$MMT_DEBUG) on $(hostname), job $JOBTAG"
    desc="${SELECT:+--select $SELECT }${SAMPLES:+--samples $SAMPLES }"
    (( FORCE )) && desc+="--force "
    (( LIST )) && desc+="(listing the samples with work) "
    (( MMT_DEBUG )) && desc+="--debug "
    (( DRY )) || mlog_start "MGScall --$MODE" "${desc% }"
    # a run that stops early (an error, scancel, the time limit) says so in OUT/mmt.log
    trap 'run_exit $?' EXIT
    trap 'exit 143' TERM INT
    resolve_inputs
    resolve_resources
    setup_tools
    TMPROOT=$(abspath "${MMT_TMP:-${TMPDIR:-/tmp}}")/mmt
    mkdir -p "$TMPROOT" || die "cannot create work space $TMPROOT"
    work_root
    DUCK_MEM="$(( MEM_GB * 6 / 10 > 1 ? MEM_GB * 6 / 10 : 1 ))GB"
    SORT_T=$(( THREADS < 4 ? THREADS : 4 ))
    SORT_MEM="$(( MEM_GB * 1024 * 3 / 10 / SORT_T > 256 ? MEM_GB * 1024 * 3 / 10 / SORT_T : 256 ))M"
    MIN_MAPQ=$MGSCALL_MIN_MAPQ CTX=$MGSCALL_CTX CLIP_MIN=$MGSCALL_CLIP_MIN
    SETTINGS="keepdup=$MGSCALL_KEEP_DUP id=$MGSCALL_MIN_IDENTITY mapq=$MGSCALL_MIN_MAPQ pct=$MGSCALL_MIN_PCT_ALIGNED bp=$MGSCALL_MIN_ALN_BP rowcov=$MGSCALL_MIN_ROW_COV ctx=$MGSCALL_CTX clip=$MGSCALL_CLIP_MIN present=$MGSCALL_PRESENT_MOL"
    log "inputs: MF4 run $MMT_MF4_RUN | GeneCat $MMT_GENECAT (Bin_$MMT_BINNER) | map $MMT_MAP"
    log "resources: $THREADS threads, $MEM_GB GB (DuckDB $DUCK_MEM, sort $SORT_T x $SORT_MEM) | work space $WROOT"
    log "settings: $SETTINGS | method $MGSCALL_METHOD"
    tool_versions | sed 's/^/tool: /'

    ensure_index
    D="$MMT_OUT/MGScall/$MODE"; STATE="$D/_state"
    mkdir -p "$STATE" "$MMT_OUT/logs/MGScall/$MODE/ipdSummary"
    migrate_layout
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
    if [[ -z $SAMPLES ]]; then plan_log; fi       # per MGS: where it is present, to do / up to date (whole-selection runs)
    if (( LIST )); then
        (( ${#todo_samples[@]} == 0 )) || printf '%s\n' "${todo_samples[@]}" >&3
        mlog END "listed ${#todo_samples[@]} sample(s) with work ($n_units MGS x sample results) - nothing processed by a listing run"
        trap - EXIT; return 0
    fi
    if (( DRY )); then say "dry run - nothing done"; trap - EXIT; return 0; fi

    write_readmes
    for S in "${todo_samples[@]}"; do
        isolated process_sample "$S"
        if (( RC == 0 )); then n_ok=$(( n_ok + 1 )); else
            failed+=("$S"); err=$(last_error)
            log "== $S FAILED (exit $RC)"
            mlog FAIL "$S: ${err:-exit $RC} - details: ${MMT_LOG#"$MMT_OUT"/}$( (( MMT_DEBUG )) && echo " ; work folder kept: $WROOT/MGScall_$MODE/$S")"
        fi
    done
    mlog END "$n_ok sample(s) processed, ${#failed[@]} failed$( (( ${#failed[@]} )) && echo " ($(id_list "${failed[@]}"))"), $n_skip already up to date ; $(( ($(date +%s) - t0) / 60 )) min"
    trap - EXIT
    if (( ${#failed[@]} )); then
        say "finished with errors after $(( ($(date +%s) - t0) / 60 )) min: ${#failed[@]} sample(s) failed: ${failed[*]} - rerun the same command to resume"
        exit 1
    fi
    say "done in $(( ($(date +%s) - t0) / 60 )) min: $n_ok sample(s) processed, $n_skip already up to date"
}

# ============================================================================ selection + plan
pick_mgs() {   # MGS_SEL (usable, selected) ; REP_FASTA[MGS] ; REP_STAMP[MGS] ; OK_MGS[MGS] ; IS_SEL[MGS]
    declare -gA REP_FASTA=() REP_STAMP=() OK_MGS=() IS_SEL=()
    declare -ga MGS_SEL=()
    local m st fa asm want=() bad=() unk=()
    declare -A seen=()
    while IFS=$'\t' read -r m st _ _ asm fa _; do
        seen[$m]=$st
        [[ $st == ok ]] || continue
        OK_MGS[$m]=1; REP_FASTA[$m]=$fa
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
                                   mlog WARN "skipped, representative not usable: ${bad[*]} (see index/report.txt)"; }
    else
        mapfile -t MGS_SEL < <(awk -F'\t' 'NR > 1 && $2 == "ok" {print $1}' "$IDX/mgs.tsv")
    fi
    (( ${#MGS_SEL[@]} )) || die "no usable MGS selected"
    for m in "${MGS_SEL[@]}"; do REP_STAMP[$m]=$(stamp "${REP_FASTA[$m]}"); IS_SEL[$m]=1; done
}

pick_samples() {   # SAMPLE_SEL ; OKS (every usable sample) ; RAW[S] CRAM[S] ASM[S] ; UNITS[S] = the selected MGS that
                   # sample can give reads to ; UMD5[assembly|MGS] = checksum of the contigs those reads come from ; ASM_MD5
    declare -gA RAW=() CRAM=() ASM=() UNITS=() UMD5=() ASM_MD5=()
    declare -ga SAMPLE_SEL=() OKS=()
    local s st raw cram asm m here mall mrep want=() bad=() unk=()
    declare -A seen=() asm_units=()
    #   mgs-rep: every MGS with a MAG in that assembly ; mag-rep: MGS whose representative comes from that assembly
    while IFS=$'\t' read -r asm m _ here mall mrep; do
        [[ -n ${IS_SEL[$m]:-} ]] || continue
        if [[ $MODE == mgs-rep ]]; then asm_units[$asm]+="$m "; UMD5[$asm|$m]=$mall
        elif [[ $here == 1 ]]; then asm_units[$asm]+="$m "; UMD5[$asm|$m]=$mrep; fi
    done < <(tail -n +2 "$IDX/units.tsv")
    while IFS=$'\t' read -r asm mall; do ASM_MD5[$asm]=$mall; done < <(tail -n +2 "$IDX/assemblies.tsv")
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
                                   mlog WARN "skipped, sample not usable: ${bad[*]} (see index/report.txt)"; }
    else
        while IFS=$'\t' read -r s st _; do   # samples with work, or with older results that may need dropping
            [[ $st == ok ]] && [[ -n ${UNITS[$s]:-} || -f $STATE/$s.tsv ]] && SAMPLE_SEL+=("$s")
        done < <(tail -n +2 "$IDX/samples.tsv")
    fi
    return 0
}

sample_fp() {   # what every result of a sample depends on (besides its MGS's representative and contigs)
    md5str "$MGSCALL_METHOD|$MODE|$SETTINGS|${CRAM[$1]}:$(stamp "${CRAM[$1]}")|${RAW[$1]}:$(stamp "${RAW[$1]}")"
}

plan_sample() {   # plan_sample S -> TODO (MGS to (re)do), CURRENT (up to date), STALE (old results to drop), N_CURRENT ;
                  # FP[MGS] = new fingerprint
    local S=$1 sfp m st fp
    declare -gA FP=()
    declare -ga TODO=() CURRENT=() STALE=()
    declare -A old_fp=() old_st=()
    N_CURRENT=0
    sfp=$(sample_fp "$S")
    if [[ -f $STATE/$S.tsv ]]; then
        while IFS=$'\t' read -r m st fp _; do old_fp[$m]=$fp; old_st[$m]=$st; done < <(tail -n +2 "$STATE/$S.tsv")
    fi
    for m in ${UNITS[$S]:-}; do
        FP[$m]="$sfp:${REP_STAMP[$m]}:${UMD5[${ASM[$S]}|$m]}"
        if (( ! FORCE )) && [[ ${old_fp[$m]:-} == "${FP[$m]}" ]] \
           && [[ ${old_st[$m]} != "done" || -s $D/$m/$S.parquet ]]; then
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
    local W LOCK="$STATE/$S.lock" t m mol pseudo n f mof repf="" ts
    ts=$(date +%s)
    lock_take "$LOCK" || die "$S: another job is working on this sample ($(cat "$LOCK/owner" 2> /dev/null)) - remove $LOCK if not"
    W="$WROOT/MGScall_$MODE/$S"               # this sample's work folder: ours while we hold the lock
    # whatever ends this sample (done, error, scancel / time limit): work folder removed (kept with --debug), lock released
    trap 'rc=$?; work_done "'"$W"'" $(( rc == 0 )); lock_release "'"$LOCK"'"' EXIT
    trap 'exit 143' TERM INT
    plan_sample "$S"                       # again, under the lock (another job may have finished some)
    if (( ${#TODO[@]} == 0 )); then
        if (( ${#STALE[@]} )); then write_state "$S" "" "$ts"; else log "== $S: up to date (another job finished it)"; fi
        return 0
    fi
    log "== $S: ${#TODO[@]} MGS to do ($MODE) | raw ${RAW[$S]} | CRAM ${CRAM[$S]}"
    local wfp
    wfp=$(md5str "${TODO[*]}|$(for m in "${TODO[@]}"; do printf '%s ' "${FP[$m]}"; done)")
    if [[ $(cat "$W/work.fp" 2> /dev/null) != "$wfp" ]] || (( FORCE )); then rm -rf "$W"; fi   # other work than last time, or --force: start clean
    [[ -d $W ]] && log "$S: resuming from the work folder of an earlier run: $W"
    mkdir -p "$W"; echo "$wfp" > "$W/work.fp"
    [[ -f $W/run.stats ]] || : > "$W/run.stats"

    # 1. read -> MAG -> MGS for this sample (shared by both modes; kept in index/reads)
    t=$(date +%s); ensure_reads "$S" "$W"; log "$S: reads index ok ($(( $(date +%s) - t )) s)"

    # 2. presence + the reads to use: per MGS to do, the molecules MATAFILER4 placed on its MAG(s) in this sample
    #    (mag-rep: on the representative only) ; fewer than MGSCALL_PRESENT_MOL = absent, not processed
    if [[ ! -f $W/select.ok ]]; then
        printf '%s\n' "${TODO[@]}" > "$W/todo.txt"
        [[ $MODE == mag-rep ]] && repf="AND is_rep"
        (( MGSCALL_KEEP_DUP )) || repf+=" AND NOT dup"
        "$DUCKDB" -bail -c "SET threads=$THREADS; SET memory_limit='$DUCK_MEM';
            CREATE TEMP TABLE s AS SELECT read, MGS FROM read_parquet('$READS_PQ')
                WHERE MGS IN (SELECT m FROM read_csv('$W/todo.txt', header=false, columns={'m':'VARCHAR'})) $repf;
            COPY (SELECT t.m AS MGS, count(s.read) AS molecules FROM read_csv('$W/todo.txt', header=false, columns={'m':'VARCHAR'}) t
                  LEFT JOIN s ON s.MGS = t.m GROUP BY t.m ORDER BY t.m) TO '$W/molecules.tsv' (HEADER, DELIMITER '\t');
            COPY (SELECT read, MGS FROM s WHERE MGS IN (SELECT MGS FROM s GROUP BY MGS HAVING count(*) >= $MGSCALL_PRESENT_MOL)
                  ORDER BY read) TO '$W/sel.tsv' (HEADER, DELIMITER '\t')"
        tail -n +2 "$W/sel.tsv" | cut -f1 > "$W/names.txt"
        printf 'selected_molecules\t%s\n' "$(wc -l < "$W/names.txt")" >> "$W/run.stats"
        log "$S: presence (molecules on the MGS's MAGs in this sample, >= $MGSCALL_PRESENT_MOL = present): $(awk -F'\t' -v k="$MGSCALL_PRESENT_MOL" \
            'NR > 1 {printf "%s %d%s  ", $1, $2, ($2 >= k ? "" : " ABSENT")}' "$W/molecules.tsv")"
        [[ -s $W/names.txt ]] || log "$S: every MGS to do is absent - raw BAM not read"
        touch "$W/select.ok"
    fi

    # 3. those reads from the raw kinetics BAM (one pass) -> one pseudo-read per strand -> one BAM per MGS
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
        rm -rf "$W/byMGS"
        mof=$(( $(ulimit -n) - 100 )); (( mof > 2000 )) && mof=2000; (( mof < 50 )) && mof=50   # BAMs open at once
        "$PYTHON" "$MMT_HOME/lib/mf4.py" split --bam "$W/bys.bam" --reads "$W/sel.tsv" --outdir "$W/byMGS" \
            --threads "$SORT_T" --max-open "$mof"
        rm -f "$W/bys.bam"
        log "$S: $found molecules extracted and split into $(( $(wc -l < "$W/byMGS/counts.tsv") - 1 )) MGS ($(( $(date +%s) - t )) s)"
        touch "$W/extract.ok"
    fi

    # 4. each MGS's reads aligned to its own representative only (MATAFILER4 already placed them)
    #    kept: primary, MAPQ >= MIN_MAPQ, identity >= MIN_IDENTITY (pbmm2 -y), no hard clip, >= MIN_PCT_ALIGNED % of the
    #    read and >= MIN_ALN_BP aligned
    if [[ -f $W/extract.ok && ! -f $W/align.ok ]]; then
        t=$(date +%s); mkdir -p "$W/aln"; : > "$W/aln/counts.tsv"
        while IFS=$'\t' read -r m mol pseudo; do
            ensure_ref "$m"
            "$PBMM2" align --preset CCS --sample "$S" -y "$MGSCALL_MIN_IDENTITY" -j "$THREADS" \
                "$(ref_dir "$m")/rep.mmi" "$W/byMGS/$m.bam" "$W/aln/$m.raw.bam"
            "$SAMTOOLS" view -u -F 2308 -q "$MGSCALL_MIN_MAPQ" "$W/aln/$m.raw.bam" \
                | "$SAMTOOLS" view -u -e "hclen==0 && (qlen-sclen)*100 >= $MGSCALL_MIN_PCT_ALIGNED*qlen && (qlen-sclen) >= $MGSCALL_MIN_ALN_BP" - \
                | "$SAMTOOLS" sort -@ "$SORT_T" -m "$SORT_MEM" -T "$W/aln/$m.srt" -o "$W/aln/$m.bam" -
            n=$("$SAMTOOLS" view -c "$W/aln/$m.bam")
            rm -f "$W/aln/$m.raw.bam"
            (( n > 0 )) || rm -f "$W/aln/$m.bam"
            printf '%s\t%s\t%s\t%s\n' "$m" "$mol" "$pseudo" "$n" >> "$W/aln/counts.tsv"
        done < <(tail -n +2 "$W/byMGS/counts.tsv")
        rm -rf "$W/byMGS"
        log "$S: aligned $(awk -F'\t' '{s += $4} END {print s + 0}' "$W/aln/counts.tsv") pseudo-reads to $(awk -F'\t' '$4 > 0' "$W/aln/counts.tsv" | wc -l) representatives ($(( $(date +%s) - t )) s)"
        touch "$W/align.ok"
    fi

    # 5. one BAM + one reference for the sample -> ipdSummary + pileup once (lib/kinetics.sh)
    if [[ -f $W/align.ok && ! -f $W/kinetics.ok ]] && awk -F'\t' '$4 > 0 {f = 1} END {exit !f}' "$W/aln/counts.tsv"; then
        t=$(date +%s)
        if [[ ! -f $W/merge.ok ]]; then
            awk -F'\t' -v d="$W/aln" '$4 > 0 {print d "/" $1 ".bam"}' "$W/aln/counts.tsv" > "$W/bams.txt"
            "$SAMTOOLS" merge -@ "$THREADS" -c -p -f -o "$W/all.bam" -b "$W/bams.txt"
            "$SAMTOOLS" index -@ "$THREADS" "$W/all.bam"
            "$PBINDEX" "$W/all.bam"
            : > "$W/ref.fa"
            while read -r f; do m=$(basename "$f" .bam); ensure_ref "$m"; cat "$(ref_dir "$m")/rep.fa" >> "$W/ref.fa"; done < "$W/bams.txt"
            "$SAMTOOLS" faidx "$W/ref.fa"
            xargs rm -f < "$W/bams.txt"            # merged: the per-MGS BAMs are no longer needed (disk)
            touch "$W/merge.ok"
        fi
        rm -rf "$W/k"; mkdir -p "$W/k"
        "$SAMTOOLS" idxstats "$W/all.bam" | awk '$3 > 0 {print $1}' > "$W/k/contigs.txt"
        kinetics_table "$W/all.bam" "$W/ref.fa" "$W/k/contigs.txt" "$W/k" "$MMT_OUT/logs/MGScall/$MODE/ipdSummary/$S.log"
        log "$S: kinetics table ready ($(( $(date +%s) - t )) s)"
        touch "$W/kinetics.ok"
    fi

    # 6. results: one parquet per MGS ; state of every MGS done in this run (work folder: removed by the EXIT trap)
    write_results "$S" "$W" "$ts"
    log "== $S: done"
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

write_results() {   # write_results S W T0 : one parquet per MGS done (each written to a temp name, checked, renamed) + state
    local S=$1 W=$2 ts=$3 m n mol pseudo aln n_ac n_out n_nomgs
    local -A MOL=() PSE=() ALN=() ROWS=()
    # molecules per MGS (presence count, absent ones included), then what the later steps kept
    if [[ -s $W/molecules.tsv ]]; then
        while IFS=$'\t' read -r m n; do MOL[$m]=$n; done < <(tail -n +2 "$W/molecules.tsv")
    fi
    if [[ -f $W/aln/counts.tsv ]]; then
        while IFS=$'\t' read -r m mol pseudo aln; do PSE[$m]=$pseudo; ALN[$m]=$aln; done < "$W/aln/counts.tsv"
    fi
    if [[ -f $W/kinetics.ok ]]; then
        local K="$W/k/k.duckdb"
        awk -F'\t' '$4 > 0 {print $1}' "$W/aln/counts.tsv" > "$W/k/mgs_done.txt"
        if [[ ! -f $W/k/written.ok ]]; then
            "$DUCKDB" -bail "$K" -c "SET threads=$THREADS; SET memory_limit='$DUCK_MEM'; SET temp_directory='$W/k/duck_tmp';
                DROP TABLE IF EXISTS ctg; DROP TABLE IF EXISTS out;
                CREATE TABLE ctg AS SELECT MGS, contig FROM read_csv('$IDX/rep_contigs.tsv.gz', delim='\t', header=true, all_varchar=true)
                    WHERE MGS IN (SELECT m FROM read_csv('$W/k/mgs_done.txt', header=false, columns={'m':'VARCHAR'}));
                CREATE TABLE out AS
                    SELECT c.MGS, k.contig, k.pos::INTEGER AS pos, CASE k.strand WHEN 0 THEN '+' ELSE '-' END AS strand, k.base,
                           k.coverage, k.score, k.ipdRatio, k.tMean, k.tErr, k.modelPrediction, k.frac, k.fracLow, k.fracUp,
                           k.call, k.identificationQv, k.pile_depth, k.n_ref, k.top_alt, k.n_alt,
                           k.ctx_alt_pct, k.ctx_alt_n, k.clip_ctx_n::USMALLINT AS clip_ctx_n
                    FROM k LEFT JOIN ctg c ON c.contig = k.contig
                    WHERE k.base IN ('A', 'C') AND k.coverage >= $MGSCALL_MIN_ROW_COV;"
            read -r n_ac n_out n_nomgs < <("$DUCKDB" -readonly "$K" -noheader -list -c "SELECT
                (SELECT count(*) FROM k WHERE base IN ('A','C') AND coverage >= $MGSCALL_MIN_ROW_COV) || ' ' ||
                (SELECT count(*) FROM out) || ' ' || (SELECT count(*) FROM out WHERE MGS IS NULL)")
            [[ $n_ac == "$n_out" && $n_nomgs == 0 ]] || die "$S: A/C rows $n_ac, table rows $n_out, rows without MGS $n_nomgs"
            "$DUCKDB" -readonly "$K" -noheader -csv -separator $'\t' -c "SELECT MGS, count(*) FROM out GROUP BY MGS ORDER BY MGS" > "$W/k/rows.tsv"
            # one sorted export, split per MGS by lib/mmtpq.py (DuckDB's partitioned write neither keeps the order inside a
            # partition nor one file per partition)
            "$DUCKDB" -bail -readonly "$K" -c "SET threads=$THREADS; SET memory_limit='$DUCK_MEM'; SET temp_directory='$W/k/duck_tmp';
                COPY (SELECT MGS, '$S' AS sample, * EXCLUDE (MGS) FROM out ORDER BY MGS, contig, pos, strand)
                TO '$W/k/export.parquet' (FORMAT parquet, COMPRESSION zstd)"
            "$PYTHON" "$MMT_HOME/lib/mmtpq.py" split --input "$W/k/export.parquet" --out "$D" --sample "$S" | sort > "$W/k/written.tsv"
            cmp -s "$W/k/written.tsv" <(sort "$W/k/rows.tsv") || die "$S: rows written per MGS differ from the kinetics table: $(diff "$W/k/written.tsv" <(sort "$W/k/rows.tsv") | head -3 | tr '\n' ' ')"
            rm -f "$W/k/export.parquet"
            touch "$W/k/written.ok"
            cat "$W/k/stats.tsv" >> "$W/run.stats" 2> /dev/null || true
        fi
        while IFS=$'\t' read -r m n; do ROWS[$m]=$n; done < "$W/k/rows.tsv"
    fi
    write_state "$S" "$W" "$ts"
}

write_state() {   # write_state S [W] [T0] : _state/<S>.tsv = old lines, minus TODO (replaced) and STALE (dropped with their
                  # result) ; one SAMPLE line in OUT/mmt.log. Uses TODO STALE FP and, from write_results, MOL PSE ALN ROWS
    local S=$1 W=${2:-} ts=${3:-} m st mol pseudo aln rows now t
    local -a out=()
    now=$(date '+%F %T')
    for m in "${STALE[@]}"; do
        rm -f "$D/$m/$S.parquet"; rmdir "$D/$m" 2> /dev/null || true
        log "$S: $m dropped (no longer usable, or no MAG of it in this sample's assembly)"
    done
    t=$(tmpname "$STATE/$S.tsv")
    {
        printf 'MGS\tstatus\tfingerprint\tmolecules\tpseudo_reads\taligned\trows\tupdated\tmmt\n'
        if [[ -f $STATE/$S.tsv ]]; then
            awk -F'\t' -v drop=" ${TODO[*]} ${STALE[*]} " 'FNR > 1 && index(drop, " " $1 " ") == 0' "$STATE/$S.tsv"
        fi
        for m in "${TODO[@]}"; do
            mol=${MOL[$m]:-0} pseudo=${PSE[$m]:-0} aln=${ALN[$m]:-0} rows=${ROWS[$m]:-0}
            if (( rows > 0 )); then st="done"
            else
                if (( mol < MGSCALL_PRESENT_MOL || mol == 0 )); then st=absent; elif (( pseudo == 0 )); then st=no_kinetics
                elif (( aln == 0 )); then st=no_aligned; else st=no_rows; fi
                rm -f "$D/$m/$S.parquet"                   # an older result for this MGS no longer holds
            fi
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$m" "$st" "${FP[$m]}" "$mol" "$pseudo" "$aln" "$rows" "$now" "$MMT_VERSION"
            if [[ $st == done ]]; then out+=("$m done ($rows sites, $mol molecules)")
            elif [[ $st == absent ]]; then out+=("$m absent ($mol molecules < $MGSCALL_PRESENT_MOL)")
            else out+=("$m $st ($mol molecules, $pseudo strand reads, $aln aligned)"); fi
        done
    } > "$t"
    publish "$t" "$STATE/$S.tsv"
    if [[ -n $W ]]; then
        t=$(tmpname "$STATE/$S.stats.tsv")
        { printf 'sample\t%s\nmode\t%s\nupdated\t%s\n' "$S" "$MODE" "$now"
          cat "$MMT_OUT/index/reads/_stats/$S.tsv" "$W/run.stats" 2> /dev/null || true; } > "$t"
        publish "$t" "$STATE/$S.stats.tsv"
    fi
    for m in "${STALE[@]}"; do out+=("$m dropped (no MAG of it in this sample any more, or no longer usable)"); done
    log "$S: state: $(awk -F'\t' -v todo=" ${TODO[*]} " 'FNR > 1 && index(todo, " " $1 " ") {n[$2]++} END {for (s in n) printf "%s %d  ", s, n[s]}' "$STATE/$S.tsv")${STALE[*]:+ dropped: ${STALE[*]}}"
    local line="" x
    for x in "${out[@]}"; do line+="${line:+ ; }$x"; done
    mlog SAMPLE "$S: $line${ts:+ ; $(( ($(date +%s) - ts + 59) / 60 )) min}"
}

write_readmes() {   # OUT/README.md (folders) and OUT/MGScall/README.md (method + columns), refreshed every run
    local t
    t=$(tmpname "$MMT_OUT/MGScall/README.md"); cp "$MMT_HOME/docs/MGScall.md" "$t"; publish "$t" "$MMT_OUT/MGScall/README.md"
    write_out_readme
}

migrate_layout() {   # MMT <= 0.2 wrote MGS=<MGS>/sample=<S>/data.parquet: rewritten once as <MGS>/<S>.parquet with the
                     # columns MGS + sample (sorted) ; fingerprints live in _state, so nothing is recomputed
    local f m S n w=0 lst done_ rc
    compgen -G "$D/MGS=*/sample=*/data.parquet" > /dev/null || return 0
    until lock_take "$D/.migrate.lock"; do
        (( w % 20 == 0 )) && say "another job is moving the results of MMT <= 0.2 into the new layout - waiting"
        sleep 15; w=$(( w + 1 ))
        (( w < 480 )) || die "still locked after 2 h ($D/.migrate.lock) - remove it if no MMT job is running"
    done
    n=0; lst=$(tmpname "$D/.migrate.list"); done_=$(tmpname "$D/.migrate.done")
    for f in "$D"/MGS=*/sample=*/data.parquet; do      # one python process for all files (there may be 10^5)
        [[ -f $f ]] || continue
        m=${f#"$D"/MGS=}; m=${m%%/*}; S=${f%/data.parquet}; S=${S##*/sample=}
        if [[ -s $D/$m/$S.parquet ]]; then rm -f "$f"; rmdir "${f%/data.parquet}" "$D/MGS=$m" 2> /dev/null || true
        else printf '%s\t%s\t%s\t%s\n' "$f" "$D/$m/$S.parquet" "$m" "$S"; fi
    done > "$lst"
    rc=0; "$PYTHON" "$MMT_HOME/lib/mmtpq.py" migrate --list "$lst" > "$done_" || rc=$?
    local bad=()
    while IFS=$'\t' read -r f r; do                    # the files written so far, also after a failure
        if [[ $r == ERROR* ]]; then mv -f "$f" "$f.unreadable"; bad+=("$f"); continue; fi   # recomputed (no new file)
        m=${f#"$D"/MGS=}; m=${m%%/*}
        rm -f "$f"; rmdir "${f%/data.parquet}" "$D/MGS=$m" 2> /dev/null || true
        n=$(( n + 1 ))
    done < "$done_"
    (( ${#bad[@]} == 0 )) || say "WARNING: ${#bad[@]} result(s) of MMT <= 0.2 could not be read (renamed *.unreadable, they are recomputed): ${bad[*]:0:3}"
    rm -f "$lst" "$done_"
    lock_release "$D/.migrate.lock"
    (( rc == 0 )) || die "moving the results of MMT <= 0.2 into the new layout failed after $n file(s) - rerun to continue"
    (( n == 0 )) || { say "moved $n result(s) of MMT <= 0.2 into the new layout ($D/<MGS>/<sample>.parquet)"
                      mlog MOVE "moved $n result(s) of MMT <= 0.2 into MGScall/$MODE/<MGS>/<sample>.parquet (not recomputed)"; }
}
