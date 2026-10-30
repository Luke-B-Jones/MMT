# shellcheck shell=bash
# MMT shared helpers, sourced by the `mmt` wrapper (which sets MMT_HOME and MMT_VERSION).
#   messages + run log | configuration | resources | tools | stamps + fingerprints | locks | ID lists | MF4 index

# ============================================================================ messages + run log
# Until start_log everything goes to the terminal. After it, stdout and stderr go to the run log, and only say / die
# also reach the terminal (saved as fd 3 / fd 4) - an sbatch .out/.err file stays a few lines long.
# The run log is kept only with -l/--log DIR (DIR/<analysis>/<name>.log, appended) ; without it it is a private file in
# the work space, removed when the run ends - after a failure its last lines are printed to the terminal first.
MMT_LOGGING=0
MMT_QUIET=0       # 1 = say writes to the log only (commands whose terminal output is data, e.g. --list-samples)
MMT_LOG=""
RUN_ENDED=0       # 1 once a run has written its END line (the EXIT trap then adds no STOP line)
log() { echo "[$(date '+%F %T')] $*"; }
say() { log "$*"; if (( MMT_LOGGING && ! MMT_QUIET )); then echo "$*" >&3; fi; }
die() {
    log "ERROR: $*" >&2
    if (( MMT_LOGGING )); then echo "ERROR: $*${MMT_LOG_DIR:+ (log: $MMT_LOG)}" >&4; fi
    exit 1
}
usage_die() {
    echo "mmt ${MMT_CMD:-}: $*  (see: mmt help ${MMT_CMD:-})" >&2
    if (( MMT_LOGGING )); then echo "mmt ${MMT_CMD:-}: $*  (see: mmt help ${MMT_CMD:-})" >&4; fi
    exit 2
}
start_log() {   # start_log ANALYSIS NAME : run log = <-l DIR>/ANALYSIS/NAME.log, else a private file under the work space
    if [[ -n ${MMT_LOG_DIR:-} ]]; then
        MMT_LOG="$(abspath "$MMT_LOG_DIR")/$1/$2.log"
    else
        MMT_LOG="$TMPROOT/logs/$1.$2.$(hostname -s).$$.log"
    fi
    mkdir -p "$(dirname "$MMT_LOG")" || usage_die "cannot create the log folder $(dirname "$MMT_LOG")"
    exec 3>&1 4>&2 >> "$MMT_LOG" 2>&1
    MMT_LOGGING=1
    log "======== mmt $MMT_CMD - MMT $(mmt_version) on $(hostname -s), job ${JOBTAG:-?}"
}
end_log() {   # end_log RC : without -l the run log is removed (after a failure its last lines go to the terminal first)
    (( MMT_LOGGING )) || return 0
    [[ -n ${MMT_LOG_DIR:-} ]] && return 0
    if (( $1 != 0 )); then { echo "---- last lines of the run log (keep whole logs with -l DIR):"; tail -n 25 "$MMT_LOG"; } >&4; fi
    if (( ${MMT_DEBUG:-0} )); then echo "run log kept (--debug): $MMT_LOG" >&4; else rm -f "$MMT_LOG"; fi
}
log_ref() { if [[ -n ${MMT_LOG_DIR:-} ]]; then printf '%s' "$MMT_LOG"; else printf 'no run log kept (run with -l DIR)'; fi; }

# ============================================================================ OUT/mmt_<date>.log (the output folder's own log)
# One human-readable log per output folder and date, appended by every run of every command: its start and end, its
# plan per MGS and what happened to every sample x MGS. Each line stands alone: "<date> <time>  <command>  <EVENT>  <text>",
# so `grep MGS.12 OUT/mmt_*.log` or `grep PDB15 OUT/mmt_*.log` gives the history of an MGS or a sample. The date in the
# name is today's, or MMT_LOG_DATE when set (a pipeline sets it once, so every job of one run writes to one file).
# Detail (every step, tool versions, timings) is in the run logs (-l DIR).
MLOG_CMD=""       # e.g. "MGScall --mgs-rep" - set by the command
MLOG_ON=0         # set by the command once OUT is known (dry runs leave it off)
MLOG_FILE=""
mlog() {   # mlog EVENT LINE [LINE...] : append to OUT/mmt_<date>.log under a lock (appends from several nodes to one NFS
           # file are not atomic) ; the lock is held for one write: one left by a killed writer is broken after ~60 s
    (( MLOG_ON )) || return 0
    local l="$MMT_OUT/.mmt.log.lock" i=0 ev=$1 h x; shift
    MLOG_FILE=${MLOG_FILE:-$MMT_OUT/mmt_${MMT_LOG_DATE:-$(date +%F)}.log}
    h="$(date '+%F %T')  $MLOG_CMD  $ev"
    until mkdir "$l" 2> /dev/null; do
        i=$(( i + 1 ))
        if (( i > 240 )); then rm -rf "$l"; i=0; fi
        sleep 0.25
    done
    for x in "$@"; do printf '%s  %s\n' "$h" "$x"; done >> "$MLOG_FILE"
    rmdir "$l" 2> /dev/null || true
}
mlog_start() {   # mlog_start "mgs-rep" ARGS... : START line (job, host, version, arguments, run log)
    MLOG_CMD=$1; MLOG_ON=1; shift
    mlog START "job ${JOBTAG:-?} on $(hostname -s) | MMT $(mmt_version) | ${*:-(no selection: all)} | run log: $(log_ref)"
}
last_error() {   # the run log's last ERROR, short: a failed command is cut to its first 80 characters
    local e re="^(line [0-9]+: ')(.*)(' failed \(exit [0-9]+\))$"
    e=$(grep -h 'ERROR: ' "$MMT_LOG" 2> /dev/null | tail -n 1 | sed 's/^\[[^]]*\] ERROR: //' || true)
    if [[ $e =~ $re ]] && (( ${#BASH_REMATCH[2]} > 80 )); then e="${BASH_REMATCH[1]}${BASH_REMATCH[2]:0:80} ...${BASH_REMATCH[3]}"; fi
    printf '%s' "${e:0:400}"
}
id_list() {   # id_list WORD... : "N: a, b, c" (all names - the log is where they are looked up)
    local IFS=,
    if (( $# )); then printf '%d: %s' "$#" "$*" | sed 's/,/, /g'; else printf '0'; fi
}
job_tag() {   # JOBTAG = <array job>_<task> | <job> | local<pid>
    JOBTAG=${SLURM_ARRAY_JOB_ID:+${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}}
    JOBTAG=${JOBTAG:-${SLURM_JOB_ID:-local$$}}
}

run_exit() {   # EXIT trap of a run: an early stop (no END line yet) is written to the mmt log with its reason ; run log
               # removed unless -l
    local rc=$1 e
    if (( rc != 0 && ! RUN_ENDED )); then
        e=$(last_error)
        if (( rc == 143 )); then e="stopped by a signal (scancel / time limit)"; fi
        mlog STOP "exit $rc: ${e:-see the run log} - $(log_ref)"
    fi
    end_log "$rc"
}

plan_log() {   # PLAN lines in the mmt log: per selected MGS where it has a MAG, what is to do, and where GeneCat's abundance
               # table has it without a MAG (needs MGS_SEL OKS UNITS P_TODO P_CUR todo_samples ; reads IDX/abundance.tsv)
    local m s n_in n_todo n_cur
    local -a lines=() nb=()
    declare -A ab=()
    while IFS=$'\t' read -r s m _; do ab[$s|$m]=1; done < <(tail -n +2 "$IDX/abundance.tsv")
    lines+=("selection: $(id_list "${MGS_SEL[@]}") ; usable samples ${#OKS[@]} ; samples with work ${#todo_samples[@]}")
    for m in "${MGS_SEL[@]}"; do
        n_in=0; nb=()
        for s in "${OKS[@]}"; do
            if [[ " ${UNITS[$s]:-} " == *" $m "* ]]; then n_in=$(( n_in + 1 )); elif [[ -n ${ab[$s|$m]:-} ]]; then nb+=("$s"); fi
        done
        n_todo=$(wc -w <<< "${P_TODO[$m]:-}"); n_cur=$(wc -w <<< "${P_CUR[$m]:-}")
        lines+=("$m: a MAG in the assembly of $n_in of ${#OKS[@]} usable samples -> to do $n_todo, up to date $n_cur ; present (GeneCat abundance) but not binned: $(id_list "${nb[@]}")")
    done
    mlog PLAN "${lines[@]}"
}

# any failing command (set -e) is reported with its line, in the log and on the terminal
on_error() { die "line $1: '$2' failed (exit $3)"; }
trap 'on_error "$LINENO" "$BASH_COMMAND" "$?"' ERR
# isolated CMD ARGS... : run CMD in a subshell that stops at its first error ; RC = its exit status ; the caller goes on.
#   (bash ignores set -e inside anything tested by if / || / &&, so the status is returned in RC instead)
isolated() {
    trap - ERR; set +e
    ( set -e; trap 'on_error "$LINENO" "$BASH_COMMAND" "$?"' ERR; "$@" )
    RC=$?
    set -e; trap 'on_error "$LINENO" "$BASH_COMMAND" "$?"' ERR
}

mmt_version() {   # e.g. "0.1.0 (git 1a2b3c4)" ; "+dirty" when the checkout has uncommitted changes
    local c d=""
    c=$(git -C "$MMT_HOME" rev-parse --short HEAD 2> /dev/null) || { echo "$MMT_VERSION (no git)"; return 0; }
    [[ -n $(git -C "$MMT_HOME" status --porcelain --untracked-files=no 2> /dev/null) ]] && d="+dirty"
    echo "$MMT_VERSION (git $c$d)"
}

# ============================================================================ configuration
load_config() {   # load_config "$@" : mmt.conf -> site.conf -> $MMT_CONFIG -> --config FILE (options come later)
    local a prev=""
    # shellcheck source=/dev/null
    source "$MMT_HOME/config/mmt.conf"
    # shellcheck source=/dev/null
    [[ -f $MMT_HOME/config/site.conf ]] && source "$MMT_HOME/config/site.conf"
    if [[ -n ${MMT_CONFIG:-} ]]; then
        [[ -f $MMT_CONFIG ]] || usage_die "MMT_CONFIG=$MMT_CONFIG: file not found"
        # shellcheck source=/dev/null
        source "$MMT_CONFIG"
    fi
    for a in "$@"; do
        if [[ $prev == --config ]]; then
            [[ -f $a ]] || usage_die "--config $a: file not found"
            # shellcheck source=/dev/null
            source "$a"
        fi
        prev=$a
    done
    return 0
}

# common_opt "$@" : one option shared by all commands at $1 ; sets COMMON_SHIFT (how many words it used, 0 = none)
common_opt() {
    COMMON_SHIFT=0
    case $1 in
        --debug|--keep-tmp) MMT_DEBUG=1; COMMON_SHIFT=1; return 0 ;;   # keep work folders (success or failure)
        -i|--mf4-run|-o|--out|--tmp|-t|--threads|--mem|--genecat|--map|--config|-l|--log|--ipd-args)
            [[ $# -ge 2 && -n $2 ]] || usage_die "$1 needs a value" ;;
        *) return 1 ;;
    esac
    case $1 in
        -i|--mf4-run) MMT_MF4_RUN=$2 ;;
        -o|--out)     MMT_OUT=$2 ;;
        --tmp)        MMT_TMP=$2 ;;
        -t|--threads) MMT_THREADS=$2 ;;
        --mem)        MMT_MEM_GB=$2 ;;
        --genecat)    MMT_GENECAT=$2 ;;
        --map)        MMT_MAP=$2 ;;
        -l|--log)     MMT_LOG_DIR=$2 ;;
        --ipd-args)   MMT_IPD_ARGS=$2 ;;
        --config)     : ;;   # already read by load_config
    esac
    COMMON_SHIFT=2
}

# ============================================================================ resources
resolve_resources() {   # THREADS, MEM_GB: option/config, else the SLURM allocation, else 4 cpus / 16 GB
    THREADS=${MMT_THREADS:-${SLURM_CPUS_PER_TASK:-4}}
    if [[ -n ${MMT_MEM_GB:-} ]]; then MEM_GB=$MMT_MEM_GB
    elif [[ -n ${SLURM_MEM_PER_NODE:-} ]]; then MEM_GB=$(( SLURM_MEM_PER_NODE / 1024 ))
    elif [[ -n ${SLURM_MEM_PER_CPU:-} ]]; then MEM_GB=$(( SLURM_MEM_PER_CPU * THREADS / 1024 ))
    else MEM_GB=16; fi
    [[ $THREADS =~ ^[0-9]+$ && $THREADS -ge 1 ]] || usage_die "threads must be a positive number (got '$THREADS')"
    [[ $MEM_GB =~ ^[0-9]+$ ]] || usage_die "memory must be whole GB (got '$MEM_GB')"
    (( MEM_GB >= 2 )) || MEM_GB=2
}

# ============================================================================ tools
need_tool() {   # need_tool VAR COMMAND LABEL : VAR = absolute path of COMMAND, or stop
    local p
    p=$(command -v -- "$2" 2> /dev/null) || die "$3 not found: '$2' - set it in config/site.conf (mmt check)"
    printf -v "$1" '%s' "$p"
}
smrt_tool() {   # smrt_tool VAR NAME : SMRT Link tool from MMT_SMRT_BIN, else PATH
    if [[ -n ${MMT_SMRT_BIN:-} ]]; then
        [[ -x $MMT_SMRT_BIN/$2 ]] || die "SMRT Link tool $2 not found in MMT_SMRT_BIN=$MMT_SMRT_BIN"
        printf -v "$1" '%s' "$MMT_SMRT_BIN/$2"
    else
        need_tool "$1" "$2" "SMRT Link tool $2"
    fi
}
setup_tools() {   # PYTHON SAMTOOLS DUCKDB IPDSUMMARY BYSTRANDIFY PBMM2 PBINDEX, with the features MMT relies on
    local h
    if [[ -n ${MMT_ENV_PREFIX:-} ]]; then
        [[ -d $MMT_ENV_PREFIX/bin ]] || die "MMT_ENV_PREFIX=$MMT_ENV_PREFIX has no bin/ folder"
        PATH="$MMT_ENV_PREFIX/bin:$PATH"
    fi
    need_tool PYTHON "$MMT_PYTHON" python
    need_tool SAMTOOLS "$MMT_SAMTOOLS" samtools
    need_tool DUCKDB "$MMT_DUCKDB" "DuckDB command-line binary"
    smrt_tool IPDSUMMARY ipdSummary
    smrt_tool BYSTRANDIFY ccs-kinetics-bystrandify
    smrt_tool PBMM2 pbmm2
    smrt_tool PBINDEX pbindex
    h=$("$SAMTOOLS" mpileup 2>&1 || true)
    [[ $h == *--reverse-del* ]] || die "samtools $SAMTOOLS has no 'mpileup --reverse-del' - needs >= 1.15"
    "$PYTHON" -c 'import pysam, pyarrow, numpy' 2> /dev/null || die "python $PYTHON cannot import pysam, pyarrow and numpy"
    read -ra IPD_EXTRA <<< "${MMT_IPD_ARGS:-}"   # handed to ipdSummary as they are (e.g. --ipdModel FILE)
    # CRAM: never let htslib download a reference (compute nodes have no internet ; MMT reads CRAMs without one)
    export REF_PATH=/nonexistent/mmt REF_CACHE=/nonexistent/mmt
}
tool_versions() {   # one line per tool (for logs and `mmt check`)
    local v
    v=$("$SAMTOOLS" --version 2>&1 || true); echo "samtools   ${v%%$'\n'*}   $SAMTOOLS"
    v=$("$PYTHON" -c 'import sys, pysam, pyarrow; print("python", sys.version.split()[0], "pysam", pysam.__version__, "pyarrow", pyarrow.__version__)' 2>&1 || true)
    echo "python     $v   $PYTHON"
    v=$("$DUCKDB" --version 2>&1 || true); echo "duckdb     $v   $DUCKDB"
    v=$("$PBMM2" --version 2>&1 || true); echo "pbmm2      ${v%%$'\n'*}   $PBMM2"
    v=$("$IPDSUMMARY" --version 2>&1 || true); echo "ipdSummary ${v%%$'\n'*}   $IPDSUMMARY"
    echo "bystrandify $BYSTRANDIFY"
    echo "pbindex    $PBINDEX"
}

# ============================================================================ stamps, fingerprints, files
stamp() { stat -L -c '%s:%Y' -- "$1" 2> /dev/null || echo missing; }   # size:mtime of a file
md5str() { printf '%s' "$1" | md5sum | cut -c1-16; }
abspath() { readlink -m -- "$1"; }
publish() { mv -f -- "$1" "$2"; }   # publish TMP FINAL : rename into place (same file system = all or nothing)
tmpname() { echo "$(dirname "$1")/.$(basename "$1").$(hostname -s).$$.tmp"; }   # private temp name next to FILE ($$: the job's
                                                                                 # shell, alive while the name is in use)
work_root() {   # WROOT = work folders of this OUT: <tmp>/mmt/<OUT folder name>_<key> (several OUTs may share one --tmp)
    local b; b=$(basename "$MMT_OUT"); b=${b//[^A-Za-z0-9._-]/_}
    WROOT="$TMPROOT/${b}_$(md5str "$MMT_OUT" | cut -c1-8)"
}
# work_done DIR OK : remove a work folder (OK = 1 success, 0 failure) unless --debug ; says where a kept one is
work_done() {
    local d=$1 ok=$2
    [[ -n $d && -d $d ]] || return 0
    if (( ${MMT_DEBUG:-0} )); then log "work folder kept (--debug): $d"; return 0; fi
    rm -rf -- "$d"
    (( ok )) || log "work folder removed after the failure (rerun with --debug to keep it): $d"
    return 0   # (shared parent folders stay: another job may be creating its folder in them right now)
}

read_ids() {   # read_ids "A,B,C" | FILE : one ID per line (a file: first word of each non-empty, non-# line)
    if [[ -f $1 ]]; then
        awk '!/^[[:space:]]*(#|$)/ {print $1}' "$1"
    else
        tr -s ', \t' '\n' <<< "$1" | awk 'NF {print $1}'
    fi
}

# ============================================================================ locks
# mkdir is atomic on every file system an HPC uses (flock is not reliable on network file systems).
# A lock left by a job that died is taken over: the owner's SLURM job is no longer in squeue, or (no SLURM)
# the owner's process is gone on this host.
lock_take() {   # lock_take DIR : 0 = we hold it, 1 = someone else does
    local d=$1 host="" pid="" job=""
    if mkdir "$d" 2> /dev/null; then
        echo "$(hostname) $BASHPID ${SLURM_JOB_ID:-none}" > "$d/owner"
        return 0
    fi
    # judging + replacing a stale lock: one process at a time (a takeover mutex left > 5 min = by a killed process)
    mkdir "$d.takeover" 2> /dev/null || { [[ -n $(find "$d.takeover" -maxdepth 0 -mmin +5 2> /dev/null) ]] && rm -rf "$d.takeover"; return 1; }
    read -r host pid job < "$d/owner" 2> /dev/null || true
    if lock_stale "$d" "$host" "$pid" "$job"; then
        log "taking over stale lock $d (owner: ${host:-?} pid ${pid:-?} job ${job:-?})"
        rm -rf "$d"
        if mkdir "$d" 2> /dev/null; then
            echo "$(hostname) $BASHPID ${SLURM_JOB_ID:-none}" > "$d/owner"
            rmdir "$d.takeover"; return 0
        fi
    fi
    rmdir "$d.takeover" 2> /dev/null; return 1
}
lock_stale() {   # lock_stale DIR HOST PID JOB
    local d=$1 host=$2 pid=$3 job=$4
    if [[ -z $host ]]; then   # owner line not written yet: only stale when the lock is older than a minute
        [[ -n $(find "$d" -maxdepth 0 -mmin +1 2> /dev/null) ]]; return
    fi
    if [[ -n $job && $job != none ]] && command -v squeue > /dev/null; then
        local o
        o=$(squeue -h -j "$job" -o %i 2>&1) && { [[ -z $o ]]; return; }   # squeue answered: stale when the job is not listed
        [[ $o == *"Invalid job id"* ]]; return                             # any other error (controller busy) = not stale
    fi
    [[ $host == "$(hostname)" ]] && ! kill -0 "$pid" 2> /dev/null
}
lock_release() { rm -rf -- "$1"; }

# ============================================================================ MF4 index (shared by all modules)
INDEX_METHOD=4   # bump when lib/mf4.py index output changes -> every OUT rebuilds its index on the next run
INDEX_FILES="samples.tsv mgs.tsv mags.tsv units.tsv assemblies.tsv contigs.tsv.gz rep_contigs.tsv.gz movies.tsv abundance.tsv inputs.tsv report.txt"
resolve_inputs() {   # absolute MMT_MAP, MMT_MF4_RUN (from the map when not given), MMT_GENECAT, MVG (MAGvsGC)
    local out rid
    [[ -n ${MMT_MAP:-} ]] || usage_die "no MATAFILER4 map: --map FILE or MMT_MAP in config/site.conf"
    [[ -f $MMT_MAP ]] || die "map not found: $MMT_MAP"
    MMT_MAP=$(abspath "$MMT_MAP")
    if [[ -z ${MMT_MF4_RUN:-} ]]; then
        out=$(awk '{sub(/\r$/, "")} $1 == "#OutPath" {print $2; exit}' "$MMT_MAP")
        rid=$(awk '{sub(/\r$/, "")} $1 == "#RunID" {print $2; exit}' "$MMT_MAP")
        [[ -n $out && -n $rid ]] || usage_die "no MATAFILER4 run folder: -i DIR (the map has no #OutPath + #RunID)"
        MMT_MF4_RUN="${out%/}/$rid"
    fi
    [[ -d $MMT_MF4_RUN ]] || die "MATAFILER4 run folder not found: $MMT_MF4_RUN"
    MMT_MF4_RUN=$(abspath "$MMT_MF4_RUN")
    [[ -n ${MMT_GENECAT:-} ]] || usage_die "no GeneCat folder: --genecat DIR or MMT_GENECAT in config/site.conf"
    MMT_GENECAT=$(abspath "$MMT_GENECAT")
    MVG="$MMT_GENECAT/Bin_$MMT_BINNER/MAGvsGC.txt.gz"
    [[ -s $MVG ]] || die "GeneCat MGS table not found: $MVG (GeneCat MGS stage finished? MMT_BINNER=$MMT_BINNER?)"
    ABUND="$MMT_GENECAT/Bin_$MMT_BINNER/Annotation/Abundance/MGS.matL0.txt"   # GeneCat's MGS abundance (presence)
}
resolve_out() {   # absolute MMT_OUT (created) ; TMPROOT = <tmp>/mmt ; WROOT = this OUT's work space there
    [[ -n ${MMT_OUT:-} ]] || usage_die "no output folder: -o DIR"
    mkdir -p "$MMT_OUT" || die "cannot create output folder $MMT_OUT"
    MMT_OUT=$(abspath "$MMT_OUT")
    TMPROOT=$(abspath "${MMT_TMP:-${TMPDIR:-/tmp}}")/mmt
    mkdir -p "$TMPROOT" || die "cannot create work space $TMPROOT"
    work_root
}

# ensure_index : <work space>/index/ (a cache: rebuilt whenever missing) describes these inputs - tables rebuilt when the
#   map, MAGvsGC, GeneCat's abundance table or settings changed: $INDEX_FILES + index.fp. Results do not depend on
#   index.fp itself but on the checksums in units.tsv / assemblies.tsv, so an edit of the map that changes nothing for a
#   sample redoes nothing. Several jobs may start at once: one builds (lock), the others wait for it.
index_current() {   # the index in IDX was built from these inputs and is complete
    local f
    [[ $(cat "$IDX/index.fp" 2> /dev/null) == "$INDEX_FP" ]] || return 1
    for f in $INDEX_FILES; do [[ -s $IDX/$f ]] || return 1; done
    # the assembly.txt and bin files it read are unchanged (re-binning changes MAG contigs without the map / MAGvsGC)
    cmp -s <(cut -f2 "$IDX/inputs.tsv") <(cut -f1 "$IDX/inputs.tsv" | while read -r f; do stamp "$f"; done)
}
ensure_index() {
    local w=0
    IDX="$WROOT/index"
    INDEX_FP=$(md5str "$INDEX_METHOD|$MMT_MF4_RUN|$MMT_GENECAT|$MMT_BINNER|$MMT_MAP|$(stamp "$MMT_MAP")|$(stamp "$MVG")|$(stamp "$ABUND")")
    index_current && return 0
    mkdir -p "$IDX"
    until lock_take "$IDX/.lock"; do
        (( w % 20 == 0 )) && log "index: another job is building it - waiting"
        sleep 15; w=$(( w + 1 ))
        index_current && return 0
        (( w < 240 )) || die "index: still locked after 1 h ($IDX/.lock) - remove it if no MMT job is running"
    done
    if index_current; then lock_release "$IDX/.lock"; return 0; fi
    log "index: building from $MMT_MF4_RUN + $MMT_GENECAT (Bin_$MMT_BINNER) + $MMT_MAP"
    rm -rf "$IDX/.build"; mkdir -p "$IDX/.build"
    if ! "$PYTHON" "$MMT_HOME/lib/mf4.py" index --mf4-run "$MMT_MF4_RUN" --genecat "$MMT_GENECAT" --map "$MMT_MAP" \
            --binner "$MMT_BINNER" --abundance "$ABUND" --out "$IDX/.build"; then
        [[ -f $IDX/.build/report.txt ]] && cp "$IDX/.build/report.txt" "$IDX/report.failed.txt"
        lock_release "$IDX/.lock"
        die "index: could not be built - see the log ($IDX/report.failed.txt)"
    fi
    local f
    for f in $INDEX_FILES; do publish "$IDX/.build/$f" "$IDX/$f"; done
    rm -f "$IDX/report.failed.txt"
    echo "$INDEX_FP" > "$IDX/.index.fp" && publish "$IDX/.index.fp" "$IDX/index.fp"
    rm -rf "$IDX/.build"
    lock_release "$IDX/.lock"
    log "index: $(grep -c . "$IDX/samples.tsv") samples (header incl.), $(grep -c . "$IDX/mgs.tsv") MGS (header incl.) - $IDX/report.txt"
    local -a il=()
    mapfile -t il < <(grep -E '^(samples in map|raw BAM read groups|MGS in MAGvsGC|member MAGs|MGS):' "$IDX/report.txt")
    mapfile -t -O "${#il[@]}" il < <(grep -E '^GeneCat abundance' "$IDX/report.txt")
    mlog INDEX "index (re)built from $MMT_MAP + $MVG - details: $IDX/report.txt" "${il[@]}"
    grep -qE '^member MAGs: ok [0-9]+$' "$IDX/report.txt" || say "WARNING: index: $(grep -m1 '^member MAGs:' "$IDX/report.txt") - see $IDX/report.txt"
}
