# shellcheck shell=bash
# MMT shared helpers, sourced by the `mmt` wrapper (which sets MMT_HOME and MMT_VERSION).
#   messages + run log | configuration | resources | tools | stamps + fingerprints | locks | ID lists | MF4 index

# ============================================================================ messages + run log
# Until start_log everything goes to the terminal. After it, stdout and stderr go to the run log, and only say / die
# also reach the terminal (saved as fd 3 / fd 4) - an sbatch .out/.err file stays a few lines long.
MMT_LOGGING=0
MMT_QUIET=0       # 1 = say writes to the log only (commands whose terminal output is data, e.g. --list-samples)
MMT_LOG=""
log() { echo "[$(date '+%F %T')] $*"; }
say() { log "$*"; if (( MMT_LOGGING && ! MMT_QUIET )); then echo "$*" >&3; fi; }
die() {
    log "ERROR: $*" >&2
    if (( MMT_LOGGING )); then echo "ERROR: $* (log: $MMT_LOG)" >&4; fi
    exit 1
}
usage_die() {
    echo "mmt ${MMT_CMD:-}: $*  (see: mmt help ${MMT_CMD:-})" >&2
    if (( MMT_LOGGING )); then echo "mmt ${MMT_CMD:-}: $*  (see: mmt help ${MMT_CMD:-})" >&4; fi
    exit 2
}
start_log() {
    MMT_LOG=$1
    mkdir -p "$(dirname "$MMT_LOG")"
    exec 3>&1 4>&2 >> "$MMT_LOG" 2>&1
    MMT_LOGGING=1
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
        -i|--mf4-run|-o|--out|--tmp|-t|--threads|--mem|--genecat|--map|--config)
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
    "$PYTHON" -c 'import pysam' 2> /dev/null || die "python $PYTHON cannot import pysam"
    # CRAM: never let htslib download a reference (compute nodes have no internet ; MMT reads CRAMs without one)
    export REF_PATH=/nonexistent/mmt REF_CACHE=/nonexistent/mmt
}
tool_versions() {   # one line per tool (for logs and `mmt check`)
    local v
    v=$("$SAMTOOLS" --version 2>&1 || true); echo "samtools   ${v%%$'\n'*}   $SAMTOOLS"
    v=$("$PYTHON" -c 'import sys, pysam; print("python", sys.version.split()[0], "pysam", pysam.__version__)' 2>&1 || true)
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
tmpname() { echo "$(dirname "$1")/.$(basename "$1").$(hostname -s).$BASHPID.tmp"; }   # private temp name next to FILE

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
INDEX_METHOD=2   # bump when lib/mf4.py index output changes -> every OUT rebuilds its index on the next run
INDEX_FILES="samples.tsv mgs.tsv mags.tsv units.tsv assemblies.tsv contigs.tsv.gz rep_contigs.tsv.gz inputs.tsv report.txt"
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
}
resolve_out() {   # absolute MMT_OUT (created)
    [[ -n ${MMT_OUT:-} ]] || usage_die "no output folder: -o DIR"
    mkdir -p "$MMT_OUT" || die "cannot create output folder $MMT_OUT"
    MMT_OUT=$(abspath "$MMT_OUT")
}

# ensure_index : OUT/index/ describes these inputs (tables rebuilt when the map, MAGvsGC or settings changed):
#   $INDEX_FILES + index.fp. Results do not depend on index.fp itself but on the checksums in units.tsv /
#   assemblies.tsv, so an edit of the map that changes nothing for a sample redoes nothing.
#   Several jobs may start at once: one builds (lock), the others wait for it.
index_current() {   # the index in IDX was built from these inputs and is complete
    local f
    [[ $(cat "$IDX/index.fp" 2> /dev/null) == "$INDEX_FP" ]] || return 1
    for f in $INDEX_FILES; do [[ -s $IDX/$f ]] || return 1; done
    # the assembly.txt and bin files it read are unchanged (re-binning changes MAG contigs without the map / MAGvsGC)
    cmp -s <(cut -f2 "$IDX/inputs.tsv") <(cut -f1 "$IDX/inputs.tsv" | while read -r f; do stamp "$f"; done)
}
ensure_index() {
    local w=0
    IDX="$MMT_OUT/index"
    INDEX_FP=$(md5str "$INDEX_METHOD|$MMT_MF4_RUN|$MMT_GENECAT|$MMT_BINNER|$MMT_MAP|$(stamp "$MMT_MAP")|$(stamp "$MVG")")
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
            --binner "$MMT_BINNER" --out "$IDX/.build"; then
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
    grep -qE '^member MAGs: ok [0-9]+$' "$IDX/report.txt" || say "WARNING: index: $(grep -m1 '^member MAGs:' "$IDX/report.txt") - see $IDX/report.txt"
}
