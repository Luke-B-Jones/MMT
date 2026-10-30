# shellcheck shell=bash
# MMT derep-ctg (alias: dereplicate-contigs, derep) - each MGS as a short list of collinear blocks ("dereplicated
# contigs"), every occurrence in every MAG down to the base (lift). Sourced by the mmt wrapper. Method, outputs and
# columns: docs/derep-ctg.md. Engine: lib/derep.py.

DEREP_METHOD=3        # bump when a code change alters results -> the next run redoes every MGS
DEREP_PREP_METHOD=2   # bump when lib/derep.py prep output changes -> prep rebuilt

module_help() {
    cat << 'EOF'
mmt derep-ctg - each MGS as a short list of collinear blocks ("dereplicated contigs")

Usage: mmt derep-ctg -o OUT [options]          (also accepted: mmt dereplicate-contigs, mmt derep)

For every MGS (GeneCat species bin): the genes of all its MAGs (MATAFILER4 prodigal genes, GeneCat 95 % catalogue)
-> gene families (catalogue genes that are variants of one gene merged by context + protein identity) -> blocks =
runs of families in the same order and orientation in the MAGs (core blocks: the persistent backbone; accessory
blocks: recurrent variable regions) -> every occurrence of every block in every MAG (instances), a representative
sequence per block (polished to the majority base of its instances), how the blocks follow each other along the
contigs (links), and every base of every instance on its block (lift: input of `mmt ctg-meth`). No reads are used.

Inputs (defaults from config/site.conf)
  -i, --mf4-run DIR   MATAFILER4 run folder (default: the map's #OutPath/#RunID)
  --genecat DIR       GeneCat folder (Bin_<binner>/MAGvsGC.txt.gz + compl.incompl.<id>.fna.clstr.idx + .prot.faa)
  --map FILE          MATAFILER4 map
Selection
  --mgs IDS           only these MGS: "MGS.12,MGS.40" or a file (one ID per line); default every MGS with a usable MAG
                      (--select is the same)
  --part I/N          the I-th of N parts of the selection (balanced by MAG count) - for SLURM arrays
Other
  --prep-only         build the shared preparation (catalogue + genes of all MGS MAGs) and stop - run once, then arrays
  --list-mgs          print the MGS that have work (one per line, nothing else) and stop
  --dry-run           show what would run and stop
  --force             redo the selected MGS even if up to date
  --debug             keep the per-MGS work files (instance alignments) in the work folder
  -t, --threads N / --mem GB / --tmp DIR / -o, --out DIR / --config FILE   as for every command
Every run is summarised in OUT/mmt.log (start, MGS done / failed, end).
Settings: DEREP_* in config/mmt.conf (documented there and in docs/derep-ctg.md); a changed setting redoes the MGS
it affects on the next run.

Output (in OUT/derep-ctg)
  <MGS>/blocks.fna.gz        one sequence per block (the dereplicated contigs)
  <MGS>/blocks.parquet       one row per block (+ its sequence)   instances.parquet  every block occurrence in every MAG
  <MGS>/lift.parquet         every base of every instance on its block (exact-match segments)
  <MGS>/block_genes.parquet  genes of each block                  links.parquet      block ends that meet on contigs
  <MGS>/genes.parquet        every gene of every MAG              families.parquet   gene families and their class
  <MGS>/contigs.parquet      block path of every contig           mags.parquet       block coverage of every MAG
  <MGS>/composition.parquet  base counts of every MAG contig in 1 kb windows (GC skew)
  <MGS>/merged_catalogue_genes.parquet, summary.tsv ; summary.tsv (all MGS) ; _prep/ (shared preparation)
A results folder of MMT <= 0.2 (OUT/dereplicate-contigs/) is moved to OUT/derep-ctg/ by the first run (its _prep kept).

Examples
  mmt derep-ctg -o OUT --mgs MGS.2345,MGS.453        # two MGS, quick
  mmt derep-ctg -o OUT                                # every MGS, one job
  SLURM: prepare once, then an array over parts, then one plain run (refreshes summary.tsv):
    sbatch -c 16 --mem=16G -t 2:00:00 --wrap 'mmt derep-ctg -o OUT --prep-only -t 16'
    sbatch --array=1-10 -c 4 --mem=8G -t 4:00:00 --wrap 'mmt derep-ctg -o OUT --part ${SLURM_ARRAY_TASK_ID}/10 --tmp $TMPDIR'
EOF
}

derep_tools() {   # PYTHON (pysam, pyarrow, edlib when merging), SAMTOOLS, MINIMAP2 (when DEREP_IDENTITY=1)
    if [[ -n ${MMT_ENV_PREFIX:-} ]]; then
        [[ -d $MMT_ENV_PREFIX/bin ]] || die "MMT_ENV_PREFIX=$MMT_ENV_PREFIX has no bin/ folder"
        PATH="$MMT_ENV_PREFIX/bin:$PATH"
    fi
    need_tool PYTHON "$MMT_PYTHON" python
    need_tool SAMTOOLS "$MMT_SAMTOOLS" samtools
    local mods="pysam, pyarrow"
    [[ $DEREP_FAM_MIN_ID != 0 ]] && mods+=", edlib"
    "$PYTHON" -c "import sys, $mods; assert sys.version_info >= (3, 10)" 2> /dev/null \
        || die "python $PYTHON: derep-ctg needs python >= 3.10 with $mods"
    MINIMAP2=""
    if [[ $DEREP_IDENTITY == 1 ]]; then
        need_tool MINIMAP2 "${MMT_MINIMAP2:-minimap2}" "minimap2 (set MMT_MINIMAP2, or DEREP_IDENTITY=0 to skip instance identity)"
    fi
}

derep_settings() {   # DEREP_JSON (for lib/derep.py), DEREP_TXT (logs, fingerprints) ; checks every value
    local k v json="" txt=""
    local -a keys=(fam_min_id fam_len_ratio dup_margin mixed_mag max_copy hub_partners min_family_mags min_edge_mags synteny
                   max_insert_bp max_skip core_prev min_block_genes min_block_bp identity)
    for k in "${keys[@]}"; do
        v=DEREP_${k^^}; v=${!v:-}
        [[ $v =~ ^(0|[1-9][0-9]*)(\.[0-9]+)?$ ]] || usage_die "DEREP_${k^^}='$v' is not a number"
        case $k in fam_min_id|fam_len_ratio|dup_margin|mixed_mag|synteny|core_prev) ;;
            *) [[ $v =~ ^[0-9]+$ ]] || usage_die "DEREP_${k^^}='$v' must be a whole number" ;; esac
        json+="${json:+,}\"$k\":$v"; txt+="${txt:+ }$k=$v"
    done
    awk -v s="$DEREP_SYNTENY" 'BEGIN {exit !(s > 0.5 && s <= 1)}' || usage_die "DEREP_SYNTENY must be > 0.5 and <= 1 (got $DEREP_SYNTENY)"
    awk -v s="$DEREP_FAM_MIN_ID" -v r="$DEREP_FAM_LEN_RATIO" -v d="$DEREP_DUP_MARGIN" -v x="$DEREP_MIXED_MAG" -v c="$DEREP_CORE_PREV" \
        'BEGIN {exit !(s <= 1 && r <= 1 && d <= 1 && x <= 1 && c <= 1)}' \
        || usage_die "DEREP_FAM_MIN_ID, _FAM_LEN_RATIO, _DUP_MARGIN, _MIXED_MAG and _CORE_PREV must be within 0..1"
    (( DEREP_MIN_EDGE_MAGS >= 1 && DEREP_MIN_FAMILY_MAGS >= 1 && DEREP_MAX_COPY >= 2 )) \
        || usage_die "DEREP_MIN_EDGE_MAGS and DEREP_MIN_FAMILY_MAGS must be >= 1, DEREP_MAX_COPY >= 2"
    DEREP_JSON="{$json}"; DEREP_TXT="$txt"
}

prep_fp() {   # PREP = D/_prep ; PREP_FP = what prep reads: MAG contigs + their assemblies (index), the assemblies'
              # contig + gene files, the catalogue (a touched map / MAGvsGC, a .fai next to an assembly: nothing to redo)
    local cidx="" f asm ctg
    PREP="$D/_prep"
    for f in "$MMT_GENECAT"/compl.incompl."$DEREP_GC_ID".fna.{clstr,ctsv}.idx{,.gz}; do [[ -s $f ]] && { cidx=$f; break; }; done
    [[ -n $cidx ]] || die "no compl.incompl.$DEREP_GC_ID.fna.clstr.idx in $MMT_GENECAT (DEREP_GC_ID?)"
    asm=$(awk -F'\t' 'NR > 1 {print $1}' "$IDX/assemblies.tsv" | while read -r a; do
        stat -L -c '%n:%s:%Y' "$a"/scaffolds.fasta.filt{,.gz} "$a"/genePred/genes{,.bac}.gff{,.gz} 2> /dev/null || true; done | md5sum | cut -c1-16)
    ctg=$( { gzip -dc "$IDX/contigs.tsv.gz" | cut -f1-4; cut -f1,5 "$IDX/samples.tsv"; } | md5sum | cut -c1-16)
    PREP_FP=$(md5str "$DEREP_PREP_METHOD|$ctg|$asm|$cidx:$(stamp "$cidx")|$(stamp "$MMT_GENECAT/compl.incompl.$DEREP_GC_ID.prot.faa")$(stamp "$MMT_GENECAT/compl.incompl.$DEREP_GC_ID.prot.faa.gz")|check=$DEREP_PREP_CHECK")
}

ensure_prep() {   # build PREP once for PREP_FP (lock: one job builds, the others wait)
    local w=0 f
    prep_current() { [[ $(cat "$PREP/prep.fp" 2> /dev/null) == "$PREP_FP" && -s $PREP/assemblies.tsv && -s $PREP/catalogue.tsv ]]; }
    prep_current && return 0
    until lock_take "$D/.prep.lock"; do
        (( w % 20 == 0 )) && say "prep: another job is building it - waiting"
        sleep 15; w=$(( w + 1 ))
        prep_current && return 0
        (( w < 960 )) || die "prep: still locked after 4 h ($D/.prep.lock) - remove it if no MMT job is running"
    done
    if prep_current; then lock_release "$D/.prep.lock"; return 0; fi
    say "prep: catalogue + genes of every MGS MAG ($THREADS threads) - once per index"
    rm -rf "$D/.prep.build"
    if ! "$PYTHON" "$MMT_HOME/lib/derep.py" prep --idx "$IDX" --genecat "$MMT_GENECAT" --gc-id "$DEREP_GC_ID" \
            --out "$D/.prep.build" --threads "$THREADS" --check "$DEREP_PREP_CHECK" --samtools "$SAMTOOLS"; then
        f=$(tail -n 1 "$MMT_LOG" | cut -c1-400)
        rm -rf "$D/.prep.build"; lock_release "$D/.prep.lock"; die "prep failed: ${f#ERROR: }"
    fi
    echo "$PREP_FP" > "$D/.prep.build/prep.fp"
    rm -rf "$D/.prep.old"; [[ -d $PREP ]] && mv "$PREP" "$D/.prep.old"
    mv "$D/.prep.build" "$PREP"; rm -rf "$D/.prep.old"
    lock_release "$D/.prep.lock"
}

module_main() {
    local IDS="" PART="" FORCE=0 DRY=0 LIST=0 PREP_ONLY=0 t0 m st fp n n_ok=0 row=0 sel desc err
    MMT_DEBUG=${MMT_DEBUG:-0}
    local -a todo=() failed=()
    while (( $# )); do
        if common_opt "$@"; then shift "$COMMON_SHIFT"; continue; fi
        case $1 in
            --mgs|--select) [[ $# -ge 2 && -n $2 ]] || usage_die "$1 needs MGS IDs or a file"; IDS=$2; shift ;;
            --part)   [[ $# -ge 2 && $2 =~ ^[0-9]+/[0-9]+$ ]] || usage_die "--part needs I/N, e.g. 3/20"; PART=$2; shift ;;
            --force)  FORCE=1 ;;
            --dry-run) DRY=1 ;;
            --list-mgs) LIST=1 ;;
            --prep-only) PREP_ONLY=1 ;;
            -h|--help) module_help; exit 0 ;;
            *) usage_die "unknown option '$1'" ;;
        esac
        shift
    done
    (( LIST )) && MMT_QUIET=1
    t0=$(date +%s)
    derep_settings
    resolve_out
    job_tag
    start_log derep-ctg "$(date +%Y%m%d-%H%M%S)_$JOBTAG"
    say "MMT $(mmt_version) derep-ctg -> $MMT_OUT/derep-ctg${MMT_LOG_DIR:+ (log: $MMT_LOG)}"
    log "command: mmt derep-ctg ${IDS:+--mgs $IDS }${PART:+--part $PART }(force=$FORCE dry=$DRY list=$LIST prep_only=$PREP_ONLY debug=$MMT_DEBUG) on $(hostname), job $JOBTAG"
    desc="${IDS:+--mgs $IDS }${PART:+--part $PART }"
    (( FORCE )) && desc+="--force "
    (( LIST )) && desc+="(listing the MGS with work) "
    (( PREP_ONLY )) && desc+="--prep-only "
    (( DRY )) || mlog_start "derep-ctg" "${desc% }"
    trap 'run_exit $?' EXIT
    trap 'exit 143' TERM INT
    resolve_inputs
    resolve_resources
    derep_tools
    log "settings: $DEREP_TXT | gc_id=$DEREP_GC_ID check=$DEREP_PREP_CHECK | method $DEREP_METHOD | $THREADS threads"
    log "tools: python $PYTHON | samtools $SAMTOOLS | minimap2 ${MINIMAP2:-(not used)}"
    ensure_index
    D="$MMT_OUT/derep-ctg"; migrate_derep; mkdir -p "$D"
    prep_fp
    if (( PREP_ONLY )); then
        ensure_prep; say "prep ready: $PREP ($(( $(grep -c . "$PREP/assemblies.tsv") - 1 )) assemblies)"
        mlog END "prep ready ($(( $(grep -c . "$PREP/assemblies.tsv") - 1 )) assemblies) ; $(( ($(date +%s) - t0) / 60 )) min"
        RUN_ENDED=1; return 0
    fi
    BASE_FP=$(md5str "$DEREP_METHOD|$PREP_FP|$DEREP_TXT")
    sel=""
    if [[ -n $IDS ]]; then
        sel="$D/.sel.$BASHPID"; read_ids "$IDS" > "$sel"
        [[ -s $sel ]] || usage_die "--mgs: no MGS IDs in '$IDS'"
    fi
    local -a pa=(--idx "$IDX" --out "$D" --base-fp "$BASE_FP")
    [[ -n $sel ]] && pa+=(--mgs "$sel")
    [[ -n $PART ]] && pa+=(--part "$PART")
    (( FORCE )) && pa+=(--force)
    local plan perr; perr=$(tmpname "$D/plan.err")
    if ! plan=$("$PYTHON" "$MMT_HOME/lib/derep.py" plan "${pa[@]}" 2> "$perr"); then
        cat "$perr"; m=$(tail -n 1 "$perr"); rm -f "$perr" "$sel"; die "plan: ${m#ERROR: }"
    fi
    cat "$perr"; rm -f "$perr"
    [[ -n $sel ]] && rm -f "$sel"
    n=0
    while IFS=$'\t' read -r m st fp _; do
        [[ -n $m ]] || continue
        if [[ $st == stale ]]; then   # (full runs only) the MGS is no longer in GeneCat's table / has no usable MAG
            (( DRY )) || { rm -rf "${D:?}/$m"; say "dropped $m (no longer an MGS with a usable MAG)"; }
            continue
        fi
        n=$(( n + 1 ))
        [[ $st == todo ]] && todo+=("$m|$fp")
    done <<< "$plan"
    if (( n == 0 )) && [[ -n $PART ]]; then say "part $PART holds no MGS - nothing to do"; mlog END "part $PART holds no MGS"; RUN_ENDED=1; return 0; fi
    (( n )) || die "no MGS selected (none with a usable MAG?)"
    say "plan: ${#todo[@]} MGS to do, $(( n - ${#todo[@]} )) up to date (selection: $n)"
    mlog PLAN "${#todo[@]} MGS to do$( (( ${#todo[@]} )) && echo " ($(id_list "${todo[@]%%|*}"))"), $(( n - ${#todo[@]} )) up to date (selection: $n)"
    if (( LIST )); then (( ${#todo[@]} == 0 )) || printf '%s\n' "${todo[@]%%|*}" >&3; mlog END "listing only"; RUN_ENDED=1; return 0; fi
    if (( DRY )); then printf '%s\n' "$plan" | awk -F'\t' 'NF {printf "  %s %s (%s MAGs)\n", $1, $2, $4}' | while read -r l; do say "$l"; done
        say "dry run - nothing done"; RUN_ENDED=1; return 0; fi
    (( ${#todo[@]} )) && ensure_prep
    for x in "${todo[@]}"; do
        m=${x%%|*}; fp=${x#*|}
        isolated run_mgs "$m" "$fp"
        if (( RC == 0 )); then n_ok=$(( n_ok + 1 )); row=0; mlog MGS "$m $(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i} NR == 2 {printf "%s: %s MAGs, %s core + %s accessory blocks, %s instances", $h["status"], $h["n_mags"], $h["n_blocks_core"], $h["n_blocks_accessory"], $h["n_instances"]}' "$D/$m/summary.tsv" 2> /dev/null)"
        else failed+=("$m"); row=$(( row + 1 )); log "== $m FAILED (exit $RC)"; err=$(last_error); mlog FAIL "$m: ${err:-exit $RC} - $(log_ref)"; fi
        (( row < 5 )) || { say "5 MGS failed in a row - stopping (a problem common to all MGS? see the log)"; break; }
    done
    "$PYTHON" "$MMT_HOME/lib/derep.py" summary --out "$D" > /dev/null || log "WARNING: summary.tsv not rebuilt"
    mlog END "$n_ok MGS done, ${#failed[@]} failed$( (( ${#failed[@]} )) && echo " ($(id_list "${failed[@]}"))") ; $(( ($(date +%s) - t0) / 60 )) min"
    RUN_ENDED=1
    if (( ${#failed[@]} )); then
        say "finished with errors after $(( ($(date +%s) - t0) / 60 )) min: ${#failed[@]} MGS failed: ${failed[*]:0:20} - rerun the same command to resume"
        exit 1
    fi
    say "done in $(( ($(date +%s) - t0) / 60 )) min: $n_ok MGS processed ; table of all MGS: $D/summary.tsv"
}

migrate_derep() {   # MMT <= 0.2: OUT/dereplicate-contigs/ with MGS=<MGS>/ folders -> OUT/derep-ctg/ (its _prep kept ; the
                    # MGS folders are dropped: DEREP_METHOD 3 redoes them anyway)
    local old="$MMT_OUT/dereplicate-contigs" w=0 d
    [[ -d $old && ! -e $D ]] || { rm -rf "$D"/MGS=* 2> /dev/null || true; return 0; }
    until lock_take "$MMT_OUT/.derep-move.lock"; do
        sleep 15; w=$(( w + 1 )); (( w < 240 )) || die "still locked after 1 h ($MMT_OUT/.derep-move.lock)"
    done
    if [[ -d $old && ! -e $D ]]; then
        if compgen -G "$old/.lock.*" > /dev/null || [[ -e $old/.prep.lock ]]; then
            lock_release "$MMT_OUT/.derep-move.lock"; die "$old has locks of running jobs - wait for them, then rerun"
        fi
        mv "$old" "$D"; say "moved $old -> $D (results of MMT <= 0.2; their MGS are redone)"
    fi
    lock_release "$MMT_OUT/.derep-move.lock"
    for d in "$D"/MGS=*; do [[ -d $d ]] && rm -rf "$d"; done
    return 0
}

run_mgs() {   # run_mgs MGS FINGERPRINT : one MGS, under a lock (two jobs never write the same MGS)
    local m=$1 fp=$2 L="$D/.lock.$1" w=0
    until lock_take "$L"; do        # another job is on it (overlapping runs): wait, then skip it if that made it current
        if (( w == 0 )); then say "$m: another job is working on it ($(cat "$L/owner" 2> /dev/null)) - waiting"
        elif (( w % 20 == 0 )); then log "$m: still waiting for the other job"; fi
        (( w < 240 )) || die "$m: still locked after 1 h ($L) - remove it if no MMT job is running"
        sleep 15; w=$(( w + 1 ))
    done
    trap 'lock_release "'"$L"'"' EXIT
    # the plan is from this run's start: another job (an array task, an early refresh run) may have done it meanwhile
    if (( ! FORCE )) && [[ -s $D/$m/summary.tsv && $(cat "$D/$m/.fp" 2> /dev/null) == "$fp" ]]; then
        log "== $m: done meanwhile by another job"; return 0
    fi
    log "== $m"
    rm -rf "$D/.$m".*.tmp "$D/.$m".*.old "$WROOT"/derep_"$m"_*          # left by a killed job (same OUT: the MGS lock is held)
    "$PYTHON" "$MMT_HOME/lib/derep.py" mgs --idx "$IDX" --prep "$PREP" --mgs "$m" --out "$D" --tmp "$WROOT" \
        --threads "$THREADS" --minimap2 "$MINIMAP2" --settings "$DEREP_JSON" --settings-txt "$DEREP_TXT" --fp "$fp" \
        --mmt-version "$(mmt_version)" --prep-fp "$PREP_FP" $( (( MMT_DEBUG )) && echo --keep-tmp) \
        || { local rc=$? e; e=$(tail -n 1 "$MMT_LOG" | cut -c1-300); e=${e#ERROR: }; die "$m (exit $rc): ${e#"$m: "}"; }
}
