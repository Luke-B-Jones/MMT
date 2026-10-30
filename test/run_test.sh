#!/usr/bin/env bash
# run_test.sh WORKDIR - end-to-end test of MMT MGScall on a synthetic MATAFILER4 run (test/mock_mf4.py). WORKDIR is wiped.
#   Needs samtools + minimap2 on PATH and:
#     MMT_TEST_PY          python with pysam (default python3)
#     MMT_TEST_IPD         a kineticsTools ipdSummary (e.g. pip install kineticsTools) ; MMT_TEST_IPD_MODEL its model
#                          (passed with MMT_IPD_ARGS="--ipdModel ...", as a user would)
#     MMT_TEST_DUCKDB      DuckDB command-line binary (default duckdb)
#   SMRT Link's bystrandify / pbmm2 / pbindex are replaced by stand-ins (test/standins.sh). Takes ~30 min on 2 cpus.
# shellcheck disable=SC2016,SC2034,SC2012  # checks are strings run by eval in t()
set -uo pipefail
W=$(readlink -m "$1"); HERE=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd); MMT="$HERE/../mmt"
PY=${MMT_TEST_PY:-python3}; IPD=${MMT_TEST_IPD:?set MMT_TEST_IPD}; DUCK=${MMT_TEST_DUCKDB:-duckdb}
MODEL=${MMT_TEST_IPD_MODEL:?set MMT_TEST_IPD_MODEL}
DUCK=$(command -v "$DUCK") || { echo "no duckdb"; exit 2; }
rm -rf "$W"; mkdir -p "$W"
NP=0 NF=0
IDXD="$W/tmp/mmt/out_$(printf '%s' "$W/out" | md5sum | cut -c1-8)/index"   # the index cache of OUT (MMT work_root)
t() { if eval "$2"; then echo "PASS  $1"; NP=$((NP + 1)); else echo "FAIL  $1"; NF=$((NF + 1)); fi; }
OUT="$W/out" LOGS="$W/logs"
run() {   # run ARGS... -> RC, TERM_TXT (terminal output), LOG (the run log written last)
    "$MMT" "$@" --config "$W/test.conf" -l "$LOGS" > "$W/term.txt" 2>&1; RC=$?
    TERM_TXT=$(cat "$W/term.txt"); LOG=$(ls -t "$LOGS"/*/*.log 2> /dev/null | head -1)
    echo "    \$ mmt $*  -> exit $RC, $(grep -c . "$W/term.txt") terminal line(s)"
}
f() { if [[ $1 == mgs-rep ]]; then echo "$OUT/mgs-rep/$2/$3.parquet"; else echo "$OUT/$1/$2/$3/$4.parquet"; fi; }
mt() { stat -c %Y "$(f "$@")" 2> /dev/null || echo none; }
ino() { stat -c %i "$(f "$@")" 2> /dev/null || echo none; }
q() { "$DUCK" -noheader -list -c "$1"; }
ML() { cat "$OUT"/mmt_*.log; }
G2A=S2.run2M2__SemiBin_3 G2B=S2.run2M2__SemiBin_5.fa.gz

echo "== mock MATAFILER4 run + stand-ins"
"$PY" "$HERE/mock_mf4.py" "$W/mock" > /dev/null || { echo "mock failed"; exit 1; }
bash "$HERE/standins.sh" "$W/bin" "$(command -v "$PY")" "$IPD"
cat > "$W/test.conf" << EOF
MMT_MAP="$W/mock/map.txt"
MMT_GENECAT="$W/mock/gc"
MMT_SMRT_BIN="$W/bin"
MMT_PYTHON="$(command -v "$PY")"
MMT_DUCKDB="$DUCK"
MMT_TMP="$W/tmp"
MMT_THREADS=2
MMT_MEM_GB=4
MMT_IPD_ARGS="--ipdModel $MODEL"
EOF

echo "== check, index (in the work space), concurrent index builds"
run check -o "$OUT"
t "check: exit 0, 2 usable MGS, GeneCat abundance read" '[[ $RC == 0 && $TERM_TXT == *"MGS: ok 2  rep_bin_not_in_bin_file 1"* && $TERM_TXT == *"GeneCat abundance: "*"3 MGS present in 4 samples"* ]]'
t "index in the work space, not in OUT" '[[ -s $IDXD/units.tsv && -s $IDXD/abundance.tsv && ! -e $OUT/index ]]'
( "$MMT" check -o "$W/out2" --config "$W/test.conf" > "$W/c1.txt" 2>&1; echo $? > "$W/c1.rc" ) &
( "$MMT" check -o "$W/out2" --config "$W/test.conf" > "$W/c2.txt" 2>&1; echo $? > "$W/c2.rc" ) &
wait
t "two index builds at once: both succeed, same index, no lock left" \
  '[[ $(cat "$W/c1.rc")$(cat "$W/c2.rc") == 00 && -z $(ls -d "$W"/tmp/mmt/out2_*/index/.lock 2> /dev/null) ]] && cmp -s "$(ls "$W"/tmp/mmt/out2_*/index/units.tsv)" "$IDXD/units.tsv"'

echo "== selection, listing, option errors"
run MGScall --mag-rep -o "$OUT" --list-samples
t "mag-rep --list-samples prints every sample with a MAG: S1 S2 S2.run2 S3" '[[ $RC == 0 && $(echo $TERM_TXT) == "S1 S2 S2.run2 S3" ]]'
run MGScall --mgs-rep -o "$OUT" --list-samples
t "mgs-rep --list-samples prints S1 S2 S2.run2 S3" '[[ $(echo $TERM_TXT) == "S1 S2 S2.run2 S3" ]]'
t "PLAN: MGS.1 present in S3 by GeneCat abundance but not binned" 'ML | grep -q "mgs-rep  PLAN  MGS.1: a MAG in the assembly of 3 of 4 usable samples -> to do 3, up to date 0 ; present (GeneCat abundance) but not binned: 1: S3"'
run MGScall --mgs-rep -o "$OUT" --select MGS.9
t "--select unknown MGS: exit 1 with a message" '[[ $RC == 1 && $TERM_TXT == *"not in MAGvsGC: MGS.9"* ]]'
run MGScall --mgs-rep -o "$OUT" --select MGS.3
t "--select unusable MGS only: warning, exit 1" '[[ $RC == 1 && $TERM_TXT == *"not usable: MGS.3"* && $TERM_TXT == *"no usable MGS"* ]]'
run MGScall -o "$OUT"
t "no mode: usage error 2" '[[ $RC == 2 ]]'
run MGScall --mgs-rep --mag-rep -o "$OUT"
t "two modes: usage error 2" '[[ $RC == 2 ]]'
run MGScall --mag-rep -o "$OUT" --samples S1,S7
t "--samples unknown: exit 1" '[[ $RC == 1 && $TERM_TXT == *"not in the map: S7"* ]]'

echo "== mgs-rep: one MGS with an injected failure, resume, then all MGS"
MMT_TEST_FAIL_IPD=S2 run MGScall --mgs-rep -o "$OUT" --select MGS.1
t "failure in S2: exit 1, others done, error on the terminal" \
  '[[ $RC == 1 && $TERM_TXT == *"ERROR"* && $TERM_TXT == *"failed: S2 "* && $(mt mgs-rep MGS.1 S1) != none && $(mt mgs-rep MGS.1 S2.run2) != none && $(mt mgs-rep MGS.1 S2) == none ]]'
t "failed sample left no lock, its work folder removed (no --debug)" \
  '[[ ! -e $OUT/mgs-rep/stats/S2.lock && -z $(ls -d "$W"/tmp/mmt/out_*/mgs-rep/S2 2> /dev/null) ]] && grep -q "work folder removed after the failure" "$LOG"'
t "ipdSummary failure: named in the run log" 'grep -q "ERROR: ipdSummary failed" "$LOG"'
t "mmt log: the failure of S2 with its reason and log, the run end with the failed sample" \
  'ML | grep -q "mgs-rep  FAIL  S2: ipdSummary failed.* - $LOGS/mgs-rep/" && ML | grep -q "mgs-rep  END  2 sample(s) processed, 1 failed (1: S2)"'
MMT_TEST_FAIL_IPD=S2 run MGScall --mgs-rep -o "$OUT" --select MGS.1 --samples S2 --debug
t "failure with --debug: work folder kept (steps up to the kinetics done), one-sample run log named after it" \
  '[[ $RC == 1 && -f $(ls -d "$W"/tmp/mmt/out_*/mgs-rep/S2 2> /dev/null | head -1)/align.ok && $LOG == "$LOGS/mgs-rep/S2.log" ]] && grep -q "work folder kept (--debug)" "$LOG"'
run MGScall --mgs-rep -o "$OUT" --select MGS.1
t "resume: only S2, from its finished steps" \
  '[[ $RC == 0 && $TERM_TXT == *"plan: 1 sample(s) with work"* && $(mt mgs-rep MGS.1 S2) != none ]] && ! grep -q "S2: .*extracted" "$LOG" && grep -q "S2: kinetics table ready" "$LOG"'
m1=$(mt mgs-rep MGS.1 S1)$(mt mgs-rep MGS.1 S2)$(mt mgs-rep MGS.1 S2.run2)
sleep 1; run MGScall --mgs-rep -o "$OUT"
t "all MGS: 4 samples, 4 results added" '[[ $RC == 0 && $TERM_TXT == *"plan: 4 sample(s) with work (4 MGS x sample results to make)"* ]]'
t "MGS.1 results of the --select run kept as they were" '[[ $m1 == "$(mt mgs-rep MGS.1 S1)$(mt mgs-rep MGS.1 S2)$(mt mgs-rep MGS.1 S2.run2)" ]]'
t "kinetics model recorded per sample (the --ipd-args model)" 'grep -q "^kinetics_model	.*SP3-C3.npz.gz (passed-in)" "$OUT/mgs-rep/stats/S1.run.tsv" && grep -q "^ipdsummary_args	--ipdModel" "$OUT/mgs-rep/stats/S1.run.tsv" && grep -q "^binding_kit	" "$OUT/mgs-rep/stats/S1.run.tsv"'
t "stats: abundance from GeneCat, MAG list, representative" '[[ $(awk -F"\t" "\$1 == \"MGS.0\" {print \$2, \$3, \$5, \$10}" "$OUT/mgs-rep/stats/S1.tsv") == "S1__SemiBin_1,S1__SemiBin_6 done 20.1 S1__SemiBin_1" ]]'

echo "== mag-rep: each MAG apart ; one-MAG samples linked to the mgs-rep table, S1 MGS.0 (2 MAGs) computed"
s0=$(date +%s); run MGScall --mag-rep -o "$OUT"
t "mag-rep: exit 0, 3 terminal lines" '[[ $RC == 0 && $(grep -c . "$W/term.txt") == 3 ]]'
t "one-MAG units: hard links of the mgs-rep files (same inode), status linked" \
  '[[ $(ino mag-rep MGS.1 S1 S1__SemiBin_4) == $(ino mgs-rep MGS.1 S1) && $(ino mag-rep MGS.0 S3 S3__SemiBin_2) == $(ino mgs-rep MGS.0 S3) && $(ino mag-rep MGS.1 S2 $G2B) == $(ino mgs-rep MGS.1 S2) ]] && [[ $(awk -F"\t" "NR > 1 && \$3 == \"linked\"" "$OUT"/mag-rep/stats/S2.tsv | wc -l) == 2 ]]'
t "S1 MGS.0: one table per MAG, computed (only S1 extracted)" \
  '[[ -s $(f mag-rep MGS.0 S1 S1__SemiBin_1) && -s $(f mag-rep MGS.0 S1 S1__SemiBin_6) && $(ino mag-rep MGS.0 S1 S1__SemiBin_1) != $(ino mgs-rep MGS.0 S1) ]] && [[ $(grep -c "extracted" "$LOG") == 1 ]] && grep -q "S1: .* extracted and split into 2 MAG" "$LOG" && grep -q "MAG layer 2" "$LOG"'
s0=$(date +%s); run MGScall --mag-rep -o "$OUT"
t "mag-rep rerun: nothing to do, < 30 s" '[[ $RC == 0 && $TERM_TXT == *"plan: 0 sample(s) with work"* && $(( $(date +%s) - s0 )) -lt 30 ]]'

echo "== mag-mag: each MAG's reads on the MAG itself"
run MGScall --mag-mag -o "$OUT"
t "mag-mag: exit 0, a table per MAG (the units of mag-rep), all computed (nothing linked)" \
  '[[ $RC == 0 && $(cd "$OUT/mag-mag" && ls -d MGS*/*/*.parquet | tr "\n" " ") == "$(cd "$OUT/mag-rep" && ls -d MGS*/*/*.parquet | tr "\n" " ")" ]] && [[ $(ino mag-mag MGS.0 S3 S3__SemiBin_2) != $(ino mag-rep MGS.0 S3 S3__SemiBin_2) && -z $(awk -F"\t" "FNR > 1 && \$3 != \"done\"" "$OUT"/mag-mag/stats/S?.tsv "$OUT"/mag-mag/stats/S2.run2.tsv) ]]'
t "mag-mag: MAG references cut from the assemblies, cached in the work space (6 MAGs: S2 and S2.run2 share theirs)" '[[ $(ls "$W"/tmp/mmt/refs_mag | wc -l) == 6 ]] && grep -q "MAG layer 2" "$LOG"'
s0=$(date +%s); run MGScall --mag-mag -o "$OUT"
t "mag-mag rerun: nothing to do, < 30 s" '[[ $RC == 0 && $TERM_TXT == *"plan: 0 sample(s) with work"* && $(( $(date +%s) - s0 )) -lt 30 ]]'

echo "== results against the truth"
DUCKDB_BIN="$DUCK" "$PY" "$HERE/check_truth.py" "$W/mock" "$OUT" "$IDXD" mgs-rep mag-rep mag-mag > "$W/truth.txt"; RC=$?
sed 's/^/    /' "$W/truth.txt"
t "truth checks" '[[ $RC == 0 ]]'

echo "== settings change, map edit, MGS membership change"
echo 'MGSCALL_MIN_ROW_COV=12' > "$W/cov12.conf"
a1=$(mt mag-rep MGS.0 S1 S1__SemiBin_1); sleep 1
run MGScall --mag-rep -o "$OUT" --samples S2 --config "$W/cov12.conf"
t "new setting: S2 recomputed (the mgs-rep table no longer matches: no link), coverage >= 12, S1 untouched" \
  '[[ $RC == 0 && $(q "SELECT min(coverage) FROM read_parquet('"'"'$(f mag-rep MGS.1 S2 $G2B)'"'"')") -ge 12 && $(ino mag-rep MGS.1 S2 $G2B) != $(ino mgs-rep MGS.1 S2) && $(mt mag-rep MGS.0 S1 S1__SemiBin_1) == "$a1" ]]'
run MGScall --mag-rep -o "$OUT"
t "back to the default setting: S2 linked again" '[[ $RC == 0 && $TERM_TXT == *"plan: 1 sample(s) with work"* && $(ino mag-rep MGS.1 S2 $G2B) == $(ino mgs-rep MGS.1 S2) ]]'
echo "# a comment added to the map" >> "$W/mock/map.txt"
run MGScall --mgs-rep -o "$OUT" --dry-run
t "map edited: index rebuilt, nothing to redo" '[[ $TERM_TXT == *"plan: 0 sample(s) with work"* ]] && grep -q "index: building" "$LOG"'
zcat "$W/mock/gc/Bin_SB/MAGvsGC.txt.gz" | grep -v '^S1__SemiBin_4' | gzip > "$W/mvg.gz" && mv "$W/mvg.gz" "$W/mock/gc/Bin_SB/MAGvsGC.txt.gz"
b1=$(mt mgs-rep MGS.0 S1); sleep 1
run MGScall --mgs-rep -o "$OUT" --dry-run
t "S1 lost its MGS.1 MAG: plan drops that result only" '[[ $TERM_TXT == *"plan S1: 0 to do, 1 up to date, to drop: MGS.1"* && $TERM_TXT == *"plan: 1 sample(s) with work (0 MGS"* ]]'
run MGScall --mgs-rep -o "$OUT"
run MGScall --mag-rep -o "$OUT"
run MGScall --mag-mag -o "$OUT"
t "S1 MGS.1 tables and state lines removed in every mode, S1 MGS.0 untouched" \
  '[[ $RC == 0 && $(mt mgs-rep MGS.1 S1) == none && ! -e $OUT/mag-rep/MGS.1/S1 && ! -e $OUT/mag-mag/MGS.1/S1 && $(mt mgs-rep MGS.0 S1) == "$b1" ]] && ! grep -q "^MGS.1" "$OUT/mgs-rep/stats/S1.tsv" "$OUT/mag-rep/stats/S1.tsv" "$OUT/mag-mag/stats/S1.tsv"'

echo "== locks"
mkdir -p "$OUT/mag-rep/stats/S1.lock"; echo "$(hostname) $$ none" > "$OUT/mag-rep/stats/S1.lock/owner"
run MGScall --mag-rep -o "$OUT" --samples S1 --force
t "sample locked by a live process: exit 1, lock kept" \
  '[[ $RC == 1 && $TERM_TXT == *"another job is working on this sample"* && -d $OUT/mag-rep/stats/S1.lock ]]'
rm -rf "$OUT/mag-rep/stats/S1.lock"
mkdir -p "$OUT/mgs-rep/stats/S3.lock"; echo "$(hostname) 4194000 none" > "$OUT/mgs-rep/stats/S3.lock/owner"
run MGScall --mgs-rep -o "$OUT" --samples S3 --force
t "lock of a dead process: taken over, sample redone, lock released" \
  '[[ $RC == 0 && ! -e $OUT/mgs-rep/stats/S3.lock ]] && grep -q "taking over stale lock" "$LOG"'

echo "== logs: without -l no run log is kept ; a failure then prints the log's last lines"
n0=$(find "$W/tmp" -name "*.log" | wc -l)
"$MMT" MGScall --mgs-rep -o "$OUT" --samples S3 --force --config "$W/test.conf" > "$W/term.txt" 2>&1; RC=$?
t "no -l: success, nothing left in the work space or OUT" '[[ $RC == 0 && $(find "$W/tmp" -name "*.log" | wc -l) == "$n0" ]] && ML | grep -q "run log: no run log kept"'
MMT_TEST_FAIL_IPD=S3 "$MMT" MGScall --mgs-rep -o "$OUT" --samples S3 --force --config "$W/test.conf" > "$W/term.txt" 2>&1; RC=$?
t "no -l, failure: the log's last lines on the terminal, log removed" '[[ $RC == 1 ]] && grep -q "last lines of the run log" "$W/term.txt" && grep -q "injected failure" "$W/term.txt" && [[ $(find "$W/tmp" -name "*.log" | wc -l) == "$n0" ]]'
run MGScall --mgs-rep -o "$OUT" --samples S3

echo "== outputs"
t "OUT holds only the analysis folders and the dated mmt log" '[[ $(ls "$OUT" | tr "\n" " ") == "mag-mag mag-rep mgs-rep mmt_$(date +%F).log " ]]'
t "each analysis folder: MGS folders + stats only" '[[ $(ls "$OUT/mgs-rep" | tr "\n" " ") == "MGS.0 MGS.1 stats " && $(ls "$OUT/mag-rep" | tr "\n" " ") == "MGS.0 MGS.1 stats " && $(ls "$OUT/mag-mag" | tr "\n" " ") == "MGS.0 MGS.1 stats " ]]'
t "no temp files left in OUT" '[[ -z $(find "$OUT" -name ".*.tmp" -o -name "*.lock") ]]'
t "work folders removed after success" '[[ -z $(find "$W/tmp/mmt" -mindepth 3 -maxdepth 3 -path "*/out_*/m??-rep/*") ]]'
t "run logs: one per sample run, named after it" '[[ -s $LOGS/mgs-rep/S2.log && -s $LOGS/mag-rep/S2.log ]] && grep -q "ipdSummary: kinetics model" "$LOGS/mgs-rep/S3.log"'

echo "== the mmt log: one line per event, every run"
t "every non-dry run has a START, every finished run an END" \
  '[[ $(ML | grep -c "  START  ") -ge 10 && $(ML | grep -c "  START  ") == $(( $(ML | grep -c "  END  ") + $(ML | grep -c "  STOP  ") )) ]]'
t "one SAMPLE line per sample processed, with every MGS (mag-rep: every MAG) outcome" \
  'ML | grep -q "mgs-rep  SAMPLE  S3: MGS.0 done ([0-9]* sites, [0-9]* molecules)" && ML | grep -q "mag-rep  SAMPLE  S1: MGS.0 S1__SemiBin_1 done ([0-9]* sites, [0-9]* molecules), S1__SemiBin_6 done" && ML | grep -q "mag-rep  SAMPLE  S3: MGS.0 S3__SemiBin_2 linked"'
t "an index build and the dropped result are logged" 'ML | grep -q "INDEX  index (re)built" && ML | grep -q "SAMPLE  S1: MGS.1 dropped"'
t "dry runs are not logged" '! ML | grep -q "dry"'
t "MMT_LOG_DATE names the log file" 'MMT_LOG_DATE=run42 "$MMT" status -o "$OUT" --config "$W/test.conf" > /dev/null 2>&1 && grep -q "status  STATUS" "$OUT/mmt_run42.log" && rm "$OUT/mmt_run42.log"'

echo "== status"
run status -o "$OUT"
t "status: per MGS and mode, planned samples and their state with names, linked, MAG tables, not binned, models" \
  '[[ $RC == 0 && $TERM_TXT == *"mgs-rep | MGS.0: planned for 4 sample(s)"* && $TERM_TXT == *"mgs-rep | MGS.0 done 4: S1, S2, S2.run2, S3"* && $TERM_TXT == *"mag-rep | MGS.1 linked 2: S2, S2.run2"* && $TERM_TXT == *"mag-rep | MGS.0: planned for 4 sample(s) "*"5 MAG table(s)"* && $TERM_TXT == *"mag-mag | MGS.0 done 4: S1, S2, S2.run2, S3"* && $TERM_TXT == *"mgs-rep | MGS.1 not binned 2: S1, S3"* && $TERM_TXT == *"kinetics model(s): "*"SP3-C3.npz.gz (passed-in)"* ]] && ML | grep -q "status  STATUS  mag-rep | MGS.0"'

echo
echo "$NP passed, $NF failed"
(( NF == 0 ))
