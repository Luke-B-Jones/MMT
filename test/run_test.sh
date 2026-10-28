#!/usr/bin/env bash
# run_test.sh WORKDIR - end-to-end test of MMT on a synthetic MATAFILER4 run (test/mock_mf4.py). WORKDIR is wiped.
#   Needs samtools + minimap2 on PATH and:
#     MMT_TEST_PY          python with pysam (default python3)
#     MMT_TEST_IPD         a kineticsTools ipdSummary (e.g. pip install kineticsTools) ; MMT_TEST_IPD_MODEL its model
#     MMT_TEST_DUCKDB      DuckDB command-line binary (default duckdb)
#   SMRT Link's bystrandify / pbmm2 / pbindex are replaced by stand-ins (test/standins.sh). Takes ~30 min on 2 cpus.
# shellcheck disable=SC2016,SC2034,SC2012  # checks are strings run by eval in t()
set -uo pipefail
W=$(readlink -m "$1"); HERE=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd); MMT="$HERE/../mmt"
PY=${MMT_TEST_PY:-python3}; IPD=${MMT_TEST_IPD:?set MMT_TEST_IPD}; DUCK=${MMT_TEST_DUCKDB:-duckdb}
DUCK=$(command -v "$DUCK") || { echo "no duckdb"; exit 2; }
rm -rf "$W"; mkdir -p "$W"
NP=0 NF=0
t() { if eval "$2"; then echo "PASS  $1"; NP=$((NP + 1)); else echo "FAIL  $1"; NF=$((NF + 1)); fi; }
OUT="$W/out"
run() {   # run ARGS... -> RC, TERM (terminal output), LOG (that run's log)
    local before; before=$(ls "$OUT"/logs/MGScall/*/*.log 2> /dev/null | sort)
    "$MMT" "$@" --config "$W/test.conf" > "$W/term.txt" 2>&1; RC=$?
    TERM_TXT=$(cat "$W/term.txt"); LOG=$(comm -13 <(echo "$before") <(ls "$OUT"/logs/MGScall/*/*.log 2> /dev/null | sort) | tail -1)
    echo "    \$ mmt $*  -> exit $RC, $(grep -c . "$W/term.txt") terminal line(s)"
}
mt() { stat -c %Y "$OUT/MGScall/$1/$2/$3.parquet" 2> /dev/null || echo none; }

echo "== mock MATAFILER4 run + stand-ins"
"$PY" "$HERE/mock_mf4.py" "$W/mock" > /dev/null || { echo "mock failed"; exit 1; }
bash "$HERE/standins.sh" "$W/bin" "$(command -v "$PY")" "$IPD" "${MMT_TEST_IPD_MODEL:-}"
cat > "$W/test.conf" << EOF
MMT_MAP="$W/mock/map.txt"
MMT_GENECAT="$W/mock/gc"
MMT_SMRT_BIN="$W/bin"
MMT_PYTHON="$(command -v "$PY")"
MMT_DUCKDB="$DUCK"
MMT_TMP="$W/tmp"
MMT_THREADS=2
MMT_MEM_GB=4
EOF

echo "== check, index, concurrent index builds"
run check -o "$OUT"
t "check: exit 0, 2 usable MGS" '[[ $RC == 0 && $TERM_TXT == *"MGS: ok 2  rep_bin_not_in_bin_file 1"* ]]'
( "$MMT" check -o "$W/out2" --config "$W/test.conf" > "$W/c1.txt" 2>&1; echo $? > "$W/c1.rc" ) &
( "$MMT" check -o "$W/out2" --config "$W/test.conf" > "$W/c2.txt" 2>&1; echo $? > "$W/c2.rc" ) &
wait
t "two index builds at once: both succeed, same index, no lock left" \
  '[[ $(cat "$W/c1.rc")$(cat "$W/c2.rc") == 00 && ! -e $W/out2/index/.lock ]] && cmp -s "$W/out2/index/units.tsv" "$OUT/index/units.tsv"'

echo "== selection, listing, option errors"
run MGScall --mag-rep -o "$OUT" --list-samples
t "mag-rep --list-samples prints only S1 S2 S2.run2" '[[ $RC == 0 && $(echo $TERM_TXT) == "S1 S2 S2.run2" ]]'
run MGScall --mgs-rep -o "$OUT" --list-samples
t "mgs-rep --list-samples prints S1 S2 S2.run2 S3" '[[ $(echo $TERM_TXT) == "S1 S2 S2.run2 S3" ]]'
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

echo "== mag-rep: full run, then nothing to do"
run MGScall --mag-rep -o "$OUT"
t "mag-rep: exit 0, 3 terminal lines" '[[ $RC == 0 && $(grep -c . "$W/term.txt") == 3 ]]'
s0=$(date +%s); run MGScall --mag-rep -o "$OUT"
t "mag-rep rerun: nothing to do, < 30 s" '[[ $RC == 0 && $TERM_TXT == *"plan: 0 sample(s) with work"* && $(( $(date +%s) - s0 )) -lt 30 ]]'

echo "== results of MMT <= 0.2 (MGS=<m>/sample=<S>/data.parquet) move into the new layout without recomputing"
o="$OUT/MGScall/mag-rep/MGS=MGS.1/sample=S2"; mkdir -p "$o"; n0=$("$DUCK" -noheader -list -c "SELECT count(*) FROM '$OUT/MGScall/mag-rep/MGS.1/S2.parquet'")
"$DUCK" -c "COPY (SELECT * EXCLUDE (MGS, sample) FROM '$OUT/MGScall/mag-rep/MGS.1/S2.parquet' ORDER BY random()) TO '$o/data.parquet' (FORMAT parquet)"
rm "$OUT/MGScall/mag-rep/MGS.1/S2.parquet"
run MGScall --mag-rep -o "$OUT"
t "old layout: moved, sorted, MGS + sample columns, same rows, nothing recomputed" \
  '[[ $RC == 0 && $TERM_TXT == *"moved 1 result(s)"* && $TERM_TXT == *"plan: 0 sample(s) with work"* && ! -e $OUT/MGScall/mag-rep/MGS=MGS.1 ]] && [[ $("$DUCK" -noheader -list -c "SELECT count(*) || '"'"' '"'"'|| min(MGS) || '"'"' '"'"' || min(sample) || '"'"' '"'"' || bool_and(pos >= coalesce(pp, 0) OR contig <> pc) FROM (SELECT *, lag(pos) OVER () pp, lag(contig) OVER () pc FROM '"'"'$OUT/MGScall/mag-rep/MGS.1/S2.parquet'"'"')") == "$n0 MGS.1 S2 true" ]]'
o="$OUT/MGScall/mag-rep/MGS=MGS.1/sample=S2.run2"; mkdir -p "$o"; echo "not a parquet file" > "$o/data.parquet"
rm "$OUT/MGScall/mag-rep/MGS.1/S2.run2.parquet"
run MGScall --mag-rep -o "$OUT"
t "old layout, a damaged file: set aside (*.unreadable), warned, that result recomputed" \
  '[[ $RC == 0 && $TERM_TXT == *"could not be read"* && -f $o/data.parquet.unreadable && -s $OUT/MGScall/mag-rep/MGS.1/S2.run2.parquet && $TERM_TXT == *"plan: 1 sample(s) with work"* ]]'
rm -rf "$OUT/MGScall/mag-rep/MGS=MGS.1"

echo "== mgs-rep: one MGS with an injected failure, resume, then all MGS"
MMT_TEST_FAIL_IPD=S2 run MGScall --mgs-rep -o "$OUT" --select MGS.1
t "failure in S2: exit 1, others done, error on the terminal" \
  '[[ $RC == 1 && $TERM_TXT == *"ERROR"* && $TERM_TXT == *"failed: S2 "* && $(mt mgs-rep MGS.1 S1) != none && $(mt mgs-rep MGS.1 S2.run2) != none && $(mt mgs-rep MGS.1 S2) == none ]]'
t "failed sample left no lock, its work folder removed (no --debug)" \
  '[[ ! -e $OUT/MGScall/mgs-rep/_state/S2.lock && -z $(ls -d "$W"/tmp/mmt/*_*/MGScall_mgs-rep/S2 2> /dev/null) ]] && grep -q "work folder removed after the failure" "$LOG"'
t "mmt.log: the failure of S2 with its reason, the run end with the failed sample" \
  'grep -q "MGScall --mgs-rep  FAIL  S2: line [0-9]*: .*failed (exit 3) - details: logs/MGScall/mgs-rep/" "$OUT/mmt.log" && grep -q "MGScall --mgs-rep  END  2 sample(s) processed, 1 failed (1: S2)" "$OUT/mmt.log"'
MMT_TEST_FAIL_IPD=S2 run MGScall --mgs-rep -o "$OUT" --select MGS.1 --samples S2 --debug
t "failure with --debug: work folder kept (steps up to the kinetics done)" \
  '[[ $RC == 1 && -f $(ls -d "$W"/tmp/mmt/*_*/MGScall_mgs-rep/S2 2> /dev/null | head -1)/align.ok ]] && grep -q "work folder kept (--debug)" "$LOG"'
run MGScall --mgs-rep -o "$OUT" --select MGS.1
t "resume: only S2, from its finished steps" \
  '[[ $RC == 0 && $TERM_TXT == *"plan: 1 sample(s) with work"* && $(mt mgs-rep MGS.1 S2) != none ]] && ! grep -q "S2: .*extracted" "$LOG" && grep -q "S2: kinetics table ready" "$LOG"'
m1=$(mt mgs-rep MGS.1 S1)$(mt mgs-rep MGS.1 S2)$(mt mgs-rep MGS.1 S2.run2)
sleep 1; run MGScall --mgs-rep -o "$OUT"
t "all MGS: 4 samples, 4 results added" '[[ $RC == 0 && $TERM_TXT == *"plan: 4 sample(s) with work (4 MGS x sample results to make)"* ]]'
t "MGS.1 results of the --select run kept as they were" '[[ $m1 == "$(mt mgs-rep MGS.1 S1)$(mt mgs-rep MGS.1 S2)$(mt mgs-rep MGS.1 S2.run2)" ]]'

echo "== results against the truth"
DUCKDB_BIN="$DUCK" "$PY" "$HERE/check_truth.py" "$W/mock" "$OUT" mag-rep mgs-rep > "$W/truth.txt"; RC=$?
sed 's/^/    /' "$W/truth.txt"
t "truth checks" '[[ $RC == 0 ]]'

echo "== settings change, map edit, MGS membership change"
echo 'MGSCALL_MIN_ROW_COV=12' > "$W/cov12.conf"
a1=$(mt mag-rep MGS.0 S1); sleep 1
run MGScall --mag-rep -o "$OUT" --samples S2 --config "$W/cov12.conf"
t "new setting: S2 redone with coverage >= 12, S1 untouched" \
  '[[ $RC == 0 && $("$DUCK" -noheader -list -c "SELECT min(coverage) FROM read_parquet('"'"'$OUT/MGScall/mag-rep/MGS.1/S2.parquet'"'"')") -ge 12 && $(mt mag-rep MGS.0 S1) == "$a1" ]]'
run MGScall --mag-rep -o "$OUT" --dry-run
t "back to the default setting: S2 is due again" '[[ $TERM_TXT == *"plan S2: 1 to do"* && $TERM_TXT == *"plan S1: 0 to do, 1 up to date"* ]]'
echo "# a comment added to the map" >> "$W/mock/map.txt"
run MGScall --mgs-rep -o "$OUT" --dry-run
t "map edited: index rebuilt, nothing to redo" '[[ $TERM_TXT == *"plan: 0 sample(s) with work"* ]] && grep -q "index: building" "$LOG"'
zcat "$W/mock/gc/Bin_SB/MAGvsGC.txt.gz" | grep -v '^S1__SemiBin_4' | gzip > "$W/mvg.gz" && mv "$W/mvg.gz" "$W/mock/gc/Bin_SB/MAGvsGC.txt.gz"
b1=$(mt mgs-rep MGS.0 S1); sleep 1
run MGScall --mgs-rep -o "$OUT" --dry-run
t "S1 lost its MGS.1 MAG: plan drops that result only" '[[ $TERM_TXT == *"plan S1: 0 to do, 1 up to date, to drop: MGS.1"* && $TERM_TXT == *"plan: 1 sample(s) with work (0 MGS"* ]]'
run MGScall --mgs-rep -o "$OUT"
t "S1 MGS.1 result and state line removed, S1 MGS.0 untouched" \
  '[[ $RC == 0 && $(mt mgs-rep MGS.1 S1) == none && $(mt mgs-rep MGS.0 S1) == "$b1" ]] && ! grep -q "^MGS.1" "$OUT/MGScall/mgs-rep/_state/S1.tsv"'

echo "== locks"
mkdir -p "$OUT/MGScall/mag-rep/_state/S1.lock"; echo "$(hostname) $$ none" > "$OUT/MGScall/mag-rep/_state/S1.lock/owner"
run MGScall --mag-rep -o "$OUT" --samples S1 --force
t "sample locked by a live process: exit 1, lock kept" \
  '[[ $RC == 1 && $TERM_TXT == *"another job is working on this sample"* && -d $OUT/MGScall/mag-rep/_state/S1.lock ]]'
rm -rf "$OUT/MGScall/mag-rep/_state/S1.lock"
mkdir -p "$OUT/MGScall/mgs-rep/_state/S3.lock"; echo "$(hostname) 4194000 none" > "$OUT/MGScall/mgs-rep/_state/S3.lock/owner"
run MGScall --mgs-rep -o "$OUT" --samples S3 --force
t "lock of a dead process: taken over, sample redone, lock released" \
  '[[ $RC == 0 && ! -e $OUT/MGScall/mgs-rep/_state/S3.lock ]] && grep -q "taking over stale lock" "$LOG"'

echo "== outputs"
t "README files written" '[[ -s $OUT/README.md && -s $OUT/MGScall/README.md ]] && grep -q "mmt.log" "$OUT/README.md"'
t "no temp files left in OUT" '[[ -z $(find "$OUT" -name ".*.tmp" -o -name "*.lock") ]]'
t "work folders removed after success" '[[ -z $(find "$W/tmp/mmt" -mindepth 3 -maxdepth 3 -path "*_*/MGScall_*/*") ]]'

echo "== OUT/mmt.log: one line per event, every run"
L="$OUT/mmt.log"
t "every non-dry run has a START, every finished run an END" \
  '[[ $(grep -c "  START  " "$L") -ge 10 && $(grep -c "  START  " "$L") == $(( $(grep -c "  END  " "$L") + $(grep -c "  STOP  " "$L") )) ]]'
t "plan per MGS: where it is present (mgs-rep) / the representative's samples (mag-rep)" \
  'grep -q "MGScall --mgs-rep  PLAN  MGS.0: a MAG in the assembly of [0-9]* of 4 usable samples -> to do" "$L" && grep -q "MGScall --mag-rep  PLAN  MGS.1: representative from the assembly of 2: S2, S2.run2" "$L"'
t "one SAMPLE line per sample processed, with every MGS outcome" \
  'grep -q "MGScall --mgs-rep  SAMPLE  S3: MGS.0 done ([0-9]* sites, [0-9]* molecules)" "$L" && grep -q "MGScall --mag-rep  SAMPLE  S1: MGS.0 done" "$L"'
t "an index build and the dropped result are logged" 'grep -q "INDEX  index (re)built" "$L" && grep -q "SAMPLE  S1: MGS.1 dropped" "$L"'
t "dry runs are not logged" '! grep -q "dry" "$L"'

echo "== presence: an MGS with fewer molecules than MGSCALL_PRESENT_MOL is absent, before the raw BAM is read"
cp -a "$OUT" "$W/outp"; echo 'MGSCALL_PRESENT_MOL=100000' > "$W/present.conf"
"$MMT" MGScall --mgs-rep -o "$W/outp" --samples S3 --config "$W/test.conf" --config "$W/present.conf" > "$W/term.txt" 2>&1; RC=$?
PL=$(ls -t "$W"/outp/logs/MGScall/mgs-rep/*.log | head -1)
t "S3 absent: state absent, old table removed, raw BAM not read" \
  '[[ $RC == 0 && $(awk -F"\t" "NR > 1 {print \$2}" "$W/outp/MGScall/mgs-rep/_state/S3.tsv" | sort -u) == absent && ! -e $W/outp/MGScall/mgs-rep/MGS.0/S3.parquet ]] && grep -q "raw BAM not read" "$PL" && ! grep -q "S3: .*extracted" "$PL"'
t "mmt.log says absent with the molecule count" 'grep -q "SAMPLE  S3: MGS.0 absent ([0-9]* molecules < 100000)" "$W/outp/mmt.log"'
"$MMT" MGScall --mgs-rep -o "$W/outp" --samples S3 --config "$W/test.conf" --config "$W/present.conf" > "$W/term.txt" 2>&1
t "absent is a result: a rerun with the same settings has nothing to do" 'grep -q "plan: 0 sample(s) with work" "$W/term.txt"'

echo "== status"
run status -o "$OUT"
t "status: per MGS and mode, planned samples and their state with names" \
  '[[ $RC == 0 && $TERM_TXT == *"MGScall --mgs-rep | MGS.0: planned for 4 sample(s)"* && $TERM_TXT == *"MGScall --mgs-rep | MGS.0 done 4: S1, S2, S2.run2, S3"* && $TERM_TXT == *"MGScall --mag-rep | MGS.1 done 2: S2, S2.run2"* ]] && grep -q "status  STATUS  MGScall --mag-rep | MGS.0" "$OUT/mmt.log"'
run status -o "$W/outp" --select MGS.0,MGS.1
t "status --select: absent counted, only states that occur" '[[ $RC == 0 && $TERM_TXT == *"-> done 3, absent 1"* && $TERM_TXT == *"MGS.0 absent 1: S3"* && $TERM_TXT != *" : "* ]]'

echo
echo "$NP passed, $NF failed"
(( NF == 0 ))
