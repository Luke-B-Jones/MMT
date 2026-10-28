#!/usr/bin/env bash
# run_ctgmeth_test.sh WORKDIR - end-to-end test of `mmt derep-ctg` + `mmt ctg-meth` on a synthetic MATAFILER4 run with
# genes, catalogue, kinetics and known conformations (test/ctgmeth_mock.py ; truth checks in test/check_ctgmeth.py).
# WORKDIR is wiped. Needs samtools + minimap2 on PATH and:
#   MMT_TEST_PY      python >= 3.10 with pysam, pyarrow, numpy, edlib, duckdb (python module, for the checks)
#   MMT_TEST_IPD     a kineticsTools ipdSummary (MMT_TEST_IPD_MODEL its model) ; MMT_TEST_DUCKDB  DuckDB CLI
# SMRT Link's bystrandify / pbmm2 / pbindex are stand-ins (test/standins.sh). ~30 min on 2 cpus.
# shellcheck disable=SC2016,SC2034,SC2012  # checks are strings run by eval in t()
set -uo pipefail
W=$(readlink -m "$1"); HERE=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd); MMT="$HERE/../mmt"
PY=$(command -v "${MMT_TEST_PY:-python3}") || { echo "no python"; exit 2; }
IPD=${MMT_TEST_IPD:?set MMT_TEST_IPD}; DUCK=$(command -v "${MMT_TEST_DUCKDB:-duckdb}") || { echo "no duckdb"; exit 2; }
rm -rf "$W"; mkdir -p "$W"
NP=0 NF=0
t() { if eval "$2"; then echo "PASS  $1"; NP=$((NP + 1)); else echo "FAIL  $1"; NF=$((NF + 1)); fi; }
OUT="$W/out"; D="$OUT/ctg-meth"
run() {   # run COMMAND ARGS... -> RC, TERM_TXT, LOG (that run's log)
    local cmd=$1 before; shift
    before=$(find "$OUT/logs" -name '20*_*.log' 2> /dev/null | sort)     # run logs of every command (not ipdSummary/<S>.log)
    "$MMT" "$cmd" "$@" --config "$W/test.conf" > "$W/term.txt" 2>&1; RC=$?
    TERM_TXT=$(cat "$W/term.txt"); LOG=$(comm -13 <(echo "$before") <(find "$OUT/logs" -name '20*_*.log' 2> /dev/null | sort) | tail -1)
    echo "    \$ mmt $cmd $*  -> exit $RC, $(grep -c . "$W/term.txt") terminal line(s)"
}
mt() { stat -c %Y "$D/$1/$2.parquet" 2> /dev/null || echo none; }
q() { "$DUCK" -noheader -list -c "$1"; }
# fracLow / fracUp left out: ipdSummary bootstraps them (its random stream depends on how the work is chunked), so they
# differ slightly between runs ; frac and everything else is deterministic
kin_cols="MGS, sample, read_set, side, dist, MAG, contig, pos, strand, base, coverage, score, ipdRatio, tMean, tErr, modelPrediction, frac, call, identificationQv, pile_depth, n_ref, top_alt, n_alt, ctx_alt_pct, ctx_alt_n, clip_ctx_n"
same_kin() {   # same_kin FILE1 FILE2 [WHERE] : identical kinetics columns (row order independent)
    [[ $(q "SELECT count(*) FROM ((SELECT $kin_cols FROM '$1' ${3:-}) EXCEPT ALL (SELECT $kin_cols FROM '$2' ${3:-}))") == 0 &&
       $(q "SELECT count(*) FROM ((SELECT $kin_cols FROM '$2' ${3:-}) EXCEPT ALL (SELECT $kin_cols FROM '$1' ${3:-}))") == 0 ]]
}

echo "== mock MATAFILER4 run (genes, catalogue, kinetics, conformations) + stand-ins"
"$PY" "$HERE/ctgmeth_mock.py" "$W/mock" > /dev/null || { echo "mock failed"; exit 1; }
bash "$HERE/standins.sh" "$W/bin" "$PY" "$IPD" "${MMT_TEST_IPD_MODEL:-}"
cat > "$W/test.conf" << EOF
MMT_MAP="$W/mock/map.txt"
MMT_GENECAT="$W/mock/gc"
MMT_SMRT_BIN="$W/bin"
MMT_PYTHON="$PY"
MMT_DUCKDB="$DUCK"
MMT_TMP="$W/tmp"
MMT_THREADS=2
MMT_MEM_GB=4
EOF

echo "== help, option and setting errors, no derep-ctg result yet"
"$MMT" help ctg-meth > "$W/h.txt" 2>&1
t "help" 'grep -q "Usage: mmt ctg-meth" "$W/h.txt"'
run ctg-meth
t "no -o: usage error 2" '[[ $RC == 2 ]]'
for bad in 'CTGMETH_MIN_MOL=2' 'CTGMETH_MIN_FRAC=1.5' 'CTGMETH_WINDOW=abc' 'CTGMETH_MIN_ROW_COV=2'; do
    echo "$bad" > "$W/bad.conf"; run ctg-meth -o "$OUT" --config "$W/bad.conf"
    t "setting $bad: usage error 2" '[[ $RC == 2 && $TERM_TXT == *CTGMETH_* ]]'
done
run ctg-meth -o "$OUT"
t "before derep-ctg: exit 1, says which MGS lack a derep-ctg result" '[[ $RC == 1 && $TERM_TXT == *"no derep-ctg result"* ]]'

echo "== derep-ctg, then ctg-meth: listing, dry run"
run derep-ctg -o "$OUT"
t "derep-ctg: exit 0" '[[ $RC == 0 && -s $OUT/derep-ctg/MGS.0/lift.parquet ]]'
run ctg-meth -o "$OUT" --list-samples
t "--list-samples prints S1 .. S5" '[[ $RC == 0 && $(echo $TERM_TXT) == "S1 S2 S3 S4 S5" ]]'
run ctg-meth -o "$OUT" --dry-run
t "--dry-run: 7 MGS x sample results, nothing written" '[[ $RC == 0 && $TERM_TXT == *"(7 MGS x sample results)"* && ! -e $D/MGS.0 ]]'

echo "== full run with an injected ipdSummary failure in S3, then resume"
MMT_TEST_FAIL_IPD=S3 run ctg-meth -o "$OUT"
t "failure in S3: exit 1, others done, S3 named" '[[ $RC == 1 && $TERM_TXT == *"failed: S3"* && $(mt MGS.0 S1) != none && $(mt MGS.0 S5) != none && $(mt MGS.0 S3) == none ]]'
t "failed sample: no lock left, work folder removed (no --debug), failure in mmt.log" \
  '[[ ! -e $D/_state/S3.lock && -z $(ls -d "$W"/tmp/mmt/*_*/ctg-meth/S3 2> /dev/null) ]] && grep -q "ctg-meth  FAIL  S3: line [0-9]*: .*(exit 3)" "$OUT/mmt.log"'
MMT_TEST_FAIL_IPD=S3 run ctg-meth -o "$OUT" --samples S3 --debug
t "failure with --debug: work folder kept for the resume" '[[ $RC == 1 && -n $(ls "$W"/tmp/mmt/*_*/ctg-meth/S3/junctions.ok 2> /dev/null) ]]'
run ctg-meth -o "$OUT"
t "resume: only S3, from its finished steps" \
  '[[ $RC == 0 && $TERM_TXT == *"plan: 1 sample(s) with work"* && $(mt MGS.0 S3) != none ]] && ! grep -q "S3: .*extracted" "$LOG" && grep -q "S3: kinetics table ready" "$LOG"'
s0=$(date +%s); run ctg-meth -o "$OUT"
t "rerun: nothing to do, < 40 s" '[[ $RC == 0 && $TERM_TXT == *"plan: 0 sample(s) with work"* && $(( $(date +%s) - s0 )) -lt 40 ]]'

echo "== results against the truth"
"$PY" "$HERE/check_ctgmeth.py" "$W/mock" "$OUT" > "$W/truth.txt"; rc=$?
grep -E "^(FAIL|all|[0-9]+ FAILED)" "$W/truth.txt" | sed 's/^/    /'
t "truth: $(grep -c ^PASS "$W/truth.txt") checks pass" '[[ $rc == 0 ]]'

echo "== JUNCTIONS=0 (S1): read set 'all' unchanged, no read-set rows ; back to the default: the same tables again"
mkdir -p "$W/keep"; cp "$D/MGS.0/S1.parquet" "$W/keep/s1.parquet"; cp "$D/MGS.0/junctions/S1.parquet" "$W/keep/j1.parquet"
echo 'CTGMETH_JUNCTIONS=0' > "$W/j0.conf"
run ctg-meth -o "$OUT" --samples S1 --config "$W/j0.conf"
t "JUNCTIONS=0: S1 recomputed" '[[ $RC == 0 && $TERM_TXT == *"(2 MGS x sample results)"* ]]'
t "JUNCTIONS=0: only 'all' rows, identical to before ; read sets still listed, none called" \
  'same_kin "$D/MGS.0/S1.parquet" "$W/keep/s1.parquet" "WHERE read_set = '"'"'all'"'"'" && [[ $(q "SELECT count(*) FROM '"'"'$D/MGS.0/S1.parquet'"'"' WHERE read_set <> '"'"'all'"'"'") == 0 && $(q "SELECT count(*) || '"'"' '"'"' || count(*) FILTER (WHERE called) FROM '"'"'$D/MGS.0/junctions/S1.parquet'"'"'") == "$(q "SELECT count(*) FROM '"'"'$W/keep/j1.parquet'"'"'") 0" ]]'
run ctg-meth -o "$OUT" --samples S1
t "default again: same sites and read sets as the first run (deterministic)" \
  '[[ $RC == 0 ]] && same_kin "$D/MGS.0/S1.parquet" "$W/keep/s1.parquet" && [[ $(q "SELECT count(*) FROM ((SELECT read_set, side, contig, pos, n_molecules, call_bp FROM '"'"'$D/MGS.0/junctions/S1.parquet'"'"') EXCEPT (SELECT read_set, side, contig, pos, n_molecules, call_bp FROM '"'"'$W/keep/j1.parquet'"'"'))") == 0 ]]'

echo "== a newer derep-ctg result: re-lift only (no kinetics), then back"
for s in S1 S2 S3 S4 S5; do cp "$D/MGS.0/$s.parquet" "$W/keep/$s.parquet"; done
echo 'DEREP_MIN_BLOCK_GENES=12' > "$W/d12.conf"
run derep-ctg -o "$OUT" --config "$W/d12.conf"
t "derep-ctg with other blocks: done" '[[ $RC == 0 ]]'
s0=$(date +%s); run ctg-meth -o "$OUT"
t "ctg-meth: every result re-lifted, no kinetics, < 60 s" \
  '[[ $RC == 0 && $TERM_TXT == *"(7 MGS x sample results)"* && $(( $(date +%s) - s0 )) -lt 60 ]] && ! grep -q "kinetics table ready" "$LOG" && grep -q "re-lifted" "$LOG"'
t "re-lift: kinetics identical, lift_fp updated, blocks from the new derep-ctg" \
  '( for s in S1 S2 S3 S4 S5; do same_kin "$D/MGS.0/$s.parquet" "$W/keep/$s.parquet" || exit 1; done ) && [[ $(q "SELECT count(DISTINCT lift_fp) FROM read_parquet('"'"'$D/MGS.0/*.parquet'"'"')") == 1 && $(q "SELECT count(DISTINCT lift_fp) FROM read_parquet('"'"'$W/keep/*.parquet'"'"')") == 1 && $(q "SELECT min(lift_fp) FROM '"'"'$D/MGS.0/S1.parquet'"'"'") != $(q "SELECT min(lift_fp) FROM '"'"'$W/keep/S1.parquet'"'"'") ]]'
"$PY" "$HERE/check_ctgmeth.py" "$W/mock" "$OUT" > "$W/truth2.txt"; rc=$?
t "re-lifted results pass every truth check (block bases, lift, read sets)" '[[ $rc == 0 ]] || { grep ^FAIL "$W/truth2.txt"; false; }'
run derep-ctg -o "$OUT"
run ctg-meth -o "$OUT"
t "derep-ctg back at the defaults: re-lift again, blocks as at first" \
  '[[ $RC == 0 && -s $LOG ]] && ! grep -q "kinetics table ready" "$LOG" && [[ $(q "SELECT count(*) FROM ((SELECT block, instance, block_pos, block_strand, ctx_dist, gene FROM '"'"'$D/MGS.0/S2.parquet'"'"') EXCEPT ALL (SELECT block, instance, block_pos, block_strand, ctx_dist, gene FROM '"'"'$W/keep/S2.parquet'"'"'))") == 0 ]]'

echo "== an MGS whose derep-ctg result disappears is skipped, its results kept"
mv "$OUT/derep-ctg/MGS.1" "$W/keep/MGS.1"; m1=$(mt MGS.1 S1)
run ctg-meth -o "$OUT"
t "warning, MGS.1 skipped, its results kept" '[[ $RC == 0 && $TERM_TXT == *"no derep-ctg result for 1 MGS (MGS.1)"* && $(mt MGS.1 S1) == "$m1" ]]'
mv "$W/keep/MGS.1" "$OUT/derep-ctg/MGS.1"

echo "== MGScall on the same folder: shares the read index, plain layout"
run MGScall --mag-rep -o "$OUT" --select MGS.0 --samples S1
t "MGScall: exit 0, read index reused, result in MGScall/mag-rep/MGS.0/S1.parquet" \
  '[[ $RC == 0 && -s $OUT/MGScall/mag-rep/MGS.0/S1.parquet && -s $LOG ]] && ! grep -q "assigning reads" "$LOG" && [[ $(q "SELECT min(MGS) || '"'"' '"'"' || min(sample) FROM '"'"'$OUT/MGScall/mag-rep/MGS.0/S1.parquet'"'"'") == "MGS.0 S1" ]]'

echo "== outputs"
t "README files, mmt.log (plan per MGS, a SAMPLE line per sample, END)" '[[ -s $D/README.md && -s $OUT/README.md ]] && grep -q "ctg-meth  PLAN  MGS.0: a MAG in the assembly of" "$OUT/mmt.log" && grep -q "ctg-meth  SAMPLE  S1: MGS.0 done (" "$OUT/mmt.log" && grep -q "ctg-meth  END  " "$OUT/mmt.log"'
t "no temp file or lock left in OUT" '[[ -z $(find "$OUT" -name ".*.tmp" -o -name "*.lock" -o -name ".*.old") ]]'
t "work folders removed after success" '[[ -z $(ls -d "$W"/tmp/mmt/*_*/ctg-meth/* 2> /dev/null) ]]'

echo "== presence: an MGS with fewer molecules than CTGMETH_PRESENT_MOL is absent (no files), raw BAM not read"
cp -a "$OUT" "$W/outp"; echo 'CTGMETH_PRESENT_MOL=100000' > "$W/present.conf"
"$MMT" ctg-meth -o "$W/outp" --samples S5 --config "$W/test.conf" --config "$W/present.conf" > "$W/term.txt" 2>&1; rc=$?
PL=$(ls -t "$W"/outp/logs/ctg-meth/*.log | head -1)
t "S5: every MGS absent, its files removed, raw BAM not read, logged" \
  '[[ $rc == 0 && $(awk -F"\t" "NR > 1 {print \$2}" "$W/outp/ctg-meth/_state/S5.tsv" | sort -u) == absent && -z $(ls "$W"/outp/ctg-meth/MGS.*/S5.parquet 2> /dev/null) ]] && grep -q "raw BAM not read" "$PL" && grep -q "SAMPLE  S5: MGS.[0-9] absent (" "$W/outp/mmt.log"'
"$MMT" ctg-meth -o "$W/outp" --samples S5 --config "$W/test.conf" --config "$W/present.conf" > "$W/term.txt" 2>&1
t "absent is a result: rerun has nothing to do" 'grep -q "plan: 0 sample(s) with work" "$W/term.txt"'
"$MMT" status -o "$W/outp" --config "$W/test.conf" > "$W/term.txt" 2>&1
t "status lists ctg-meth per MGS" 'grep -q "ctg-meth | MGS.0: planned for" "$W/term.txt" && grep -q "ctg-meth | MGS.0 absent 1: S5" "$W/term.txt" && grep -q "derep-ctg | done" "$W/term.txt"'
t "index: movies.tsv + depth bins of every MAG contig" \
  '[[ $(tail -n +2 "$OUT/index/movies.tsv" | wc -l) == 5 && $(q "SELECT count(DISTINCT contig) FROM read_parquet('"'"'$OUT/index/depth/*.parquet'"'"')") -ge 9 ]]'

echo
echo "$NP passed, $NF failed"
(( NF == 0 ))
