#!/usr/bin/env bash
# run_motexport_test.sh WORKDIR - end-to-end test of `mmt export-mot` on a synthetic MATAFILER4 run
#   (test/motexport_mock.py). WORKDIR is wiped. Needs samtools on PATH and python with pysam, pyarrow, numpy
#   (MMT_TEST_PY, default python3). With MOT_BIN, a temporary six-column expansion is checked by older `mot plan`;
#   the shipped schema-2 table itself requires an updated MOT reader. ~1 min.
set -uo pipefail
W=$(readlink -m "$1"); HERE=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd); MMT="$HERE/../mmt"
PY=${MMT_TEST_PY:-python3}
rm -rf "$W"; mkdir -p "$W"
NP=0 NF=0
t() { if eval "$2"; then echo "PASS  $1"; NP=$((NP + 1)); else echo "FAIL  $1"; NF=$((NF + 1)); fi; }
"$PY" "$HERE/motexport_mock.py" "$W/mock" > /dev/null || { echo "mock failed"; exit 1; }
M=$W/mock
printf 'MMT_PYTHON="%s"\n' "$(command -v "$PY")" > "$W/test.conf"
run() {   # run OUT.tsv ARGS... -> RC
    mkdir -p "$(dirname "$1")"
    "$MMT" export-mot -i "$M/mf4/run" -gc "$M/gc" -o "$1" --config "$W/test.conf" -t 2 "${@:2}" > "$1.term" 2>&1; RC=$?
    echo "    \$ mmt export-mot ... -o $(basename "$1") ${*:2}  -> exit $RC"
}
mot_plan() {   # mot_plan TABLE -> 0 when mot accepts it
    [[ -n ${MOT_BIN:-} ]] || return 0
    "$PY" - "$1" "$1.legacy.tsv" <<'PY'
import csv, json, sys
with open(sys.argv[1]) as src, open(sys.argv[2], 'w') as out:
    fields = ['sample', 'bam', 'ref', 'role', 'reads', 'group']
    w = csv.DictWriter(out, fields, delimiter='\t', extrasaction='ignore', lineterminator='\n'); w.writeheader()
    for r in csv.DictReader(src, delimiter='\t'):
        if r['mot_status'] == 'ready':
            for b, reads in zip(json.loads(r['bam']), json.loads(r['reads'])):
                w.writerow(dict(r, bam=b, reads=reads))
PY
    "$MOT_BIN" plan -i "$1.legacy.tsv" -q > "$1.motplan" 2>&1
}

echo "== whole cohort (map from GeneCat's records): PARTIAL expected"
run "$W/full/mot_input.tsv"
t "exit 3 (PARTIAL)" '(( RC == 3 ))'
t "terminal report ends with Status: PARTIAL" 'tail -1 "$W/full/mot_input.tsv.term" | grep -qx "Status: PARTIAL"'
t "map taken from GeneCat's GCmaps.ori" 'grep -q "GCmaps.ori" "$W/full/mot_input.tsv.term"'
"$PY" "$HERE/check_motexport.py" "$M" "$W/full/mot_input.tsv" full | sed 's/^/  /'; t "checks (full)" '(( PIPESTATUS[0] == 0 ))'
if [[ -n ${MOT_BIN:-} ]]; then
    t "L: mot plan accepts the PARTIAL table (only validated rows are written)" 'mot_plan "$W/full/mot_input.tsv"'
fi

echo "== rerun: nothing changed -> every sample and assembly reused, identical outputs"
cp "$W/full/mot_input.tsv" "$W/first.tsv"
find "$W/full/mot_input.motlib" -type f \( -name '*.reads.txt.gz' -o -name '*.fa.gz' \) -exec md5sum {} + | sort > "$W/first.md5"
run "$W/full/mot_input.tsv"
t "rerun: same exit, every sample reused" '(( RC == 3 )) && (( $(grep -c "(reused)" "$W/full/mot_input.tsv.term") == 10 ))'
t "rerun: identical table and files" 'cmp -s "$W/first.tsv" "$W/full/mot_input.tsv" && cmp -s "$W/first.md5" <(find "$W/full/mot_input.motlib" -type f \( -name "*.reads.txt.gz" -o -name "*.fa.gz" \) -exec md5sum {} + | sort)'
touch "$M/mf4/run/S2/mapping/S2-smd.cram"
run "$W/full/mot_input.tsv"
t "a changed CRAM redoes that sample only" '(( $(grep -c "(reused)" "$W/full/mot_input.tsv.term") == 9 )) && grep -q "S2: .*s)$" "$W/full/mot_input.tsv.term"'

echo "== clean map (S1-S5): COMPLETE expected"
run "$W/clean/mot_input.tsv" --map "$M/map_clean.txt"
t "exit 0 (COMPLETE)" '(( RC == 0 ))'
"$PY" "$HERE/check_motexport.py" "$M" "$W/clean/mot_input.tsv" clean | sed 's/^/  /'; t "checks (clean)" '(( PIPESTATUS[0] == 0 ))'
if [[ -n ${MOT_BIN:-} ]]; then
    t "L: mot plan accepts the table (exit 0)" 'mot_plan "$W/clean/mot_input.tsv"'
    # mot's read-name loader on every selection (ZmwSet::load via `mot plan -b ... -r ...` is per sample; the table
    # check above already loaded each reads file once)
fi

echo "== --exclude-qcfail"
run "$W/noqc/mot_input.tsv" --exclude-qcfail
t "exit 3" '(( RC == 3 ))'
"$PY" "$HERE/check_motexport.py" "$M" "$W/noqc/mot_input.tsv" noqc | sed 's/^/  /'; t "checks (noqc)" '(( PIPESTATUS[0] == 0 ))'

echo "== a broken table is refused by mot (the test can fail)"
if [[ -n ${MOT_BIN:-} ]]; then
    "$PY" - "$W/clean/mot_input.tsv" "$W/bad.tsv" <<'PY'
import csv, json, sys
with open(sys.argv[1]) as src:
    rows = list(csv.DictReader(src, delimiter='\t')); fields = list(rows[0])
active = [r for r in rows if r['mot_status'] == 'ready']
active[1]['bam'], active[1]['reads'] = active[0]['bam'], active[0]['reads']
with open(sys.argv[2], 'w') as out:
    w = csv.DictWriter(out, fields, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
PY
    t "L (negative): two samples selecting the same molecules of one BAM -> mot plan refuses" '[[ -s $W/bad.tsv ]] && ! mot_plan "$W/bad.tsv" && grep -q "same molecules" "$W/bad.tsv.motplan"'
fi
echo "== missing coverage stays NA and PARTIAL; molecule selections remain usable"
mv "$M/mf4/run/S1/mapping/S1-smd.bam.coverage.gz" "$M/S1.depth.gz"
run "$W/missing/mot_input.tsv" --map "$M/map_clean.txt"
t "missing depth -> PARTIAL" '(( RC == 3 ))'
"$PY" - "$W/missing/mot_input.tsv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
r = next(r for r in rows if r['group'] == 'S1' and r['MAG'] == 'S1__SemiBin_1')
assert r['normalized_abundance'] == 'NA' and r['abundance_status'] == 'unavailable'
assert r['mot_status'] == 'ready'
PY
t "missing depth never becomes zero" '(( $? == 0 ))'
printf 'sample\tcoverage\nS1\tS1.depth.gz\n' > "$M/coverage.tsv"
run "$W/override/mot_input.tsv" --map "$M/map_clean.txt" --coverage-map "$M/coverage.tsv"
t "explicit relative coverage override restores COMPLETE" '(( RC == 0 ))'
"$PY" -m unittest discover -s "$HERE" -p test_motabundance.py -q
t "analytical abundance tests" '(( $? == 0 ))'
echo "$NP passed, $NF failed"
exit $(( NF > 0 ))
