#!/usr/bin/env bash
# run_motexport_test.sh WORKDIR - end-to-end test of `mmt export-mot` on a synthetic MATAFILER4 run
#   (test/motexport_mock.py). WORKDIR is wiped. Needs samtools on PATH and python with pysam, pyarrow, numpy
#   (MMT_TEST_PY, default python3). With MOT_BIN (a mot binary) the tables are also checked by `mot plan`, mot's own
#   validation of an input table (files, BAM headers, read-name lists, molecules in two samples). ~1 min.
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
    "$MOT_BIN" plan -i "$1" -q > "$1.motplan" 2>&1
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
    t "L: mot plan lists every mot sample" '(( $(grep -c . <(cut -f1 "$W/clean/mot_input.tsv" | tail -n +2 | sort -u)) == 8 ))'
    # mot's read-name loader on every selection (ZmwSet::load via `mot plan -b ... -r ...` is per sample; the table
    # check above already loaded each reads file once)
fi

echo "== --exclude-qcfail"
run "$W/noqc/mot_input.tsv" --exclude-qcfail
t "exit 3" '(( RC == 3 ))'
"$PY" "$HERE/check_motexport.py" "$M" "$W/noqc/mot_input.tsv" noqc | sed 's/^/  /'; t "checks (noqc)" '(( PIPESTATUS[0] == 0 ))'

echo "== a broken table is refused by mot (the test can fail)"
if [[ -n ${MOT_BIN:-} ]]; then
    # row 3 gets row 2's BAM and a copy of row 2's selection (another file, the same molecules)
    cp "$(awk -F'\t' 'NR == 2 {print $5}' "$W/clean/mot_input.tsv")" "$W/copy.reads.txt.gz"
    awk -F'\t' -v OFS='\t' -v c="$W/copy.reads.txt.gz" 'NR == 2 {b = $2} NR == 3 {$2 = b; $5 = c} 1' \
        "$W/clean/mot_input.tsv" > "$W/bad.tsv"
    t "L (negative): two samples selecting the same molecules of one BAM -> mot plan refuses" '[[ -s $W/bad.tsv ]] && ! mot_plan "$W/bad.tsv" && grep -q "same molecules" "$W/bad.tsv.motplan"'
fi
echo "$NP passed, $NF failed"
exit $(( NF > 0 ))
