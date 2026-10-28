#!/usr/bin/env bash
# run_derep_test.sh WORKDIR - end-to-end test of `mmt derep-ctg` on a synthetic MATAFILER4 + GeneCat run
# with known blocks (test/derep_mock.py, truth checks in test/check_derep.py). WORKDIR is wiped.
#   Needs samtools + minimap2 on PATH and MMT_TEST_PY = a python >= 3.10 with pysam, pyarrow, pandas, edlib
#   (default python3). Takes ~1 min.
# shellcheck disable=SC2016,SC2034  # checks are strings run by eval in t()
set -uo pipefail
W=$(readlink -m "$1"); HERE=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd); MMT="$HERE/../mmt"
PY=$(command -v "${MMT_TEST_PY:-python3}") || { echo "no python"; exit 2; }
rm -rf "$W"; mkdir -p "$W"
NP=0 NF=0
t() { if eval "$2"; then echo "PASS  $1"; NP=$((NP + 1)); else echo "FAIL  $1"; NF=$((NF + 1)); fi; }
OUT="$W/out"; D="$OUT/derep-ctg"
run() {   # run ARGS... -> RC, TERM_TXT (terminal output)
    "$MMT" derep-ctg "$@" --config "$W/test.conf" > "$W/term.txt" 2>&1; RC=$?
    TERM_TXT=$(cat "$W/term.txt")
    echo "    \$ mmt derep-ctg $*  -> exit $RC, $(grep -c . "$W/term.txt") terminal line(s)"
}
conf() { printf '%s\n' "$@" > "$W/extra.conf"; }        # conf 'DEREP_X=1' ... : settings for the next run (--config extra)
mt() { stat -c %Y "$D/$1/summary.tsv" 2> /dev/null || echo none; }
same_results() {   # same_results DIR1 DIR2 : identical block / instance / gene tables for every MGS
    "$PY" - "$1" "$2" << 'EOF'
import sys, glob, os, pyarrow.parquet as pq
a, b = sys.argv[1:3]
ls = lambda x: sorted(os.path.basename(d) for d in glob.glob(f"{x}/*") if os.path.isfile(f"{d}/.fp"))
ms = ls(a)
assert ms and ms == ls(b), "different MGS"
for m in ms:
    for f in ("blocks", "block_genes", "instances", "genes", "links", "families", "lift", "composition", "contig_seqs"):
        x, y = pq.read_table(f"{a}/{m}/{f}.parquet"), pq.read_table(f"{b}/{m}/{f}.parquet")
        assert x.equals(y), f"{m}/{f} differs"
EOF
}

echo "== mock MATAFILER4 run + GeneCat"
"$PY" "$HERE/derep_mock.py" "$W/mock" > /dev/null || { echo "mock failed"; exit 1; }
cat > "$W/test.conf" << EOF
MMT_MAP="$W/mock/map.txt"
MMT_GENECAT="$W/mock/gc"
MMT_PYTHON="$PY"
MMT_TMP="$W/tmp"
MMT_THREADS=2
MMT_MEM_GB=4
EOF
conf ''
GC="$W/mock/gc"; CIDX="$GC/compl.incompl.95.fna.clstr.idx"; MVG="$GC/Bin_SB/MAGvsGC.txt.gz"; MAP="$W/mock/map.txt"
asm_of() { awk -F'\t' -v s="$1" '$1 == s {print $5}' "$OUT/index/samples.tsv"; }   # assembly folder of a sample

echo "== help and option errors"
"$MMT" help dereplicate-contigs > "$W/h.txt" 2>&1
t "help (old name accepted)" 'grep -q "Usage: mmt derep-ctg" "$W/h.txt"'
run
t "no -o: usage error 2" '[[ $RC == 2 ]]'
run -o "$OUT" --bogus
t "unknown option: usage error 2" '[[ $RC == 2 && $TERM_TXT == *"unknown option"* ]]'
run -o "$OUT" --part x
t "--part x: usage error 2" '[[ $RC == 2 ]]'
for bad in 'DEREP_SYNTENY=0.5' 'DEREP_MAX_SKIP=1.5' 'DEREP_FAM_MIN_ID=abc' 'DEREP_FAM_MIN_ID=1.2' 'DEREP_MAX_COPY=1' 'DEREP_MIN_EDGE_MAGS=0'; do
    conf "$bad"; run -o "$OUT" --config "$W/extra.conf"
    t "setting $bad: usage error 2" '[[ $RC == 2 && $TERM_TXT == *DEREP_* ]]'
done
conf ''

echo "== listing, dry run, prep only"
run -o "$OUT" --list-mgs
t "--list-mgs prints exactly the 3 MGS" '[[ $RC == 0 && $(echo $TERM_TXT) == "MGS.0 MGS.1 MGS.2" ]]'
run -o "$OUT" --dry-run
t "--list-mgs / --dry-run: nothing done, not even the prep" '[[ $RC == 0 && $TERM_TXT == *"dry run"* && ! -e $D/MGS.0 && ! -e $D/_prep ]]'
run -o "$OUT" --prep-only -t 4
t "--prep-only: prep ready, no MGS done" '[[ $RC == 0 && -s $D/_prep/assemblies.tsv && ! -e $D/MGS.0 ]]'
run -o "$OUT" --mgs MGS.9
t "--mgs unknown MGS: exit 1 with message" '[[ $RC == 1 && $TERM_TXT == *"not in MAGvsGC: MGS.9"* ]]'
run -o "$OUT" --part 3/2
t "--part 3/2: exit 1" '[[ $RC == 1 ]]'
run -o "$OUT" --part 4/5 --dry-run
t "--part 4/5 of 3 MGS: an empty part is not an error" '[[ $RC == 0 && $TERM_TXT == *"holds no MGS"* ]]'
run -o "$OUT" --mgs "MGS.1 MGS.2" --dry-run
t "--mgs with spaces: both MGS" '[[ $RC == 0 && $TERM_TXT == *"(selection: 2)"* ]]'
run -o "$OUT" -t abc --dry-run
t "-t abc: usage error 2, message on the terminal" '[[ $RC == 2 && $TERM_TXT == *"threads must be a positive number"* ]]'
conf 'DEREP_MAX_SKIP=03'; run -o "$OUT" --config "$W/extra.conf" --dry-run; conf ''
t "DEREP_MAX_SKIP=03: usage error 2" '[[ $RC == 2 ]]'

echo "== SLURM-like: two parts at once on a fresh folder (shared index + prep build), then on the prepared one"
O3="$W/out3"
( "$MMT" derep-ctg -o "$O3" --part 1/2 --config "$W/test.conf" > "$W/q1.txt" 2>&1; echo $? > "$W/q1.rc" ) &
( "$MMT" derep-ctg -o "$O3" --part 2/2 --config "$W/test.conf" > "$W/q2.txt" 2>&1; echo $? > "$W/q2.rc" ) &
wait
t "two parts from scratch: both exit 0, prep built once, all 3 MGS" \
  '[[ $(cat "$W/q1.rc")$(cat "$W/q2.rc") == 00 && $(cat "$W"/q?.txt | grep -c "prep: catalogue") == 1 && $(ls -d "$O3"/derep-ctg/MGS.* | wc -l) == 3 ]]'

( "$MMT" derep-ctg -o "$OUT" --part 1/2 --config "$W/test.conf" > "$W/p1.txt" 2>&1; echo $? > "$W/p1.rc" ) &
( "$MMT" derep-ctg -o "$OUT" --part 2/2 --config "$W/test.conf" > "$W/p2.txt" 2>&1; echo $? > "$W/p2.rc" ) &
wait
t "two parts: both exit 0, all 3 MGS done, no lock / temp left" \
  '[[ $(cat "$W/p1.rc")$(cat "$W/p2.rc") == 00 && -s $D/MGS.0/blocks.parquet && -s $D/MGS.1/blocks.parquet && -s $D/MGS.2/blocks.parquet ]] && [[ -z $(ls -d "$D"/.lock* "$D"/.MGS* "$D"/.prep* "$W"/tmp/mmt/derep_* 2> /dev/null) ]]'
t "parts are balanced: MGS.0 (20 MAGs) alone in part 1" 'grep -q "plan: 1 MGS to do" "$W/p1.txt" && grep -q "plan: 2 MGS to do" "$W/p2.txt"'
run -o "$OUT"
t "rerun: nothing to do, summary of all 3 MGS" '[[ $RC == 0 && $TERM_TXT == *"plan: 0 MGS to do, 3 up to date"* && $(grep -c . "$D/summary.tsv") == 4 ]]'
t "README.md and mmt.log written (START, MGS done, END)" '[[ -s $D/README.md && -s $OUT/README.md ]] && grep -q "derep-ctg  MGS  MGS.0 done: 20 MAGs" "$OUT/mmt.log" && grep -q "derep-ctg  END  " "$OUT/mmt.log" && [[ $(grep -c "derep-ctg  START" "$OUT/mmt.log") == $(( $(grep -c "derep-ctg  END" "$OUT/mmt.log") + $(grep -c "derep-ctg  STOP" "$OUT/mmt.log") )) ]]'

echo "== results against the truth"
"$PY" "$HERE/check_derep.py" "$W/mock" "$D" > "$W/truth.txt"; rc=$?
grep -E "^(FAIL|all|[0-9]+ FAILED)" "$W/truth.txt" | sed 's/^/    /'
t "truth: $(grep -c ^PASS "$W/truth.txt") checks pass" '[[ $rc == 0 ]]'

echo "== determinism: a fresh output folder (relative -o, one job) gives identical tables"
( cd "$W" && "$MMT" derep -o out_b --config "$W/test.conf" > "$W/b.txt" 2>&1 ); rc=$?
t "fresh run: exit 0, same tables (also as the concurrent parts run)" '[[ $rc == 0 ]] && same_results "$D" "$W/out_b/derep-ctg" && same_results "$O3/derep-ctg" "$D"'

echo "== selection, force"
m0=$(mt MGS.0) m1=$(mt MGS.1) m2=$(mt MGS.2); sleep 1
printf '# comment\nMGS.1 anything\n\n' > "$W/ids.txt"
run -o "$OUT" --mgs "$W/ids.txt" --force
t "--mgs FILE --force: only MGS.1 redone" '[[ $RC == 0 && $(mt MGS.0) == "$m0" && $(mt MGS.1) != "$m1" && $(mt MGS.2) == "$m2" ]]'

echo "== settings change -> all redone ; DEREP_IDENTITY=0 needs no minimap2 ; DEREP_FAM_MIN_ID=0 merges nothing"
conf 'DEREP_IDENTITY=0' 'DEREP_FAM_MIN_ID=0' 'MMT_MINIMAP2=/nonexistent/minimap2'
run -o "$OUT" --config "$W/extra.conf"
t "settings change: 3 MGS redone" '[[ $RC == 0 && $TERM_TXT == *"plan: 3 MGS to do"* ]]'
t "identity 0: identities and lift of representatives only ; fam_min_id 0: no merges" '"$PY" -c "
import pandas as pd; d=\"$D/MGS.0\"; i = pd.read_parquet(f\"{d}/instances.parquet\")
assert i[~i.is_rep].identity.isna().all() and set(pd.read_parquet(f\"{d}/lift.parquet\").instance) == set(i[i.is_rep].instance)
assert len(pd.read_parquet(f\"{d}/merged_catalogue_genes.parquet\")) == 0"'
t "settings recorded in summary.tsv" 'grep -q "identity=0" "$D/MGS.0/summary.tsv"'

echo "== failure in one MGS (minimap2 fails for MGS.1), resume"
mkdir -p "$W/bin"; printf '#!/bin/bash\n[[ "$*" == *derep_MGS.1_* ]] && { echo "injected failure" >&2; exit 3; }\nexec %s "$@"\n' "$(command -v minimap2)" > "$W/bin/mm2"
chmod +x "$W/bin/mm2"
conf "MMT_MINIMAP2=$W/bin/mm2"
run -o "$OUT" --config "$W/extra.conf"
t "failure: exit 1, MGS.1 named, others done" '[[ $RC == 1 && $TERM_TXT == *"1 MGS failed: MGS.1"* ]] && grep -q "identity=1" "$D/MGS.0/summary.tsv" && grep -q "identity=0" "$D/MGS.1/summary.tsv"'
t "failed MGS: previous result kept, no lock / temp left" '[[ -s $D/MGS.1/blocks.parquet && -z $(ls -d "$D"/.lock* "$D"/.MGS* "$W"/tmp/mmt/derep_* 2> /dev/null) ]]'
conf ''
run -o "$OUT"
t "resume: only MGS.1 redone" '[[ $RC == 0 && $TERM_TXT == *"plan: 1 MGS to do, 2 up to date"* ]] && grep -q "identity=1" "$D/MGS.1/summary.tsv"'
t "back at the defaults: same tables as the first run" 'same_results "$D" "$W/out_b/derep-ctg"'

echo "== locks"
sleep 25 & SP=$!
mkdir "$D/.lock.MGS.2"; echo "$(hostname) $SP none" > "$D/.lock.MGS.2/owner"
s0=$(date +%s); run -o "$OUT" --mgs MGS.2 --force
t "MGS locked by a live job: waits for it, then takes the lock over" \
  '[[ $RC == 0 && $(( $(date +%s) - s0 )) -ge 15 ]] && grep -q "another job is working on it" "$OUT"/logs/derep-ctg/*.log'
wait "$SP" 2> /dev/null
sleep 0 & DP=$!; wait "$DP"
mkdir "$D/.lock.MGS.2"; echo "$(hostname) $DP none" > "$D/.lock.MGS.2/owner"
run -o "$OUT" --mgs MGS.2 --force
t "lock of a dead job is taken over at once" '[[ $RC == 0 && ! -e $D/.lock.MGS.2 && ! -e $D/.lock.MGS.2.takeover ]]'

echo "== input changes"
cp "$CIDX" "$W/clstr.idx.bak"
gzip "$CIDX"
run -o "$OUT"
t "catalogue index gzipped: prep rebuilt, all redone, same tables" \
  '[[ $RC == 0 && $TERM_TXT == *"prep: catalogue"* && $TERM_TXT == *"plan: 3 MGS to do"* ]] && same_results "$D" "$W/out_b/derep-ctg"'
rm "$CIDX.gz"; mv "$W/clstr.idx.bak" "${CIDX%.clstr.idx}.ctsv.idx"
run -o "$OUT"
t ".ctsv.idx naming: prep rebuilt, same tables" '[[ $RC == 0 && $TERM_TXT == *"prep: catalogue"* ]] && same_results "$D" "$W/out_b/derep-ctg"'
mv "${CIDX%.clstr.idx}.ctsv.idx" "$CIDX"
A=$(asm_of P03)
gzip "$A/scaffolds.fasta.filt"
run -o "$OUT"
t "an assembly gzipped (not bgzip): prep rebuilt, same tables" \
  '[[ $RC == 0 && $TERM_TXT == *"prep: catalogue"* ]] && same_results "$D" "$W/out_b/derep-ctg"'
B=$(asm_of P04); PROT="$GC/compl.incompl.95.prot.faa"; cp "$PROT" "$W/prot.bak"
bgzip "$B/scaffolds.fasta.filt" && bgzip "$PROT"
run -o "$OUT"
t "an assembly + the catalogue proteins bgzipped: same tables, no index written next to them" \
  '[[ $RC == 0 ]] && same_results "$D" "$W/out_b/derep-ctg" && [[ -z $(ls "$B"/*.gzi "$B"/*.fai "$GC"/*.gzi "$GC"/*.fai 2> /dev/null) ]]'
rm "$PROT.gz"; cp "$W/prot.bak" "$PROT"
"$PY" - "$PROT" << 'EOF'
import sys   # every catalogue id gets the next one's protein: the catalogue no longer describes the genes
R = open(sys.argv[1]).read().split(">")[1:]
open(sys.argv[1], "w").write("".join(">" + a.split("\n", 1)[0] + "\n" + b.split("\n", 1)[1] for a, b in zip(R, R[1:] + R[:1])))
EOF
run -o "$OUT"
t "catalogue proteins that do not belong to the genes: prep stops and says so" '[[ $RC == 1 && $TERM_TXT == *"do not resemble the protein of their catalogue gene"* ]]'
cp "$W/prot.bak" "$PROT"
GFF=$(ls "$A"/genePred/genes.gff*); cp "$GFF" "$W/gff.bak"
"$PY" - "$GFF" "$OUT/index/contigs.tsv.gz" << 'EOF'
import gzip, sys
f, idx = sys.argv[1:3]
op = (lambda p, m: gzip.open(p, m + "t")) if f.endswith(".gz") else open
mag = {l.split("\t")[1] for l in gzip.open(idx, "rt")}
lines = list(op(f, "r")); i = next(i for i, l in enumerate(lines) if l.split("\t")[0] in mag and "\tCDS\t" in l and "partial=00" in l)
del lines[i]; op(f, "w").writelines(lines)
EOF
cp -r "$D/_prep" "$W/prep.before"
run -o "$OUT"
t "GFF lost a catalogue gene: prep fails with a clear message, old prep kept, no lock" \
  '[[ $RC == 1 && $TERM_TXT == *"prep failed"* && ! -e $D/.prep.lock && ! -e $D/.prep.build ]] && grep -q "catalogue member genes not in" "$OUT"/logs/derep-ctg/*.log && diff -rq "$D/_prep" "$W/prep.before" > /dev/null'
cp "$W/gff.bak" "$GFF"
run -o "$OUT"
t "GFF restored: runs again" '[[ $RC == 0 ]]'
zcat "$MVG" | awk -F'\t' '$2 != "MGS.2"' | gzip > "$W/mvg.gz" && cp "$MVG" "$W/mvg.bak" && cp "$W/mvg.gz" "$MVG"
run -o "$OUT" --mgs MGS.1
t "MGS.2 gone from MAGvsGC, selection run: MGS.2 folder kept" '[[ $RC == 0 && -d $D/MGS.2 ]]'
run -o "$OUT"
t "full run: MGS.2 dropped, summary has 2 MGS" '[[ $RC == 0 && $TERM_TXT == *"dropped MGS.2"* && ! -e $D/MGS.2 && $(grep -c . "$D/summary.tsv") == 3 ]]'
cp "$W/mvg.bak" "$MVG"

echo "== what redoes what"
run -o "$OUT"
t "inputs restored: all redone" '[[ $RC == 0 && $TERM_TXT == *"plan: 3 MGS to do"* ]]'
cp "$MAP" "$W/map.bak"
touch "$MVG"; awk -F'\t' -v OFS='\t' '/^#SmplID/ {print $0, "Participant"; next} /^#/ {print; next} {print $0, "X" NR}' "$W/map.bak" > "$MAP"
run -o "$OUT"
t "MAGvsGC touched + a new map column: nothing redone, prep kept" '[[ $RC == 0 && $TERM_TXT == *"plan: 0 MGS to do"* && $TERM_TXT != *"prep: catalogue"* ]]'
samtools faidx "$(asm_of P05)/scaffolds.fasta.filt"
run -o "$OUT"
t "a .fai written next to an assembly: nothing redone" '[[ $RC == 0 && $TERM_TXT == *"plan: 0 MGS to do"* && $TERM_TXT != *"prep: catalogue"* ]]'
zcat "$W/mvg.bak" | awk -F'\t' -v OFS='\t' '$1 == "P01__SemiBin_2" {$7 = 50} 1' | gzip > "$MVG"
m0=$(mt MGS.0) m1=$(mt MGS.1); sleep 1
run -o "$OUT"
t "completeness of one MGS.1 MAG changed: only MGS.1 redone" '[[ $RC == 0 && $TERM_TXT == *"plan: 1 MGS to do"* && $(mt MGS.0) == "$m0" && $(mt MGS.1) != "$m1" ]]'
BF="$(asm_of P02)/Binning/SB/P02"; cp "$BF" "$W/bin.bak"
awk -F'\t' -v OFS='\t' '$2 == "SemiBin_2" && !done {$2 = "0"; done = 1} 1' "$W/bin.bak" > "$BF"
run -o "$OUT"
t "a bin file changed (re-binning): index + prep rebuilt, MGS redone" '[[ $RC == 0 && $TERM_TXT == *"prep: catalogue"* && $TERM_TXT == *"plan: 3 MGS to do"* ]]'
cp "$W/map.bak" "$MAP"; cp "$W/mvg.bak" "$MVG"; cp "$W/bin.bak" "$BF"

echo "== a MAG of three co-binned strains: P14 + the MGS.0 contigs of P12 and P13"
"$PY" "$HERE/derep_mock.py" --mix "$W/mock" P14 P12 P13 || echo "mix failed"
run -o "$OUT"
"$PY" "$HERE/check_derep.py" "$W/mock" "$D" > "$W/truth2.txt"; rc=$?
grep -E "^FAIL" "$W/truth2.txt" | sed 's/^/    /'
t "mixed MAG: flagged, its copies make no linker, every truth check still passes" \
  '[[ $RC == 0 && $rc == 0 ]] && "$PY" -c "import pandas as pd; m = pd.read_parquet(\"$D/MGS.0/mags.parquet\"); assert m.set_index(\"MAG\").mixed.to_dict() == {g: g == \"P14__SemiBin_1\" for g in m.MAG}"'

echo "== a results folder of MMT <= 0.2 (OUT/dereplicate-contigs, MGS=<MGS> folders) is moved, its prep reused"
O4="$W/out4"; mkdir -p "$O4"; cp -r "$OUT/index" "$O4/"; cp -r "$D" "$O4/dereplicate-contigs"
for d in "$O4"/dereplicate-contigs/MGS.*; do mv "$d" "$O4/dereplicate-contigs/MGS=${d##*/}"; done
"$MMT" derep-ctg -o "$O4" --config "$W/test.conf" > "$W/term.txt" 2>&1; RC=$?; TERM_TXT=$(cat "$W/term.txt")
t "old folder moved, prep reused, MGS redone in the new layout, no MGS= folder left" \
  '[[ $RC == 0 && $TERM_TXT == *"moved "* && $TERM_TXT != *"prep: catalogue"* && $TERM_TXT == *"plan: 3 MGS to do"* && ! -e $O4/dereplicate-contigs && -s $O4/derep-ctg/MGS.0/lift.parquet && -z $(ls -d "$O4"/derep-ctg/MGS=* 2> /dev/null) ]]'

echo
echo "$NP passed, $NF failed"
(( NF == 0 ))
