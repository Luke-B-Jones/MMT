#!/usr/bin/env python3
"""check_ctgmeth.py MOCKDIR OUT - compare `mmt ctg-meth` results on test/ctgmeth_mock.py with the truth and check every
table against the assemblies themselves. Prints PASS / FAIL lines ; exit 1 on any FAIL."""
import collections, glob, json, os, sys
import duckdb, pyarrow.parquet as pq

M, OUT = sys.argv[1:3]
T = json.load(open(f"{M}/truth.json"))
D, DR = f"{OUT}/ctg-meth", f"{OUT}/derep-ctg"
nf = 0
def t(name, ok, detail=""):
    global nf
    nf += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}{'' if ok else '   ' + str(detail)[:400]}")
con = duckdb.connect()
q = lambda sql: con.execute(sql).fetchall()
COMP = str.maketrans("ACGT", "TGCA")
CTG = T["contigs"]
TOL = 20

# ---------------------------------------------------------------- files, schemas, columns
want = {("MGS.0", s) for s in ("S1", "S2", "S3", "S4", "S5")} | {("MGS.1", "S1"), ("MGS.1", "S2")}
got = {(f.split("/")[-2], os.path.basename(f)[:-8]) for f in glob.glob(f"{D}/MGS.*/*.parquet")}
t("sites files: one per MGS x sample with a MAG", got == want, sorted(got ^ want))
for kind in ("", "junctions/", "edges/", "reads/"):
    fs = sorted(glob.glob(f"{D}/MGS.*/{kind}*.parquet"))
    sch = {str(pq.read_schema(f)) for f in fs}
    t(f"{kind or 'sites/'}: {len(fs)} files, one schema", len(fs) == len(want) and len(sch) == 1, len(sch))
bad = q(f"""SELECT count(*) FROM read_parquet('{D}/MGS.*/*.parquet', filename = true)
           WHERE filename NOT LIKE '%/' || MGS || '/' || sample || '.parquet'""")[0][0]
t("sites: MGS and sample columns name the file", bad == 0, bad)
S = q(f"SELECT MGS, sample, read_set, side, dist, MAG, contig, pos, strand, base, coverage, ipdRatio, block, instance, block_pos, "
      f"block_strand, ctx_dist, gene FROM read_parquet('{D}/MGS.*/*.parquet')")
t("sites: rows", len(S) > 100000, len(S))

# ---------------------------------------------------------------- every row = its contig's base ; on its block's base
mag_of = {}
for f in glob.glob(f"{DR}/MGS.*/instances.parquet"):
    for m, c in q(f"SELECT DISTINCT MAG, contig FROM '{f}'"): mag_of[c] = m
nb = nlift = nbad_lift = 0
bseq = {}
for f in glob.glob(f"{DR}/MGS.*/blocks.parquet"):
    for b, s in q(f"SELECT block, seq FROM '{f}'"): bseq[(f.split("/")[-2], b)] = s
for mgs, smp, rs, side, dist, mag, c, pos, st, base, cov, ipd, blk, ins, bp, bst, cd, g in S:
    cb = CTG[c][pos - 1] if st == "+" else CTG[c][pos - 1].translate(COMP)
    nb += cb != base or base not in "AC" or cov < 3
    if blk is not None:
        nlift += 1
        b = bseq[(mgs, blk)][bp - 1]
        nbad_lift += (b if bst == "+" else b.translate(COMP)) != base
t("sites: every row is an A or C of its contig on its strand, coverage >= 3", nb == 0, nb)
t(f"sites: every lifted row ({nlift}) has its base at block_pos on block_strand of blocks.parquet seq", nlift > 50000 and nbad_lift == 0, nbad_lift)
junk = [r for r in S if r[6].endswith("__C99")]
t("sites: no row on the unbinned contigs", not junk, len(junk))
t("sites: MGS.1 rows only on the Y contigs, MGS.0 never there",
  all(r[6].endswith("__C9") for r in S if r[0] == "MGS.1") and not any(r[6].endswith("__C9") for r in S if r[0] == "MGS.0"))
# ctx_dist = distance to the nearest lift segment end (spot check against lift.parquet)
seg = collections.defaultdict(list)
for f in glob.glob(f"{DR}/MGS.*/lift.parquet"):
    for c, a, b, i in q(f"SELECT contig, c_start, c_end, instance FROM '{f}'"): seg[c].append((a, b, i))
chk = bad = 0
for r in S[::97]:
    c, pos, cd, ins = r[6], r[7], r[16], r[13]
    hit = [(a, b, i) for a, b, i in seg[c] if a <= pos <= b]
    if r[12] is None: bad += bool(hit); continue
    chk += 1
    bad += not (len(hit) == 1 and hit[0][2] == ins and cd == min(pos - hit[0][0], hit[0][1] - pos, 255))
t(f"sites: block / instance / ctx_dist agree with lift.parquet ({chk} rows checked)", chk > 300 and bad == 0, bad)

# ---------------------------------------------------------------- methylation truth (read set 'all', away from events)
def motif(c):
    s, x = set(), CTG[c]
    i = x.find("GATC")
    while i >= 0: s.add((i + 2, "+")); s.add((i + 3, "-")); i = x.find("GATC", i + 1)
    return s
MOT = {c: motif(c) for c in CTG}
near = collections.defaultdict(list)
for e in T["events"]:
    for p in e["bp"]: near[(e["sample"], e["contig"])].append(p)
for smp in T["samples"]:
    on = [r[11] for r in S if r[1] == smp and r[2] == "all" and r[9] == "A" and r[10] >= 5 and (r[7], r[8]) in MOT[r[6]]
          and all(abs(r[7] - p) > 8000 for p in near[(smp, r[6])])]
    off = [r[11] for r in S if r[1] == smp and r[2] == "all" and r[9] == "A" and r[10] >= 5 and (r[7], r[8]) not in MOT[r[6]]]
    mo, mf = sum(on) / max(1, len(on)), sum(off) / max(1, len(off))
    t(f"{smp}: 'all' GATC A ipdRatio {mo:.2f} ({len(on)} sites) vs other A {mf:.2f}", len(on) > 30 and mo > 2 and mf < 1.4)

# ---------------------------------------------------------------- read sets = the planted events
J = q(f"""SELECT sample, MGS, MAG, read_set, side, kind, event, contig, pos, bp_side, partner_contig, partner_pos, n_molecules,
                 locus_molecules, frac, called, call_bp, n_used, n_rows, flags, molecules, anchor_block, anchor_block_pos, locus_ref,
                 locus, anchor_instance
          FROM read_parquet('{D}/MGS.*/junctions/*.parquet')""")
t("junctions: REF rows have locus_ref = their molecules, ALT rows locus_ref <= locus_molecules - n_molecules",
  all(r[23] == r[12] for r in J if r[5] == "ref") and all(r[23] is None or r[23] <= r[13] - r[12] for r in J if r[5] == "alt"))
JS = collections.defaultdict(list)
for r in J: JS[(r[0], r[3])].append(r)
def sets_at(smp, c, p, kind=None, called=None):
    out = set()
    for r in J:
        if r[0] == smp and r[7] == c and abs(r[8] - (p if r[9] == "R" else p + 1)) <= TOL and (kind is None or r[5] == kind) \
                and (called is None or r[15] == called):
            out.add(r[3])
    return out
def ipd(smp, mgs, rsets, c):
    v = [r[11] for r in S if r[1] == smp and r[0] == mgs and r[2] in rsets and r[9] == "A" and (r[7], r[8]) in MOT[c] and r[4] >= 30]
    return sum(v) / max(1, len(v)), len(v)
planted = set()
for e in T["events"]:
    smp, c, k = e["sample"], e["contig"], e["kind"]
    if k in ("inversion", "deletion", "misassembly"):
        alt = [sets_at(smp, c, p, "alt") for p in e["bp"]]
        ref = [sets_at(smp, c, p, "ref", True) for p in e["bp"]]
        a_all = set().union(*alt); planted |= {(smp, x) for x in a_all}
        ev = {r[6] for x in a_all for r in JS[(smp, x)]}
        want_ev = {"inversion"} if k != "deletion" else {"deletion"}
        t(f"{smp} {k} at {c}:{e['bp']}: ALT read set(s) at both breakpoints, event {sorted(ev)}", all(alt) and ev == want_ev, (alt, ev))
        if k == "misassembly":
            fl = {f for x in a_all for r in JS[(smp, x)] for f in (r[19] or "").split(",")}
            fr = min(1 - (r[23] or 0) / r[13] for x in a_all for r in JS[(smp, x)])     # not-as-assembled share
            t(f"{smp} misassembly: {fr:.2f} of the molecules not as assembled at both breakpoints, flagged assembly_suspect, "
              f"no REF read set called", fr >= 0.9 and "assembly_suspect" in fl and not any(ref), fl)
            continue
        t(f"{smp} {k}: REF read set called at both breakpoints", all(ref), ref)
        called_alt = {x for x in a_all if any(r[15] for r in JS[(smp, x)])}
        t(f"{smp} {k}: ALT read sets called", len(called_alt) >= 1, a_all)
        if e.get("short"):
            fl = {f for x in a_all for r in JS[(smp, x)] for f in (r[19] or "").split(",")}
            t(f"{smp} short inversion (aligned through): inner side not used for methylation, sets linked",
              any(f.endswith("_not_aligned") for f in fl) and any(f.startswith("linked:") for f in fl), fl)
        if "alt_meth" in e:
            ra, na = ipd(smp, e["MGS"], called_alt, c)
            rr, nr = ipd(smp, e["MGS"], set().union(*ref), c)
            if e["alt_meth"] == 0.0:
                t(f"{smp} {k}: GATC ipdRatio of ALT reads {ra:.2f} ({na}) low, of REF reads {rr:.2f} ({nr}) high",
                  na >= 5 and nr >= 5 and ra < 1.5 and rr > 2.2)
            else:
                t(f"{smp} {k}: GATC ipdRatio of ALT reads {ra:.2f} ({na}) and REF reads {rr:.2f} ({nr}) both high",
                  na >= 5 and nr >= 5 and ra > 1.8 and rr > 2.2)
    elif k == "insertion":
        s_ = sets_at(smp, c, e["bp"][0])
        planted |= {(smp, x) for x in s_}
        ev = {r[6] for x in s_ for r in JS[(smp, x)]}
        t(f"{smp} insertion at {c}:{e['bp'][0]}: listed (as {sorted(ev)})", ev & {"insertion", "unplaced"}, ev)
    elif k in ("circular", "contig_join"):
        s_ = {x for x in sets_at(smp, c, e["bp"][0], "alt") if JS[(smp, x)][0][6] == k}
        planted |= {(smp, x) for x in s_}
        called = any(r[15] for x in s_ for r in JS[(smp, x)])
        t(f"{smp} {k} at {c}:{e['bp'][0]}: one read set, called", len(s_) == 1 and called, s_)
        if k == "contig_join":
            t(f"{smp} contig_join partner = {e['partner']}", all(r[10] in (c, e["partner"]) for x in s_ for r in JS[(smp, x)]))
for smp in ("S4", "S5"):
    s_ = {r[3] for r in J if r[0] == smp and r[6] == "circular"}
    planted |= {(smp, x) for x in s_}
false = {(r[0], r[3], r[6], r[7], r[8]) for r in J if r[5] == "alt" and r[12] >= 5 and (r[0], r[3]) not in planted}
t("no ALT read set of >= 5 molecules that was not planted", not false, sorted(false)[:5])

# ---------------------------------------------------------------- read set rows, molecules, reads table, edges
bad = [(r[0], r[3], r[4]) for r in J if r[15] and r[18] != sum(1 for x in S if x[1] == r[0] and x[2] == r[3] and x[3] == r[4])]
t("junctions: n_rows = the rows of each read set side", not bad, bad[:5])
bad = [x for x in S if x[2] != "all" and not (0 <= x[4] <= next(r[16] for r in JS[(x[1], x[2])] if r[4] == x[3]))]
t("sites: read-set rows lie within their callable distance", not bad, len(bad))
t("junctions: called sides have >= 5 molecules used and a callable distance",
  all(r[17] >= 5 and r[16] > 0 for r in J if r[15]))
R_ = q(f"SELECT sample, molecule, read_sets, in_all, MAG FROM read_parquet('{D}/MGS.*/reads/*.parquet')")
rmol = collections.defaultdict(set)
for smp, m, rs, ia, mag in R_:
    for x in rs or []: rmol[(smp, x)].add(m)
bad = [(r[0], r[3]) for r in J if not set(r[20]) <= rmol[(r[0], r[3])]]
t("reads table: every molecule listed for a read set carries it in read_sets", not bad, bad[:5])
# ALT and REF of one breakpoint pair on (sample, MAG, contig, locus) ; REF rows: locus = the gap after pos (R) / before (L)
refk = {(r[0], r[2], r[7], r[24]) for r in J if r[5] == "ref"}
bad = [(r[0], r[3], r[4]) for r in J if r[5] == "alt" and r[23] is not None and (r[0], r[2], r[7], r[24]) not in refk]
t("junctions: every ALT side whose breakpoint has a REF count pairs with that REF set on (sample, MAG, contig, locus)", not bad, bad[:5])
t("junctions: locus = pos (R) / pos - 1 (L) within the clustering tolerance",
  all(r[24] is not None and abs(r[24] - (r[8] if r[9] == "R" else r[8] - 1)) <= TOL for r in J if r[5] != "unplaced"))
# anchor block position lies on the anchor instance's lift (islands too)
lift = collections.defaultdict(list)
for f in glob.glob(f"{DR}/MGS.*/lift.parquet"):
    for c, a_, b_, i in q(f"SELECT contig, c_start, c_end, instance FROM '{f}'"): lift[c].append((a_, b_, i))
bad = [(r[0], r[3], r[4]) for r in J if r[22] is not None and not any(a_ <= r[8] <= b_ and i == r[25] for a_, b_, i in lift[r[7]])]
t("junctions: anchor_block_pos only where the breakpoint end lies on the anchor instance's lift", not bad, bad[:5])
E = q(f"SELECT sample, kind, contig, pos, n_span, n_event, n_runoff, n_end FROM read_parquet('{D}/MGS.*/edges/*.parquet')")
pl = [r for r in E if r[1] == "contig_end" and r[2].endswith("__C2") and r[0] in ("S1", "S4", "S5")]
t("edges: the circular plasmid's contig ends: reads run off both", pl and all(r[6] >= 5 for r in pl), pl)
inner = [r for r in E if r[1] == "instance_end" and 2000 < r[3] < len(CTG[r[2]]) - 2000]
t(f"edges: instance ends inside contigs are spanned by molecules, or face a breakpoint ({len(inner)})",
  inner and all(r[4] > 0 or r[5] > 0 for r in inner) and sum(r[4] > 0 for r in inner if r[5] == 0) >= 0.9 * sum(r[5] == 0 for r in inner),
  [r for r in inner if r[4] == 0][:5])
altmol = collections.defaultdict(set)                 # (sample, contig) -> molecules of the ALT sets with an end there
alt_of = collections.defaultdict(set)
for r in J:
    if r[5] == "alt": alt_of[(r[0], r[7])].add(r[3])
for (smp, c), rs in alt_of.items():
    for x in rs: altmol[(smp, c)] |= rmol[(smp, x)]
bad = [r[:6] for r in E if r[1] == "instance_end" and r[5] > len(altmol[(r[0], r[2])])]
t("edges: n_event counts each molecule once", not bad, bad[:5])
st = {}
for f in glob.glob(f"{D}/stats/*.tsv"):
    if f.endswith(".run.tsv"): continue
    for l in open(f).read().splitlines()[1:]:
        x = l.split("\t"); st[(x[0], os.path.basename(f)[:-4])] = x[1]
t("state: every result done", st == {k: "done" for k in want}, st)
t("co-assembly: S4 and S5 rows on the same MAG", {r[5] for r in S if r[1] == "S4" and r[0] == "MGS.0"} == {r[5] for r in S if r[1] == "S5" and r[0] == "MGS.0"} == {"S5M2__SemiBin_1"})
print(f"\n{'all passed' if nf == 0 else f'{nf} FAILED'}")
sys.exit(nf > 0)
