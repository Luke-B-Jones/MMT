#!/usr/bin/env python3
"""check_derep.py MOCKDIR OUT/derep-ctg [CORE_PREV] - compare `mmt derep-ctg` results with test/derep_mock.py truth, and
check the base-level frame (lift) against the assemblies themselves. Prints PASS/FAIL lines, exit 1 on any FAIL."""
import bisect, collections, gzip, json, os, sys
import pandas as pd

M, D = sys.argv[1], sys.argv[2]
T = json.load(open(f"{M}/truth.json"))
nf = 0
def t(name, ok, detail=""):
    global nf
    nf += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}{'' if ok else '   ' + str(detail)}")

def blocks(m):
    bg = pd.read_parquet(f"{D}/{m}/block_genes.parquet")
    b = pd.read_parquet(f"{D}/{m}/blocks.parquet")
    fams = {k: list(v.sort_values("slot").family) for k, v in bg.groupby("block")}
    return b, fams, bg

def same_order(found, want):   # want (truth family ids, after allele merge) equals found, up to reversal
    return found == want or found == want[::-1]

fam = pd.read_parquet(f"{D}/MGS.0/families.parquet")
f_of = {int(c): int(r.family) for _, r in fam.iterrows() for c in r.catalogue}
canon = lambda xs: [f_of.get(x, x) for x in xs]
padj = dict(zip(fam.family, fam.prevalence_adj))
CORE = float(sys.argv[3]) if len(sys.argv) > 3 else 0.7   # DEREP_CORE_PREV of the run: rarer genes are bridged, not members
core_of = lambda xs: [x for x in canon(xs) if padj.get(x, 0) >= CORE]
b0, f0, bg0 = blocks("MGS.0")
core = {k: v for k, v in f0.items() if b0.set_index("block").tier[k] == "core"}
acc = {k: v for k, v in f0.items() if b0.set_index("block").tier[k] == "accessory"}
for name, want in T["segments"]["core"].items():
    w = core_of(want)
    hit = [k for k, v in core.items() if same_order(v, w)]
    t(f"MGS.0 core segment {name} ({len(w)} of its {len(want)} genes >= core prevalence) is exactly one core block", len(hit) == 1,
      [(k, len(v), len(set(v) & set(w))) for k, v in core.items() if set(v) & set(w)])
t("MGS.0 has exactly 7 core blocks", len(core) == 7, len(core))
for name, want in T["segments"]["accessory"].items():
    w = canon(want)
    hit = [k for k, v in acc.items() if same_order(v, w)]
    t(f"MGS.0 accessory {name} ({len(w)} genes) is exactly one accessory block", len(hit) == 1,
      [(k, len(v), len(set(v) & set(w))) for k, v in f0.items() if set(v) & set(w)])
t("MGS.0 has exactly 4 accessory blocks", len(acc) == 4, len(acc))
a, b = T["allele"]
t("divergent allele (own catalogue gene) merged into one family", f_of[a] == f_of[b], (f_of.get(a), f_of.get(b)))
x, y = T["not_merged"]
t("replacement gene R not merged with the gene it replaces", f_of[x] != f_of[y])
cls = fam.set_index("family").cls
t("transposase (3 copies per MAG) is a linker", cls[f_of[T["linker"]]] == "linker", cls.get(f_of[T["linker"]]))
t("contaminant genes (one MAG) are rare", all(cls[f_of[c]] == "rare" for c in T["rare"] if c in f_of))
t("replaced / replacing genes stay out of blocks", all(pd.isna(fam.set_index("family").block[f_of[c]]) for c in T["not_merged"]))
# instances: every MAG whose contig holds a whole core segment (all genes complete, no deletion, no big insert)
# must have one complete instance with exactly the first-gene-start .. last-gene-end span
inst = pd.read_parquet(f"{D}/MGS.0/instances.parquet")
genes = pd.read_parquet(f"{D}/MGS.0/genes.parquet")
seg_names = {}
for name, want in T["segments"]["core"].items():
    for c in want: seg_names[c] = name
cat_name = {v: k for k, v in T["cats"].items()}
n_exp = n_ok = 0; bad = []
bid_of = {name: [k for k, v in core.items() if same_order(v, core_of(w))][0] for name, w in T["segments"]["core"].items()
          if [k for k, v in core.items() if same_order(v, core_of(w))]}
for mag, rows in T["genes"].items():
    if T["mags"][mag]["MGS"] != "MGS.0": continue
    fl = T["mags"][mag]["flags"]
    for name, want in T["segments"]["core"].items():
        if name not in bid_of: continue
        wn = {cat_name[c].split("|")[0] for c in want if padj.get(f_of.get(c, c), 0) >= CORE}
        on = [r for r in rows if r[0].rstrip("*") in wn]
        if len(on) != len(wn) or len({r[1] for r in on}) != 1 or any(r[4] != "00" for r in on): continue
        if fl.get("phage") and ((name == "S1" and fl["phage"] == "A") or (name == "S6" and fl["phage"] == "B") or (name == "S3" and fl["phage"] == "C")): continue
        n_exp += 1
        s0, e0 = min(r[2] for r in on), max(r[3] for r in on)
        hit = inst[(inst.MAG == mag) & (inst.block == bid_of[name]) & (inst.contig == on[0][1])]
        ok = len(hit) == 1 and bool(hit.complete.iloc[0]) and hit.c_start.iloc[0] == s0 and hit.c_end.iloc[0] == e0
        n_ok += ok
        if not ok: bad.append((mag, name, s0, e0, hit[["c_start", "c_end", "complete", "n_present"]].values.tolist()))
t(f"complete core instances found with exact spans ({n_ok}/{n_exp})", n_exp > 50 and n_ok == n_exp, bad[:3])
# phage strains: the host block is split in two around the 30-kb prophage
for mag, rows in T["genes"].items():
    fl = T["mags"][mag].get("flags") or {}
    if fl.get("phage"):
        host = {"A": "S1", "B": "S6", "C": "S3"}[fl["phage"]]
        n = ((inst.MAG == mag) & (inst.block == bid_of.get(host, "?"))).sum()
        t(f"{mag}: host block {host} split around its prophage (2 instances)", n >= 2, n)
pb = [k for k, v in acc.items() if same_order(v, canon(T["segments"]["accessory"]["P"]))]
if pb:
    pi = inst[inst.block == pb[0]]
    t("prophage block: 3 MAGs, nested in nothing (between two host pieces)", pi.MAG.nunique() == 3 and (pi.host_instance < 0).all(), pi[["MAG", "host_instance"]].values.tolist())
    t("prophage families are shared with MGS.1 (shared_mgs = 2)", b0.set_index("block").max_mgs[pb[0]] == 2)
# I1 (6 kb, below max insert) sits inside S2 instances: nested, slot = 12 (after the 13th gene) or 12 counted from the other end
i1 = [k for k, v in acc.items() if same_order(v, canon(T["segments"]["accessory"]["I1"]))]
if i1:
    ii = inst[inst.block == i1[0]]
    hosts = inst.set_index("instance").block.reindex(ii.host_instance[ii.host_instance >= 0]).tolist()
    t("I1 instances nested in S2 instances", len(hosts) >= 5 and set(hosts) == {bid_of.get("S2")}, hosts)
    hp = set(ii.host_slot[ii.host_instance >= 0])
    nS2 = len(core[bid_of["S2"]])
    t("I1 insertion slot is consistent (after S2 position 12 or its mirror)", hp <= {12, nS2 - 14}, hp)
# links: inversion of S3 and translocation of S5 show as alternative block-end adjacencies
L = pd.read_parquet(f"{D}/MGS.0/links.parquet")
pairs = {frozenset((r.block_a, r.block_b)): r.n_mags for r in L[L.terminal].itertuples()}
S = collections.defaultdict(lambda: "?", bid_of)
t("S2-S3 adjacency seen (normal + inverted)", frozenset((S["S2"], S["S3"])) in pairs)
t("S1-S5 adjacency (translocation) seen", frozenset((S["S1"], S["S5"])) in pairs)
two = L[L.terminal & (((L.block_a == S["S2"]) & (L.block_b == S["S3"])) | ((L.block_a == S["S3"]) & (L.block_b == S["S2"])))]
t("S2-S3 joined at two different S3 ends (inversion)", two.shape[0] >= 2, two[["block_a", "end_a", "block_b", "end_b", "n_mags"]].values.tolist())
# block sequences: one per block, = blocks.parquet seq ; polished = the representative instance's sequence with
# substitutions only (n_polished of them)
fa = {}
with gzip.open(f"{D}/MGS.0/blocks.fna.gz", "rt") as fh:
    for l in fh:
        if l.startswith(">"): k = l[1:].split()[0]; fa[k] = []
        else: fa[k].append(l.strip())
t("blocks.fna.gz: one sequence per block, lengths as in blocks.parquet, = blocks.parquet seq",
  {k: len("".join(v)) for k, v in fa.items()} == dict(zip(b0.block, b0.length)) and {k: "".join(v) for k, v in fa.items()} == dict(zip(b0.block, b0.seq)))
t("block coverage of MAG bp >= 0.95 (MGS.0)", pd.read_parquet(f"{D}/MGS.0/mags.parquet").bp_blocks.sum() / pd.read_parquet(f"{D}/MGS.0/mags.parquet").bp.sum() >= 0.95)
t("median instance identity close to 1 - strain SNP rate", 0.97 < inst.identity.median() < 0.999, inst.identity.median())
t("genes table: every gene once, instance set for every block gene",
  genes.gene.is_unique and genes[genes.block.notna()].instance.notna().all())
# MGS.1 and MGS.2
b1, f1, _ = blocks("MGS.1")
t("MGS.1: 2 core blocks", (b1.tier == "core").sum() == 2, b1[["block", "tier", "n_genes"]].values.tolist())
b2, f2, _ = blocks("MGS.2")
t("MGS.2 (one MAG): every block from its two contigs, all instances complete",
  len(b2) >= 1 and pd.read_parquet(f"{D}/MGS.2/instances.parquet").complete.all())

# ---------------------------------------------------------------- base-level frame (lift) against the assemblies
COMP = str.maketrans("ACGTN", "TGCAN")
rc = lambda q: q.translate(COMP)[::-1]
idx = os.environ.get("MMT_TEST_IDX") or os.path.join(os.path.dirname(os.path.abspath(D)), "index")   # MMT >= 0.5: a cache in the work space
asm = {r.contig: r.assembly for r in pd.read_csv(f"{idx}/contigs.tsv.gz", sep="\t").itertuples()}
seqs = {}
def contig(c):
    if c not in seqs:
        a = asm[c]; f = next(p for p in (f"{a}/scaffolds.fasta.filt", f"{a}/scaffolds.fasta.filt.gz") if os.path.exists(p))
        cur = None
        for l in (gzip.open(f, "rt") if f.endswith(".gz") else open(f)):
            if l[0] == ">": cur = l[1:].split()[0]; seqs[cur] = []
            else: seqs[cur].append(l.strip())
        for k in list(seqs):
            if isinstance(seqs[k], list): seqs[k] = "".join(seqs[k]).upper()
    return seqs[c]
for m in ("MGS.0", "MGS.1", "MGS.2"):
    lf = pd.read_parquet(f"{D}/{m}/lift.parquet"); bl = pd.read_parquet(f"{D}/{m}/blocks.parquet").set_index("block")
    ins = pd.read_parquet(f"{D}/{m}/instances.parquet").set_index("instance")
    seq = bl.seq.to_dict()
    bad = sum((contig(r.contig)[r.c_start - 1:r.c_end] if r.ori > 0 else rc(contig(r.contig)[r.c_start - 1:r.c_end]))
              != seq[r.block][r.b_start - 1:r.b_start - 1 + r.c_end - r.c_start + 1] for r in lf.itertuples())
    t(f"{m} lift: every base of all {len(lf)} segments equals its block base (both orientations: "
      f"{(lf.ori > 0).sum()} +, {(lf.ori < 0).sum()} -)", bad == 0 and len(lf) > 0, bad)
    ov = 0
    for c, g in lf.sort_values(["contig", "c_start"]).groupby("contig"):
        ov += int((g.c_start.values[1:] <= g.c_end.values[:-1]).sum())
    t(f"{m} lift: segments never overlap on a contig", ov == 0, ov)
    dup = 0
    for i, g in lf.groupby("instance"):
        iv = sorted(zip(g.b_start, g.b_start + g.c_end - g.c_start))
        dup += sum(b <= a_end for (_, a_end), (b, _) in zip(iv, iv[1:]))
    t(f"{m} lift: a block base maps to at most one base of each instance", dup == 0, dup)
    inside = (lf.join(ins[["contig", "c_start", "c_end"]], on="instance", rsuffix="_i")
              .eval("c_start >= c_start_i and c_end <= c_end_i and contig == contig_i")).all()
    t(f"{m} lift: segments lie inside their instance", bool(inside))
    reps = ins[ins.is_rep]
    npol_ok = all(len(seq[r.block]) == r.c_end - r.c_start + 1 and
                  sum(a != b for a, b in zip(seq[r.block], (lambda q: q if r.strand == "+" else rc(q))(contig(r.contig)[r.c_start - 1:r.c_end])))
                  == bl.n_polished[r.block] for r in reps.itertuples())
    t(f"{m} polished block = representative instance + n_polished substitutions", npol_ok)
    # gene frame: each block gene copy's lifted bases fall inside its block_genes seq_start..seq_end
    bgm = pd.read_parquet(f"{D}/{m}/block_genes.parquet").set_index(["block", "slot"])
    gn = pd.read_parquet(f"{D}/{m}/genes.parquet")
    gn = gn[gn.block.notna() & gn.slot.notna()]
    seg_by = {c: (g.c_start.values, g) for c, g in lf.sort_values("c_start").groupby("contig")}
    tot = ok_ = 0
    for r in gn.itertuples():
        w = bgm.loc[(r.block, int(r.slot))]
        if pd.isna(w.seq_start) or r.contig not in seg_by: continue
        st, g = seg_by[r.contig]
        k = max(0, bisect.bisect_right(st, r.c_start) - 1)
        for x in g.iloc[k:].itertuples():
            if x.c_start > r.c_end: break
            a_, b_ = max(x.c_start, r.c_start), min(x.c_end, r.c_end)
            if a_ > b_ or x.block != r.block: continue
            for p in range(a_, b_ + 1, 7):
                bp = x.b_start + (p - x.c_start) if x.ori > 0 else x.b_start + (x.c_end - p)
                tot += 1; ok_ += w.seq_start <= bp <= w.seq_end
    t(f"{m} gene frame: {ok_}/{tot} lifted bases of block-gene copies fall inside their gene's seq_start..seq_end",
      tot > 0 and ok_ >= 0.95 * tot, (ok_, tot))
    cp = pd.read_parquet(f"{D}/{m}/composition.parquet")
    cl = pd.read_parquet(f"{D}/{m}/contigs.parquet").set_index("contig").length
    cs = pd.read_parquet(f"{D}/{m}/contig_seqs.parquet")
    t(f"{m} contig_seqs: every MAG contig once, its assembly sequence and length",
      sorted(cs.contig) == sorted(cl.index) and all(r.seq == contig(r.contig) and r.length == len(r.seq) == cl[r.contig] for r in cs.itertuples()))
    sums = cp.groupby("contig").apply(lambda g: (int((g.a + g.c + g.g + g.t).sum()), int(g.w_end.max())), include_groups=False)
    t(f"{m} composition: every MAG contig, counts add up to its length",
      set(sums.index) == set(cl.index) and all(v == (cl[c], cl[c]) for c, v in sums.items()))
    nbad = sum(1 for r in ins.itertuples() for k, i in ((r.nb5_key, r.nb5_instance), (r.nb3_key, r.nb3_instance))
               if pd.notna(i) and not k.startswith("host:") and not k.startswith(ins.block[int(i)] + ":"))
    t(f"{m} instances: neighbour keys name the neighbour's block", nbad == 0, nbad)
    t(f"{m} instances: n_mag_copies / dup_in_mag consistent",
      bool((ins.groupby(["block", "MAG"]).size().rename("n").to_frame().join(ins.set_index(["block", "MAG"]).n_mag_copies)
            .eval("n == n_mag_copies")).all()) and not ins[ins.n_mag_copies == 1].dup_in_mag.any())
print(f"\n{'all passed' if nf == 0 else f'{nf} FAILED'}")
sys.exit(nf > 0)
