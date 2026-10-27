#!/usr/bin/env python3
"""Synthetic MATAFILER4 run + GeneCat gene catalogue with a known block structure, for `mmt dereplicate-contigs`.

usage: derep_mock.py OUTDIR [SEED]          OUTDIR is wiped ; truth in OUTDIR/truth.json ; needs samtools on PATH

Species X = MGS.0, 20 MAGs (strains) from 18 samples (two strains binned in one sample P05; a co-assembly G1 of
P17 + P18 holding two strains). Backbone S1 S2 S3 S4 S5 S6 (145 genes). Variation:
  S3 inverted (30 %) ; S5 moved between S1 and S2 (20 %)                         -> block boundaries
  I1  6-gene island inside S2 (40 %, ~6 kb  <= max insert)                       -> S2 stays whole, I1 accessory block
  I2 20-gene island inside S4 (35 %, ~22 kb > max insert)                        -> S4 splits in S4a | S4b, I2 accessory
  P  25-gene prophage at two different places (2 strains each) and in 2 MAGs of species Y
  PL 12-gene plasmid contig (50 %) ; S2 gene 5: a 9 %-divergent allele (own catalogue gene, 50 %) -> merged
  S6 gene 12 replaced by an unrelated gene R (50 %) -> not merged, both bridged ; S2 gene 20 absent in 25 % and S5
  genes 3-4 absent in 25 % -> blocks stay whole ; transposase T, 3 copies per MAG ->
  linker ; single-gene deletions ; 0.3-1.2 % strain SNPs ; contigs broken at random (genes cut at contig ends become
  prodigal partial genes), some contigs reverse-complemented, one contig dropped in some MAGs ; a contaminant contig
  (species W) in one MAG ; one unbinned contig per assembly. Species Y = MGS.1 (4 MAGs), Z = MGS.2 (1 MAG).

usage: derep_mock.py --mix OUTDIR DEST DONOR...   the MGS.0 MAG of sample DEST also gets copies of the MGS.0 contigs
                                                  (sequence, genes, bins, catalogue members) of the DONOR samples:
                                                  a MAG of co-binned strains
"""
import gzip, json, os, random, shutil, sys


def mix(H, dest, donors):
    run = f"{H}/mf4/mock"
    asm = lambda s: open(f"{run}/{s}/assemblies/metag/assembly.txt").read().strip()
    def fasta(p):
        d, k = {}, None
        for l in (gzip.open(p, "rt") if p.endswith(".gz") else open(p)):
            if l[0] == ">": k = l[1:].split()[0]; d[k] = []
            else: d[k].append(l.strip())
        return {k: "".join(v) for k, v in d.items()}
    dst, new, fa, gff, fna, bins = asm(dest), {}, [], [], [], []
    for s in donors:
        a = asm(s); seq = fasta(f"{a}/scaffolds.fasta.filt")
        for l in open(f"{a}/Binning/SB/{s}").read().splitlines()[1:]:
            c, b = l.split("\t")
            if b == "SemiBin_1":
                new[c] = f"{dest}__C{900 + len(new)}_L={len(seq[c])}="
                fa.append(f">{new[c]}\n{seq[c]}\n"); bins.append(f"{new[c]}\tSemiBin_1\n")
        gff += ["\t".join([new[x[0]]] + x[1:]) for x in (l.split("\t") for l in gzip.open(f"{a}/genePred/genes.gff.gz", "rt")) if x[0] in new]
        fna += [f">{new[c]}_{n}\n{q}\n" for g, q in fasta(f"{a}/genePred/genes.shrtHD.fna.gz").items()
                for c, _, n in [g.rpartition("_")] if c in new]
    open(f"{dst}/scaffolds.fasta.filt", "a").write("".join(fa)); open(f"{dst}/Binning/SB/{dest}", "a").write("".join(bins))
    for f, add in (("genes.gff.gz", gff), ("genes.shrtHD.fna.gz", fna)):
        t = gzip.open(f"{dst}/genePred/{f}", "rt").read()
        with gzip.open(f"{dst}/genePred/{f}", "wt") as o: o.write(t + "".join(add))
    cidx = f"{H}/gc/compl.incompl.95.fna.clstr.idx"
    out = []
    for l in open(cidx):
        l = l.rstrip("\n")
        out.append(l + "".join(f",>{new[c]}_{n}" for m in l.partition("\t")[2].split(",") for c, _, n in [m[1:].rpartition("_")] if c in new) + "\n")
    open(cidx, "w").writelines(out)


if sys.argv[1] == "--mix": mix(os.path.abspath(sys.argv[2]), sys.argv[3], sys.argv[4:]); sys.exit()
H = os.path.abspath(sys.argv[1]); SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 5
shutil.rmtree(H, ignore_errors=True); os.makedirs(H)
R = random.Random(SEED)
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]
GC = {a + b + c: aa for (a, b, c), aa in zip(((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
      "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
def translate(s): return "".join(GC.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
def gene(): return "ATG" + "".join(R.choice(CODONS) for _ in range(R.randint(200, 480))) + "TAA"
def mut(s, r): return "".join(R.choice("ACGT".replace(b, "")) if R.random() < r else b for b in s)
def spacer(): return "".join(R.choice("ACGT") for _ in range(R.randint(80, 300)))

# ---------------------------------------------------------------- gene families (truth names) with base sequences
FAM = {}                                         # name -> dict(seq, strand)
def fam(name, strand=None): FAM[name] = dict(seq=gene(), strand=strand or R.choice((1, -1))); return name
SEG = {k: [fam(f"{k}g{i}") for i in range(n)] for k, n in (("S1", 30), ("S2", 25), ("S3", 20), ("S4", 30), ("S5", 15), ("S6", 25))}
I1 = [fam(f"I1g{i}") for i in range(6)]; I2 = [fam(f"I2g{i}") for i in range(20)]; P = [fam(f"Pg{i}") for i in range(25)]
PL = [fam(f"PLg{i}") for i in range(12)]; fam("R"); fam("T")
YS = {k: [fam(f"{k}g{i}") for i in range(n)] for k, n in (("Y1", 20), ("Y2", 15))}
ZS = [fam(f"Zg{i}") for i in range(18)]; W = [fam(f"Wg{i}") for i in range(8)]
ALLELE = mut(FAM["S2g5"]["seq"][3:-3], 0.09); ALLELE = "ATG" + ALLELE + "TAA"
SPACE = {}                                       # consensus spacer after a gene, per (gene, next gene)
def sp(a, b): return SPACE.setdefault((a, b), spacer())

# ---------------------------------------------------------------- strains of X
NX = 20
flags = [dict(inv=i < 6, trans=6 <= i < 10, I1=i % 2 == 1 and i < 16, I2=i in (0, 2, 4, 8, 10, 12, 14),
              phage={16: "A", 17: "C", 18: "B"}.get(i), PL=i % 2 == 0, allele=i % 2 == 0,
              R=i >= 10, snp=R.uniform(0.003, 0.012)) for i in range(NX)]
flags[3]["del"] = ["S4g7"]; flags[12]["del"] = ["S1g20"]
for i in (0, 1, 2, 3, 16): flags[i].setdefault("del", []).append("S2g20")          # soft-core gene: absent in 25 %
for i in (10, 11, 12, 13, 14): flags[i].setdefault("del", []).extend(["S5g3", "S5g4"])   # 2-gene deletion in 25 %


def x_layout(f):
    """ordered (gene name, strand, sequence) of strain f's chromosome"""
    order = [("S1", 1), ("S5", 1), ("S2", 1), ("S3", 1), ("S4", 1), ("S6", 1)] if f["trans"] else \
            [("S1", 1), ("S2", 1), ("S3", 1), ("S4", 1), ("S5", 1), ("S6", 1)]
    out = []
    for seg, _ in order:
        genes = list(SEG[seg])
        if seg == "S2" and f["I1"]: genes = genes[:13] + I1 + genes[13:]
        if seg == "S4" and f["I2"]: genes = genes[:16] + I2 + genes[16:]
        if seg == "S1" and f["phage"] == "A": genes = genes[:11] + P + genes[11:]
        if seg == "S6" and f["phage"] == "B": genes = genes[:6] + P + genes[6:]
        if seg == "S3" and f["phage"] == "C": genes = genes[:9] + P + genes[9:]
        if seg == "S6" and f["R"]: genes = [("R" if g == "S6g12" else g) for g in genes]
        genes = [g for g in genes if g not in f.get("del", [])]
        items = [(g + "*", FAM[g]["strand"], ALLELE) if (g == "S2g5" and f["allele"]) else (g, FAM[g]["strand"], FAM[g]["seq"]) for g in genes]
        if seg == "S3" and f["inv"]: items = [(g, -s, q) for g, s, q in reversed(items)]
        out += items
    return out


def assemble(items, name_gap=sp):
    """items -> sequence + gene coordinates (1-based, gene strand)"""
    seq, genes, pos = [], [], 0
    for k, (g, s, q) in enumerate(items):
        if k:
            gap = name_gap(items[k - 1][0], g); seq.append(gap); pos += len(gap)
        gq = q if s > 0 else rc(q)
        seq.append(gq); genes.append((g, pos + 1, pos + len(gq), s)); pos += len(gq)
    return "".join(seq), genes


def strainify(seq, genes, rate, transposons, R_):
    """SNPs (codon-safe: never create or destroy a stop inside genes ; good enough: mutate 3rd codon positions and
    spacers) + transposon copies inserted between genes"""
    s = list(seq)
    for i in range(len(s)):
        if R_.random() < rate: s[i] = R_.choice("ACGT".replace(s[i], ""))
    seq = "".join(s)
    # keep coding sequences valid: re-impose start/stop and remove new internal stops by restoring the original codon
    for g, a, b, st in genes:
        q = seq[a - 1:b]; q = q if st > 0 else rc(q)
        orig = ALLELE if g == "S2g5*" else FAM[g]["seq"] if g in FAM else q
        cod = [q[i:i + 3] for i in range(0, len(q), 3)]
        cod = [c if (c not in ("TAA", "TAG", "TGA") or i == len(cod) - 1) else (orig[3 * i:3 * i + 3] if len(orig) > 3 * i + 2 else "GCT")
               for i, c in enumerate(cod)]
        cod[0], cod[-1] = "ATG", "TAA"
        q = "".join(cod); q = q if st > 0 else rc(q)
        seq = seq[:a - 1] + q + seq[b:]
    for _ in range(transposons):                    # insert T after a random gene (between genes)
        k = R_.randrange(len(genes) - 1); a = genes[k][2]
        t = FAM["T"]["seq"] if FAM["T"]["strand"] > 0 else rc(FAM["T"]["seq"])
        ins = "ACGTTGCA" + t + "TTGCAACG"
        seq = seq[:a] + ins + seq[a:]
        genes = [(g, x + (len(ins) if x > a else 0), y + (len(ins) if x > a else 0), st) for g, x, y, st in genes]
        genes.insert(k + 1, ("T", a + 9, a + 8 + len(t), FAM["T"]["strand"]))
    return seq, genes


def break_up(seq, genes, n_breaks, R_, drop=False, min_len=4000):
    """cut into contigs ; genes cut by a contig end become partial ; some contigs reverse-complemented"""
    L = len(seq); cuts = sorted(R_.sample(range(min_len, L - min_len), n_breaks)) if n_breaks else []
    cuts = [c for i, c in enumerate(cuts) if i == 0 or c - cuts[i - 1] >= min_len]
    bounds = [0] + cuts + [L]
    ctgs = []
    for a, b in zip(bounds, bounds[1:]):
        gs = []
        for g, x, y, st in genes:
            if y < a + 1 or x > b: continue
            xa, yb = max(x, a + 1), min(y, b)
            part = ("1" if xa > x else "0") + ("1" if yb < y else "0")     # prodigal: left / right contig edge
            if yb - xa + 1 < 90: continue
            # partial genes: prodigal reports them in frame to the contig end ; keep the cut coordinates
            gs.append((g, xa - a, yb - a, st, part, (yb - xa + 1) / (y - x + 1)))
        ctgs.append([seq[a:b], gs])
    if drop and len(ctgs) > 2: ctgs.pop(R_.randrange(1, len(ctgs) - 1))
    out = []
    for q, gs in ctgs:
        if R_.random() < 0.4:                       # reverse-complement the contig
            Lq = len(q)
            gs = [(g, Lq - y + 1, Lq - x + 1, -st, part[::-1], fr) for g, x, y, st, part, fr in gs]
            q = rc(q)
        out.append((q, sorted(gs, key=lambda t: (t[1], t[2]))))
    return out


# ---------------------------------------------------------------- samples, assemblies, MAGs
RUN = os.path.join(H, "mf4", "mock")
samples = [f"P{i:02d}" for i in range(19)]                  # P00..P18
X_ASM = ["P00", "P01", "P02", "P03", "P05", "P05", "P04", "P06", "P07", "P08", "P09", "P10", "P11", "P12", "P13", "P14",
         "P15", "P16", "G1", "G1"]                       # strain -> assembly (P05: two strains ; G1: co-assembly P17 + P18)
ASM = {}                                                    # assembly id -> dict(dir, prefix, binfile_sample, contigs=[(name, seq, genes, bin)])
for s in samples[:17]:
    ASM[s] = dict(dir=f"{RUN}/{s}/assemblies/metag", prefix=s, bins=s, ctgs=[])
ASM["G1"] = dict(dir=f"{RUN}/AssmblGrp_G1/metag", prefix="P18M2", bins="P18", ctgs=[])
sample_asm = {s: s for s in samples[:17]}; sample_asm.update({"P17": "G1", "P18": "G1"})
truth = dict(MGS={}, mags={}, segments={}, cats={})
mags = []                                                   # (MAG name, MGS, assembly id, completeness, contamination)
bin_n = {}


def add_mag(asm_id, mgs, contigs_, compl, cont, tag):
    A = ASM[asm_id]; bin_n[asm_id] = bin_n.get(asm_id, 0) + 1
    b = f"SemiBin_{bin_n[asm_id]}"; mag = f"{A['prefix']}__{b}"
    for q, gs in contigs_:
        k = len(A["ctgs"]) + 1
        A["ctgs"].append((f"{A['prefix']}__C{k}_L={len(q)}=", q, gs, b))
    mags.append((mag, mgs, asm_id, compl, cont)); truth["mags"][mag] = dict(MGS=mgs, tag=tag)
    return mag


for i, f in enumerate(flags):
    items = x_layout(f)
    seq, genes = assemble(items)
    seq, genes = strainify(seq, genes, f["snp"], 3, R)
    drop = i in (2, 9, 15)
    ctgs = break_up(seq, genes, R.randint(1, 4), R, drop=drop)
    if f["PL"]:
        pq_, pg = assemble([(g, FAM[g]["strand"], FAM[g]["seq"]) for g in PL]); pq_, pg = strainify(pq_, pg, f["snp"], 0, R)
        ctgs += break_up(pq_, pg, 0, R)
    if i == 7:
        wq, wg = assemble([(g, FAM[g]["strand"], FAM[g]["seq"]) for g in W]); ctgs += break_up(wq, wg, 0, R)
    tot = sum(len(q) for q, _ in ctgs)
    mag = add_mag(X_ASM[i], "MGS.0", ctgs,
                  round(100 - (8 if drop else R.uniform(0, 3)), 2), round(R.uniform(0.2, 2.5) + (3 if i == 7 else 0), 2), f"X{i}")
    truth["mags"][mag].update(flags=f, snp=f["snp"])
for j, a in enumerate(("P00", "P01", "P02", "P03")):        # species Y, the prophage in two of them
    items = [(g, FAM[g]["strand"], FAM[g]["seq"]) for g in YS["Y1"]]
    if j < 2: items = items[:8] + [(g, FAM[g]["strand"], FAM[g]["seq"]) for g in P] + items[8:]
    items += [(g, FAM[g]["strand"], FAM[g]["seq"]) for g in YS["Y2"]]
    seq, genes = assemble(items); seq, genes = strainify(seq, genes, 0.004, 0, R)
    add_mag(a, "MGS.1", break_up(seq, genes, 1, R), 97.0, 1.0, f"Y{j}")
seq, genes = assemble([(g, FAM[g]["strand"], FAM[g]["seq"]) for g in ZS])
add_mag("P04", "MGS.2", break_up(seq, genes, 1, R), 92.0, 0.5, "Z0")
for a, A in ASM.items():                                    # one unbinned contig per assembly
    uq, ug = assemble([(f"U_{a}_{k}", 1, gene()) for k in range(5)], name_gap=lambda x, y: spacer())
    A["ctgs"].append((f"{A['prefix']}__C{len(A['ctgs']) + 1}_L={len(uq)}=", uq, [(g, x, y, s, "00", 1.0) for g, x, y, s in ug], None))

# ---------------------------------------------------------------- catalogue: cluster = (family, allele) ; partial
# pieces < 90 % of their gene = own cluster ; unbinned genes = own clusters
clusters, members, ids = {}, {}, []
def cid(key):
    if key not in clusters: clusters[key] = len(clusters) + 1; members[clusters[key]] = []
    return clusters[key]
gene_cat = {}
for a, A in sorted(ASM.items()):
    for name, q, gs, b in A["ctgs"]:
        for n, (g, x, y, st, part, fr) in enumerate(gs, 1):
            gq = q[x - 1:y] if st > 0 else rc(q[x - 1:y])
            if g.startswith("U_"): key = ("U", g)
            elif part != "00" and (y - x + 1) < 100: key = ("frag", name, n)   # GeneCat clusterCov 0.9 of the SHORTER gene: fragments join
            elif g == "S2g5*": key = ("S2g5", "allele")
            else: key = (g, "base")
            c = cid(key); members[c].append((f"{name}_{n}", gq, part)); gene_cat[(name, n)] = c
truth["cats"] = {f"{k[0]}|{k[1]}" if len(k) == 2 else f"{k[0]}|{k[1]}|{k[2]}": v for k, v in clusters.items()}
GCD = os.path.join(H, "gc"); os.makedirs(f"{GCD}/Bin_SB/Genomes/MGS_ctg", exist_ok=True)
with open(f"{GCD}/compl.incompl.95.fna", "w") as o, open(f"{GCD}/compl.incompl.95.prot.faa", "w") as p, \
        open(f"{GCD}/compl.incompl.95.fna.clstr.idx", "w") as ix, open(f"{GCD}/Matrix.genes2rows.txt", "w") as g2r:
    ix.write("#Gene\tmembers\n"); g2r.write("#idx\tcluster\trep\n")
    for c in sorted(members, key=lambda c: R.random()):     # file order is not catalogue order (as GeneCat)
        ms = members[c]; rep = max(ms, key=lambda m: (m[2] == "00", len(m[1])))
        o.write(f">{c}\n{rep[1]}\n"); p.write(f">{c}\n{translate(rep[1])}\n")
        ms = ms[:]; R.shuffle(ms)
        ix.write(f"{c}\t" + ",".join(">" + m[0] for m in ms) + "\n"); g2r.write(f"{c}\t>cl{c}\t>{rep[0]}\n")

# ---------------------------------------------------------------- MATAFILER4 assemblies, gene predictions, bins
for a, A in ASM.items():
    d = A["dir"]; os.makedirs(f"{d}/genePred", exist_ok=True); os.makedirs(f"{d}/Binning/SB", exist_ok=True)
    with open(f"{d}/scaffolds.fasta.filt", "w") as o:
        for name, q, _, _ in A["ctgs"]: o.write(f">{name}\n" + "\n".join(q[i:i + 100] for i in range(0, len(q), 100)) + "\n")
    with gzip.open(f"{d}/genePred/genes.gff.gz", "wt") as gf, gzip.open(f"{d}/genePred/genes.shrtHD.fna.gz", "wt") as fn:
        gf.write("##gff-version  3\n")
        for k, (name, q, gs, _) in enumerate(A["ctgs"], 1):
            gf.write(f"# Sequence Data: seqnum={k};seqlen={len(q)};seqhdr=\"{name}\"\n# Model Data: version=Prodigal.v2.6.3;run_type=Metagenomic;\n")
            for n, (g, x, y, st, part, fr) in enumerate(gs, 1):
                gf.write(f"{name}\tProdigal_v2.6.3\tCDS\t{x}\t{y}\t99.9\t{'+' if st > 0 else '-'}\t0\tID={k}_{n};partial={part};start_type=ATG;"
                         f"rbs_motif=None;rbs_spacer=None;gc_cont=0.500;conf=99.99;score=99.9;cscore=90.0;sscore=9.9;rscore=0.0;uscore=0.0;tscore=0.0;\n")
                gq = q[x - 1:y] if st > 0 else rc(q[x - 1:y])
                fn.write(f">{name}_{n}\n{gq}\n")
    with open(f"{d}/Binning/SB/{A['bins']}", "w") as o:
        o.write("Sequence ID\tBin\n" + "".join(f"{name}\t{b}\n" for name, _, _, b in A["ctgs"] if b))
for s in samples:
    sd = f"{RUN}/{s}"; os.makedirs(f"{sd}/assemblies/metag", exist_ok=True); os.makedirs(f"{sd}/mapping", exist_ok=True)
    open(f"{sd}/assemblies/metag/assembly.txt", "w").write(ASM[sample_asm[s]]["dir"] + "/\n")
    open(f"{sd}/mapping/{s}-smd.cram", "w").close()
raw = os.path.join(H, "raw"); os.makedirs(raw)
for s in samples: open(f"{raw}/m84_{s}.hifi_reads.bam", "w").close()
with open(f"{H}/map.txt", "w") as o:
    o.write(f"#OutPath\t{H}/mf4/\n#RunID\tmock\n#DirPath\t{raw}/\n#SmplID\tSmplPrefix\tAssmblGrps\n")
    for s in samples: o.write(f"{s}\tm84_{s}.hifi_reads\t{'G1' if s in ('P17', 'P18') else ''}\n")
# MAGvsGC: rows grouped by MGS, representative starred ; GeneCat representative FASTAs
rep_of = {"MGS.0": [m for m, mg, *_ in mags if truth["mags"][m]["tag"] == "X13"][0], "MGS.1": [m for m, mg, *_ in mags if mg == "MGS.1"][0],
          "MGS.2": [m for m, mg, *_ in mags if mg == "MGS.2"][0]}
with gzip.open(f"{GCD}/Bin_SB/MAGvsGC.txt.gz", "wt") as o:
    o.write("MAG\tMGS\tRepresentative4MGS\tMatch2MGS\tUniqueness\tAssociatedMGS\tCompleteness\tContamination\tN50\n")
    for mgs in ("MGS.0", "MGS.1", "MGS.2"):
        for m, mg, a, comp, cont in mags:
            if mg == mgs: o.write(f"{m}\t{mg}\t{'*' if rep_of[mgs] == m else ''}\t1\t1\t\t{comp}\t{cont}\t50000\n")
for mgs, m in rep_of.items():
    a = [x for x in mags if x[0] == m][0][2]; b = m.split("__")[-1]
    with gzip.open(f"{GCD}/Bin_SB/Genomes/MGS_ctg/{mgs}.ctgs.{m}.fna.gz", "wt") as o:
        for name, q, _, bb in ASM[a]["ctgs"]:
            if bb == b: o.write(f">{name}\n{q}\n")

# ---------------------------------------------------------------- truth
def fid(g): return clusters.get((g, "base"))
def fids(gs): return [x for x in map(fid, gs) if x is not None]
seg = {k: fids(v) for k, v in SEG.items()}
seg["S4a"], seg["S4b"] = seg["S4"][:16], seg["S4"][16:]
seg["S6"] = [x for x in seg["S6"] if x != fid("S6g12")]
truth["segments"] = dict(core={k: seg[k] for k in ("S1", "S2", "S3", "S4a", "S4b", "S5", "S6")},
                         accessory={"I1": fids(I1), "I2": fids(I2), "P": fids(P), "PL": fids(PL)},
                         MGS1={"Y1": fids(YS["Y1"]), "Y2": fids(YS["Y2"])}, MGS2={"Z": fids(ZS)})
truth["allele"] = [clusters[("S2g5", "base")], clusters[("S2g5", "allele")]]
truth["not_merged"] = [fid("S6g12"), fid("R")]
truth["linker"] = fid("T"); truth["rare"] = fids(W); truth["phage"] = fids(P)
truth["rep"] = rep_of
# per MAG: where each segment sits (contig, first gene start, last gene end, all genes complete, gene count)
loc = {}
for a, A in ASM.items():
    for name, q, gs, b in A["ctgs"]:
        if not b: continue
        mag = f"{A['prefix']}__{b}"
        for n, (g, x, y, st, part, fr) in enumerate(gs, 1):
            loc.setdefault(mag, []).append((g, name, x, y, part))
truth["genes"] = loc
json.dump(truth, open(f"{H}/truth.json", "w"))
print(H)
