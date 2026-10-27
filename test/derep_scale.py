#!/usr/bin/env python3
"""Large synthetic MATAFILER4 run for timing / memory tests of `mmt dereplicate-contigs` (no exact truth).

usage: derep_scale.py OUTDIR [N_MAGS=150] [N_GENES=2500] [SEED=1]

One MGS of N_MAGS strains (each ~N_GENES genes, ~1 kb each) with 8 inversion polymorphisms (10-40 %), 40 islands
(5-40 genes, 5-60 %), single-gene presence/absence, a 3-copy insertion element, 1 % SNPs, 5-30 contig breaks per MAG,
3 % of genes with a second catalogue allele ; plus one 3-MAG MGS. Writes the same MATAFILER4 / GeneCat files as
derep_mock.py (assemblies, prodigal GFF + gene FASTA, bins, catalogue index + proteins, MAGvsGC, map).
"""
import gzip, os, shutil, sys
import numpy as np

H = os.path.abspath(sys.argv[1]); NM = int(sys.argv[2]) if len(sys.argv) > 2 else 150
NG = int(sys.argv[3]) if len(sys.argv) > 3 else 2500; rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
shutil.rmtree(H, ignore_errors=True); os.makedirs(H)
B = np.frombuffer(b"ACGT", np.uint8); C = np.zeros(256, np.uint8); C[list(b"ACGT")] = list(b"TGCA")
STOPS = {b"TAA", b"TAG", b"TGA"}
codons = np.array([[a, b, c] for a in b"ACGT" for b in b"ACGT" for c in b"ACGT" if bytes([a, b, c]) not in STOPS], np.uint8)
def gene(): return np.concatenate([np.frombuffer(b"ATG", np.uint8), codons[rng.integers(0, len(codons), rng.integers(100, 500))].ravel(), np.frombuffer(b"TAA", np.uint8)])
def rc(a): return C[a[::-1]]
GCODE = {bytes(c): aa for c, aa in zip((bytes([x, y, z]) for x in b"TCAG" for y in b"TCAG" for z in b"TCAG"),
         "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
def tr(a): b = a.tobytes(); return "".join(GCODE.get(b[i:i + 3], "X") for i in range(0, len(b) - 2, 3))

fams = [dict(seq=gene(), st=int(rng.choice((1, -1)))) for _ in range(NG)]
islands = []
for n in rng.integers(5, 41, 40):                       # each island: its own new gene families
    ids = list(range(len(fams), len(fams) + int(n))); fams += [dict(seq=gene(), st=int(rng.choice((1, -1)))) for _ in ids]
    islands.append((int(rng.integers(0, NG - 1)), ids, float(rng.uniform(0.05, 0.6))))
T = len(fams); fams.append(dict(seq=gene(), st=1))
inv = [(int(a), int(a + rng.integers(20, 200)), float(rng.uniform(0.1, 0.4))) for a in rng.integers(0, NG - 250, 8)]
allele = {int(f): None for f in rng.choice(NG, int(0.03 * NG), replace=False)}
for f in allele:
    s = fams[f]["seq"].copy(); m = rng.random(len(s)) < 0.08; m[:3] = m[-3:] = False
    s[m] = B[rng.integers(0, 4, m.sum())]
    for i in range(3, len(s) - 3, 3):
        if s[i:i + 3].tobytes() in STOPS: s[i:i + 3] = fams[f]["seq"][i:i + 3]
    allele[f] = s
SPACE = {}
RUN = f"{H}/mf4/run"; GCD = f"{H}/gc"; os.makedirs(f"{GCD}/Bin_SB/Genomes/MGS_ctg")
cat = {}; members = {}
def cid(k):
    if k not in cat: cat[k] = len(cat) + 1; members[cat[k]] = []
    return cat[k]
mags, maplines = [], []
for i in range(NM + 3):
    sname = f"S{i:03d}"; d = f"{RUN}/{sname}/assemblies/metag"; os.makedirs(f"{d}/genePred"); os.makedirs(f"{d}/Binning/SB")
    os.makedirs(f"{RUN}/{sname}/mapping"); open(f"{RUN}/{sname}/mapping/{sname}-smd.cram", "w").close()
    open(f"{RUN}/{sname}/assemblies/metag/assembly.txt", "w").write(d + "/\n"); maplines.append(sname)
    small = i >= NM
    order = [(g, 1) for g in (range(60) if small else range(NG))]                   # (family, orientation)
    if not small:
        for a, b, fr in inv:                                                         # inversion: reversed, strands flipped
            if rng.random() < fr: order[a:b] = [(g, -o) for g, o in reversed(order[a:b])]
        for site, ids, fr in islands:
            if rng.random() < fr:
                k = next(j for j, (g, _) in enumerate(order) if g == site); order[k + 1:k + 1] = [(x, 1) for x in ids]
        order = [(g, o) for g, o in order if not (g < NG and rng.random() < 0.01)]   # single-gene absence
        for _ in range(3): order.insert(int(rng.integers(1, len(order))), (T, 1))
    parts, genes, pos, prev = [], [], 0, None
    for g, o in order:
        if (prev, g) not in SPACE and (g, prev) not in SPACE: SPACE[(prev, g)] = B[rng.integers(0, 4, int(rng.integers(50, 250)))]
        sp = SPACE[(prev, g)] if (prev, g) in SPACE else rc(SPACE[(g, prev)])   # shared intergenic sequence ; inverted = rc
        parts.append(sp); pos += len(sp); prev = g
        q = allele[g] if (g in allele and i % 2) else fams[g]["seq"]
        st = fams[g]["st"] * o
        parts.append(q if st > 0 else rc(q)); genes.append((g, pos + 1, pos + len(q), st, g in allele and i % 2 == 1)); pos += len(q)
    seq = np.concatenate(parts)
    m = rng.random(len(seq)) < 0.01; seq[m] = B[rng.integers(0, 4, m.sum())]
    for g, a, b, st, al in genes:                                                    # keep start / stop codons
        if st > 0: seq[a - 1:a + 2] = list(b"ATG"); seq[b - 3:b] = list(b"TAA")
        else: seq[a - 1:a + 2] = list(b"TTA"); seq[b - 3:b] = list(b"CAT")
    cuts = np.sort(rng.choice(np.arange(5000, len(seq) - 5000), 3 if small else int(rng.integers(5, 31)), replace=False))
    bounds = [0] + [int(c) for c in cuts] + [len(seq)]
    ctgs = []
    for a, b in zip(bounds, bounds[1:]):
        if b - a < 3000: continue
        gs = []
        for g, x, y, st, al in genes:
            xa, yb = max(x, a + 1), min(y, b)
            if yb - xa + 1 < 100: continue
            gs.append((g, xa - a, yb - a, st, ("1" if xa > x else "0") + ("1" if yb < y else "0"), al))
        ctgs.append((seq[a:b], gs))
    with open(f"{d}/scaffolds.fasta.filt", "wb") as fa, gzip.open(f"{d}/genePred/genes.gff.gz", "wt") as gf, \
            gzip.open(f"{d}/genePred/genes.shrtHD.fna.gz", "wb") as fn, open(f"{d}/Binning/SB/{sname}", "w") as bf:
        bf.write("Sequence ID\tBin\n")
        for k, (q, gs) in enumerate(ctgs, 1):
            name = f"{sname}__C{k}_L={len(q)}="
            fa.write(f">{name}\n".encode() + q.tobytes() + b"\n"); bf.write(f"{name}\tSemiBin_1\n")
            gf.write(f"# Sequence Data: seqnum={k};seqlen={len(q)};seqhdr=\"{name}\"\n")
            for n, (g, x, y, st, part, al) in enumerate(sorted(gs, key=lambda t: t[1]), 1):
                gf.write(f"{name}\tProdigal_v2.6.3\tCDS\t{x}\t{y}\t9\t{'+' if st > 0 else '-'}\t0\tID={k}_{n};partial={part};start_type=ATG;\n")
                gq = q[x - 1:y] if st > 0 else rc(q[x - 1:y])
                fn.write(f">{name}_{n}\n".encode() + gq.tobytes() + b"\n")
                c = cid(("s", g)) if small else cid((g, al))
                members[c].append((f"{name}_{n}", gq if part == "00" else None))
    mags.append((f"{sname}__SemiBin_1", "MGS.1" if small else "MGS.0"))
with open(f"{GCD}/compl.incompl.95.fna", "wb") as o, open(f"{GCD}/compl.incompl.95.prot.faa", "w") as p, \
        open(f"{GCD}/compl.incompl.95.fna.clstr.idx", "w") as ix:
    ix.write("#Gene\tmembers\n")
    for c, ms in members.items():
        rep = next((q for _, q in ms if q is not None), None)
        if rep is None: rep = np.frombuffer(b"ATGAAATAA", np.uint8)
        o.write(f">{c}\n".encode() + rep.tobytes() + b"\n"); p.write(f">{c}\n{tr(rep)}\n")
        ix.write(f"{c}\t" + ",".join(">" + n for n, _ in ms) + "\n")
with gzip.open(f"{GCD}/Bin_SB/MAGvsGC.txt.gz", "wt") as o:
    o.write("MAG\tMGS\tRepresentative4MGS\tCompleteness\tContamination\n")
    for k, (m, g) in enumerate(mags): o.write(f"{m}\t{g}\t{'*' if k in (0, NM) else ''}\t95\t1\n")
for k in (0, NM):
    m, g = mags[k]; s = m.split("__")[0]
    with open(f"{RUN}/{s}/assemblies/metag/scaffolds.fasta.filt", "rb") as fi, gzip.open(f"{GCD}/Bin_SB/Genomes/MGS_ctg/{g}.ctgs.{m}.fna.gz", "wb") as fo:
        fo.write(fi.read())
os.makedirs(f"{H}/raw")
for s in maplines: open(f"{H}/raw/m_{s}.hifi_reads.bam", "w").close()
with open(f"{H}/map.txt", "w") as o:
    o.write(f"#OutPath\t{H}/mf4/\n#RunID\trun\n#DirPath\t{H}/raw/\n#SmplID\tSmplPrefix\n" + "".join(f"{s}\tm_{s}.hifi_reads\n" for s in maplines))
print(H, len(cat), "catalogue genes")
