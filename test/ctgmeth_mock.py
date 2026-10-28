#!/usr/bin/env python3
"""Synthetic MATAFILER4 run + GeneCat (genes, catalogue, bins, MAGvsGC) WITH PacBio kinetics and known conformations,
for `mmt derep-ctg` + `mmt ctg-meth` end to end.

usage: ctgmeth_mock.py OUTDIR          (OUTDIR is wiped ; needs pysam, samtools, minimap2 ; truth in OUTDIR/truth.json)

Species X = MGS.0 (chromosome ~ 48 kb: segments A B C D of 10 genes ; an intergenic 1.2 kb region INV inside A, a
180-bp invertible promoter IVT with 25-bp inverted repeats inside C), species Y = MGS.1 (14 genes). Methylation: m6A at
the A of every GATC (both strands), methylated fraction 0.9 - except where a conformation says otherwise.
  S1   MAG = ancestral X + circular plasmid PL ; reads: INV inverted in 40 % of molecules (those are unmethylated),
       IVT flipped in 50 %, PL circular ; Y
  S2   MAG = X with a 2.5 kb element MGE (2 genes) between d5 and d6 ; reads: MGE absent in 40 % (unmethylated), a novel
       3 kb insertion INS between b4 and b5 in 30 % (methylated) ; Y (low coverage)
  S3   MAG = X with C inverted, assembled in two contigs (break between d4 and d5) ; reads span the break
  S4 + S5 co-assembled (G1): MAG = A B D C (C moved after D) + PL, its contig has b6-b7 inverted although no read is
       (misassembly: ALT ~ 100 %)
Every sample's own strain carries 0.3 % SNPs (synonymous / intergenic). CRAM = minimap2 map-hifi to the sample's own
assembly, 3 % QC-fail (0x200), 2 % duplicate (0x400), no reference embedded.
"""
import array, gzip, json, os, random, shutil, subprocess, sys
import pysam

H = os.path.abspath(sys.argv[1]); shutil.rmtree(H, ignore_errors=True); os.makedirs(H)
R = random.Random(11)
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]
def rnd(n): return "".join(R.choice("ACGT") for _ in range(n))
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]
STOP = {"TAA", "TAG", "TGA"}
GC = {a + b + c: aa for (a, b, c), aa in zip(((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
      "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
def translate(s): return "".join(GC.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
def gene(): return "ATG" + "".join(R.choice(CODONS) for _ in range(R.randint(180, 330))) + "TAA"

FAM = {}
def fam(name): FAM[name] = dict(seq=gene(), strand=R.choice((1, -1))); return name
SEG = {k: [fam(f"{k}{i}") for i in range(10)] for k in "abcd"}
PL = [fam(f"p{i}") for i in range(5)]
MGE = [fam("m0"), fam("m1")]
Y = [fam(f"y{i}") for i in range(14)]
SP = {}
def spacer(a, b):              # the same gene pair read the other way round has the reverse-complement spacer
    return rc(SP[(b, a)]) if (b, a) in SP else SP.setdefault((a, b), rnd(R.randint(150, 300)))
INV, IR, IVT, INS = rnd(1200), rnd(25), rnd(130), rnd(3000)
IVT_EL = IR + IVT + rc(IR)                      # 180 bp, invertible between its inverted repeats


def assemble(items):
    """items: ('g', name, strand) | ('n', name, seq) -> sequence, genes [(name, start, end, strand)], marks {name: (start, end)}"""
    seq, genes, marks, pos, prev = [], [], {}, 0, None
    for it in items:
        if prev is not None:
            sp = spacer(prev, it[1]); seq.append(sp); pos += len(sp)
        if it[0] == "g":
            q = FAM[it[1]]["seq"] if it[2] > 0 else rc(FAM[it[1]]["seq"])
            genes.append((it[1], pos + 1, pos + len(q), it[2]))
        else:
            q = it[2]; marks[it[1]] = (pos + 1, pos + len(q))
        seq.append(q); pos += len(q); prev = it[1]
    return "".join(seq), genes, marks


def seg(k, rev=False):
    g = [("g", n, FAM[n]["strand"]) for n in SEG[k]]
    if k == "a": g = g[:5] + [("n", "INV", INV)] + g[5:]
    if k == "c": g = g[:4] + [("n", "IVT", IVT_EL)] + g[4:]
    return [(t, n, (-x if rev else x) if t == "g" else (rc(x) if rev else x)) for t, n, x in (g[::-1] if rev else g)]


def strainify(seq, genes, rate):
    """SNPs at intergenic and third codon positions that create no stop codon"""
    s, cod = list(seq), {}
    for name, a, b, st in genes:
        for i in range(a - 1, b):
            off = (i - (a - 1)) if st > 0 else (b - 1 - i)
            cod[i] = (off % 3, name, a, b, st)
    for i in range(len(s)):
        if R.random() >= rate: continue
        c = cod.get(i)
        if c and c[0] != 2: continue
        old = s[i]; s[i] = R.choice("ACGT".replace(old, ""))
        if c:
            _, _, a, b, st = c
            off = (i - (a - 1)) if st > 0 else (b - 1 - i)
            k0 = off - off % 3
            q = "".join(s[a - 1:b]); q = q if st > 0 else rc(q)
            if off < len(q) - 3 and q[k0:k0 + 3] in STOP: s[i] = old
    return "".join(s)


# ---------------------------------------------------------------- MAGs (assemblies) and the molecules of each sample
X_ANC = seg("a") + seg("b") + seg("c") + seg("d")
LAYOUT = {
    "S1": X_ANC,
    "S2": seg("a") + seg("b") + seg("c") + seg("d")[:6] + [("g", "m0", 1), ("g", "m1", -1)] + seg("d")[6:],
    "S3": seg("a") + seg("b") + seg("c", rev=True) + seg("d"),
    "G1": seg("a") + seg("b") + seg("d") + seg("c"),
}
PLAS = [("g", n, FAM[n]["strand"]) for n in PL]
YL = [("g", n, FAM[n]["strand"]) for n in Y]
truth = dict(motif="GATC", frac=0.9, mags={}, events=[], samples=["S1", "S2", "S3", "S4", "S5"])
asm = {}        # assembly -> dict(dir, prefix, bins, contigs [(name, seq, genes, bin)])
RUN = os.path.join(H, "mf4", "mock")
for a, prefix, binsmp in (("S1", "S1", "S1"), ("S2", "S2", "S2"), ("S3", "S3", "S3"), ("G1", "S5M2", "S5")):
    d = f"{RUN}/{a}/assemblies/metag" if a != "G1" else f"{RUN}/AssmblGrp_G1/metag"
    asm[a] = dict(dir=d, prefix=prefix, bins=binsmp, ctgs=[])
strain = {}
for a in ("S1", "S2", "S3", "G1"):
    q, genes, marks = assemble(LAYOUT[a])
    q = strainify(q, genes, 0.003)
    strain[a] = (q, genes, marks)
# what the reads are made from: per sample a list of (molecule sequence, circular, weight, methylated fraction, tag)
def flip(q, a, b): return q[:a - 1] + rc(q[a - 1:b]) + q[b:]
mol = {}
q1, g1, m1 = strain["S1"]
inv_a, inv_b = m1["INV"]; iv_a, iv_b = m1["IVT"]; iv_a += 25; iv_b -= 25      # the part between the inverted repeats
pl_q, pl_g, _ = assemble(PLAS); pl_q = strainify(pl_q, pl_g, 0.003)
y_q, y_g, _ = assemble(YL)
mol["S1"] = [(q1, False, 0.3, 0.9, "ref"), (flip(q1, inv_a, inv_b), False, 0.2, 0.0, "inv"), (flip(q1, iv_a, iv_b), False, 0.25, 0.9, "ivt"),
             (flip(flip(q1, inv_a, inv_b), iv_a, iv_b), False, 0.15, 0.0, "inv+ivt"), (pl_q, True, 0.1, 0.9, "plasmid"),
             (y_q, False, 0.25, 0.9, "Y")]
q2, g2, m2 = strain["S2"]
mge = (min(a for n, a, b, s in g2 if n == "m0") - 120, max(b for n, a, b, s in g2 if n == "m1") + 120)   # element incl. its spacers
mge_del = q2[:mge[0] - 1] + q2[mge[1]:]
ins_at = max(b for n, a, b, s in g2 if n == "b4") + 60
mol["S2"] = [(q2, False, 0.45, 0.9, "ref"), (mge_del, False, 0.3, 0.0, "mge_del"), (q2[:ins_at] + INS + q2[ins_at:], False, 0.25, 0.9, "ins"),
             (y_q, False, 0.08, 0.9, "Y")]
q3, g3, m3 = strain["S3"]
mol["S3"] = [(q3, False, 1.0, 0.9, "ref")]
qg, gg, mg = strain["G1"]
mis = (min(a for n, a, b, s in gg if n == "b6") - 50, max(b for n, a, b, s in gg if n == "b7") + 50)
mol["S4"] = [(qg, False, 0.85, 0.9, "ref"), (pl_q, True, 0.15, 0.9, "plasmid")]
mol["S5"] = mol["S4"]

def contigs_of(a):
    """the MAG contigs of an assembly: [(name, seq, genes, bin)] ; S3 broken between d4 and d5 ; G1 misassembled"""
    q, genes, marks = strain[a]
    P = asm[a]["prefix"]
    if a == "S3":
        cut = max(b for n, a_, b, s in genes if n == "d4") + 80
        return [(f"{P}__C1", q[:cut], [g for g in genes if g[2] <= cut], "SemiBin_1"),
                (f"{P}__C2", q[cut:], [(n, x - cut, y - cut, s) for n, x, y, s in genes if x > cut], "SemiBin_1")]
    if a == "G1":
        a_, b_ = mis
        qm = flip(q, a_, b_)
        gm = [(n, (a_ + b_ - y) if a_ <= x and y <= b_ else x, (a_ + b_ - x) if a_ <= x and y <= b_ else y, -s if a_ <= x and y <= b_ else s)
              for n, x, y, s in genes]
        return [(f"{P}__C1", qm, sorted(gm, key=lambda g: g[1]), "SemiBin_1"), (f"{P}__C2", pl_q, pl_g, "SemiBin_1")]
    out = [(f"{P}__C1", q, genes, "SemiBin_1")]
    if a == "S1": out.append((f"{P}__C2", pl_q, pl_g, "SemiBin_1"))
    return out

for a in asm:
    asm[a]["ctgs"] = contigs_of(a)
for a in ("S1", "S2"):
    asm[a]["ctgs"].append((f"{asm[a]['prefix']}__C9", y_q, y_g, "SemiBin_2"))
for a in asm:                                     # one unbinned junk contig per assembly
    asm[a]["ctgs"].append((f"{asm[a]['prefix']}__C99", rnd(5000), [], None))

# truth events in MAG contig coordinates (1-based; breakpoint = gap after the given base)
P1, P2, P3, PG = (asm[a]["prefix"] for a in ("S1", "S2", "S3", "G1"))
truth["events"] = [
    dict(sample="S1", MGS="MGS.0", kind="inversion", contig=f"{P1}__C1", bp=[inv_a - 1, inv_b], alt=0.35, alt_meth=0.0),
    dict(sample="S1", MGS="MGS.0", kind="inversion", contig=f"{P1}__C1", bp=[iv_a - 1, iv_b], alt=0.4, alt_meth=0.9, short=True),
    dict(sample="S1", MGS="MGS.0", kind="circular", contig=f"{P1}__C2", bp=[len(pl_q)]),
    dict(sample="S2", MGS="MGS.0", kind="deletion", contig=f"{P2}__C1", bp=[mge[0] - 1, mge[1]], alt=0.4, alt_meth=0.0),
    dict(sample="S2", MGS="MGS.0", kind="insertion", contig=f"{P2}__C1", bp=[ins_at], alt=0.25, insert=len(INS)),
    dict(sample="S3", MGS="MGS.0", kind="contig_join", contig=f"{P3}__C1", bp=[len(asm["S3"]["ctgs"][0][1])], partner=f"{P3}__C2"),
    dict(sample="S4", MGS="MGS.0", kind="misassembly", contig=f"{PG}__C1", bp=[mis[0] - 1, mis[1]]),
    dict(sample="S5", MGS="MGS.0", kind="misassembly", contig=f"{PG}__C1", bp=[mis[0] - 1, mis[1]]),
    dict(sample="S4", MGS="MGS.0", kind="circular", contig=f"{PG}__C2", bp=[len(pl_q)]),
]

# ---------------------------------------------------------------- catalogue (cluster = gene family), files
clusters, members = {}, {}
for a in sorted(asm):
    for name, q, genes, b in asm[a]["ctgs"]:
        for n, (g, x, y, st) in enumerate(genes, 1):
            c = clusters.setdefault(g, len(clusters) + 1)
            gq = q[x - 1:y] if st > 0 else rc(q[x - 1:y])
            members.setdefault(c, []).append((f"{name}_{n}", gq))
GCD = os.path.join(H, "gc"); os.makedirs(f"{GCD}/Bin_SB/Genomes/MGS_ctg", exist_ok=True)
with open(f"{GCD}/compl.incompl.95.prot.faa", "w") as p, open(f"{GCD}/compl.incompl.95.fna.clstr.idx", "w") as ix:
    ix.write("#Gene\tmembers\n")
    for c in sorted(members, key=lambda c: R.random()):
        p.write(f">{c}\n{translate(members[c][0][1])}\n")
        ix.write(f"{c}\t" + ",".join(">" + m for m, _ in members[c]) + "\n")
for a, A in asm.items():
    d = A["dir"]; os.makedirs(f"{d}/genePred", exist_ok=True); os.makedirs(f"{d}/Binning/SB", exist_ok=True)
    with open(f"{d}/scaffolds.fasta.filt", "w") as o:
        for name, q, _, _ in A["ctgs"]: o.write(f">{name}\n" + "\n".join(q[i:i + 100] for i in range(0, len(q), 100)) + "\n")
    with gzip.open(f"{d}/genePred/genes.gff.gz", "wt") as gf, gzip.open(f"{d}/genePred/genes.shrtHD.fna.gz", "wt") as fn:
        gf.write("##gff-version  3\n")
        for k, (name, q, genes, _) in enumerate(A["ctgs"], 1):
            for n, (g, x, y, st) in enumerate(genes, 1):
                gf.write(f"{name}\tProdigal_v2.6.3\tCDS\t{x}\t{y}\t99.9\t{'+' if st > 0 else '-'}\t0\tID={k}_{n};partial=00;start_type=ATG;\n")
                fn.write(f">{name}_{n}\n{q[x - 1:y] if st > 0 else rc(q[x - 1:y])}\n")
    with open(f"{d}/Binning/SB/{A['bins']}", "w") as o:
        o.write("Sequence ID\tBin\n" + "".join(f"{name}\t{b}\n" for name, _, _, b in A["ctgs"] if b))
sample_asm = {"S1": "S1", "S2": "S2", "S3": "S3", "S4": "G1", "S5": "G1"}
for s, a in sample_asm.items():
    sd = f"{RUN}/{s}"; os.makedirs(f"{sd}/assemblies/metag", exist_ok=True); os.makedirs(f"{sd}/mapping", exist_ok=True)
    open(f"{sd}/assemblies/metag/assembly.txt", "w").write(asm[a]["dir"] + "/\n")
mags = [(f"{asm[a]['prefix']}__SemiBin_1", "MGS.0", a) for a in ("S1", "S2", "S3", "G1")] + \
       [(f"{asm[a]['prefix']}__SemiBin_2", "MGS.1", a) for a in ("S1", "S2")]
with gzip.open(f"{GCD}/Bin_SB/MAGvsGC.txt.gz", "wt") as o:
    o.write("MAG\tMGS\tRepresentative4MGS\tCompleteness\tContamination\tN50\n")
    for m, g, a in mags: o.write(f"{m}\t{g}\t{'*' if m in ('S1__SemiBin_1', 'S1__SemiBin_2') else ''}\t98.0\t0.5\t50000\n")
for m, g, a in mags:
    if m.startswith("S1__"):
        b = m.split("__")[1]
        with gzip.open(f"{GCD}/Bin_SB/Genomes/MGS_ctg/{g}.ctgs.{m}.fna.gz", "wt") as o:
            for name, q, _, bb in asm[a]["ctgs"]:
                if bb == b: o.write(f">{name}\n{q}\n")
truth["mags"] = {m: dict(MGS=g, assembly=a) for m, g, a in mags}
truth["contigs"] = {name: q for A in asm.values() for name, q, _, _ in A["ctgs"]}


# ---------------------------------------------------------------- reads with kinetics, CRAMs
def kin(read, fr):
    """ip per base of `read`'s template strand: high at the A of GATC on the template (fraction fr of the sites)"""
    L = len(read); ip = [R.randint(8, 14) for _ in range(L)]; t = rc(read)
    i = t.find("GATC")
    while i >= 0:
        if R.random() < fr: ip[L - 2 - i] = R.randint(45, 60)
        i = t.find("GATC", i + 1)
    return array.array("B", ip)

COV = {"S1": 25, "S2": 25, "S3": 18, "S4": 12, "S5": 12}          # chromosome ; plasmid 30x, Y 15x (S1) / 6x (S2)
TAGCOV = {("S1", "plasmid"): 30, ("S4", "plasmid"): 30, ("S5", "plasmid"): 30, ("S1", "Y"): 15, ("S2", "Y"): 6}
raw_dir = os.path.join(H, "raw"); os.makedirs(raw_dir)
maplines = [f"#OutPath\t{H}/mf4/", "#RunID\tmock", f"#DirPath\t{raw_dir}/", "#SmplID\tSmplPrefix\tAssmblGrps"]
truth["reads"], truth["flag"], truth["tag"] = {}, {}, {}
for s in truth["samples"]:
    hdr = {"HD": {"VN": "1.6", "SO": "unknown", "pb": "5.0.0"},
           "RG": [{"ID": "1a2b3c4d", "PL": "PACBIO", "PU": f"m84_{s}", "PM": "REVIO",
                   "DS": "READTYPE=CCS;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-894-200;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000"}]}
    raw = os.path.join(raw_dir, f"m84_{s}.hifi_reads.bam"); z = 0
    with pysam.AlignmentFile(raw, "wb", header=hdr) as o:
        for src, circ, w, fr, tag in mol[s]:
            n = int(len(src) * TAGCOV[(s, tag)] / 6000) if (s, tag) in TAGCOV else int(len(src) * COV[s] * w / 6000)
            for _ in range(n):
                L = R.randint(4000, 8000)
                if circ:
                    st = R.randrange(len(src)); seq = (src + src)[st:st + min(L, len(src) - 1)]
                else:
                    L = min(L, len(src)); st = R.randint(0, len(src) - L); seq = src[st:st + L]
                seq = "".join(R.choice("ACGT".replace(b, "")) if R.random() < 0.001 else b for b in seq)
                if R.random() < .5: seq = rc(seq)
                z += 1; q = f"m84_{s}/{z * 7 + 1}/ccs"
                truth["tag"][q] = tag
                a = pysam.AlignedSegment(o.header); a.query_name = q; a.flag = 4; a.query_sequence = seq
                a.set_tag("RG", "1a2b3c4d"); a.set_tag("zm", z * 7 + 1); a.set_tag("rq", 0.999); a.set_tag("np", 12)
                a.set_tag("fi", kin(seq, fr)); a.set_tag("ri", kin(rc(seq), fr))
                o.write(a)
    maplines.append(f"{s}\tm84_{s}.hifi_reads\t{'G1' if s in ('S4', 'S5') else ''}")
    ref = os.path.join(asm[sample_asm[s]]["dir"], "scaffolds.fasta.filt")
    cram = os.path.join(RUN, s, "mapping", f"{s}-smd.cram")
    sam = subprocess.run(f"samtools fastq {raw} 2>/dev/null | minimap2 -ax map-hifi -t 2 {ref} - 2>/dev/null",
                         shell=True, check=True, capture_output=True, text=True).stdout.splitlines()
    out = []
    for line in sam:
        if not line.startswith("@"):
            f = line.split("\t"); fl = int(f[1])
            if not fl & 0x904:
                r = R.random()
                if r < 0.03: fl |= 0x200
                elif r < 0.05: fl |= 0x400
                truth["flag"][f[0]] = fl
            f[1] = str(fl); line = "\t".join(f)
        out.append(line)
    subprocess.run(f"samtools sort -O cram --reference {ref} -o {cram} -", shell=True, check=True, input="\n".join(out) + "\n", text=True)
open(os.path.join(H, "map.txt"), "w").write("\n".join(maplines) + "\n")
json.dump(truth, open(os.path.join(H, "truth.json"), "w"))
print("mock:", H, {s: sum(1 for k in truth["tag"] if k.startswith(f"m84_{s}/")) for s in truth["samples"]})
