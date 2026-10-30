#!/usr/bin/env python3
"""A small synthetic MATAFILER4 run + GeneCat MGS stage with PacBio kinetics, for testing MMT end to end.

usage: mock_mf4.py OUTDIR        (OUTDIR is wiped; needs pysam, samtools, minimap2)

Species / MAGs (truth in OUTDIR/truth.json)
  A   2 contigs; m6A on the template A of GAAG, methylated fraction per sample FR[sample]
  A2  strain of A in S3: 0.3 % SNPs + a 4 kb island inside contig 1
  A3  second strain of A in S1 (1 % SNPs on contig 2), its own bin -> --mgs-rep uses it, --mag-rep does not
  B   2 contigs; m6A on GATC (both strands), fraction 0.9 everywhere
  C   binned, in no MGS (its reads must never be used) ; J = unbinned junk contig
Samples: S1 (A + A3 + B + J) ; S2 + S2.run2 co-assembled (group G2: aliases S2M1, S2.run2M2) (A + B) ; S3 (A2 + C)
MGS.0  rep S1__SemiBin_1 (A) ; S1__SemiBin_6 (A3) ; S2.run2M2__SemiBin_3 (A) ; S3__SemiBin_2 (A2) ; S9__SemiBin_1 (no sample)
MGS.1  rep S2.run2M2__SemiBin_5.fa.gz (B; MAG name ends .fa.gz like real SemiBin names) ; S1__SemiBin_4 (B) ; canopy row
MGS.3  rep S3__SemiBin_7 - not in the bin file (unusable MGS)
GeneCat abundance (Bin_SB/Annotation/Abundance/MGS.matL0.txt): MGS.1 also in S3, where no MAG of it was binned
CRAM (MATAFILER4 style): minimap2 map-hifi to the sample's own assembly, CRAM without embedded reference,
  3 % of primary alignments flagged QC-fail (0x200), 2 % duplicate (0x400), unmapped reads kept.
"""
import array, gzip, json, os, random, shutil, subprocess, sys
import pysam

H = os.path.abspath(sys.argv[1]); shutil.rmtree(H, ignore_errors=True)
random.seed(7)
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]
def rnd(n): return "".join(random.choice("ACGT") for _ in range(n))
def mut(s, r): return "".join(random.choice("ACGT".replace(b, "")) if random.random() < r else b for b in s)

A = [rnd(36000), rnd(24000)]
ISLAND = rnd(4000); ISL_AT = 18000
A2 = [mut(A[0][:ISL_AT], 0.003) + ISLAND + mut(A[0][ISL_AT:], 0.003), mut(A[1], 0.003)]
A3 = [mut(A[1], 0.01)]
B = [rnd(30000), rnd(18000)]
C = [rnd(20000)]
J = [rnd(8000)]
FR = {"S1": 0.9, "S2": 0.5, "S2.run2": 0.5, "S3": 0.9}
MOTIF = {"A": "GAAG", "A2": "GAAG", "A3": "GAAG", "B": "GATC", "C": "GAAG", "J": "GAAG"}
DEPTH = {"A": 20, "A2": 20, "A3": 10, "B": 20, "C": 15, "J": 15}

RUN = os.path.join(H, "mf4", "mock")                 # MATAFILER4 run folder = #OutPath + #RunID
asm = {   # assembly -> (dir, [(contig, seq, species, bin)])
    "S1": (f"{RUN}/S1/assemblies/metag", [("S1__C1", A[0], "A", "SemiBin_1"), ("S1__C2", A[1], "A", "SemiBin_1"),
                                         ("S1__C6", A3[0], "A3", "SemiBin_6"),
                                         ("S1__C3", B[0], "B", "SemiBin_4"), ("S1__C4", B[1], "B", "SemiBin_4"),
                                         ("S1__C5", J[0], "J", "0")]),
    "G2": (f"{RUN}/AssmblGrp_G2/metag", [("S2.run2M2__C1", A[0], "A", "SemiBin_3"), ("S2.run2M2__C2", A[1], "A", "SemiBin_3"),
                                         ("S2.run2M2__C3", B[0], "B", "SemiBin_5.fa.gz"),
                                         ("S2.run2M2__C4", B[1], "B", "SemiBin_5.fa.gz")]),
    "S3": (f"{RUN}/S3/assemblies/metag", [("S3__C1", A2[0], "A2", "SemiBin_2"), ("S3__C2", A2[1], "A2", "SemiBin_2"),
                                         ("S3__C3", C[0], "C", "SemiBin_9")]),
}
bin_file = {"S1": "S1", "G2": "S2.run2", "S3": "S3"}   # GeneCat reads Binning/SB/<SmplID> (alias resolved)
sample_asm = {"S1": "S1", "S2": "G2", "S2.run2": "G2", "S3": "S3"}
for a, (d, ctgs) in asm.items():
    os.makedirs(os.path.join(d, "Binning", "SB"), exist_ok=True)
    with open(os.path.join(d, "scaffolds.fasta.filt"), "w") as o:
        for c, s, _, _ in ctgs: o.write(f">{c} len={len(s)}\n{s}\n")
    with open(os.path.join(d, "Binning", "SB", bin_file[a]), "w") as o:
        o.write("Sequence ID\tBin\n" + "".join(f"{c}\t{b}\n" for c, _, _, b in ctgs))
for s, a in sample_asm.items():
    sd = os.path.join(RUN, s); os.makedirs(os.path.join(sd, "assemblies", "metag"), exist_ok=True)
    os.makedirs(os.path.join(sd, "mapping"), exist_ok=True)
    open(os.path.join(sd, "assemblies", "metag", "assembly.txt"), "w").write(asm[a][0] + "/\n")

# ---- GeneCat: MAGvsGC (block per MGS, rows in GeneCat order) + GeneCat's representative FASTAs
gc = os.path.join(H, "gc", "Bin_SB"); os.makedirs(os.path.join(gc, "Genomes", "MGS_ctg"), exist_ok=True)
mvg = [("S2.run2M2__SemiBin_3", "MGS.0", ""), ("S1__SemiBin_1", "MGS.0", "*"), ("S1__SemiBin_6", "MGS.0", ""),
       ("S3__SemiBin_2", "MGS.0", ""), ("S9__SemiBin_1", "MGS.0", ""),
       ("Cano__77", "MGS.1", ""), ("S2.run2M2__SemiBin_5.fa.gz", "MGS.1", "*"), ("S1__SemiBin_4", "MGS.1", ""),
       ("S3__SemiBin_7", "MGS.3", "*"), ("S3__SemiBin_9", "?", "")]
with gzip.open(os.path.join(gc, "MAGvsGC.txt.gz"), "wt") as o:
    o.write("MAG\tMGS\tRepresentative4MGS\tCompleteness\tContamination\tN50\tother\n")
    for mag, m, st in mvg: o.write(f"{mag}\t{m}\t{st}\t95.1\t1.2\t5000\tg1,g2\n")
# GeneCat's MGS abundance (R-style matrix: header one field short) ; MGS.1 present in S3 without a MAG there
os.makedirs(os.path.join(gc, "Annotation", "Abundance"), exist_ok=True)
with open(os.path.join(gc, "Annotation", "Abundance", "MGS.matL0.txt"), "w") as o:
    o.write("S1\tS2\tS2.run2\tS3\n" "MGS.0\t20.1\t18.2\t18.4\t19.9\n" "MGS.1\t20.3\t20.0\t19.8\t0.4\n" "MGS.3\t0\t0\t0\t14.0\n")
for m, mag, a, b in (("MGS.0", "S1__SemiBin_1", "S1", "SemiBin_1"), ("MGS.1", "S2.run2M2__SemiBin_5.fa.gz", "G2", "SemiBin_5.fa.gz")):
    name = f"{m}.ctgs.{mag}"; name = name[:-3] if name.endswith(".gz") else name
    name = (name if name.endswith((".fa", ".fna", ".fasta")) else name + ".fna") + ".gz"
    with gzip.open(os.path.join(gc, "Genomes", "MGS_ctg", name), "wt") as o:
        for c, s, _, bb in asm[a][1]:
            if bb == b: o.write(f">{c}\n{s}\n")

# ---- raw kinetics BAMs + MATAFILER-style CRAMs
def kin(read, motif, fr):
    L = len(read); ip = [random.randint(8, 14) for _ in range(L)]; t = rc(read)
    i = t.find(motif)
    while i >= 0:
        if random.random() < fr: ip[L - 2 - i] = random.randint(45, 60)   # template A at motif position 2
        i = t.find(motif, i + 1)
    return array.array("B", ip)

comm = {"S1": [("A", A), ("A3", A3), ("B", B), ("J", J)], "S2": [("A", A), ("B", B)],
        "S2.run2": [("A", A), ("B", B)], "S3": [("A2", A2), ("C", C)]}
raw_dir = os.path.join(H, "raw"); os.makedirs(raw_dir, exist_ok=True)
maplines = [f"#OutPath\t{H}/mf4/", "#RunID\tmock", f"#DirPath\t{raw_dir}/", "#SmplID\tSmplPrefix\tAssmblGrps"]
truth = {"FR": FR, "reads": {}, "flag": {}}
for s, sp in comm.items():
    hdr = {"HD": {"VN": "1.6", "SO": "unknown", "pb": "5.0.0"},
           "RG": [{"ID": "1a2b3c4d", "PL": "PACBIO", "PU": f"m84_{s}", "PM": "REVIO",
                   "DS": "READTYPE=CCS;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-894-200;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000"}]}
    raw = os.path.join(raw_dir, f"m84_{s}.hifi_reads.bam"); z = 0
    with pysam.AlignmentFile(raw, "wb", header=hdr) as o:
        for name, ctgs in sp:
            tot = sum(map(len, ctgs)); n = int(tot * DEPTH[name] / 5000)
            fr = 0.9 if name == "B" else FR[s]
            for _ in range(n):
                g = random.choices(ctgs, [len(c) for c in ctgs])[0]; L = min(random.randint(3000, 7000), len(g))
                st = random.randint(0, len(g) - L); seq = mut(g[st:st + L], 0.001)
                if random.random() < .5: seq = rc(seq)
                z += 1; q = f"m84_{s}/{z * 7 + 1}/ccs"
                truth["reads"][q] = name
                a = pysam.AlignedSegment(o.header); a.query_name = q; a.flag = 4; a.query_sequence = seq
                a.set_tag("RG", "1a2b3c4d"); a.set_tag("zm", z * 7 + 1); a.set_tag("rq", 0.999); a.set_tag("np", 12)
                a.set_tag("fi", kin(seq, MOTIF[name], fr)); a.set_tag("ri", kin(rc(seq), MOTIF[name], fr))
                o.write(a)
    maplines.append(f"{s}\tm84_{s}.hifi_reads\t{'G2' if s.startswith('S2') else ''}")
    ref = os.path.join(asm[sample_asm[s]][0], "scaffolds.fasta.filt")
    cram = os.path.join(RUN, s, "mapping", f"{s}-smd.cram")
    sam = subprocess.run(f"samtools fastq {raw} 2>/dev/null | minimap2 -ax map-hifi -t 2 {ref} - 2>/dev/null",
                         shell=True, check=True, capture_output=True, text=True).stdout.splitlines()
    out = []
    for line in sam:
        if not line.startswith("@"):
            f = line.split("\t"); fl = int(f[1])
            if not fl & 0x904:
                r = random.random()
                if r < 0.03: fl |= 0x200
                elif r < 0.05: fl |= 0x400
                truth["flag"][f[0]] = fl
            f[1] = str(fl); line = "\t".join(f)
        out.append(line)
    subprocess.run(f"samtools sort -O cram --reference {ref} -o {cram} -", shell=True, check=True,
                   input="\n".join(out) + "\n", text=True)
open(os.path.join(H, "map.txt"), "w").write("\n".join(maplines) + "\n")
truth["seq"] = {"A": A, "A2": A2, "A3": A3, "B": B}; truth["island"] = [ISL_AT, len(ISLAND)]
json.dump(truth, open(os.path.join(H, "truth.json"), "w"))
print("mock MATAFILER4 run:", H)
