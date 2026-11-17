#!/usr/bin/env python3
"""A small MATAFILER4 run + GeneCat MGS stage for testing `mmt export-mot`, with exact truth (truth.json).

usage: motexport_mock.py OUTDIR        (OUTDIR is wiped; needs pysam + samtools)

Samples (raw PacBio BAMs with kinetics unless said; CRAMs written as MATAFILER4 does: sample reads -> own assembly)
  S1   test A + F: one BAM, MAG S1__SemiBin_1 (1 contig, the MGS.0 representative); CRAM primaries flagged 0x200 / 0x400
       (kept, selected); a low-quality bin (not in MAGvsGC) and an unbinned contig -> off-MAG ; input_raw.txt names a
       moved BAM (missing) -> resolved through the map
  S2   tests B + C: MAG S2__SemiBin_1 = 3 contigs ; MAG S2__SemiBin_2 = one circular contig (scaffolds.fasta.circ),
       MGS "?" (no MGS) ; GeneCat FASTA of S2__SemiBin_1 (MGS.1 representative)
  S3   tests D + G + H: three BAMs (three movies, input_raw.txt lists them) ; MAGs S3__SemiBin_1/2/3: a molecule split
       over MAG 2 and 3 (supplementary), one at MAPQ 0, one MAG + unbinned contig, primary + supplementary on one MAG
       (listed once), by-strand names (/ccs/fwd + /ccs/rev) of one molecule
  S4 S5 test E: co-assembly group G (aliases S4M1, S5M2), one MAG S5M2__SemiBin_1 ; each sample's own CRAM and BAM
  S6   test I: PacBio CRAM, no raw BAM anywhere (map folder empty, input_raw.txt names a missing file)
  S7   test J: Illumina-only (fastq inputs, Illumina read names in the CRAM)
  S8   test J: ONT BAM (no kinetics, UUID read names)
  S9   PacBio BAM without kinetics
  S10  test K: CRAM made against an older assembly: contig S10__C1 was polished afterwards (+12 bp)
"""
import array, gzip, json, os, random, shutil, subprocess, sys
import pysam

H = os.path.abspath(sys.argv[1]); shutil.rmtree(H, ignore_errors=True)
R = random.Random(11)
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]
def rnd(n): return "".join(R.choice("ACGT") for _ in range(n))

RUN = f"{H}/mf4/run"; RAW = f"{H}/raw"; GC = f"{H}/gc"
os.makedirs(RAW); os.makedirs(f"{GC}/Bin_SB/Genomes/MGS_ctg"); os.makedirs(f"{GC}/LOGandSUB")
truth = {"selected": {}, "status": {}, "rows": [], "mags": {}, "unresolved_samples": ["S6", "S9", "S10"],
         "no_rows": ["S6", "S7", "S8", "S9"]}
CHEM = "BINDINGKIT=103-496-900;SEQUENCINGKIT=103-496-700"


def movie(s, k=1):
    return f"m84{sum(map(ord, s)) % 900 + 100:03d}_251101_{R.randint(100000, 999999)}_s{k}"


# ---------------------------------------------------------------- assemblies: name -> (dir, [(contig, seq, bin)])
asm = {}
def assembly(key, d, ctgs, binfile, circ=()):
    os.makedirs(f"{d}/Binning/SB", exist_ok=True)
    with open(f"{d}/scaffolds.fasta.filt", "w") as o:
        for c, s, _ in ctgs:
            o.write(f">{c} len={len(s)}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    with open(f"{d}/Binning/SB/{binfile}", "w") as o:
        o.write("Sequence ID\tBin\n" + "".join(f"{c}\t{b}\n" for c, _, b in ctgs if b))
    if circ:
        with open(f"{d}/scaffolds.fasta.circ", "w") as o:
            for c, s, _ in ctgs:
                if c in circ:
                    o.write(f">{c}\n{s}\n")
    asm[key] = (d, ctgs)

assembly("S1", f"{RUN}/S1/assemblies/metag", [("S1__C1", rnd(12000), "SemiBin_1"), ("S1__C2", rnd(6000), "SemiBin_9"),
                                              ("S1__C3", rnd(5000), "")], "S1")
assembly("S2", f"{RUN}/S2/assemblies/metag", [("S2__C1", rnd(7000), "SemiBin_1"), ("S2__C2", rnd(5000), "SemiBin_1"),
                                              ("S2__C3", rnd(4000), "SemiBin_1"), ("S2__C4", rnd(9000), "SemiBin_2")],
         "S2", circ=("S2__C4",))
SH = rnd(3000)                                                       # a region MAG 2 and MAG 3 share
assembly("S3", f"{RUN}/S3/assemblies/metag", [("S3__C1", rnd(8000), "SemiBin_1"), ("S3__C2", rnd(6000) + SH, "SemiBin_2"),
                                              ("S3__C3", SH + rnd(6000), "SemiBin_3"), ("S3__C4", rnd(4000), "")], "S3")
assembly("G", f"{RUN}/AssmblGrp_G/metag", [("S5M2__C1", rnd(9000), "SemiBin_1")], "S5")
assembly("S6", f"{RUN}/S6/assemblies/metag", [("S6__C1", rnd(6000), "SemiBin_1")], "S6")
assembly("S7", f"{RUN}/S7/assemblies/metag", [("S7__C1", rnd(6000), "SemiBin_1")], "S7")
assembly("S8", f"{RUN}/S8/assemblies/metag", [("S8__C1", rnd(6000), "SemiBin_1")], "S8")
assembly("S9", f"{RUN}/S9/assemblies/metag", [("S9__C1", rnd(6000), "SemiBin_1")], "S9")
OLD10 = [("S10__C1", rnd(7000), "SemiBin_1"), ("S10__C2", rnd(5000), "SemiBin_1")]
NEW10 = [(OLD10[0][0], OLD10[0][1][:3500] + "ACGTACGTACGT" + OLD10[0][1][3500:], "SemiBin_1"), OLD10[1]]
assembly("S10", f"{RUN}/S10/assemblies/metag", NEW10, "S10")
sample_asm = {"S1": "S1", "S2": "S2", "S3": "S3", "S4": "G", "S5": "G", "S6": "S6", "S7": "S7", "S8": "S8", "S9": "S9",
              "S10": "S10"}
for s, a in sample_asm.items():
    os.makedirs(f"{RUN}/{s}/mapping", exist_ok=True); os.makedirs(f"{RUN}/{s}/assemblies/metag", exist_ok=True)
    open(f"{RUN}/{s}/assemblies/metag/assembly.txt", "w").write(asm[a][0] + "/\n")
seq = {c: s for a in asm.values() for c, s, _ in a[1]}

# ---------------------------------------------------------------- GeneCat
mvg = [("S1__SemiBin_1", "MGS.0", "*"), ("S2__SemiBin_1", "MGS.1", "*"), ("S2__SemiBin_2", "?", ""),
       ("Cano__5", "MGS.0", ""), ("S3__SemiBin_1", "MGS.0", ""), ("S3__SemiBin_2", "MGS.2", "*"),
       ("S3__SemiBin_3", "MGS.3", "*"), ("S5M2__SemiBin_1", "MGS.1", ""), ("S6__SemiBin_1", "MGS.4", "*"),
       ("S7__SemiBin_1", "MGS.4", ""), ("S8__SemiBin_1", "MGS.4", ""), ("S9__SemiBin_1", "MGS.4", ""),
       ("S10__SemiBin_1", "MGS.5", "*")]
with gzip.open(f"{GC}/Bin_SB/MAGvsGC.txt.gz", "wt") as o:
    o.write("MAG\tMGS\tRepresentative4MGS\tCompleteness\tContamination\n")
    for m, g, st in mvg:
        o.write(f"{m}\t{g}\t{st}\t91.0\t1.5\n")
for m, g, a, b in (("S1__SemiBin_1", "MGS.0", "S1", "SemiBin_1"), ("S2__SemiBin_1", "MGS.1", "S2", "SemiBin_1")):
    with gzip.open(f"{GC}/Bin_SB/Genomes/MGS_ctg/{g}.ctgs.{m}.fna.gz", "wt") as o:
        for c, s, bb in asm[a][1]:
            if bb == b:
                o.write(f">{c}\n{s}\n")
for m, _, _ in mvg:
    if not m.startswith("Cano"):
        a = sample_asm[m.split("__")[0].replace("S5M2", "S5")]
        truth["mags"][m] = [c for c, _, b in asm[a][1] if b == m.split("__")[1]]

# ---------------------------------------------------------------- raw BAMs
def write_bam(path, mv, names, kinetics=True, readtype="CCS"):
    hdr = {"HD": {"VN": "1.6", "SO": "unknown", "pb": "5.0.0"},
           "RG": [{"ID": "rg" + mv[-6:], "PL": "PACBIO", "PU": mv, "PM": "REVIO",
                   "DS": f"READTYPE={readtype};Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;{CHEM};BASECALLERVERSION=5.0.0;"
                         "FRAMERATEHZ=100.000000"}]}
    with pysam.AlignmentFile(path, "wb", header=hdr) as o:
        for q in names:
            a = pysam.AlignedSegment(o.header); a.query_name = q; a.flag = 4; a.query_sequence = rnd(300)
            a.query_qualities = pysam.qualitystring_to_array("I" * 300)
            a.set_tag("RG", hdr["RG"][0]["ID"])
            if kinetics:
                a.set_tag("fi", array.array("B", [R.randint(5, 40) for _ in range(300)]))
                a.set_tag("ri", array.array("B", [R.randint(5, 40) for _ in range(300)]))
                a.set_tag("fn", 12); a.set_tag("rn", 12)
            o.write(a)

# ---------------------------------------------------------------- CRAM records
class Cram:
    def __init__(self, s, ref_ctgs):
        self.s, self.ref, self.recs, self.n = s, ref_ctgs, [], 0
    def aln(self, q, ctg, pos, L, flag=0, mapq=60, clip=None):
        s = seq[ctg] if ctg in seq else dict(self.ref)[ctg]
        s = dict(self.ref).get(ctg, s)
        sq = s[pos - 1:pos - 1 + L]
        cig = f"{L}M"
        if clip:                                                     # (left, right) hard clips of a split read
            cig = (f"{clip[0]}H" if clip[0] else "") + cig + (f"{clip[1]}H" if clip[1] else "")
        self.recs.append((q, flag, ctg, pos, mapq, cig, sq))
    def write(self):
        d = f"{RUN}/{self.s}/mapping"
        ref = f"{d}/.ref.fa"
        with open(ref, "w") as o:
            for c, s in self.ref:
                o.write(f">{c}\n{s}\n")
        pysam.faidx(ref)
        sam = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(f"@SQ\tSN:{c}\tLN:{len(s)}\n" for c, s in self.ref)
        for q, flag, c, pos, mq, cig, sq in self.recs:
            sam += f"{q}\t{flag}\t{c}\t{pos}\t{mq}\t{cig}\t*\t0\t0\t{sq}\t*\n"
        subprocess.run(["samtools", "sort", "-O", "cram", "--reference", ref, "-o", f"{d}/{self.s}-smd.cram", "-"],
                       input=sam, text=True, check=True)
        os.remove(ref); os.remove(ref + ".fai")
        open(f"{d}/{self.s}-smd.cram.retained.size", "w").write("1000\n")

def ctgs(a): return [(c, s) for c, s, _ in asm[a][1]]
def sel(s, mag, bam, names):
    truth["selected"].setdefault(f"{s}|{mag}|{bam}", []).extend(names)
def stat(s, name, st):
    truth["status"][f"{s}|{name}"] = st

maplines = ["#SmplID\tSmplPrefix\tReads\tSeqTech\tAssmblGrps\tNote", f"#OutPath\t{H}/mf4/", "#RunID\trun"]

# S1: one BAM; QC-flagged primaries; off-MAG reads; input_raw.txt names a moved BAM
mv = movie("S1"); os.makedirs(f"{RAW}/c1")
cr = Cram("S1", ctgs("S1")); names = []
for z in range(1, 41):
    q = f"{mv}/{z}/ccs"; names.append(q)
    if z <= 30:
        fl = 0x200 if z in (3, 4) else 0x400 if z == 5 else 0
        cr.aln(q, "S1__C1", 100 + z * 300, 1500, flag=fl)
        stat("S1", f"{mv}/{z}", "unique"); sel("S1", "S1__SemiBin_1", "S1.hifi_reads", [f"{mv}/{z}"])
    elif z <= 35:
        cr.aln(q, "S1__C2", 100 + z * 50, 1500)                     # low-quality bin: not a final MAG
    else:
        cr.aln(q, "S1__C3", 100 + z * 50, 1500)                     # unbinned
write_bam(f"{RAW}/c1/S1.hifi_reads.bam", mv, names); cr.write()
open(f"{RUN}/S1/input_raw.txt", "w").write(f"{RAW}/old_location/S1.hifi_reads.bam")
maplines += [f"#DirPath\t{RAW}/c1/", "S1\tS1.hifi_reads\t\tPB\t\t"]

# S2: multi-contig MAG + circular MAG
mv = movie("S2"); cr = Cram("S2", ctgs("S2")); names = []
for z, (c, mag) in enumerate([("S2__C1", 1), ("S2__C2", 1), ("S2__C3", 1), ("S2__C4", 2)] * 6, 1):
    q = f"{mv}/{z}/ccs"; names.append(q); cr.aln(q, c, 200 + z * 40, 1200)
    sel("S2", f"S2__SemiBin_{mag}", "S2.hifi_reads", [f"{mv}/{z}"])
os.makedirs(f"{RAW}/c2"); write_bam(f"{RAW}/c2/S2.hifi_reads.bam", mv, names); cr.write()
maplines += [f"#DirPath\t{RAW}/c2/", "S2\tS2.hifi_reads\t\tPB\t\t"]

# S3: three BAMs ; competition cases
mvs = [movie("S3", k) for k in (1, 2, 3)]; cr = Cram("S3", ctgs("S3")); per = {m: [] for m in mvs}
os.makedirs(f"{RAW}/c3")
z = 0
for i, m in enumerate(mvs):
    lab = f"S3.run{i + 1}.hifi_reads"
    for _ in range(8):                                               # unique on MAG 1 (D)
        z += 1; q = f"{m}/{z}/ccs"; per[m].append(q); cr.aln(q, "S3__C1", 100 + z * 50, 1500)
        sel("S3", "S3__SemiBin_1", lab, [f"{m}/{z}"]); stat("S3", f"{m}/{z}", "unique")
m = mvs[0]; lab = "S3.run1.hifi_reads"
z += 1; q = f"{m}/{z}/ccs"; per[m].append(q)                         # G: split over MAG 2 and MAG 3
cr.aln(q, "S3__C2", 1000, 1500, clip=(0, 800)); cr.aln(q, "S3__C3", 4000, 800, flag=0x800, clip=(1500, 0))
stat("S3", f"{m}/{z}", "split_across_mags")
z += 1; q = f"{m}/{z}/ccs"; per[m].append(q)                         # G: inside the shared region, MAPQ 0
cr.aln(q, "S3__C2", 6100, 2500, mapq=0); stat("S3", f"{m}/{z}", "ambiguous_low_mapq")
z += 1; q = f"{m}/{z}/ccs"; per[m].append(q)                         # MAG 2 + unbinned contig
cr.aln(q, "S3__C2", 2000, 1500, clip=(0, 600)); cr.aln(q, "S3__C4", 500, 600, flag=0x800, clip=(1500, 0))
stat("S3", f"{m}/{z}", "split_mag_unbinned")
z += 1; q = f"{m}/{z}/ccs"; per[m].append(q)                         # H: primary + supplementary on MAG 3
cr.aln(q, "S3__C3", 4000, 1500, clip=(0, 700)); cr.aln(q, "S3__C3", 7000, 700, flag=0x800, clip=(1500, 0))
stat("S3", f"{m}/{z}", "unique"); sel("S3", "S3__SemiBin_3", lab, [f"{m}/{z}"])
z += 1                                                               # H: by-strand records of one molecule
for st_ in ("fwd", "rev"):
    q = f"{m}/{z}/ccs/{st_}"; per[m].append(q); cr.aln(q, "S3__C2", 1500, 1200, flag=0 if st_ == "fwd" else 16)
stat("S3", f"{m}/{z}", "unique"); sel("S3", "S3__SemiBin_2", lab, [f"{m}/{z}"])
for i, m in enumerate(mvs):
    write_bam(f"{RAW}/c3/S3.run{i + 1}.hifi_reads.bam", m, per[m])
cr.write()
open(f"{RUN}/S3/input_raw.txt", "w").write(";".join(f"{RAW}/c3/S3.run{i}.hifi_reads.bam" for i in (1, 2, 3)))
maplines += [f"#DirPath\t{RAW}/c3/", "S3\tS3.run\t\tPB\t\t"]

# S4, S5: co-assembly G ; MAG S5M2__SemiBin_1 ; each sample its own BAM + CRAM
os.makedirs(f"{RAW}/c4")
for s in ("S4", "S5"):
    mv = movie(s); cr = Cram(s, ctgs("G")); names = []
    for z in range(1, 13):
        q = f"{mv}/{z}/ccs"; names.append(q); cr.aln(q, "S5M2__C1", 100 + z * 400, 1500)
        sel(s, "S5M2__SemiBin_1", f"{s}.hifi_reads", [f"{mv}/{z}"])
    write_bam(f"{RAW}/c4/{s}.hifi_reads.bam", mv, names); cr.write()
maplines += [f"#DirPath\t{RAW}/c4/", "S4\tS4.hifi_reads\t\tPB\tG\t", "S5\tS5.hifi_reads\t\tPB\tG\t"]

# S6: PacBio CRAM, no raw BAM
mv = movie("S6"); cr = Cram("S6", ctgs("S6"))
for z in range(1, 6):
    cr.aln(f"{mv}/{z}/ccs", "S6__C1", 100 + z * 300, 1500)
cr.write(); os.makedirs(f"{RAW}/c6")
open(f"{RUN}/S6/input_raw.txt", "w").write(f"{RAW}/gone/S6.hifi_reads.bam")
maplines += [f"#DirPath\t{RAW}/c6/", "S6\tS6.hifi_reads\t\tPB\t\t"]

# S7: Illumina only
cr = Cram("S7", ctgs("S7")); os.makedirs(f"{RAW}/c7")
for z in range(1, 6):
    cr.aln(f"A00123:8:H5KTJ:1:1101:{1000 + z}:2000", "S7__C1", 100 + z * 300, 150)
cr.write()
for r in (1, 2):
    with gzip.open(f"{RAW}/c7/S7_R{r}.fastq.gz", "wt") as o:
        o.write("@x\nACGT\n+\nIIII\n")
open(f"{RUN}/S7/input_raw.txt", "w").write(f"{RAW}/c7/S7_R1.fastq.gz;{RAW}/c7/S7_R2.fastq.gz")
maplines += [f"#DirPath\t{RAW}/c7/", "S7\tS7_\t\tIllumina\t\t"]

# S8: ONT BAM (no kinetics, UUID names) ; S9: PacBio BAM without kinetics
os.makedirs(f"{RAW}/c8")
cr = Cram("S8", ctgs("S8")); names = []
for z in range(1, 6):
    q = f"0b4c{z:04d}-1f2e-4a5b-9c8d-{z:012d}"; names.append(q); cr.aln(q, "S8__C1", 100 + z * 300, 1500)
cr.write(); write_bam(f"{RAW}/c8/S8.ont.bam", "ONTRUN", names, kinetics=False, readtype="")
mv = movie("S9"); cr = Cram("S9", ctgs("S9")); names = []
for z in range(1, 6):
    q = f"{mv}/{z}/ccs"; names.append(q); cr.aln(q, "S9__C1", 100 + z * 300, 1500)
cr.write(); write_bam(f"{RAW}/c8/S9.hifi_reads.bam", mv, names, kinetics=False)
maplines += [f"#DirPath\t{RAW}/c8/", "S8\tS8.ont\t\tONT\t\t", "S9\tS9.hifi_reads\t\tPB\t\t"]

# S10: CRAM against the old assembly (S10__C1 polished afterwards, +12 bp)
mv = movie("S10"); cr = Cram("S10", [(c, s) for c, s, _ in OLD10]); names = []
for z in range(1, 11):
    q = f"{mv}/{z}/ccs"; names.append(q)
    c = "S10__C1" if z <= 5 else "S10__C2"
    cr.aln(q, c, 100 + z * 300, 1500)
    if c == "S10__C2":
        sel("S10", "S10__SemiBin_1", "S10.hifi_reads", [f"{mv}/{z}"])
    stat("S10", f"{mv}/{z}", "reference_mismatch" if c == "S10__C1" else "unique")
os.makedirs(f"{RAW}/c10"); write_bam(f"{RAW}/c10/S10.hifi_reads.bam", mv, names); cr.write()
maplines += [f"#DirPath\t{RAW}/c10/", "S10\tS10.hifi_reads\t\tPB\t\t"]

# ---------------------------------------------------------------- maps: the user's + GeneCat's records of it
open(f"{H}/map.txt", "w").write("\n".join(maplines) + "\n")
shutil.copy(f"{H}/map.txt", f"{GC}/LOGandSUB/map.0.txt")
open(f"{GC}/LOGandSUB/inmap.txt", "w").write("map.0.txt\n")
open(f"{GC}/LOGandSUB/GCmaps.ori", "w").write(f"{H}/map.txt")
clean = [l for l in maplines if not any(l.startswith(x + "\t") for x in ("S6", "S7", "S8", "S9", "S10"))]
open(f"{H}/map_clean.txt", "w").write("\n".join(clean) + "\n")
json.dump(truth, open(f"{H}/truth.json", "w"), indent=1)
print("mock:", H)
