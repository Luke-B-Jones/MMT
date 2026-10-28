#!/usr/bin/env python3
"""MMT ctg-meth engine (method and every column: docs/ctg-meth.md).

  refs       MAG reference FASTAs from one assembly (one pass over it)
  junctions  pbmm2 alignments of a sample's pseudo-reads to their own MAGs (name-sorted) -> pieces (reads table),
             breakpoints, read sets (REF / ALT conformations), windows + the window BAM for ipdSummary
  write      sorted kinetics export + derep-ctg lift -> per MGS: <S>.parquet, junctions/<S>.parquet, edges/<S>.parquet,
             reads/<S>.parquet
  relift     the same files re-annotated from a newer derep-ctg result (kinetics untouched)
  derep_mags the MAGs of derep-ctg results (is a result current for the index?)

Coordinates are 1-based and inclusive on contigs, blocks and reads. A breakpoint end (BE) = (contig, pos, side): side R
= the sequence ending at pos is joined, at its right end, to something else ; side L = the sequence starting at pos is
joined at its left end. A breakpoint (locus) = the gap after base g of a contig (R at g, L at g + 1).
"""
import argparse, bisect, collections, gzip, json, os, re, sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(q): return q.translate(COMP)[::-1]
def die(msg): sys.exit(f"ERROR: {msg}")
def log(msg): print(msg, file=sys.stderr, flush=True)
def opener(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)
M_OPS, Q_OPS, R_OPS = {0, 7, 8}, {0, 1, 4, 7, 8}, {0, 2, 3, 7, 8}     # aligned / query-consuming / ref-consuming


# ============================================================================================== refs
def cmd_refs(a):
    """one pass over the assembly FASTA: the contigs of each MAG -> OUT/<n>/mag.fa (upper case) ; --contigs: contig TAB
    MAG TAB n (n = folder name)"""
    want = {}
    for l in open(a.contigs):
        c, m, n = l.rstrip("\n").split("\t")
        want[c] = n
    outs, cur, got = {}, None, set()
    with opener(a.fasta) as fh:
        for l in fh:
            if l[0] == ">":
                c = l[1:].split()[0]
                cur = None
                if c in want:
                    n = want[c]
                    if n not in outs: os.makedirs(os.path.join(a.out, n), exist_ok=True); outs[n] = open(os.path.join(a.out, n, "mag.fa"), "w")
                    cur = outs[n]; cur.write(f">{c}\n"); got.add(c)
            elif cur:
                cur.write(l.strip().upper() + "\n")
    for o in outs.values(): o.close()
    miss = set(want) - got
    if miss: die(f"{len(miss)} MAG contigs not in {a.fasta} (e.g. {next(iter(miss))})")


# ============================================================================================== junctions
class Piece:
    """an aligned stretch of a read without a large event inside: contig interval [r0, r1) (0-based), read interval
    [q0, q1) in the pseudo-read's own orientation, strand +1 / -1, ..."""
    __slots__ = ("c", "r0", "r1", "q0", "q1", "st", "mapq", "supp", "div", "aln", "syn", "rec", "all")
    def __init__(s, c, r0, r1, q0, q1, st, mapq, supp, div, aln, syn, rec, all_):
        s.c, s.r0, s.r1, s.q0, s.q1, s.st, s.mapq, s.supp, s.div, s.aln, s.syn, s.rec, s.all = \
            c, r0, r1, q0, q1, st, mapq, supp, div, aln, syn, rec, all_
    def exit_be(s): return (s.c, s.r1, "R") if s.st > 0 else (s.c, s.r0 + 1, "L")     # where the read leaves the piece
    def entry_be(s): return (s.c, s.r0 + 1, "L") if s.st > 0 else (s.c, s.r1, "R")    # where the read enters it


RE_BIG = re.compile(r"(\d{3,})[IDN]")


def record_pieces(r, ri, S, fasta, edlib, stats):
    """one alignment record -> pieces (split at CIGAR D / I >= min_sv and at dense difference regions whose read
    segment is the reverse complement of the reference = short inversions, which aligners pass straight through)"""
    cig = r.cigartuples
    L = r.infer_read_length()
    lc = (cig[0][1] if cig[0][0] in (4, 5) else 0) + (cig[1][1] if len(cig) > 1 and cig[0][0] == 5 and cig[1][0] == 4 else 0)
    hl = cig[0][1] if cig[0][0] == 5 else 0
    st, c = (-1 if r.is_reverse else 1), r.reference_name
    hard = cig[0][0] == 5 or cig[-1][0] == 5
    ok_all = not r.is_supplementary and r.mapping_quality >= S["min_mapq"] and not hard
    cst = r.get_cigar_stats()[0]
    if cst[8] + cst[1] + cst[2] < S["dense_n"] and not any(int(n) >= S["min_sv"] for n in RE_BIG.findall(r.cigarstring)):
        aln = int(cst[0] + cst[7] + cst[8])                     # fast path: no event possible in this record
        qa = r.query_alignment_length                            # = samtools qlen - sclen (the 'all' filter)
        ok_all = ok_all and qa * 100 >= S["min_pct_aligned"] * L and qa >= S["min_aln_bp"]
        q0, q1 = lc, L - ((cig[-1][1] if cig[-1][0] in (4, 5) else 0) + (cig[-2][1] if len(cig) > 1 and cig[-1][0] == 5 and cig[-2][0] == 4 else 0))
        if st < 0: q0, q1 = L - q1, L - q0
        return [Piece(c, r.reference_start, r.reference_end, q0, q1, st, r.mapping_quality, r.is_supplementary,
                      int(cst[8] + cst[1] + cst[2]), aln, False, ri, ok_all)], L
    # walk: segments between large events ; difference events for the dense-region scan
    segs, ev, rp, qp = [], [], r.reference_start, lc
    s_r, s_q, div, aln = rp, qp, 0, 0
    for op, n in cig:
        if op in (4, 5): continue
        if op in M_OPS:
            if op == 8: ev.append((rp, n)); div += n
            rp += n; qp += n; aln += n
        elif op == 1:
            if n >= S["min_sv"]:
                segs.append((s_r, rp, s_q, qp, div, aln, "ins", n)); s_r, s_q, div, aln = rp, qp + n, 0, 0
            else: ev.append((rp, n)); div += n
            qp += n
        elif op in (2, 3):
            if n >= S["min_sv"]:
                segs.append((s_r, rp, s_q, qp, div, aln, "del", n)); s_r, s_q, div, aln = rp + n, qp, 0, 0
            else: ev.append((rp, n)); div += n
            rp += n
    segs.append((s_r, rp, s_q, qp, div, aln, None, 0))
    qa = r.query_alignment_length                                # = samtools qlen - sclen (the 'all' filter)
    ok_all = ok_all and qa * 100 >= S["min_pct_aligned"] * L and qa >= S["min_aln_bp"]
    # dense regions -> inversion check -> split
    inv = []
    if len(ev) and sum(n for _, n in ev) >= S["dense_n"]:
        i = tot = 0
        marks = []
        for j, (p, n) in enumerate(ev):
            tot += n
            while p - ev[i][0] >= S["dense_w"]: tot -= ev[i][1]; i += 1
            if tot >= S["dense_n"]: marks.append([ev[i][0], p + max(1, n)])
        reg = []
        for m in marks:
            if reg and m[0] <= reg[-1][1]: reg[-1][1] = max(reg[-1][1], m[1])
            else: reg.append(m)
        for ra, rb in reg:
            if not 40 <= rb - ra <= 5000: continue
            qa, qb = ref2q(cig, r.reference_start, lc, ra), ref2q(cig, r.reference_start, lc, rb)
            if qa is None or qb is None or qb - qa < 30: continue
            seq = r.query_sequence
            rs = seq[qa - hl:qb - hl].upper() if seq else ""
            ref = fasta.fetch(c, ra, rb).upper()
            if not rs or not ref: continue
            d = edlib.align(rs, rc(ref), mode="NW", task="distance")["editDistance"]
            if 1 - d / max(len(rs), len(ref)) >= S["inv_id"]:
                inv.append((ra, rb, qa, qb)); stats["dense_inversions"] += 1
            else:
                stats["dense_other"] += 1
    out = []
    def add(r0, r1, q0, q1, dv, al, syn, s_):
        if r1 <= r0: return
        q0o, q1o = (q0, q1) if st > 0 else (L - q1, L - q0)        # -> the pseudo-read's own orientation
        out.append(Piece(c, r0, r1, q0o, q1o, s_, r.mapping_quality, r.is_supplementary, dv, al, syn, ri, ok_all))
    for s_r, e_r, s_q, e_q, dv, al, _, _ in segs:
        cuts = [x for x in inv if s_r <= x[0] and x[1] <= e_r]
        if not cuts: add(s_r, e_r, s_q, e_q, dv, al, False, st); continue
        pr, pq = s_r, s_q
        for ra, rb, qa, qb in cuts:
            frac = (ra - pr) / max(1, e_r - s_r)
            add(pr, ra, pq, qa, round(dv * frac), ra - pr, False, st)
            add(ra, rb, qa, qb, 0, rb - ra, True, -st)                   # the inverted segment: opposite strand
            pr, pq = rb, qb
        add(pr, e_r, pq, e_q, round(dv * (e_r - pr) / max(1, e_r - s_r)), e_r - pr, False, st)
    return out, L


def ref2q(cig, rstart, lc, target):
    """query position (alignment orientation, hard clips counted) of reference position `target` (0-based)"""
    rp, qp = rstart, lc
    for op, n in cig:
        if op in (4, 5): continue
        if op in M_OPS:
            if rp + n > target: return qp + (target - rp)
            rp += n; qp += n
        elif op == 1: qp += n
        elif op in (2, 3):
            if rp + n > target: return qp
            rp += n
    return qp if rp == target else None


def read_chain(recs, S, fasta, edlib, stats):
    """all records of one pseudo-read -> (pieces in read order after merging non-events, read length)"""
    ps, L = [], 0
    for i, r in enumerate(recs):
        p, L = record_pieces(r, i, S, fasta, edlib, stats)
        ps += p
    ps = [p for p in ps if p.mapq >= S["min_mapq"] and (p.syn or p.aln >= S["min_piece"])]
    ps.sort(key=lambda p: (p.q0, p.q1))
    out = []
    for p in ps:                     # same contig + strand, continuing within min_sv on both axes: one piece
        if out:
            q = out[-1]
            gr = (p.r0 - q.r1) if p.st > 0 else (q.r0 - p.r1)
            gq = p.q0 - q.q1
            if q.c == p.c and q.st == p.st and not q.syn and not p.syn and -S["min_sv"] < gr < S["min_sv"] and -S["min_sv"] < gq < S["min_sv"]:
                q.r0, q.r1, q.q1 = min(q.r0, p.r0), max(q.r1, p.r1), max(q.q1, p.q1)
                q.div += p.div + abs(gr) + max(0, gq); q.aln += p.aln; q.mapq = min(q.mapq, p.mapq); q.all = q.all and p.all
                continue
        out.append(p)
    return out, L


def observations(ps, L, S):
    """consecutive pieces -> (BE_a, BE_b | None, insert / overhang bp, ext_a, ext_b, syn_a, syn_b) ; ext = aligned bp of the
    piece next to that end (how far the read reaches into that side) ; syn = that piece is a short inversion the aligner
    passed through (its bases are aligned the wrong way round: never a methylation window)"""
    ob = []
    if not ps: return ob
    for p, n in zip(ps, ps[1:]):
        ob.append((p.exit_be(), n.entry_be(), max(0, n.q0 - p.q1), p.r1 - p.r0, n.r1 - n.r0, p.syn, n.syn))
    if ps[0].q0 >= S["overhang"]: ob.append((ps[0].entry_be(), None, ps[0].q0, ps[0].r1 - ps[0].r0, 0, ps[0].syn, False))
    if L - ps[-1].q1 >= S["overhang"]: ob.append((ps[-1].exit_be(), None, L - ps[-1].q1, ps[-1].r1 - ps[-1].r0, 0, ps[-1].syn, False))
    return ob


def molecules(path):
    """name-sorted BAM -> (molecule, {'fwd': [records], 'rev': [records]})"""
    import pysam
    cur, grp = None, {}
    with pysam.AlignmentFile(path, "rb", check_sq=False) as f:
        for r in f.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary: continue
            m, _, t = r.query_name.rpartition("/")
            if t not in ("fwd", "rev"): m, t = r.query_name, "fwd"
            if m != cur:
                if cur is not None: yield cur, grp
                cur, grp = m, {}
            grp.setdefault(t, []).append(r)
    if cur is not None: yield cur, grp


def be_locus(be):   # BE -> breakpoint (contig, g): the gap after base g
    c, p, sd = be
    return (c, p) if sd == "R" else (c, p - 1)


def cluster(vals, tol):
    """sorted positions -> list of clusters (lists of indices into vals) ; a gap > tol starts a new cluster"""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    out = []
    for i in order:
        if out and vals[i] - vals[out[-1][-1]] <= tol: out[-1].append(i)
        else: out.append([i])
    return out


def cmd_junctions(a):
    import edlib, pysam, numpy as np
    import pyarrow as pa, pyarrow.parquet as pq
    S = json.loads(a.settings)
    os.makedirs(a.out, exist_ok=True)
    mags = [l.rstrip("\n").split("\t") for l in open(a.mags) if l.strip()]      # MAG, MGS, ref fasta, name-sorted raw BAM, unaligned BAM
    stats = collections.Counter()
    obs, mol_mag, clen = [], {}, {}
    piece_cols, piece_tabs = collections.defaultdict(list), []
    mag_mgs = {m: g for m, g, *_ in mags}
    for mag, mgs, ref, raw, _ in mags:
        fasta = pysam.FastaFile(ref)
        for c, n in zip(fasta.references, fasta.lengths): clen[c] = n
        for mol, grp in molecules(raw):
            stats["molecules"] += 1
            chains = {}
            for t in ("fwd", "rev"):
                if t in grp:
                    stats["pseudo_reads"] += 1
                    chains[t] = read_chain(grp[t], S, fasta, edlib, stats)
            t = "fwd" if chains.get("fwd", ([], 0))[0] else "rev"
            ps, L = chains.get(t, ([], 0))
            if not ps: stats["molecules_without_piece"] += 1; continue
            if "fwd" in chains and "rev" in chains and chains["fwd"][0] and chains["rev"][0]:
                stats["both_strands"] += 1
                stats["both_strands_same_pieces"] += len(chains["fwd"][0]) == len(chains["rev"][0])
            mol_mag[mol] = mag
            for o in observations(ps, L, S): obs.append((mol,) + o)
            # pieces in molecule orientation (the fwd pseudo-read = the molecule ; a rev one is its reverse complement)
            seq = ps if t == "fwd" else ps[::-1]
            for k, p in enumerate(seq):
                q0, q1, s_ = (p.q0, p.q1, p.st) if t == "fwd" else (L - p.q1, L - p.q0, -p.st)
                for col, v in (("MAG", mag), ("molecule", mol), ("piece", k), ("contig", p.c), ("r_start", p.r0 + 1), ("r_end", p.r1),
                               ("strand", "+" if s_ > 0 else "-"), ("q_start", q0 + 1), ("q_end", q1), ("q_len", L), ("mapq", p.mapq),
                               ("supplementary", p.supp), ("synthetic", p.syn), ("div", p.div), ("in_all", p.all)):
                    piece_cols[col].append(v)
            # events at each piece's ends (contig orientation): junction (next / previous piece), overhang, end
            n = len(seq)
            for k, p in enumerate(seq):
                qs_ev = "junction" if k > 0 else ("overhang" if (piece_cols["q_start"][-n + k] - 1) >= S["overhang"] else "end")
                qe_ev = "junction" if k < n - 1 else ("overhang" if L - piece_cols["q_end"][-n + k] >= S["overhang"] else "end")
                s_ = 1 if piece_cols["strand"][-n + k] == "+" else -1
                piece_cols["ev_start"].append(qs_ev if s_ > 0 else qe_ev); piece_cols["ev_end"].append(qe_ev if s_ > 0 else qs_ev)
        fasta.close()
        if piece_cols: piece_tabs.append(pa.table(dict(piece_cols))); piece_cols.clear()
    stats["observations"] = len(obs)
    log(f"junctions: {stats['molecules']} molecules, {len(obs)} breakpoint observations")

    # ---------------------------------------------------------------- read sets
    sets, loci_bes = build_sets(obs, S, clen, stats)
    # REF: molecules spanning a breakpoint contiguously (>= anchor bp on both sides), at breakpoints with a large ALT set
    P = pa.concat_tables(piece_tabs).combine_chunks() if piece_tabs else None
    need = sorted({lc for s_ in sets if s_["n"] >= S["min_mol"] for lc in s_["loci"] if lc})
    ref_sets = ref_molecules(P, need, sets, S, clen, stats)
    allsets = finalize(sets, ref_sets, S, stats)
    # ---------------------------------------------------------------- windows
    wins = make_windows(allsets, S, clen, mags, stats)
    # window records (both pseudo-reads of every molecule used), ip orientation check, divergence per window
    write_window_bam(wins, mags, S, a.out, stats)
    # ---------------------------------------------------------------- outputs
    member = collections.defaultdict(list)
    for s_ in allsets:
        for m in s_["mols"]: member[m].append(s_["id"])
    if P is not None:
        mol = P.column("molecule").to_pylist()
        P = P.append_column("read_sets", pa.array([",".join(member.get(m, [])) for m in mol], pa.string()))
        P = P.append_column("MGS", pa.array([mag_mgs[m] for m in P.column("MAG").to_pylist()], pa.string()))
        pq.write_table(P, f"{a.out}/pieces.parquet", compression="zstd")
    else:
        open(f"{a.out}/pieces.none", "w").close()
    ctg_mag = {}
    for mag, mgs, ref, raw, _ in mags:
        for l in open(ref + ".fai"): ctg_mag[l.split("\t", 1)[0]] = mag
    with open(f"{a.out}/sets.tsv", "w") as o:
        cols = ["read_set", "side", "MGS", "MAG", "kind", "event", "contig", "pos", "bp_side", "partner_kind", "partner_contig",
                "partner_pos", "partner_side", "insert_bp", "event_bp", "locus", "n_molecules", "locus_molecules", "locus_ref", "frac", "called", "call_bp",
                "n_used", "div_kb", "flags", "molecules"]
        o.write("\t".join(cols) + "\n")
        for s_ in allsets:
            for row in set_rows(s_, mag_mgs, ctg_mag):
                o.write("\t".join("" if row.get(k) is None else str(row[k]) for k in cols) + "\n")
    with open(f"{a.out}/windows.tsv", "w") as o:
        o.write("wname\tread_set\tside\tMGS\tMAG\tcontig\tw_c_start\tw_c_end\tprefix\tbp\tbp_side\tD\n")
        for w in wins:
            o.write(f"{w['name']}\t{w['set']}\t{w['side']}\t{mag_mgs[w['MAG']]}\t{w['MAG']}\t{w['contig']}\t{w['a']}\t{w['b']}\t"
                    f"{len(w['pre'])}\t{w['bp']}\t{w['bp_side']}\t{w['D']}\n")
    with open(f"{a.out}/stats.tsv", "w") as o:
        for k in sorted(stats): o.write(f"junction_{k}\t{stats[k]}\n")


def build_sets(obs, S, clen, stats):
    """observations -> ALT read sets: BEs clustered per (contig, side) within tol ; a set = (BE cluster, partner cluster
    | unplaced, insert length class) with its molecules"""
    bes = []
    for i, (m, a_, b_, ins, ea, eb, sa, sb) in enumerate(obs):
        bes.append((a_, i, 0))
        if b_: bes.append((b_, i, 1))
    by = collections.defaultdict(list)
    for be, i, k in bes: by[(be[0], be[2])].append((be[1], i, k))
    cid = {}                                  # (obs index, end) -> cluster key (contig, pos, side)
    for (c, sd), v in by.items():
        for cl in cluster([x[0] for x in v], S["tol"]):
            pos = sorted(v[j][0] for j in cl)[len(cl) // 2]
            for j in cl: cid[(v[j][1], v[j][2])] = (c, pos, sd)
    groups = collections.defaultdict(list)
    for i, (m, a_, b_, ins, ea, eb, sa, sb) in enumerate(obs):
        ka, kb = cid[(i, 0)], (cid[(i, 1)] if b_ else None)
        if kb is not None and kb < ka: ka, kb, ea, eb, sa, sb = kb, ka, eb, ea, sb, sa
        groups[(ka, kb)].append((m, ins, ea, eb, sa, sb))
    sets = []
    for (ka, kb), v in groups.items():
        v.sort(key=lambda x: x[1])
        cls = [[v[0]]]
        for x in v[1:]:                       # insert length classes (a different element at the same place)
            if kb is not None and x[1] - cls[-1][-1][1] > max(S["tol"], 0.1 * cls[-1][-1][1]) and max(x[1], cls[-1][-1][1]) >= S["min_sv"]: cls.append([x])
            else: cls[-1].append(x)
        for cl in cls:
            mols, syn = {}, [0, 0]
            for m, ins, ea, eb, sa, sb in cl:  # one entry per molecule (its longest reach on each side)
                o = mols.get(m)
                mols[m] = (max(ins, o[0]), max(ea, o[1]), max(eb, o[2])) if o else (ins, ea, eb)
                syn[0] += sa; syn[1] += sb
            ev, size = classify(ka, kb, sorted(x[0] for x in mols.values())[len(mols) // 2], S, clen)
            if ev is None: stats["observations_small"] += len(mols); continue
            if len(mols) < 2 and ev not in ("circular", "contig_join"):
                stats["singleton_observations"] += 1; continue
            sets.append(dict(kind="alt", a=ka, b=kb, event=ev, size=size, ins=sorted(x[0] for x in mols.values())[len(mols) // 2],
                             mols=mols, n=len(mols), loci=[be_locus(ka), be_locus(kb) if kb else None],
                             inner=[syn[0] * 2 > len(cl), syn[1] * 2 > len(cl)]))
    # breakpoints within tol bp are one (an insertion's two ends overlap by its microhomology)
    lv = collections.defaultdict(set)
    for s_ in sets:
        for lc in s_["loci"]:
            if lc: lv[lc[0]].add(lc[1])
    rep_ = {}
    for c, v in lv.items():
        v = sorted(v)
        for cl in cluster(v, S["tol"]):
            m = sorted(v[i] for i in cl)[(len(cl) - 1) // 2]
            for i in cl: rep_[(c, v[i])] = (c, m)
    for s_ in sets: s_["loci"] = [rep_[lc] if lc else None for lc in s_["loci"]]
    return sets, None


def classify(ka, kb, ins, S, clen):
    """event of a set ; None = below min_sv (not a large rearrangement)"""
    if kb is None: return "unplaced", None
    (c1, p1, s1), (c2, p2, s2) = ka, kb
    endish = lambda c, p, s: (s == "R" and p >= clen[c] - S["end_tol"]) or (s == "L" and p <= S["end_tol"])
    if c1 == c2 and s1 != s2:
        (pr, _), (pl, _) = ((p1, s1), (p2, s2)) if s1 == "R" else ((p2, s2), (p1, s1))
        if pr >= clen[c1] - S["end_tol"] and pl <= S["end_tol"]: return "circular", None
        if pl > pr:
            d = pl - pr - 1
            if d >= S["min_sv"] and ins >= S["min_sv"]: return "replacement", max(d, ins)
            if d >= S["min_sv"]: return "deletion", d
            if ins >= S["min_sv"]: return "insertion", ins
            return None, None
        size = pr - pl + 1                          # the read goes back: a tandem duplication of [pl, pr] ...
        if size >= S["min_sv"]: return "duplication", size
        return ("insertion", ins) if ins >= S["min_sv"] else (None, None)   # ... or a few bp of microhomology at an insertion
    if c1 == c2:
        return "inversion", abs(p1 - p2)
    if endish(c1, p1, s1) and endish(c2, p2, s2): return "contig_join", None
    return "translocation", None


def ref_molecules(P, need, sets, S, clen, stats):
    """REF read set per breakpoint in `need`: molecules whose piece covers [g - anchor + 1, g + anchor] ; molecules that
    are also ALT at that breakpoint are dropped from both"""
    import numpy as np
    out = []
    if P is None or not need: return out
    alt_at = collections.defaultdict(set)
    for s_ in sets:
        for lc in s_["loci"]:
            if lc: alt_at[lc] |= set(s_["mols"])
    c_all = np.array(P.column("contig").to_pylist(), dtype=object)
    r0_all, r1_all = np.array(P.column("r_start").to_pylist()) - 1, np.array(P.column("r_end").to_pylist())
    mol_all = P.column("molecule").to_pylist()
    idx = collections.defaultdict(list)
    for i, c in enumerate(c_all): idx[c].append(i)
    A = S["anchor"]
    for c, g in need:
        if g <= 0 or g >= clen[c]: continue                      # contig end: no continuation possible
        ii = np.array(idx.get(c, []), dtype=np.int64)
        if not len(ii): continue
        hit = ii[(r0_all[ii] <= g - A) & (r1_all[ii] >= g + A)]
        mols = {}
        for i in hit:
            m = mol_all[i]
            if m in alt_at[(c, g)]: stats["ref_alt_conflicts"] += 1; continue
            mols[m] = (0, g - int(r0_all[i]), int(r1_all[i]) - g)
        for m in alt_at[(c, g)] & set(mols): mols.pop(m)
        out.append(dict(kind="ref", a=(c, g, "R"), b=(c, g + 1, "L"), event="ref", size=None, ins=0, mols=mols, n=len(mols),
                        loci=[(c, g), (c, g)], inner=[False, False]))
    return out


def finalize(sets, refs, S, stats):
    """ids, locus molecule counts, frac, called flags, links between sets, flags"""
    import itertools
    at = collections.defaultdict(set)                             # breakpoint -> every molecule REF or ALT there
    for s_ in sets + refs:
        for lc in s_["loci"]:
            if lc: at[lc] |= set(s_["mols"])
    ref_n = {s_["loci"][0]: s_["n"] for s_ in refs}
    for s_ in sets + refs:
        s_["locus_n"] = [len(at[lc]) if lc else None for lc in s_["loci"]]
        s_["locus_ref"] = [ref_n.get(lc) if lc else None for lc in s_["loci"]]   # None: REF not counted there
        fr = [s_["n"] / n for n in s_["locus_n"] if n]
        s_["frac"] = round(max(fr), 4) if fr else None
    # the reads deny the assembled conformation: at every breakpoint of the set, >= 90 % of the molecules are ALT
    # (any ALT set there: the two junctions of an inversion share their breakpoints) and fewer than min_mol are REF
    deny = lambda lc: lc in at and ref_n.get(lc, 0) < S["min_mol"] and len(at[lc]) - ref_n.get(lc, 0) >= 0.9 * len(at[lc])
    alt_loci = set()                                               # breakpoints with a real alternative conformation
    for s_ in sets:
        placed = s_["b"] is not None
        real = s_["n"] >= S["min_mol"] and (s_["frac"] or 0) >= S["min_frac"]
        s_["called"] = bool(S["junctions"]) and placed and real
        if real: alt_loci |= {lc for lc in s_["loci"] if lc}
        fl = []
        if s_["event"] not in ("circular", "contig_join", "unplaced") and s_["n"] >= S["min_mol"] and all(deny(lc) for lc in s_["loci"] if lc):
            fl.append("assembly_suspect")
        s_["flags"] = fl
    for s_ in refs:
        s_["called"] = bool(S["junctions"]) and s_["n"] >= S["min_mol"] and s_["loci"][0] in alt_loci
        s_["flags"] = []
    sets.sort(key=lambda s_: (s_["a"], s_["b"] or ("",), s_["ins"]))
    refs.sort(key=lambda s_: s_["a"])
    for i, s_ in enumerate(sets): s_["id"] = f"J{i + 1:04d}"
    for i, s_ in enumerate(refs): s_["id"] = f"R{i + 1:04d}"
    of = collections.defaultdict(list)                            # molecule -> sets (of >= 2 molecules)
    for i, s_ in enumerate(sets):
        if s_["n"] >= 2:
            for m in s_["mols"]: of[m].append(i)
    shared = collections.Counter(p for v in of.values() for p in itertools.combinations(v, 2))
    for (i, j), sh in sorted(shared.items()):                     # inversion pairs etc.: >= 80 % shared
        x, y = sets[i], sets[j]
        if sh >= 0.8 * min(x["n"], y["n"]):
            x["flags"].append(f"linked:{y['id']}"); y["flags"].append(f"linked:{x['id']}")
    for s_ in sets:     # a read folding back on itself (R-R / L-L a few hundred bp apart, no partner junction): artefact?
        if s_["event"] == "inversion" and abs(s_["a"][1] - s_["b"][1]) < 300 and not any(f.startswith("linked:") for f in s_["flags"]):
            s_["flags"].append("fold_back")
    stats["sets_alt"] = len(sets); stats["sets_ref"] = len(refs)
    stats["sets_called"] = sum(s_["called"] for s_ in sets + refs)
    return refs + sets


def make_windows(allsets, S, clen, mags, stats):
    """per called set and side: callable distance D (largest n(d) * d with n(d) >= min_mol molecules reaching d + M
    into that side), the molecules used, the window = flank + junction-aware continuation (K bases of what follows the
    breakpoint in that conformation: partner sequence ; an insert as N)"""
    import pysam
    M, K = S["ctx"] + 5, S["ctx"] + 5
    fa = {}
    for mag, mgs, ref, raw, un in mags:
        f = pysam.FastaFile(ref)
        for c in f.references: fa[c] = (ref, mag)
        f.close()
    open_ = {}
    def seq(c, a_, b_):            # 1-based inclusive, clipped to the contig
        a_, b_ = max(1, a_), min(clen[c], b_)
        if b_ < a_: return ""
        p = fa[c][0]
        if p not in open_: open_[p] = pysam.FastaFile(p)
        return open_[p].fetch(c, a_ - 1, b_).upper()
    def cont_after(be, ins):       # what the read reads after leaving through BE's partner side, K bases
        c, p, sd = be
        q = seq(c, p, p + K - 1) if sd == "L" else rc(seq(c, p - K + 1, p))
        return ("N" * min(ins, K) + q)[:K]
    wins, n = [], 0
    for s_ in allsets:
        s_["side"] = {}
        if not s_["called"]: continue
        sides = [(1, s_["a"], 1), (2, s_["b"], 2)]
        for side, be, k in sides:
            c, p, sd = be
            if s_["inner"][side - 1]:                 # the aligner passed this inverted stretch straight through
                s_["side"][side] = dict(D=None, mols=[]); s_["flags"].append(f"side{side}_not_aligned"); continue
            ext = sorted((v[k] for v in s_["mols"].values()), reverse=True)
            best = (0, 0)
            for i, e in enumerate(ext):            # n(d) = i + 1 molecules reach >= e ; d = e - M
                d = min(S["window"], e - M)
                if i + 1 >= S["min_mol"] and d > 0 and (i + 1) * d > best[0] * best[1]: best = (i + 1, d)
            if best[0] == 0: s_["side"][side] = dict(D=None, mols=[]); continue
            D = best[1]
            used = sorted(m for m, v in s_["mols"].items() if v[k] >= D + M)
            if s_["kind"] == "ref":
                partner = s_["b"] if side == 1 else s_["a"]
                ins = 0
            else:
                partner = s_["b"] if side == 1 else s_["a"]
                ins = s_["ins"]
            if sd == "R":
                a_, b_ = max(1, p - D - M + 1), p
                pre, suf = "", cont_after(partner, ins)
            else:
                a_, b_ = p, min(clen[c], p + D + M - 1)
                pre, suf = rc(cont_after(partner, ins)), ""
            n += 1
            w = dict(name=f"w{n:06d}", set=s_["id"], side=side, MAG=fa[c][1], contig=c, a=a_, b=b_, pre=pre, suf=suf, bp=p,
                     bp_side=sd, D=D, mols=used, seq=pre + seq(c, a_, b_) + suf)
            s_["side"][side] = dict(D=D, mols=used, w=w)
            wins.append(w)
    for f in open_.values(): f.close()
    stats["windows"] = len(wins)
    names = set(clen)
    if any(w["name"] in names for w in wins): die("a window name equals a MAG contig name")
    return wins


def clip_record(r, a_, b_, off, full):
    """aligned record -> its part on contig [a_, b_] (1-based) as a record on a window contig whose base 1 is contig
    base a_ - off ; parts outside -> soft clip ; hard clips restored from `full` (seq, qual, tags of the unaligned
    pseudo-read, native orientation) ; None when < 50 aligned bases remain. Returns (cigar, ref_start0, seq, qual, n_diff, n_aln)"""
    cig = r.cigartuples
    hl = cig[0][1] if cig[0][0] == 5 else 0
    hr = cig[-1][1] if cig[-1][0] == 5 else 0
    if hl or hr:
        if full is None: return None
        s0, q0 = full
        sq = s0 if not r.is_reverse else rc(s0)
        ql = (q0 if not r.is_reverse else q0[::-1]) if q0 is not None else None
        if r.query_sequence and sq[hl:len(sq) - hr if hr else None] != r.query_sequence: return None
    else:
        sq, ql = r.query_sequence, r.query_qualities
    ops = [(op, n) for op, n in cig if op != 5]
    if hl: ops = [(4, hl)] + ops
    if hr: ops = ops + [(4, hr)]
    lo, hi = a_ - 1, b_                                  # 0-based [lo, hi)
    new, rp, qp, start, nd, na = [], r.reference_start, 0, None, 0, 0
    def push(op, n):
        if n <= 0: return
        if new and new[-1][0] == op: new[-1] = (op, new[-1][1] + n)
        else: new.append((op, n))
    for op, n in ops:
        if op == 4: push(4, n); qp += n; continue
        if op in M_OPS:
            # split the run at lo / hi
            s_, e_ = rp, rp + n
            x0, x1 = max(s_, lo), min(e_, hi)
            if x1 <= x0:
                push(4, n)
            else:
                push(4, x0 - s_)
                if start is None: start = x0
                push(op, x1 - x0); na += x1 - x0
                if op == 8: nd += x1 - x0
                push(4, e_ - x1)
            rp += n; qp += n
        elif op == 1:
            inside = start is not None and lo < rp < hi and (new and new[-1][0] in M_OPS)
            push(1 if inside else 4, n); qp += n
            if inside: nd += n
        elif op in (2, 3):
            if start is not None and lo < rp and rp + n < hi and new and new[-1][0] in M_OPS:
                push(2, n); nd += n
            rp += n
    # trailing / leading I or D next to a soft clip -> soft clip / drop
    def tidy(cl):
        out = []
        for op, n in cl:
            if op == 1 and (not out or out[-1][0] == 4): op = 4
            if op == 2 and (not out or out[-1][0] == 4): continue
            if out and out[-1][0] == op: out[-1] = (op, out[-1][1] + n)
            else: out.append((op, n))
        return out
    new = tidy(new)
    new = tidy(new[::-1])[::-1]
    new = [(4 if op == 1 else op, n) if k in (0, len(new) - 1) and op == 1 else (op, n) for k, (op, n) in enumerate(new)]
    if na < 50 or start is None: return None
    return new, start - (a_ - 1) + off, sq, ql, nd, na


def write_window_bam(wins, mags, S, out, stats):
    """window records: both pseudo-reads of every molecule a window uses, one record each (the one covering the anchor
    base next to the breakpoint, most overlap with the flank), clipped to the window, on the window contig.
    Hard-clipped records (supplementary) are rebuilt from the unaligned pseudo-read: full sequence, kinetics in native
    orientation (PacBio convention - checked here on reverse-strand primary records)."""
    import pysam
    need = collections.defaultdict(list)                  # molecule -> windows
    for w in wins:
        for m in w["mols"]: need[m].append(w)
    by_mag = collections.defaultdict(set)
    for w in wins: by_mag[w["MAG"]] |= set(w["mols"])
    raws = [(mag, raw, un) for mag, mgs, ref, raw, un in mags]

    def chosen(recs, w):
        c, p = w["contig"], w["bp"]
        cands = [r for r in recs if r.reference_name == c and r.reference_start < p <= r.reference_end and r.mapping_quality >= S["min_mapq"]]
        if not cands: return None
        return max(cands, key=lambda r: (min(r.reference_end, w["b"]) - max(r.reference_start, w["a"] - 1), -r.is_supplementary))
    hard = lambda r: r.cigartuples[0][0] == 5 or r.cigartuples[-1][0] == 5
    # pass 1: which chosen records are hard-clipped (their full read is needed) ; reverse primaries for the ip check
    full_need, chk, hdr = set(), {}, None
    for mag, raw, un in raws:
        with pysam.AlignmentFile(raw, "rb", check_sq=False) as f:
            if hdr is None: hdr = f.header.to_dict()
            for mol, grp in molecules(raw):
                ws = need.get(mol)
                for t, recs in grp.items():
                    if len(chk) < 100:
                        for r in recs:
                            if r.is_reverse and not r.is_supplementary and r.has_tag("ip") and not hard(r): chk[r.query_name] = list(r.get_tag("ip"))
                    for w in ws or ():
                        r = chosen(recs, w)
                        if r is not None and hard(r): full_need.add(r.query_name)
    # pass 2: unaligned pseudo-reads
    full, ok_, bad = {}, 0, 0
    for mag, raw, un in raws:
        with pysam.AlignmentFile(un, "rb", check_sq=False) as f:
            for r in f.fetch(until_eof=True):
                q = r.query_name
                if q in chk:
                    ip = list(r.get_tag("ip")) if r.has_tag("ip") else None
                    if ip == chk[q]: ok_ += 1
                    elif ip is not None and ip[::-1] == chk[q]: bad += 1
                if q in full_need: full[q] = (r.query_sequence, r.query_qualities, r.get_tags())
    stats["ip_orientation_checked"] = ok_ + bad
    if bad: die(f"{bad} of {ok_ + bad} reverse-strand records carry their kinetics reversed (expected native orientation, "
                "as in PacBio BAM) - window records would be wrong")
    if len(full) < len(full_need): stats["window_records_without_full_read"] = len(full_need) - len(full)
    hdr["SQ"] = [{"SN": w["name"], "LN": len(w["seq"])} for w in wins] or [{"SN": "w000000", "LN": 1}]
    if "HD" in hdr: hdr["HD"]["SO"] = "unknown"
    tid = {w["name"]: i for i, w in enumerate(wins)}
    with open(f"{out}/windows.fa", "w") as o:
        for w in wins: o.write(f">{w['name']}\n{w['seq']}\n")
    drop = {"SA", "MD", "NM", "ms", "AS", "nn", "tp", "cm", "s1", "s2", "de", "rl", "cg", "cs"}
    n_rec, n_chk = collections.Counter(), 0
    # pass 3: the window records
    with pysam.AlignmentFile(f"{out}/win.bam", "wb", header=hdr) as o:
        for mag, raw, un in raws:
            if mag not in by_mag: continue
            for mol, grp in molecules(raw):
                ws = need.get(mol)
                if not ws: continue
                for t, recs in grp.items():
                    for w in ws:
                        r = chosen(recs, w)
                        if r is None: continue
                        fu = full.get(r.query_name) if hard(r) else None
                        if hard(r) and fu is None: continue
                        x = clip_record(r, w["a"], w["b"], len(w["pre"]), (fu[0], fu[1]) if fu else None)
                        if x is None: continue
                        cig, start, sq, ql, nd, na = x
                        a = pysam.AlignedSegment(o.header)
                        a.query_name, a.flag = r.query_name, (16 if r.is_reverse else 0)
                        a.reference_id, a.reference_start, a.mapping_quality = tid[w["name"]], start, r.mapping_quality
                        a.cigartuples, a.query_sequence = cig, sq
                        if ql is not None: a.query_qualities = ql
                        a.set_tags([(k, v) for k, v in (fu[2] if fu else r.get_tags()) if k not in drop])
                        if fu is not None or n_chk < 1000:
                            n_chk += fu is None
                            if not check_pairs(r, a, w): die(f"window record of {r.query_name} does not follow its alignment (bug)")
                        o.write(a)
                        n_rec[(w["set"], w["side"])] += 1
                        w["nd"] = w.get("nd", 0) + nd; w["na"] = w.get("na", 0) + na
    stats["window_records"] = sum(n_rec.values())
    for w in wins: w["n_records"] = n_rec[(w["set"], w["side"])]


def check_pairs(r, a, w):
    """every aligned (read base, window base) pair of the new record is an aligned pair of the original"""
    off = w["a"] - 1 - len(w["pre"])
    orig = set(r.get_aligned_pairs(matches_only=True))
    hl = r.cigartuples[0][1] if r.cigartuples[0][0] == 5 else 0
    return all((q - hl, p + off) in orig for q, p in a.get_aligned_pairs(matches_only=True))


def set_rows(s_, mag_mgs, ctg_mag):
    """junction table rows (one per side) of a read set, before lift"""
    rows = []
    mag = ctg_mag[s_["a"][0]]
    sides = [(1, s_["a"], s_["b"])] + ([(2, s_["b"], s_["a"])] if s_["b"] is not None else [])
    for side, be, pt in sides:
        info = s_.get("side", {}).get(side, {})
        w = info.get("w")
        rows.append(dict(read_set=s_["id"], side=side, MGS=mag_mgs.get(mag), MAG=mag, kind=s_["kind"], event=s_["event"],
                         contig=be[0], pos=be[1], bp_side=be[2],
                         partner_kind=("contig" if s_["kind"] == "ref" else "placed") if pt is not None else "unplaced",
                         partner_contig=pt[0] if pt else None, partner_pos=pt[1] if pt else None, partner_side=pt[2] if pt else None,
                         insert_bp=s_["ins"], event_bp=s_["size"],
                         locus=s_["loci"][side - 1][1] if len(s_["loci"]) >= side and s_["loci"][side - 1] else None, n_molecules=s_["n"],
                         locus_molecules=s_["locus_n"][side - 1] if len(s_["locus_n"]) >= side else None,
                         locus_ref=s_["locus_ref"][side - 1] if len(s_["locus_ref"]) >= side else None,
                         frac=s_["frac"], called=int(bool(s_["called"] and w)), call_bp=info.get("D"), n_used=len(info.get("mols", [])),
                         div_kb=round(1000 * w["nd"] / w["na"], 3) if w and w.get("na") else None, flags=",".join(s_["flags"]),
                         molecules=",".join(info["mols"] if w else (sorted(s_["mols"]) if side == 1 else []))))
    return rows


# ============================================================================================== write / relift
LIFT_COLS = ["block", "instance", "block_pos", "block_strand", "ctx_dist", "gene", "up_gene", "up_dist"]
PROMOTER = 300        # up_gene: a gene start at most this far downstream (in the gene's orientation)


class Lift:
    """derep-ctg layer of one MGS for one sample's MAGs: lift segments, block sequences, genes, instances"""
    def __init__(s, d, mags, fp=None):
        import numpy as np, pyarrow.parquet as pq
        s.np, s.ok = np, os.path.isfile(f"{d}/lift.parquet")
        mags = sorted(mags)
        rfp = lambda: (open(f"{d}/.fp").read().strip() if os.path.isfile(f"{d}/.fp") else None)
        if fp and rfp() != fp: die(f"the derep-ctg result in {d} changed since this run started (rerun)")
        if mags:
            if not s.ok: die(f"no derep-ctg result in {d} (lift.parquet) - run mmt derep-ctg for it")
            have = set(pq.read_table(f"{d}/mags.parquet", columns=["MAG"]).column("MAG").to_pylist())
            if set(mags) - have:
                die(f"the derep-ctg result in {d} lacks MAG(s) {' '.join(sorted(set(mags) - have)[:3])} - rerun mmt derep-ctg for it")
        rd = lambda f, cols: pq.read_table(f"{d}/{f}", columns=cols, filters=[("MAG", "in", mags)]).to_pydict() if s.ok and mags else {k: [] for k in cols}
        L = rd("lift.parquet", ["contig", "c_start", "c_end", "b_start", "ori", "instance", "block"])
        G = rd("genes.parquet", ["contig", "gene", "c_start", "c_end", "strand"])
        I = rd("instances.parquet", ["instance", "MAG", "block", "contig", "c_start", "c_end", "strand", "host_instance", "nb5_key", "nb3_key",
                                     "nb5_gap", "nb3_gap"])
        s.I = I
        B = pq.read_table(f"{d}/blocks.parquet", columns=["block", "seq"]).to_pydict() if s.ok else {"block": [], "seq": []}
        s.bidx = {b: i for i, b in enumerate(B["block"])}
        s.seq = np.frombuffer("".join(B["seq"]).encode(), np.uint8)
        s.boff = np.cumsum([0] + [len(q) for q in B["seq"]])[:-1].astype(np.int64)
        s.blocks = np.array(B["block"], dtype=object)
        s.code = {}
        for c in L["contig"] + G["contig"] + I["contig"]: s.code.setdefault(c, len(s.code))
        key = lambda cs, ps: (np.array([s.code[c] for c in cs], np.int64) << 32) + np.array(ps, np.int64)
        o = np.argsort(key(L["contig"], L["c_start"]), kind="stable")
        s.lk = key(L["contig"], L["c_start"])[o]
        s.lc = np.array([s.code[c] for c in L["contig"]], np.int64)[o]
        s.le, s.lb = np.array(L["c_end"], np.int64)[o], np.array(L["b_start"], np.int64)[o]
        s.lo, s.li = np.array(L["ori"], np.int64)[o], np.array(L["instance"], np.int64)[o]
        s.lbl = np.array([s.bidx[b] for b in L["block"]], np.int64)[o]
        o = np.argsort(key(G["contig"], G["c_start"]), kind="stable")
        s.gk, s.gc = key(G["contig"], G["c_start"])[o], np.array([s.code[c] for c in G["contig"]], np.int64)[o]
        s.ge, s.gn = np.array(G["c_end"], np.int64)[o], np.array(G["gene"], dtype=object)[o]
        gs = np.array(G["strand"], dtype=object)[o]
        plus, minus = gs == "+", gs == "-"
        s.pk, s.pn = s.gk[plus], s.gn[plus]                                   # + genes by start
        mk = (s.gc[minus] << 32) + s.ge[minus]; om = np.argsort(mk, kind="stable")
        s.mk, s.mn = mk[om], s.gn[minus][om]                                   # - genes by end
        if fp and rfp() != fp: die(f"the derep-ctg result in {d} was replaced while it was read (rerun)")

    def rows(s, t):
        """lift columns of the rows of arrow table t (contig, pos, strand, base) -> dict of arrow arrays ; stops when a
        lifted base differs from its block base (exact-match segments: it never does unless positions are wrong)"""
        import pyarrow as pa, pyarrow.compute as pc
        np, n = s.np, t.num_rows
        typ = dict(block=pa.string(), instance=pa.int32(), block_pos=pa.int32(), block_strand=pa.string(), ctx_dist=pa.uint8(),
                   gene=pa.string(), up_gene=pa.string(), up_dist=pa.int16())
        out = {k: pa.nulls(n, typ[k]) for k in LIFT_COLS}
        if not n or not s.code: return out
        d = pc.dictionary_encode(t.column("contig")).combine_chunks()
        code = np.array([s.code.get(c, -1) for c in d.dictionary.to_pylist()], np.int64)[d.indices.to_numpy(zero_copy_only=False)]
        pos = t.column("pos").to_numpy().astype(np.int64)
        plus = pc.equal(t.column("strand"), "+").to_numpy(zero_copy_only=False)
        isA = pc.equal(t.column("base"), "A").to_numpy(zero_copy_only=False)
        key = (np.maximum(code, 0) << 32) + pos
        dict_arr = lambda idx, ok, names: pa.DictionaryArray.from_arrays(pa.array(np.where(ok, idx, 0), mask=~ok), pa.array(names)).cast(pa.string())
        if len(s.lk):
            i = np.searchsorted(s.lk, key, side="right") - 1
            j = np.maximum(i, 0)
            ok = (i >= 0) & (code >= 0) & (s.lc[j] == code) & (pos <= s.le[j])
            c0 = s.lk[j] & 0xFFFFFFFF
            ori = s.lo[j]
            bp = np.where(ori > 0, s.lb[j] + (pos - c0), s.lb[j] + (s.le[j] - pos))
            bplus = np.where(ori > 0, plus, ~plus)
            if ok.any():
                b = s.seq[np.where(ok, s.boff[s.lbl[j]] + bp - 1, 0)]
                want = np.where(bplus, np.where(isA, 65, 67), np.where(isA, 84, 71))      # A/C on the block strand
                bad = int((ok & (b != want)).sum())
                if bad: die(f"{bad} of {int(ok.sum())} lifted rows have another base than their block (positions shifted? bug)")
            out["block"] = dict_arr(s.lbl[j], ok, list(s.blocks) or [""])
            out["instance"] = pa.array(s.li[j].astype(np.int32), mask=~ok)
            out["block_pos"] = pa.array(bp.astype(np.int32), mask=~ok)
            out["block_strand"] = dict_arr(np.where(bplus, 0, 1), ok, ["+", "-"])
            out["ctx_dist"] = pa.array(np.minimum(np.minimum(pos - c0, s.le[j] - pos), 255).clip(0).astype(np.uint8), mask=~ok)
        if len(s.gk):
            g = np.searchsorted(s.gk, key, side="right") - 1
            hit, gi = np.zeros(n, bool), np.zeros(n, np.int64)
            for back in (0, 1, 2):                             # overlapping genes: the latest start holding pos
                h = np.maximum(g - back, 0)
                m = (~hit) & (g - back >= 0) & (code >= 0) & (s.gc[h] == code) & (pos <= s.ge[h])
                gi[m] = h[m]; hit |= m
            out["gene"] = dict_arr(gi, hit, list(s.gn))
            big = np.int64(10 ** 9)
            dp, dm = np.full(n, big), np.full(n, big)
            u = v = np.zeros(n, np.int64)
            if len(s.pk):                                      # + gene starting 1..PROMOTER bp downstream
                u = np.minimum(np.searchsorted(s.pk, key, side="right"), len(s.pk) - 1)
                okp = (code >= 0) & ((s.pk[u] >> 32) == code) & ((s.pk[u] & 0xFFFFFFFF) > pos)
                dp = np.where(okp, (s.pk[u] & 0xFFFFFFFF) - pos, big)
            if len(s.mk):                                      # - gene ending 1..PROMOTER bp upstream (its start is its end)
                v = np.maximum(np.searchsorted(s.mk, key, side="left") - 1, 0)
                okm = (code >= 0) & ((s.mk[v] >> 32) == code) & ((s.mk[v] & 0xFFFFFFFF) < pos)
                dm = np.where(okm, pos - (s.mk[v] & 0xFFFFFFFF), big)
            usep = (dp <= dm) & (dp <= PROMOTER)
            usem = ~usep & (dm <= PROMOTER)
            names = list(s.pn) + list(s.mn)
            idx = np.where(usep, u, np.where(usem, len(s.pn) + v, 0))
            out["up_gene"] = dict_arr(idx, usep | usem, names or [""])
            out["up_dist"] = pa.array(np.where(usep, dp, np.where(usem, dm, 0)).astype(np.int16), mask=~(usep | usem))
        return out

    def anchor(s, contig, pos, side):
        """the instance a breakpoint end's anchor runs into: the smallest one holding pos (an island inside its host
        included), else the nearest one on the anchor side (a top-level one on a tie) -> (instance, block, block_pos of
        pos on that instance (NULL off its lift), block end 5|3 facing the breakpoint, bp from the breakpoint to that end:
        > 0 = the breakpoint lies outside the instance, <= 0 = inside)"""
        if contig is None or not s.ok: return (None,) * 5
        pos, best = int(pos), None
        for k in range(len(s.I["instance"])):
            if s.I["contig"][k] != contig: continue
            a_, b_ = s.I["c_start"][k], s.I["c_end"][k]
            if side == "R" and a_ <= pos: d = pos - b_
            elif side == "L" and b_ >= pos: d = a_ - pos
            else: continue
            key = (0, b_ - a_, 0) if d <= 0 else (1, d, s.I["host_instance"][k] not in (-1, None))
            if best is None or key < best[0]: best = (key, k, d)
        if best is None: return (None,) * 5
        _, k, d = best
        plus = s.I["strand"][k] == "+"
        end = (3 if plus else 5) if side == "R" else (5 if plus else 3)
        return s.I["instance"][k], s.I["block"][k], s.block_pos(contig, pos, s.I["instance"][k]), end, d

    def block_pos(s, contig, pos, inst=None):
        """position on the block of the lift segment holding contig:pos (None off the lift, or on another instance)"""
        np, pos = s.np, int(pos)
        c = s.code.get(contig)
        if c is None or not len(s.lk): return None
        i = int(np.searchsorted(s.lk, (c << 32) + pos, side="right")) - 1
        if i < 0 or s.lc[i] != c or pos > s.le[i] or (inst is not None and s.li[i] != inst): return None
        c0 = int(s.lk[i] & 0xFFFFFFFF)
        return int(s.lb[i] + (pos - c0) if s.lo[i] > 0 else s.lb[i] + (s.le[i] - pos))


def sites_table(t, L, kin_fp, lift_fp):
    """kinetics rows (+ read_set, side, dist) -> final sites table: lift columns, fingerprints, sorted: 'all' rows on
    blocks by block, instance, block_pos, block strand ; 'all' rows off blocks by MAG, contig, pos, strand ; read sets by
    read_set, side, dist, strand"""
    import numpy as np, pyarrow as pa, pyarrow.compute as pc, mmtpq
    n = t.num_rows
    sch = mmtpq.schema("ctg_sites")
    lf = L.rows(t)
    cols = {c: t.column(c) for c in t.column_names}
    for k in LIFT_COLS: cols[k] = lf[k]
    cols["kin_fp"] = pa.array([kin_fp] * n, pa.string()); cols["lift_fp"] = pa.array([lift_fp] * n, pa.string())
    tt = pa.table({f.name: cols[f.name] for f in sch}).cast(sch)
    if not n: return tt
    def codes(col):
        d = pc.dictionary_encode(pc.fill_null(col, "")).combine_chunks()
        rank = np.argsort(np.argsort(np.array(d.dictionary.to_pylist(), dtype=object)))
        return rank[d.indices.to_numpy(zero_copy_only=False)]
    num = lambda col: pc.fill_null(col.cast(pa.int64()), -1).to_numpy(zero_copy_only=False)
    al = pc.equal(tt.column("read_set"), "all").to_numpy(zero_copy_only=False)
    lifted = pc.is_valid(tt.column("block")).to_numpy(zero_copy_only=False)
    grp = np.where(al, np.where(lifted, 0, 1), 2)
    c1 = np.where(grp == 0, codes(tt.column("block")), np.where(grp == 1, codes(tt.column("MAG")), codes(tt.column("read_set"))))
    c2 = np.where(grp == 0, num(tt.column("instance")), np.where(grp == 1, codes(tt.column("contig")), num(tt.column("side"))))
    c3 = np.where(grp == 0, num(tt.column("block_pos")), np.where(grp == 1, num(tt.column("pos")), num(tt.column("dist"))))
    st = np.where(grp == 0, pc.equal(pc.fill_null(tt.column("block_strand"), "+"), "-").to_numpy(zero_copy_only=False),
                  pc.equal(tt.column("strand"), "-").to_numpy(zero_copy_only=False)).astype(np.int8)
    return tt.take(pa.array(np.lexsort((st, c3, c2, c1, grp))))


def junction_table(sets, L, sites, kin_fp, lift_fp, sample):
    """sets.tsv rows of one MGS -> junction table with anchors on blocks and the number of site rows per read set side"""
    import pyarrow as pa, pyarrow.compute as pc
    cnt = collections.Counter()
    if sites is not None and sites.num_rows:
        x = sites.select(["read_set", "side"]).filter(pc.not_equal(sites.column("read_set"), "all"))
        if x.num_rows:
            g = x.group_by(["read_set", "side"]).aggregate([("read_set", "count")])
            cnt.update({(r_, s_): n for r_, s_, n in zip(*(g.column(k).to_pylist() for k in ("read_set", "side", "read_set_count")))})
    cols = collections.defaultdict(list)
    I = lambda x: None if x in ("", None) else int(x)
    F = lambda x: None if x in ("", None) else float(x)
    for r in sets:
        an = L.anchor(r["contig"], int(r["pos"]), r["bp_side"])
        pk = r["partner_contig"] or None
        pa_ = L.anchor(pk, int(r["partner_pos"]), r["partner_side"]) if pk else (None,) * 5
        row = dict(MGS=r["MGS"], sample=sample, MAG=r["MAG"], read_set=r["read_set"], side=int(r["side"]), kind=r["kind"], event=r["event"],
                   contig=r["contig"], pos=int(r["pos"]), bp_side=r["bp_side"], partner_kind=r["partner_kind"], partner_contig=pk,
                   partner_pos=I(r["partner_pos"]), partner_side=r["partner_side"] or None, insert_bp=I(r["insert_bp"]), event_bp=I(r["event_bp"]),
                   locus=I(r.get("locus")),
                   n_molecules=int(r["n_molecules"]), locus_molecules=I(r["locus_molecules"]), locus_ref=I(r.get("locus_ref")),
                   frac=F(r["frac"]), called=r["called"] == "1",
                   call_bp=I(r["call_bp"]), n_used=int(r["n_used"]), n_rows=cnt[(r["read_set"], int(r["side"]))], div_kb=F(r["div_kb"]),
                   flags=r["flags"] or None, molecules=[m for m in r["molecules"].split(",") if m],
                   anchor_instance=an[0], anchor_block=an[1], anchor_block_pos=an[2], anchor_end=an[3], anchor_end_dist=an[4],
                   partner_instance=pa_[0], partner_block=pa_[1], partner_block_pos=pa_[2], partner_end=pa_[3], partner_end_dist=pa_[4],
                   kin_fp=kin_fp, lift_fp=lift_fp)
        for k, v in row.items(): cols[k].append(v)
    return pa.table({f.name: pa.array(cols.get(f.name, []), type=f.type) for f in JUNCTION_SCHEMA()})


def JUNCTION_SCHEMA():
    import pyarrow as pa
    S_, I8, I16, I32, F, B = pa.string(), pa.int8(), pa.int16(), pa.int32(), pa.float32(), pa.bool_()
    return pa.schema([("MGS", S_), ("sample", S_), ("MAG", S_), ("read_set", S_), ("side", I8), ("kind", S_), ("event", S_),
                      ("contig", S_), ("pos", I32), ("bp_side", S_), ("partner_kind", S_), ("partner_contig", S_),
                      ("partner_pos", I32), ("partner_side", S_), ("insert_bp", I32), ("event_bp", I32), ("locus", I32), ("n_molecules", I32),
                      ("locus_molecules", I32), ("locus_ref", I32), ("frac", F), ("called", B), ("call_bp", I32), ("n_used", I32), ("n_rows", I32),
                      ("div_kb", F), ("flags", S_), ("molecules", pa.list_(S_)), ("anchor_instance", I32), ("anchor_block", S_),
                      ("anchor_block_pos", I32), ("anchor_end", I8), ("anchor_end_dist", I32), ("partner_instance", I32),
                      ("partner_block", S_), ("partner_block_pos", I32), ("partner_end", I8), ("partner_end_dist", I32),
                      ("kin_fp", S_), ("lift_fp", S_)])


def READS_SCHEMA():
    import pyarrow as pa
    S_ = pa.string()
    return pa.schema([("MGS", S_), ("sample", S_), ("MAG", S_), ("molecule", S_), ("piece", pa.int16()), ("contig", S_),
                      ("r_start", pa.int32()), ("r_end", pa.int32()), ("strand", S_), ("q_start", pa.int32()), ("q_end", pa.int32()),
                      ("q_len", pa.int32()), ("mapq", pa.uint8()), ("supplementary", pa.bool_()), ("synthetic", pa.bool_()),
                      ("ev_start", S_), ("ev_end", S_), ("div", pa.int32()), ("in_all", pa.bool_()), ("read_sets", pa.list_(S_))])


def EDGES_SCHEMA():
    import pyarrow as pa
    S_, I32 = pa.string(), pa.int32()
    return pa.schema([("MGS", S_), ("sample", S_), ("MAG", S_), ("kind", S_), ("contig", S_), ("pos", I32), ("bp_side", S_),
                      ("instance", I32), ("block", S_), ("block_end", pa.int8()), ("nb_key", S_), ("gap", I32), ("n_span", I32),
                      ("n_event", I32), ("n_runoff", I32), ("n_end", I32), ("lift_fp", S_)])


def edges_table(reads, junc, L, mgs, sample, lift_fp, S):
    """every instance end and contig end of the sample's MAGs: molecules spanning it contiguously (>= anchor bp each
    side), molecules with a breakpoint in its gap (listed read sets facing the same way), run-off / ending reads at
    contig ends"""
    import numpy as np, pyarrow as pa, pyarrow.compute as pc
    A, T, E = S["anchor"], S["tol"], S["end_tol"]
    byc, mset = {}, collections.defaultdict(set)          # contig -> (r_start, r_end, molecule code) ; read set -> molecules
    if reads.num_rows:
        mol = pc.dictionary_encode(reads.column("molecule")).combine_chunks()
        code = mol.indices.to_numpy(zero_copy_only=False)
        rs_, re_ = (reads.column(k).to_numpy().astype(np.int64) for k in ("r_start", "r_end"))
        ctg = pc.dictionary_encode(reads.column("contig")).combine_chunks()
        ci, names = ctg.indices.to_numpy(zero_copy_only=False), ctg.dictionary.to_pylist()
        for j, c in enumerate(names):
            w = ci == j
            byc[c] = (rs_[w], re_[w], code[w])
        for m, v in zip(mol.indices.to_pylist(), reads.column("read_sets").to_pylist()):
            for x in v or (): mset[x].add(m)
    J = collections.defaultdict(list)                      # (contig, bp_side) -> [(pos, read set)] of ALT rows
    if junc.num_rows:
        for c, p, sd, k, r_ in zip(*(junc.column(x).to_pylist() for x in ("contig", "pos", "bp_side", "kind", "read_set"))):
            if k == "alt": J[(c, sd)].append((p, r_))
    cols = collections.defaultdict(list)
    def add(**kw):
        for f in EDGES_SCHEMA():
            cols[f.name].append(kw.get(f.name))
    def span(c, lo, hi):          # molecules with one piece covering [lo, hi]
        if c not in byc: return 0
        a0, a1, m = byc[c]
        return int(np.unique(m[(a0 <= lo) & (a1 >= hi)]).size)
    I = L.I
    for k in range(len(I["instance"])):
        c, a_, b_, plus = I["contig"][k], I["c_start"][k], I["c_end"][k], I["strand"][k] == "+"
        for side, e in (("R", b_), ("L", a_)):
            five = (side == "L") == plus
            key, gap = (I["nb5_key"][k], I["nb5_gap"][k]) if five else (I["nb3_key"][k], I["nb3_gap"][k])
            g = gap if gap is not None else 0
            lo, hi = (e - A + 1, e + A) if side == "R" else (e - A, e + A - 1)
            lo_, hi_ = (min(e, e + g), max(e, e + g)) if side == "R" else (min(e - g, e), max(e - g, e))   # g < 0: overlap
            hit = {r_ for p, r_ in J.get((c, side), ()) if lo_ - T <= p <= hi_ + T}
            ev = len(set().union(*(mset[r_] for r_ in hit))) if hit else 0
            add(MGS=mgs, sample=sample, MAG=I["MAG"][k], kind="instance_end", contig=c, pos=e, bp_side=side, instance=I["instance"][k],
                block=I["block"][k], block_end=5 if five else 3, nb_key=key, gap=gap, n_span=span(c, lo, hi), n_event=ev, lift_fp=lift_fp)
    return pa.table({f.name: pa.array(cols.get(f.name, []), type=f.type) for f in EDGES_SCHEMA()})


def contig_end_rows(reads, clen, mag_of, mgs, sample, lift_fp, S):
    """contig ends: reads running off them (unaligned >= overhang bp, or continuing on another piece) / ending there"""
    import pyarrow as pa
    T, E = S["tol"], S["end_tol"]
    out = collections.defaultdict(lambda: [0, 0])
    if reads.num_rows:
        for c, a_, b_, es, ee in zip(*(reads.column(x).to_pylist() for x in ("contig", "r_start", "r_end", "ev_start", "ev_end"))):
            if c not in clen: continue
            if a_ <= E:
                if es in ("overhang", "junction"): out[(c, "L")][0] += 1
                elif a_ <= T + 1 and es == "end": out[(c, "L")][1] += 1
            if b_ >= clen[c] - E:
                if ee in ("overhang", "junction"): out[(c, "R")][0] += 1
                elif b_ >= clen[c] - T and ee == "end": out[(c, "R")][1] += 1
    cols = collections.defaultdict(list)
    for c in sorted(clen):
        for sd, pos in (("L", 1), ("R", clen[c])):
            ro, en = out[(c, sd)]
            for f, v in (("MGS", mgs), ("sample", sample), ("MAG", mag_of[c]), ("kind", "contig_end"), ("contig", c), ("pos", pos),
                         ("bp_side", "L" if sd == "L" else "R"), ("instance", None), ("block", None), ("block_end", None), ("nb_key", None),
                         ("gap", None), ("n_span", None), ("n_event", None), ("n_runoff", ro), ("n_end", en), ("lift_fp", lift_fp)):
                cols[f].append(v)
    return pa.table({f.name: pa.array(cols.get(f.name, []), type=f.type) for f in EDGES_SCHEMA()})


def read_sets_tsv(path):
    if not os.path.isfile(path): return []
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(h, l.rstrip("\n").split("\t"))) for l in fh]


def reads_table(P, mgs, sample):
    """pieces of one MGS -> reads table (read_sets: list of the read sets the molecule belongs to)"""
    import pyarrow as pa, pyarrow.compute as pc
    sch = READS_SCHEMA()
    if P is None: return pa.table({f.name: pa.array([], type=f.type) for f in sch})
    t = P.filter(pc.equal(P.column("MGS"), mgs))
    rs = t.column("read_sets").combine_chunks()
    rs = pc.split_pattern(pc.if_else(pc.equal(rs, ""), pa.scalar(None, pa.string()), rs), ",")
    rs = pc.fill_null(rs, pa.scalar([], pa.list_(pa.string())))
    cols = {f.name: t.column(f.name) for f in sch if f.name in t.column_names}
    cols["sample"] = pa.array([sample] * t.num_rows, pa.string()); cols["read_sets"] = rs
    return pa.table({f.name: cols[f.name] for f in sch}).cast(sch).sort_by([("molecule", "ascending"), ("piece", "ascending")])


def mgs_files(d, mgs, sample):
    return dict(sites=f"{d}/{mgs}/{sample}.parquet", junctions=f"{d}/{mgs}/junctions/{sample}.parquet",
                edges=f"{d}/{mgs}/edges/{sample}.parquet", reads=f"{d}/{mgs}/reads/{sample}.parquet")


def write_mgs(d, mgs, sample, sites, junc, edges, reads):
    import mmtpq
    f = mgs_files(d, mgs, sample)
    for k in ("junctions", "edges", "reads"): os.makedirs(os.path.dirname(f[k]), exist_ok=True)
    mmtpq.write(f["reads"], reads, READS_SCHEMA()); mmtpq.write(f["edges"], edges, EDGES_SCHEMA())
    mmtpq.write(f["junctions"], junc, JUNCTION_SCHEMA())
    return mmtpq.write(f["sites"], sites, "ctg_sites")          # last: the sites file marks the MGS as written


def cmd_write(a):
    """sorted kinetics export + junction outputs + derep-ctg -> the four files of every listed MGS (empty files for an
    MGS without rows: every glob over a folder keeps one schema)"""
    import pyarrow as pa, pyarrow.parquet as pq, mmtpq
    S = json.loads(a.settings)
    todo = [l.rstrip("\n").split("\t") for l in open(a.mgs_list) if l.strip()]       # MGS, kin_fp, lift_fp, MAGs, derep .fp
    P = pq.read_table(f"{a.jx}/pieces.parquet") if os.path.isfile(f"{a.jx}/pieces.parquet") else None
    sets = collections.defaultdict(list)
    for r in read_sets_tsv(f"{a.jx}/sets.tsv"): sets[r["MGS"]].append(r)
    ctg = {}
    for l in open(a.contigs):                                  # contig, MAG, MGS, length of the sample's MAG contigs
        c, m, g, n = l.rstrip("\n").split("\t")
        ctg[c] = (m, g, int(n))
    base = [f for f in mmtpq.schema("ctg_sites") if f.name not in LIFT_COLS + ["kin_fp", "lift_fp"]]
    have = set()
    if a.export and os.path.isfile(a.export):
        have = set(pq.read_table(a.export, columns=["MGS"]).column("MGS").unique().to_pylist())
    if have - {m for m, *_ in todo}: die(f"rows of MGS not in the list: {sorted(have - {m for m, *_ in todo})[:5]}")
    for mgs, kin_fp, lift_fp, mags, dfp in todo:
        L = Lift(os.path.join(a.derep, mgs), mags.split(","), dfp)
        t = pq.read_table(a.export, filters=[("MGS", "=", mgs)]) if mgs in have else pa.table({f.name: pa.array([], type=f.type) for f in base})
        st = sites_table(t, L, kin_fp, lift_fp)
        jt = junction_table(sets.get(mgs, []), L, st, kin_fp, lift_fp, a.sample)
        rd = reads_table(P, mgs, a.sample)
        ce = contig_end_rows(rd, {c: v[2] for c, v in ctg.items() if v[1] == mgs}, {c: v[0] for c, v in ctg.items()}, mgs, a.sample, lift_fp, S)
        ed = pa.concat_tables([edges_table(rd, jt, L, mgs, a.sample, lift_fp, S), ce])
        n = write_mgs(a.out, mgs, a.sample, st, jt, ed, rd)
        rs_ = jt.column("read_set").to_pylist()
        print(f"{mgs}\t{n}\t{len(set(rs_))}\t{len({r_ for r_, c in zip(rs_, jt.column('called').to_pylist()) if c})}", flush=True)


def cmd_relift(a):
    """existing files of one MGS x sample re-annotated with the current derep-ctg result (no kinetics)"""
    import pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc
    S = json.loads(a.settings)
    f = mgs_files(a.out, a.mgs, a.sample)
    old = pq.read_table(f["sites"])
    mags = sorted(set(pq.read_table(f["reads"], columns=["MAG"]).column("MAG").to_pylist()) |
                  set(old.column("MAG").to_pylist()) | set(a.mags.split(",") if a.mags else []))
    L = Lift(os.path.join(a.derep, a.mgs), mags, a.derep_fp)
    kin = old.column("kin_fp")[0].as_py() if old.num_rows else a.kin_fp
    st = sites_table(old, L, kin or a.kin_fp, a.lift_fp)
    oj = pq.read_table(f["junctions"])
    sets = [{k: ("" if v is None else ("1" if v is True else "0" if v is False else (",".join(v) if isinstance(v, list) else str(v))))
             for k, v in r.items()} for r in oj.to_pylist()]
    jt = junction_table(sets, L, st, kin or a.kin_fp, a.lift_fp, a.sample)
    rd = pq.read_table(f["reads"])
    oe = pq.read_table(f["edges"])
    ce = oe.filter(pc.equal(oe.column("kind"), "contig_end"))
    ce = ce.set_column(ce.schema.get_field_index("lift_fp"), "lift_fp", pa.array([a.lift_fp] * ce.num_rows, pa.string()))
    ed = pa.concat_tables([edges_table(rd, jt, L, a.mgs, a.sample, a.lift_fp, S), ce.cast(EDGES_SCHEMA())])
    print(write_mgs(a.out, a.mgs, a.sample, st, jt, ed, rd))


def cmd_derep_mags(a):
    """MGS TAB MAG of the derep-ctg results of the listed MGS (an unreadable result lists nothing: not current)"""
    import pyarrow.parquet as pq
    for m in (l.strip() for l in open(a.mgs)):
        if not m: continue
        try: mags = pq.read_table(f"{a.derep}/{m}/mags.parquet", columns=["MAG"]).column("MAG").to_pylist()
        except Exception as e: log(f"WARNING: {m}: derep-ctg result not readable ({e})"); continue
        sys.stdout.write("".join(f"{m}\t{g}\n" for g in mags))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("refs"); p.add_argument("--fasta", required=True); p.add_argument("--contigs", required=True); p.add_argument("--out", required=True)
    p = sp.add_parser("junctions"); p.add_argument("--mags", required=True); p.add_argument("--settings", required=True); p.add_argument("--out", required=True)
    p = sp.add_parser("write"); p.add_argument("--export"); p.add_argument("--jx", required=True); p.add_argument("--derep", required=True)
    p.add_argument("--out", required=True); p.add_argument("--sample", required=True); p.add_argument("--mgs-list", required=True)
    p.add_argument("--settings", required=True); p.add_argument("--contigs", required=True)
    p = sp.add_parser("relift"); p.add_argument("--out", required=True); p.add_argument("--sample", required=True)
    p.add_argument("--mgs", required=True); p.add_argument("--derep", required=True); p.add_argument("--kin-fp", default="")
    p.add_argument("--lift-fp", required=True); p.add_argument("--settings", required=True); p.add_argument("--mags", default="")
    p.add_argument("--derep-fp", default="")
    p = sp.add_parser("derep_mags"); p.add_argument("--derep", required=True); p.add_argument("--mgs", required=True)
    a = ap.parse_args()
    {"refs": cmd_refs, "junctions": cmd_junctions, "write": cmd_write, "relift": cmd_relift, "derep_mags": cmd_derep_mags}[a.cmd](a)


if __name__ == "__main__":
    main()
