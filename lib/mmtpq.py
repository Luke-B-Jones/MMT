#!/usr/bin/env python3
"""MMT result tables: fixed schemas and one parquet writer for every result file.

  split    a sorted kinetics export (MGScall) -> one file per MGS: OUT/<MGS>/<sample>.parquet
  migrate  MGScall files of the old layout (MGS=<m>/sample=<S>/data.parquet) -> the new one (a list, one process)

Every result file has a fixed schema (also when empty: a glob read takes its schema from the first file), MGS and
sample as columns, zstd, dictionary encoding for everything except positions, which are delta-encoded (DELTA_BINARY_
PACKED: ~3 B/row less than a dictionary on positions), row groups of 256k rows. Files are written to a temp name and
renamed. DuckDB >= 1.1 and pyarrow read them. Positions are signed (int32): differences of positions never overflow.
"""
import argparse, os, sys
import pyarrow as pa, pyarrow.parquet as pq

S, F, I8, U8, U16, I16, I32, B = pa.string(), pa.float32(), pa.int8(), pa.uint8(), pa.uint16(), pa.int16(), pa.int32(), pa.bool_()
KIN = [("coverage", I32), ("score", I32), ("ipdRatio", F), ("tMean", F), ("tErr", F), ("modelPrediction", F), ("frac", F),
       ("fracLow", F), ("fracUp", F), ("call", S), ("identificationQv", F), ("pile_depth", I32), ("n_ref", I32),
       ("top_alt", S), ("n_alt", I32), ("ctx_alt_pct", U8), ("ctx_alt_n", I32), ("clip_ctx_n", U16)]
SCHEMA = {
    "mgscall": [("MGS", S), ("sample", S), ("contig", S), ("pos", I32), ("strand", S), ("base", S)] + KIN,
    "ctg_sites": [("MGS", S), ("sample", S), ("read_set", S), ("side", I8), ("dist", I32), ("MAG", S), ("contig", S),
                  ("pos", I32), ("strand", S), ("base", S), ("block", S), ("instance", I32), ("block_pos", I32),
                  ("block_strand", S), ("ctx_dist", U8), ("gene", S), ("up_gene", S), ("up_dist", I16)] + KIN
                 + [("kin_fp", S), ("lift_fp", S)],
}
DELTA = {"pos", "block_pos", "r_start", "r_end", "c_start", "c_end", "b_start", "q_start", "q_end", "w_start", "w_end", "start", "end"}
ROWS = 262144


def schema(name):
    return pa.schema(SCHEMA[name]) if isinstance(name, str) else name


def cast(t, sch):
    """table -> exactly this schema (column order, types) ; missing columns are an error"""
    miss = [f.name for f in sch if f.name not in t.column_names]
    if miss: raise SystemExit(f"ERROR: table lacks columns {miss}")
    return pa.table([t.column(f.name).cast(f.type) for f in sch], schema=sch)


class Writer:
    """streaming writer of one result file: temp name, schema, encodings ; close() renames into place"""
    def __init__(s, path, sch):
        s.path, s.sch, s.n = path, schema(sch), 0
        s.tmp = os.path.join(os.path.dirname(path) or ".", f".{os.path.basename(path)}.{os.uname()[1]}.{os.getpid()}.tmp")
        d = [f.name for f in s.sch if f.name in DELTA and pa.types.is_integer(f.type)]
        kw = dict(compression="zstd", write_statistics=True)
        try:
            s.w = pq.ParquetWriter(s.tmp, s.sch, use_dictionary=[f.name for f in s.sch if f.name not in d],
                                   column_encoding={c: "DELTA_BINARY_PACKED" for c in d}, **kw)
        except (TypeError, ValueError, NotImplementedError, pa.ArrowNotImplementedError):    # old pyarrow: no column_encoding
            s.w = pq.ParquetWriter(s.tmp, s.sch, **kw)

    def write(s, t):
        t = cast(t, s.sch)
        if t.num_rows:
            s.w.write_table(t, row_group_size=ROWS); s.n += t.num_rows

    def close(s):
        s.w.close(); os.replace(s.tmp, s.path)
        return s.n

    def abort(s):
        try: s.w.close()
        finally:
            if os.path.exists(s.tmp): os.remove(s.tmp)


def write(path, t, sch):
    w = Writer(path, sch)
    try: w.write(t)
    except BaseException: w.abort(); raise
    return w.close()


# ============================================================================ MGScall
def cmd_split(a):
    """sorted export (MGS, contig, pos, strand) -> D/<MGS>/<sample>.parquet ; prints 'MGS<TAB>rows' per file"""
    import numpy as np, pyarrow.compute as pc
    f = pq.ParquetFile(a.input)
    cur, w, done = None, None, set()
    try:
        for b in f.iter_batches(batch_size=ROWS):
            t = pa.Table.from_batches([b])
            if t.column("MGS").null_count: raise SystemExit("ERROR: export has rows without MGS")
            d = pc.dictionary_encode(t.column("MGS")).combine_chunks()
            codes, names = d.indices.to_numpy(), d.dictionary.to_pylist()
            cut = np.flatnonzero(np.diff(codes)) + 1
            for i, j in zip(np.r_[0, cut], np.r_[cut, len(codes)]):
                m = names[codes[i]]
                if m != cur:
                    if w: print(f"{cur}\t{w.close()}", flush=True); done.add(cur); w = None
                    if m in done: raise SystemExit(f"ERROR: export not grouped by MGS ({m} again)")
                    cur = m
                    os.makedirs(os.path.join(a.out, cur), exist_ok=True)
                    w = Writer(os.path.join(a.out, cur, f"{a.sample}.parquet"), "mgscall")
                w.write(t.slice(int(i), int(j - i)))
        if w: print(f"{cur}\t{w.close()}", flush=True); w = None
    except BaseException:
        if w: w.abort()
        raise


def cmd_migrate(a):
    """old-layout MGScall files (no MGS / sample columns, maybe unsorted) -> new files, sorted, with both columns ;
    --list: old TAB new TAB MGS TAB sample per line ; prints 'old<TAB>rows' per file written (the caller removes it)"""
    for l in open(a.list):
        if not l.strip(): continue
        old, new, mgs, sample = l.rstrip("\n").split("\t")
        try: t = pq.read_table(old)
        except Exception as e:                     # a damaged old file: reported, the caller sets it aside
            print(f"{old}\tERROR {type(e).__name__}: {str(e)[:200]}".replace("\n", " "), flush=True); continue
        n = t.num_rows
        for k, v in (("MGS", mgs), ("sample", sample)):
            if k not in t.column_names: t = t.append_column(k, pa.array([v] * n, S))
        t = t.sort_by([("contig", "ascending"), ("pos", "ascending"), ("strand", "ascending")])
        os.makedirs(os.path.dirname(new) or ".", exist_ok=True)
        if write(new, t, "mgscall") != n: raise SystemExit(f"ERROR: {old}: row count changed")
        print(f"{old}\t{n}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("split"); p.add_argument("--input", required=True); p.add_argument("--out", required=True)
    p.add_argument("--sample", required=True)
    p = sp.add_parser("migrate"); p.add_argument("--list", required=True)
    a = ap.parse_args()
    {"split": cmd_split, "migrate": cmd_migrate}[a.cmd](a)


if __name__ == "__main__":
    main()
