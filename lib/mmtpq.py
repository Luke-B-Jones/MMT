#!/usr/bin/env python3
"""MMT result tables: fixed schemas and one parquet writer for every result file.

  split    a sorted kinetics export (MGScall) -> one file per MGS: OUT/<MGS>/<sample>.parquet, or the path --paths
           gives for that MGS (MGS TAB relative path)

Every result file has a fixed schema (also when empty: a glob read takes its schema from the first file), MGS, MAG
and sample as columns (constant columns cost < 0.1 bit per row), zstd, dictionary encoding for everything except
positions, which are delta-encoded (DELTA_BINARY_PACKED: ~3 B/row less than a dictionary on positions), row groups of
256k rows, and the column conventions as file metadata (key "mmt": pyarrow schema.metadata, duckdb
parquet_kv_metadata()). Files are written to a temp name and renamed. DuckDB >= 1.1 and pyarrow read them. Positions are
signed (int32): differences of positions never overflow.
"""
import argparse, json, os, sys
import pyarrow as pa, pyarrow.parquet as pq

S, F, I8, U8, U16, I16, I32, B = pa.string(), pa.float32(), pa.int8(), pa.uint8(), pa.uint16(), pa.int16(), pa.int32(), pa.bool_()
KIN = [("coverage", I32), ("score", I32), ("ipdRatio", F), ("tMean", F), ("tErr", F), ("modelPrediction", F), ("frac", F),
       ("fracLow", F), ("fracUp", F), ("call", S), ("identificationQv", F), ("pile_depth", I32), ("n_ref", I32),
       ("top_alt", S), ("n_alt", I32), ("ctx_alt_pct", U8), ("ctx_alt_n", I32), ("clip_ctx_n", U16)]
SCHEMA = {
    "mgscall": [("MGS", S), ("MAG", S), ("sample", S), ("contig", S), ("pos", I32), ("strand", S), ("base", S)] + KIN,
    "ctg_sites": [("MGS", S), ("sample", S), ("read_set", S), ("side", I8), ("dist", I32), ("MAG", S), ("contig", S),
                  ("pos", I32), ("strand", S), ("base", S), ("block", S), ("instance", I32), ("block_pos", I32),
                  ("block_strand", S), ("ctx_dist", U8), ("gene", S), ("up_gene", S), ("up_dist", I16)] + KIN
                 + [("kin_fp", S), ("lift_fp", S)],
}
META = {   # column conventions, stored in every file (key "mmt")
    "mgscall": {
        "MAG": "MAG(s) whose reads made this table (mgs-rep: every MAG of the MGS in the sample, comma-separated)",
        "contig": "reference contig, its MATAFILER4 name (mgs-rep / mag-rep: of the representative MAG ; mag-mag: of the MAG itself)",
        "pos": "1-based position on the contig",
        "strand": "+ / - : strand of the base whose kinetics are reported",
        "base": "base on that strand (- strand: complement of the contig base)",
        "rows": "every A and C with >= 3 kinetic reads (ipdSummary writes no row below 3) ; G/T only where ipdSummary "
                "flags a kinetic signal (call modified_base) ; no row = < 3 reads",
        "coverage": "kinetic reads (molecules) on this strand at this position",
        "frac": "methylated fraction, only on rows called m6A / m4C with coverage > 10 (ipdSummary)",
        "pile_depth/n_ref/top_alt/n_alt": "reads at the position on this strand / carrying the contig base / most common other "
                                          "observation / its reads (same reads as the kinetics)",
        "ctx_alt_pct/ctx_alt_n/clip_ctx_n": "within +-ctx bp on this strand: highest % (and count) of reads differing from "
                                            "the contig ; reads with a clip point",
    },
}
DELTA = {"pos", "block_pos", "r_start", "r_end", "c_start", "c_end", "b_start", "q_start", "q_end", "w_start", "w_end", "start", "end"}
ROWS = 262144


def schema(name):
    if not isinstance(name, str):
        return name
    return pa.schema(SCHEMA[name], metadata={"mmt": json.dumps(META[name])} if name in META else None)


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
    """sorted export (MGS, contig, pos, strand) -> D/<MGS>/<sample>.parquet (or --paths) ; prints 'MGS<TAB>rows' per file"""
    import numpy as np, pyarrow.compute as pc
    paths = dict(l.rstrip("\n").split("\t")[:2] for l in open(a.paths) if l.strip()) if a.paths else {}
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
                    path = os.path.join(a.out, paths.get(cur, os.path.join(cur, f"{a.sample}.parquet")))
                    os.makedirs(os.path.dirname(path), exist_ok=True)
                    w = Writer(path, "mgscall")
                w.write(t.slice(int(i), int(j - i)))
        if w: print(f"{cur}\t{w.close()}", flush=True); w = None
    except BaseException:
        if w: w.abort()
        raise


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("split"); p.add_argument("--input", required=True); p.add_argument("--out", required=True)
    p.add_argument("--sample", required=True); p.add_argument("--paths", default="")
    a = ap.parse_args()
    {"split": cmd_split}[a.cmd](a)


if __name__ == "__main__":
    main()
