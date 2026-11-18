"""Sample-specific MF4 depth -> MAG coverage per Gbp contributing to assembly depth.

No MGS score, BAM byte size, read-selection count or co-assembly pooled depth is used.
Missing/unvalidated evidence is unavailable, never an inferred zero.
"""
import csv
import functools
import gzip
import hashlib
import json
import math
import os

METHOD = 1
UNIT = "fold_coverage_per_Gbp_assembly_mapped"
METRICS = ["normalized_abundance", "mag_mean_depth", "mag_breadth_1x", "mag_depth_sd",
           "mag_mapped_bases", "sample_mapped_bases", "mag_reference_bp", "sample_reference_bp",
           "mapped_base_fraction"]


def identity(path):
    st = os.stat(path)
    return [os.path.realpath(path), st.st_size, st.st_mtime_ns]


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def coverage_overrides(path, samples):
    """TSV sample/coverage; paths relative to the TSV. Optional format: bedgraph/percontig."""
    if not path:
        return {}
    out = {}
    with open(path) as fh:
        reader = csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")), delimiter="\t")
        if not {"sample", "coverage"} <= set(reader.fieldnames or []):
            raise ValueError(f"{path}: expected sample and coverage columns")
        for row in reader:
            s = row["sample"].strip()
            if s not in samples or s in out:
                raise ValueError(f"{path}: unknown or repeated sample {s!r}")
            p = row["coverage"].strip()
            if not p:
                raise ValueError(f"{path}: empty coverage path for {s}")
            p = os.path.abspath(os.path.join(os.path.dirname(path), p))
            fmt = row.get("format", "").strip()
            if fmt not in ("", "bedgraph", "percontig"):
                raise ValueError(f"{path}: {s}: format must be bedgraph or percontig")
            if not os.path.isfile(p):
                raise ValueError(f"{path}: {s}: coverage file does not exist: {p}")
            out[s] = dict(path=p, format=fmt, discovery="--coverage-map")
    return out


def find_coverage(run, sample, overrides):
    """Bounded discovery in this sample only. Never borrow assembly-level or support depth."""
    if sample in overrides:
        return overrides[sample]
    sd = os.path.join(run, sample)
    for ext in ("bam", "cram"):
        base = os.path.join(sd, "mapping", f"{sample}-smd.{ext}.coverage")
        for p in (base + ".gz", base):
            if os.path.isfile(p):
                return dict(path=p, format="bedgraph", discovery="sample mapping")
    # rdCover's summaries are sample-specific even for a shared assembly.
    base = os.path.join(sd, "assemblies", "metag", "ContigStats", "Coverage.percontig")
    for p in (base + ".gz", base):
        if os.path.isfile(p):
            return dict(path=p, format="percontig", discovery="sample ContigStats")
    return dict(path="", format="", discovery="missing")


@functools.lru_cache(maxsize=4)
def _lengths(path, signature):
    lengths = {}
    name, length = None, 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                if name in lengths:
                    raise ValueError(f"{path}: duplicate FASTA contig {name}")
                length = 0
            else:
                if name is None and line.strip():
                    raise ValueError(f"{path}: sequence before first FASTA header")
                length += len(line.strip())
    if name is not None:
        lengths[name] = length
    if not lengths or any(L <= 0 for L in lengths.values()):
        raise ValueError(f"{path}: empty assembly or zero-length contig")
    return lengths


def read_depth(path, fmt, lengths):
    """Contig -> (depth integral, covered bp or None, squared-depth integral or None).

    Full reference coverage is required, including zero-depth intervals: MF4 uses
    samtools depth -aa. This detects truncated/stale tracks, not just matching MAG names.
    """
    totals, ends = {}, {}
    previous, closed = None, set()
    with open_text(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            f = line.rstrip("\r\n").split("\t")
            if not fmt:
                fmt = "bedgraph" if len(f) == 4 else "percontig" if len(f) == 2 else ""
            if not fmt or len(f) != (4 if fmt == "bedgraph" else 2):
                raise ValueError(f"{path}:{lineno}: expected four-column bedGraph or two-column contig mean depth")
            c = f[0]
            if c not in lengths:
                raise ValueError(f"{path}:{lineno}: contig absent from sample assembly: {c}")
            if fmt == "percontig":
                depth = float(f[1])
                if c in totals or not math.isfinite(depth) or depth < 0:
                    raise ValueError(f"{path}:{lineno}: duplicate contig or invalid depth")
                totals[c] = (depth * lengths[c], None, None)
                continue
            start, end, depth = int(f[1]), int(f[2]), float(f[3])
            if (start != ends.get(c, 0) or end <= start or end > lengths[c]
                    or not math.isfinite(depth) or depth < 0 or not depth.is_integer()):
                raise ValueError(f"{path}:{lineno}: invalid depth, gap/overlap, or reference length mismatch")
            if c != previous:
                if c in closed:
                    raise ValueError(f"{path}:{lineno}: contig intervals are not contiguous")
                if previous is not None:
                    closed.add(previous)
                previous = c
            n, d = end - start, int(depth)
            bp, covered, squares = totals.get(c, (0, 0, 0))
            totals[c] = (bp + n * d, covered + (n if d else 0), squares + n * d * d)
            ends[c] = end
    missing = set(lengths) - set(totals)
    if missing or (fmt == "bedgraph" and any(ends[c] != L for c, L in lengths.items() if c in ends)):
        raise ValueError(f"{path}: coverage does not span the complete assembly including zero-depth contigs "
                         f"({len(missing)} contigs missing); sparse tracks are not assumed complete")
    return totals, fmt


def quantify(job):
    """One sample; cached independently of molecule assignments and their MAPQ settings."""
    source, fasta, mags = job["source"], job["fasta"], job["mags"]
    result = dict(sample=job["sample"], source=source, status="unavailable", reason="", mags={})
    if not source["path"]:
        result["reason"] = "no sample-specific MF4 depth track or contig-mean table; supply --coverage-map"
        return result
    if not fasta:
        result["reason"] = "sample assembly FASTA unavailable"
        return result
    try:
        inputs = [METHOD, identity(source["path"]), source, identity(fasta), mags]
        fp = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()
        cache = job["cache"]
        if os.path.isfile(cache):
            with open(cache) as fh:
                saved = json.load(fh)
            if saved.get("fingerprint") == fp:
                return saved["result"]
        lengths = _lengths(fasta, tuple(identity(fasta)))
        depth, fmt = read_depth(source["path"], source["format"], lengths)
        total = math.fsum(v[0] for v in depth.values())
        result.update(status="ok" if total > 0 else "zero_sample_depth", sample_mapped_bases=total,
                      sample_reference_bp=sum(lengths.values()), source=dict(source, format=fmt))
        owners = {}
        for m, cs in mags.items():
            for c in set(cs):
                owners.setdefault(c, []).append(m)
        for m, cs in mags.items():
            cs = set(cs)
            r = {k: None for k in METRICS}
            r.update(status="ok", reason="", sample_mapped_bases=total,
                     sample_reference_bp=result["sample_reference_bp"])
            if not cs or cs - set(lengths):
                r.update(status="reference_mismatch", reason="MAG contigs absent from sample assembly")
            elif any(len(owners[c]) > 1 for c in cs):
                r.update(status="shared_contigs", reason="contig assigned to more than one final MAG")
            else:
                L = sum(lengths[c] for c in cs)
                bp = math.fsum(depth[c][0] for c in cs)
                mean = bp / L
                r.update(mag_reference_bp=L, mag_mapped_bases=bp, mag_mean_depth=mean,
                         normalized_abundance=mean * 1e9 / total if total else None,
                         mapped_base_fraction=bp / total if total else None)
                if fmt == "bedgraph":
                    r.update(mag_breadth_1x=sum(depth[c][1] for c in cs) / L,
                             mag_depth_sd=math.sqrt(max(0, math.fsum(depth[c][2] for c in cs) / L - mean * mean)))
                if not total:
                    r.update(status="zero_sample_depth", reason="sample-wide depth integral is zero")
            result["mags"][m] = r
        if identity(source["path"]) != inputs[1] or identity(fasta) != inputs[3]:
            raise ValueError("coverage or assembly changed during quantification; rerun after MF4 finishes")
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        temp = f"{cache}.{os.getpid()}.tmp"
        with open(temp, "w") as fh:
            json.dump(dict(fingerprint=fp, result=result), fh)
        os.replace(temp, cache)
    except (OSError, ValueError, EOFError, OverflowError) as e:
        result.update(status="invalid_coverage", reason=str(e))
        result["mags"] = {}
    return result


def pair_metrics(sample, mag, evidence):
    """Foreign assembly pairs are unprofiled, not zero; absent MGS has no effect on MAG depth."""
    empty = {k: None for k in METRICS}
    if not sample["asm"] or not mag["asm"]:
        return dict(empty, abundance_status="reference_unavailable", abundance_reason="sample or MAG assembly unavailable")
    if os.path.realpath(sample["asm"]) != os.path.realpath(mag["asm"]):
        return dict(empty, abundance_status="not_profiled", abundance_reason="sample was not mapped to this MAG's assembly")
    if mag["status"] != "ok":
        return dict(empty, abundance_status="MAG_unavailable", abundance_reason=mag["status"])
    if not evidence:
        return dict(empty, abundance_status="unavailable", abundance_reason="no coverage evidence")
    r = evidence.get("mags", {}).get(mag["id"])
    if not r:
        return dict(empty, abundance_status=evidence["status"], abundance_reason=evidence["reason"])
    return dict({k: r[k] for k in METRICS}, abundance_status=r["status"], abundance_reason=r["reason"])
