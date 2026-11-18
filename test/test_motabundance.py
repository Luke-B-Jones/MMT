"""Analytical abundance truth, discovery boundaries, incomplete evidence and cache invalidation."""
import gzip
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import motabundance as a
import motexport


class AbundanceTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.fa = self.root / "assembly.fa"
        self.fa.write_text(">a\n" + "A" * 100 + "\n>b\n" + "C" * 300 + "\n>unbinned\n" + "G" * 600 + "\n")
        self.track = self.root / "depth.gz"
        self.write_track("a\t0\t100\t10\nb\t0\t300\t2\nunbinned\t0\t600\t5\n")

    def tearDown(self):
        self.tmp.cleanup()

    def write_track(self, text):
        with gzip.open(self.track, "wt") as fh:
            fh.write(text)

    def quantify(self, path=None, fmt="bedgraph", mags=None):
        return a.quantify(dict(sample="S", source=dict(path=str(path or self.track), format=fmt),
                               fasta=str(self.fa), mags=mags or {"MAG": ["a", "b"]},
                               cache=str(self.root / "cache" / "S.json")))

    def test_length_weighted_truth_and_unbinned_denominator(self):
        r = self.quantify(); m = r["mags"]["MAG"]
        self.assertEqual(r["status"], "ok")
        self.assertEqual(r["sample_mapped_bases"], 4600)
        self.assertEqual(m["mag_reference_bp"], 400)
        self.assertEqual(m["mag_mapped_bases"], 1600)
        self.assertEqual(m["mag_mean_depth"], 4)
        self.assertAlmostEqual(m["normalized_abundance"], 4e9 / 4600)
        self.assertAlmostEqual(m["mapped_base_fraction"], 1600 / 4600)
        self.assertEqual(m["mag_breadth_1x"], 1)
        self.assertAlmostEqual(m["mag_depth_sd"], 12 ** .5)

    def test_depth_rescaling_is_invariant_and_cache_refreshes(self):
        first = self.quantify()["mags"]["MAG"]["normalized_abundance"]
        self.write_track("a\t0\t100\t20\nb\t0\t300\t4\nunbinned\t0\t600\t10\n")
        r = self.quantify()
        self.assertEqual(r["sample_mapped_bases"], 9200)
        self.assertEqual(r["mags"]["MAG"]["normalized_abundance"], first)

    def test_zero_MAG_is_distinct_from_zero_sample(self):
        self.write_track("a\t0\t100\t0\nb\t0\t300\t0\nunbinned\t0\t600\t5\n")
        m = self.quantify()["mags"]["MAG"]
        self.assertEqual(m["normalized_abundance"], 0)
        self.assertEqual(m["mag_breadth_1x"], 0)
        self.write_track("a\t0\t100\t0\nb\t0\t300\t0\nunbinned\t0\t600\t0\n")
        r = self.quantify()
        self.assertEqual(r["status"], "zero_sample_depth")
        self.assertIsNone(r["mags"]["MAG"]["normalized_abundance"])

    def test_bad_or_incomplete_tracks_never_yield_a_number(self):
        for text in ("a\t0\t100\t10\nb\t0\t300\t2\n",  # missing unbinned denominator
                     "a\t0\t100\t10\nb\t0\t301\t2\nunbinned\t0\t600\t5\n",  # stale reference
                     "a\t1\t100\t10\nb\t0\t300\t2\nunbinned\t0\t600\t5\n",  # sparse
                     "a\t0\t100\tnan\nb\t0\t300\t2\nunbinned\t0\t600\t5\n"):
            with self.subTest(text=text):
                self.write_track(text)
                r = self.quantify()
                self.assertEqual(r["status"], "invalid_coverage")
                self.assertFalse(r["mags"])

    def test_contig_mean_fallback_has_no_invented_breadth(self):
        p = self.root / "means"
        p.write_text("a\t10\nb\t2\nunbinned\t5\n")
        m = self.quantify(p, "percontig")["mags"]["MAG"]
        self.assertAlmostEqual(m["normalized_abundance"], 4e9 / 4600)
        self.assertIsNone(m["mag_breadth_1x"])
        self.assertIsNone(m["mag_depth_sd"])

    def test_shared_contigs_are_not_double_counted(self):
        r = self.quantify(mags={"M1": ["a"], "M2": ["a", "b"]})
        self.assertTrue(all(m["status"] == "shared_contigs" for m in r["mags"].values()))
        self.assertTrue(all(m["normalized_abundance"] is None for m in r["mags"].values()))

    def test_sample_specific_discovery_and_explicit_override(self):
        d = self.root / "run" / "S" / "mapping"; d.mkdir(parents=True)
        (d / "OTHER-smd.bam.coverage.gz").touch()
        (d / "S.sup-smd.bam.coverage.gz").touch()
        self.assertEqual(a.find_coverage(str(self.root / "run"), "S", {})["path"], "")
        p = d / "S-smd.cram.coverage.gz"; p.touch()
        self.assertEqual(a.find_coverage(str(self.root / "run"), "S", {})["path"], str(p))
        override = self.root / "overrides.tsv"
        override.write_text("sample\tcoverage\nS\tdepth.gz\n")
        got = a.coverage_overrides(str(override), {"S"})
        self.assertEqual(got["S"]["path"], str(self.track))
        self.assertEqual(a.find_coverage(str(self.root / "run"), "S", got)["discovery"], "--coverage-map")
        override.write_text("sample\tcoverage\nUNKNOWN\tdepth.gz\n")
        with self.assertRaises(ValueError):
            a.coverage_overrides(str(override), {"S"})

    def test_foreign_assembly_pair_is_not_zero(self):
        r = a.pair_metrics({"asm": "/one"}, {"asm": "/two", "status": "ok"}, None)
        self.assertEqual(r["abundance_status"], "not_profiled")
        self.assertIsNone(r["normalized_abundance"])

    def test_multiple_binners_require_explicit_choice_and_plain_table_supported(self):
        for b in ("SB", "MB2"):
            d = self.root / f"Bin_{b}"; d.mkdir()
            (d / "MAGvsGC.txt").write_text("MAG\tMGS\n")
        with self.assertRaises(SystemExit):
            motexport.find_binner(str(self.root), "")
        self.assertEqual(motexport.find_binner(str(self.root), "MB2"), "MB2")

    def test_compressed_bin_membership(self):
        import mf4
        p = self.root / "Binning" / "SB"
        p.mkdir(parents=True)
        with gzip.open(p / "S.gz", "wt") as fh:
            fh.write("Sequence ID\tBin\na\tbin1\nb\tbin1\n")
        got, paths = mf4.bin_contigs({"S__bin1": "MGS.1"}, {}, {"S": str(self.root)}, "SB")
        self.assertEqual(got["S__bin1"]["contigs"], ["a", "b"])
        self.assertTrue(paths[0].endswith(".gz"))

    def test_coassembled_samples_use_independent_depth(self):
        other = self.root / "other.depth"
        other.write_text("a\t0\t100\t10\nb\t0\t300\t2\nunbinned\t0\t600\t50\n")
        first = self.quantify()["mags"]["MAG"]["normalized_abundance"]
        second = self.quantify(other)["mags"]["MAG"]["normalized_abundance"]
        self.assertAlmostEqual(second, 4e9 / 31600)
        self.assertNotEqual(first, second)


if __name__ == "__main__":
    unittest.main()
