"""BED parsing, merging, intersection and the BED/VCF coordinate conversion."""
import gzip
import os
import tempfile
import unittest

from gtt.errors import BedFormatError
from gtt.regions import RegionSet, build_comparison_regions, intersect_intervals, merge_intervals

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
BED = os.path.join(FIX, "synthetic_regions.bed")             # chrT 0-30
BED_NARROW = os.path.join(FIX, "synthetic_regions_narrow.bed")   # chrT 0-15, chrS 0-10
BED_NOCHR = os.path.join(FIX, "synthetic_regions_nochr.bed")     # T 0-30
BED_SPLIT = os.path.join(FIX, "synthetic_regions_split.bed")     # 0-10, 10-20, 18-30


class TestIntervalArithmetic(unittest.TestCase):
    def test_merge_joins_touching_and_overlapping(self):
        self.assertEqual(merge_intervals([(0, 10), (10, 20), (18, 30)]), [(0, 30)])
        self.assertEqual(merge_intervals([(50, 60), (0, 10)]), [(0, 10), (50, 60)])
        self.assertEqual(merge_intervals([(5, 5), (0, 3)]), [(0, 3)])   # empty span dropped

    def test_intersection(self):
        self.assertEqual(intersect_intervals([(0, 30)], [(10, 15), (20, 40)]),
                         [(10, 15), (20, 30)])
        self.assertEqual(intersect_intervals([(0, 10)], [(10, 20)]), [])   # touching, not overlapping


class TestParsing(unittest.TestCase):
    def test_comment_track_and_browser_lines_are_skipped(self):
        rs = RegionSet.from_bed(BED)
        self.assertEqual(rs.n_intervals, 1)
        self.assertEqual(rs.total_bases, 30)

    def test_intervals_are_merged_on_load(self):
        rs = RegionSet.from_bed(BED_SPLIT)
        self.assertEqual(rs.n_intervals, 1)
        self.assertEqual(rs.total_bases, 30)

    def test_gzip_bed(self):
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "r.bed.gz")
            with gzip.open(path, "wt") as fh:
                fh.write("chrT\t0\t10\n")
            self.assertEqual(RegionSet.from_bed(path).total_bases, 10)

    def test_whitespace_separated_bed_is_accepted(self):
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "r.bed")
            with open(path, "w", newline="\n") as fh:
                fh.write("chrT 0 10\n")
            self.assertEqual(RegionSet.from_bed(path).total_bases, 10)

    def test_malformed_beds_are_refused(self):
        with tempfile.TemporaryDirectory() as td:
            def write(text):
                path = os.path.join(td, "bad.bed")
                with open(path, "w", newline="\n") as fh:
                    fh.write(text)
                return path
            for text in ("chrT\t10\n", "chrT\tx\t10\n", "chrT\t20\t10\n", "chrT\t-5\t10\n"):
                with self.subTest(text=text):
                    with self.assertRaises(BedFormatError):
                        RegionSet.from_bed(write(text))
            with self.assertRaises(BedFormatError):
                RegionSet.from_bed(write("# only a comment\n"))
            with self.assertRaises(BedFormatError):
                RegionSet.from_bed(os.path.join(td, "absent.bed"))


class TestMembership(unittest.TestCase):
    """BED is 0-based half-open, VCF is 1-based inclusive: chrT 0-30 means POS 1..30."""

    def setUp(self):
        self.rs = RegionSet.from_bed(BED)

    def test_boundaries(self):
        self.assertTrue(self.rs.contains("chrT", 1))
        self.assertTrue(self.rs.contains("chrT", 30))
        self.assertFalse(self.rs.contains("chrT", 31))
        self.assertFalse(self.rs.contains("chrT", 0))

    def test_unknown_contig_is_outside(self):
        self.assertFalse(self.rs.contains("chrS", 5))
        self.assertFalse(self.rs.has_contig("chrS"))

    def test_chr_prefix_is_ignored_in_both_directions(self):
        self.assertTrue(self.rs.contains("T", 5))
        nochr = RegionSet.from_bed(BED_NOCHR)
        self.assertTrue(nochr.contains("chrT", 5))
        self.assertEqual(nochr.total_bases, self.rs.total_bases)

    def test_span_mode_requires_the_whole_ref_allele_inside(self):
        self.assertTrue(self.rs.contains_span("chrT", 28, 3))   # 28,29,30 all inside
        self.assertFalse(self.rs.contains_span("chrT", 29, 3))  # 29,30,31 - spills over
        # The same position passes the default start test, which looks at POS only.
        self.assertTrue(self.rs.contains("chrT", 29))


class TestCombination(unittest.TestCase):
    def test_single_file(self):
        rs, prov = build_comparison_regions(regions=BED)
        self.assertEqual(rs.total_bases, 30)
        self.assertIn("regions", prov)
        self.assertEqual(prov["combination"], "synthetic_regions.bed")

    def test_two_files_are_intersected(self):
        rs, prov = build_comparison_regions(regions_baseline=BED, regions_followup=BED_NARROW)
        self.assertEqual(rs.total_bases, 15)        # chrT 0-30 AND chrT 0-15
        self.assertEqual(rs.contigs, ["T"])         # chrS is in only one file
        self.assertFalse(rs.contains("chrS", 5))
        self.assertIn("intersection of", prov["combination"])

    def test_no_files_means_no_restriction(self):
        rs, prov = build_comparison_regions()
        self.assertIsNone(rs)
        self.assertEqual(prov, {})


if __name__ == "__main__":
    unittest.main()
