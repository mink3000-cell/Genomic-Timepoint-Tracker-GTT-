"""End-to-end comparison with --reference and --regions, plus v0.1.0 parity.

The fixture pair ``synthetic_repeat_baseline.vcf`` / ``synthetic_repeat_followup.vcf``
writes the same three indels at different anchors, so the counts below change in
a fully predictable way when a reference is supplied. Every expected number was
computed by hand from the fixture sequence before this test was run; see
``tests/fixtures/README.md``.
"""
import io
import json
import os
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

from gtt.cli import main
from gtt.compare import compare_files
from gtt.report import check_wording, render_markdown, write_report

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
BASE = os.path.join(FIX, "synthetic_repeat_baseline.vcf")
FOLLOW = os.path.join(FIX, "synthetic_repeat_followup.vcf")
REF = os.path.join(FIX, "synthetic_reference.fa")
BED = os.path.join(FIX, "synthetic_regions.bed")             # chrT 1-30
BED_NARROW = os.path.join(FIX, "synthetic_regions_narrow.bed")   # chrT 1-15, chrS 1-10
BED_EDGE = os.path.join(FIX, "synthetic_regions_edge.bed")       # chrT 1-20

# Fixtures used for the v0.1.0 parity check.
OLD_BASE = os.path.join(FIX, "synthetic_baseline.vcf")
OLD_FOLLOW = os.path.join(FIX, "synthetic_followup.vcf")

#: Exactly the keys v0.1.0 wrote. New keys must appear only when a new option is used.
V010_SUMMARY_KEYS = {
    "gtt_version", "baseline", "followup", "compatibility", "options",
    "baseline_stats", "followup_stats", "counts", "baseline_only_contexts",
    "followup_only_contexts", "contigs_compared", "contigs_baseline_only",
    "contigs_followup_only", "output_files",
}
V010_COUNT_KEYS = {
    "reported_in_both", "baseline_only", "followup_only",
    "genotype_differences_among_shared", "reported_alleles_baseline",
    "reported_alleles_followup", "jaccard", "normalization_level",
}
V010_STATS_KEYS = {
    "records", "alleles_emitted", "multiallelic_records", "skipped_star",
    "skipped_no_alt", "skipped_symbolic", "skipped_not_pass", "symbolic_emitted",
    "not_pass_records", "no_call_records", "hom_ref_records", "n_contigs_seen",
}


def variants(rows):
    return sorted(r["variant"] for r in rows)


class TestTrimOnlyBaseline(unittest.TestCase):
    """Without --reference the two anchorings stay apart - the v0.1.0 limitation."""

    @classmethod
    def setUpClass(cls):
        cls.res = compare_files(BASE, FOLLOW, collect=True)

    def test_counts(self):
        self.assertEqual(self.res.normalization_level, "trim_only")
        self.assertEqual(self.res.n_shared, 1)
        self.assertEqual(self.res.n_baseline_only, 5)
        self.assertEqual(self.res.n_followup_only, 4)

    def test_the_same_indel_appears_in_both_difference_tables(self):
        self.assertIn("chrT-17-AA-A", variants(self.res.baseline_only_rows))
        self.assertIn("chrT-11-AA-A", variants(self.res.followup_only_rows))


class TestLeftAligned(unittest.TestCase):
    """With --reference the equivalent spellings collapse onto one key."""

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.TemporaryDirectory()
        cls.res = compare_files(BASE, FOLLOW, outdir=cls.dir.name, reference=REF,
                                collect=True, emit_shared=True)
        cls.md = render_markdown(cls.res)

    @classmethod
    def tearDownClass(cls):
        cls.dir.cleanup()

    def test_counts(self):
        self.assertEqual(self.res.normalization_level, "left_aligned")
        self.assertEqual(self.res.n_shared, 4)
        self.assertEqual(self.res.n_baseline_only, 2)
        self.assertEqual(self.res.n_followup_only, 1)

    def test_shared_keys_are_the_left_aligned_ones(self):
        self.assertEqual(variants(self.res.shared_rows),
                         ["chrS-1-AA-A", "chrT-10-CA-C", "chrT-20-TCA-T", "chrT-4-T-G"])

    def test_only_the_real_differences_remain(self):
        # chrT-5-C-T stays baseline-only because its REF does not match the reference.
        self.assertEqual(variants(self.res.baseline_only_rows), ["chrT-35-G-A", "chrT-5-C-T"])
        self.assertEqual(variants(self.res.followup_only_rows), ["chrT-37-C-T"])

    def test_per_allele_outcomes_are_counted(self):
        b = self.res.baseline_stats
        self.assertEqual(b.left_aligned_shifted, 2)
        self.assertEqual(b.left_aligned_unchanged, 2)
        self.assertEqual(b.ref_mismatches, 1)
        self.assertEqual(b.truncated_at_contig_start, 1)
        self.assertEqual(self.res.followup_stats.ref_mismatches, 0)
        self.assertEqual(self.res.followup_stats.truncated_at_contig_start, 1)

    def test_ref_mismatches_are_written_out(self):
        path = os.path.join(self.dir.name, "ref_mismatches.tsv")
        self.assertTrue(os.path.exists(path))
        with open(path, encoding="utf-8") as fh:
            lines = fh.read().splitlines()
        self.assertEqual(len(lines), 2)                       # header + one mismatch
        self.assertIn("chrT", lines[1])
        self.assertEqual(self.res.ref_mismatch_rows[0]["ref_in_file"], "C")
        self.assertEqual(self.res.ref_mismatch_rows[0]["ref_in_fasta"], "A")

    def test_report_states_the_normalization_level_and_the_mismatch(self):
        self.assertIn("left_aligned", self.md)
        self.assertIn("Reference check and left-alignment", self.md)
        self.assertIn("Tan,", self.md)
        self.assertIn("1 REF allele disagrees with the reference sequence", self.md)
        body = self.md.split("## 5. How to read these numbers")[0]
        self.assertEqual(check_wording(body), [])

    def test_summary_json_carries_the_reference_block(self):
        paths = write_report(self.res, self.dir.name)
        with open(paths["summary_json"], encoding="utf-8") as fh:
            summary = json.load(fh)
        self.assertEqual(summary["counts"]["normalization_level"], "left_aligned")
        self.assertEqual(summary["counts"]["ref_mismatches"], 1)
        self.assertEqual(len(summary["ref_mismatches"]), 1)


class TestRegions(unittest.TestCase):
    def test_variants_outside_the_regions_are_not_interrogated(self):
        res = compare_files(BASE, FOLLOW, reference=REF, regions=BED, collect=True)
        self.assertEqual(res.n_shared, 3)
        self.assertEqual(res.n_baseline_only, 1)          # the REF-mismatch allele
        self.assertEqual(res.n_followup_only, 0)
        self.assertEqual(res.n_not_interrogated_baseline, 2)
        self.assertEqual(res.n_not_interrogated_followup, 2)

    def test_not_interrogated_is_not_counted_as_absent(self):
        res = compare_files(BASE, FOLLOW, reference=REF, regions=BED, collect=True)
        outside = variants(res.not_interrogated_rows)
        self.assertIn("chrT-35-G-A", outside)             # would be baseline-only without --regions
        self.assertNotIn("chrT-35-G-A", variants(res.baseline_only_rows))
        for row in res.not_interrogated_rows:
            self.assertIn("not a statement that the variant is present or absent",
                          row["statement"])

    def test_two_bed_files_are_intersected(self):
        res = compare_files(BASE, FOLLOW, reference=REF, regions_baseline=BED,
                            regions_followup=BED_NARROW, collect=True)
        self.assertEqual(res.region_provenance["comparison_regions"]["total_bases"], 15)
        self.assertEqual(res.n_shared, 2)                 # chrT-20-TCA-T now falls outside
        self.assertEqual(res.n_not_interrogated_baseline, 3)
        self.assertEqual(res.n_not_interrogated_followup, 3)

    def test_region_mode_span_is_stricter_at_an_edge(self):
        start = compare_files(BASE, FOLLOW, reference=REF, regions=BED_EDGE)
        span = compare_files(BASE, FOLLOW, reference=REF, regions=BED_EDGE, region_mode="span")
        self.assertEqual(start.n_shared, 3)               # POS 20 is inside
        self.assertEqual(span.n_shared, 2)                # REF 'TCA' runs to 22, outside
        self.assertEqual(span.n_not_interrogated_baseline, 3)

    def test_table_is_written_only_when_asked(self):
        with tempfile.TemporaryDirectory() as td:
            compare_files(BASE, FOLLOW, outdir=td, regions=BED)
            self.assertFalse(os.path.exists(os.path.join(td, "not_interrogated.tsv")))
        with tempfile.TemporaryDirectory() as td:
            compare_files(BASE, FOLLOW, outdir=td, regions=BED, emit_not_interrogated=True)
            path = os.path.join(td, "not_interrogated.tsv")
            self.assertTrue(os.path.exists(path))
            with open(path, encoding="utf-8") as fh:
                self.assertEqual(len(fh.read().splitlines()), 5)   # header + 4 alleles

    def test_report_section_present(self):
        res = compare_files(BASE, FOLLOW, regions=BED)
        md = render_markdown(res)
        self.assertIn("Region restriction (what was interrogated)", md)
        self.assertIn("not interrogated", md)
        self.assertEqual(check_wording(md.split("## 5. How to read these numbers")[0]), [])


class TestV010Parity(unittest.TestCase):
    """Without the new options, v0.2.0 must produce v0.1.0's output."""

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.TemporaryDirectory()
        cls.res = compare_files(OLD_BASE, OLD_FOLLOW, outdir=cls.dir.name, emit_shared=True)
        cls.paths = write_report(cls.res, cls.dir.name)
        with open(cls.paths["summary_json"], encoding="utf-8") as fh:
            cls.summary = json.load(fh)

    @classmethod
    def tearDownClass(cls):
        cls.dir.cleanup()

    def test_counts_unchanged(self):
        self.assertEqual(self.res.n_shared, 7)
        self.assertEqual(self.res.n_baseline_only, 4)
        self.assertEqual(self.res.n_followup_only, 3)
        self.assertEqual(self.res.n_genotype_differences, 2)
        self.assertEqual(self.res.normalization_level, "trim_only")

    def test_no_new_output_files(self):
        self.assertEqual(sorted(os.listdir(self.dir.name)),
                         ["baseline_only.tsv", "comparison_report.md", "comparison_summary.json",
                          "followup_only.tsv", "genotype_differences.tsv", "shared.tsv"])

    def test_summary_keys_are_exactly_the_v010_ones(self):
        self.assertEqual(set(self.summary), V010_SUMMARY_KEYS)
        self.assertEqual(set(self.summary["counts"]), V010_COUNT_KEYS)
        self.assertEqual(set(self.summary["baseline_stats"]), V010_STATS_KEYS)
        self.assertEqual(set(self.summary["followup_stats"]), V010_STATS_KEYS)

    def test_report_keeps_the_trim_only_caveat(self):
        with open(self.paths["report_md"], encoding="utf-8") as fh:
            md = fh.read()
        self.assertIn("Full left-alignment requires", md)
        self.assertNotIn("Region restriction", md)
        self.assertNotIn("Reference check", md)


class TestCliPhase1(unittest.TestCase):
    def test_compare_with_reference_and_regions(self):
        with tempfile.TemporaryDirectory() as td:
            buf = io.StringIO()
            with redirect_stdout(buf), redirect_stderr(io.StringIO()):
                rc = main(["compare", BASE, FOLLOW, "-o", td, "--reference", REF,
                           "--regions", BED, "--emit-not-interrogated"])
            self.assertEqual(rc, 0)
            out = json.loads(buf.getvalue())
            self.assertEqual(out["normalization_level"], "left_aligned")
            self.assertEqual(out["reported_in_both"], 3)
            self.assertEqual(out["not_interrogated_baseline"], 2)
            self.assertTrue(os.path.exists(os.path.join(td, "not_interrogated.tsv")))
            self.assertTrue(os.path.exists(os.path.join(td, "ref_mismatches.tsv")))

    def test_faidx_command(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = main(["faidx", REF, "--fetch", "chrT:11-18"])
        self.assertEqual(rc, 0)
        self.assertIn("AAAAAAAA", buf.getvalue())

    def test_regions_command_reports_the_intersection(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = main(["regions", BED, BED_NARROW])
        self.assertEqual(rc, 0)
        info = json.loads(buf.getvalue())
        self.assertEqual(info["intersection"]["total_bases"], 15)

    def test_missing_reference_is_a_clean_error(self):
        with tempfile.TemporaryDirectory() as td:
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                rc = main(["compare", BASE, FOLLOW, "-o", td,
                           "--reference", os.path.join(td, "nope.fa")])
            self.assertEqual(rc, 2)


if __name__ == "__main__":
    unittest.main()
