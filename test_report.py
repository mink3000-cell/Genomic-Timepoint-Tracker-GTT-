"""Report rendering and vocabulary enforcement."""
import json
import os
import tempfile
import unittest

from gtt.compare import compare_files
from gtt.report import check_wording, render_markdown, write_report

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
BASE = os.path.join(FIX, "synthetic_baseline.vcf")
FOLLOW = os.path.join(FIX, "synthetic_followup.vcf")


class TestReport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = compare_files(BASE, FOLLOW, collect=True)
        cls.md = render_markdown(cls.res)

    def test_required_conservative_phrases(self):
        self.assertIn("reported in the baseline file but not reported in the follow-up file", self.md)
        self.assertIn("reported in the follow-up file but not reported in the baseline file", self.md)
        self.assertIn("it does not and cannot establish that a variant was gained or lost", self.md)

    def test_no_forbidden_interpretive_terms(self):
        # The caveat text names the concepts only in the permitted negated forms.
        body = self.md.split("## 5. How to read these numbers")[0]
        self.assertEqual(check_wording(body), [])

    def test_synthetic_banner_present(self):
        self.assertIn("Synthetic test data", self.md)
        self.assertIn("not a biological observation", self.md.lower())

    def test_build_and_normalization_are_explicit(self):
        self.assertIn("GRCh38", self.md)
        self.assertIn("trim_only", self.md)
        self.assertIn("Full left-alignment requires", self.md)

    def test_counts_rendered(self):
        self.assertIn("| reported in both files | 7 |", self.md)

    def test_write_report_emits_md_and_json(self):
        with tempfile.TemporaryDirectory() as td:
            paths = write_report(self.res, td)
            self.assertTrue(os.path.exists(paths["report_md"]))
            with open(paths["summary_json"], encoding="utf-8") as fh:
                summary = json.load(fh)
            self.assertEqual(summary["counts"]["reported_in_both"], 7)
            self.assertEqual(summary["counts"]["baseline_only"], 4)
            self.assertEqual(summary["compatibility"]["baseline_build"], "GRCh38")
            self.assertTrue(summary["baseline"]["synthetic"])


if __name__ == "__main__":
    unittest.main()
