"""End-to-end CLI tests."""
import io
import json
import os
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr

from gtt.cli import main

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
BASE = os.path.join(FIX, "synthetic_baseline.vcf")
FOLLOW = os.path.join(FIX, "synthetic_followup.vcf")


def run(argv):
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        code = main(argv)
    return code, out.getvalue(), err.getvalue()


class TestCli(unittest.TestCase):
    def test_inspect_with_count(self):
        code, out, _ = run(["inspect", BASE, "--count"])
        self.assertEqual(code, 0)
        info = json.loads(out)
        self.assertEqual(info["inferred_build"], "GRCh38")
        self.assertEqual(info["stats"]["records"], 10)
        self.assertEqual(info["stats"]["alleles_emitted"], 11)
        self.assertTrue(info["synthetic"])

    def test_compare_writes_outputs(self):
        with tempfile.TemporaryDirectory() as td:
            code, out, err = run(["compare", BASE, FOLLOW, "-o", td])
            self.assertEqual(code, 0)
            counts = json.loads(out)
            self.assertEqual(counts["reported_in_both"], 7)
            self.assertEqual(counts["baseline_only"], 4)
            self.assertEqual(counts["followup_only"], 3)
            self.assertIn("comparison_report.md", err)
            for f in ("comparison_report.md", "comparison_summary.json",
                      "baseline_only.tsv", "followup_only.tsv"):
                self.assertTrue(os.path.exists(os.path.join(td, f)), f)
            self.assertFalse(os.path.exists(os.path.join(td, "shared.tsv")))

    def test_build_mismatch_exits_nonzero(self):
        with tempfile.TemporaryDirectory() as td:
            code, _, err = run(["compare", BASE,
                                os.path.join(FIX, "synthetic_followup_grch37.vcf"), "-o", td])
            self.assertEqual(code, 2)
            self.assertIn("BuildMismatchError", err)


if __name__ == "__main__":
    unittest.main()
