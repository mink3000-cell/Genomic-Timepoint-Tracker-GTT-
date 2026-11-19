"""Comparison-engine tests against fixtures with hand-computed expected results."""
import os
import tempfile
import unittest

from gtt.compare import (
    CTX_FILTERED,
    CTX_NO_CALL,
    CTX_NO_RECORD,
    CTX_OTHER_ALLELE,
    compare_files,
)
from gtt.errors import BuildMismatchError, UnsortedInputError

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
BASE = os.path.join(FIX, "synthetic_baseline.vcf")
FOLLOW = os.path.join(FIX, "synthetic_followup.vcf")

EXPECTED_SHARED = {
    "chr1-100000-A-G", "chr1-200000-C-T", "chr1-300000-CTT-C", "chr1-800000-A-G",
    "chr1-900000-T-C", "chr2-500000-G-A", "chr2-800000-G-A",
}
EXPECTED_BASELINE_ONLY = {
    "chr1-400000-T-A", "chr1-450000-G-C", "chr2-500000-G-T", "chrX-1000000-A-T",
}
EXPECTED_FOLLOWUP_ONLY = {"chr1-450000-G-T", "chr1-600000-C-G", "chr2-700000-A-AT"}


class TestDefaultComparison(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = compare_files(BASE, FOLLOW, collect=True, emit_shared=True)

    def _members(self, rows):
        return {r["variant"] for r in rows}

    def test_counts(self):
        self.assertEqual(self.res.n_shared, 7)
        self.assertEqual(self.res.n_baseline_only, 4)
        self.assertEqual(self.res.n_followup_only, 3)

    def test_membership(self):
        self.assertEqual(self._members(self.res.shared_rows), EXPECTED_SHARED)
        self.assertEqual(self._members(self.res.baseline_only_rows), EXPECTED_BASELINE_ONLY)
        self.assertEqual(self._members(self.res.followup_only_rows), EXPECTED_FOLLOWUP_ONLY)

    def test_sets_are_disjoint_and_complete(self):
        s, b, f = (self._members(self.res.shared_rows),
                   self._members(self.res.baseline_only_rows),
                   self._members(self.res.followup_only_rows))
        self.assertEqual(s & b, set())
        self.assertEqual(s & f, set())
        self.assertEqual(self.res.n_reported_baseline, len(s | b))
        self.assertEqual(self.res.n_reported_followup, len(s | f))

    def test_normalization_matches_equivalent_deletion(self):
        row = next(r for r in self.res.shared_rows if r["variant"] == "chr1-300000-CTT-C")
        self.assertEqual(row["baseline_as_written"], "chr1:300000:CTT>CTT".replace("CTT>CTT", "CTT>C"))
        self.assertEqual(row["followup_as_written"], "chr1:300000:CTTA>CA")

    def test_genotype_differences_do_not_change_membership(self):
        self.assertEqual(self.res.n_genotype_differences, 2)
        diffs = {r["variant"] for r in self.res.genotype_difference_rows}
        self.assertEqual(diffs, {"chr1-200000-C-T", "chr1-900000-T-C"})
        self.assertTrue(diffs <= EXPECTED_SHARED)

    def test_multiallelic_genotype_projected_per_allele(self):
        row = next(r for r in self.res.shared_rows if r["variant"] == "chr2-500000-G-A")
        self.assertEqual(row["baseline_genotype"], "1/2")
        self.assertEqual(row["followup_genotype"], "0/1")
        self.assertEqual(row["baseline_allele_copies"], 1)
        self.assertEqual(row["followup_allele_copies"], 1)
        self.assertTrue(row["genotype_agrees"])

    def test_contexts(self):
        ctx = {r["variant"]: r["other_file_context"] for r in self.res.baseline_only_rows}
        self.assertEqual(ctx["chr1-400000-T-A"], CTX_NO_RECORD)
        self.assertEqual(ctx["chr1-450000-G-C"], CTX_OTHER_ALLELE)
        self.assertEqual(ctx["chr2-500000-G-T"], CTX_OTHER_ALLELE)
        self.assertEqual(ctx["chrX-1000000-A-T"], CTX_NO_RECORD)
        fctx = {r["variant"]: r["other_file_context"] for r in self.res.followup_only_rows}
        self.assertEqual(fctx["chr1-450000-G-T"], CTX_OTHER_ALLELE)
        self.assertEqual(fctx["chr1-600000-C-G"], CTX_NO_RECORD)

    def test_statements_are_conservative(self):
        for r in self.res.baseline_only_rows + self.res.followup_only_rows:
            self.assertIn("not reported in the", r["statement"])
            self.assertIn("not evidence of a biological gain or loss", r["statement"])
            for banned in ("lost", "newly acquired", "de novo", "disappeared"):
                self.assertNotIn(banned, r["statement"].lower())

    def test_contig_only_in_one_file_is_tracked(self):
        self.assertIn("chrX", self.res.contigs_baseline_only)
        self.assertEqual(self.res.contigs_followup_only, [])
        self.assertEqual(sorted(self.res.contigs_compared), ["chr1", "chr2"])

    def test_summary_jaccard(self):
        self.assertAlmostEqual(self.res.summary()["jaccard"], 7 / 14)


class TestStrictnessOptions(unittest.TestCase):
    def test_pass_only_moves_filtered_allele_to_baseline_only(self):
        res = compare_files(BASE, FOLLOW, collect=True, emit_shared=True, pass_only=True)
        self.assertEqual(res.n_shared, 6)
        self.assertEqual(res.n_baseline_only, 5)
        self.assertEqual(res.n_followup_only, 3)
        row = next(r for r in res.baseline_only_rows if r["variant"] == "chr1-800000-A-G")
        self.assertEqual(row["other_file_context"], CTX_FILTERED)
        self.assertEqual(row["context_detail"], "LowQual")
        self.assertIn("filtered", row["statement"])

    def test_require_carried_moves_no_call_to_baseline_only(self):
        res = compare_files(BASE, FOLLOW, collect=True, emit_shared=True, require_carried=True)
        self.assertEqual(res.n_shared, 6)
        self.assertEqual(res.n_baseline_only, 5)
        row = next(r for r in res.baseline_only_rows if r["variant"] == "chr1-900000-T-C")
        self.assertEqual(row["other_file_context"], CTX_NO_CALL)
        self.assertEqual(res.n_genotype_differences, 1)

    def test_keep_star_adds_an_allele_to_both(self):
        res = compare_files(BASE, FOLLOW, collect=True, emit_shared=True, skip_star=False)
        self.assertEqual(res.n_shared, 8)
        self.assertIn("chr2-800000-G-*", {r["variant"] for r in res.shared_rows})


class TestIdentityAndOrdering(unittest.TestCase):
    def test_self_comparison_yields_no_differences(self):
        res = compare_files(BASE, BASE, collect=True, emit_shared=True)
        self.assertEqual((res.n_baseline_only, res.n_followup_only), (0, 0))
        self.assertEqual(res.n_shared, 11)
        self.assertEqual(res.n_genotype_differences, 0)
        self.assertEqual(res.summary()["jaccard"], 1.0)

    def test_reversing_inputs_swaps_the_difference_sets(self):
        fwd = compare_files(BASE, FOLLOW, collect=True)
        rev = compare_files(FOLLOW, BASE, collect=True)
        self.assertEqual(fwd.n_shared, rev.n_shared)
        self.assertEqual({r["variant"] for r in fwd.baseline_only_rows},
                         {r["variant"] for r in rev.followup_only_rows})
        self.assertEqual({r["variant"] for r in fwd.followup_only_rows},
                         {r["variant"] for r in rev.baseline_only_rows})

    def test_contig_order_difference_is_tolerated(self):
        # Follow-up lists chr2 before chr1; results must be unchanged.
        import shutil
        with tempfile.TemporaryDirectory() as td:
            lines = open(FOLLOW, encoding="utf-8").read().splitlines()
            head = [l for l in lines if l.startswith("#")]
            body = [l for l in lines if not l.startswith("#")]
            reordered = [l for l in body if l.startswith("chr2")] + \
                        [l for l in body if not l.startswith("chr2")]
            path = os.path.join(td, "reordered.vcf")
            with open(path, "w", newline="\n", encoding="utf-8") as fh:
                fh.write("\n".join(head + reordered) + "\n")
            res = compare_files(BASE, path, collect=True)
            self.assertEqual((res.n_shared, res.n_baseline_only, res.n_followup_only), (7, 4, 3))

    def test_interleaved_contigs_rejected(self):
        with self.assertRaises(UnsortedInputError):
            compare_files(BASE, os.path.join(FIX, "synthetic_interleaved.vcf"))

    def test_build_mismatch_refused_before_comparing(self):
        with self.assertRaises(BuildMismatchError):
            compare_files(BASE, os.path.join(FIX, "synthetic_followup_grch37.vcf"))

    def test_chr_naming_harmonize(self):
        res = compare_files(BASE, os.path.join(FIX, "synthetic_followup_nochr.vcf"),
                            chr_naming="harmonize", collect=True)
        self.assertEqual({r["variant"] for r in res.shared_rows}, {"1-100000-A-G"})
        self.assertEqual({r["variant"] for r in res.followup_only_rows}, {"1-600000-C-G"})


class TestOutputFiles(unittest.TestCase):
    def test_tsv_files_written(self):
        with tempfile.TemporaryDirectory() as td:
            res = compare_files(BASE, FOLLOW, outdir=td, emit_shared=True)
            for name in ("baseline_only", "followup_only", "shared", "genotype_differences"):
                path = os.path.join(td, f"{name}.tsv")
                self.assertTrue(os.path.exists(path), name)
                with open(path, encoding="utf-8") as fh:
                    rows = fh.read().strip().split("\n")
                self.assertGreaterEqual(len(rows), 2)
            with open(os.path.join(td, "baseline_only.tsv"), encoding="utf-8") as fh:
                self.assertEqual(len(fh.read().strip().split("\n")), 1 + res.n_baseline_only)


if __name__ == "__main__":
    unittest.main()
