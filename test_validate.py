"""Reference-build and header compatibility tests."""
import os
import unittest

from gtt.errors import BuildMismatchError
from gtt.validate import check_compatibility, infer_build, require_compatible, validate_header
from gtt.vcfio import read_header

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


def H(name):
    return read_header(os.path.join(FIX, name))


class TestBuildInference(unittest.TestCase):
    def test_grch38_from_contig_length(self):
        build, evidence = infer_build(H("synthetic_baseline.vcf"))
        self.assertEqual(build, "GRCh38")
        self.assertIn("248956422", evidence)

    def test_grch37_from_contig_length(self):
        build, _ = infer_build(H("synthetic_followup_grch37.vcf"))
        self.assertEqual(build, "GRCh37")

    def test_unprefixed_names_still_fingerprint(self):
        build, _ = infer_build(H("synthetic_followup_nochr.vcf"))
        self.assertEqual(build, "GRCh38")


class TestCompatibility(unittest.TestCase):
    def test_matching_inputs_are_compatible(self):
        rep = check_compatibility(H("synthetic_baseline.vcf"), H("synthetic_followup.vcf"))
        self.assertTrue(rep.compatible)
        self.assertEqual(rep.baseline_build, "GRCh38")
        self.assertEqual(len(rep.shared_contigs), 3)
        self.assertIn("SYNTHETIC_SUBJECT", rep.sample_note)

    def test_build_mismatch_is_refused(self):
        rep = check_compatibility(H("synthetic_baseline.vcf"), H("synthetic_followup_grch37.vcf"))
        self.assertFalse(rep.compatible)
        self.assertIn("build mismatch", rep.refusal_reason)
        self.assertIn("does not lift over", rep.refusal_reason)
        with self.assertRaises(BuildMismatchError):
            require_compatible(rep)

    def test_build_mismatch_override(self):
        rep = check_compatibility(H("synthetic_baseline.vcf"), H("synthetic_followup_grch37.vcf"),
                                  allow_build_mismatch=True)
        self.assertTrue(rep.compatible)
        self.assertTrue(any("OVERRIDE" in w for w in rep.warnings))

    def test_chr_naming_strict_refuses(self):
        rep = check_compatibility(H("synthetic_baseline.vcf"), H("synthetic_followup_nochr.vcf"))
        self.assertTrue(rep.chr_naming_mismatch)
        self.assertFalse(rep.compatible)
        self.assertIn("--chr-naming harmonize", rep.refusal_reason)

    def test_chr_naming_harmonize_allows(self):
        rep = check_compatibility(H("synthetic_baseline.vcf"), H("synthetic_followup_nochr.vcf"),
                                  chr_naming="harmonize")
        self.assertTrue(rep.compatible)
        self.assertEqual(len(rep.shared_contigs), 3)

    def test_synthetic_fixtures_are_flagged(self):
        warnings = validate_header(H("synthetic_baseline.vcf"))
        self.assertTrue(any("SYNTHETIC TEST FIXTURE" in w for w in warnings))


if __name__ == "__main__":
    unittest.main()
