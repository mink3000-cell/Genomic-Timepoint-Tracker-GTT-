"""Optional cross-check of GTT's left-alignment against `bcftools norm`.

Skipped cleanly when bcftools is not installed. bcftools is a validation tool
only: GTT itself never calls it at runtime and never depends on it.
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "scripts"))

from crosscheck_bcftools import crosscheck, have_bcftools   # noqa: E402

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
REF = os.path.join(FIX, "synthetic_reference.fa")
BASE = os.path.join(FIX, "synthetic_repeat_baseline.vcf")
FOLLOW = os.path.join(FIX, "synthetic_repeat_followup.vcf")


@unittest.skipUnless(have_bcftools(), "bcftools is not installed (optional validation only)")
class TestAgainstBcftools(unittest.TestCase):
    def test_followup_fixture_matches_bcftools_norm(self):
        # The follow-up fixture has no REF mismatch, so every allele is comparable.
        res = crosscheck(FOLLOW, REF)
        self.assertTrue(res["agree"],
                        f"GTT only: {res['only_gtt']}; bcftools only: {res['only_bcftools']}")

    def test_baseline_fixture_differs_only_at_the_ref_mismatch(self):
        # chrT:5 C>T disagrees with the reference. GTT keeps it as written and
        # counts it; bcftools refuses to normalize it too, but may drop it.
        res = crosscheck(BASE, REF)
        unexpected = [k for k in res["only_gtt"] if k[1] != 5]
        self.assertEqual(unexpected, [], f"unexpected GTT-only keys: {unexpected}")


if __name__ == "__main__":
    unittest.main()
