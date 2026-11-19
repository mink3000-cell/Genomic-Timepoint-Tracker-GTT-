"""Normalization unit tests."""
import unittest

from gtt.normalize import classify_variant, is_symbolic, normalize_allele


class TestNormalizeAllele(unittest.TestCase):
    def test_snv_unchanged(self):
        self.assertEqual(normalize_allele(100, "A", "G"), (100, "A", "G"))

    def test_lowercase_is_uppercased(self):
        self.assertEqual(normalize_allele(100, "a", "g"), (100, "A", "G"))

    def test_common_suffix_trimmed(self):
        # CTTA>CA and CTT>C describe the same 2 bp deletion.
        self.assertEqual(normalize_allele(100, "CTTA", "CA"), (100, "CTT", "C"))
        self.assertEqual(normalize_allele(100, "CTT", "C"), (100, "CTT", "C"))

    def test_suffix_is_trimmed_before_prefix(self):
        # Right-then-left trimming (the vt/bcftools convention): CTT>CT loses the
        # shared trailing T, giving CT>C at the original POS - not TT>T at POS+1.
        self.assertEqual(normalize_allele(100, "CTT", "CT"), (100, "CT", "C"))

    def test_common_prefix_shifts_pos(self):
        # No shared suffix here, so the shared prefix is trimmed and POS advances.
        self.assertEqual(normalize_allele(100, "CAT", "CG"), (101, "AT", "G"))

    def test_padded_snv_trims_to_snv(self):
        self.assertEqual(normalize_allele(100, "AGT", "CGT"), (100, "A", "C"))
        self.assertEqual(normalize_allele(100, "GTA", "GTC"), (102, "A", "C"))

    def test_insertion_already_minimal(self):
        self.assertEqual(normalize_allele(602113, "T", "TGCCCA"), (602113, "T", "TGCCCA"))

    def test_at_least_one_base_retained(self):
        pos, ref, alt = normalize_allele(100, "AAAA", "AAAA")
        self.assertTrue(len(ref) >= 1 and len(alt) >= 1)

    def test_symbolic_and_placeholder_alleles_untouched(self):
        self.assertEqual(normalize_allele(100, "N", "<DEL>"), (100, "N", "<DEL>"))
        self.assertEqual(normalize_allele(100, "CTT", "*"), (100, "CTT", "*"))
        self.assertEqual(normalize_allele(100, "CTT", "."), (100, "CTT", "."))
        self.assertTrue(is_symbolic("<DUP:TANDEM>"))
        self.assertTrue(is_symbolic("C[chr2:321682["))
        self.assertFalse(is_symbolic("ACGT"))

    def test_idempotent(self):
        once = normalize_allele(100, "CTTA", "CA")
        twice = normalize_allele(*once)
        self.assertEqual(once, twice)


class TestClassify(unittest.TestCase):
    def test_classes(self):
        self.assertEqual(classify_variant("A", "G"), "snv")
        self.assertEqual(classify_variant("AC", "GT"), "mnv")
        self.assertEqual(classify_variant("T", "TGCCCA"), "insertion")
        self.assertEqual(classify_variant("CTT", "C"), "deletion")
        self.assertEqual(classify_variant("CTT", "GA"), "indel")
        self.assertEqual(classify_variant("N", "<DEL>"), "symbolic")
        self.assertEqual(classify_variant("CTT", "*"), "upstream_deletion")


if __name__ == "__main__":
    unittest.main()
