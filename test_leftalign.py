"""Left-alignment against a reference, with hand-computed expected answers.

The fixture contig ``chrT`` is::

    pos  1         11        21        31
         ACGTACGTAC AAAAAAAA GT CACACACACACA TTGGCCAA
                    ^ 8-base A homopolymer, anchored by the C at position 10
                               ^ six CA units, anchored by the T at position 20

Every expectation below was worked out on paper from that sequence before the
code was run. ``chrS`` is ``AAAAACCCCC``: a deletion inside its leading A run
cannot be anchored at all, because there is no base before position 1.
"""
import os
import unittest

from gtt.fasta import IndexedFasta
from gtt.normalize import (
    LA_AT_CONTIG_START,
    LA_CONTIG_ABSENT,
    LA_NOT_APPLICABLE,
    LA_REF_MISMATCH,
    LA_SHIFT_LIMIT,
    LA_SHIFTED,
    LA_UNCHANGED,
    left_align_allele,
    normalize_allele,
)

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
REF = os.path.join(FIX, "synthetic_reference.fa")


class TestLeftAlign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fa = IndexedFasta(REF)

    @classmethod
    def tearDownClass(cls):
        cls.fa.close()

    def align(self, contig, pos, ref, alt):
        """Trim, then left-align - the order GTT uses on every allele."""
        p, r, a = normalize_allele(pos, ref, alt)
        return left_align_allele(self.fa, contig, p, r, a)

    # -- homopolymer ------------------------------------------------------
    def test_homopolymer_deletion_shifts_to_the_anchor_base(self):
        # Deleting one A from the run at 11-18, written at three anchors.
        for pos in (11, 13, 17):
            with self.subTest(pos=pos):
                p, r, a, status, seen = self.align("chrT", pos, "AA", "A")
                self.assertEqual((p, r, a), (10, "CA", "C"))
                self.assertEqual(status, LA_SHIFTED)
                self.assertIsNone(seen)

    def test_homopolymer_insertion_shifts_to_the_anchor_base(self):
        for pos in (11, 15, 18):
            with self.subTest(pos=pos):
                p, r, a, status, _ = self.align("chrT", pos, "A", "AA")
                self.assertEqual((p, r, a), (10, "C", "CA"))
                self.assertEqual(status, LA_SHIFTED)

    def test_two_spellings_of_one_deletion_become_one_key(self):
        left = self.align("chrT", 11, "AA", "A")[:3]
        right = self.align("chrT", 17, "AA", "A")[:3]
        self.assertEqual(left, right)

    # -- dinucleotide repeat ---------------------------------------------
    def test_dinucleotide_repeat_deletion(self):
        # One CA unit removed from the repeat at 21-32, at two anchors.
        for pos in (24, 30):
            with self.subTest(pos=pos):
                p, r, a, status, _ = self.align("chrT", pos, "ACA", "A")
                self.assertEqual((p, r, a), (20, "TCA", "T"))
                self.assertEqual(status, LA_SHIFTED)

    def test_dinucleotide_repeat_insertion(self):
        p, r, a, status, _ = self.align("chrT", 32, "A", "ACA")
        self.assertEqual((p, r, a), (20, "T", "TCA"))
        self.assertEqual(status, LA_SHIFTED)

    # -- things that must not move ---------------------------------------
    def test_snv_is_untouched(self):
        p, r, a, status, _ = self.align("chrT", 4, "T", "G")
        self.assertEqual((p, r, a, status), (4, "T", "G", LA_UNCHANGED))

    def test_already_left_aligned_deletion_is_unchanged(self):
        p, r, a, status, _ = self.align("chrT", 10, "CA", "C")
        self.assertEqual((p, r, a, status), (10, "CA", "C", LA_UNCHANGED))

    def test_insertion_outside_a_repeat_is_unchanged(self):
        p, r, a, status, _ = self.align("chrT", 33, "T", "TGG")
        self.assertEqual((p, r, a, status), (33, "T", "TGG", LA_UNCHANGED))

    def test_symbolic_and_star_alleles_are_not_applicable(self):
        for alt in ("<DEL>", "*", "."):
            with self.subTest(alt=alt):
                self.assertEqual(self.align("chrT", 10, "C", alt)[3], LA_NOT_APPLICABLE)

    # -- the honest failure modes ----------------------------------------
    def test_ref_mismatch_is_reported_not_corrected(self):
        # chrT:5 is A in the reference; this allele claims C.
        p, r, a, status, seen = self.align("chrT", 5, "C", "T")
        self.assertEqual((p, r, a), (5, "C", "T"))
        self.assertEqual(status, LA_REF_MISMATCH)
        self.assertEqual(seen, "A")

    def test_ref_past_the_end_of_the_contig_is_a_mismatch(self):
        p, r, a, status, seen = self.align("chrS", 10, "CC", "C")
        self.assertEqual(status, LA_REF_MISMATCH)
        self.assertIsNone(seen)

    def test_deletion_at_a_contig_start_stops_at_base_one(self):
        p, r, a, status, _ = self.align("chrS", 3, "AA", "A")
        self.assertEqual((p, r, a), (1, "AA", "A"))
        self.assertEqual(status, LA_AT_CONTIG_START)

    def test_contig_absent_from_the_reference(self):
        p, r, a, status, _ = self.align("chrZ", 10, "CA", "C")
        self.assertEqual((p, r, a, status), (10, "CA", "C", LA_CONTIG_ABSENT))

    def test_shift_limit_is_reported(self):
        p, r, a = normalize_allele(17, "AA", "A")
        out = left_align_allele(self.fa, "chrT", p, r, a, max_shift=3)
        self.assertEqual(out[3], LA_SHIFT_LIMIT)
        self.assertGreater(out[0], 10)      # stopped partway, and says so

    def test_lowercase_reference_still_matches(self):
        # chrL is soft-masked in the file; the comparison is case-insensitive.
        self.assertEqual(self.align("chrL", 5, "A", "G")[3], LA_UNCHANGED)


if __name__ == "__main__":
    unittest.main()
