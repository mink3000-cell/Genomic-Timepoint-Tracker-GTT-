"""Header parsing, record decomposition and genotype projection."""
import os
import unittest

from gtt.errors import VcfFormatError
from gtt.model import GT_HEMI, GT_HET, GT_HOM_ALT, GT_HOM_REF, GT_NO_CALL
from gtt.vcfio import VcfReader, _parse_gt, read_header

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


class TestHeader(unittest.TestCase):
    def test_synthetic_header_parsed(self):
        h = read_header(os.path.join(FIX, "synthetic_baseline.vcf"))
        self.assertEqual(h.fileformat, "VCFv4.2")
        self.assertEqual(h.samples, ["SYNTHETIC_SUBJECT"])
        self.assertEqual(h.contigs["chr1"], 248956422)
        self.assertTrue(h.chr_prefixed)
        self.assertTrue(h.synthetic)
        self.assertIn("LowQual", h.filters)
        self.assertIn("GT", h.formats)

    def test_unprefixed_contigs_detected(self):
        h = read_header(os.path.join(FIX, "synthetic_followup_nochr.vcf"))
        self.assertFalse(h.chr_prefixed)
        self.assertIn("1", h.contigs)


class TestMalformed(unittest.TestCase):
    def _expect_error(self, name, needle):
        with self.assertRaises(VcfFormatError) as ctx:
            r = VcfReader(os.path.join(FIX, name))
            for _ in r:
                pass
        self.assertIn(needle, str(ctx.exception).lower())

    def test_missing_fileformat(self):
        self._expect_error("malformed_no_fileformat.vcf", "fileformat")

    def test_missing_chrom_line(self):
        self._expect_error("malformed_no_chrom_line.vcf", "#chrom")

    def test_data_before_header(self):
        self._expect_error("malformed_data_before_header.vcf", "before the #chrom")

    def test_non_integer_pos(self):
        self._expect_error("malformed_bad_pos.vcf", "not an integer")

    def test_invalid_ref_alphabet(self):
        self._expect_error("malformed_bad_ref.vcf", "plain [acgtn]")

    def test_column_count_mismatch(self):
        self._expect_error("malformed_column_count.vcf", "header declares")

    def test_missing_file(self):
        with self.assertRaises(VcfFormatError):
            read_header(os.path.join(FIX, "does_not_exist.vcf"))


class TestGenotypeProjection(unittest.TestCase):
    def test_statuses(self):
        self.assertEqual(_parse_gt("0/1", 1)[0], GT_HET)
        self.assertEqual(_parse_gt("1/1", 1)[0], GT_HOM_ALT)
        self.assertEqual(_parse_gt("0/0", 1)[0], GT_HOM_REF)
        self.assertEqual(_parse_gt("./.", 1)[0], GT_NO_CALL)
        self.assertEqual(_parse_gt("1", 1)[0], GT_HEMI)
        self.assertEqual(_parse_gt("1|0", 1)[0], GT_HET)
        self.assertTrue(_parse_gt("1|0", 1)[2])

    def test_multiallelic_projection(self):
        # GT 1/2: one copy of ALT 1 and one of ALT 2 - each is het for its own allele.
        self.assertEqual(_parse_gt("1/2", 1), (GT_HET, 1, False))
        self.assertEqual(_parse_gt("1/2", 2), (GT_HET, 1, False))
        self.assertEqual(_parse_gt("1/2", 3)[0], "other_alt_only")
        self.assertEqual(_parse_gt("2/2", 1)[0], "other_alt_only")


class TestDecomposition(unittest.TestCase):
    def setUp(self):
        self.calls = list(VcfReader(os.path.join(FIX, "synthetic_baseline.vcf")))

    def test_one_row_per_alt_with_star_skipped(self):
        # 10 records; one has 2 real ALTs, one has an ALT plus '*' (skipped).
        self.assertEqual(len(self.calls), 11)
        keys = {(c.contig, c.pos, c.ref, c.alt) for c in self.calls}
        self.assertIn(("chr2", 500000, "G", "A"), keys)
        self.assertIn(("chr2", 500000, "G", "T"), keys)
        self.assertNotIn(("chr2", 800000, "G", "*"), keys)

    def test_stats(self):
        r = VcfReader(os.path.join(FIX, "synthetic_baseline.vcf"))
        list(r)
        self.assertEqual(r.stats.records, 10)
        self.assertEqual(r.stats.alleles_emitted, 11)
        self.assertEqual(r.stats.multiallelic_records, 2)
        self.assertEqual(r.stats.skipped_star, 1)
        self.assertEqual(r.stats.contigs_seen, ["chr1", "chr2", "chrX"])

    def test_keep_star_option(self):
        r = VcfReader(os.path.join(FIX, "synthetic_baseline.vcf"), skip_star=False)
        calls = list(r)
        self.assertEqual(len(calls), 12)

    def test_original_representation_preserved(self):
        r = VcfReader(os.path.join(FIX, "synthetic_followup.vcf"))
        by_key = {c.key: c for c in r}
        trimmed = by_key[("chr1", 300000, "CTT", "C")]
        self.assertEqual((trimmed.orig_ref, trimmed.orig_alt), ("CTTA", "CA"))
        self.assertEqual(trimmed.var_type, "deletion")

    def test_haploid_and_filters(self):
        by_key = {c.key: c for c in self.calls}
        self.assertEqual(by_key[("chrX", 1000000, "A", "T")].gt_status, GT_HEMI)
        self.assertTrue(all(c.is_pass for c in self.calls))

    def test_unknown_sample_rejected(self):
        with self.assertRaises(VcfFormatError):
            VcfReader(os.path.join(FIX, "synthetic_baseline.vcf"), sample="NO_SUCH_SAMPLE")


if __name__ == "__main__":
    unittest.main()
