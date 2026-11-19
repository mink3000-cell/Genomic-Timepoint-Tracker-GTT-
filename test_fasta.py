"""Indexed FASTA access: index building, index parsing, random access.

Expected sequences are hand-written here, not read from the fixture, so a
corrupted fixture fails the test instead of silently redefining it.
"""
import gzip
import os
import shutil
import tempfile
import unittest

from gtt.errors import FastaError
from gtt.fasta import IndexedFasta, build_fai, read_fai, write_fai

FIX = os.path.join(os.path.dirname(__file__), "fixtures")
REF = os.path.join(FIX, "synthetic_reference.fa")

# The three fixture contigs, written out by hand.
CHR_T = "ACGTACGTAC" + "AAAAAAAA" + "GT" + "CACACACACACA" + "TTGGCCAA"   # 40 bp
CHR_S = "AAAAACCCCC"                                                      # 10 bp
CHR_L = "ACGTACGTACGTACGTACGT"                                            # 20 bp, lowercase in file


class TestBuildIndex(unittest.TestCase):
    def test_record_fields(self):
        recs = {r.name: r for r in build_fai(REF)}
        self.assertEqual(set(recs), {"chrT", "chrS", "chrL"})
        self.assertEqual(recs["chrT"].length, 40)
        self.assertEqual(recs["chrS"].length, 10)
        self.assertEqual(recs["chrL"].length, 20)
        # linewidth = linebases + the newline byte(s) actually used in the file
        for rec in recs.values():
            self.assertGreaterEqual(rec.linewidth, rec.linebases + 1)

    def test_offset_points_at_the_first_base(self):
        recs = {r.name: r for r in build_fai(REF)}
        with open(REF, "rb") as fh:
            for name, expected in (("chrT", CHR_T), ("chrS", CHR_S)):
                fh.seek(recs[name].offset)
                self.assertEqual(fh.read(len(expected[:10])).decode(), expected[:10])

    def test_roundtrip_through_a_written_index(self):
        with tempfile.TemporaryDirectory() as td:
            fai = os.path.join(td, "x.fai")
            recs = build_fai(REF)
            write_fai(recs, fai)
            again = read_fai(fai)
            self.assertEqual([(r.name, r.length, r.offset, r.linebases, r.linewidth) for r in recs],
                             [(r.name, r.length, r.offset, r.linebases, r.linewidth) for r in again])

    def test_index_is_created_when_missing(self):
        with tempfile.TemporaryDirectory() as td:
            copy = os.path.join(td, "ref.fa")
            shutil.copy(REF, copy)
            self.assertFalse(os.path.exists(copy + ".fai"))
            with IndexedFasta(copy) as fa:
                self.assertTrue(fa.index_created)
                self.assertTrue(fa.index_written)
                self.assertEqual(fa.fetch("chrT", 0, 40), CHR_T)
            self.assertTrue(os.path.exists(copy + ".fai"))
            # Second open reuses the file instead of rebuilding it.
            with IndexedFasta(copy) as fa2:
                self.assertFalse(fa2.index_created)

    def test_crlf_fasta_is_indexed_correctly(self):
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "crlf.fa")
            with open(path, "wb") as fh:
                fh.write(b">c1 synthetic\r\nACGTACGTAC\r\nTTTTGGGGCC\r\nAC\r\n")
            with IndexedFasta(path) as fa:
                self.assertEqual(fa.records["c1"].linewidth, 12)
                self.assertEqual(fa.fetch("c1", 0, 22), "ACGTACGTACTTTTGGGGCCAC")
                self.assertEqual(fa.fetch("c1", 8, 12), "ACTT")

    def test_gzip_fasta_is_refused_with_advice(self):
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "r.fa.gz")
            with gzip.open(path, "wb") as fh:
                fh.write(b">c1\nACGT\n")
            with self.assertRaises(FastaError) as cm:
                build_fai(path)
            self.assertIn("Decompress", str(cm.exception))

    def test_ragged_lines_are_refused(self):
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "ragged.fa")
            with open(path, "w", newline="\n") as fh:
                fh.write(">c1\nACGTAC\nAC\nACGTAC\n")   # short line in the middle
            with self.assertRaises(FastaError):
                build_fai(path)

    def test_empty_and_missing_inputs_are_refused(self):
        with tempfile.TemporaryDirectory() as td:
            empty = os.path.join(td, "e.fa")
            open(empty, "w").close()
            with self.assertRaises(FastaError):
                build_fai(empty)
            with self.assertRaises(FastaError):
                build_fai(os.path.join(td, "nope.fa"))


class TestFetch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fa = IndexedFasta(REF)

    @classmethod
    def tearDownClass(cls):
        cls.fa.close()

    def test_whole_contigs(self):
        self.assertEqual(self.fa.fetch("chrT", 0, 40), CHR_T)
        self.assertEqual(self.fa.fetch("chrS", 0, 10), CHR_S)

    def test_fetch_crosses_a_line_boundary(self):
        # The fixture wraps at 20 bases, so 15..25 spans two lines.
        self.assertEqual(self.fa.fetch("chrT", 15, 25), CHR_T[15:25])

    def test_soft_masked_sequence_is_uppercased(self):
        self.assertEqual(self.fa.fetch("chrL", 0, 20), CHR_L)

    def test_single_base_is_one_based(self):
        self.assertEqual(self.fa.base("chrT", 1), "A")
        self.assertEqual(self.fa.base("chrT", 11), "A")
        self.assertEqual(self.fa.base("chrT", 10), "C")
        self.assertEqual(self.fa.base("chrT", 40), "A")

    def test_empty_slice(self):
        self.assertEqual(self.fa.fetch("chrT", 5, 5), "")

    def test_contig_name_aliases(self):
        self.assertEqual(self.fa.resolve("chrT"), "chrT")
        self.assertEqual(self.fa.resolve("T"), "chrT")
        self.assertIsNone(self.fa.resolve("chrZ"))
        self.assertEqual(self.fa.length("T"), 40)

    def test_out_of_range_is_refused(self):
        with self.assertRaises(FastaError):
            self.fa.fetch("chrT", 0, 41)
        with self.assertRaises(FastaError):
            self.fa.fetch("chrT", -1, 5)
        with self.assertRaises(FastaError):
            self.fa.fetch("chrZ", 0, 5)

    def test_index_that_runs_past_the_end_of_the_file_is_detected(self):
        # An index claiming more sequence than the file holds is caught when the
        # read comes up short. (An index that is merely wrong in the middle of a
        # file cannot be caught this way - hence the mtime warning in IndexedFasta.)
        with tempfile.TemporaryDirectory() as td:
            copy = os.path.join(td, "ref.fa")
            shutil.copy(REF, copy)
            recs = build_fai(copy)
            recs[-1].length += 50                      # last contig only
            write_fai(recs, copy + ".fai")
            with IndexedFasta(copy) as fa:
                with self.assertRaises(FastaError) as cm:
                    fa.fetch("chrL", 40, 60)
                self.assertIn("stale", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
