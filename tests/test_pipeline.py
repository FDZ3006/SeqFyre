"""
test_pipeline.py — pengujian unit pipeline SeqFyre.

Jalankan dari root proyek:
    python -m unittest discover -s tests -v

Pengujian sengaja memakai parser fallback (use_biopython=False) agar
hasilnya deterministik dan independen dari ketersediaan BioPython.
"""

import os
import sys
import unittest

# Pastikan paket 'seqfyre' dapat diimpor saat tes dijalankan dari mana pun.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from seqfyre import Analyzer, Parser, ParserError  # noqa: E402
from seqfyre.models import SequenceRecord  # noqa: E402

FASTA = (
    ">seq1 Thermus thermophilus HB8 16S\n"
    "GGGGCCCCATAT\n"          # G4 C4 A2 T2  -> GC = 8/12 = 66.67%
    ">seq2 Escherichia coli 16S\n"
    "AAAATTTTGGCC\n"          # A4 T4 G2 C2  -> GC = 4/12 = 33.33%
)

FASTQ = (
    "@seq1 Thermus thermophilus\n"
    "GGGGCCCCATAT\n"
    "+\n"
    "IIIIIIIIIIII\n"
)


class TestSequenceRecord(unittest.TestCase):
    def setUp(self):
        self.rec = SequenceRecord("NR_001.1", "NR_001.1 Thermus thermophilus HB8", "GGGGCCCCATAT")

    def test_length(self):
        self.assertEqual(self.rec.length, 12)

    def test_frequency_is_dict(self):
        freq = self.rec.nucleotide_frequency()
        self.assertIsInstance(freq, dict)
        self.assertEqual(freq["G"], 4)
        self.assertEqual(freq["C"], 4)
        self.assertEqual(freq["A"], 2)
        self.assertEqual(freq["T"], 2)

    def test_gc_content(self):
        self.assertAlmostEqual(self.rec.gc_content, 66.6667, places=3)

    def test_organism_extraction(self):
        self.assertEqual(self.rec.organism, "Thermus thermophilus")

    def test_n_excluded_from_gc(self):
        rec = SequenceRecord("x", "x org name", "GGCCNNNN")  # GC=4 valid=4 -> 100%
        self.assertEqual(rec.nucleotide_frequency()["N"], 4)
        self.assertAlmostEqual(rec.gc_content, 100.0)


class TestParser(unittest.TestCase):
    def setUp(self):
        self.parser = Parser(use_biopython=False)  # paksa jalur fallback

    def test_parse_fasta(self):
        recs = self.parser.parse_upload("x.fasta", FASTA.encode())
        self.assertEqual(len(recs), 2)
        self.assertEqual(recs[0].id, "seq1")

    def test_parse_fastq(self):
        recs = self.parser.parse_upload("x.fastq", FASTQ.encode())
        self.assertEqual(len(recs), 1)
        self.assertEqual(recs[0].sequence, "GGGGCCCCATAT")

    def test_empty_raises(self):
        with self.assertRaises(ParserError):
            self.parser.parse_upload("x.fasta", b"   ")


class TestAnalyzer(unittest.TestCase):
    def setUp(self):
        recs = Parser(use_biopython=False).parse_upload("x.fasta", FASTA.encode())
        self.an = Analyzer(recs)

    def test_sort_descending(self):
        ordered = self.an.sort_by_gc(descending=True)
        self.assertEqual(ordered[0].id, "seq1")  # GC tertinggi di atas
        self.assertGreater(ordered[0].gc_content, ordered[1].gc_content)

    def test_top_n_clamped(self):
        # Minta 10 padahal cuma 2 sekuens -> tetap kembalikan 2.
        self.assertEqual(len(self.an.top_n(10)), 2)

    def test_csv_has_header_and_rows(self):
        csv = self.an.to_csv_string()
        self.assertIn("Accession_ID", csv)
        self.assertIn("seq1", csv)

    def test_plots_generated(self):
        plots = self.an.plots_base64()
        self.assertEqual(set(plots), {"histogram_gc", "barchart_gc", "composition", "scatter_len_gc"})
        for uri in plots.values():
            self.assertTrue(uri.startswith("data:image/png;base64,"))

    def test_summary_counts(self):
        s = self.an.summary()
        self.assertEqual(s["total_sequences"], 2)


if __name__ == "__main__":
    unittest.main()
