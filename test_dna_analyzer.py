import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dna_analyzer import (
    clean_sequence, validate_sequence, InvalidSequenceError,
    nucleotide_frequency, gc_content, complement, reverse_complement,
    transcribe, translate, find_motif, parse_fasta,
)


class TestDNAAnalyzer(unittest.TestCase):

    def test_clean_sequence(self):
        self.assertEqual(clean_sequence(" atg c\n"), "ATGC")

    def test_validate_sequence_valid(self):
        self.assertTrue(validate_sequence("ATGC"))

    def test_validate_sequence_invalid(self):
        with self.assertRaises(InvalidSequenceError):
            validate_sequence("ATXZ")

    def test_validate_sequence_empty(self):
        with self.assertRaises(InvalidSequenceError):
            validate_sequence("")

    def test_nucleotide_frequency(self):
        self.assertEqual(nucleotide_frequency("AATTGGCC"), {"A": 2, "T": 2, "G": 2, "C": 2})

    def test_gc_content(self):
        self.assertEqual(gc_content("GGCC"), 100.0)
        self.assertEqual(gc_content("AATT"), 0.0)

    def test_complement(self):
        self.assertEqual(complement("ATGC"), "TACG")

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")

    def test_transcribe(self):
        self.assertEqual(transcribe("ATGC"), "AUGC")

    def test_translate_stops_at_stop_codon(self):
        # ATG GCC TAA -> M A * (stop) -> "MA"
        self.assertEqual(translate("ATGGCCTAA"), "MA")

    def test_find_motif(self):
        self.assertEqual(find_motif("ATGATGATG", "ATG"), [0, 3, 6])
        self.assertEqual(find_motif("ATGC", "TTT"), [])

    def test_parse_fasta(self):
        path = os.path.join(os.path.dirname(__file__), "..", "sample_data", "sample.fasta")
        records = parse_fasta(path)
        self.assertIn("sample_gene_1 example DNA sequence", records)
        self.assertTrue(records["sample_gene_1 example DNA sequence"].startswith("ATG"))


if __name__ == "__main__":
    unittest.main()
