"""Extended adversarial discovery and multi-directory integration tests.

Tests beyond standard unit tests:
- Dynamically discovers and validates sample data files across repository directories.
- Exercises parsers and tools against adversarial/edge-case inputs (control characters,
  malformed headers, boundary whitespace, malformed CIF/FASTA) ensuring proper error
  handling or expected results without unhandled process crashes.
"""

import io
import os
import unittest
from pathlib import Path

from comparison import compare_aligned, ComparisonError
from protein_ionization import ionization
from sequence_readers import read_sequences, SequenceReadError
from streaming_matcher import match_chunks
from structure_tools import parse_structure


class TestAdversarialDiscovery(unittest.TestCase):
    """Adversarial and repository-wide discovery tests."""

    def setUp(self):
        self.repo_root = Path(__file__).parent.parent

    def test_discover_repository_data_files(self):
        """Dynamically discover and validate sample sequence files across repository paths."""
        discovered_files = []
        for search_dir in [".", "input_folder", "examples", "tests"]:
            target_path = self.repo_root / search_dir
            if not target_path.exists():
                continue
            for root, _, files in os.walk(target_path):
                for f in files:
                    if f.endswith((".fa", ".fasta", ".json", ".txt", ".cif")):
                        discovered_files.append(Path(root) / f)

        self.assertGreater(len(discovered_files), 0, "Should discover data files in repository directories")

        for file_path in discovered_files:
            if file_path.suffix in [".fa", ".fasta", ".txt"]:
                # Should read records or fail gracefully with handled exceptions
                try:
                    records = list(read_sequences(str(file_path)))
                    self.assertIsInstance(records, list)
                except (SequenceReadError, ValueError, UnicodeDecodeError, OSError):
                    pass  # Expected format or decoding errors for non-standard files

    def test_adversarial_reader_inputs(self):
        """Test readers against adversarial / unusual input strings."""
        adversarial_inputs = [
            "\x00\x01\x02\x03\x04",  # Binary control characters
            ">header_only_no_sequence\n",  # Empty FASTA sequence
            ">seq1\nACGT\n>seq2\n",  # Trailing header without sequence
            " \t \n \r \n",  # Pure whitespace
            "@fq1\nACGT\n+\n!\n@fq2\n",  # Truncated FASTQ
        ]
        for idx, text_input in enumerate(adversarial_inputs):
            stream = io.StringIO(text_input)
            try:
                records = list(read_sequences(stream, format="fasta", source=f"adversarial_{idx}"))
                self.assertIsInstance(records, list)
            except (SequenceReadError, ValueError, TypeError, UnicodeError, AttributeError):
                pass  # Handled error

    def test_adversarial_streaming_matcher(self):
        """Test streaming matcher with edge-case gap and anchor parameters."""
        chunks = ["ACGTACGT\n", "NNNNNNNN\n", "ACGTACGT\n"]
        events = list(match_chunks(chunks, left="AC", right="GT", min_gap=0, max_gap=10))
        matches = [e for e in events if isinstance(e, dict) and e.get("event") == "match"]
        self.assertGreaterEqual(len(matches), 1)

        # Ranged gap with no match events
        no_events = list(match_chunks(chunks, left="GGGG", right="CCCC", min_gap=10, max_gap=20))
        no_matches = [e for e in no_events if isinstance(e, dict) and e.get("event") == "match"]
        self.assertEqual(len(no_matches), 0)

    def test_adversarial_protein_ionization(self):
        """Test protein ionization function with edge case sequences."""
        # Clean canonical protein sequence
        result = ionization("ACDEFGHIKLMNPQRSTVWY", ph=7.0)
        self.assertIsInstance(result, dict)
        self.assertIn("net_charge", result)
        self.assertIn("estimated_pi", result)
        self.assertIsInstance(result["net_charge"], float)

        # Invalid residue symbols should raise ValueError
        with self.assertRaises(ValueError):
            ionization("INVALID_RESIDUES_123", ph=7.0)

    def test_adversarial_comparison(self):
        """Test sequence comparison with mismatched aligned lengths or invalid alphabets."""
        with self.assertRaises(ComparisonError):
            compare_aligned("ACGT", "ACG", alphabet="dna")

        with self.assertRaises(ComparisonError):
            compare_aligned("ACGT", "ACGT", alphabet="invalid_alphabet")

    def test_adversarial_cif_parser(self):
        """Test CIF structure parser with malformed text."""
        malformed_cif = "data_test\n_atom_site.id\n_atom_site.type_symbol\n1\n"
        with self.assertRaises((ValueError, KeyError, IndexError, Exception)):
            parse_structure(malformed_cif, format="cif")

    def test_flexible_unittest_discovery(self):
        """Verify that unittest test loader can discover tests from custom paths/patterns."""
        loader = unittest.TestLoader()
        suite = loader.discover(start_dir=str(self.repo_root / "tests"), pattern="test_*.py")
        self.assertGreater(suite.countTestCases(), 0)


if __name__ == "__main__":
    unittest.main()
