"""Adversarial contracts and repository-wide sequence fixture integration.

Readers validate record structure, not biological alphabets. Each case states
its output or exact rejection; unexpected implementation exceptions propagate.
The custom character map below is a caller-side simulation, not a reader codec.
"""

import hashlib
import io
import os
import unittest
from pathlib import Path

from comparison import compare_aligned, ComparisonError
from protein_ionization import ionization
from sequence_readers import read_sequences, SequenceReadError
from streaming_matcher import match_chunks
from structure_tools import Atom, parse_structure


# Residue hashes/lengths are recorded independently in tests/fixtures/SOURCE.md.
KNOWN_FASTA = {
    "tests/fixtures/P00698_mature.fasta": (
        "sp|P00698|LYSC_CHICK", 129,
        "eaf1f8d7f5bf78bbefe1100d912a09893361c57dda057f634adb5c5af6dfe79c"),
    "tests/fixtures/P37840.fasta": (
        "sp|P37840|SYUA_HUMAN", 140,
        "9a06387610edd099466d614941c42a433f444feb6a8f25375bb0b6ad8ec9c6f4"),
    "tests/fixtures/P37840_bom.fasta": (
        "sp|P37840|SYUA_HUMAN", 140,
        "9a06387610edd099466d614941c42a433f444feb6a8f25375bb0b6ad8ec9c6f4"),
}
CIF_HEADER = """data_adversarial
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
"""
CIF_ROW = "1 C CA ALA A 1 0 1 2\n"


class TestAdversarialDiscovery(unittest.TestCase):
    """Adversarial and repository-wide discovery tests."""

    def setUp(self):
        self.repo_root = Path(__file__).resolve().parent.parent

    def test_discover_repository_data_files(self):
        """Walk the repository once and require known fixtures to parse correctly."""
        sequence_suffixes = {
            ".fa", ".fas", ".fasta", ".fna", ".ffn", ".faa", ".frn",
            ".fq", ".fastq", ".gb", ".gbk", ".genbank", ".gp",
        }
        # Skip generated/tool trees, so local environments and CI reports do not
        # become sample inputs. Root, input_folder/, examples/ and tests/ are
        # all covered by this single traversal without duplicate visits.
        excluded_dirs = {
            ".git", ".venv", "venv", "__pycache__", ".pytest_cache",
            "node_modules", "reports", "build", "dist",
        }
        discovered = set()
        for root, dirs, files in os.walk(self.repo_root):
            dirs[:] = sorted(name for name in dirs if name not in excluded_dirs)
            for name in files:
                path = Path(root) / name
                suffixes = path.suffixes
                if suffixes and suffixes[-1].lower() in {".gz", ".bz2", ".xz", ".lzma"}:
                    suffixes = suffixes[:-1]
                # Structure parsing takes decoded text; compressed sequence
                # files alone use read_sequences' compression contract.
                if path.suffix.lower() == ".cif" or (suffixes and suffixes[-1].lower() in sequence_suffixes):
                    discovered.add(path.relative_to(self.repo_root).as_posix())
        # JSON configs, requirements.txt and golden ROI prose are not sequence
        # fixtures: those ambiguous formats require explicit fixture options.
        self.assertTrue(set(KNOWN_FASTA).issubset(discovered),
                        f"Missing fixtures: {set(KNOWN_FASTA) - discovered}")
        self.assertGreaterEqual(len(discovered), 3)
        parsed = {}
        for relative in sorted(discovered):
            with self.subTest(path=relative):
                path = self.repo_root / relative
                if path.suffix.lower() == ".cif":
                    structure = parse_structure(path.read_text(encoding="utf-8"), format="cif")
                    self.assertGreater(len(structure["atoms"]), 0)
                    continue
                records = list(read_sequences(path))
                self.assertGreater(len(records), 0, f"No records in {relative}")
                for index, record in enumerate(records, 1):
                    self.assertTrue(record.id)
                    self.assertTrue(record.sequence)
                    self.assertEqual(record.source, str(path))
                    self.assertEqual(record.provenance.record_index, index)
                parsed[relative] = records
                if relative in KNOWN_FASTA:
                    identifier, length, digest = KNOWN_FASTA[relative]
                    self.assertEqual(len(records), 1)
                    self.assertEqual(records[0].id, identifier)
                    self.assertEqual(len(records[0].sequence), length)
                    self.assertEqual(hashlib.sha256(records[0].sequence.encode("ascii")).hexdigest(), digest)
                    self.assertEqual(records[0].provenance.format, "fasta")
                    self.assertIsNone(records[0].quality)
        self.assertTrue(set(KNOWN_FASTA).issubset(parsed))
        self.assertEqual(parsed["tests/fixtures/P37840.fasta"][0].sequence,
                         parsed["tests/fixtures/P37840_bom.fasta"][0].sequence)

    def test_adversarial_reader_inputs(self):
        """Structural errors carry the documented type, reason, source and line."""
        cases = [
            ("\x00\x01\x02\x03\x04", "fasta", "sequence before FASTA header", 1),
            (">header_only_no_sequence\n", "fasta", "empty FASTA sequence", 1),
            (">\nAC\n", "fasta", "empty FASTA header", 1),
            (" >id\nAC\n", "fasta", "sequence before FASTA header", 1),
            ("@fq1\nACGT\n+\n!\n", "fastq", "FASTQ quality-length mismatch", 1),
            ("@fq1\nACGT\n+\n!\n@fq2\n", "fastq", "FASTQ quality-length mismatch", 1),
            ("@fq1\nACGT\n", "fastq", "missing FASTQ + separator", 1),
            ("@fq1\nA C\n+\n!!!\n", "fastq", "empty/whitespace FASTQ sequence line", 2),
            ("@fq1\nAC\n+\n!\x00\n", "fastq", "invalid FASTQ quality (expected ASCII 33..126)", 4),
            ("@fq1\nAC\n+other\n!!\n", "fastq", "FASTQ + header does not match @ header", 3),
        ]
        for index, (text, format, reason, line) in enumerate(cases):
            source = f"adversarial_{index}"
            with self.subTest(text=text, format=format):
                with self.assertRaises(SequenceReadError) as caught:
                    list(read_sequences(io.StringIO(text), format=format, source=source))
                error = caught.exception
                self.assertIs(type(error), SequenceReadError)
                self.assertEqual((error.reason, error.source, error.line), (reason, source, line))
                self.assertEqual(str(error), f"{source}: line {line}: {reason}")

    def test_adversarial_reader_valid_outputs(self):
        for format in ("fasta", "fastq", "text"):
            for text in ("", " \t \n \r \n"):
                with self.subTest(format=format, text=text):
                    self.assertEqual(list(read_sequences(io.StringIO(text), format=format)), [])
        # The reader preserves symbols, including NUL and stop markers, after
        # structural validation. Alphabet validation belongs to downstream tools.
        text = "\ufeff;comment\n>id description\n a C\x00P* \n;ignored\n"
        row, = read_sequences(io.StringIO(text), format="fasta", source="valid")
        self.assertEqual((row.id, row.description, row.sequence, row.quality, row.source),
                         ("id", "description", "aC\x00P*", None, "valid"))
        self.assertEqual((row.provenance.start_line, row.provenance.end_line), (2, 4))
        row, = read_sequences(io.StringIO("@id\nACGT\n+\n@+!~"), format="fastq")
        self.assertEqual((row.id, row.sequence, row.quality), ("id", "ACGT", "@+!~"))
        self.assertEqual((row.provenance.start_line, row.provenance.end_line), (1, 4))

    def test_reader_yields_valid_prefix_before_trailing_error(self):
        records = read_sequences(io.StringIO(">seq1\nACGT\n>seq2\n"), format="fasta", source="trailing")
        first = next(records)
        self.assertEqual((first.id, first.sequence, first.provenance.record_index), ("seq1", "ACGT", 1))
        with self.assertRaises(SequenceReadError) as caught:
            next(records)
        self.assertEqual((caught.exception.reason, caught.exception.source, caught.exception.line),
                         ("empty FASTA sequence", "trailing", 3))

    def test_custom_character_table_encoding_simulation(self):
        """Caller translation stops at either marker, including at payload edges."""
        custom_charmap = {
            "ア": "A", "イ": "C", "ウ": "G", "エ": "T",
            "オ": "N", "カ": "A", "キ": "C", "ク": "G",
            "\x50": None, "\x00": None,
        }
        cases = [
            ("アイウエオカキク\x50余剰ペロード", "ACGTNACG"),
            ("アイウエオカキク\x00アイ", "ACGTNACG"),
            ("アイ", "AC"),  # EOF is also a valid payload boundary.
            ("ア?イ", "AC"),  # Unmapped caller-side symbols are omitted.
            ("", ""),
        ]
        for terminator in ("\x50", "\x00"):
            cases.extend([(terminator + "アイ", ""), ("アイ" + terminator, "AC"),
                          ("ア" + terminator + "イ", "A"),
                          ("ア" + terminator + terminator + "イ", "A")])
        for payload, expected in cases:
            with self.subTest(payload=payload):
                symbols = []
                for char in payload:
                    mapped = custom_charmap.get(char, "")
                    if mapped is None:
                        break
                    symbols.append(mapped)
                translated = "".join(symbols)
                self.assertEqual(translated, expected)
                stream = io.StringIO(f">custom_charmap_record\n{translated}\n")
                if not expected:
                    # A terminator-first/empty payload is not a valid FASTA
                    # record: conversion must not silently invent residues.
                    with self.assertRaises(SequenceReadError) as caught:
                        list(read_sequences(stream, format="fasta", source="charmap_test"))
                    self.assertEqual((caught.exception.reason, caught.exception.source, caught.exception.line),
                                     ("empty FASTA sequence", "charmap_test", 1))
                else:
                    row, = read_sequences(stream, format="fasta", source="charmap_test")
                    self.assertEqual((row.id, row.sequence, row.source, row.quality),
                                     ("custom_charmap_record", expected, "charmap_test", None))
                    self.assertEqual((row.provenance.start_line, row.provenance.end_line), (1, 2))

    def test_adversarial_streaming_matcher(self):
        """A match survives every chunk split, including zero-length gaps."""
        text = "ACNNGT"
        for split in range(len(text) + 1):
            with self.subTest(split=split):
                events = list(match_chunks([text[:split], text[split:]], left="AC", right="GT",
                                           min_gap=2, max_gap=2))
                self.assertEqual([e["event"] for e in events], ["match", "complete"])
                self.assertEqual(events[0]["match"]["span"], [1, 6])
                self.assertEqual(events[0]["match"]["middle"], "NN")
                self.assertEqual(events[0]["match"]["gap"], 2)
                self.assertEqual(events[1]["match_count"], 1)
        events = list(match_chunks(["AC", "GT"], left="AC", right="GT", min_gap=0, max_gap=0))
        self.assertEqual(events[0]["match"]["gap"], 0)
        self.assertEqual(events[-1]["match_count"], 1)
        no_events = list(match_chunks([text], left="GGGG", right="CCCC", min_gap=10, max_gap=20))
        self.assertEqual([e["event"] for e in no_events], ["complete"])
        self.assertEqual(no_events[0]["match_count"], 0)

    def test_adversarial_protein_ionization(self):
        # For alanine only the free termini contribute. Their pKa values are
        # 7.59 and 3.55, so equal positive/negative charge occurs at 5.57.
        result = ionization("A", ph=5.57)
        self.assertEqual(result["schema"], "protein-ionization/1")
        self.assertAlmostEqual(result["net_charge"], 0.0)
        self.assertAlmostEqual(result["estimated_pi"], 5.57)
        with self.assertRaisesRegex(ValueError, "^sequence must contain ASCII letters and whitespace only$"):
            ionization("INVALID_RESIDUES_123", ph=7.0)
        with self.assertRaisesRegex(ValueError, "^ionization requires a nonempty canonical protein sequence$"):
            ionization("AX", ph=7.0)

    def test_adversarial_comparison(self):
        result = compare_aligned("ACGT", "ACGA", alphabet="dna")
        self.assertEqual(result["percent_identity"], 75.0)
        self.assertEqual(result["counts"]["matches"], 3)
        self.assertEqual(result["counts"]["mismatches"], 1)
        with self.assertRaisesRegex(ComparisonError, "aligned sequences must have equal length"):
            compare_aligned("ACGT", "ACG", alphabet="dna")
        with self.assertRaisesRegex(ComparisonError, "^alphabet must be one of dna, rna, protein$"):
            compare_aligned("ACGT", "ACGT", alphabet="invalid_alphabet")

    def test_adversarial_cif_parser(self):
        cases = [
            ("data_test\n_atom_site.id\n_atom_site.type_symbol\n1\n", "missing CIF scalar value"),
            (CIF_HEADER + "1 C\n", "incomplete CIF loop row"),
            (CIF_HEADER + "1 C\nstop_\n", "incomplete CIF loop row"),
            (CIF_HEADER + "1 C CA ALA A 1 nan 1 2\n", "nonfinite atom coordinates"),
            ("data_test\n_title 'unclosed\n", "unterminated CIF quoted value"),
            ("data_test\n_title\n;unclosed\n", "unterminated CIF multiline value"),
            ("data_test\n_title\n;value\n;extra\n", "text after CIF multiline terminator"),
            ("data_test\n", "structure has no supported atomic coordinates"),
        ]
        for text, reason in cases:
            with self.subTest(text=text):
                with self.assertRaises(ValueError) as caught:
                    parse_structure(text, format="cif")
                self.assertIs(type(caught.exception), ValueError)
                self.assertEqual(str(caught.exception), reason)

    def test_cif_valid_boundaries_and_terminators(self):
        # A complete loop row may end at EOF, # comment, or stop_. Quoted
        # delimiter-looking tokens are values, not loop boundaries.
        for ending in ("", "#\n", "stop_\n"):
            text = CIF_HEADER + CIF_ROW + ending
            with self.subTest(ending=ending):
                structure = parse_structure(text, format="cif")
                self.assertEqual(structure["atoms"], [Atom("1", "CA", "C", "ALA", "A", "1", "1", "", 0., 1., 2.)])
                self.assertEqual(structure["format"], "cif")
                self.assertEqual(structure["source_sha256"], hashlib.sha256(text.encode()).hexdigest())
        text = "data_adversarial\n_title\n;metadata\n;\n" + CIF_HEADER.split("\n", 1)[1]
        text += CIF_ROW.replace(" CA ", " 'stop_' ") + "stop_\n"
        self.assertEqual(parse_structure(text)["atoms"][0].name, "stop_")

    def test_flexible_unittest_discovery(self):
        loader = unittest.TestLoader()
        suite = loader.discover(start_dir=str(self.repo_root / "tests"), pattern="test_*.py")
        self.assertEqual(loader.errors, [])
        self.assertGreater(suite.countTestCases(), 0)


if __name__ == "__main__":
    unittest.main()
