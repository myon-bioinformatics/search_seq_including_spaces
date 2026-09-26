import hashlib
import pathlib
import unittest

from protein.parse import Prepared, SequenceError, prepare, read_fasta

FIXTURES = pathlib.Path(__file__).parent / "fixtures"
FIXTURE_SHA256 = {
    "P00698_mature.fasta": ("eaf1f8d7f5bf78bbefe1100d912a09893361c57dda057f634adb5c5af6dfe79c", 129),
    "P37840.fasta": ("9a06387610edd099466d614941c42a433f444feb6a8f25375bb0b6ad8ec9c6f4", 140),
}


class TestReadFasta(unittest.TestCase):
    def test_multiple_records_and_line_wrapping(self):
        text = ">one\nACDE\nFGHI\n>two desc\nKLMN\n"
        self.assertEqual(read_fasta(text), [("one", "ACDEFGHI"), ("two desc", "KLMN")])

    def test_crlf_blank_lines_and_comments(self):
        text = ">one\r\n\r\n; a comment\r\nAC DE\r\n\r\n"
        self.assertEqual(read_fasta(text), [("one", "AC DE")])

    def test_headerless_sequence_is_one_record(self):
        self.assertEqual(read_fasta("ACDE\nFG\n"), [("", "ACDEFG")])

    def test_data_before_first_header_is_an_error(self):
        with self.assertRaises(SequenceError):
            read_fasta("ACDE\n>one\nFG\n")

    def test_empty_text(self):
        self.assertEqual(read_fasta(""), [])
        self.assertEqual(read_fasta("\n\n; only comments\n"), [])

    def test_header_only_record_is_returned_unvalidated(self):
        self.assertEqual(read_fasta(">empty\n"), [("empty", "")])


class TestPrepare(unittest.TestCase):
    def test_uppercases_and_drops_whitespace_keeping_input_positions(self):
        p = prepare("h", "ac d\te\nf")
        self.assertIsInstance(p, Prepared)
        self.assertEqual(p.residues, "ACDEF")
        self.assertEqual(p.input_positions, (1, 2, 3, 4, 5))
        self.assertEqual(p.classes, ("canonical",) * 5)

    def test_terminal_stop_is_stripped_and_recorded(self):
        p = prepare("h", "ACD*")
        self.assertEqual(p.residues, "ACD")
        self.assertTrue(p.stripped_terminal_stop)
        self.assertFalse(prepare("h", "ACD").stripped_terminal_stop)

    def test_internal_stop_is_an_error_by_default(self):
        with self.assertRaises(SequenceError) as ctx:
            prepare("h", "AC*DE")
        self.assertEqual(ctx.exception.position, 3)
        self.assertEqual(ctx.exception.symbol, "*")

    def test_internal_stop_kept_as_marker_when_allowed(self):
        p = prepare("h", "AC*DE*", allow_internal_stop=True)
        self.assertEqual(p.residues, "AC*DE")
        self.assertEqual(p.classes[2], "stop")
        self.assertTrue(p.stripped_terminal_stop)

    def test_gap_is_an_error_by_default(self):
        with self.assertRaises(SequenceError) as ctx:
            prepare("h", "AC-DE")
        self.assertEqual(ctx.exception.position, 3)

    def test_strip_gaps_keeps_input_positions(self):
        p = prepare("h", "A--CD-E", strip_gaps=True)
        self.assertEqual(p.residues, "ACDE")
        self.assertEqual(p.input_positions, (1, 4, 5, 7))
        self.assertEqual(p.stripped_gaps, 3)

    def test_invalid_symbol_reports_input_position(self):
        for raw, pos, sym in (("ACD1E", 4, "1"), ("AC Ω", 3, "Ω"), ("a.c", 2, ".")):
            with self.assertRaises(SequenceError) as ctx:
                prepare("rec", raw)
            self.assertEqual((ctx.exception.position, ctx.exception.symbol), (pos, sym))
            self.assertEqual(ctx.exception.header, "rec")
            self.assertIn("rec", str(ctx.exception))

    def test_empty_inputs_are_errors(self):
        for raw in ("", "   \n", "*"):
            with self.assertRaises(SequenceError):
                prepare("h", raw)
        with self.assertRaises(SequenceError):
            prepare("h", "---", strip_gaps=True)
        with self.assertRaises(SequenceError):
            prepare("h", "**", allow_internal_stop=True)

    def test_ambiguous_and_rare_residues_are_kept_with_their_class(self):
        p = prepare("h", "ABZJXUO")
        self.assertEqual(p.residues, "ABZJXUO")
        self.assertEqual(p.classes, ("canonical", "ambiguous", "ambiguous", "ambiguous",
                                     "ambiguous", "rare-but-defined", "rare-but-defined"))

    def test_prepared_is_immutable(self):
        p = prepare("h", "ACD")
        with self.assertRaises(AttributeError):
            p.residues = "X"  # type: ignore[misc]


class TestFixtures(unittest.TestCase):
    def test_fixtures_match_recorded_hashes(self):
        for name, (sha, length) in FIXTURE_SHA256.items():
            records = read_fasta((FIXTURES / name).read_text(encoding="utf-8"))
            self.assertEqual(len(records), 1, name)
            p = prepare(*records[0])
            self.assertEqual(len(p.residues), length, name)
            self.assertEqual(hashlib.sha256(p.residues.encode()).hexdigest(), sha, name)
            self.assertEqual(set(p.classes), {"canonical"}, name)
            self.assertEqual(p.input_positions, tuple(range(1, length + 1)), name)

    def test_source_md_lists_every_fixture_hash(self):
        source = (FIXTURES / "SOURCE.md").read_text(encoding="utf-8")
        for name, (sha, _) in FIXTURE_SHA256.items():
            self.assertIn(name, source)
            self.assertIn(sha, source)


if __name__ == "__main__":
    unittest.main()
