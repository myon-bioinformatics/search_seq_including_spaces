import pathlib
import unittest

from protein.features import protein_profiles
from protein.parse import prepare, read_fasta
from protein.render import render_ascii
from protein.roi import regions_of_interest

FIXTURES = pathlib.Path(__file__).parent / "fixtures"
MARGIN = 8


def load(name):
    return prepare(*read_fasta((FIXTURES / name).read_text(encoding="utf-8"))[0])


LYSOZYME = load("P00698_mature.fasta")
SYNUCLEIN = load("P37840.fasta")
# 1161 residues: a length whose last ruler label ("1161") once got cut to "1".
LONG = prepare("long", LYSOZYME.residues * 9)


def result(seq, **kwargs):
    return regions_of_interest(protein_profiles(seq), **kwargs)


def alignment_errors(text, regions, width):
    """Ruler labels on their residue column, '['/']' on region ends, lines fit."""
    errors = []
    lines = text.splitlines()
    starts = {r["start"] for r in regions}
    ends = {r["end"] for r in regions}
    seq_rows = [k for k, line in enumerate(lines) if line.startswith("SEQ     ")]
    for block, k in enumerate(seq_rows):
        first = block * width + 1
        ruler, roi = lines[k - 2], lines[k + 3]
        if not roi.startswith("ROI"):
            errors.append(f"line {k + 4}: no ROI track")
        col = MARGIN
        while col < len(ruler):
            if ruler[col] != " ":
                label = ruler[col:].split(" ")[0]
                if int(label) != first + col - MARGIN:
                    errors.append(f"ruler label {label} at residue {first + col - MARGIN}")
                if col - MARGIN + len(label) > len(lines[k]) - MARGIN:
                    errors.append(f"ruler label {label} runs past the block")
                col += len(label)
            else:
                col += 1
        for col, ch in enumerate(roi[MARGIN:]):
            position = first + col
            if ch == "[" and position not in starts:
                errors.append(f"'[' at {position}")
            if ch == "]" and position not in ends:
                errors.append(f"']' at {position}")
        for position in starts:
            if first <= position < first + width and roi[MARGIN + position - first:][:1] != "[":
                errors.append(f"no '[' at {position}")
    errors += [f"line {k + 1} too long" for k, line in enumerate(lines) if len(line) > width + MARGIN]
    return errors


class TestLayout(unittest.TestCase):
    def test_alignment_at_several_widths(self):
        for seq in (LYSOZYME, SYNUCLEIN, LONG):
            doc = result(seq)
            for width in (20, 40, 60, 80):
                text = render_ascii(doc, width=width)
                self.assertEqual(alignment_errors(text, doc["regions"], width), [],
                                 (len(seq.residues), width))

    def test_alignment_at_every_width_for_alpha_synuclein(self):
        # Covers regions that end exactly at, or one past, a block boundary.
        doc = result(SYNUCLEIN)
        for width in range(20, 141):
            text = render_ascii(doc, width=width)
            self.assertEqual(alignment_errors(text, doc["regions"], width), [], width)

    def test_long_sequence_has_regions_and_a_complete_last_label(self):
        doc = result(LONG)
        self.assertGreater(len(doc["regions"]), 2)
        text = render_ascii(doc, width=60)
        self.assertIn("1141", text)  # the last label that fits (1141-1161 block)

    def test_plain_ascii_without_colour(self):
        for seq in (LYSOZYME, SYNUCLEIN):
            text = render_ascii(result(seq))
            self.assertNotIn("\x1b", text)
            self.assertTrue(text.isascii())
            self.assertTrue(text.endswith("\n"))

    def test_non_ascii_header_is_escaped(self):
        text = render_ascii(result(prepare("β-シヌクレイン", SYNUCLEIN.residues)))
        self.assertTrue(text.isascii())
        self.assertIn("\\u03b2", text)


class TestContent(unittest.TestCase):
    def test_header_has_disclaimer_legend_and_summary(self):
        lines = render_ascii(result(SYNUCLEIN)).splitlines()
        header = [line for line in lines if line.startswith("#")]
        self.assertEqual(header[0], "# ROI: heuristic, sequence-only; not a functional/binding site")
        self.assertIn("# CHG9  + net >= +3   - net <= -3   . otherwise   ? n/a", header)
        self.assertIn("# HYD9  # >= 1.5   + 0 to 1.5   . < 0   ? n/a", header)
        self.assertIn("# sp|P37840|SYUA_HUMAN len=140 ROIs=4 coverage=40%", header)

    def test_tracks_match_the_design_example(self):
        # Block 61-120 of alpha-synuclein as printed in the design document.
        lines = render_ascii(result(SYNUCLEIN)).splitlines()
        k = lines.index("SEQ     EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP")
        self.assertEqual(lines[k - 2], "        61        71        81        91        101       111")
        self.assertEqual(lines[k + 1],
                         "CHG9    ..............................................--..----------")
        self.assertEqual(lines[k + 2],
                         "HYD9    .....+#++#+###++++.++..+++++++#++....................+++....")
        self.assertEqual(lines[k + 3],
                         "ROI      [===========]                                [=============")
        self.assertEqual(lines[k + 4], "         ROI-03                                       ROI-04 ->")

    def test_region_frame(self):
        text = render_ascii(result(SYNUCLEIN))
        self.assertIn("ROI-04  107-136  len=30  class=heuristic  mode=simple\n"
                      "  context  ...GKNEEG[APQEGILEDMPVDPDNEAYEMPSEEGYQDY]EPEA\n"
                      "  reason   charge_cluster  ##########  28/30 residues\n"
                      "           rule |net charge, w=9| >= 3 (tool-default)\n"
                      "  detail   net charge -10 (simplified): D/E=10, K/R=0\n"
                      "  edges    +/-4 (window 9); right edge truncated at 136\n", text)

    def test_user_thresholds_and_warnings_are_shown(self):
        text = render_ascii(result("E" * 30, thresholds={"charge_cluster": 4}))
        self.assertIn("thresholds=user-set charge_cluster=4", text)
        self.assertIn("+ net >= +4", text)
        self.assertIn("# WARNING: ROI coverage", text)

    def test_no_regions(self):
        text = render_ascii(result("ACDEFGHIKLMNPQRSTVWY"))
        self.assertIn("ROIs=0", text)
        self.assertNotIn("[", text.split("\n\n", 1)[1])

    def test_close_region_labels_go_to_extra_rows(self):
        from protein.render import _label_rows
        rows = _label_rows([(0, "ROI-01"), (3, "ROI-02"), (8, "ROI-03"), (13, "ROI-04 ->")], 15)
        self.assertEqual(rows, ["ROI-01  ROI-03 ",
                                "   ROI-02      ",
                                "      ROI-04 ->"])  # moved left to stay in the block

    def test_truncated_region_label_is_marked(self):
        from protein.render import _label_rows
        self.assertEqual(_label_rows([(0, "<- ROI-04")], 7), ["<- ROI~"])


class TestArguments(unittest.TestCase):
    def test_bad_width_and_charset(self):
        doc = result(LYSOZYME)
        for width in (19, 0, True, 60.0):
            with self.assertRaises(ValueError):
                render_ascii(doc, width=width)
        with self.assertRaises(ValueError):
            render_ascii(doc, charset="unicode")
