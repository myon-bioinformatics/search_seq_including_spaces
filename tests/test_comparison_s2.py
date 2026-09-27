import unittest

from comparison import ComparisonError
from comparison_s2 import render_ascii_diff, sliding_identity


class TestSlidingIdentity(unittest.TestCase):
    def test_windows_match_hand_calculation(self):
        doc = sliding_identity(
            "AAAACCCC",
            "AATACACC",
            alphabet="dna",
            window=4,
            step=2,
        )
        self.assertEqual(
            [(w["start"], w["end"], w["percent_identity"]) for w in doc["windows"]],
            [(1, 4, 75.0), (3, 6, 50.0), (5, 8, 75.0)],
        )
        self.assertEqual(
            [w["divergence_percent"] for w in doc["windows"]],
            [25.0, 50.0, 25.0],
        )

    def test_coverage_respects_excluded_ambiguity(self):
        doc = sliding_identity(
            "AANNAA",
            "AAGGAA",
            alphabet="dna",
            window=3,
            step=3,
        )
        self.assertEqual(doc["windows"][0]["coverage_percent"], 33.333333)
        self.assertEqual(doc["windows"][1]["coverage_percent"], 66.666667)

    def test_compatible_metric_stays_separate(self):
        doc = sliding_identity(
            "RRYN",
            "AGCA",
            alphabet="dna",
            ambiguity_policy="compatible",
            window=4,
        )
        window = doc["windows"][0]
        self.assertEqual(window["percent_identity"], 0.0)
        self.assertEqual(window["percent_compatible"], 100.0)
        self.assertEqual(window["divergence_percent"], 100.0)

    def test_gap_mismatch_policy(self):
        doc = sliding_identity(
            "AC-GT",
            "ACTGT",
            alphabet="dna",
            gap_policy="mismatch",
            window=5,
        )
        window = doc["windows"][0]
        self.assertEqual(window["counts"]["gap_mismatches"], 1)
        self.assertEqual(window["percent_identity"], 80.0)
        self.assertEqual(window["coverage_percent"], 100.0)

    def test_zero_comparable_window_has_null_identity_and_divergence(self):
        doc = sliding_identity("NN--", "NN--", alphabet="dna", window=4)
        window = doc["windows"][0]
        self.assertIsNone(window["percent_identity"])
        self.assertIsNone(window["divergence_percent"])
        self.assertEqual(window["coverage_percent"], 0.0)

    def test_region_and_step_use_aligned_coordinates(self):
        doc = sliding_identity(
            "AAAACCCC",
            "AAAACCCC",
            alphabet="dna",
            window=2,
            step=2,
            start=3,
            end=8,
        )
        self.assertEqual(
            [(w["start"], w["end"]) for w in doc["windows"]],
            [(3, 4), (5, 6), (7, 8)],
        )

    def test_window_larger_than_region_is_error(self):
        with self.assertRaisesRegex(ComparisonError, "exceeds selected"):
            sliding_identity("ACGT", "ACGT", alphabet="dna", window=5)

    def test_bad_window_or_step_is_error(self):
        for kwargs in ({"window": 0}, {"window": True}, {"window": 2, "step": 0}):
            with self.assertRaises(ComparisonError):
                sliding_identity("ACGT", "ACGT", alphabet="dna", **kwargs)


class TestAsciiDiff(unittest.TestCase):
    def test_markers_keep_exact_compatible_mismatch_and_excluded_distinct(self):
        text = render_ascii_diff(
            "ARN-C",
            "AGC-T",
            alphabet="dna",
            ambiguity_policy="compatible",
            gap_policy="mismatch",
            width=10,
        )
        self.assertIn("A ARN-C", text)
        self.assertIn("B AGC-T", text)
        self.assertIn("|~~ .", text)

    def test_gap_mismatch_marker(self):
        text = render_ascii_diff(
            "AC-GT",
            "ACTGT",
            alphabet="dna",
            gap_policy="mismatch",
            width=10,
        )
        self.assertIn("||^||", text)

    def test_excluded_ambiguity_marker(self):
        text = render_ascii_diff(
            "ANAA",
            "AGAA",
            alphabet="dna",
            width=10,
        )
        self.assertIn("|?||", text)

    def test_width_validation(self):
        with self.assertRaisesRegex(ComparisonError, "width"):
            render_ascii_diff("ACGT", "ACGT", alphabet="dna", width=9)
