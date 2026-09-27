import json
import unittest

from comparison import ComparisonError, compare_aligned


class TestIdentity(unittest.TestCase):
    def test_identical_dna(self):
        doc = compare_aligned("ACGTACGT", "ACGTACGT", alphabet="dna")
        self.assertEqual(doc["counts"]["matches"], 8)
        self.assertEqual(doc["counts"]["mismatches"], 0)
        self.assertEqual(doc["percent_identity"], 100.0)
        self.assertEqual(doc["coverage_percent"], 100.0)

    def test_one_mismatch(self):
        doc = compare_aligned("ACGTACGT", "ACGTTCGT", alphabet="dna")
        self.assertEqual((doc["counts"]["matches"], doc["counts"]["mismatches"]), (7, 1))
        self.assertEqual(doc["percent_identity"], 87.5)

    def test_all_mismatch(self):
        doc = compare_aligned("AAAA", "CCCC", alphabet="dna")
        self.assertEqual(doc["percent_identity"], 0.0)

    def test_lowercase_and_whitespace_are_normalized(self):
        doc = compare_aligned("ac gt\n", "ACGT", alphabet="dna")
        self.assertEqual(doc["percent_identity"], 100.0)

    def test_region_is_one_based_inclusive(self):
        doc = compare_aligned("AAAACCCC", "AATACCCC", alphabet="dna", start=2, end=4)
        self.assertEqual(doc["region"], {"start": 2, "end": 4, "columns": 3})
        self.assertEqual(doc["percent_identity"], 66.666667)

    def test_protein_exact_identity(self):
        doc = compare_aligned("ACDEFG", "ACNEFG", alphabet="protein")
        self.assertEqual(doc["percent_identity"], 83.333333)


class TestGaps(unittest.TestCase):
    def test_exclude_any_gap_by_default(self):
        doc = compare_aligned("AC-GT", "ACTGT", alphabet="dna")
        self.assertEqual(doc["counts"]["gap_columns"], 1)
        self.assertEqual(doc["counts"]["excluded_positions"], 1)
        self.assertEqual(doc["counts"]["compared_positions"], 4)
        self.assertEqual(doc["percent_identity"], 100.0)
        self.assertEqual(doc["coverage_percent"], 80.0)

    def test_one_sided_gap_can_be_mismatch(self):
        doc = compare_aligned("AC-GT", "ACTGT", alphabet="dna", gap_policy="mismatch")
        self.assertEqual(doc["counts"]["gap_mismatches"], 1)
        self.assertEqual(doc["counts"]["compared_positions"], 5)
        self.assertEqual(doc["counts"]["mismatches"], 1)
        self.assertEqual(doc["percent_identity"], 80.0)
        self.assertEqual(doc["coverage_percent"], 100.0)

    def test_gap_gap_is_always_excluded(self):
        for policy in ("exclude", "mismatch"):
            doc = compare_aligned("AC--GT", "AC--GT", alphabet="dna", gap_policy=policy)
            self.assertEqual(doc["counts"]["gap_columns"], 2)
            self.assertEqual(doc["counts"]["excluded_positions"], 2)
            self.assertEqual(doc["percent_identity"], 100.0)


class TestAmbiguity(unittest.TestCase):
    def test_exclude_is_default(self):
        doc = compare_aligned("ANR", "AGA", alphabet="dna")
        self.assertEqual(doc["counts"]["ambiguous_columns"], 2)
        self.assertEqual(doc["counts"]["compared_positions"], 1)
        self.assertEqual(doc["percent_identity"], 100.0)
        self.assertEqual(doc["coverage_percent"], 33.333333)

    def test_strict_compares_literal_symbols(self):
        doc = compare_aligned("ANN", "ANA", alphabet="dna", ambiguity_policy="strict")
        self.assertEqual(doc["counts"]["matches"], 2)
        self.assertEqual(doc["counts"]["mismatches"], 1)
        self.assertEqual(doc["percent_identity"], 66.666667)
        self.assertIsNone(doc["percent_compatible"])

    def test_compatible_reports_separate_metric(self):
        doc = compare_aligned("RRYN", "AGCA", alphabet="dna", ambiguity_policy="compatible")
        self.assertEqual(doc["counts"]["matches"], 0)
        self.assertEqual(doc["percent_identity"], 0.0)
        self.assertEqual(doc["counts"]["compatible_matches"], 4)
        self.assertEqual(doc["percent_compatible"], 100.0)

    def test_rna_uses_u_in_ambiguity_sets(self):
        doc = compare_aligned("Y", "U", alphabet="rna", ambiguity_policy="compatible")
        self.assertEqual(doc["percent_compatible"], 100.0)

    def test_protein_ambiguity_is_explicit(self):
        doc = compare_aligned("BJZX", "DEIA", alphabet="protein",
                              ambiguity_policy="compatible")
        self.assertEqual(doc["percent_compatible"], 50.0)


class TestErrors(unittest.TestCase):
    def test_unequal_lengths_are_not_zipped_silently(self):
        with self.assertRaisesRegex(ComparisonError, "equal length"):
            compare_aligned("ACGT", "ACG", alphabet="dna")

    def test_empty_is_rejected(self):
        with self.assertRaisesRegex(ComparisonError, "must not be empty"):
            compare_aligned("", "", alphabet="dna")

    def test_invalid_symbol_reports_position(self):
        with self.assertRaisesRegex(ComparisonError, "position 3"):
            compare_aligned("ACQG", "ACAG", alphabet="dna")

    def test_non_ascii_is_rejected_before_case_conversion(self):
        with self.assertRaisesRegex(ComparisonError, "non-ASCII"):
            compare_aligned("AıGT", "AIGT", alphabet="protein")

    def test_bad_region(self):
        for kwargs in ({"start": 0}, {"start": True}, {"start": 3, "end": 2},
                       {"end": 99}):
            with self.assertRaises(ComparisonError, msg=kwargs):
                compare_aligned("ACGT", "ACGT", alphabet="dna", **kwargs)

    def test_bad_policies(self):
        for kwargs in ({"alphabet": "auto"}, {"gap_policy": "ignore"},
                       {"ambiguity_policy": "mean"}):
            with self.assertRaises(ComparisonError, msg=kwargs):
                compare_aligned("ACGT", "ACGT", alphabet=kwargs.pop("alphabet", "dna"),
                                **kwargs)


class TestDocument(unittest.TestCase):
    def test_provenance_and_schema(self):
        doc = compare_aligned("ACGT", "ACGT", alphabet="dna")
        self.assertEqual(doc["schema"], "sequence-comparison/v1")
        self.assertEqual(doc["coordinate_convention"], "1-based-inclusive")
        self.assertEqual(doc["alignment"], "pre-aligned-required")
        self.assertEqual(doc["provenance"]["gap_policy"], "exclude")
        self.assertEqual(doc["provenance"]["ambiguity_policy"], "exclude")
        self.assertEqual(len(doc["provenance"]["sequence_a_sha256"]), 64)

    def test_deterministic_json(self):
        left = compare_aligned("ACGT", "ACGT", alphabet="dna")
        right = compare_aligned("ACGT", "ACGT", alphabet="dna")
        self.assertEqual(json.dumps(left, sort_keys=True), json.dumps(right, sort_keys=True))

    def test_zero_comparable_positions_uses_null_percentages(self):
        doc = compare_aligned("NN--", "NN--", alphabet="dna")
        self.assertEqual(doc["counts"]["compared_positions"], 0)
        self.assertIsNone(doc["percent_identity"])
        self.assertEqual(doc["coverage_percent"], 0.0)
