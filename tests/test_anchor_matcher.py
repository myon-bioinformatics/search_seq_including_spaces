import json
import itertools
import random
import unittest

from anchor_matcher import find_matches, normalize_sequence, reverse_complement


class AnchorMatcherTests(unittest.TestCase):
    def test_normalization_coordinates_and_exact_gap(self):
        raw = " \taC\r\nN t G\n"
        normalized = normalize_sequence(raw)
        self.assertEqual(normalized.sequence, "ACNTG")
        self.assertEqual(normalized.source_positions, (3, 4, 7, 9, 11))
        match = find_matches(raw, "a c", "t g", min_gap=1)["matches"][0]
        self.assertEqual(match["span"], [1, 5])
        self.assertEqual(match["source_span"], [3, 11])
        self.assertEqual(match["middle"], "N")
        self.assertEqual(match["middle_source_span"], [7, 7])
        self.assertEqual(json.loads(json.dumps(match)), match)

    def test_overlapping_pairs_no_duplicate_front_back_scan(self):
        matches = find_matches("AAAA", "A", "A", min_gap=0, max_gap=1)["matches"]
        self.assertEqual([(m["span"], m["gap"]) for m in matches],
                         [([1, 2], 0), ([1, 3], 1), ([2, 3], 0), ([2, 4], 1), ([3, 4], 0)])
        self.assertIsNone(matches[0]["middle_span"])
        self.assertIsNone(matches[0]["middle_source_span"])
        self.assertEqual(matches[0]["middle"], "")

    def test_gap_clamped_to_available_sequence(self):
        self.assertEqual(len(find_matches("AC", "A", "C", max_gap=10**9)["matches"]), 1)
        self.assertEqual(find_matches("AC", "A", "C", min_gap=10**9)["matches"], [])
        self.assertEqual(find_matches("", "A", "C")["matches"], [])

    def test_literal_N_is_not_removed_or_a_wildcard(self):
        self.assertEqual(find_matches("ANT", "A", "T", min_gap=1)["matches"][0]["middle"], "N")
        self.assertEqual(find_matches("AGT", "AN", "T")["matches"], [])
        self.assertEqual(len(find_matches("AGT", "AN", "T", iupac=True)["matches"]), 1)
        self.assertEqual(find_matches("ACT", "AR", "T", iupac=True)["matches"], [])
        self.assertEqual(len(find_matches("ART", "AG", "T", iupac=True)["matches"]), 1)

    def test_reverse_complement_ambiguity_sets(self):
        codes = "ACGTRYSWKMBDHVN"
        self.assertEqual(reverse_complement(codes), "NBDHVKMWSRYACGT")
        self.assertEqual(reverse_complement(reverse_complement(codes)), codes)
        sets = {"R":"AG", "Y":"CT", "S":"CG", "W":"AT", "K":"GT", "M":"AC",
                "B":"CGT", "D":"AGT", "H":"ACT", "V":"ACG", "N":"ACGT"}
        for code, bases in sets.items():
            opposite = reverse_complement(code)
            self.assertEqual(set(reverse_complement(bases)), set(sets.get(opposite, opposite)))

    def test_reverse_coordinates_in_original_reference(self):
        raw = " t\nTaC g "
        result = find_matches(raw, "CG", "AA", min_gap=1, strand="-")
        self.assertEqual(result["matches"], [{"strand":"-", "span":[1,5], "source_span":[2,8],
            "left_span":[4,5], "middle_span":[3,3], "right_span":[1,2],
            "middle_source_span":[5,5], "left":"CG", "middle":"T", "right":"AA", "gap":1}])

    def test_palindrome_strands_are_distinct_observations(self):
        matches = find_matches("AT", "A", "T", strand="both")["matches"]
        self.assertEqual([m["strand"] for m in matches], ["+", "-"])
        self.assertEqual([m["span"] for m in matches], [[1,2], [1,2]])

    def test_literal_protein_letters(self):
        self.assertEqual(find_matches("MQKLV", "MQ", "LV", min_gap=1)["matches"][0]["middle"], "K")

    def test_invalid_inputs_and_no_silent_truncation(self):
        for text in ["A-C", ">id\nAC", "A1C", "ß", "A\0C"]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                normalize_sequence(text)
        with self.assertRaises(TypeError):
            normalize_sequence(None)
        for options in [{"min_gap":True}, {"max_gap":-1}, {"max_gap":1.5},
                        {"min_gap":2,"max_gap":1}, {"max_matches":0},
                        {"iupac":1}, {"strand":"unknown"}]:
            with self.subTest(options=options), self.assertRaises(ValueError):
                find_matches("AC", "A", "C", **options)
        for left, right in [("", "A"), ("A", " \n")]:
            with self.assertRaises(ValueError):
                find_matches("AC", left, right)
        with self.assertRaises(ValueError):
            reverse_complement("AU")
        with self.assertRaises(ValueError):
            find_matches("AU", "A", "U", iupac=True)
        with self.assertRaises(ValueError):
            find_matches("AAAA", "A", "A", max_matches=2)
        self.assertEqual(len(find_matches("AAAA", "A", "A", max_matches=3)["matches"]),3)

    def test_literal_matches_against_independent_bruteforce(self):
        rng = random.Random(314)
        for _ in range(100):
            sequence = "".join(rng.choice("ACGT") for _ in range(12))
            left, right = rng.choice("ACGT"), rng.choice("ACGT")
            expected = []
            for start, stop in itertools.product(range(12), range(12)):
                gap = stop - start - 1
                if 0 <= gap <= 3 and sequence[start] == left and sequence[stop] == right:
                    expected.append(([start+1,stop+1], sequence[start+1:stop]))
            actual = find_matches(sequence, left, right, max_gap=3)["matches"]
            self.assertEqual([(m["span"],m["middle"]) for m in actual],expected)
