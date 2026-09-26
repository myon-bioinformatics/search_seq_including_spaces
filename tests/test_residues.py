import unittest

from protein import residues as R

# Expected values restated from the cited sources (test data, not a second
# property table used by the code).
KYTE_DOOLITTLE_1982 = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
IMGT_DONOR_ONLY = set("RKW")
IMGT_ACCEPTOR_ONLY = set("DE")
IMGT_DONOR_AND_ACCEPTOR = set("NQHSTY")
IMGT_NONE = set("ACGILMFPV")


class TestTableShape(unittest.TestCase):
    def test_22_unique_entries(self):
        self.assertEqual(len(R.RESIDUES), 22)
        self.assertEqual(set(R.RESIDUES), R.CANONICAL | R.RARE_BUT_DEFINED)

    def test_keys_match_codes(self):
        for key, row in R.RESIDUES.items():
            self.assertEqual(len(key), 1)
            self.assertTrue(key.isupper())
            self.assertEqual(row.code, key)

    def test_canonical_is_20(self):
        self.assertEqual(len(R.CANONICAL), 20)
        self.assertEqual(R.RARE_BUT_DEFINED, frozenset("UO"))
        self.assertFalse(R.CANONICAL & R.RARE_BUT_DEFINED)

    def test_table_is_read_only(self):
        with self.assertRaises(TypeError):
            R.RESIDUES["A"] = None  # type: ignore[index]
        with self.assertRaises(AttributeError):
            R.RESIDUES["A"].kd = 0.0  # type: ignore[misc]

    def test_every_field_has_a_source(self):
        for key in ("kd", "mw_residue_avg", "hbond", "charge_simplified", "water_avg"):
            self.assertIn(key, R.SOURCES)
        self.assertTrue(R.TABLE_VERSION)


class TestValues(unittest.TestCase):
    def test_kyte_doolittle(self):
        for code, value in KYTE_DOOLITTLE_1982.items():
            self.assertEqual(R.RESIDUES[code].kd, value, code)

    def test_undefined_values_are_none_for_rare_residues(self):
        for code in "UO":
            row = R.RESIDUES[code]
            self.assertIsNone(row.kd)
            self.assertIsNone(row.charge_simplified)
            self.assertIsNone(row.hbond_donor)
            self.assertIsNone(row.hbond_acceptor)

    def test_simplified_charge(self):
        charged = {c: R.RESIDUES[c].charge_simplified for c in R.CANONICAL}
        self.assertEqual({c for c, q in charged.items() if q == -1}, set("DE"))
        self.assertEqual({c for c, q in charged.items() if q == 1}, set("KR"))
        self.assertEqual(charged["H"], 0)

    def test_aromatic_and_sulfur_flags(self):
        self.assertEqual({c for c, r in R.RESIDUES.items() if r.aromatic}, set("FWY"))
        self.assertEqual({c for c, r in R.RESIDUES.items() if r.sulfur}, set("CM"))
        self.assertFalse(R.RESIDUES["U"].sulfur)  # selenium, not sulfur

    def test_imgt_hbond_classes(self):
        for c in R.CANONICAL:
            row = R.RESIDUES[c]
            got = (row.hbond_donor, row.hbond_acceptor)
            if c in IMGT_DONOR_ONLY:
                self.assertEqual(got, (True, False), c)
            elif c in IMGT_ACCEPTOR_ONLY:
                self.assertEqual(got, (False, True), c)
            elif c in IMGT_DONOR_AND_ACCEPTOR:
                self.assertEqual(got, (True, True), c)
            else:
                self.assertIn(c, IMGT_NONE)
                self.assertEqual(got, (False, False), c)

    def test_masses_positive_and_ordered(self):
        masses = {c: r.mw_residue_avg for c, r in R.RESIDUES.items()}
        self.assertTrue(all(m > 0 for m in masses.values()))
        self.assertEqual(min(masses, key=masses.get), "G")
        self.assertEqual(max(masses, key=masses.get), "O")
        self.assertEqual(masses["I"], masses["L"])

    def test_masses_reproduce_uniprot_sequence_masses(self):
        # UniProt sequence masses (Da, rounded) recorded in tests/fixtures/SOURCE.md.
        asyn = ("MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK"
                "EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP"
                "DNEAYEMPSEEGYQDYEPEA")
        mass = sum(R.RESIDUES[c].mw_residue_avg for c in asyn) + R.WATER_AVG
        self.assertEqual(round(mass), 14460)

    def test_water_is_the_documented_derived_value(self):
        self.assertAlmostEqual(R.WATER_AVG, 2 * 1.00794 + 15.9994, places=5)
        self.assertIn("1.00794", R.SOURCES["water_avg"])


class TestAmbiguity(unittest.TestCase):
    def test_candidate_sets(self):
        self.assertEqual(R.AMBIGUOUS["B"], frozenset("DN"))
        self.assertEqual(R.AMBIGUOUS["Z"], frozenset("EQ"))
        self.assertEqual(R.AMBIGUOUS["J"], frozenset("IL"))
        self.assertEqual(R.AMBIGUOUS["X"], R.CANONICAL)

    def test_candidates_are_canonical_and_rare_are_not_ambiguous(self):
        for code, candidates in R.AMBIGUOUS.items():
            self.assertTrue(candidates <= R.CANONICAL, code)
            self.assertNotIn(code, R.RESIDUES)
        self.assertFalse(set(R.AMBIGUOUS) & R.RARE_BUT_DEFINED)


class TestSymbolClass(unittest.TestCase):
    def test_each_class(self):
        cases = {
            "A": "canonical", "w": "canonical",
            "U": "rare-but-defined", "o": "rare-but-defined",
            "B": "ambiguous", "z": "ambiguous", "J": "ambiguous", "X": "ambiguous",
            "*": "stop", "-": "gap",
            "1": "invalid", ".": "invalid", "Ω": "invalid", " ": "invalid",
        }
        for ch, expected in cases.items():
            self.assertEqual(R.symbol_class(ch), expected, repr(ch))
            self.assertIn(expected, R.SYMBOL_CLASSES)

    def test_rejects_non_single_characters(self):
        for bad in ("", "AC", None, 1):
            with self.assertRaises(ValueError):
                R.symbol_class(bad)  # type: ignore[arg-type]

    def test_non_ascii_is_invalid_even_if_it_upper_cases_to_a_residue(self):
        # Regression: str.upper() maps these onto canonical letters.
        self.assertEqual("ı".upper(), "I")   # dotless i
        self.assertEqual("ſ".upper(), "S")   # long s
        for ch in ("ı", "ſ", "K", "Ａ", "é"):
            self.assertEqual(R.symbol_class(ch), "invalid", repr(ch))

    def test_ascii_canonical_residues_unchanged(self):
        for code in R.CANONICAL:
            self.assertEqual(R.symbol_class(code), "canonical")
            self.assertEqual(R.symbol_class(code.lower()), "canonical")
