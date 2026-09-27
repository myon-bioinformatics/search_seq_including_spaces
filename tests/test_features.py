import json
import math
import pathlib
import random
import unittest
from collections import Counter
from fractions import Fraction

from protein import residues as R
from protein.features import (entropy_values, protein_features, protein_profiles,
                              window_values)
from protein.parse import SequenceError, prepare, read_fasta
from protein.provenance import EVIDENCE_CLASSES, FIELD_CLASSES, ROW_KEYS

FIXTURES = pathlib.Path(__file__).parent / "fixtures"
# UniProt P00698 residues 1-18 (signal peptide); precursor = this + mature fixture.
LYSOZYME_SIGNAL_PEPTIDE = "MRSLLILVLCFLPLAALG"


def load(name):
    return prepare(*read_fasta((FIXTURES / name).read_text(encoding="utf-8"))[0])


def feature(doc, name):
    return next(f for f in doc["features"] if f["name"] == name)


def value(seq, name, **kwargs):
    return feature(protein_features(seq, **kwargs), name)["value"]


def column(seq, name, **kwargs):
    return [row[name] for row in protein_profiles(seq, **kwargs)["residues"]]


class TestHandCalculated(unittest.TestCase):
    def test_molecular_weight(self):
        # 2 x Gly residue mass + water = 2 x 57.0519 + 18.01528 = 132.11908
        self.assertEqual(value("GG", "molecular_weight_avg"), 132.1191)
        # Ala + Trp + water = 71.0788 + 186.2132 + 18.01528 = 275.30728
        self.assertEqual(value("AW", "molecular_weight_avg"), 275.3073)

    def test_net_charge(self):
        self.assertEqual(value("DEKR", "net_charge_simplified"), 0)
        self.assertEqual(value("DDK", "net_charge_simplified"), -1)
        self.assertEqual(value("KRH", "net_charge_simplified"), 2)  # His counts 0

    def test_length_and_composition(self):
        self.assertEqual(value("AAC", "length"), 3)
        composition = value("AAC", "composition")
        self.assertEqual(len(composition), 26)  # 20 canonical + U, O + B, Z, J, X
        self.assertEqual({k: v for k, v in composition.items() if v}, {"A": 2, "C": 1})

    def test_fractions(self):
        self.assertEqual(value("FWYA", "fraction_aromatic"), 0.75)
        self.assertEqual(value("CMAA", "fraction_sulfur"), 0.5)
        self.assertEqual(value("PGAA", "fraction_pro_gly"), 0.5)
        self.assertEqual(value("AAA", "fraction_aromatic"), 0.0)

    def test_positions_are_input_positions(self):
        p = prepare("h", "AC-CGPW", strip_gaps=True)
        self.assertEqual(value(p, "cys_positions"), [2, 4])
        self.assertEqual(value(p, "gly_positions"), [5])
        self.assertEqual(value(p, "pro_positions"), [6])
        self.assertEqual(value(p, "aromatic_positions"), [7])

    def test_window_values_by_hand(self):
        # K K D A with window 3: position 2 = (K, K, D), position 3 = (K, D, A).
        self.assertEqual(column("KKDA", "charge_window", window=3), [None, 1.0, 0.0, None])
        # (-3.9 - 3.9 - 3.5) / 3 = -3.7667 and (-3.9 - 3.5 + 1.8) / 3 = -1.8667
        self.assertEqual(column("KKDA", "hydropathy_window", window=3),
                         [None, -3.7667, -1.8667, None])
        self.assertEqual(column("KKDA", "hydropathy_window", window=1),
                         [-3.9, -3.9, -3.5, 1.8])

    def test_entropy_by_hand(self):
        # Window 4 at 0-based i covers i-2 .. i+1: only i = 2 fits in 4 residues.
        self.assertEqual(column("ACDE", "entropy_window", entropy_window=4),
                         [None, None, 2.0, None])
        # Written as 0.0, never -0.0 (-(1 * log2 1) is -0.0 before rounding).
        self.assertEqual(json.dumps(column("AAAA", "entropy_window", entropy_window=4)),
                         "[null, null, 0.0, null]")
        self.assertEqual(column("AACC", "entropy_window", entropy_window=4),
                         [None, None, 1.0, None])

    def test_edges_are_null(self):
        rows = protein_profiles("A" * 20)["residues"]
        filled = [r["position"] for r in rows if r["hydropathy_window"] is not None]
        self.assertEqual(filled, list(range(5, 17)))  # window 9: 4 residues each end
        self.assertEqual(protein_profiles("A" * 5)["residues"][2]["hydropathy_window"], None)


class TestFixtureValues(unittest.TestCase):
    def test_net_charge(self):
        self.assertEqual(value(load("P00698_mature.fasta"), "net_charge_simplified"), 8)
        self.assertEqual(value(load("P37840.fasta"), "net_charge_simplified"), -9)

    def test_molecular_weight_matches_uniprot(self):
        # UniProt sequence masses (average, Da): P37840 14,460; P00698 precursor 16,239.
        self.assertEqual(round(value(load("P37840.fasta"), "molecular_weight_avg")), 14460)
        precursor = LYSOZYME_SIGNAL_PEPTIDE + load("P00698_mature.fasta").residues
        self.assertEqual(len(precursor), 147)
        self.assertEqual(round(value(precursor, "molecular_weight_avg")), 16239)

    def test_mature_lysozyme_mass_regression(self):
        self.assertEqual(value(load("P00698_mature.fasta"), "molecular_weight_avg"), 14313.1393)


def naive_window(values, window, scale):
    """Direct per-position recomputation with Fractions (reference)."""
    half, out = window // 2, []
    for i in range(len(values)):
        if i < half or i + half >= len(values):
            out.append(None)
            continue
        used = [Fraction(v) for v in values[i - half:i + half + 1] if v is not None]
        if len(used) < (window + 1) // 2:
            out.append(None)
            continue
        mean = sum(used, Fraction(0)) / len(used)
        result = round(float(mean * (window if scale else 1)), 4)
        out.append(0.0 if result == 0 else result)
    return out


def naive_entropy(residues, classes, window):
    out = []
    for i in range(len(residues)):
        start = i - window // 2
        if start < 0 or start + window > len(residues):
            out.append(None)
            continue
        known = [c for c, cls in zip(residues[start:start + window], classes[start:start + window])
                 if cls in ("canonical", "rare-but-defined")]
        if len(known) < (window + 1) // 2:
            out.append(None)
            continue
        n = len(known)
        h = math.fsum(-(c / n) * math.log2(c / n) for c in Counter(known).values())
        out.append(round(h, 4) + 0.0)
    return out


class TestRunningSumsAgainstNaive(unittest.TestCase):
    def random_sequence(self, rng, n):
        alphabet = "ACDEFGHIKLMNPQRSTVWY" * 4 + "BZJXUO"
        return "".join(rng.choice(alphabet) for _ in range(n))

    def test_profiles_match_direct_recomputation(self):
        rng = random.Random(20260926)
        for trial in range(30):
            seq = self.random_sequence(rng, rng.randint(1, 120))
            p = prepare("r", seq)
            for policy in ("exclude", "mean-of-candidates"):
                for window in (1, 3, 9, 15):
                    doc = protein_profiles(p, window=window, ambiguity=policy)
                    rows = doc["residues"]
                    hyd = [None if r["hydropathy"] is None else self._exact_hydropathy(r, policy)
                           for r in rows]
                    charge = [r["charge"] for r in rows]
                    self.assertEqual([r["hydropathy_window"] for r in rows],
                                     naive_window(hyd, window, scale=False), (trial, seq, window))
                    self.assertEqual([r["charge_window"] for r in rows],
                                     naive_window(charge, window, scale=True), (trial, seq, window))
            for ew in (1, 4, 12):
                self.assertEqual(entropy_values(p.residues, p.classes, ew),
                                 naive_entropy(p.residues, p.classes, ew), (trial, seq, ew))

    @staticmethod
    def _exact_hydropathy(row, policy):
        code = row["residue"]
        if row["class"] == "ambiguous":
            candidates = R.AMBIGUOUS[code]
            return sum((Fraction(str(R.RESIDUES[c].kd)) for c in candidates),
                       Fraction(0)) / len(candidates)
        return Fraction(str(R.RESIDUES[code].kd))

    def test_window_values_accepts_ints_and_fractions(self):
        value_of = {"a": 1, "n": None, "t": Fraction(1, 3)}
        self.assertEqual(window_values("ant", value_of, 3), [None, 0.6667, None])
        self.assertEqual(window_values("nna", value_of, 3), [None, None, None])
        self.assertEqual(window_values("", value_of, 3), [])
        # A tiny negative mean rounds to 0.0, written without a minus sign.
        self.assertEqual(json.dumps(window_values("m", {"m": Fraction(-1, 100000)}, 1)), "[0.0]")


class TestAmbiguityPolicy(unittest.TestCase):
    def test_exclude_is_the_default_and_records_positions(self):
        doc = protein_features(prepare("h", "AC D-XBK", strip_gaps=True))
        self.assertEqual(doc["record"]["ambiguity"],
                         {"ambiguity_policy": "exclude", "ambiguous_count": 2,
                          "excluded_positions": [5, 6]})
        self.assertIsNone(feature(doc, "molecular_weight_avg")["value"])
        self.assertEqual(doc["provenance"]["params"]["ambiguity"], "exclude")

    def test_exclude_nulls_ambiguous_values_and_skips_them_in_windows(self):
        rows = protein_profiles("AAAAXAAAA")["residues"]
        self.assertIsNone(rows[4]["hydropathy"])
        self.assertIsNone(rows[4]["charge"])
        self.assertEqual(rows[4]["hydropathy_window"], 1.8)  # mean of the 8 Ala

    def test_window_needs_at_least_half_usable_residues(self):
        # Window 9 needs 5 usable residues.
        self.assertEqual(column("AAAAAXXXX", "hydropathy_window")[4], 1.8)
        self.assertIsNone(column("AAAAXXXXX", "hydropathy_window")[4])
        self.assertIsNone(column("AAAAUUUUU", "hydropathy_window")[4])  # U has no KD value

    def test_mean_of_candidates(self):
        rows = protein_profiles("BZJX", window=1, ambiguity="mean-of-candidates")["residues"]
        self.assertEqual([r["hydropathy"] for r in rows], [-3.5, -3.5, 4.15, -0.49])
        self.assertEqual([r["charge"] for r in rows], [None] * 4)  # never a candidate mean
        doc = protein_features("GB", ambiguity="mean-of-candidates")
        # Gly + mean(Asp, Asn) + water = 57.0519 + 114.5962 + 18.01528
        self.assertEqual(feature(doc, "molecular_weight_avg")["value"], 189.6634)
        self.assertEqual(doc["record"]["ambiguity"]["excluded_positions"], [])
        self.assertEqual(value("XA", "fraction_aromatic", ambiguity="mean-of-candidates"), 0.075)
        self.assertEqual(value("XA", "fraction_aromatic"), 0.0)

    def test_net_charge_range_is_policy_independent(self):
        for policy in ("exclude", "mean-of-candidates"):
            f = feature(protein_features("KBZXJ", ambiguity=policy), "net_charge_simplified")
            # K +1; B [-1, 0]; Z [-1, 0]; X [-1, +1]; J 0
            self.assertEqual(f["range"], [-2, 2])
            self.assertIsNone(f["value"])
        f = feature(protein_features("KJ"), "net_charge_simplified")
        self.assertEqual((f["value"], f["range"]), (1, [1, 1]))

    def test_error_policy(self):
        with self.assertRaises(SequenceError) as ctx:
            protein_profiles(prepare("rec", "AC XD"), ambiguity="error")
        self.assertEqual((ctx.exception.position, ctx.exception.symbol), (3, "X"))
        self.assertEqual(value("ACD", "length", ambiguity="error"), 3)

    def test_unknown_policy(self):
        for call in (protein_features, protein_profiles):
            with self.assertRaises(ValueError):
                call("ACD", ambiguity="mean")


class TestRareResiduesAndStops(unittest.TestCase):
    def test_u_and_o(self):
        doc = protein_features("AUOK")
        self.assertEqual(feature(doc, "net_charge_simplified")["unmodeled_positions"], [2, 3])
        self.assertEqual(feature(doc, "net_charge_simplified")["value"], 1)
        # 71.0788 + 150.0388 + 237.3018 + 128.1741 + 18.01528
        self.assertEqual(feature(doc, "molecular_weight_avg")["value"], 604.6088)
        rows = protein_profiles("AUOK", window=1)["residues"]
        self.assertEqual([r["hydropathy"] for r in rows], [1.8, None, None, -3.9])
        self.assertEqual([r["charge"] for r in rows], [0, None, None, 1])
        self.assertEqual(doc["record"]["ambiguity"]["ambiguous_count"], 0)

    def test_internal_stop_markers(self):
        p = prepare("h", "AC*DE*", allow_internal_stop=True)
        doc = protein_features(p)
        self.assertEqual(feature(doc, "length")["value"], 4)
        self.assertIsNone(feature(doc, "molecular_weight_avg")["value"])
        self.assertEqual(doc["record"]["residue_classes"]["stop"], 1)
        self.assertTrue(doc["record"]["modifications"]["stripped_terminal_stop"])
        stop_row = protein_profiles(p, window=1)["residues"][2]
        self.assertEqual((stop_row["residue"], stop_row["class"], stop_row["hydropathy"],
                          stop_row["charge"], stop_row["hydropathy_window"]),
                         ("*", "stop", None, None, None))


class TestArguments(unittest.TestCase):
    def test_bad_windows(self):
        for window in (0, 2, 8, -1, True, 9.0, "9"):
            with self.assertRaises(ValueError, msg=repr(window)):
                protein_profiles("ACDEFGHIK", window=window)
        for window in (0, -3, 2.5):
            with self.assertRaises(ValueError):
                protein_profiles("ACDEFGHIK", entropy_window=window)

    def test_unsupported_options(self):
        for kwargs in ({"scale": "eisenberg"}, {"charge_model": "ph-aware"}, {"edge": "shrink"}):
            with self.assertRaises(ValueError):
                protein_profiles("ACDEFGHIK", **kwargs)

    def test_input_types(self):
        with self.assertRaises(TypeError):
            protein_features(b"ACD")  # type: ignore[arg-type]
        with self.assertRaises(SequenceError):
            protein_features("AC1D")
        self.assertEqual(protein_features("acd")["record"]["length"], 3)


class TestDocument(unittest.TestCase):
    def documents(self):
        for name in ("P00698_mature.fasta", "P37840.fasta"):
            p = load(name)
            yield protein_features(p)
            yield protein_profiles(p)
        yield protein_features("ABZJXUO*", ambiguity="mean-of-candidates")
        yield protein_profiles(prepare("h", "AC*DUX", allow_internal_stop=True))

    def test_top_level_keys(self):
        for doc in self.documents():
            self.assertEqual(doc["schema"], "protein-hints/1")
            self.assertLessEqual({"schema", "record", "provenance", "field_classes",
                                  "limitations"}, set(doc))
            self.assertEqual(set(doc["record"]["ambiguity"]),
                             {"ambiguity_policy", "ambiguous_count", "excluded_positions"})

    def test_every_value_field_has_an_evidence_class(self):
        for doc in self.documents():
            for name, cls in doc["field_classes"].items():
                self.assertEqual(FIELD_CLASSES[name], cls)
                self.assertIn(cls, EVIDENCE_CLASSES)
            for f in doc.get("features", []):
                self.assertIn(f["name"], doc["field_classes"])
                self.assertEqual(f["evidence_class"], FIELD_CLASSES[f["name"]])
            for row in doc.get("residues", []):
                for key in row:
                    self.assertTrue(key in ROW_KEYS or key in doc["field_classes"], key)

    def test_output_is_deterministic_and_json_safe(self):
        for make in (lambda: protein_features(load("P37840.fasta")),
                     lambda: protein_profiles(load("P37840.fasta"))):
            first = json.dumps(make(), sort_keys=True, allow_nan=False)
            self.assertEqual(first, json.dumps(make(), sort_keys=True, allow_nan=False))
            self.assertNotIn("-0.0,", first)

    def test_record_block(self):
        doc = protein_features(load("P37840.fasta"))
        record = doc["record"]
        self.assertEqual(record["id"], "sp|P37840|SYUA_HUMAN")
        self.assertEqual(record["length"], 140)
        self.assertEqual(record["sequence_sha256"],
                         "9a06387610edd099466d614941c42a433f444feb6a8f25375bb0b6ad8ec9c6f4")
        self.assertEqual(record["residue_classes"],
                         {"canonical": 140, "rare-but-defined": 0, "ambiguous": 0, "stop": 0})

    def test_profile_rows(self):
        p = prepare("h", "A-CD", strip_gaps=True)
        rows = protein_profiles(p, window=1)["residues"]
        self.assertEqual([(r["position"], r["input_position"]) for r in rows],
                         [(1, 1), (2, 3), (3, 4)])
        self.assertEqual(set(rows[0]), set(ROW_KEYS) | {"hydropathy", "charge",
                                                        "hydropathy_window", "charge_window",
                                                        "entropy_window"})
