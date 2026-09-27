import json
import pathlib
import unittest

from protein import roi as ROI
from protein.features import protein_profiles
from protein.parse import prepare, read_fasta
from protein.provenance import FIELD_CLASSES
from protein.roi import DEFAULT_THRESHOLDS, LABEL, regions_of_interest

FIXTURES = pathlib.Path(__file__).parent / "fixtures"


def load(name):
    return prepare(*read_fasta((FIXTURES / name).read_text(encoding="utf-8"))[0])


LYSOZYME = load("P00698_mature.fasta")
SYNUCLEIN = load("P37840.fasta")


def detect(seq, profile_kwargs=None, **kwargs):
    return regions_of_interest(protein_profiles(seq, **(profile_kwargs or {})), **kwargs)


def spans(doc):
    return [(r["start"], r["end"]) for r in doc["regions"]]


def summary(doc):
    """(ROI count, coverage in whole percent, longest ROI) as in the design document."""
    lengths = [r["length"] for r in doc["regions"]]
    return (len(lengths), round(doc["summary"]["roi_coverage"] * 100), max(lengths, default=0))


class TestMeasuredDefaults(unittest.TestCase):
    """[measured] ROI threshold sensitivity (design document, 2026-09-26)."""

    def test_lysozyme(self):
        doc = detect(LYSOZYME)
        self.assertEqual(spans(doc), [(75, 80), (111, 117)])
        self.assertEqual(doc["summary"]["roi_coverage"], 0.1008)

    def test_alpha_synuclein(self):
        doc = detect(SYNUCLEIN)
        self.assertEqual(spans(doc), [(14, 20), (47, 52), (62, 74), (107, 136)])
        self.assertEqual(doc["summary"]["roi_coverage"], 0.4)
        self.assertEqual(doc["summary"]["warnings"], [])

    def test_alpha_synuclein_c_terminal_region(self):
        region = detect(SYNUCLEIN)["regions"][3]
        self.assertEqual(region["region_id"], "ROI-04")
        self.assertEqual(region["triggers"], [{
            "feature": "charge_cluster", "residues_triggered": 28, "threshold": ">=3",
            "threshold_source": "tool-default", "rule": "|net charge, w=9| >= 3"}])
        self.assertEqual(region["net_charge_simplified"],
                         {"value": -10, "negative": 10, "positive": 0, "not_counted": 0})
        self.assertEqual(region["edge_truncated"], "end")  # 137-140 cannot be evaluated
        self.assertEqual(region["boundary_uncertainty"], 4)
        self.assertEqual(region["context"], {"before": "GKNEEG",
                                             "region": "APQEGILEDMPVDPDNEAYEMPSEEGYQDY",
                                             "after": "EPEA"})

    def test_sweep_table(self):
        # (changed value, profile kwargs, roi kwargs, lysozyme, alpha-synuclein)
        sweep = [
            ({}, {"thresholds": {"charge_cluster": 2}}, (2, 14, 10), (5, 49, 31)),
            ({}, {}, (2, 10, 7), (4, 40, 30)),
            ({}, {"thresholds": {"charge_cluster": 4}}, (1, 5, 6), (4, 24, 13)),
            ({}, {"thresholds": {"charge_cluster": 5}}, (1, 5, 6), (3, 19, 13)),
            ({}, {"thresholds": {"hydropathy_peak": 1.2}}, (2, 12, 8), (5, 48, 30)),
            ({}, {"thresholds": {"hydropathy_peak": 1.8}}, (2, 10, 7), (3, 36, 30)),
            ({"window": 7}, {}, (2, 11, 7), (5, 37, 13)),
            ({"window": 11}, {}, (1, 5, 7), (3, 36, 30)),
            ({}, {"merge_gap": 0}, (1, 5, 6), (2, 23, 26)),
            ({}, {"merge_gap": 4}, (2, 10, 7), (3, 73, 66)),
            ({}, {"min_len": 3}, (3, 12, 7), (5, 42, 30)),
            ({}, {"min_len": 7}, (1, 5, 7), (3, 36, 30)),
        ]
        for profile_kwargs, kwargs, lysozyme, synuclein in sweep:
            label = (profile_kwargs, kwargs)
            self.assertEqual(summary(detect(LYSOZYME, profile_kwargs, **kwargs)), lysozyme, label)
            self.assertEqual(summary(detect(SYNUCLEIN, profile_kwargs, **kwargs)), synuclein, label)


class TestRuns(unittest.TestCase):
    def runs(self, pattern, gap=2, min_len=1):
        mask = [c in "#" for c in pattern]
        stops = [c == "*" for c in pattern]
        return ROI._runs(mask, stops, gap, min_len)

    def test_gap_merge_and_min_len(self):
        self.assertEqual(self.runs("##..##"), [(0, 5)])        # gap 2 joined
        self.assertEqual(self.runs("##...##"), [(0, 1), (5, 6)])  # gap 3 splits
        self.assertEqual(self.runs("##..##", gap=1), [(0, 1), (4, 5)])
        self.assertEqual(self.runs(".###..#.", min_len=5), [(1, 6)])
        self.assertEqual(self.runs(".###....#", min_len=4), [])
        self.assertEqual(self.runs("........"), [])

    def test_prefix_mean_needs_half_of_the_range(self):
        prefix = ROI._Prefix([1, None, None, None, 3])
        self.assertIsNone(prefix.mean(0, 5))            # 2 of 5 usable
        self.assertEqual(prefix.mean(3, 5), 3)          # 1 of 2 usable
        self.assertIsNone(prefix.mean(-1, 2))           # past the start
        self.assertIsNone(prefix.mean(4, 6))            # past the end

    def test_stop_marker_ends_a_run(self):
        self.assertEqual(self.runs("##*##"), [(0, 1), (3, 4)])
        self.assertEqual(self.runs("##.*.##"), [(0, 1), (5, 6)])


class TestBehaviour(unittest.TestCase):
    def test_nulls_split_regions(self):
        # Windows centred in the X run have < 5 usable residues -> null -> no trigger.
        doc = detect("E" * 20 + "X" * 9 + "E" * 20)
        self.assertEqual(len(doc["regions"]), 2)
        middle = 25  # centre of the X run
        self.assertFalse(any(r["start"] <= middle <= r["end"] for r in doc["regions"]))

    def test_count_triggers_use_literal_known_residue_counts(self):
        # Two real aromatics with three excluded X residues must not be
        # extrapolated to an apparent count of three in a 9-residue window.
        doc = detect("A" * 20 + "AFXXXWAAA" + "A" * 20)
        features = {t["feature"] for r in doc["regions"] for t in r["triggers"]}
        self.assertNotIn("aromatic_cluster", features)

    def test_low_complexity_only_region_marks_start_truncation(self):
        doc = detect("Q" * 20 + "ACDEFGHIKLMNPQRSTVWY",
                     thresholds={"charge_cluster": 99, "charge_transition": 99,
                                 "hydropathy_peak": 99, "hydropathy_gradient": 99,
                                 "cys_cluster": 99, "aromatic_cluster": 99,
                                 "pro_gly_cluster": 99})
        region = doc["regions"][0]
        self.assertEqual(region["triggers"][0]["feature"], "low_complexity")
        self.assertIn(region["edge_truncated"], ("start", "both"))

    def test_charge_transition_needs_opposite_signs(self):
        def features(seq):
            return {t["feature"] for r in detect(seq)["regions"] for t in r["triggers"]}
        # Centre A: left KKKKK (+5), right KAAAA (+1): same sign, difference 4.
        self.assertNotIn("charge_transition", features("A" * 5 + "KKKKK" + "A" + "KAAAA" + "A" * 5))
        # Centre A: left EEEEE (-5), right KKKKK (+5).
        self.assertIn("charge_transition", features("A" * 5 + "EEEEE" + "A" + "KKKKK" + "A" * 5))

    def test_region_away_from_the_start_is_not_truncated_there(self):
        region = detect("A" * 7 + "E" * 20)["regions"][0]
        self.assertEqual((region["start"], region["edge_truncated"]), (6, "end"))

    def test_edges_and_coverage_warning(self):
        doc = detect("E" * 30)
        self.assertEqual(spans(doc), [(5, 26)])
        self.assertEqual(doc["regions"][0]["edge_truncated"], "both")
        self.assertEqual(len(doc["summary"]["warnings"]), 1)
        self.assertIn("50%", doc["summary"]["warnings"][0])

    def test_no_regions(self):
        doc = detect("ACDEFGHIKLMNPQRSTVWY")
        self.assertEqual((doc["regions"], doc["summary"]["roi_count"]), ([], 0))

    def test_stop_marker_splits_regions(self):
        doc = detect(prepare("h", "E" * 15 + "*" + "E" * 15, allow_internal_stop=True))
        self.assertEqual(len(doc["regions"]), 2)
        self.assertNotIn("*", "".join(r["context"]["region"] for r in doc["regions"]))

    def test_input_positions_are_reported(self):
        p = prepare("h", "-" * 3 + "E" * 20, strip_gaps=True)
        region = detect(p)["regions"][0]
        self.assertEqual((region["start"], region["input_start"]), (5, 8))

    def test_threshold_change_changes_result_and_provenance(self):
        default, changed = detect(SYNUCLEIN), detect(SYNUCLEIN, thresholds={"charge_cluster": 4})
        self.assertNotEqual(spans(default), spans(changed))
        self.assertEqual(default["provenance"]["roi"]["thresholds"]["charge_cluster"], 3)
        self.assertEqual(changed["provenance"]["roi"]["thresholds"]["charge_cluster"], 4)
        self.assertEqual(changed["provenance"]["roi"]["threshold_sources"]["charge_cluster"], "user")
        self.assertEqual(default["provenance"]["roi"]["threshold_sources"]["charge_cluster"],
                         "tool-default")

    def test_input_document_is_not_modified(self):
        profiles = protein_profiles(SYNUCLEIN)
        before = json.dumps(profiles, sort_keys=True)
        regions_of_interest(profiles)
        self.assertEqual(json.dumps(profiles, sort_keys=True), before)


class TestDocument(unittest.TestCase):
    def test_regions_are_heuristic_and_explained(self):
        for seq in (LYSOZYME, SYNUCLEIN):
            doc = detect(seq)
            self.assertEqual(doc["field_classes"]["regions"], FIELD_CLASSES["regions"])
            self.assertEqual(FIELD_CLASSES["regions"], "heuristic")
            for region in doc["regions"]:
                self.assertEqual(region["evidence_class"], "heuristic")
                self.assertEqual(region["label"], LABEL)
                self.assertIsNone(region["score"])
                self.assertTrue(region["triggers"])  # never a region without reasons
                self.assertTrue(region["context"]["region"])

    def test_never_called_a_hotspot_or_site_prediction(self):
        text = json.dumps(detect(SYNUCLEIN)).lower()
        self.assertNotIn("hotspot", text)
        self.assertIn("not a functional-site prediction", text)

    def test_provenance_records_the_method(self):
        roi = detect(LYSOZYME)["provenance"]["roi"]
        self.assertEqual(roi["roi_method"], "simple-threshold/1")
        self.assertEqual((roi["merge_gap"], roi["min_len"], roi["half_window"]), (2, 5, 5))
        self.assertEqual(set(roi["thresholds"]), set(DEFAULT_THRESHOLDS))

    def test_deterministic(self):
        self.assertEqual(json.dumps(detect(SYNUCLEIN), sort_keys=True, allow_nan=False),
                         json.dumps(detect(SYNUCLEIN), sort_keys=True, allow_nan=False))


class TestArguments(unittest.TestCase):
    def test_bad_arguments(self):
        profiles = protein_profiles("ACDEFGHIK")
        for kwargs in ({"thresholds": {"unknown": 1}},
                       {"thresholds": {"charge_cluster": "3"}},
                       {"thresholds": {"charge_cluster": float("nan")}},
                       {"thresholds": {"charge_cluster": True}},
                       {"thresholds": {"charge_cluster": -1}},
                       {"thresholds": {"low_complexity": -0.1}},
                       {"merge_gap": -1}, {"min_len": 0}, {"min_len": True},
                       {"half_window": 0}):
            with self.assertRaises(ValueError, msg=kwargs):
                regions_of_interest(profiles, **kwargs)

    def test_needs_a_profile_document(self):
        from protein.features import protein_features
        for doc in ({}, protein_features("ACD"), None):
            with self.assertRaises(ValueError):
                regions_of_interest(doc)  # type: ignore[arg-type]
