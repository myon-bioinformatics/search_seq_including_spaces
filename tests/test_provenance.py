import hashlib
import sys
import unittest

from protein import provenance as P


class TestConstants(unittest.TestCase):
    def test_schema_and_version(self):
        self.assertEqual(P.SCHEMA, "protein-hints/1")
        self.assertRegex(P.TOOL_VERSION, r"^\d+\.\d+\.\d+$")

    def test_five_evidence_classes(self):
        self.assertEqual(P.EVIDENCE_CLASSES, (
            "sequence-derived", "heuristic", "predicted",
            "structure-derived-experimental", "structure-derived-predicted"))

    def test_field_classes_use_known_classes_and_are_read_only(self):
        self.assertTrue(P.FIELD_CLASSES)
        for name, cls in P.FIELD_CLASSES.items():
            self.assertIn(cls, P.EVIDENCE_CLASSES, name)
        with self.assertRaises(TypeError):
            P.FIELD_CLASSES["new"] = "heuristic"  # type: ignore[index]

    def test_row_keys_are_not_value_fields(self):
        self.assertFalse(set(P.ROW_KEYS) & set(P.FIELD_CLASSES))


class TestHelpers(unittest.TestCase):
    def test_sequence_sha256(self):
        self.assertEqual(P.sequence_sha256("ACD"), hashlib.sha256(b"ACD").hexdigest())
        with self.assertRaises(UnicodeEncodeError):
            P.sequence_sha256("AÇD")

    def test_field_classes_subset_is_sorted_and_strict(self):
        got = P.field_classes(["length", "charge", "length"])
        self.assertEqual(list(got), ["charge", "length"])
        self.assertEqual(got["charge"], "sequence-derived")
        with self.assertRaises(KeyError):
            P.field_classes(["not_registered"])

    def test_provenance_block(self):
        params = {"window": 9}
        block = P.provenance(table_version="aa-props/1", params=params)
        self.assertEqual(set(block), {"tool_version", "python", "table_version", "params",
                                      "extras"})
        self.assertEqual(block["python"], "{}.{}.{}".format(*sys.version_info[:3]))
        self.assertEqual(block["table_version"], "aa-props/1")
        params["window"] = 3  # the block keeps its own copy
        self.assertEqual(block["params"], {"window": 9})
