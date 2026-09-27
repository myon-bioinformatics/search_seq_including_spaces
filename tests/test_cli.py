import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parent.parent
TOOL = ROOT / "sequence_tool.py"
FIXTURES = ROOT / "tests" / "fixtures"
GOLDEN = ROOT / "tests" / "golden"
PYTHON_VERSION = "{}.{}.{}".format(*sys.version_info[:3])
FIXTURE_NAMES = ("P00698_mature", "P37840")


def run(*args, stdin=None, cwd=ROOT):
    """Run the CLI without site-packages; return (exit code, stdout bytes, stderr text)."""
    result = subprocess.run([sys.executable, "-S", str(TOOL), *map(str, args)],
                            input=stdin, capture_output=True, cwd=cwd, timeout=60)
    return result.returncode, result.stdout, result.stderr.decode("utf-8")


class TestGoldenOutputs(unittest.TestCase):
    def test_features_json_matches_golden(self):
        for name in FIXTURE_NAMES:
            code, out, err = run("protein", "features", FIXTURES / f"{name}.fasta")
            self.assertEqual((code, err), (0, ""), name)
            normalized = out.replace(f'"python": "{PYTHON_VERSION}"'.encode(),
                                     b'"python": "<python>"')
            self.assertEqual(normalized, (GOLDEN / f"{name}.features.json").read_bytes(), name)

    def test_profile_tsv_matches_golden(self):
        for name in FIXTURE_NAMES:
            code, out, err = run("protein", "profile", FIXTURES / f"{name}.fasta")
            self.assertEqual((code, err), (0, ""), name)
            self.assertEqual(out, (GOLDEN / f"{name}.profile.tsv").read_bytes(), name)

    def test_same_input_gives_same_bytes(self):
        for args in (("protein", "features"), ("protein", "profile", "--json")):
            first = run(*args, FIXTURES / "P37840.fasta")
            self.assertEqual(first, run(*args, FIXTURES / "P37840.fasta"))
            self.assertEqual(first[0], 0)

    def test_bom_and_plain_fixture_give_the_same_output(self):
        plain = run("protein", "profile", "--jsonl", FIXTURES / "P37840.fasta")
        bom = run("protein", "profile", "--jsonl", FIXTURES / "P37840_bom.fasta")
        self.assertEqual(plain, bom)


class TestFormats(unittest.TestCase):
    TWO_RECORDS = b">one first\nACDEFGHIKL\n>two\nKKKK*\n"

    def test_jsonl_writes_one_document_per_record(self):
        code, out, _ = run("protein", "features", "--jsonl", "-", stdin=self.TWO_RECORDS)
        self.assertEqual(code, 0)
        docs = [json.loads(line) for line in out.decode("utf-8").splitlines()]
        self.assertEqual([d["record"]["id"] for d in docs], ["one", "two"])
        self.assertTrue(docs[1]["record"]["modifications"]["stripped_terminal_stop"])

    def test_json_needs_exactly_one_record(self):
        code, out, err = run("protein", "features", "--json", "-", stdin=self.TWO_RECORDS)
        self.assertEqual((code, out), (2, b""))
        self.assertIn("--jsonl", err)

    def test_profile_json_and_options_reach_the_document(self):
        code, out, _ = run("protein", "profile", "--json", "--window", "5",
                           "--entropy-window", "6", "--ambiguity", "mean-of-candidates",
                           "-", stdin=b">x\nACDEFGHIKLBZ\n")
        self.assertEqual(code, 0)
        doc = json.loads(out)
        self.assertEqual(doc["provenance"]["params"]["window"], 5)
        self.assertEqual(doc["provenance"]["params"]["entropy_window"], 6)
        self.assertEqual(doc["record"]["ambiguity"]["ambiguity_policy"], "mean-of-candidates")
        self.assertEqual(len(doc["residues"]), 12)

    def test_tsv_uses_na_for_missing_values(self):
        code, out, _ = run("protein", "profile", "-", stdin=b">x\nACX\n")
        self.assertEqual(code, 0)
        lines = out.decode("utf-8").splitlines()
        self.assertEqual(lines[0].split("\t")[:4], ["id", "position", "input_position", "residue"])
        self.assertEqual(lines[3].split("\t")[:7], ["x", "3", "3", "X", "ambiguous", "NA", "NA"])

    def test_non_ascii_header_is_written_as_utf8(self):
        code, out, _ = run("protein", "features", "-", stdin=">β-シヌクレイン\nACD\n".encode())
        self.assertEqual(code, 0)
        self.assertEqual(json.loads(out.decode("utf-8"))["record"]["id"], "β-シヌクレイン")

    def test_runs_from_another_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, out, _ = run("fasta", "validate", FIXTURES / "P37840.fasta", cwd=tmp)
        self.assertEqual(code, 0)
        self.assertIn(b"sp|P37840|SYUA_HUMAN\t140\t140\t0\t0\t0\tfalse\t0", out)


class TestFastaValidate(unittest.TestCase):
    def test_reports_every_record_and_fails_on_any_error(self):
        code, out, err = run("fasta", "validate", "-", stdin=b">a\nACXD\n>b\nAC1D\n>c\nKK\n")
        self.assertEqual(code, 2)
        rows = out.decode("utf-8").splitlines()
        self.assertEqual([r.split("\t")[0] for r in rows], ["id", "a", "c"])
        self.assertIn("position 3", err)
        self.assertIn("'b'", err)

    def test_json_report(self):
        code, out, _ = run("fasta", "validate", "--json", "--strip-gaps", "-",
                           stdin=b">a\nAC-D\n>b\nA*C\n")
        self.assertEqual(code, 2)
        report = json.loads(out)
        self.assertFalse(report["valid"])
        self.assertEqual(report["records"][0]["modifications"]["stripped_gaps"], 1)
        self.assertEqual(report["records"][1]["error"],
                         {"reason": "internal stop (use allow_internal_stop=True to keep it)",
                          "position": 2, "symbol": "*"})

    def test_internal_stop_can_be_allowed(self):
        code, out, _ = run("fasta", "validate", "--allow-internal-stop", "-", stdin=b">a\nA*C\n")
        self.assertEqual(code, 0)
        self.assertIn(b"a\t2\t2\t0\t0\t1\tfalse\t0", out)


class TestInputErrors(unittest.TestCase):
    def assert_input_error(self, *args, stdin=None, message=""):
        code, out, err = run(*args, stdin=stdin)
        self.assertEqual((code, out), (2, b""), err)
        self.assertIn(message, err)
        self.assertNotIn("Traceback", err)

    def test_invalid_symbol(self):
        self.assert_input_error("protein", "features", "-", stdin=b">r\nAC1D\n",
                                message="position 3")

    def test_ambiguity_error_policy(self):
        self.assert_input_error("protein", "profile", "--ambiguity", "error", "-",
                                stdin=b">r\nACXD\n", message="symbol 'X'")

    def test_missing_file_and_not_utf8(self):
        self.assert_input_error("protein", "features", ROOT / "no-such-file.fasta",
                                message="cannot read")
        self.assert_input_error("protein", "features", "-", stdin=b">r\nAC\xff\n",
                                message="not UTF-8")

    def test_no_records(self):
        self.assert_input_error("fasta", "validate", "-", stdin=b"\n; comment only\n",
                                message="no FASTA records")

    def test_bad_options(self):
        self.assert_input_error("protein", "profile", "--window", "8", "-", stdin=b"ACD",
                                message="odd")
        self.assert_input_error("protein", "profile", "--window", "x", "-", stdin=b"ACD",
                                message="integer")
        self.assert_input_error("protein", "features", "--ambiguity", "mean", "-",
                                stdin=b"ACD", message="invalid choice")
        self.assert_input_error("protein", "features", "--json", "--jsonl", "-",
                                stdin=b"ACD", message="not allowed")
        self.assert_input_error(message="required")
