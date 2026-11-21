"""sequence_tool.py CLI: protein/compare golden outputs and errors, the sequence
subcommands (match, candidates, analysis, preprocess) and fetch with offline I/O."""
import hashlib, io, json, pathlib, subprocess, sys, tempfile, unittest
from contextlib import redirect_stderr, redirect_stdout; from pathlib import Path; from unittest.mock import Mock, patch

from seqlib.acquisition import Client
from sequence_tool import main

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

    def test_roi_ascii_matches_golden(self):
        # Completion check for PR A-C: the [measured] ROIs, reproduced byte for byte.
        for name in FIXTURE_NAMES:
            code, out, err = run("protein", "roi", "--ascii", FIXTURES / f"{name}.fasta")
            self.assertEqual((code, err), (0, ""), name)
            self.assertEqual(out, (GOLDEN / f"{name}.roi.txt").read_bytes(), name)

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


class TestCompare(unittest.TestCase):
    def test_aligned_identity_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = pathlib.Path(tmp) / "a.fa"
            b = pathlib.Path(tmp) / "b.fa"
            a.write_text(">a\nACGTACGT\n", encoding="utf-8")
            b.write_text(">b\nACGTTCGT\n", encoding="utf-8")
            code, out, err = run("compare", a, b, "--alphabet", "dna")
        self.assertEqual((code, err), (0, ""))
        doc = json.loads(out)
        self.assertEqual(doc["percent_identity"], 87.5)
        self.assertEqual(doc["coverage_percent"], 100.0)
        self.assertEqual(doc["records"], {"a": {"id": "a"}, "b": {"id": "b"}})

    def test_compare_region_and_gap_policy(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = pathlib.Path(tmp) / "a.fa"
            b = pathlib.Path(tmp) / "b.fa"
            a.write_text(">a\nAC-GTA\n", encoding="utf-8")
            b.write_text(">b\nACTGCA\n", encoding="utf-8")
            code, out, _ = run("compare", a, b, "--alphabet", "dna",
                               "--gap-policy", "mismatch", "--start", "2", "--end", "5")
        self.assertEqual(code, 0)
        doc = json.loads(out)
        self.assertEqual(doc["region"], {"start": 2, "end": 5, "columns": 4})
        self.assertEqual(doc["counts"]["gap_mismatches"], 1)
        self.assertEqual(doc["counts"]["compared_positions"], 4)

    def test_compare_rejects_multiple_records_and_unequal_alignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = pathlib.Path(tmp) / "a.fa"
            b = pathlib.Path(tmp) / "b.fa"
            a.write_text(">a\nACGT\n>a2\nACGT\n", encoding="utf-8")
            b.write_text(">b\nACGT\n", encoding="utf-8")
            code, out, err = run("compare", a, b, "--alphabet", "dna")
            self.assertEqual((code, out), (2, b""))
            self.assertIn("exactly one record", err)
            a.write_text(">a\nACG\n", encoding="utf-8")
            code, out, err = run("compare", a, b, "--alphabet", "dna")
        self.assertEqual((code, out), (2, b""))
        self.assertIn("equal length", err)


class TestCompareS2(unittest.TestCase):
    def test_compare_sliding_windows_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = pathlib.Path(tmp) / "a.fa"
            b = pathlib.Path(tmp) / "b.fa"
            a.write_text(">a\nAAAACCCC\n", encoding="utf-8")
            b.write_text(">b\nAATACACC\n", encoding="utf-8")
            code, out, err = run(
                "compare", a, b, "--alphabet", "dna", "--window", "4", "--step", "2"
            )
        self.assertEqual((code, err), (0, ""))
        doc = json.loads(out)
        self.assertEqual(
            [(w["start"], w["end"], w["percent_identity"]) for w in doc["windows"]],
            [(1, 4, 75.0), (3, 6, 50.0), (5, 8, 75.0)],
        )

    def test_compare_ascii_and_double_stdin_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = pathlib.Path(tmp) / "a.fa"
            b = pathlib.Path(tmp) / "b.fa"
            a.write_text(">a\nAC-GT\n", encoding="utf-8")
            b.write_text(">b\nACTGT\n", encoding="utf-8")
            code, out, err = run(
                "compare", a, b, "--alphabet", "dna", "--gap-policy", "mismatch", "--ascii"
            )
        self.assertEqual((code, err), (0, ""))
        text = out.decode("utf-8")
        self.assertIn("records: a=a b=b  alphabet=dna  gap=mismatch  ambiguity=exclude", text)
        self.assertIn("||^||", text)

        code, out, err = run("compare", "-", "-", "--alphabet", "dna", stdin=b"ACGT\n")
        self.assertEqual((code, out), (2, b""))
        self.assertIn("cannot read both inputs from stdin", err)

    def test_compare_rejects_incompatible_s2_options(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = pathlib.Path(tmp) / "a.fa"
            b = pathlib.Path(tmp) / "b.fa"
            a.write_text(">a\nACGT\n", encoding="utf-8")
            b.write_text(">b\nACGT\n", encoding="utf-8")
            cases = (
                (("--ascii", "--window", "3"), "--ascii and --window"),
                (("--ascii", "--step", "2"), "--step is only meaningful with --window"),
                (("--width", "20"), "--width is only meaningful with --ascii"),
                (("--step", "2"), "--step is only meaningful with --window"),
            )
            for extra, message in cases:
                code, out, err = run("compare", a, b, "--alphabet", "dna", *extra)
                self.assertEqual((code, out), (2, b""), extra)
                self.assertIn(message, err)
                self.assertNotIn("Traceback", err)


class TestRoi(unittest.TestCase):
    def test_json_regions_and_threshold_override(self):
        code, out, _ = run("protein", "roi", FIXTURES / "P37840.fasta")
        self.assertEqual(code, 0)
        doc = json.loads(out)
        self.assertEqual([(r["start"], r["end"]) for r in doc["regions"]],
                         [(14, 20), (47, 52), (62, 74), (107, 136)])
        code, out, _ = run("protein", "roi", "--threshold", "charge_cluster=4",
                           "--threshold", "hydropathy_peak=1.8", "--merge-gap", "1",
                           "--min-len", "4", FIXTURES / "P37840.fasta")
        roi = json.loads(out)["provenance"]["roi"]
        self.assertEqual((roi["thresholds"]["charge_cluster"], roi["thresholds"]["hydropathy_peak"],
                          roi["merge_gap"], roi["min_len"]), (4, 1.8, 1, 4))

    def test_ascii_for_several_records(self):
        code, out, _ = run("protein", "roi", "--ascii", "--width", "20", "-",
                           stdin=b">a\n" + b"E" * 30 + b"\n>b\nACDEFGHIKL\n")
        self.assertEqual(code, 0)
        text = out.decode("ascii")
        self.assertEqual(text.count("# ROI: heuristic"), 2)
        self.assertTrue(all(len(line) <= 28 for line in text.splitlines()))

    def test_bad_roi_options(self):
        for args in (("--threshold", "unknown=1"), ("--threshold", "charge_cluster"),
                     ("--threshold", "charge_cluster=inf"), ("--width", "19"),
                     ("--merge-gap", "-1"), ("--min-len", "0")):
            code, out, err = run("protein", "roi", *args, FIXTURES / "P37840.fasta")
            self.assertEqual((code, out), (2, b""), args)
            self.assertNotIn("Traceback", err)


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


class SequenceCommandsCliTests(unittest.TestCase):
    def test_match_jsonl_pack_stdin_and_errors(self):
        code,out,err=run('match','-','--input-format','fasta','--left','AC','--right','GT','--min-gap',1,'--jsonl',stdin=b'>x\nACNGT\n>y\nTTTT\n')
        self.assertEqual((code,err),(0,''))
        docs=[json.loads(line) for line in out.splitlines()]
        self.assertEqual(docs[0]['result']['matches'][0]['middle'],'N')
        self.assertEqual(docs[1]['result']['matches'],[])
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'patterns.json'
            p.write_text(json.dumps({'schema':'sequence-patterns/1','patterns':[{'id':'p','left':'A','right':'T'}]}))
            code,out,err=run('match','-','--input-format','jsonl','--patterns',p,'--sequence-path','/payload/seq',stdin=b'{"payload":{"seq":"AT"}}\n')
            self.assertEqual((code,err),(0,''))
            self.assertEqual(json.loads(out)['pattern_id'],'p')
            code,out,err=run('match','-','--input-format','fasta','--patterns',p,'--left','A',stdin=b'')
            self.assertEqual((code,out),(2,b''))
        for opts in (('--min-gap','-1'),('--max-matches','0')):
            code,out,err=run('match','-','--input-format','fasta','--left','A','--right','T',*opts,stdin=b'')
            self.assertEqual((code,out),(2,b''))
            self.assertNotIn('Traceback',err)

    def test_match_shortlist_counts_and_stable_ties(self):
        args = ('match', '-', '--input-format', 'fasta', '--left', 'A', '--right', 'T', '--top', 2)
        sample = b'>first\nAT\n>second\nATAT\n>third\nAT\n>empty\nCC\n'
        code, out, err = run(*args, stdin=sample)
        self.assertEqual((code, err), (0, ''))
        report = json.loads(out)
        self.assertEqual((report['scanned_pairs'], report['matching_pairs'], report['omitted_matching_pairs']), (4, 3, 1))
        self.assertEqual([(r['record_id'], r['match_count'], r['first_span']) for r in report['results']],
                         [('second', 2, [1, 2]), ('first', 1, [1, 2])])
        self.assertNotIn('review', report['results'][0])
        self.assertEqual(json.loads(run(*args, stdin=b'>none\nCC\n')[1])['results'], [])
        for extra in (('--jsonl',), ('--chunk-size', 3), ('--objective', 'review')):
            self.assertEqual(run(*args, *extra, stdin=sample)[0], 2)

    def test_jev_review_reranks_shortlist_and_retains_uncertain_hits(self):
        from seqlib.jev_adapter import JevClassifier
        from sequence_tool import _build_parser
        payloads = []
        answers = iter([('no', .99), ('yes', .6), None, ('unknown', .95), ('yes', .85), ('yes', .98)])
        def transport(payload, **kwargs):
            payloads.append(payload['state'])
            answer = next(answers)
            if answer is None:
                raise TimeoutError('private transport detail')
            choice, confidence = answer
            scores = {k: float(k == choice) for k in ('yes', 'no', 'unknown')}
            return json.dumps({'answers': {'needs_deeper_analysis': {
                'type': 'choice', 'choice': choice, 'confidence': confidence, 'probabilities': scores}}}).encode()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'private.fa'
            path.write_text(''.join(f'>r{i} private description\nAT\n' for i in range(7)))
            args = _build_parser().parse_args(['match', str(path), '--left', 'A', '--right', 'T',
                                               '--top', '6', '--jev-model', 'test', '--objective', 'review',
                                               '--min-confidence', '0.8'])
            classifier = JevClassifier(model='test', transport=transport)
            with patch('seqlib.sequence_commands.JevClassifier', return_value=classifier):
                output, code = args.handler(args)
        self.assertEqual(code, 0)
        rows = json.loads(output)['results']
        self.assertEqual([r['record_id'] for r in rows], ['r5', 'r4', 'r1', 'r2', 'r3', 'r0'])
        self.assertEqual(rows[3]['review'], {'decision': 'unknown', 'confidence': None, 'reason': 'classifier-failed'})
        self.assertEqual(len(payloads), 6)
        self.assertTrue(all(set(item) == {'kind', 'objective', 'pattern_id', 'match_count', 'record_length'}
                            for item in payloads))

    def test_chunk_jsonl_and_summary_cli(self):
        code,out,err=run('match','-','--input-format','fasta','--left','A','--right','T','--chunk-size','4','--jsonl',stdin=b'>x\nCCCCAT\n')
        self.assertEqual((code,err),(0,''))
        self.assertEqual(json.loads(out.splitlines()[-1])['match_count'],1)
        with tempfile.TemporaryDirectory() as d:
            a,b=Path(d)/'a.fa',Path(d)/'b.fa'
            a.write_text('>a\nAACC');b.write_text('>b\nAATT')
            code,out,err=run('compare',a,b,'--alphabet','dna','--window','2','--step','2','--summary','--top','1')
            self.assertEqual((code,err),(0,''))
            self.assertEqual(json.loads(out)['omitted_items'],1)

    def test_analysis_candidates_and_preprocess_receipt(self):
        for args,data,schema in (
            (('candidates','--mode','pairwise'),b'>x\nNNN\n','sequence-candidates/1'),
            (('nucleotide-features',),b'>x\nACGTN\n','nucleotide-features/1'),
            (('rna-fold',),b'>x\nGGGAAACCC\n','rna-nussinov/1'),
            (('protein-ionization',),b'>x\nINGAR\n','protein-ionization/1'),
            (('preset','--kind','poly-a'),b'>x\nAAAAAA\n','sequence-marker-clues/1'),
            (('translate',),b'>x\nATGGCC\n','nucleotide-translation/1')):
            code,out,err=run(*args,'-','--input-format','fasta',stdin=data)
            self.assertEqual((code,err),(0,''))
            self.assertEqual(json.loads(out)['schema'],schema)
        with tempfile.TemporaryDirectory() as d:
            receipt=Path(d)/'receipt.json'
            args=('preprocess','-','--input-format','fasta','--export-format','mega','--receipt',receipt)
            code,out,err=run(*args,stdin=b'>a\nACGT\n>b\nAC-T\n')
            self.assertEqual((code,err),(0,''))
            self.assertTrue(out.startswith(b'#MEGA'))
            self.assertEqual(json.loads(receipt.read_text())['output_records'],2)
            before=receipt.read_bytes()
            code,out,err=run(*args,stdin=b'>a\nAAAA\n')
            self.assertEqual((code,out),(2,b''))
            self.assertEqual(receipt.read_bytes(),before)


class FetchCliTests(unittest.TestCase):
    def invoke(self, args, payload=b'', error=None):
        transport = Mock(return_value=payload, side_effect=error)
        def client(**kwargs):
            return Client(transport=transport, retries=0, **kwargs)
        buffer = io.BytesIO()
        stdout = io.TextIOWrapper(buffer, encoding='utf-8')
        stderr = io.StringIO()
        with patch('seqlib.sequence_commands.Client', side_effect=client), \
                redirect_stdout(stdout), redirect_stderr(stderr):
            try:
                code = main(['fetch', *args])
            except SystemExit as exc:
                code = exc.code
        stdout.flush()
        return code, buffer.getvalue(), stderr.getvalue(), transport

    def test_structure_cif_pdb_and_default(self):
        for fmt in (None, 'cif', 'pdb'):
            with self.subTest(format=fmt):
                expected = fmt or 'cif'
                payload = b'data_1CRN\n#\n' if expected == 'cif' else b'HEADER    TEST\nEND\n'
                args = ['--operation', 'records', '--database', 'structure',
                        '--ids', '1crn', '--timeout', '2', '--max-bytes', '1000']
                if fmt:
                    args += ['--format', fmt]
                code, out, err, transport = self.invoke(args, payload)
                self.assertEqual((code, err), (0, ''))
                result = json.loads(out)
                url = f'https://files.rcsb.org/download/1CRN.{expected}'
                transport.assert_called_once_with(url, timeout=2, max_bytes=1000)
                self.assertEqual(result['schema'], 'structure-download/1')
                self.assertEqual((result['id'], result['format']), ('1CRN', expected))
                self.assertEqual(result['text'], payload.decode())
                self.assertEqual(result['source'], url)
                self.assertEqual(result['sha256'], hashlib.sha256(payload).hexdigest())
                self.assertTrue(result['mmdb_lookup'].endswith('1CRN'))

    def test_invalid_records_fail_before_transport(self):
        cases = [
            ('structure', []),
            ('structure', ['--ids', '1CRN', '2ABC']),
            ('structure', ['--ids', '../secret']),
            ('structure', ['--ids', '1CRN', '--query', 'x']),
            *[('structure', ['--ids', '1CRN', '--format', fmt])
              for fmt in ('fasta', 'genbank', 'xml')],
            *[(db, ['--ids', 'NM_1.2', '--format', fmt])
              for db in ('nuccore', 'protein') for fmt in ('cif', 'pdb')],
            ('pubmed', ['--ids', '123', '--format', 'fasta']),
        ]
        for db, extra in cases:
            with self.subTest(database=db, args=extra):
                code, out, err, transport = self.invoke(
                    ['--operation', 'records', '--database', db, *extra])
                self.assertEqual((code, out), (2, b''))
                self.assertIn('error:', err)
                self.assertNotIn('Traceback', err)
                transport.assert_not_called()

    def test_sequence_and_pubmed_records_unchanged(self):
        for db in ('nuccore', 'protein'):
            for fmt in (None, 'fasta', 'genbank'):
                with self.subTest(database=db, format=fmt):
                    payload = (b'LOCUS       TEST 4 bp DNA\nACCESSION   TEST\n'
                               b'ORIGIN\n        1 acgt\n//\n'
                               if fmt == 'genbank' else b'>TEST\nACGT\n')
                    args = ['--operation', 'records', '--database', db, '--ids', 'TEST']
                    if fmt:
                        args += ['--format', fmt]
                    code, out, err, transport = self.invoke(args, payload)
                    self.assertEqual((code, err), (0, ''))
                    result = json.loads(out)
                    self.assertEqual(result['schema'], 'ncbi-sequences/1')
                    self.assertEqual(result['records'][0]['sequence'],
                                     'acgt' if fmt == 'genbank' else 'ACGT')
                    rettype = ('gb' if db == 'nuccore' else 'gp') if fmt == 'genbank' else 'fasta'
                    self.assertIn('rettype=' + rettype, transport.call_args.args[0])
        code, out, err, transport = self.invoke(
            ['--operation', 'records', '--database', 'pubmed', '--ids', '123'],
            b'<PubmedArticleSet/>')
        self.assertEqual((code, err), (0, ''))
        self.assertEqual(json.loads(out)['articles'], [])
        self.assertIn('efetch.fcgi?db=pubmed', transport.call_args.args[0])

    def test_structure_search_summary_keep_ncbi_dispatch(self):
        for operation, extra, payload in (
            ('search', ['--query', 'crambin'],
             b'{"esearchresult":{"count":"1","idlist":["123"]}}'),
            ('summary', ['--ids', '123'], b'{"result":{"uids":["123"]}}'),
        ):
            with self.subTest(operation=operation):
                args = ['--operation', operation, '--database', 'structure', *extra]
                code, out, err, transport = self.invoke(args, payload)
                self.assertEqual((code, err), (0, ''))
                self.assertEqual(json.loads(out)['database'], 'structure')
                self.assertIn(f'/e{operation}.fcgi?db=structure', transport.call_args.args[0])
                code, out, err, transport = self.invoke(args + ['--format', 'cif'])
                self.assertEqual((code, out), (2, b''))
                transport.assert_not_called()

    def test_structure_transport_failure_and_size_limit(self):
        args = ['--operation', 'records', '--database', 'structure', '--ids', '1CRN']
        code, out, err, transport = self.invoke(args, error=TimeoutError('private detail'))
        self.assertEqual((code, out), (2, b''))
        self.assertIn('fetch failed (TimeoutError)', err)
        self.assertNotIn('private detail', err)
        self.assertEqual(transport.call_count, 1)
        code, out, err, transport = self.invoke(args + ['--max-bytes', '2'], b'abc')
        self.assertEqual((code, out), (2, b''))
        self.assertIn('oversized', err)
