"""Reorder a saved SSIS shortlist with a local NLI judge, or evaluate fixed toy searches."""
import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


class LocalJudge:
    def __init__(self, directory):
        import onnxruntime as ort
        from tokenizers import Tokenizer

        config = json.loads((directory / "config.json").read_text())
        self.labels = {int(k): v.lower() for k, v in config["id2label"].items()}
        if set(self.labels.values()) != {"entailment", "neutral", "contradiction"}:
            raise ValueError("model must have three named NLI labels")
        self.tokenizer = Tokenizer.from_file(str(directory / "tokenizer.json"))
        self.session = ort.InferenceSession(str(directory / "model.onnx"), providers=["CPUExecutionProvider"])
        self.names = {i.name for i in self.session.get_inputs()}

    def score(self, evidence, objective):
        import numpy as np

        pair = self.tokenizer.encode(evidence, objective)
        if len(pair.ids) > 512:
            raise ValueError("evidence and objective exceed the 512-token model window")
        inputs = {"input_ids": [pair.ids], "attention_mask": [pair.attention_mask],
                  "token_type_ids": [pair.type_ids]}
        logits = self.session.run(None, {k: np.asarray(inputs[k], dtype=np.int64) for k in self.names})[0][0]
        weights = np.exp(logits - logits.max())
        probabilities = weights / weights.sum()
        if not np.isfinite(probabilities).all():
            raise ValueError("model returned nonfinite probabilities")
        return {self.labels[i]: float(p) for i, p in enumerate(probabilities)}


def review(report, objective, judge):
    if report.get("schema") != "match-shortlist/1" or not objective.strip():
        raise ValueError("a match --top report and a nonempty objective are required")
    rows = []
    for row in report["results"]:
        count, length = row["match_count"], row["record_length"]
        if type(count) is not int or type(length) is not int or min(count, length) < 1:
            raise ValueError("expected positive integer count and record length")
        density = 1000 * count / length
        evidence = (f"This record contains {count} matches in {length} characters. "
                    f"There are {density:.3f} matches per 1000 characters.")
        scores = judge.score(evidence, objective)
        rows.append({**row, "matches_per_1000": density, "evidence": evidence, "nli": scores})
    # Stable sort retains input priority on ties. Neutral and contradiction stay visible.
    rows.sort(key=lambda r: -r["nli"]["entailment"])
    return {"schema": "ssis-local-review/1", "objective": objective,
            "omitted_matching_pairs": report["omitted_matching_pairs"], "results": rows}


def evaluate(judge):
    fixtures = json.loads((HERE / "cases.json").read_text())
    trials = []
    for case in fixtures:
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "toy.fa"
            text = []
            for record, count, length in case["records"]:
                assert 0 < 2 * count <= length
                text.append(f">{record}\n" + "AT" * count + "C" * (length - 2 * count))
            fasta.write_text("\n".join(text) + "\n")
            output = subprocess.check_output([sys.executable, "-S", str(ROOT / "sequence_tool.py"),
                                              "match", str(fasta), "--left", "A", "--right", "T",
                                              "--top", str(len(text))])
        shortlist = json.loads(output)
        for objective, expected in case["objectives"]:
            result = review(shortlist, objective, judge)
            rows = result["results"]
            for row in rows:
                row.pop("source")  # temporary path varies across runs
            orders = {"nli": rows, "count": shortlist["results"],
                      "density": sorted(shortlist["results"], key=lambda r: -r["match_count"] / r["record_length"]),
                      "length": sorted(shortlist["results"], key=lambda r: r["record_length"])}
            ranks = {k: [r["record_id"] for r in ordered].index(expected) + 1 for k, ordered in orders.items()}
            trials.append({"case": case["name"], "objective": objective, "expected": expected,
                           "rank": ranks, "results": rows})
    metrics = {k: {"first": sum(t["rank"][k] == 1 for t in trials),
                   "mean_reciprocal_rank": math.fsum(1 / t["rank"][k] for t in trials) / len(trials)}
               for k in trials[0]["rank"]}
    return {"objectives": len(trials), "metrics": metrics, "trials": trials}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", required=True, type=Path)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("evaluate")
    run = sub.add_parser("review")
    run.add_argument("shortlist", type=Path)
    run.add_argument("--objective", required=True)
    args = parser.parse_args()
    judge = LocalJudge(args.model_dir)
    result = (evaluate(judge) if args.command == "evaluate"
              else review(json.loads(args.shortlist.read_text()), args.objective, judge))
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
