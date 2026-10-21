"""Independent basic alignment arithmetic checks across corruption pipelines.

This sweep checks internal consistency, not full biological correctness.
In particular, coherent deletion columns can pass these checks even when
missing read coverage should not have been represented as a deletion.

Usage:
    .venv/bin/python tools/metadata_biology_audit/coordinate_sweep.py
    .venv/bin/python tools/metadata_biology_audit/coordinate_sweep.py --output sweep.json
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re
from typing import Callable

import GenAIRR as ga


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def check_record(record: dict) -> list[str]:
    """Check record slices and arithmetic without using GenAIRR's validator."""
    errors = []
    sequence = (
        _reverse_complement(record["sequence"])
        if record["rev_comp"]
        else record["sequence"]
    )
    alignment = record["sequence_alignment"]
    germline = record["germline_alignment"]
    if alignment.replace("-", "").upper() != sequence.upper():
        errors.append("alignment_not_query")
    if len(alignment) != len(germline):
        errors.append("alignment_lengths")

    for segment in "vdj":
        if record[f"{segment}_sequence_start"] is None:
            continue
        query_start = record[f"{segment}_sequence_start"]
        query_end = record[f"{segment}_sequence_end"]
        aln_start = record[f"{segment}_alignment_start"]
        aln_end = record[f"{segment}_alignment_end"]
        if aln_start is None:
            continue
        observed = sequence[query_start:query_end].upper()
        aligned = alignment[aln_start:aln_end].replace("-", "").upper()
        if observed != aligned:
            errors.append(f"{segment}_slice")

        operators = [
            (int(count), operator)
            for count, operator in re.findall(
                r"(\d+)([MIDNSHP=X])", record[f"{segment}_cigar"]
            )
        ]
        query_consumption = sum(n for n, op in operators if op in "MI=X")
        reference_consumption = sum(n for n, op in operators if op in "MD=X")
        alignment_consumption = sum(n for n, op in operators if op in "MID=X")
        if query_consumption != query_end - query_start:
            errors.append(f"{segment}_cigar_query_span")
        reference_span = (
            record[f"{segment}_germline_end"]
            - record[f"{segment}_germline_start"]
        )
        if reference_consumption != reference_span:
            errors.append(f"{segment}_cigar_germline_span")
        if alignment_consumption != aln_end - aln_start:
            errors.append(f"{segment}_cigar_alignment_span")
        if aln_end > aln_start:
            matches = sum(
                query_base == reference_base and query_base not in "-N"
                for query_base, reference_base in zip(
                    alignment[aln_start:aln_end].upper(),
                    germline[aln_start:aln_end].upper(),
                )
            )
            expected_identity = matches / (aln_end - aln_start)
            identity = record[f"{segment}_identity"]
            if identity is None or abs(expected_identity - identity) > 1e-12:
                errors.append(f"{segment}_identity")

    if record["read_layout"]:
        if record["r1_sequence"] != record["sequence"][
            record["r1_start"]:record["r1_end"]
        ]:
            errors.append("r1_window")
        if record["r2_sequence"] != _reverse_complement(
            record["sequence"][record["r2_start"]:record["r2_end"]]
        ):
            errors.append("r2_window")
    return errors


def _profiles() -> dict[str, Callable]:
    return {
        "baseline": lambda experiment: experiment,
        "substitutions": lambda experiment: (
            experiment.mutate(count=12)
            .pcr_amplify(count=4)
            .sequencing_errors(count=4)
            .ambiguous_base_calls(count=4)
        ),
        "indels": lambda experiment: experiment.polymerase_indels(count=10),
        "end_loss": lambda experiment: (
            experiment.end_loss_5prime(length=20).end_loss_3prime(length=12)
        ),
        "paired_end": lambda experiment: experiment.paired_end(
            r1_length=75, r2_length=75, insert_size=150
        ),
        "all_forward": lambda experiment: (
            experiment.mutate(count=12)
            .polymerase_indels(count=10)
            .pcr_amplify(count=4)
            .sequencing_errors(count=4)
            .ambiguous_base_calls(count=4)
            .end_loss_5prime(length=20)
            .end_loss_3prime(length=12)
            .paired_end(r1_length=75, r2_length=75, insert_size=150)
        ),
    }


def run_sweep(n: int = 100, seed: int = 0) -> dict:
    """Inspect ``12 * n`` records and retain up to three failures per profile."""
    if n <= 0:
        raise ValueError("n must be positive")
    results = {}
    for chain in ("human_igh", "human_igk"):
        for name, profile in _profiles().items():
            experiment = profile(ga.Experiment.on(chain).recombine())
            records = experiment.run_records(n=n, seed=seed).records
            errors: Counter = Counter()
            examples = []
            records_with_errors = 0
            for index, record in enumerate(records):
                issues = check_record(record)
                errors.update(issues)
                if issues:
                    records_with_errors += 1
                    if len(examples) < 3:
                        examples.append(
                            {"index": index, "errors": issues, "record": record}
                        )
            results[f"{chain}/{name}"] = {
                "n": len(records),
                "records_with_errors": records_with_errors,
                "errors": dict(errors),
                "examples": examples,
            }
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=100, help="records per profile")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path, help="optional detailed JSON output")
    args = parser.parse_args()
    results = run_sweep(n=args.n, seed=args.seed)
    for name, result in results.items():
        print(
            f"{name}: {result['n']} records; "
            f"{result['records_with_errors']} with errors"
        )
    total = sum(result["n"] for result in results.values())
    failed = sum(result["records_with_errors"] for result in results.values())
    print(f"Total: {total} records; {failed} with errors.")
    if args.output:
        args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        print(f"Detailed output: {args.output}")
    return int(failed > 0)


if __name__ == "__main__":
    raise SystemExit(main())
