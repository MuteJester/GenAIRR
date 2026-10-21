"""Deterministic probes for the coordinate and alignment biology audit.

Run from the repository root with its freshly built Rust extension:

    .venv/bin/python tools/metadata_biology_audit/coordinates.py

These are diagnostic probes, not regression tests that require discrepancies to
remain present. ``reproduced`` becomes false when the corresponding output
matches the independent expectation. Partial end-loss cases are explicitly
classified as semantic concerns, not definite AIRR violations.
The short custom reference sequences
isolate alignment mechanics; they are not full biological receptor references.
"""

from __future__ import annotations

import json
import re
from typing import Any

import GenAIRR as ga
from GenAIRR import _engine as ge


def _vj_reference() -> Any:
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("IGKV1*01", "IGKV1", b"AAACCCTGT", anchor=6)
    cfg.add_j_allele("IGKJ1*01", "IGKJ1", b"TTTGGG", anchor=0)
    return cfg


def _vdj_reference() -> Any:
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("IGHV1*01", "IGHV1", b"AAACCCTGT", anchor=6)
    cfg.add_d_allele("IGHD1*01", "IGHD1", b"AACGTA")
    cfg.add_j_allele("IGHJ1*01", "IGHJ1", b"TGGGGG", anchor=0)
    return cfg


def _experiment(cfg: Any) -> ga.Experiment:
    lengths = {"np1_lengths": [(3, 1.0)]}
    if cfg.chain_type == "vdj":
        lengths["np2_lengths"] = [(3, 1.0)]
    return (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine(**lengths)
        .trim(enabled=False)
    )


def _record(experiment: ga.Experiment, cfg: Any) -> tuple[Any, dict]:
    outcome = experiment.run(n=1, seed=0)[0]
    return outcome, outcome.airr_record(cfg)


def _interval(record: dict, field: str) -> list:
    return [record[f"{field}_start"], record[f"{field}_end"]]


def _cigar_consumption(cigar: str, operators: str) -> int:
    return sum(
        int(count)
        for count, operator in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if operator in operators
    )


def _finding(
    identifier: str,
    title: str,
    expected: dict,
    observed: dict,
    outcome: Any,
    cfg: Any,
    **data: Any,
) -> dict:
    return {
        "id": identifier,
        "title": title,
        "expected": expected,
        "observed": observed,
        "reproduced": observed != expected,
        "validator_issues": outcome.validate_record(cfg),
        "data": data,
    }


def run_checks() -> list[dict]:
    """Return six independently specified coordinate/alignment findings."""
    findings = []

    cfg = _vj_reference()
    _, forward = _record(_experiment(cfg), cfg)
    outcome, reverse = _record(
        _experiment(cfg).random_strand_orientation(prob=1), cfg
    )

    def biological_annotations(record: dict) -> dict:
        return {
            "junction": record["junction"],
            "junction_aa": record["junction_aa"],
            "v_sequence_interval": _interval(record, "v_sequence"),
            "j_sequence_interval": _interval(record, "j_sequence"),
        }

    findings.append(
        _finding(
            "coordinate_reverse_orientation",
            "Reverse orientation changes coding-strand annotations",
            biological_annotations(forward),
            biological_annotations(reverse),
            outcome,
            cfg,
            sequence=reverse["sequence"],
            sequence_alignment=reverse["sequence_alignment"],
            rev_comp=reverse["rev_comp"],
            expectation=(
                "AIRR rev_comp=true annotations use the reverse complement "
                "of sequence, so the coding-strand annotations must agree "
                "with the otherwise identical forward record."
            ),
            standard="https://docs.airr-community.org/en/stable/datarep/rearrangements.html",
            source=["engine_rs/src/airr_record/sequence.rs:18"],
        )
    )

    cfg = _vdj_reference()
    outcome, record = _record(_experiment(cfg).invert_d(prob=1), cfg)
    d_start, d_end = _interval(record, "d_alignment")
    findings.append(
        _finding(
            "coordinate_d_inversion",
            "Inverted D alignment fabricates deletions and mismatches",
            {
                "d_reference_consumption": 6,
                "d_germline_span": 6,
                "d_identity": 1.0,
                "d_deletion_bases": 0,
                "d_germline_alignment": "TACGTT",
            },
            {
                "d_reference_consumption": _cigar_consumption(
                    record["d_cigar"], "MD=X"
                ),
                "d_germline_span": (
                    record["d_germline_end"] - record["d_germline_start"]
                ),
                "d_identity": record["d_identity"],
                "d_deletion_bases": _cigar_consumption(record["d_cigar"], "D"),
                "d_germline_alignment": record["germline_alignment"][d_start:d_end],
            },
            outcome,
            cfg,
            d_cigar=record["d_cigar"],
            d_inverted=record["d_inverted"],
            n_indels=record["n_indels"],
            expectation=(
                "AACGTA reverse-complements to TACGTT. This untrimmed, "
                "unmutated six-base segment must align to its oriented "
                "reference without deletions."
            ),
            source=[
                "engine_rs/src/airr_record/segment_projection.rs:176",
                "engine_rs/src/airr_record/walk/mod.rs:217",
            ],
        )
    )

    cfg = _vj_reference()
    plan = ge.PassPlan()
    for segment in ("V", "J"):
        plan.push_sample_allele(segment, cfg)
    plan.push_assemble("V")
    plan.push_p_addition("V_3", [(3, 1.0)])
    plan.push_generate_np("NP1", [(3, 1.0)])
    plan.push_p_addition("J_5", [(2, 1.0)])
    plan.push_assemble("J")
    outcome = ge.run(plan, 0, refdata=cfg, allow_curatable_refdata=True)
    record = outcome.airr_record(cfg)
    findings.append(
        _finding(
            "coordinate_p_alignment_order",
            "P nucleotides move to the end of the alignment",
            {"ungapped_sequence_alignment": record["sequence"].upper()},
            {
                "ungapped_sequence_alignment": record["sequence_alignment"]
                .replace("-", "")
                .upper()
            },
            outcome,
            cfg,
            sequence=record["sequence"],
            sequence_alignment=record["sequence_alignment"],
            p_v_3_length=record["p_v_3_length"],
            p_j_5_length=record["p_j_5_length"],
            expectation=(
                "An alignment must preserve query-base order. The P bases "
                "must remain between their adjacent coding regions."
            ),
            source=[
                "engine_rs/src/airr_record/walk/mod.rs:115",
                "engine_rs/src/passes/p_addition.rs:293",
            ],
        )
    )

    for side, segment, germline_interval in ((5, "v", [3, 9]), (3, "j", [0, 3])):
        cfg = _vj_reference()
        experiment = _experiment(cfg)
        experiment = (
            experiment.end_loss_5prime(length=3)
            if side == 5
            else experiment.end_loss_3prime(length=3)
        )
        outcome, record = _record(experiment, cfg)
        findings.append(
            _finding(
                f"coordinate_end_loss_{side}",
                f"{side}' read loss is represented as a sequence deletion",
                {
                    "germline_interval": germline_interval,
                    "identity": 1.0,
                    "cigar_deletion_bases": 0,
                },
                {
                    "germline_interval": _interval(record, f"{segment}_germline"),
                    "identity": record[f"{segment}_identity"],
                    "cigar_deletion_bases": _cigar_consumption(
                        record[f"{segment}_cigar"], "D"
                    ),
                },
                outcome,
                cfg,
                segment=segment.upper(),
                cigar=record[f"{segment}_cigar"],
                sequence_alignment=record["sequence_alignment"],
                realized_end_loss=record[f"end_loss_{side}_length"],
                expectation=(
                    "Read/primer loss removes coverage of terminal germline "
                    "bases; the observed unchanged bases retain identity "
                    "1.0 against their covered reference interval."
                ),
                source=[
                    "engine_rs/src/airr_record/segment_projection.rs:176",
                    "engine_rs/src/airr_record/walk/mod.rs:217",
                ],
            )
        )

    cfg = _vj_reference()
    _, baseline = _record(_experiment(cfg), cfg)
    outcome, record = _record(_experiment(cfg).end_loss_5prime(length=999), cfg)
    findings.append(
        _finding(
            "coordinate_empty_pool_provenance",
            "Complete read loss discards the realized end-loss count",
            {"sequence_length": 0, "end_loss_5_length": len(baseline["sequence"])},
            {
                "sequence_length": record["sequence_length"],
                "end_loss_5_length": record["end_loss_5_length"],
            },
            outcome,
            cfg,
            expectation=(
                "A sampled end loss larger than the read is clamped to its "
                "length. Empty observation must retain that realized count."
            ),
            source=["engine_rs/src/airr_record/builder.rs:75"],
        )
    )
    for finding in findings:
        if finding["id"] in ("coordinate_end_loss_5", "coordinate_end_loss_3"):
            finding["classification"] = (
                "alignment convention requiring a downstream contract decision"
            )
            finding["data"]["qualification"] = (
                "GenAIRR intentionally represents missing read ends as terminal "
                "deletions. AIRR allows tool-dependent identity definitions; "
                "this difference from a covered-interval alignment is a semantic "
                "concern, not independently sufficient evidence of invalid AIRR."
            )
    return findings


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2))
