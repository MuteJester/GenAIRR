"""Reproduce metadata projection/export defects without changing the engine."""

import csv
import tempfile
from pathlib import Path

import GenAIRR as ga


def run_checks():
    findings = []

    def repertoire(with_post):
        experiment = (
            ga.Experiment.on("human_igh")
            .recombine()
            .clonal_repertoire(n_clones=1, max_size=1)
        )
        if with_post:
            experiment.sequencing_errors(count=0)
        compiled = experiment.compile()
        result = compiled.run_records(seed=42)
        return result, compiled.refdata

    before, _ = repertoire(False)
    after, refdata = repertoire(True)
    fields = (
        "d_trim_3", "j_trim_5", "d_cigar", "j_cigar",
        "d_identity", "j_identity",
    )
    findings.append({
        "id": "clonal_founder_trims",
        "title": "A no-op post-fork step loses founder trim metadata",
        "expected": "Identical molecules retain identical trims and alignments.",
        "observed": {
            "sequences_equal": before[0]["sequence"] == after[0]["sequence"],
            "without_post": {k: before[0][k] for k in fields},
            "with_noop_post": {k: after[0][k] for k in fields},
        },
        "reproduced": (
            before[0]["sequence"] == after[0]["sequence"]
            and any(before[0][k] != after[0][k] for k in fields)
        ),
        "validator_issues": after.validate_records(refdata).failures,
    })

    compiled = ga.Experiment.on("human_igh").recombine().compile()
    result = compiled.run_records(n=1, seed=42)
    record = result[0]
    segment_offsets = []
    for segment in ("v", "d", "j"):
        query_start = record[f"{segment}_sequence_start"]
        germline_start = record[f"{segment}_germline_start"]
        expected_prefix = (
            (f"{query_start}S" if query_start else "")
            + (f"{germline_start}N" if germline_start else "")
        )
        cigar = record[f"{segment}_cigar"]
        if expected_prefix and not cigar.startswith(expected_prefix):
            segment_offsets.append({
                "segment": segment,
                "cigar": cigar,
                "query_start_0based": query_start,
                "germline_start_0based": germline_start,
                "required_leading_operations": expected_prefix,
            })
    findings.append({
        "id": "airr_cigar_offsets",
        "title": "CIGAR lacks the leading AIRR query/reference offsets",
        "expected": "Nonzero query/reference starts have leading S/N operators.",
        "observed": segment_offsets,
        "reproduced": bool(segment_offsets),
        "validator_issues": result.validate_records(compiled.refdata).failures,
    })

    with tempfile.TemporaryDirectory(prefix="genairr-metadata-audit-") as folder:
        exported = {}
        strict_offset_errors = []
        for strict in (None, False, True):
            path = Path(folder) / f"strict_{strict}.tsv"
            if strict is None:
                result.to_tsv(str(path))
            else:
                result.to_tsv(str(path), airr_strict=strict)
            with path.open(newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            exported["default" if strict is None else str(strict)] = {
                key: row[key]
                for key in ("v_sequence_start", "rev_comp", "productive", "stop_codon")
            }
            if strict is True:
                for segment in ("v", "d", "j"):
                    query_start = record[f"{segment}_sequence_start"]
                    reference_start = record[f"{segment}_germline_start"]
                    prefix = (f"{query_start}S" if query_start else "")
                    prefix += f"{reference_start}N" if reference_start else ""
                    if not row[f"{segment}_cigar"].startswith(prefix):
                        strict_offset_errors.append(segment)
    # Keep the native-record finding visible while distinguishing the AIRR
    # interchange boundary, which now supplies the missing placement offsets.
    findings[-1]["strict_tsv_offset_errors"] = strict_offset_errors
    findings.append({
        "id": "airr_tsv_serialization",
        "title": "Default TSV coordinates and even strict TSV booleans violate AIRR",
        "expected": "AIRR TSV starts are 1-based and booleans are T/F.",
        "observed": exported,
        "reproduced": (
            exported["default"]["v_sequence_start"] != "1"
            or exported["True"]["rev_comp"] not in ("T", "F")
        ),
        "validator_issues": "Record validator does not inspect serialized TSV.",
    })

    result[0]["sequence"] = "AAA"
    result[0]["productive"] = "not a boolean"
    result[0]["junction_start"] = 999999
    corrupted_report = result.validate_records(compiled.refdata)
    stamp_rejected = False
    stamped = None
    try:
        stamped = (
            ga.Experiment.on("human_igh").recombine()
            .with_metadata(sequence="AAA", productive="wrong")
            .run_records(n=1, seed=42, validate_records=True)
        )
    except ValueError:
        stamp_rejected = True
    findings.append({
        "id": "validator_ignores_records",
        "title": "Validation checks a rebuilt record instead of the emitted dictionary",
        "expected": "A malformed emitted record is rejected.",
        "observed": {
            "invalid_record_validation_ok": corrupted_report.ok,
            "stamped_sequence": stamped[0]["sequence"] if stamped is not None else None,
            "stamped_productive": stamped[0]["productive"] if stamped is not None else None,
            "metadata_collision_rejected": stamp_rejected,
            "run_records_validate_true_returned": stamped is not None,
        },
        "reproduced": corrupted_report.ok or not stamp_rejected,
        "validator_issues": corrupted_report.failures,
    })
    return findings
