#!/usr/bin/env python3
"""Reproduce the bundled catalogue quarantine registry using GenAIRR only.

This is a structural/annotation audit, not a source curation or biological
certification. It never rewrites archived cartridges or guesses corrections.
"""
import argparse
from collections import Counter
import hashlib
import json
from importlib.resources import files

from GenAIRR import data
from GenAIRR._refdata_resolver import dataconfig_to_refdata


def audit_all():
    entries = {}
    for name in data.list_configs():
        cfg = getattr(data, name)
        ref = dataconfig_to_refdata(cfg)
        issues = ref.validate()
        counts = Counter(i["kind"] for i in issues)
        artifact = files("GenAIRR.data").joinpath("builtin_dataconfigs", name + ".pkl")
        entries[name] = {
            "artifact_sha256": "sha256:" + hashlib.sha256(artifact.read_bytes()).hexdigest(),
            "dataconfig_checksum": cfg.schema_sha256,
            "status": "quarantined" if issues else "structurally_valid_evidence_unresolved",
            "issue_counts": dict(sorted(counts.items())),
            "fatal_issue_count": sum(i["severity"] == "fatal" for i in issues),
            "curatable_issue_count": sum(i["severity"] == "curatable" for i in issues),
            "source_evidence": "authoritative_source_snapshot_not_bundled",
            "coding_boundary_evidence": "unknown",
            "landmark_evidence": "legacy_coordinates_not_independently_verified",
            "correction": "none; source-backed curation required",
        }
    return {
        "schema_version": 1,
        "audit_method": "genairr_native_structural_and_anchor_validation_v1",
        "scope": "Provided catalogue bytes and declared rules; no empirical calibration or genomic-boundary certification",
        "quarantine_policy": "Default strict compilation rejects every reported issue. allow_curatable_refdata() explicitly permits curatable diagnostic catalogues only; fatal structural errors remain blocked.",
        "cartridges": entries,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    with open(args.output, "w") as f:
        json.dump(audit_all(), f, indent=2, sort_keys=True)
        f.write("\n")

if __name__ == "__main__":
    main()
