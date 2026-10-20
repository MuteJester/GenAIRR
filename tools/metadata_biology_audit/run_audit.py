#!/usr/bin/env python3
"""Run the October 2026 biological metadata audit reproductions.

These are diagnostic probes, separate from the normal regression suite.
A reproduced discrepancy is reported in JSON; exit status is zero when the probes
ran successfully, not when the metadata is biologically correct.

Usage:
    .venv/bin/python tools/metadata_biology_audit/run_audit.py --output audit.json
"""

import argparse
import importlib
import json
import subprocess
import sys
from pathlib import Path

import GenAIRR
from GenAIRR import _engine


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--section", choices=("all", "productivity", "coordinates", "germline", "integration"),
        default="all",
    )
    args = parser.parse_args()
    sections = (
        ("productivity", "coordinates", "germline", "integration")
        if args.section == "all" else (args.section,)
    )
    checks = []
    for section in sections:
        checks.extend(importlib.import_module(section).run_checks())
    for check in checks:
        check.setdefault("classification", "confirmed correctness defect")
    root = Path(__file__).resolve().parents[2]
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, check=True,
    ).stdout.strip()
    branch = subprocess.run(
        ["git", "branch", "--show-current"], cwd=root,
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    tracked_changes = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"], cwd=root,
        capture_output=True, text=True, check=True,
    ).stdout
    report = {
        "audit_date": "2026-10-20",
        "commit": revision,
        "branch": branch,
        "tracked_source_modified": bool(tracked_changes.strip()),
        "python": sys.version,
        "package": GenAIRR.__file__,
        "extension": _engine.__file__,
        "checks_run": len(checks),
        "discrepancies_reproduced": sum(bool(check["reproduced"]) for check in checks),
        "confirmed_defect_probes_reproduced": sum(
            bool(check["reproduced"])
            and check["classification"] == "confirmed correctness defect"
            for check in checks
        ),
        "checks": checks,
    }
    text = json.dumps(report, indent=2, default=str)
    if args.output:
        args.output.write_text(text + "\n")
        print(
            f'{report["discrepancies_reproduced"]}/{len(checks)} discrepancies reproduced '
            f'({report["confirmed_defect_probes_reproduced"]} confirmed-defect probes); '
            f'{args.output}'
        )
    else:
        print(text)


if __name__ == "__main__":
    main()
