"""Export the standalone duplex projection without flattening its ancestry."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from ._result_export import _DEFAULT_COLUMN_ORDER, _to_airr_strict


def export_duplex_projection(projection, path, *, format="json", airr_strict=True):
    """Write the full JSON envelope, or its AIRR record as one TSV row.

    JSON preserves native zero-based coordinates. Strict TSV converts only
    AIRR record coordinates; the JSON ancestry cell keeps its declared native
    coordinate spaces. Unknown scalars become empty cells, never zeroes.
    """
    if not isinstance(projection, dict) or projection.get("model") != "duplex_airr_projection_v1":
        raise ValueError("unsupported duplex projection schema")
    if not isinstance(projection.get("record"), dict) or not isinstance(projection.get("paths"), dict):
        raise ValueError("duplex projection requires record and paths")
    if format not in ("json", "tsv"):
        raise ValueError("format must be 'json' or 'tsv'")
    destination = Path(path)
    if format == "json":
        # Serialize before opening so invalid values cannot truncate a file.
        payload = json.dumps(projection, sort_keys=True, allow_nan=False, indent=2) + "\n"
        destination.write_text(payload, encoding="utf-8")
        return
    source = projection["record"]
    row = _to_airr_strict(source) if airr_strict else dict(source)
    if airr_strict:
        # Canonical D bounds remain forward coordinates; its CIGAR traverses
        # the oriented reference. Asymmetric inversion needs that offset.
        for segment in ("v", "d", "j"):
            aligned = projection["paths"].get(segment, {})
            cigar = aligned.get("cigar")
            q = aligned.get("query_start")
            g = aligned.get("oriented_reference_start")
            if cigar is not None and q is not None and g is not None:
                row[f"{segment}_cigar"] = (f"{q}S" if q else "") + (f"{g}N" if g else "") + cigar
    columns = [key for key in _DEFAULT_COLUMN_ORDER if key in row]
    columns.extend(sorted(set(row).difference(columns)))
    for key, value in row.items():
        if isinstance(value, bool):
            row[key] = "T" if value else "F"
        elif isinstance(value, (dict, list)):
            row[key] = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
