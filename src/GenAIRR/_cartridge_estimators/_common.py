"""Shared module-level helpers for the cartridge estimators (behavior-
preserving). Moved verbatim from the original single-file module."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

_NP_CANONICAL_BASES = ("A", "C", "G", "T")


def _split_tie_set(raw: str) -> List[str]:
    """Parse a comma-separated AIRR tie-set call into a clean list.

    Empty / whitespace-only entries drop out. The returned list
    preserves insertion order — important for the
    ``ambiguous="truth_first"`` policy which keeps only the first
    entry."""
    return [token.strip() for token in raw.split(",") if token.strip()]


def _allele_names_for_segment(
    buckets: Dict[str, List[Any]],
) -> "set[str]":
    """Flatten a builder's per-gene allele bucket dict into a set
    of canonical allele names. Used by
    :meth:`ReferenceCartridgeBuilder.estimate_allele_usage` to
    decide which AIRR-record allele names are "known" to the
    cartridge under construction."""
    out: "set[str]" = set()
    for allele_list in buckets.values():
        for allele in allele_list:
            name = getattr(allele, "name", None)
            if isinstance(name, str) and name:
                out.add(name)
    return out


def _load_rearrangements(
    source: Any,
) -> "tuple[List[Dict[str, Any]], str]":
    """Normalise a ``rearrangements`` argument into a
    ``(records, source_label)`` tuple.

    Accepted inputs:

    - ``list[dict]`` — used verbatim; source label
      ``"records:N"`` carrying the row count.
    - path-like (``str`` / ``Path``) — parsed via
      :class:`csv.DictReader` with ``delimiter='\\t'`` (AIRR-C TSV
      convention); source label is the path.
    - open text file handle — parsed via :class:`csv.DictReader`;
      source label is ``"file:N"``.

    Raises :class:`TypeError` for unrecognised shapes."""
    import csv

    if isinstance(source, list):
        return list(source), f"records:{len(source)}"
    if hasattr(source, "read"):
        rows = list(csv.DictReader(source, delimiter="\t"))
        return rows, f"file:{len(rows)}"
    if isinstance(source, (str, Path)):
        path = Path(source)
        with open(path, "r", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return rows, str(path)
    raise TypeError(
        f"rearrangements must be a list[dict], path, or open text "
        f"file, got {type(source).__name__}"
    )
