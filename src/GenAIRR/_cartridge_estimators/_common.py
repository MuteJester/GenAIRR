"""Shared module-level helpers for the cartridge estimators (behavior-
preserving). Moved verbatim from the original single-file module."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

_NP_CANONICAL_BASES = ("A", "C", "G", "T")


def _reject_nonzero_p_provenance(
    records: List[Dict[str, Any]], *, estimator: str,
) -> None:
    """Prevent combined AIRR N/P measurements from becoming an N-only model."""
    columns = ("p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length")
    for row_index, row in enumerate(records):
        for column in columns:
            raw = row.get(column)
            if raw is None or raw == "":
                continue
            try:
                value = float(raw)
            except (TypeError, ValueError):
                # Validation of malformed P fields belongs to the P estimator.
                # This guard recognizes explicit numeric nonzero provenance.
                continue
            if value != 0:
                raise ValueError(
                    f"{estimator}: row {row_index} has nonzero {column}={raw!r}; "
                    "AIRR NP fields combine N and P bases and cannot safely "
                    "train separate N and P models. N-only provenance is "
                    "required; subtracting P counts after sequence "
                    "transformations is not sufficient."
                )


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
