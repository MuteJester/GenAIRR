"""``SimulationResult`` — a list-like wrapper around a batch of AIRR
records produced by :meth:`Experiment.run_records`.

The wrapper behaves like a sequence of dicts (``len`` /
``__getitem__`` / iteration), and adds export helpers so users can
write the batch out as AIRR TSV, FASTA, or a pandas DataFrame
without touching the underlying engine objects.

Records are plain ``dict[str, object]`` — see
:func:`._airr_record.outcome_to_airr_record` for the field set. The
underlying ``Outcome`` objects are kept on the wrapper as
``.outcomes`` for advanced consumers who need access to the trace
or revision history.
"""
from __future__ import annotations

import csv
from typing import Any, Dict, Iterator, List, Optional, Sequence, Union, overload


# Canonical AIRR-style column order used by ``to_csv`` / ``to_tsv``.
# Anchors the header row to a stable shape across runs and keeps
# the output diff-friendly.
# Coordinate fields whose values are 0-based half-open (Python
# convention) by default. ``airr_strict=True`` exports add 1 to each
# *_start field to convert to the AIRR-spec 1-based-inclusive form;
# the *_end fields keep their existing value because Python's
# half-open `end` already equals the 1-based-inclusive `end`.
_AIRR_STRICT_START_FIELDS = (
    "v_sequence_start",
    "d_sequence_start",
    "j_sequence_start",
    "v_alignment_start",
    "d_alignment_start",
    "j_alignment_start",
    "v_germline_start",
    "d_germline_start",
    "j_germline_start",
    "junction_start",
)


def _to_airr_strict(rec: Dict[str, Any]) -> Dict[str, Any]:
    """Return a copy of ``rec`` with every coord ``*_start`` field
    converted from 0-based half-open to 1-based-inclusive. ``*_end``
    fields are unchanged because the Python half-open end already
    equals the 1-based-inclusive end (e.g. ``[0, 5)`` ↔ ``[1, 5]``).

    ``None`` start values are left as-is. Unset / missing keys are
    silently ignored.
    """
    converted = dict(rec)
    for field in _AIRR_STRICT_START_FIELDS:
        val = converted.get(field)
        if isinstance(val, int):
            converted[field] = val + 1
    return converted


_DEFAULT_COLUMN_ORDER = [
    # AIRR metadata
    "sequence_id",
    "sequence",
    "sequence_aa",
    "sequence_alignment",
    "germline_alignment",
    "germline_alignment_d_mask",
    "sequence_length",
    "rev_comp",
    "locus",
    # Calls
    "v_call",
    "v_cigar",
    "v_score",
    "v_identity",
    "v_support",
    "v_sequence_start",
    "v_sequence_end",
    "v_alignment_start",
    "v_alignment_end",
    "v_germline_start",
    "v_germline_end",
    "v_trim_5",
    "v_trim_3",
    "d_call",
    "d_cigar",
    "d_score",
    "d_identity",
    "d_support",
    "d_sequence_start",
    "d_sequence_end",
    "d_alignment_start",
    "d_alignment_end",
    "d_germline_start",
    "d_germline_end",
    "d_trim_5",
    "d_trim_3",
    "j_call",
    "j_cigar",
    "j_score",
    "j_identity",
    "j_support",
    "j_sequence_start",
    "j_sequence_end",
    "j_alignment_start",
    "j_alignment_end",
    "j_germline_start",
    "j_germline_end",
    "j_trim_5",
    "j_trim_3",
    "c_call",
    # Junction
    "junction",
    "junction_aa",
    "junction_start",
    "junction_end",
    "junction_length",
    # NP regions
    "np1",
    "np1_aa",
    "np1_length",
    "np2",
    "np2_aa",
    "np2_length",
    # Functionality
    "productive",
    "vj_in_frame",
    "stop_codon",
    # SHM + corruption (non-AIRR; GenAIRR additions)
    "n_mutations",
    "mutation_rate",
    "n_pcr_errors",
    "n_quality_errors",
    "n_indels",
    "is_contaminant",
]


class SimulationResult:
    """List-like wrapper around a batch of AIRR records.

    ``result[i]`` returns the i-th record dict; ``len(result)`` is
    the number of records; iteration yields records in order.

    The original ``Outcome`` objects (with their full trace +
    revision history) are kept on ``.outcomes`` for advanced
    inspection — most users won't need them.
    """

    __slots__ = ("_records", "_outcomes")

    def __init__(
        self,
        records: Sequence[Dict[str, Any]],
        outcomes: Optional[Sequence] = None,
    ) -> None:
        self._records: List[Dict[str, Any]] = list(records)
        # ``outcomes`` is optional: callers that built records by
        # other means (e.g. round-tripping a TSV) don't have the
        # underlying Outcome objects available.
        self._outcomes: Optional[List] = (
            list(outcomes) if outcomes is not None else None
        )

    @classmethod
    def from_outcomes(
        cls,
        outcomes: Sequence,
        refdata: Any,
        *,
        id_prefix: str = "seq",
    ) -> "SimulationResult":
        """Build a :class:`SimulationResult` from a list of Rust
        ``Outcome`` objects + the refdata they ran against.

        Each record's ``sequence_id`` is set to ``f"{id_prefix}{i}"``
        (e.g. ``seq0``, ``seq1``, …) so AIRR-format consumers see a
        unique per-row identifier out of the box.
        """
        from ._airr_record import outcome_to_airr_record

        records = [
            outcome_to_airr_record(
                o, refdata, sequence_id=f"{id_prefix}{i}"
            )
            for i, o in enumerate(outcomes)
        ]
        return cls(records, outcomes=outcomes)

    # ── list-like access ────────────────────────────────────────────

    @property
    def records(self) -> List[Dict[str, Any]]:
        """The underlying list of record dicts. Mutation through this
        view propagates back into the result."""
        return self._records

    @property
    def outcomes(self) -> Optional[List]:
        """The underlying list of ``Outcome`` objects, or ``None``
        when this :class:`SimulationResult` was built from records
        directly (e.g. loaded from a TSV)."""
        return self._outcomes

    def __len__(self) -> int:
        return len(self._records)

    def __iter__(self) -> Iterator[Dict[str, Any]]:
        return iter(self._records)

    @overload
    def __getitem__(self, key: int) -> Dict[str, Any]: ...
    @overload
    def __getitem__(self, key: slice) -> List[Dict[str, Any]]: ...
    def __getitem__(
        self, key: Union[int, slice]
    ) -> Union[Dict[str, Any], List[Dict[str, Any]]]:
        return self._records[key]

    def __repr__(self) -> str:
        return f"<SimulationResult n={len(self._records)}>"

    # ── exports ─────────────────────────────────────────────────────

    def to_dataframe(self, *, airr_strict: bool = False):
        """Return a :class:`pandas.DataFrame` with one row per record.

        ``airr_strict=True`` converts all 0-based half-open coord
        ``*_start`` fields to the AIRR-spec 1-based-inclusive form
        (``*_end`` fields are unchanged). Useful when handing the
        DataFrame off to AIRR-strict downstream tooling.

        Raises ``ImportError`` if pandas isn't installed (pandas is
        an optional extra: ``pip install GenAIRR[all]``).
        """
        try:
            import pandas as pd
        except ImportError as exc:  # pragma: no cover — depends on env
            raise ImportError(
                "to_dataframe() requires pandas. Install with "
                "`pip install GenAIRR[all]` or `pip install pandas`."
            ) from exc

        if not self._records:
            return pd.DataFrame(columns=_DEFAULT_COLUMN_ORDER)
        records = (
            [_to_airr_strict(r) for r in self._records]
            if airr_strict
            else self._records
        )
        return pd.DataFrame(records, columns=self._column_order())

    def to_tsv(self, path: str, *, airr_strict: bool = False) -> None:
        """Write the records as AIRR-style TSV (tab-separated). The
        header row uses :data:`_DEFAULT_COLUMN_ORDER`.

        ``airr_strict=True`` converts coord ``*_start`` fields to
        1-based-inclusive (AIRR spec).
        """
        self._write_delimited(path, "\t", airr_strict=airr_strict)

    def to_csv(self, path: str, *, airr_strict: bool = False) -> None:
        """Write the records as comma-separated values. Convenience
        alongside :meth:`to_tsv` — most analysis tooling prefers TSV
        for AIRR data.

        ``airr_strict=True`` converts coord ``*_start`` fields to
        1-based-inclusive (AIRR spec).
        """
        self._write_delimited(path, ",", airr_strict=airr_strict)

    def to_fasta(self, path: str, *, prefix: str = "seq") -> None:
        """Write the assembled sequences as FASTA. Each record gets
        a header of the form ``">{prefix}{i}|v_call=...|j_call=..."``.
        """
        with open(path, "w", encoding="utf-8") as fh:
            for i, rec in enumerate(self._records):
                seq = rec.get("sequence", "")
                v_call = rec.get("v_call") or ""
                j_call = rec.get("j_call") or ""
                fh.write(f">{prefix}{i}|v_call={v_call}|j_call={j_call}\n")
                fh.write(f"{seq}\n")

    def to_fastq(
        self,
        path: str,
        *,
        quality: str = "illumina",
        prefix: str = "seq",
        **quality_kwargs,
    ) -> None:
        """Write the assembled sequences as FASTQ (Phase 12.E / G1).

        Each record produces:

        ::

            @{prefix}{i}|v_call=...|j_call=...
            <sequence (uppercase)>
            +
            <Phred+33 quality string>

        Parameters:
            path: output file path.
            quality: name of the quality model — ``"illumina"``
                (smoothed trapezoid, default) or ``"constant"``
                (single Q value across the read).
            prefix: per-read header prefix (default ``"seq"``).
            **quality_kwargs: forwarded to the quality model
                constructor. ``ConstantQualityModel`` accepts
                ``q`` (default 30), ``low_q`` (default 10),
                ``n_q`` (default 2). ``IlluminaQualityModel``
                accepts ``peak_q``, ``start_q``, ``end_q``,
                ``ramp_len``, ``tail_len``, ``low_q``, ``n_q``.

        FASTQ uppercases the sequence bases — GenAIRR's lowercase
        corruption-marker convention is preserved by routing
        lowercase positions to ``low_q`` in the quality string,
        the standard FASTQ way of conveying low-confidence bases.
        """
        from ._qmodel import phred_to_ascii, resolve_quality_model

        model = resolve_quality_model(quality, **quality_kwargs)
        with open(path, "w", encoding="utf-8") as fh:
            for i, rec in enumerate(self._records):
                seq = rec.get("sequence", "")
                v_call = rec.get("v_call") or ""
                j_call = rec.get("j_call") or ""
                q_array = model.quality_array(seq)
                if len(q_array) != len(seq):
                    raise RuntimeError(
                        f"quality model returned {len(q_array)} scores for "
                        f"{len(seq)}-base sequence"
                    )
                q_string = phred_to_ascii(q_array)
                fh.write(f"@{prefix}{i}|v_call={v_call}|j_call={j_call}\n")
                fh.write(f"{seq.upper()}\n")
                fh.write("+\n")
                fh.write(f"{q_string}\n")

    # ── internals ───────────────────────────────────────────────────

    def _column_order(self) -> List[str]:
        """Pick the column order for tabular exports. Starts with
        the canonical default order; appends any extra columns the
        records have (e.g. when callers add custom fields)."""
        seen = set(_DEFAULT_COLUMN_ORDER)
        extras: List[str] = []
        for rec in self._records:
            for key in rec:
                if key not in seen:
                    seen.add(key)
                    extras.append(key)
        return _DEFAULT_COLUMN_ORDER + extras

    def _write_delimited(
        self, path: str, delimiter: str, *, airr_strict: bool = False
    ) -> None:
        columns = self._column_order()
        with open(path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=columns,
                delimiter=delimiter,
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            for rec in self._records:
                source = _to_airr_strict(rec) if airr_strict else rec
                # Replace ``None`` with empty string so CSV columns
                # don't carry literal ``"None"`` strings.
                row = {k: ("" if v is None else v) for k, v in source.items()}
                writer.writerow(row)
