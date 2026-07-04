"""Serialization/export surface for SimulationResult, extracted verbatim as
a mixin (behavior-preserving). SimulationResult inherits these methods; self
resolves via MRO.
"""
from __future__ import annotations

import csv
from typing import Any, Dict, List


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
    # Per-segment indel counters (docs/indel_provenance_audit.md
    # §6.2). NP1/NP2 indels are excluded from these but still
    # count toward `n_indels`.
    "n_v_indels",
    "n_d_indels",
    "n_j_indels",
    # Per-segment SHM mutation counters
    # (docs/mutation_provenance_audit.md). Aggregated by walking
    # `outcome.events()`, filtering to the mutate.{uniform,s5f}
    # passes only, and bucketing each `BaseChanged` by carried
    # segment. NP1+NP2 roll into `n_np_mutations`. Sum equals
    # `n_mutations` by construction.
    "n_v_mutations",
    "n_d_mutations",
    "n_j_mutations",
    "n_np_mutations",
    # V-subregion SHM partition
    # (docs/v_subregion_mutation_counters_audit.md). Six fields
    # that partition `n_v_mutations` by the assigned V allele's
    # IMGT subregion intervals: five canonical labels plus
    # `n_v_unannotated_mutations` for V events that can't be
    # attributed (missing assignment, empty annotations, V-side
    # CDR3 stretch, or indel-inserted V bases). Aggregated in the
    # same `outcome.events()` walk as the per-segment counters,
    # using the same `mutate.{uniform,s5f}` pass-name filter.
    # On bundled human OGRDB cartridges the unannotated bucket is
    # 0 on every record under the canonical pass order.
    "n_fwr1_mutations",
    "n_cdr1_mutations",
    "n_fwr2_mutations",
    "n_cdr2_mutations",
    "n_fwr3_mutations",
    "n_v_unannotated_mutations",
    # Observation-stage length loss (EndLossPass / primer_trim_*).
    # Distinct from recombination-stage v_trim_*/j_trim_*. See
    # docs/primer_trim_end_loss_audit.md §6.1.
    "end_loss_5_length",
    "end_loss_3_length",
    "is_contaminant",
    # D inversion provenance (V(D)J inversion event). True when the
    # simulation committed the D allele in reverse-complement
    # orientation; false for VJ chains, VDJ chains without
    # `Experiment.invert_d(...)`, and inversion decisions that
    # landed on the forward branch. See
    # `docs/d_inversion_design.md` §6.3.
    "d_inverted",
    # Receptor revision provenance (Slice E of the receptor-revision
    # roadmap). `receptor_revision_applied` is True iff the
    # `ReceptorRevisionPass` fired with applied=True; `original_v_call`
    # carries the V allele name the recombine pass originally committed
    # (empty when no revision happened). `v_call` continues to report
    # the post-revision identity. See
    # `docs/receptor_revision_design.md` §7.
    "receptor_revision_applied",
    "original_v_call",
    # Paired-end / read layout (Slice A of the paired-end roadmap).
    # All eight fields default to empty / None / 0 / empty under
    # the legacy single-molecule projection; the Slice B/C projection
    # layer populates them. The order here matches the Rust struct
    # layout in `engine_rs/src/airr_record/record.rs`. See
    # `docs/paired_end_design.md` §10.
    "read_layout",
    "r1_sequence",
    "r2_sequence",
    "r1_start",
    "r1_end",
    "r2_start",
    "r2_end",
    "insert_size",
]


class _ResultExport:
    """Serialization/export mixin for :class:`SimulationResult`.

    Behavior-preserving extraction: every method resolves ``self``
    (``self._records`` etc.) via the MRO of the concrete
    :class:`SimulationResult` that inherits this mixin.
    """

    __slots__ = ()

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
        """Write the assembled sequences as FASTQ.

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

    def to_paired_fastq(
        self,
        r1_path: str,
        r2_path: str,
        *,
        quality: str = "illumina",
        overwrite: bool = False,
        **quality_kwargs,
    ) -> None:
        """Write the per-record paired-end reads as two FASTQ files.

        Each AIRR record contributes one R1 record (to ``r1_path``)
        and one R2 record (to ``r2_path``):

        ::

            R1 file:                R2 file:
            @{sequence_id}/1        @{sequence_id}/2
            <r1_sequence upper>     <r2_sequence upper>
            +                       +
            <Phred+33 quality>      <Phred+33 quality>

        Read names use the AIRR record's own ``sequence_id`` with the
        canonical Illumina-portable ``/1`` / ``/2`` suffix (older
        convention but universally accepted by BWA / STAR /
        samtools / Picard; the seven-field colon-separated full
        Illumina header doesn't have a GenAIRR analogue — no flow
        cell, no lane, no index — and is out of scope here per
        `docs/fastq_export_design.md` §5).

        ``r2_sequence`` is **already** the reverse complement of
        ``sequence[r2_start:r2_end]`` at projection time (the AIRR
        validator's `PairedEndWindowMismatch { side: R2 }` enforces
        the invariant); this writer outputs it verbatim. Applying a
        second RC would corrupt the read.

        Parameters:
            r1_path: output path for the R1 FASTQ file.
            r2_path: output path for the R2 FASTQ file.
            quality: name of the quality model — ``"illumina"``
                (smoothed trapezoid, default) or ``"constant"``
                (single Q value across the read). Same vocabulary
                as :meth:`to_fastq`. The model is consulted
                independently for R1 and R2; both reads get their
                own quality string starting from position 0 (this
                is the correct Illumina-style behaviour — each
                read has its own ramp-up and tail).
            overwrite: when ``False`` (default) the writer raises
                ``FileExistsError`` if either output path already
                exists. Set ``True`` to allow overwriting.
            **quality_kwargs: forwarded to the quality model
                constructor — same surface as :meth:`to_fastq`.

        Raises:
            FileExistsError: when ``overwrite=False`` and either
                output path already exists.
            ValueError: when a record's ``read_layout`` is not
                ``"paired_end"`` (the experiment hasn't run
                ``.paired_end(...)``), or when ``r1_sequence`` /
                ``r2_sequence`` is empty.
            RuntimeError: when the quality model produces a
                quality array whose length disagrees with the
                read sequence (same shape as :meth:`to_fastq`).

        FASTQ uppercases the read bases — GenAIRR's lowercase
        corruption-marker convention is preserved by routing
        lowercase positions to the model's ``low_q`` parameter,
        same as :meth:`to_fastq`.
        """
        import os

        from ._qmodel import phred_to_ascii, resolve_quality_model

        # ── 1. Output-path overwrite guard. ──────────────────
        if not overwrite:
            for label, path in (("r1_path", r1_path), ("r2_path", r2_path)):
                if os.path.exists(path):
                    raise FileExistsError(
                        f"to_paired_fastq: {label}={path!r} already exists "
                        "and overwrite=False. Pass overwrite=True to "
                        "replace it."
                    )

        # ── 2. Resolve the quality model once. ───────────────
        model = resolve_quality_model(quality, **quality_kwargs)

        # ── 3. Walk records, validating layout + writing. ────
        with open(r1_path, "w", encoding="utf-8") as r1_fh, open(
            r2_path, "w", encoding="utf-8"
        ) as r2_fh:
            for i, rec in enumerate(self._records):
                sequence_id = rec.get("sequence_id") or f"seq{i}"
                # 3a. Read-layout guard.
                read_layout = rec.get("read_layout", "")
                if read_layout != "paired_end":
                    raise ValueError(
                        f"to_paired_fastq: record {i} "
                        f"(sequence_id={sequence_id!r}) has "
                        f"read_layout={read_layout!r} — paired-end FASTQ "
                        "export requires read_layout='paired_end'. Run "
                        "Experiment.paired_end(r1_length=…, "
                        "insert_size=…) on the experiment before "
                        "exporting."
                    )
                r1_seq = rec.get("r1_sequence") or ""
                r2_seq = rec.get("r2_sequence") or ""
                # 3b. Empty-window guard. Belt-and-suspenders —
                # the projection kernel shouldn't produce empty
                # windows on a paired-layout record, but a downstream
                # consumer that hand-edited the record dict would
                # otherwise silently emit a zero-length FASTQ read.
                if not r1_seq:
                    raise ValueError(
                        f"to_paired_fastq: record {i} "
                        f"(sequence_id={sequence_id!r}) has empty "
                        "r1_sequence despite read_layout='paired_end'"
                    )
                if not r2_seq:
                    raise ValueError(
                        f"to_paired_fastq: record {i} "
                        f"(sequence_id={sequence_id!r}) has empty "
                        "r2_sequence despite read_layout='paired_end'"
                    )
                # 3c. Quality strings. Each read is scored
                # independently — Illumina-style ramp shape resets
                # per read, which is the correct biological model
                # (R1 and R2 don't share a per-base quality
                # profile).
                q_r1 = model.quality_array(r1_seq)
                if len(q_r1) != len(r1_seq):
                    raise RuntimeError(
                        f"to_paired_fastq: quality model returned "
                        f"{len(q_r1)} scores for {len(r1_seq)}-base R1 "
                        f"sequence (record {i})"
                    )
                q_r2 = model.quality_array(r2_seq)
                if len(q_r2) != len(r2_seq):
                    raise RuntimeError(
                        f"to_paired_fastq: quality model returned "
                        f"{len(q_r2)} scores for {len(r2_seq)}-base R2 "
                        f"sequence (record {i})"
                    )
                # 3d. Write the two 4-line records.
                r1_fh.write(f"@{sequence_id}/1\n")
                r1_fh.write(f"{r1_seq.upper()}\n")
                r1_fh.write("+\n")
                r1_fh.write(f"{phred_to_ascii(q_r1)}\n")
                r2_fh.write(f"@{sequence_id}/2\n")
                r2_fh.write(f"{r2_seq.upper()}\n")
                r2_fh.write("+\n")
                r2_fh.write(f"{phred_to_ascii(q_r2)}\n")

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
