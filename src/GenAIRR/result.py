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


def _allele_name_or_empty(refdata: Any, segment: str, allele_id: Optional[int]) -> str:
    """Look up an allele name from refdata by id; return ``""`` when
    the id is None or the lookup fails (defensive).
    """
    if allele_id is None:
        return ""
    try:
        if segment == "V":
            return refdata.v_allele(int(allele_id)).name
        if segment == "D":
            return refdata.d_allele(int(allele_id)).name
        if segment == "J":
            return refdata.j_allele(int(allele_id)).name
    except Exception:
        return ""
    return ""


def _inject_truth_columns(outcome: Any, refdata: Any, record: Dict[str, Any]) -> None:
    """Append `truth_v_call` / `truth_d_call` / `truth_j_call`
    columns to ``record`` from the originally-sampled allele ids
    stored in the simulation's `assignments`. Distinct from
    `v_call` / `d_call` / `j_call`, which are evidence-driven and
    can change under heavy SHM.
    """
    sim = outcome.final_simulation()
    record["truth_v_call"] = _allele_name_or_empty(refdata, "V", sim.v_allele_id())
    record["truth_d_call"] = _allele_name_or_empty(refdata, "D", sim.d_allele_id())
    record["truth_j_call"] = _allele_name_or_empty(refdata, "J", sim.j_allele_id())


class ValidationReport:
    """Aggregate report from :meth:`SimulationResult.validate_records`.

    Carries the result of the **public AIRR output correctness** gate
    over a batch — the downstream contract that says "every projected
    record is internally consistent with its outcome."

    Attributes:
      ``count``    — total records checked.
      ``failures`` — list of failing records, each as a dict with
                     ``record_index``, ``sequence_id``, and ``issues``
                     (the list of dicts returned by
                     :meth:`Outcome.validate_record`).
      ``ok``       — True iff every record passed (``failures`` is
                     empty).

    The report is truthy iff ``ok`` is True, so ``assert report``
    works as a one-line CI guard.

    Sibling gate: :meth:`Outcome.check_live_call_cache_parity` —
    internal cache-correctness check on the state that feeds
    projection. If both fail on the same batch, fix parity first
    (a stale cache leaks into projection); see
    :meth:`SimulationResult.validate_records` docstring for the
    full troubleshooting rule.
    """

    __slots__ = ("count", "failures")

    def __init__(self, count: int, failures: List[Dict[str, Any]]) -> None:
        self.count = count
        self.failures: List[Dict[str, Any]] = failures

    @property
    def ok(self) -> bool:
        """True iff every record validated without issues."""
        return not self.failures

    def __bool__(self) -> bool:
        return self.ok

    def __len__(self) -> int:
        """Number of failing records (not total). Use ``.count`` for
        the total."""
        return len(self.failures)

    def __repr__(self) -> str:
        if self.ok:
            return f"<ValidationReport ok=True count={self.count}>"
        return (
            f"<ValidationReport ok=False count={self.count} "
            f"failures={len(self.failures)}>"
        )

    def summary(self) -> Dict[str, int]:
        """Histogram of issue kinds across all failures. Useful for
        ``print(report.summary())`` to see what's wrong at a glance."""
        counts: Dict[str, int] = {}
        for failure in self.failures:
            for issue in failure["issues"]:
                kind = issue.get("kind", "Unknown")
                counts[kind] = counts.get(kind, 0) + 1
        return counts


class FamilyValidationReport:
    """Aggregate report from :meth:`SimulationResult.validate_families`.

    Carries the result of the **clonal-family consistency** gate
    over a batch — the downstream contract that says "every group of
    records sharing a ``clone_id`` agrees on the recombination-time
    truth fields the engine puts on each descendant by construction."

    Attributes:
      ``count``               — total records inspected.
      ``family_count``        — number of distinct ``clone_id`` groups
                                (``0`` when the batch carries no
                                clonal records).
      ``members_per_family``  — mapping ``{clone_id: n_members}``;
                                empty when ``family_count == 0``.
      ``failures``            — list of failing-group dicts. Each
                                carries ``clone_id``,
                                ``record_indices``, ``issue_kind``,
                                and ``values`` (the divergent values
                                observed in the group, sorted).
      ``ok``                  — ``True`` iff ``failures`` is empty.

    The report is truthy iff ``ok`` is True so ``assert report``
    works as a one-line CI guard.

    Sibling gate: :class:`ValidationReport` — the per-record AIRR
    projection postcondition check. Family validation runs *after*
    that gate passes; a family-divergence failure means the records
    are each internally consistent with their own outcome, but they
    disagree across the clonal group on a field that should be
    invariant by construction (heavy-SHM tie-set widening on the
    live-call ``v_call`` etc. is intentionally NOT a divergence —
    only ``truth_*_call`` invariance is enforced).

    The validator is a strict subset of the audit's §6 spec
    (see ``docs/clonal_family_design.md``): pre-SHM junction
    invariance, mutation-distance distribution, parent-trace
    reconstruction, and lineage topology are deliberately deferred
    to later slices.
    """

    __slots__ = ("count", "family_count", "members_per_family", "failures")

    def __init__(
        self,
        count: int,
        family_count: int,
        members_per_family: Dict[Any, int],
        failures: List[Dict[str, Any]],
    ) -> None:
        self.count = count
        self.family_count = family_count
        self.members_per_family: Dict[Any, int] = dict(members_per_family)
        self.failures: List[Dict[str, Any]] = failures

    @property
    def ok(self) -> bool:
        """True iff every family group passed every invariant check."""
        return not self.failures

    def __bool__(self) -> bool:
        return self.ok

    def __len__(self) -> int:
        """Number of failing groups (not total). Use ``.count`` for
        the total record count and ``.family_count`` for the total
        group count."""
        return len(self.failures)

    def __repr__(self) -> str:
        if self.ok:
            return (
                f"<FamilyValidationReport ok=True count={self.count} "
                f"family_count={self.family_count}>"
            )
        return (
            f"<FamilyValidationReport ok=False count={self.count} "
            f"family_count={self.family_count} "
            f"failures={len(self.failures)}>"
        )

    def summary(self) -> Dict[str, int]:
        """Histogram of issue kinds across all failing groups. Useful
        for ``print(report.summary())`` to see what's wrong at a
        glance."""
        counts: Dict[str, int] = {}
        for failure in self.failures:
            kind = failure.get("issue_kind", "Unknown")
            counts[kind] = counts.get(kind, 0) + 1
        return counts


# ── family-layer invariants ──────────────────────────────────────────
#
# (Field, issue_kind) pairs the family validator enforces. Each pair
# names a recombination-time provenance field that should be
# identical across every descendant of a clone *when the field is
# present*. We deliberately skip the field for a group when it is
# absent from every record in the group (e.g. ``expose_provenance``
# was not enabled), and require it to be present on every record once
# any record in the group carries it.
_FAMILY_TRUTH_INVARIANTS: tuple = (
    ("truth_v_call", "TruthVCallDiverges"),
    ("truth_d_call", "TruthDCallDiverges"),
    ("truth_j_call", "TruthJCallDiverges"),
)


class SimulationResult:
    """List-like wrapper around a batch of AIRR records.

    ``result[i]`` returns the i-th record dict; ``len(result)`` is
    the number of records; iteration yields records in order.

    The original ``Outcome`` objects (with their full trace +
    revision history) are kept on ``.outcomes`` for advanced
    inspection — most users won't need them.
    """

    __slots__ = ("_records", "_outcomes", "_parents")

    def __init__(
        self,
        records: Sequence[Dict[str, Any]],
        outcomes: Optional[Sequence] = None,
        parents: Optional[Sequence] = None,
    ) -> None:
        self._records: List[Dict[str, Any]] = list(records)
        # ``outcomes`` is optional: callers that built records by
        # other means (e.g. round-tripping a TSV) don't have the
        # underlying Outcome objects available.
        self._outcomes: Optional[List] = (
            list(outcomes) if outcomes is not None else None
        )
        # ``parents`` is the per-clone parent ``Outcome`` list; only
        # populated by ``CompiledClonalExperiment.run_records``.
        # ``None`` for non-clonal results — keep the absence
        # symmetric with the non-clonal record dict shape (no
        # ``clone_id`` / ``parent_id`` columns).
        self._parents: Optional[List] = (
            list(parents) if parents is not None else None
        )

    @classmethod
    def from_outcomes(
        cls,
        outcomes: Sequence,
        refdata: Any,
        *,
        id_prefix: str = "seq",
        expose_provenance: bool = False,
    ) -> "SimulationResult":
        """Build a :class:`SimulationResult` from a list of Rust
        ``Outcome`` objects + the refdata they ran against.

        Each record's ``sequence_id`` is set to ``f"{id_prefix}{i}"``
        (e.g. ``seq0``, ``seq1``, …) so AIRR-format consumers see a
        unique per-row identifier out of the box.

        ``expose_provenance=True`` adds ``truth_v_call``,
        ``truth_d_call``, ``truth_j_call`` columns containing the
        *originally-sampled* allele names — distinct from the
        evidence-driven ``v_call`` / ``d_call`` / ``j_call`` fields,
        which reflect what an aligner would see. Pair them at the
        Python level to compute aligner-vs-truth accuracy without a
        side truth file.
        """
        from ._airr_record import outcome_to_airr_record

        records = [
            outcome_to_airr_record(
                o, refdata, sequence_id=f"{id_prefix}{i}"
            )
            for i, o in enumerate(outcomes)
        ]
        if expose_provenance:
            for outcome, rec in zip(outcomes, records):
                _inject_truth_columns(outcome, refdata, rec)
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

    @property
    def parents(self) -> Optional[List]:
        """Per-clone parent ``Outcome`` objects for clonal results;
        ``None`` for non-clonal results and for results built from
        records directly.

        ``parents[c]`` is the recombination ancestor of clone ``c``:
        every descendant record with
        ``record["clone_id"] == record["parent_id"] == c`` was
        produced by running the post-fork plan from this parent's
        :meth:`final_simulation`.

        The parent ``Outcome`` carries the pre-fork addressed-choice
        ``.trace()``, the pre-fork ``.events()`` ledger, the
        per-revision IR history (``.revision(i)``), and the final
        assembled IR (``.final_simulation()``). Use these for
        replay, lineage analysis, or building a parent-aware family
        validator (Slice 3+ scope).

        The flat ``.outcomes`` list continues to carry **only the
        descendant outcomes** (one entry per AIRR record); parents
        live exclusively here. ``len(.parents)`` equals the clonal
        pipeline's ``n_clones``; ``len(.outcomes)`` equals
        ``n_clones * per_clone``.
        """
        return self._parents

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

    # ── validation ──────────────────────────────────────────────────

    def validate_records(self, refdata: Any) -> "ValidationReport":
        """**Public AIRR output correctness check.**

        Run the postcondition validator over every record in this
        result and return a :class:`ValidationReport`. This is the
        gate a downstream consumer cares about: "is each projected
        AIRR record internally consistent with the engine state that
        produced it?"

        Each record is re-derived independently from its original
        ``Outcome`` (trace + event ledger + final ``Simulation``)
        and compared against the projected dict. A record passes
        when ``outcome.validate_record(refdata, sequence_id=...)``
        returns an empty issue list. Failures collect the
        ``record_index``, ``sequence_id``, and the issue dicts.

        **Companion check** — for engine-side integrity, see
        :meth:`Outcome.check_live_call_cache_parity` (returns the
        cached-vs-fresh divergence on the live-call cache that
        *feeds* projection).

        **Troubleshooting rule** — if a CI run has both this
        validator AND the parity harness failing on the same batch,
        fix the parity divergence FIRST: a stale cache can leak
        into projection and produce spurious validator failures.
        Once parity is green, rerun the validator; remaining
        failures point at a real projection-layer bug.

        ``refdata`` must be the same :class:`RefDataConfig` the
        outcomes were produced against; passing a different refdata
        will misreport mismatches against an unrelated allele pool.

        Raises ``RuntimeError`` when this result was built without
        attached outcomes (e.g. loaded from a TSV); the validator
        needs the engine state, not just the projected record.
        """
        if self._outcomes is None:
            raise RuntimeError(
                "validate_records requires the original Outcome objects "
                "(with their trace + event ledger), but this "
                "SimulationResult was built from records only (e.g. "
                "loaded from TSV). Re-run the simulation with "
                "Experiment.run_records() to get a result with attached "
                "outcomes."
            )
        failures: List[Dict[str, Any]] = []
        for i, (outcome, record) in enumerate(zip(self._outcomes, self._records)):
            sequence_id = str(record.get("sequence_id", ""))
            issues = outcome.validate_record(refdata, sequence_id=sequence_id)
            if issues:
                failures.append(
                    {
                        "record_index": i,
                        "sequence_id": sequence_id,
                        "issues": issues,
                    }
                )
        return ValidationReport(count=len(self._outcomes), failures=failures)

    def validate_families(
        self, refdata: Any = None
    ) -> "FamilyValidationReport":
        """**Clonal-family consistency check** — a strict subset of
        the audit's §6 family-layer invariants
        (``docs/clonal_family_design.md``).

        Groups records by ``clone_id`` and asserts the recombination-
        time truth fields agree across every descendant of a clone.
        ``refdata`` is reserved for forward compatibility with the
        deeper family-layer checks (pre-SHM junction, mutation-
        distance distribution) the audit's later slices add; this
        slice's invariants are all dict-only and ignore ``refdata``.

        **Currently enforced invariants:**

        - ``truth_v_call`` constant within each ``clone_id``, **when
          present**. Skipped silently for clones whose records were
          projected without ``expose_provenance=True``.
        - ``truth_d_call`` same.
        - ``truth_j_call`` same.
        - ``clone_id`` is present on every record once any record
          in the batch carries one (a batch that mixes clonal and
          non-clonal records raises ``CloneIdMissing``).

        **Non-clonal results return ok with ``family_count == 0``.**
        This makes ``result.validate_families()`` a safe no-op on a
        flat batch — the call site does not need to branch on the
        result's clonal-ness.

        **Records-only results work.** Unlike
        :meth:`validate_records`, this validator does not require
        the underlying ``Outcome`` objects — every check is on
        record-dict fields — so a ``SimulationResult`` loaded from
        TSV can still be family-validated.

        **Not enforced yet** (deferred per the audit's §14
        out-of-scope list): mutation-distance distribution, pre-SHM
        junction invariance, parent-trace reconstruction, lineage
        topology, ``original_v_call`` / ``d_inverted`` invariance.

        Returns a :class:`FamilyValidationReport` carrying
        ``count``, ``family_count``, ``members_per_family``, and
        ``failures``.
        """
        del refdata  # forward-compat placeholder, see docstring.

        records = self._records
        total = len(records)

        # Identify whether any record carries a non-null ``clone_id``.
        # A batch that has no clonal records returns an ok no-op
        # report; a batch where SOME records carry a clone_id and
        # others don't is structurally broken (mixed clonal/flat
        # outputs) and surfaces a ``CloneIdMissing`` failure.
        any_clonal = any(
            "clone_id" in r and r["clone_id"] is not None for r in records
        )
        if not any_clonal:
            return FamilyValidationReport(
                count=total,
                family_count=0,
                members_per_family={},
                failures=[],
            )

        failures: List[Dict[str, Any]] = []
        missing_clone_id: List[int] = [
            i for i, r in enumerate(records)
            if r.get("clone_id") is None
        ]
        if missing_clone_id:
            failures.append(
                {
                    "clone_id": None,
                    "record_indices": missing_clone_id,
                    "issue_kind": "CloneIdMissing",
                    "values": [],
                }
            )

        # Group records by clone_id (ignoring records missing the tag
        # — they already surfaced above).
        by_clone: Dict[Any, List[int]] = {}
        for i, r in enumerate(records):
            cid = r.get("clone_id")
            if cid is None:
                continue
            by_clone.setdefault(cid, []).append(i)

        for cid, indices in by_clone.items():
            members = [records[i] for i in indices]
            for field, kind in _FAMILY_TRUTH_INVARIANTS:
                # Skip the field for this group if no record carries
                # it (consumer opted out of expose_provenance). Once
                # any record in the group carries it, every record
                # in the group must — the consistency check still
                # runs but the "missing on some" surface would be a
                # different bug; we treat ``None`` and absent as the
                # same "no opinion" signal here.
                present = [
                    rec.get(field) for rec in members
                    if field in rec and rec.get(field) is not None
                ]
                if not present:
                    continue
                distinct = sorted({str(v) for v in present})
                if len(distinct) > 1:
                    failures.append(
                        {
                            "clone_id": cid,
                            "record_indices": list(indices),
                            "issue_kind": kind,
                            "values": distinct,
                        }
                    )

        members_per_family = {cid: len(idxs) for cid, idxs in by_clone.items()}
        return FamilyValidationReport(
            count=total,
            family_count=len(by_clone),
            members_per_family=members_per_family,
            failures=failures,
        )

    def validate_families_with_parents(
        self, refdata: Any = None
    ) -> "FamilyValidationReport":
        """**Parent-aware clonal-family validator** — the deeper
        diagnostic that compares every descendant against its
        actual parent ``Outcome`` (Slice 3 of the clonal-family
        audit; see ``docs/clonal_parent_outcome_design.md`` §6).

        Sibling of :meth:`validate_families`. That validator is
        record-only (groups by ``clone_id``, compares truth fields
        across siblings); this one **requires the parent outcomes
        to be available on the result** and compares each
        descendant against its parent directly. Use this when you
        want to confirm "the descendants reflect the recombination
        ancestor they came from," not just "siblings agree with
        each other."

        **Currently enforced invariants** — all derived from
        record-vs-parent comparison only:

        - Structural:
          - ``ParentsMissing`` when records carry ``clone_id`` /
            ``parent_id`` but ``self.parents`` is ``None``.
          - ``ParentIdMissing`` for records without a non-null
            ``parent_id`` in a result that has parents available.
          - ``ParentIdOutOfRange`` when ``record["parent_id"]``
            is not in ``[0, len(self.parents))``.
        - Truth-allele consistency (requires ``refdata``):
          - ``ParentTruthVCallMismatch`` /
            ``ParentTruthDCallMismatch`` /
            ``ParentTruthJCallMismatch`` — descendant's
            ``truth_*_call`` (from ``expose_provenance=True``)
            disagrees with the parent's projected truth allele.
        - Provenance consistency (no ``refdata`` needed):
          - ``ParentDInvertedMismatch`` — descendant's
            ``d_inverted`` disagrees with the parent's. D inversion
            is a pre-fork decision, so divergence indicates a
            structural bug.
          - ``ParentOriginalVCallMismatch`` — descendant's
            ``original_v_call`` (receptor-revision provenance)
            disagrees with the parent's. Same reasoning.

        **Without ``refdata``:** only the structural checks
        (``ParentsMissing`` / ``ParentIdMissing`` /
        ``ParentIdOutOfRange``) run. All value comparisons —
        truth alleles, ``d_inverted``, ``original_v_call`` —
        require projecting the parent ``Outcome`` to an AIRR
        record, which today goes through the Rust projector and
        needs ``refdata``. Slice 3 deliberately stays Python-only;
        a lighter-weight refdata-free parent accessor for
        ``d_inverted`` etc. is deferred until a Rust slice surfaces
        one.

        **Skipped silently for fields not present on descendants:**
        if ``expose_provenance`` was off, the descendants don't
        carry ``truth_*_call`` and those checks are no-ops.

        **Non-clonal results return ok with ``family_count=0``** —
        same safe no-op shape as :meth:`validate_families`. Slice
        3 deliberately does NOT raise "not clonal" here.

        **Not enforced yet** (deferred per audit §6, §14):

        - **Pre-SHM junction invariance.** The descendant's
          ``junction`` AIRR field is post-SHM; the pre-SHM junction
          lives only inside the parent's IR. A proper check would
          require a parent-derived ``junction_pre_shm`` field on
          either records or a future ``FamilyRecord`` projection
          (Slice 4+). This validator does NOT compare
          ``descendant.junction`` against any parent-derived
          value today.
        - **Mutation-distance distribution.** Comparing the
          parent's assembled sequence to each descendant's
          post-SHM sequence to verify SHM mass is plausible.
          Requires projecting the parent's pool to a sequence
          string — out of scope for this slice.
        - **Plan-split pre-fork pass enumeration.** "Parent should
          not carry descendant-only observation fields like PCR
          / paired-end / quality errors" is pinned at the
          contract-test level (the pre-fork plan's pass names)
          rather than enforced at runtime here.

        **Not wired into ``validate_records=True``.** This is an
        explicit deeper diagnostic surface. The
        ``validate_records=True`` gate continues to run only the
        per-record postcondition validator and the field-only
        :meth:`validate_families`. Callers who want parent-aware
        checks invoke this method explicitly.

        Returns a :class:`FamilyValidationReport`. Failure dicts
        carry ``clone_id``, ``parent_id``, ``record_indices``,
        ``issue_kind``, ``parent_value``, and ``child_values``
        (the latter two are ``None`` / ``[]`` for structural
        failures that don't compare values).
        """
        records = self._records
        total = len(records)

        any_clonal = any(
            "clone_id" in r and r["clone_id"] is not None for r in records
        )
        any_parent_id = any(
            "parent_id" in r and r["parent_id"] is not None for r in records
        )

        # Non-clonal: ok no-op (matches validate_families).
        if not any_clonal and not any_parent_id:
            return FamilyValidationReport(
                count=total,
                family_count=0,
                members_per_family={},
                failures=[],
            )

        failures: List[Dict[str, Any]] = []
        parents = self._parents

        # Parents missing entirely — surface and bail (no comparable
        # parent state to run further checks against).
        if parents is None:
            failures.append(
                {
                    "clone_id": None,
                    "parent_id": None,
                    "record_indices": [
                        i for i, r in enumerate(records)
                        if r.get("clone_id") is not None
                        or r.get("parent_id") is not None
                    ],
                    "issue_kind": "ParentsMissing",
                    "parent_value": None,
                    "child_values": [],
                }
            )
            # Compute aggregation by clone_id so the report still
            # carries useful structure even though no comparison
            # ran.
            by_clone: Dict[Any, List[int]] = {}
            for i, r in enumerate(records):
                cid = r.get("clone_id")
                if cid is None:
                    continue
                by_clone.setdefault(cid, []).append(i)
            return FamilyValidationReport(
                count=total,
                family_count=len(by_clone),
                members_per_family={
                    cid: len(idxs) for cid, idxs in by_clone.items()
                },
                failures=failures,
            )

        # Parents are present. Group descendants by parent_id and
        # surface structural failures (missing / out of range) per
        # record.
        n_parents = len(parents)
        missing_parent_id: List[int] = []
        out_of_range: Dict[Any, List[int]] = {}
        by_parent: Dict[int, List[int]] = {}
        for i, r in enumerate(records):
            pid = r.get("parent_id")
            if pid is None:
                # Records without a parent_id in a clonal result are
                # structurally broken — surface them and skip per-
                # parent invariant checks for those records.
                if r.get("clone_id") is not None:
                    missing_parent_id.append(i)
                continue
            if not isinstance(pid, int) or isinstance(pid, bool):
                out_of_range.setdefault(pid, []).append(i)
                continue
            if pid < 0 or pid >= n_parents:
                out_of_range.setdefault(pid, []).append(i)
                continue
            by_parent.setdefault(pid, []).append(i)

        if missing_parent_id:
            failures.append(
                {
                    "clone_id": None,
                    "parent_id": None,
                    "record_indices": missing_parent_id,
                    "issue_kind": "ParentIdMissing",
                    "parent_value": None,
                    "child_values": [],
                }
            )
        for pid, idxs in out_of_range.items():
            failures.append(
                {
                    "clone_id": None,
                    "parent_id": pid,
                    "record_indices": idxs,
                    "issue_kind": "ParentIdOutOfRange",
                    "parent_value": pid,
                    "child_values": [],
                }
            )

        # Build parent projections lazily — only for parent ids that
        # have descendants pointing at them. The Rust AIRR projector
        # ``outcome_to_airr_record`` requires refdata; without it we
        # can't materialize the parent's projected fields and
        # therefore can't run any value-comparison check. Structural
        # checks (missing / out-of-range parent_id) already ran
        # above. Document this in the method's contract: refdata-
        # less calls run structural checks only.
        parent_projections: Dict[int, Dict[str, Any]] = {}
        if by_parent and refdata is not None:
            from ._airr_record import outcome_to_airr_record

            for pid in by_parent:
                proj = outcome_to_airr_record(
                    parents[pid],
                    refdata,
                    sequence_id=f"parent{pid}",
                )
                _inject_truth_columns(parents[pid], refdata, proj)
                parent_projections[pid] = proj

        # The (field, issue_kind) pairs the parent-aware validator
        # enforces when ``refdata`` is available. All of these are
        # comparisons against parent-projected values; refdata is
        # required to materialize them.
        FIELD_CHECKS = [
            ("truth_v_call", "ParentTruthVCallMismatch"),
            ("truth_d_call", "ParentTruthDCallMismatch"),
            ("truth_j_call", "ParentTruthJCallMismatch"),
            ("d_inverted", "ParentDInvertedMismatch"),
            ("original_v_call", "ParentOriginalVCallMismatch"),
        ]

        for pid, indices in by_parent.items():
            parent_proj = parent_projections.get(pid)
            if parent_proj is None:
                # Refdata was None — value comparisons skipped.
                continue
            members = [records[i] for i in indices]
            # Use the first member's clone_id for the failure
            # group; per-clone_id divergence within a parent_id
            # group would have surfaced via the family validator
            # already.
            clone_id_for_group = members[0].get("clone_id", pid)
            for field, kind in FIELD_CHECKS:
                parent_val = parent_proj.get(field)
                # Skip when parent has no opinion on the field
                # (e.g. ``original_v_call`` is empty when receptor
                # revision didn't run; comparing "" to "" would
                # always pass but comparing "" to a real allele
                # name would surface a real bug — so we only skip
                # when the field is genuinely missing).
                if parent_val is None:
                    continue
                divergent_indices: List[int] = []
                divergent_values: List[Any] = []
                for idx, rec in zip(indices, members):
                    if field not in rec:
                        continue
                    child_val = rec[field]
                    if child_val is None:
                        continue
                    if child_val != parent_val:
                        divergent_indices.append(idx)
                        if child_val not in divergent_values:
                            divergent_values.append(child_val)
                if divergent_indices:
                    failures.append(
                        {
                            "clone_id": clone_id_for_group,
                            "parent_id": pid,
                            "record_indices": divergent_indices,
                            "issue_kind": kind,
                            "parent_value": parent_val,
                            "child_values": sorted(
                                divergent_values, key=lambda v: str(v)
                            ),
                        }
                    )

        return FamilyValidationReport(
            count=total,
            family_count=len(by_parent),
            members_per_family={
                pid: len(idxs) for pid, idxs in by_parent.items()
            },
            failures=failures,
        )

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


class SimulationResultWithLineages(SimulationResult):
    """A :class:`SimulationResult` that also carries per-clone lineage trees.

    Produced by :meth:`CompiledLineageExperiment.run_records`. Adds a
    ``.lineage_trees`` property that exposes the raw
    :class:`~GenAIRR._engine.LineageTree` objects (one per clone) for
    ground-truth export via ``.to_newick()``, ``.to_fasta()``, and
    ``.to_node_table_tsv()``.
    """

    __slots__ = ("_lineage_trees",)

    def __init__(
        self,
        records: "Sequence[Dict[str, Any]]",
        outcomes: "Optional[Sequence]" = None,
        parents: "Optional[Sequence]" = None,
        lineage_trees: "Optional[Sequence]" = None,
    ) -> None:
        super().__init__(records, outcomes, parents)
        self._lineage_trees: "Optional[List]" = (
            list(lineage_trees) if lineage_trees is not None else None
        )

    @property
    def lineage_trees(self) -> "Optional[List]":
        """Per-clone :class:`~GenAIRR._engine.LineageTree` objects, or ``None``.

        Each tree supports ``.validate()``, ``.to_newick()``,
        ``.to_fasta()``, and ``.to_node_table_tsv()`` for ground-truth
        export and downstream phylogenetic analysis.
        """
        return self._lineage_trees
