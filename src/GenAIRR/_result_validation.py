"""Record/family validation surface for SimulationResult: the public
ValidationReport / FamilyValidationReport result types plus the validate_*
methods extracted verbatim as a mixin (behavior-preserving).
"""
from __future__ import annotations

from typing import Any, Dict, List

from ._result_common import _inject_truth_columns


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


class _ResultValidation:
    """Record/family validation mixin for :class:`SimulationResult`.

    Behavior-preserving extraction: every method resolves ``self``
    (``self._records`` / ``self._outcomes`` / ``self._parents``) via
    the MRO of the concrete :class:`SimulationResult` that inherits
    this mixin.
    """

    __slots__ = ()

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
