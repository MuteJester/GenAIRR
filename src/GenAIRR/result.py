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

from typing import Any, Dict, Iterator, List, Optional, Sequence, Union, overload

from ._result_common import _inject_truth_columns
from ._result_export import (
    _ResultExport,
    _DEFAULT_COLUMN_ORDER,
    _to_airr_strict,
)
from ._result_validation import (
    _ResultValidation,
    ValidationReport,
    FamilyValidationReport,
)


class SimulationResult(_ResultExport, _ResultValidation):
    """List-like wrapper around a batch of AIRR records.

    ``result[i]`` returns the i-th record dict; ``len(result)`` is
    the number of records; iteration yields records in order.

    The original ``Outcome`` objects (with their full trace +
    revision history) are kept on ``.outcomes`` for advanced
    inspection — most users won't need them.
    """

    __slots__ = ("_records", "_outcomes", "_parents", "_genotypes")

    def __init__(
        self,
        records: Sequence[Dict[str, Any]],
        outcomes: Optional[Sequence] = None,
        parents: Optional[Sequence] = None,
    ) -> None:
        self._records: List[Dict[str, Any]] = list(records)
        # Per-subject ground-truth ``Genotype`` objects, populated by
        # ``CompiledExperiment.run_records`` when a genotype is attached.
        # ``None`` for non-genotype results.
        self._genotypes: Optional[List] = None
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
    def genotypes(self) -> Optional[List]:
        """Per-subject ground-truth ``Genotype`` objects when the
        experiment had a genotype attached, else ``None``."""
        return self._genotypes

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
