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

from ._result_common import _as_native_refdata, _inject_truth_columns
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

    Results produced by a simulation are bound to the native reference
    they ran against (``.reference_content_hash``); validating them
    against a different reference raises ``ValueError``.
    """

    __slots__ = ("_records", "_outcomes", "_parents", "_genotypes", "_reference_content_hash")

    def __init__(
        self,
        records: Sequence[Dict[str, Any]],
        outcomes: Optional[Sequence] = None,
        parents: Optional[Sequence] = None,
        *,
        reference_content_hash: Optional[str] = None,
    ) -> None:
        if reference_content_hash is not None and not isinstance(reference_content_hash, str):
            raise TypeError("reference_content_hash must be a string or None")
        self._reference_content_hash: Optional[str] = reference_content_hash
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
        # populated by ``CompiledClonalExperiment.run_records`` and
        # duplex repertoire/lineage runs.
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
        ``truth_d_call``, and ``truth_j_call`` compatibility aliases
        for the standard generating calls. Post-mutation inference
        is always available separately as ``recalled_*_call``.

        ``refdata`` may also be a :class:`GenAIRR.Reference`. The result
        is bound to the reference's content hash.
        """
        from ._airr_record import outcome_to_airr_record

        refdata = _as_native_refdata(refdata)
        records = [
            outcome_to_airr_record(
                o, refdata, sequence_id=f"{id_prefix}{i}"
            )
            for i, o in enumerate(outcomes)
        ]
        if expose_provenance:
            for outcome, rec in zip(outcomes, records):
                _inject_truth_columns(outcome, refdata, rec)
        return cls(records, outcomes=outcomes, reference_content_hash=refdata.content_hash())

    # ── list-like access ────────────────────────────────────────────

    @property
    def reference_content_hash(self) -> Optional[str]:
        """Content hash of the native reference these records were generated
        against (for a portable reference, its ``execution_hash``), or
        ``None`` when the result was built from records alone."""
        return getattr(self, "_reference_content_hash", None)

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
        """Per-clone founder ``Outcome`` objects for fixed-size clonal and
        duplex repertoire/lineage results;
        ``None`` for non-clonal results and for results built from
        records directly.

        ``parents[record["parent_id"]]`` is the record's recombination
        ancestor. The dense parent index and ``clone_id`` have separate
        meanings; consumers should not assume they are equal.

        The parent ``Outcome`` carries the pre-fork addressed-choice
        ``.trace()``, the pre-fork ``.events()`` ledger, the
        per-revision IR history (``.revision(i)``), and the final
        assembled IR (``.final_simulation()``). Use these for
        replay, lineage analysis, or building a parent-aware family
        validator (Slice 3+ scope).

        The flat ``.outcomes`` list continues to carry **only the
        descendant outcomes** (one entry per AIRR record). Parent outcomes
        live here. Duplex lineage runs retain accepted extinct founders as
        well, even when their trees produce no observed records.
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

    def save_family_replay(self, path, *, replace=False):
        """Save a sealed, captured duplex lineage run for self-contained replay.

        Run with ``capture_family_replay=True`` first. Changed records or
        annotations reject publication; existing destinations require the
        explicit ``replace=True`` option.
        """
        from ._family_replay import _bundle_from_result, FamilyReplayError

        if getattr(self, "_family_replay_capture", None) is None:
            raise FamilyReplayError("missing_capture", "result", "run a supported lineage with capture_family_replay=True")
        _bundle_from_result(self).write_to(path, replace=replace)


class SimulationResultWithLineages(SimulationResult):
    """A :class:`SimulationResult` that also carries per-clone lineage trees.

    Produced by :meth:`CompiledLineageExperiment.run_records`. Adds a
    ``.lineage_trees`` property that exposes the raw
    :class:`~GenAIRR._engine.LineageTree` objects (one per clone) for
    ground-truth export via ``.to_newick()``, ``.to_fasta()``, and
    ``.to_node_table_tsv()``.
    """

    __slots__ = (
        "_lineage_trees", "_families", "_lineage_record_origins",
        "_lineage_origin_trees", "_lineage_origin_nodes",
        "_family_replay_capture", "_family_replay_seal", "_family_replay_refdata",
    )

    def __init__(
        self,
        records: "Sequence[Dict[str, Any]]",
        outcomes: "Optional[Sequence]" = None,
        parents: "Optional[Sequence]" = None,
        lineage_trees: "Optional[Sequence]" = None,
        families: "Optional[Sequence]" = None,
        lineage_record_origins: "Optional[Sequence]" = None,
        *,
        reference_content_hash: "Optional[str]" = None,
    ) -> None:
        super().__init__(records, outcomes, parents, reference_content_hash=reference_content_hash)
        self._family_replay_capture = None
        self._family_replay_seal = None
        self._family_replay_refdata = None
        self._lineage_trees: "Optional[List]" = (
            list(lineage_trees) if lineage_trees is not None else None
        )
        self._families: "Optional[List]" = (
            list(families) if families is not None else None
        )
        # Immutable origin addressing lets validation detect edited namespace
        # tags even when different cells happen to have identical DNA/history.
        self._lineage_record_origins = (
            tuple(tuple(origin) for origin in lineage_record_origins)
            if lineage_record_origins is not None else None
        )
        self._lineage_origin_trees = (
            tuple(self._lineage_trees or ())
            if lineage_record_origins is not None else None
        )
        if self._lineage_record_origins is None:
            self._lineage_origin_nodes = None
        else:
            # Read each native tree once, then retain only the observed nodes
            # needed to check record metadata without cloning a whole arena
            # for every record during validation.
            native_nodes = [tree.nodes() for tree in self._lineage_origin_trees]
            self._lineage_origin_nodes = tuple(
                native_nodes[parent_id][node_id]
                for parent_id, _, node_id in self._lineage_record_origins
            )

    @property
    def refdata(self):
        """Exact native references for captured/replayed duplex families, else None."""
        return self._family_replay_refdata

    @property
    def families(self) -> "Optional[List]":
        """Native family outcomes for duplex runs, aligned with lineage trees.

        ``families[i].replay_node(node_id, refdata)`` replays a retained cell's
        actual mutation traces from its founder without new random choices.
        This is an in-memory replay interface; it is not a durable family
        file format. Descendant repair IDs describe inherited founder
        ancestry. Insertion IDs are local to a branch, identified by its
        clone and lineage-node fields.
        """
        return self._families

    @property
    def lineage_trees(self) -> "Optional[List]":
        """Per-clone :class:`~GenAIRR._engine.LineageTree` objects, or ``None``.

        Each tree supports ``.validate()``, ``.to_newick()``,
        ``.to_fasta()``, and ``.to_node_table_tsv()`` for ground-truth
        export and downstream phylogenetic analysis.
        """
        return self._lineage_trees
