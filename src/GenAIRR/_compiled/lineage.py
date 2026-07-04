"""Runtime wrapper for the BCR affinity-maturation lineage experiment.

:class:`CompiledLineageExperiment` grows one lineage tree per clone via
the Rust ``simulate_family_outcomes`` kernel and returns per-observed-node
AIRR records with lineage metadata. It is intentionally thin: it holds a
reference to the compiled engine simulators and source context so record
export can render a faithful narrative. The execution loop delegates
straight into the engine; the wrapper exists for ergonomics and to keep
the public surface stable across engine refactors.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, TYPE_CHECKING

from GenAIRR import _engine  # private Rust extension submodule

from .._pipeline_ir import _LineageForkStep

if TYPE_CHECKING:
    from ..dataconfig import DataConfig
    from ..result import SimulationResultWithLineages


class CompiledLineageExperiment:
    """A compiled experiment that grows BCR affinity-maturation lineage trees.

    Wraps a pre-fork :class:`GenAIRR._engine.CompiledSimulator` (founder
    recombination) and a :class:`~GenAIRR._pipeline_ir._LineageForkStep`
    (lineage parameters). :meth:`run_records` grows one lineage tree per
    clone via the Rust ``simulate_family_outcomes`` kernel and returns a
    :class:`~GenAIRR.result.SimulationResultWithLineages` whose records are
    per-observed-node AIRR dicts with lineage metadata.
    """

    __slots__ = (
        "_pre",
        "_step",
        "_refdata",
        "_post",
        "_post_steps",
        "_dataconfig",
        "_metadata",
    )

    def __init__(
        self,
        pre_simulator: "_engine.CompiledSimulator",
        step: "_LineageForkStep",
        refdata: "_engine.RefDataConfig",
        post_simulator: Optional["_engine.CompiledSimulator"] = None,
        post_steps: Sequence[Any] = (),
        dataconfig: Optional["DataConfig"] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        self._pre = pre_simulator
        self._step = step
        self._refdata = refdata
        # Per-observed-cell library-prep / sequencing corruption plan.
        # ``None`` keeps the pristine-read path (no artefacts applied).
        self._post = post_simulator
        self._post_steps: Tuple[Any, ...] = tuple(post_steps)
        self._dataconfig = dataconfig
        self._metadata = dict(metadata) if metadata else {}

    @property
    def refdata(self) -> "_engine.RefDataConfig":
        return self._refdata

    def run_records(
        self,
        *,
        seed: int = 0,
        strict: bool = False,
        expose_provenance: bool = False,
        validate_records: bool = False,
    ) -> "SimulationResultWithLineages":
        """Grow lineage trees and return per-observed-node AIRR records.

        Each clone is seeded at ``seed + clone_idx * 1_000_000`` so
        independent clones are reproducible and non-overlapping.

        The returned :class:`SimulationResultWithLineages` carries the
        per-record ``Outcome`` that produced each AIRR record on its
        ``.outcomes`` attribute (index-aligned with ``.records``). On the
        pristine-read path these are the per-observed-node SHM
        ``Outcome`` objects (``fam.observed_outcomes()``); on the
        corruption path they are the per-node merged
        recombination+SHM+artefact ``Outcome`` objects. Since each node
        ``Outcome`` is self-consistent, ``result.validate_records(refdata)``
        passes.

        ``expose_provenance=True`` appends ``truth_v_call`` /
        ``truth_d_call`` / ``truth_j_call`` columns to each record from
        the founder allele assignments carried on the per-record
        ``Outcome``.

        ``validate_records=True`` runs
        :meth:`SimulationResult.validate_records` (per-record
        postcondition) and :meth:`SimulationResult.validate_families`
        (clonal-family consistency by ``clone_id``) on the freshly built
        batch, raising the matching validation error on any failure.
        """
        from GenAIRR import _engine as _eng
        from .._airr_record import outcome_to_airr_record
        from .._s5f_loader import load_builtin_s5f_kernel
        from ..result import SimulationResultWithLineages, _inject_truth_columns

        step = self._step
        mutability, substitution = load_builtin_s5f_kernel(step.s5f_model)

        records: List[Dict[str, Any]] = []
        # Per-record source ``Outcome`` objects, index-aligned with
        # ``records`` so ``result.validate_records`` can re-derive each
        # record from the engine state that produced it.
        outcomes: List["_engine.Outcome"] = []
        lineage_trees = []

        # Bound on founder-survival retries. Sampling draws from the LIVING
        # final-generation population, so an extinct founder (drew 0 offspring)
        # yields zero observed cells. With ``allow_extinction=False`` (default)
        # we condition each clone on survival by re-growing the family with a
        # fresh deterministic sub-seed until it survives or we exhaust the
        # bound. The sub-seed offset (a prime times the attempt index) keeps the
        # retries deterministic — same top-level seed => same result.
        _MAX_SURVIVAL_ATTEMPTS = 50
        _SUBSEED_STRIDE = 7_919  # prime; keeps retry sub-seeds well-separated
        allow_extinction = getattr(step, "allow_extinction", False)

        for clone_idx in range(step.n_clones):
            clone_seed = int(seed) + clone_idx * 1_000_000
            fam = None
            tree = None
            for attempt in range(_MAX_SURVIVAL_ATTEMPTS):
                family_seed = clone_seed + attempt * _SUBSEED_STRIDE
                # _pre.run() returns a single Outcome (not a list) — mirror
                # CompiledClonalExperiment which calls self._pre.run(seed=...).
                founder = self._pre.run(seed=family_seed, strict=strict)
                candidate = _eng.simulate_family_outcomes(
                    founder,
                    self._refdata,
                    mutability,
                    substitution,
                    step.rate,
                    step.lambda_base,
                    0.0,      # lambda_mut: positional slot 6 (inert; hardcoded)
                    step.max_generations,
                    step.n_max,
                    step.n_sample,
                    family_seed,
                    selection_strength=step.selection_strength,
                    beta=step.beta,
                    target_aa=step.target_aa,
                    mature_substitutions=step.mature_substitutions,
                )
                survived = len(candidate.observed_outcomes()) > 0
                if survived or allow_extinction:
                    fam = candidate
                    tree = candidate.tree()
                    break
                # extinct + survival required => retry with next sub-seed
            else:
                # Exhausted the retry bound without a surviving family.
                raise ValueError(
                    f"clonal_lineage: clone {clone_idx} went extinct on every "
                    f"one of {_MAX_SURVIVAL_ATTEMPTS} survival attempts (each "
                    "founder drew 0 offspring). Increase lambda_base (offspring "
                    "Poisson mean) and/or max_generations so families reliably "
                    "survive, or pass allow_extinction=True to accept extinct "
                    "clones (producing fewer than n_clones families)."
                )

            if fam is None:
                # allow_extinction=True and the final attempt was extinct:
                # skip this clone (it contributes no observed cells/records).
                continue

            lineage_trees.append(tree)
            if len(fam.observed_outcomes()) == 0:
                # allow_extinction=True and this family is extinct: keep its
                # (empty) tree for export parity but emit no records.
                continue

            if self._post is None:
                # Pristine-read path: project each observed node's
                # synthesized SHM Outcome straight to an AIRR record.
                # ``observed_outcomes()`` is index-aligned with both the
                # observed-node list and ``airr_records()``.
                observed_nodes = [n for n in tree.nodes() if n.observed]
                recs = fam.airr_records(self._refdata)
                observed_outcomes = fam.observed_outcomes()
                for node, rec, node_outcome in zip(
                    observed_nodes, recs, observed_outcomes
                ):
                    self._stamp_lineage_metadata(rec, clone_idx, node)
                    if expose_provenance and node_outcome is not None:
                        _inject_truth_columns(node_outcome, self._refdata, rec)
                    records.append(rec)
                    outcomes.append(node_outcome)
                continue

            # Corruption path: per observed node, run the library-prep /
            # sequencing corruption plan FROM the node's post-SHM
            # simulation, then merge the founder-recombination + SHM
            # provenance with the corruption trace / events so the AIRR
            # record reports trims, v/d/j, SHM counts AND artefact
            # counters (n_quality_errors, n_pcr_errors, n_indels, …).
            node_outcomes = fam.node_outcomes()
            for node, base_outcome in zip(tree.nodes(), node_outcomes):
                if base_outcome is None:
                    continue
                node_seed = clone_seed + 1 + node.id
                corruption_outcome = self._post.run_from(
                    base_outcome.final_simulation(), node_seed, strict=strict
                )
                merged = _eng.merge_lineage_corruption(
                    base_outcome, corruption_outcome
                )
                rec = outcome_to_airr_record(
                    merged,
                    self._refdata,
                    sequence_id=f"clone{clone_idx}_node{node.id}",
                )
                self._stamp_lineage_metadata(rec, clone_idx, node)
                if expose_provenance:
                    _inject_truth_columns(merged, self._refdata, rec)
                records.append(rec)
                outcomes.append(merged)

        result = SimulationResultWithLineages(
            records, outcomes=outcomes, lineage_trees=lineage_trees
        )
        if validate_records:
            from .._validation import (
                _raise_on_family_validation_failure,
                _raise_on_validation_failure,
            )

            _raise_on_validation_failure(result.validate_records(self._refdata))
            _raise_on_family_validation_failure(result.validate_families())
        return result

    @staticmethod
    def _stamp_lineage_metadata(rec: Dict[str, Any], clone_idx: int, node) -> None:
        """Stamp clone + lineage-node provenance onto an AIRR record."""
        rec["clone_id"] = clone_idx
        rec["lineage_node_id"] = node.id
        rec["lineage_parent_id"] = (
            node.parent_id if node.parent_id is not None else -1
        )
        rec["lineage_generation"] = node.generation
        rec["lineage_abundance"] = node.abundance
        # AIRR-standard abundance field that abundance-aware tools (Change-O /
        # SCOPer / dowser) read. Mirror lineage_abundance so the genotype-
        # collapsed observed cell carries its represented count both ways.
        rec["duplicate_count"] = node.abundance
        rec["lineage_affinity"] = node.affinity
        rec["sequence_id"] = f"clone{clone_idx}_node{node.id}"

    def __repr__(self) -> str:
        return (
            f"<CompiledLineageExperiment n_clones={self._step.n_clones} "
            f"max_gen={self._step.max_generations} chain={self._refdata.chain_type}>"
        )
