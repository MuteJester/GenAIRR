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
        "_contracts",
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
        contracts: Optional["_engine.ContractSet"] = None,
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
        self._contracts = contracts

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
        capture_family_replay: bool = False,
        family_replay_limits=None,
        _family_replay_request=None,
    ) -> "SimulationResultWithLineages":
        """Grow lineage trees and return per-observed-node AIRR records.

        Each clone is seeded at ``seed + clone_idx * 1_000_000`` so
        independent clones are reproducible and non-overlapping. Declared
        productivity contracts constrain each division's mutation support,
        as well as the founder and any post-lineage observation passes.

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

        ``capture_family_replay=True`` records and seals a self-contained duplex
        family history for :meth:`SimulationResult.save_family_replay`. It
        validates all native evidence even when ``validate_records`` is false.
        The optional ``family_replay_limits`` bounds capture work and storage.
        """
        from GenAIRR import _engine as _eng
        from .._airr_record import outcome_to_airr_record
        from .._s5f_loader import load_builtin_s5f_kernel
        from ..result import SimulationResultWithLineages, _inject_truth_columns

        if type(capture_family_replay) is not bool:
            raise TypeError("capture_family_replay must be Boolean")
        capture = None
        if capture_family_replay:
            from .._family_replay import _capture_request, _start_capture

            request = (_family_replay_request if _family_replay_request is not None
                       else _capture_request({}, family_replay_limits))
            capture = _start_capture(self, request, strict, expose_provenance, seed)
        elif family_replay_limits is not None or _family_replay_request is not None:
            raise ValueError("family replay limits/context require capture_family_replay=True")

        step = self._step
        mutability, substitution = load_builtin_s5f_kernel(step.s5f_model)

        records: List[Dict[str, Any]] = []
        # Per-record source ``Outcome`` objects, index-aligned with
        # ``records`` so ``result.validate_records`` can re-derive each
        # record from the engine state that produced it.
        outcomes: List["_engine.Outcome"] = []
        lineage_trees = []
        duplex = "duplex.recombine" in self._pre.pass_names()
        parents = [] if duplex else None
        families = [] if duplex else None
        # Collected on every path, not just duplex. These are the immutable
        # (family, clone, node) addresses validation resolves a record back to
        # its tree node through; gating them on duplex meant that on the
        # legacy path no lineage field was validated at all, and all eight of
        # them could be edited with validate_records still returning ok.
        lineage_record_origins = []

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
            # Repair is a once-per-founder choice. Conditioning duplex growth
            # on survival must not reweight the source history by redrawing it.
            if capture is not None:
                from .._family_replay import _native, _truth_aliases

                founder = _native(f"clones[{clone_idx}].source", self._pre.run_for_family_capture,
                                  clone_seed, strict, capture.budget)
            else:
                founder = self._pre.run(seed=clone_seed, strict=strict) if duplex else None
            captured_clone = capture.source(self._pre, founder, clone_idx) if capture is not None else None
            fam = None
            tree = None
            for attempt in range(_MAX_SURVIVAL_ATTEMPTS):
                family_seed = clone_seed + attempt * _SUBSEED_STRIDE
                # _pre.run() returns a single Outcome (not a list) — mirror
                # CompiledClonalExperiment which calls self._pre.run(seed=...).
                if not duplex:
                    founder = self._pre.run(seed=family_seed, strict=strict)
                capture_kwargs = ({"capture_family_replay": True, "budget": capture.budget}
                                  if capture is not None else {})

                def grow_candidate():
                    return _eng.simulate_family_outcomes(
                        founder,
                        self._refdata,
                        mutability,
                        substitution,
                        step.rate,
                        step.lambda_base,
                        0.0,  # lambda_mut: positional slot 6 (inert; hardcoded)
                        step.max_generations,
                        step.n_max,
                        step.n_sample,
                        family_seed,
                        selection_strength=step.selection_strength,
                        beta=step.beta,
                        target_aa=step.target_aa,
                        mature_substitutions=step.mature_substitutions,
                        respect=self._contracts,
                        strict=strict,
                        **capture_kwargs,
                    )
                candidate = (_native(f"clones[{clone_idx}].attempts[{attempt}]", grow_candidate)
                             if capture is not None else grow_candidate())
                if capture is not None:
                    capture.attempt(captured_clone, candidate)
                survived = (_native(f"clones[{clone_idx}].attempts[{attempt}]", _eng._family_replay_has_observed,
                                    candidate, budget=capture.budget)
                            if capture is not None else len(candidate.observed_outcomes()) > 0)
                if survived or allow_extinction:
                    fam = candidate
                    tree = candidate.tree()
                    if capture is not None:
                        captured_clone["accepted_attempt"] = attempt
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
            parent_id = len(parents) if parents is not None else None
            if parents is not None:
                parents.append(founder)
                families.append(fam)
            has_records = survived if capture is not None else len(fam.observed_outcomes()) > 0
            if not has_records:
                # allow_extinction=True and this family is extinct: keep its
                # (empty) tree for export parity but emit no records.
                continue

            if self._post is None:
                # Pristine-read path: project each observed node's native
                # SHM Outcome straight to an AIRR record.
                # ``observed_outcomes()`` is index-aligned with both the
                # observed-node list and ``airr_records()``.
                observed_nodes = [n for n in tree.nodes() if n.observed]
                if capture is not None:
                    observed_outcomes = _native(f"clones[{clone_idx}].nodes", _eng._family_replay_observed_outcomes,
                                                fam, budget=capture.budget)
                    recs = (_native(f"clones[{clone_idx}].observations[{node.id}]", _eng._family_replay_airr_record,
                                    outcome, self._refdata, f"clone{clone_idx}_node{node.id}", budget=capture.budget)
                            for node, outcome in zip(observed_nodes, observed_outcomes))
                else:
                    recs = fam.airr_records(self._refdata)
                    observed_outcomes = fam.observed_outcomes()
                for node, rec, node_outcome in zip(
                    observed_nodes, recs, observed_outcomes
                ):
                    self._stamp_lineage_metadata(rec, clone_idx, node)
                    if capture is not None:
                        capture.observation(captured_clone, node.id)
                    if duplex:
                        rec["parent_id"] = parent_id
                        rec["lineage_mutations_from_parent"] = node.mutations_from_parent
                    lineage_record_origins.append(
                        (parent_id if parent_id is not None else clone_idx,
                         clone_idx, node.id))
                    if expose_provenance and node_outcome is not None:
                        if capture is not None:
                            _truth_aliases(rec)
                        else:
                            _inject_truth_columns(node_outcome, self._refdata, rec)
                    if capture is not None:
                        capture.record(rec)
                    records.append(rec)
                    outcomes.append(node_outcome)
                continue

            # Corruption path: per observed node, run the library-prep /
            # sequencing corruption plan FROM the node's post-SHM
            # simulation, then merge the founder-recombination + SHM
            # provenance with the corruption trace / events so the AIRR
            # record reports trims, v/d/j, SHM counts AND artefact
            # counters (n_quality_errors, n_pcr_errors, n_indels, …).
            node_outcomes = (_native(f"clones[{clone_idx}].nodes", _eng._family_replay_node_outcomes,
                                     fam, budget=capture.budget)
                             if capture is not None else fam.node_outcomes())
            for node, base_outcome in zip(tree.nodes(), node_outcomes):
                if base_outcome is None:
                    continue
                node_seed = clone_seed + 1 + node.id
                if capture is not None:
                    initial = _native(f"clones[{clone_idx}].observations[{node.id}]",
                                      _eng._family_replay_final_simulation, base_outcome, budget=capture.budget)
                    corruption_outcome = _native(
                        f"clones[{clone_idx}].observations[{node.id}]", self._post.run_for_family_capture,
                        node_seed, strict, capture.budget, initial,
                    )
                else:
                    corruption_outcome = self._post.run_from(
                        base_outcome.final_simulation(), node_seed, strict=strict
                    )
                if capture is not None:
                    capture.observation(captured_clone, node.id, self._post, corruption_outcome)
                if capture is not None:
                    merged = _native(f"clones[{clone_idx}].observations[{node.id}]", _eng._family_replay_merge_outcomes,
                                     base_outcome, corruption_outcome, budget=capture.budget)
                    rec = _native(f"clones[{clone_idx}].observations[{node.id}]", _eng._family_replay_airr_record,
                                  merged, self._refdata, f"clone{clone_idx}_node{node.id}", budget=capture.budget)
                else:
                    merged = _eng.merge_lineage_corruption(base_outcome, corruption_outcome)
                    rec = outcome_to_airr_record(merged, self._refdata, sequence_id=f"clone{clone_idx}_node{node.id}")
                self._stamp_lineage_metadata(rec, clone_idx, node)
                if duplex:
                    rec["parent_id"] = parent_id
                    rec["lineage_mutations_from_parent"] = node.mutations_from_parent
                lineage_record_origins.append(
                    (parent_id if parent_id is not None else clone_idx,
                     clone_idx, node.id))
                if expose_provenance:
                    if capture is not None:
                        _truth_aliases(rec)
                    else:
                        _inject_truth_columns(merged, self._refdata, rec)
                if capture is not None:
                    capture.record(rec)
                records.append(rec)
                outcomes.append(merged)

        result = SimulationResultWithLineages(
            records, outcomes=outcomes, parents=parents,
            lineage_trees=lineage_trees, families=families,
            lineage_record_origins=lineage_record_origins,
            reference_content_hash=self._refdata.content_hash(),
        )
        if capture is not None:
            result._family_replay_capture = capture
            result._family_replay_refdata = self._refdata
            if not capture.request.deferred:
                from .._family_replay import _finish_capture

                _finish_capture(result)
        if validate_records:
            from .._validation import (
                _raise_on_family_validation_failure,
                _raise_on_validation_failure,
            )

            validation_kwargs = ({"_family_replay_budget": capture.budget} if capture is not None else {})
            _raise_on_validation_failure(result.validate_records(self._refdata, **validation_kwargs))
            _raise_on_family_validation_failure(result.validate_families())
            if duplex:
                _raise_on_family_validation_failure(
                    result.validate_families_with_parents(self._refdata, **validation_kwargs)
                )
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
