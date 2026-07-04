from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional

from .._compiled import (
    CompiledClonalExperiment,
    CompiledLineageExperiment,
    CompiledRepertoireExperiment,
)


class _RunMixin:
    __slots__ = ()

    def run_cohort(
        self,
        genotypes,
        *,
        n_per_subject: int = 1,
        seed: int = 0,
        counts=None,
        strict: bool = False,
        expose_provenance: bool = False,
        validate_records: bool = False,
        allow_curatable_refdata: Optional[bool] = None,
    ) -> "CohortResult":
        """Run a cohort: N subjects, each with its own diploid genotype, in one
        call. A Python loop around the single-subject genotype path — each
        subject is compiled and run independently, records are tagged with
        ``subject_id`` and given a namespaced ``sequence_id``, and the per-subject
        ``SimulationResult`` (with its own refdata) is collected into a
        :class:`~GenAIRR.cohort.CohortResult`.

        ``n_per_subject`` applies to every subject; ``counts`` (a parallel
        sequence, same length as ``genotypes``) overrides it per subject and may
        contain ``0`` (that subject appears with zero records). Subject IDs are
        taken from each genotype, or auto-assigned ``subject_0..N-1`` when all are
        unset; mixed/duplicate IDs raise. Per-subject sub-seeds are derived
        deterministically from ``seed``.

        Mutually exclusive with :meth:`with_genotype`, :meth:`restrict_alleles`,
        and ``recombine(*_allele_weights=...)`` (the genotype owns allele
        expression). :meth:`receptor_revision` is supported (each subject's
        replacement V is restricted to its carried alleles on the drawn
        chromosome); clonal forks are not supported in this release.
        """
        import copy as _copy
        import random as _random

        from ..cohort import (
            CohortResult,
            CohortSubjectResult,
            _resolve_counts,
            _resolve_subject_ids,
        )
        from ..genotype import Genotype
        from ..result import SimulationResult

        gts = list(genotypes)
        if not gts:
            raise ValueError("run_cohort: genotypes must be a non-empty sequence")
        for g in gts:
            if not isinstance(g, Genotype):
                raise TypeError(
                    f"run_cohort: every element must be a Genotype, got "
                    f"{type(g).__name__}")

        # Mutual exclusions — mirror with_genotype / compile.
        if self._genotype is not None:
            raise ValueError(
                "run_cohort() and with_genotype() are mutually exclusive")
        if any(v is not None for v in self._locks.values()):
            raise ValueError(
                "run_cohort() and restrict_alleles() are mutually exclusive")
        if self._user_allele_weights_set:
            raise ValueError(
                "run_cohort() and recombine(*_allele_weights=...) are mutually "
                "exclusive: the genotype owns allele expression")
        # receptor_revision is supported per subject (each subject's compile
        # builds a genotype-aware revision pass restricted to that subject's
        # carried alleles on the drawn chromosome) — no rejection here.
        if self._has_clonal_fork():
            raise ValueError(
                "run_cohort() is not supported together with expand_clones() / "
                "clonal_lineage() / clonal_repertoire() in this release")
        # with_metadata() is applied per subject below, but it must not overwrite
        # cohort-owned columns.
        _cohort_owned = {"subject_id", "sequence_id", "haplotype"}
        _md_collision = _cohort_owned & set(self._metadata)
        if _md_collision:
            raise ValueError(
                f"run_cohort: with_metadata keys {sorted(_md_collision)} are "
                f"cohort-owned (reserved); rename them")

        # Cartridge-hash check (same as with_genotype).
        live_hash = self._refdata.content_hash()
        for g in gts:
            if g._source_hash != live_hash:
                raise ValueError(
                    "run_cohort: a genotype was built against a different cartridge "
                    f"(content hash {g._source_hash!r} != experiment {live_hash!r})")

        resolved_counts = _resolve_counts(len(gts), n_per_subject, counts)
        subject_ids = _resolve_subject_ids([g.subject_id for g in gts])
        base_rng = _random.Random(seed)

        subjects = []
        for g, sid, count in zip(gts, subject_ids, resolved_counts):
            sub_seed = base_rng.getrandbits(63)
            snap = g._snapshot()
            snap.subject_id = sid
            # Clone the uncompiled experiment; never mutate self.
            exp_i = _copy.copy(self)
            exp_i._genotype = snap
            compiled = exp_i.compile(allow_curatable_refdata=allow_curatable_refdata)
            refdata_i = compiled.refdata
            if count == 0:
                # run() rejects n < 1; build an empty result and stamp the
                # genotype manually (run_records would otherwise do it).
                res = SimulationResult.from_outcomes(
                    [], refdata_i, expose_provenance=expose_provenance)
                res._genotypes = [snap]
            else:
                res = compiled.run_records(
                    n=count, seed=sub_seed, strict=strict,
                    expose_provenance=expose_provenance,
                    validate_records=validate_records)
            # Apply with_metadata() per subject (parity with run_records), then
            # namespace sequence_id so combined export never collides. Metadata
            # is stamped first; the cohort-owned sequence_id rewrite wins (and a
            # collision was already rejected up front).
            if self._metadata:
                for rec in res.records:
                    for key, value in self._metadata.items():
                        rec[key] = value
            for rec in res.records:
                rec["sequence_id"] = f"{sid}_{rec.get('sequence_id', '')}"
            subjects.append(CohortSubjectResult(
                subject_id=sid, genotype=snap, result=res, refdata=refdata_i,
                seed=sub_seed, count=count))

        return CohortResult(subjects)

    def run_records(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
        expose_provenance: bool = False,
        allow_curatable_refdata: Optional[bool] = None,
        validate_records: bool = False,
    ) -> "SimulationResult":
        """Compile and run, then return the batch as a
        :class:`SimulationResult` ready for ``.to_csv`` / ``.to_fasta``
        / ``.to_dataframe`` export.

        For non-clonal experiments ``n`` defaults to 1. For clonal
        experiments (when the pipeline contains :meth:`with_clonal
        _structure`) ``n`` defaults to ``n_clones * size`` and may
        be omitted; passing ``n`` explicitly is allowed only if it
        matches that product.

        ``expose_provenance=True`` appends ``truth_v_call``,
        ``truth_d_call``, ``truth_j_call`` columns containing the
        originally-sampled allele names — distinct from the
        evidence-driven ``v_call`` / ``d_call`` / ``j_call`` fields
        an aligner would produce. Useful for benchmarking aligners
        against ground truth without keeping a side truth file.

        ``strict`` semantics match :meth:`run` — strict-mode applies
        only to **fresh sampling**. Trace replay
        (:meth:`CompiledExperiment.replay_from_trace_file`) consumes
        recorded sentinel values verbatim, so a permissive trace
        replays cleanly even with ``strict=True``. See
        ``docs/productive_failure_mode_audit.md`` §5.

        ``validate_records=True`` runs
        :meth:`SimulationResult.validate_records` on the freshly
        built batch before returning. If any record fails the
        postcondition validator the call raises
        :class:`GenAIRR._validation.RecordValidationFailedError`
        (a :class:`RuntimeError` subclass) carrying a
        machine-greppable summary of the failures. The check costs
        roughly one outcome-side re-derivation per record, so it
        defaults to ``False``; flip it on in CI or when chasing a
        suspected projection bug. The validator runs **before**
        any ``with_metadata`` stamps are applied, matching the
        order :meth:`SimulationResult.validate_records` would see
        on a separate post-hoc call (metadata columns are
        per-batch annotations, not engine-derived fields).

        Returns a :class:`SimulationResult`; clonal records carry
        an integer ``clone_id`` field per row.
        """
        compiled = self.compile(allow_curatable_refdata=allow_curatable_refdata)
        if isinstance(compiled, CompiledClonalExperiment):
            result = compiled.run_records(
                n=n,
                seed=seed,
                strict=strict,
                expose_provenance=expose_provenance,
                validate_records=validate_records,
            )
        elif isinstance(compiled, CompiledRepertoireExperiment):
            if n is not None:
                raise ValueError(
                    "The 'n' parameter is not supported for clonal_repertoire "
                    "experiments. The number of records depends on the per-clone "
                    "sizes drawn from the heavy-tailed distribution and the "
                    "read-collapse, not a fixed product."
                )
            result = compiled.run_records(
                seed=seed,
                strict=strict,
                expose_provenance=expose_provenance,
                validate_records=validate_records,
            )
        elif isinstance(compiled, CompiledLineageExperiment):
            if n is not None:
                raise ValueError(
                    "The 'n' parameter is not supported for clonal_lineage experiments. "
                    "The number of observed records depends on the lineage trees "
                    "grown from n_clones / n_sample / selection, not a fixed product."
                )
            result = compiled.run_records(
                seed=seed,
                strict=strict,
                expose_provenance=expose_provenance,
                validate_records=validate_records,
            )
        else:
            result = compiled.run_records(
                n=1 if n is None else n,
                seed=seed,
                strict=strict,
                expose_provenance=expose_provenance,
                validate_records=validate_records,
            )
        if self._metadata:
            for rec in result.records:
                for key, value in self._metadata.items():
                    rec[key] = value
        return result

    def run(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
        allow_curatable_refdata: Optional[bool] = None,
    ) -> List["_engine.Outcome"]:
        """Compile and run this experiment ``n`` times.

        Equivalent to
        ``self.compile().run(n=n, seed=seed, strict=strict)``.
        Returns a list of :class:`GenAIRR._engine.Outcome` objects in
        clone-major order for clonal experiments.

        Attach :meth:`productive_only` (or any future constraint
        method) to the chain to require admissible records; the
        runtime filters NP base draws, length samples, and mutation
        / contamination substitutions in real time so the resulting
        sequences satisfy the bundle by construction.

        Statically impossible contract configurations fail during
        ``compile()`` with ``ValueError``. For runtime residue
        — i.e., empty admissible support emerging dynamically at
        sample time — ``strict=False`` (default) lets a pass consume
        the slot as its explicit no-op / sentinel; ``strict=True``
        raises :class:`GenAIRR._engine.StrictSamplingError` instead.

        Note the two error paths use **different exception classes**:
        ``ValueError`` for compile-time preconditions,
        ``StrictSamplingError`` (subclass of ``Exception``, NOT of
        ``ValueError``) for runtime empty-support. A bare
        ``except ValueError:`` will not catch the runtime case. See
        ``docs/productive_failure_mode_audit.md`` §6.1.

        ``strict`` only governs **fresh sampling**. Trace replay
        (:meth:`CompiledExperiment.replay_from_trace_file`) consumes
        recorded values verbatim; a permissive-recorded sentinel
        trace replays cleanly even with ``strict=True``. To re-execute
        a trace under strict-fresh semantics, call
        ``simulator.run(seed=<original_seed>, strict=True)`` instead.

        **Output-correctness validation is on** :meth:`run_records`
        **only.** This method returns raw ``Outcome`` objects, which
        have no projected AIRR record to validate; pass
        ``validate_records=True`` to :meth:`run_records` to opt into
        the post-build check (which raises
        :class:`GenAIRR._validation.RecordValidationFailedError` on
        any failure). For an outcome-by-outcome post-hoc check
        without re-running, build a :class:`SimulationResult` via
        :meth:`SimulationResult.from_outcomes` and call
        :meth:`SimulationResult.validate_records` on it.
        """
        compiled = self.compile(allow_curatable_refdata=allow_curatable_refdata)
        if isinstance(compiled, CompiledClonalExperiment):
            return compiled.run(n=n, seed=seed, strict=strict)
        return compiled.run(n=1 if n is None else n, seed=seed, strict=strict)

    def stream(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
    ) -> Iterator["_engine.Outcome"]:
        """Compile and lazily yield :class:`GenAIRR._engine.Outcome`
        objects. See :meth:`CompiledExperiment.stream` for full
        semantics."""
        return self.compile().stream(n=n, seed=seed, strict=strict)

    def stream_records(
        self,
        *,
        n: Optional[int] = None,
        seed: int = 0,
        strict: bool = False,
        id_prefix: str = "seq",
    ) -> Iterator[Dict[str, Any]]:
        """Compile and lazily yield AIRR-format record dicts. See
        :meth:`CompiledExperiment.stream_records`."""
        return self.compile().stream_records(
            n=n,
            seed=seed,
            strict=strict,
            id_prefix=id_prefix,
        )
