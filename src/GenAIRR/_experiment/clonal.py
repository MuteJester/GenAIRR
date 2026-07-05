from __future__ import annotations

import warnings
from typing import Optional

from .._pipeline_ir import (
    _ClonalForkStep,
    _LineageForkStep,
    _RepertoireForkStep,
)
from .._step_validation import _descendant_phase_step_classifier


class _ClonalMixin:
    __slots__ = ()

    def expand_clones(
        self,
        *,
        n_clones: int,
        per_clone: int,
    ) -> "Experiment":
        """Expand the pipeline into clonal lineages.

        Marks the per-clone / per-descendant boundary in the chain:
        steps appended *before* this call run **once per clone** —
        typically just :meth:`recombine`, which establishes the
        parent V/D/J + trim + NP + assembled IR for the clonal
        family. Steps appended *after* this call run **once per
        read inside the family** — typically :meth:`mutate` and the
        library-prep / sequencing-stage steps, which introduce
        per-read divergence within the clone.

        Concrete shape::

            exp = (Experiment.on("human_igh")
                   .recombine()
                   .expand_clones(n_clones=10, per_clone=20)
                   .mutate(rate=0.05)
                   .pcr_amplify(count=2))
            result = exp.run_records(seed=0)
            # 10 clones × 20 descendants = 200 records.
            # Each record carries a ``clone_id`` integer in [0, 10).

        ``n`` can be omitted from :meth:`run_records` for a clonal
        experiment — the runtime expands ``n_clones * per_clone``
        records automatically. Passing ``n`` is allowed only when
        ``n == n_clones * per_clone``.

        .. deprecated::
            Use :meth:`clonal_lineage` for BCR affinity-maturation
            trees, or :meth:`clonal_repertoire` for TCR / flat-BCR
            abundance repertoires with clone-size distributions.
            ``expand_clones`` remains supported for fixed-size flat
            star expansion.

        Constraints:
        - Both ``n_clones`` and ``per_clone`` must be positive ints.
        - At most one expansion per pipeline; calling this method
          twice raises ``ValueError``.

        Implementation note: the runtime forks the parent's IR
        (final ``Simulation`` after the pre-fork plan) into
        descendants by running the post-fork plan from that IR
        with distinct seeds. Within a clone, every descendant
        shares the same recombination provenance (V allele, trim,
        NP bases) and only diverges through the post-fork passes.
        """
        warnings.warn(
            "Experiment.expand_clones() is deprecated. Use "
            "Experiment.clonal_lineage() for BCR affinity-maturation trees "
            "(internal SHM, lineage_trees, no per_clone) or "
            "Experiment.clonal_repertoire() for TCR / flat-BCR abundance "
            "repertoires with clone-size distributions and duplicate_count. "
            "Neither is a drop-in replacement for fixed per_clone output; "
            "expand_clones() remains supported for legacy fixed-size flat "
            "star expansion.",
            DeprecationWarning,
            stacklevel=2,
        )
        if not isinstance(n_clones, int) or isinstance(n_clones, bool) or n_clones < 1:
            raise ValueError(
                f"n_clones must be a positive int, got {n_clones!r}"
            )
        if not isinstance(per_clone, int) or isinstance(per_clone, bool) or per_clone < 1:
            raise ValueError(f"per_clone must be a positive int, got {per_clone!r}")
        if self._has_clonal_fork():
            raise ValueError(
                "expand_clones() / clonal_lineage() / clonal_repertoire() "
                "can only be called once per pipeline"
            )
        # Unified descendant-phase ordering guard. Scan the appended
        # step list for any step that came from a descendant-phase
        # DSL method (see ``_descendant_phase_step_classifier``).
        # Pre-fork placement of these steps either silently misreports
        # the AIRR field (trace-sourced fields don't survive the
        # parent→descendant boundary — Bugs C / E / F) or collapses
        # descendant diversity (the pass runs once on the parent IR
        # and every clone member inherits an identical effect).
        # Either failure mode is a clonal-semantics violation, so
        # reject at the DSL boundary with a uniform message that
        # names the offending step and the fix.
        for step in self._steps:
            offending_method = _descendant_phase_step_classifier(step)
            if offending_method is None:
                continue
            if offending_method == "mutate":
                detail = (
                    "SHM is descendant-specific in GenAIRR's current "
                    "clonal model"
                )
            else:
                detail = (
                    "it is descendant-specific and must be sampled "
                    "independently for each clone member"
                )
            raise ValueError(
                f"{offending_method} must be called after "
                f"expand_clones(); {detail}. Move "
                f"{offending_method}(...) after expand_clones(...)."
            )
        self._steps.append(_ClonalForkStep(n_clones=n_clones, size=per_clone))
        return self

    def clonal_repertoire(
        self,
        *,
        n_clones: int,
        size_distribution: str = "power_law",
        exponent: float = 2.0,
        mu: float = 1.0,
        sigma: float = 1.0,
        max_size: int = 1000,
        unexpanded_fraction: float = 0.0,
    ) -> "Experiment":
        """Expand the pipeline into a non-tree clonal **repertoire**.

        This is the non-tree clonal model (contrast with
        :meth:`clonal_lineage`, which grows real affinity-maturation
        lineage *trees*). It generalizes the deprecated
        :meth:`expand_clones`: instead of a fixed ``per_clone`` size,
        each clone draws a size from a heavy-tailed distribution (with
        an unexpanded-singleton fraction). That many reads pass through
        the post-fork library-prep / sequencing passes, and identical
        reads are genotype-collapsed into AIRR records carrying a
        standard ``duplicate_count`` that records abundance.

        Like :meth:`expand_clones`, this call marks the per-clone /
        per-read boundary: steps appended *before* it run **once per
        clone** (typically just :meth:`recombine`, establishing the
        clonal V/D/J + trim + NP backbone); steps appended *after* it
        run **once per read** (the library-prep / sequencing passes
        that introduce per-read divergence).

        Concrete shape::

            exp = (Experiment.on("human_igh")
                   .recombine()
                   .clonal_repertoire(n_clones=200, max_size=500,
                                      unexpanded_fraction=0.3)
                   .sequencing_errors(rate=0.005))
            result = exp.run_records(seed=0)
            # Each record carries a `clone_id` and a `duplicate_count`.

        Parameters:
        - ``n_clones`` — number of clones (positive int).
        - ``size_distribution`` — ``"power_law"`` or ``"lognormal"``.
        - ``exponent`` — power-law exponent (``> 0``; used when
          ``size_distribution="power_law"``).
        - ``mu`` / ``sigma`` — lognormal parameters (``sigma >= 0``;
          used when ``size_distribution="lognormal"``).
        - ``max_size`` — clamps the largest clone. Because the total
          number of reads simulated is roughly the **sum** of the drawn
          sizes when post-fork passes are present, keep ``max_size``
          modest to bound runtime.
        - ``unexpanded_fraction`` — fraction of clones forced to size 1
          (unexpanded singletons), in ``[0, 1]``.

        TCR works out of the box (no SHM): the reads diverge only
        through the post-fork sequencing passes. Note that ``mutate``
        after this call is rejected on TCR by :meth:`mutate`'s own TCR
        guard; on BCR an optional post-fork ``mutate`` adds flat SHM.

        **No-corruption shortcut:** a clone with no post-fork passes
        emits identical copies, so it collapses to a single record whose
        ``duplicate_count`` equals the drawn size.

        Constraints:
        - At most one fork per pipeline — calling this when an
          :meth:`expand_clones`, :meth:`clonal_lineage`, or
          :meth:`clonal_repertoire` fork is already present raises
          ``ValueError``.
        - The same descendant-phase ordering guard as
          :meth:`expand_clones` applies: a descendant-phase step (e.g.
          ``.mutate()``) appended *before* this call is rejected;
          :meth:`recombine` is fine.
        """
        if not isinstance(n_clones, int) or isinstance(n_clones, bool) or n_clones < 1:
            raise ValueError(
                f"n_clones must be a positive int, got {n_clones!r}"
            )
        if size_distribution not in ("power_law", "lognormal"):
            raise ValueError(
                "size_distribution must be 'power_law' or 'lognormal', got "
                f"{size_distribution!r}"
            )
        if not (isinstance(exponent, (int, float)) and exponent > 0):
            raise ValueError(f"exponent must be > 0, got {exponent!r}")
        if not (isinstance(sigma, (int, float)) and sigma >= 0):
            raise ValueError(f"sigma must be >= 0, got {sigma!r}")
        if not isinstance(max_size, int) or isinstance(max_size, bool) or max_size < 1:
            raise ValueError(f"max_size must be a positive int, got {max_size!r}")
        if not (
            isinstance(unexpanded_fraction, (int, float))
            and 0.0 <= unexpanded_fraction <= 1.0
        ):
            raise ValueError(
                "unexpanded_fraction must be in [0, 1], got "
                f"{unexpanded_fraction!r}"
            )
        if any(
            isinstance(s, (_ClonalForkStep, _RepertoireForkStep, _LineageForkStep))
            for s in self._steps
        ):
            raise ValueError(
                "clonal_repertoire() / expand_clones() / clonal_lineage() "
                "can only be called once per pipeline"
            )
        # Same descendant-phase ordering guard as expand_clones: a
        # descendant-phase step (mutate / corruption / paired_end)
        # appended before the fork is rejected.
        for step in self._steps:
            offending_method = _descendant_phase_step_classifier(step)
            if offending_method is None:
                continue
            if offending_method == "mutate":
                detail = (
                    "SHM is descendant-specific in GenAIRR's current "
                    "clonal model"
                )
            else:
                detail = (
                    "it is descendant-specific and must be sampled "
                    "independently for each read"
                )
            raise ValueError(
                f"{offending_method} must be called after "
                f"clonal_repertoire(); {detail}. Move "
                f"{offending_method}(...) after clonal_repertoire(...)."
            )
        self._steps.append(
            _RepertoireForkStep(
                n_clones=n_clones,
                size_distribution=size_distribution,
                exponent=float(exponent),
                mu=float(mu),
                sigma=float(sigma),
                max_size=max_size,
                unexpanded_fraction=float(unexpanded_fraction),
            )
        )
        return self

    def clonal_lineage(
        self,
        *,
        n_clones: int,
        max_generations: int = 10,
        n_max: int = 1000,
        n_sample: int = 50,
        rate: float = 0.05,
        lambda_base: float = 1.5,
        selection_strength: float = 0.0,
        beta: float = 1.0,
        target_aa: Optional[str] = None,
        mature_substitutions: int = 5,
        s5f_model: str = "hh_s5f",
        allow_extinction: bool = False,
    ) -> "Experiment":
        """Grow BCR lineage trees (neutral by default; set ``selection_strength > 0``
        and optionally ``target_aa`` to enable affinity maturation).

        Each clone gets its own lineage tree produced by the Rust
        ``simulate_family_outcomes`` kernel. The returned
        :class:`~GenAIRR.result.SimulationResultWithLineages` carries:

        - ``.records`` — one AIRR dict per *observed* (genotype-collapsed) cell,
          tagged with ``clone_id``, ``lineage_node_id``, ``lineage_parent_id``,
          ``lineage_generation``, ``lineage_abundance``, and
          ``lineage_affinity``. Because identical genotypes are collapsed before
          sampling, the number of records per clone is ≤ ``n_sample``; the
          ``lineage_abundance`` field accounts for how many sampled cells were
          represented by each observed record. Mutation counts (``n_mutations``,
          ``n_v_mutations``, …) are pool-derived and self-consistent.
        - ``.lineage_trees`` — one :class:`~GenAIRR._engine.LineageTree`
          per clone for ground-truth export (Newick, FASTA, node table TSV).

        Parameters
        ----------
        n_clones:
            Number of independent clonal lineages to grow.
        max_generations:
            Maximum depth of the lineage tree (≤ 1000).
        n_max:
            Per-generation LIVING-population carrying capacity: the live
            population per generation is capped at this (the tree can contain
            more total nodes across generations). It is NOT a hard cap on the
            total number of cells per clone.
        n_sample:
            Number of cells to sample as observed leaves. Records returned
            per clone are ≤ ``n_sample`` because identical genotypes are
            collapsed (duplicates are counted in ``lineage_abundance``).
        rate:
            Per-base SHM rate for within-lineage mutations.
        lambda_base:
            Poisson mean for offspring count at affinity 0.
        selection_strength:
            Selection pressure; ``0.0`` = neutral drift (fitness is 1.0 for
            every cell). This disables selection but does not force
            ``lineage_affinity`` to 0 when a target sequence is supplied.
            Set ``> 0`` to enable affinity maturation; combine with
            ``target_aa`` for a fixed sequence target.
        beta:
            Scaling factor for the affinity term in ``exp(−beta·distance)``.
        target_aa:
            Amino-acid sequence of the full receptor used to compute
            per-cell affinity via a BLOSUM62-weighted distance (compared
            position-wise against the cell's translated receptor; only
            the overlapping prefix is scored when lengths differ). Must be
            a non-empty string of standard amino-acid letters
            (``ACDEFGHIKLMNPQRSTVWY``). When ``None``, an auto target is
            generated from the founder by applying ``mature_substitutions``
            random residue changes whenever selection is enabled. In fully
            neutral mode (``selection_strength=0`` and ``target_aa=None``),
            no affinity model is built and ``lineage_affinity`` is 0.
        mature_substitutions:
            Number of amino-acid substitutions used to build the auto
            target (when ``target_aa`` is ``None``).
        s5f_model:
            Bundled S5F kernel name for within-lineage mutation context
            (``"hh_s5f"``, ``"hkl_s5f"``, …).
        allow_extinction:
            Sampling draws from the LIVING final-generation population, so a
            founder that draws 0 offspring goes extinct and yields zero
            observed cells/records. With ``allow_extinction=False`` (default)
            each requested clone is conditioned on survival: an extinct family
            is retried with a fresh deterministic sub-seed (up to a bounded
            number of attempts) so you reliably get ``n_clones`` families. With
            ``allow_extinction=True`` extinction is accepted and the extinct
            clone is skipped, producing fewer families than ``n_clones``.

        **BCR-only guard:** ``clonal_lineage`` applies S5F somatic
        hypermutation, which is a B-cell process. Calling it on a TCR-configured
        experiment raises ``ValueError`` (immunoglobulin / BCR loci only). TCR
        clone-size primitives exist in the engine but are not yet exposed as a
        DSL workflow.
        """
        import math
        import warnings

        # --- BCR-only guard (mirror mutate()'s TCR rejection) ---
        # clonal_lineage applies S5F somatic hypermutation, a B-cell process,
        # so it must reject TCR loci. ``_is_tcr_refdata`` inspects the first V
        # allele name prefix (TR* => TCR, IG* => BCR) on the already-bound
        # refdata, exactly as mutate() does. Firing here, at call time and
        # before compile(), guarantees the clear BCR-only message instead of a
        # downstream cartridge / compile error.
        if self._is_tcr_refdata():
            locus = self._refdata.v_allele(0).name if self._refdata.v_pool_size() else "?"
            raise ValueError(
                "clonal_lineage models B-cell somatic hypermutation and "
                "supports immunoglobulin (BCR) loci only; the locus "
                f"'{locus}' is a TCR locus. (TCR clone-size simulation is not "
                "yet exposed in the DSL.)"
            )
        # --- allow_extinction ---
        if not isinstance(allow_extinction, bool):
            raise ValueError(
                f"allow_extinction must be a bool, got {allow_extinction!r}"
            )

        # --- n_clones ---
        if isinstance(n_clones, bool) or not isinstance(n_clones, int) or n_clones < 1:
            raise ValueError(f"n_clones must be a positive int, got {n_clones!r}")
        # --- max_generations ---
        if (
            isinstance(max_generations, bool)
            or not isinstance(max_generations, int)
            or max_generations < 1
            or max_generations > 1000
        ):
            raise ValueError(
                f"max_generations must be a positive int <= 1000, got {max_generations!r}"
            )
        # --- n_max ---
        if isinstance(n_max, bool) or not isinstance(n_max, int) or n_max < 1:
            raise ValueError(f"n_max must be a positive int, got {n_max!r}")
        # --- n_sample ---
        if isinstance(n_sample, bool) or not isinstance(n_sample, int) or n_sample < 1:
            raise ValueError(f"n_sample must be a positive int, got {n_sample!r}")
        # --- rate ---
        if not isinstance(rate, (int, float)) or not math.isfinite(rate) or rate < 0 or rate > 1:
            raise ValueError(f"rate must be a float in [0, 1], got {rate!r}")
        # --- lambda_base ---
        if not isinstance(lambda_base, (int, float)) or not math.isfinite(lambda_base) or lambda_base < 0:
            raise ValueError(f"lambda_base must be a finite non-negative float, got {lambda_base!r}")
        # --- beta ---
        if not isinstance(beta, (int, float)) or not math.isfinite(beta) or beta < 0:
            raise ValueError(f"beta must be a finite non-negative float, got {beta!r}")
        # --- selection_strength ---
        if not isinstance(selection_strength, (int, float)) or not math.isfinite(selection_strength) or selection_strength < 0:
            raise ValueError(
                f"selection_strength must be a finite non-negative float, got {selection_strength!r}"
            )
        # --- mature_substitutions ---
        if (
            isinstance(mature_substitutions, bool)
            or not isinstance(mature_substitutions, int)
            or mature_substitutions < 0
        ):
            raise ValueError(
                f"mature_substitutions must be a non-negative int, got {mature_substitutions!r}"
            )
        # --- target_aa ---
        _VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
        if target_aa is not None:
            if not isinstance(target_aa, str) or len(target_aa) == 0:
                raise ValueError(
                    "target_aa must be a non-empty amino-acid string "
                    "(letters from ACDEFGHIKLMNPQRSTVWY)"
                )
            target_aa = target_aa.upper()
            invalid = set(target_aa) - _VALID_AA
            if invalid:
                raise ValueError(
                    f"target_aa contains invalid characters {sorted(invalid)!r}; "
                    "only standard amino-acid letters (ACDEFGHIKLMNPQRSTVWY) are allowed"
                )
            if len(target_aa) < 30:
                warnings.warn(
                    f"target_aa has length {len(target_aa)}, which is shorter than a typical "
                    "receptor sequence (~300+ aa). If this is an epitope sequence rather than "
                    "the full receptor, affinity scoring will be based only on the overlapping "
                    "prefix — consider supplying the full translated receptor instead.",
                    UserWarning,
                    stacklevel=2,
                )
        # --- s5f_model (validate at call time, not at run time) ---
        from GenAIRR._s5f_loader import _BUILTIN_S5F_MODELS
        _s5f_key = s5f_model.lower().strip()
        if _s5f_key not in _BUILTIN_S5F_MODELS:
            avail = ", ".join(f'"{k}"' for k in sorted(_BUILTIN_S5F_MODELS))
            raise ValueError(
                f"Unknown s5f_model {s5f_model!r}. Available: {avail}"
            )
        # --- reject duplicate fork steps ---
        for s in self._steps:
            if isinstance(s, (_ClonalForkStep, _RepertoireForkStep, _LineageForkStep)):
                raise ValueError(
                    "clonal_lineage() / expand_clones() / clonal_repertoire() "
                    "can only be called once per pipeline"
                )
        # --- descendant-phase guard (same as expand_clones) ---
        for step in self._steps:
            offending_method = _descendant_phase_step_classifier(step)
            if offending_method is None:
                continue
            raise ValueError(
                f"{offending_method} must be called after "
                f"clonal_lineage(); it is descendant-specific and must be sampled "
                f"independently for each clone member. Move "
                f"{offending_method}(...) after clonal_lineage(...)."
            )
        self._steps.append(
            _LineageForkStep(
                n_clones=n_clones,
                max_generations=max_generations,
                n_max=n_max,
                n_sample=n_sample,
                rate=rate,
                lambda_base=lambda_base,
                selection_strength=selection_strength,
                beta=beta,
                target_aa=target_aa,
                mature_substitutions=mature_substitutions,
                s5f_model=s5f_model,
                allow_extinction=allow_extinction,
            )
        )
        return self
