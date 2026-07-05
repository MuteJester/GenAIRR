from __future__ import annotations

import warnings
from typing import Dict, Iterable, Optional, Tuple

from GenAIRR import _engine

from .._normalize import (
    _normalize_count,
    _normalize_lengths,
    _to_immutable_byte_pair_matrix,
    _to_immutable_byte_pairs,
    _to_immutable_pairs,
)
from .._pipeline_ir import (
    _DEFAULT_NP_LENGTHS,
    _InvertDStep,
    _PairedEndStep,
    _RecombineStep,
    _ReceptorRevisionStep,
)


class _RecombinationMixin:
    __slots__ = ()

    def recombine(
        self,
        *,
        np1_lengths: Optional[Iterable[Tuple[int, float]]] = None,
        np2_lengths: Optional[Iterable[Tuple[int, float]]] = None,
        v_allele_weights: Optional[Dict[str, float]] = None,
        d_allele_weights: Optional[Dict[str, float]] = None,
        j_allele_weights: Optional[Dict[str, float]] = None,
    ) -> "Experiment":
        """Append a standard V(D)J recombination step.

        Compiles to:
        - **VJ:** sample V → sample J → (trim V_3, J_5) → assemble V
          → generate NP1 → assemble J.
        - **VDJ:** sample V → sample D → sample J → (trim V_3, D_5,
          D_3, J_5) → assemble V → generate NP1 → assemble D →
          generate NP2 → assemble J.

        ``np1_lengths`` / ``np2_lengths`` default to the species'
        empirical NP-length distributions (from
        ``DataConfig.NP_lengths``) when the experiment is bound to
        a DataConfig. For raw-RefDataConfig experiments where no
        empirical data is available, both fall back to the uniform
        ``[(0, 1.0), ..., (6, 1.0)]`` distribution and emit a
        :class:`UserWarning` so the caller knows the synthetic
        default is being used. Pass an explicit iterable of
        ``(length, weight)`` tuples to override the default.
        Passing ``np2_lengths`` on a VJ chain raises ``ValueError``
        (VJ chains have no NP2 region — there's no D segment to
        bracket).

        **Exonuclease trim** is enabled by default and uses the
        empirical per-segment trim distributions from the bound
        ``DataConfig`` when available. To disable trim, or supply
        custom trim distributions, call :meth:`trim` *after*
        :meth:`recombine` in the chain (before any mutation /
        corruption step). On a raw ``RefDataConfig`` without trim
        data, recombine emits a :class:`UserWarning` and falls back
        to a no-op trim.

        ``v_allele_weights`` / ``d_allele_weights`` /
        ``j_allele_weights`` — optional ``{allele_name:
        weight}`` dicts that bias allele sampling. Listed alleles
        get the supplied positive weight; unlisted alleles default
        to 1.0, so e.g. ``v_allele_weights={"IGHV3-23*01": 100}``
        boosts that allele while keeping every other V allele
        possible at 1/100th the rate. Mutually exclusive with the
        per-segment :meth:`restrict_alleles` restriction. Raises
        ``ValueError`` for unknown allele names or non-positive
        weights.
        """
        # Explicit allele weights conflict with an attached genotype:
        # the genotype owns allele presence + within-gene expression, and
        # the phased lowering ignores recombine-step weights. Reject the
        # combination instead of silently dropping the weights.
        if any(
            w is not None
            for w in (v_allele_weights, d_allele_weights, j_allele_weights)
        ):
            self._user_allele_weights_set = True
            if self._genotype is not None:
                raise ValueError(
                    "recombine(*_allele_weights=...) and with_genotype() are "
                    "mutually exclusive: the genotype owns allele expression"
                )

        # VJ chains have no NP2 region — surface user mistakes loudly
        # instead of silently dropping the argument.
        if np2_lengths is not None and self._refdata.chain_type != "vdj":
            raise ValueError(
                f"np2_lengths is only valid for VDJ chains; the bound "
                f"refdata is {self._refdata.chain_type!r} (no D segment, "
                f"no NP2 region). Drop the np2_lengths kwarg or bind a "
                f"VDJ refdata."
            )

        # Trim is default-on; .trim() may later disable or replace it.
        trim = True
        defaults = self._recombine_defaults() if trim or self._dataconfig else None

        # Raw-RefDataConfig path: there's no DataConfig backing this
        # experiment, so empirical NP lengths and trim distributions
        # don't exist. We fall back to a uniform NP-length default —
        # historically silent — and surface a warning so the synthetic
        # default isn't mistaken for real biology. The trim warning is
        # emitted at compile() time instead, after .trim() has had a
        # chance to disable trim explicitly.
        if self._dataconfig is None:
            if np1_lengths is None or (
                self._refdata.chain_type == "vdj" and np2_lengths is None
            ):
                warnings.warn(
                    "Experiment bound to a raw RefDataConfig with no empirical "
                    "NP-length distribution; falling back to uniform "
                    "[(0, 1.0), ..., (6, 1.0)]. Pass np1_lengths "
                    "(and np2_lengths for VDJ) explicitly to silence this.",
                    UserWarning,
                    stacklevel=2,
                )

        np1 = (
            _normalize_lengths(np1_lengths)
            if np1_lengths is not None
            else self._default_np_lengths(defaults, "np1")
        )
        np2 = (
            _normalize_lengths(np2_lengths)
            if np2_lengths is not None
            else self._default_np_lengths(defaults, "np2")
        )

        if trim and defaults is not None:
            trim_v_3 = _to_immutable_pairs(defaults.get("trim_v_3"))
            trim_d_5 = _to_immutable_pairs(defaults.get("trim_d_5"))
            trim_d_3 = _to_immutable_pairs(defaults.get("trim_d_3"))
            trim_j_5 = _to_immutable_pairs(defaults.get("trim_j_5"))
        else:
            trim_v_3 = trim_d_5 = trim_d_3 = trim_j_5 = None

        # Cartridge-owned NP base distributions (Slice — Typed NP
        # base model). Resolution lives in
        # ``_dataconfig_extract.extract_recombine_defaults``:
        # ``typed reference_models.np_bases → uniform`` (the
        # legacy auto-lift of ``NP_first_bases`` / ``NP_transitions``
        # is deliberately deferred). ``None`` here means the
        # pre-slice ``UniformBase`` applies.
        if defaults is not None:
            np1_base_pairs = _to_immutable_byte_pairs(defaults.get("np1_bases"))
            np2_base_pairs = _to_immutable_byte_pairs(defaults.get("np2_bases"))
            np1_markov_transitions = _to_immutable_byte_pair_matrix(
                defaults.get("np1_markov_transitions")
            )
            np2_markov_transitions = _to_immutable_byte_pair_matrix(
                defaults.get("np2_markov_transitions")
            )
            p_v_3_lengths = _to_immutable_pairs(
                defaults.get("p_v_3_lengths")
            )
            p_d_5_lengths = _to_immutable_pairs(
                defaults.get("p_d_5_lengths")
            )
            p_d_3_lengths = _to_immutable_pairs(
                defaults.get("p_d_3_lengths")
            )
            p_j_5_lengths = _to_immutable_pairs(
                defaults.get("p_j_5_lengths")
            )
        else:
            np1_base_pairs = None
            np2_base_pairs = None
            np1_markov_transitions = None
            np2_markov_transitions = None
            p_v_3_lengths = None
            p_d_5_lengths = None
            p_d_3_lengths = None
            p_j_5_lengths = None

        # Precedence (Slice — Allele Usage Estimation v1):
        #
        #   1. explicit `*_allele_weights=` kwarg
        #   2. cartridge `reference_models.allele_usage`
        #   3. uniform (legacy default)
        #
        # The kwarg path stays load-bearing for ad-hoc bias —
        # only when the kwarg is omitted do we fall through to
        # the typed cartridge plane (or uniform). Legacy
        # `gene_use_dict` is NOT consulted.
        if v_allele_weights is None and defaults is not None:
            v_allele_weights = defaults.get("allele_usage_v")
        if d_allele_weights is None and defaults is not None:
            d_allele_weights = defaults.get("allele_usage_d")
        if j_allele_weights is None and defaults is not None:
            j_allele_weights = defaults.get("allele_usage_j")
        weights_v = self._resolve_allele_weights("V", v_allele_weights)
        weights_d = self._resolve_allele_weights("D", d_allele_weights)
        weights_j = self._resolve_allele_weights("J", j_allele_weights)

        step = _RecombineStep(
            np1_lengths=np1,
            np2_lengths=np2,
            trim_v_3=trim_v_3,
            trim_d_5=trim_d_5,
            trim_d_3=trim_d_3,
            trim_j_5=trim_j_5,
            weights_v=weights_v,
            weights_d=weights_d,
            weights_j=weights_j,
            np1_base_pairs=np1_base_pairs,
            np2_base_pairs=np2_base_pairs,
            np1_markov_transitions=np1_markov_transitions,
            np2_markov_transitions=np2_markov_transitions,
            p_v_3_lengths=p_v_3_lengths,
            p_d_5_lengths=p_d_5_lengths,
            p_d_3_lengths=p_d_3_lengths,
            p_j_5_lengths=p_j_5_lengths,
        )
        self._steps.append(step)
        return self

    def trim(
        self,
        *,
        enabled: bool = True,
        v_3: Optional[Iterable[Tuple[int, float]]] = None,
        d_5: Optional[Iterable[Tuple[int, float]]] = None,
        d_3: Optional[Iterable[Tuple[int, float]]] = None,
        j_5: Optional[Iterable[Tuple[int, float]]] = None,
    ) -> "Experiment":
        """Configure exonuclease trim on the preceding :meth:`recombine`
        step.

        Trim is a per-segment exonuclease step that lives biologically
        *inside* V(D)J recombination (between allele selection and NP
        insertion). It is on by default with empirical distributions
        sourced from the bound ``DataConfig``. Use this method only
        when you need to override that default:

        - ``trim(enabled=False)`` — disable trim entirely. Equivalent to
          recombining against the raw allele endpoints.
        - ``trim(v_3=..., d_5=..., d_3=..., j_5=...)`` — supply custom
          per-segment trim-length distributions as ``(length, weight)``
          iterables. Omitted segments keep their empirical defaults
          (or no-op on raw RefDataConfig).

        **Position in the chain.** ``trim()`` is configuration applied
        to the most recent ``recombine()`` step. It must appear after
        ``recombine()`` and before any mutation / corruption step.
        Calling it before ``recombine()`` raises ``ValueError``;
        calling it after a mutation step raises ``ValueError``.

        Raises ``ValueError`` if no ``recombine()`` step has been
        appended yet, or if the chain already has steps that would
        biologically follow recombination (mutate, corrupt_*).
        """
        from dataclasses import replace as _replace

        # Find the most recent recombine step.
        rec_idx: Optional[int] = None
        for i in range(len(self._steps) - 1, -1, -1):
            if isinstance(self._steps[i], _RecombineStep):
                rec_idx = i
                break
        if rec_idx is None:
            raise ValueError(
                "trim() must be called after recombine(); no recombine "
                "step is on this Experiment yet."
            )

        # Anything appended after the recombine step is wrong-ordered
        # for trim configuration.
        for j in range(rec_idx + 1, len(self._steps)):
            offending = type(self._steps[j]).__name__
            raise ValueError(
                f"trim() must be called immediately after recombine() — "
                f"before any mutation / corruption / clonal-fork step. "
                f"Found a {offending!r} between the latest recombine() "
                f"and this trim() call."
            )

        prior: _RecombineStep = self._steps[rec_idx]  # type: ignore[assignment]
        if not enabled:
            # Disable all trim slots, keep NP + weights untouched.
            new_step = _replace(
                prior,
                trim_v_3=None,
                trim_d_5=None,
                trim_d_3=None,
                trim_j_5=None,
                trim_overridden=True,
            )
            self._steps[rec_idx] = new_step
            return self

        # Override individual distributions; pass-through on None.
        def _resolve(
            current: Optional[Tuple[Tuple[int, float], ...]],
            override: Optional[Iterable[Tuple[int, float]]],
        ) -> Optional[Tuple[Tuple[int, float], ...]]:
            if override is None:
                return current
            return _normalize_lengths(override)

        new_step = _replace(
            prior,
            trim_v_3=_resolve(prior.trim_v_3, v_3),
            trim_d_5=_resolve(prior.trim_d_5, d_5),
            trim_d_3=_resolve(prior.trim_d_3, d_3),
            trim_j_5=_resolve(prior.trim_j_5, j_5),
            trim_overridden=True,
        )
        self._steps[rec_idx] = new_step
        return self

    def invert_d(self, *, prob: float = 0.05) -> "Experiment":
        """Append a D-segment inversion step.

        Models V(D)J inversion: with probability ``prob`` the sampled
        D allele is committed in reverse-complement orientation
        instead of forward. Biologically, the RSS heptamers around D
        can pair head-to-head and the D segment flips before joining;
        prevalence is low (~1–5 %) but real.

        Engine path: a Bool is recorded under
        ``sample_allele.d.inverted`` and the
        :class:`~GenAIRR._engine.InvertDPass` commits
        ``ReverseComplement`` on the D :class:`AlleleInstance`
        between sampling and assembly. The
        :class:`~GenAIRR._engine.AssembleSegmentPass`(D) (Slice B)
        consumes the orientation flag and emits the reverse-
        complemented D slice into the pool. Trace replay re-fires
        the same orientation decision deterministically.

        **VDJ chains only.** VJ chains have no D pool; calling this
        method on a VJ experiment raises ``ValueError``.

        **At most once per experiment.** Calling :meth:`invert_d`
        twice raises ``ValueError`` — v1 picks a single per-pipeline
        inversion probability rather than supporting last-one-wins
        semantics (which would be silent for an over-eager builder).

        **Position in the chain.** Append after :meth:`recombine`
        (and any :meth:`trim` override). Calling :meth:`invert_d`
        before :meth:`recombine` raises ``ValueError`` at compile
        time — the lowering needs the recombine sequence already
        materialised in the engine plan so the explicit
        ``before(invert_d, assemble.d)`` schedule edge can fire.

        The DSL does **not** expose the per-record orientation in the
        AIRR record yet — that's the Slice E follow-up. End-to-end
        observability today is via the trace
        (``sample_allele.d.inverted``) and the pool bytes.

        Returns ``self`` so the call chains fluently.
        """
        if self.chain_type != "vdj":
            raise ValueError(
                f"invert_d is only valid for VDJ chains (current chain_type={self.chain_type!r})"
            )
        if any(isinstance(s, _InvertDStep) for s in self._steps):
            raise ValueError(
                "invert_d already configured on this experiment; v1 accepts at "
                "most one inversion step per pipeline. Build a fresh Experiment "
                "if you need a different probability."
            )
        # Ordering guard — D inversion is a recombination-time
        # decision and must be inherited by every clone descendant.
        # When placed post-fork, the lowering silently drops the
        # inversion probability (no recombine step in the post-fork
        # half to consume it), producing records with d_inverted=False
        # even at prob=1.0. Reject at the DSL boundary instead.
        if self._has_clonal_fork():
            raise ValueError(
                "invert_d must be called before the clonal fork; D "
                "inversion is a recombination-time decision and must "
                "be inherited by all clone descendants. Move the "
                "invert_d(...) call before clonal_lineage(...), "
                "clonal_repertoire(...), or expand_clones(...)."
            )
        if not isinstance(prob, (int, float)):
            raise ValueError(
                f"invert_d prob must be a number in [0.0, 1.0], got {type(prob).__name__}"
            )
        prob_f = float(prob)
        # NaN check FIRST — NaN fails every `<=` comparison silently
        # (`(0.0 <= nan)` is False), which would otherwise surface as
        # a misleading "out of [0.0, 1.0]" message. Explicit NaN
        # rejection gives the user the specific reason.
        if prob_f != prob_f:
            raise ValueError("invert_d prob must be a number, got NaN")
        if not (0.0 <= prob_f <= 1.0):
            raise ValueError(
                f"invert_d prob must be in [0.0, 1.0], got {prob_f}"
            )
        self._steps.append(_InvertDStep(prob=prob_f))
        return self

    def receptor_revision(self, *, prob: float = 0.05, same_haplotype: bool = True) -> "Experiment":
        """Append a receptor-revision step.

        Models post-recombination V-segment replacement: with
        probability ``prob`` the V slot is reassigned to a different
        germline V allele and the V slice in the pool is rewritten.
        Biologically, receptor revision is a B-cell tolerance
        mechanism that lets a B cell escape autoreactivity by
        replacing its V segment via secondary VDJ-recombination-like
        rearrangement on the already-assembled receptor.

        Engine path: a Bool is recorded under
        ``receptor_revision.applied`` for every simulation. On
        ``true``, the
        :class:`~GenAIRR._engine.ReceptorRevisionPass` additionally
        records the replacement allele id at
        ``receptor_revision.v_allele`` and the derived 3' trim at
        ``receptor_revision.v_trim_3``, then commits
        ``AssignmentChanged`` + ``TrimChanged`` + ``SegmentReplaced``
        against V through a single
        :class:`~GenAIRR._engine.SimulationBuilder`. Slice C's
        same-length retained constraint
        (``allele.len() - trim_3 == old_v_len``, 5' trim fixed at 0)
        keeps downstream pool positions stable; the
        :class:`~GenAIRR._engine.LiveCallRefreshHook` (Slice B)
        reacts to ``SegmentReplaced`` with an AllStructural-
        equivalent V/D/J re-walk.

        **VDJ chains only.** Receptor revision is heavy-chain v1;
        calling this method on a VJ experiment raises
        ``ValueError``.

        **At most once per experiment.** Calling
        :meth:`receptor_revision` twice raises ``ValueError`` —
        last-one-wins semantics would silently override an
        over-eager builder.

        **Position in the chain.** Appended after :meth:`recombine`;
        the lowering inlines the engine ``push_receptor_revision``
        call immediately after ``push_assemble("J")`` so the pass
        sees the fully-assembled V/D/J/NP pool. Subsequent
        :meth:`mutate` / corruption passes lower after this step in
        the plan, giving the canonical "recombine → revise →
        mutate/corrupt" order the design doc §2 requires.

        AIRR records expose ``receptor_revision_applied`` (bool) and
        ``original_v_call`` (the pre-revision V; empty when no
        revision applied). ``v_call`` / ``truth_v_call`` are the
        **post**-revision V. The three ``receptor_revision.*`` trace
        records above remain available for replay.

        Returns ``self`` so the call chains fluently.
        """
        if self.chain_type != "vdj":
            raise ValueError(
                "receptor_revision is only valid for VDJ chains "
                f"(current chain_type={self.chain_type!r})"
            )
        if any(isinstance(s, _ReceptorRevisionStep) for s in self._steps):
            raise ValueError(
                "receptor_revision already configured on this experiment; "
                "v1 accepts at most one revision step per pipeline. Build "
                "a fresh Experiment if you need a different probability."
            )
        # Ordering guard — receptor revision is a recombination/
        # ancestor-time decision and must be inherited by every clone
        # descendant. When placed post-fork, the lowering silently
        # drops the revision probability (no recombine step in the
        # post-fork half to consume it), producing records with
        # receptor_revision_applied=False and empty original_v_call
        # even at prob=1.0. Reject at the DSL boundary instead.
        if self._has_clonal_fork():
            raise ValueError(
                "receptor_revision must be called before "
                "the clonal fork; receptor revision is a "
                "recombination-time decision and must be inherited by "
                "all clone descendants. Move the "
                "receptor_revision(...) call before clonal_lineage(...), "
                "clonal_repertoire(...), or expand_clones(...)."
            )
        if not isinstance(prob, (int, float)):
            raise ValueError(
                "receptor_revision prob must be a number in [0.0, 1.0], "
                f"got {type(prob).__name__}"
            )
        prob_f = float(prob)
        # NaN first — see the matching `invert_d` rationale.
        if prob_f != prob_f:
            raise ValueError("receptor_revision prob must be a number, got NaN")
        if not (0.0 <= prob_f <= 1.0):
            raise ValueError(
                f"receptor_revision prob must be in [0.0, 1.0], got {prob_f}"
            )
        if not isinstance(same_haplotype, bool):
            raise ValueError(
                "receptor_revision same_haplotype must be a bool, got "
                f"{same_haplotype!r}"
            )
        self._steps.append(
            _ReceptorRevisionStep(prob=prob_f, same_haplotype=same_haplotype)
        )
        return self

    def paired_end(
        self,
        *,
        r1_length,
        r2_length=None,
        insert_size,
    ) -> "Experiment":
        """Append a paired-end / read-layout step.

        Models the Illumina paired-end read layout: each fragment
        produces R1 (forward from the 5' adapter) and R2
        (reverse-complemented from the 3' adapter) windows over the
        final projected molecule, plus an *insert size* that locates
        R2's 3' end. The DSL exposes three integer distributions:

        - ``r1_length`` — required.
        - ``r2_length`` — defaults to ``r1_length`` when ``None``.
          Many Illumina libraries do run asymmetric (R2 quality
          drops faster); the explicit shape lets callers opt in.
        - ``insert_size`` — required.

        Each accepts the same three shapes the rest of the DSL
        already uses for length-like distributions:

        - ``int`` — fixed value.
        - ``(low, high)`` — uniform integer in the closed
          interval ``[low, high]``.
        - ``[(value, weight), …]`` — explicit empirical
          distribution.

        Engine path: a trace-only
        :class:`~GenAIRR._engine.PairedEndSamplingPass` records
        three Ints at ``paired_end.r1_length`` /
        ``paired_end.r2_length`` / ``paired_end.insert_size``;
        the AIRR builder reads them back at projection time and
        populates the eight ``read_layout`` / ``r1_sequence`` /
        ``r2_sequence`` / ``r1_start`` / ``r1_end`` / ``r2_start`` /
        ``r2_end`` / ``insert_size`` fields via the Slice B
        projection kernel. ``rec.sequence`` is the only
        coordinate space — end-loss and rev-comp projections have
        already finalised the molecule by the time paired-end
        windows are drawn (design doc §6 / §7).

        **Both VDJ and VJ chains supported.** Paired-end is a
        sequencing-stage observable, not a biology mechanism;
        it makes sense on every chain.

        **At most once per experiment.** Calling
        :meth:`paired_end` twice raises ``ValueError`` —
        last-one-wins semantics would silently override an
        over-eager builder.

        **Position in the chain.** The compile pre-pass extracts
        the step and pushes the engine pass at the **end** of the
        plan, after every IR-mutating / corruption / orientation
        step. Even though the pass is trace-only, recording the
        choices last keeps the trace order aligned with the
        biological/readout order (recombine → mutation →
        corruption → end-loss → paired-end).

        Returns ``self`` so the call chains fluently.
        """
        if any(isinstance(s, _PairedEndStep) for s in self._steps):
            raise ValueError(
                "paired_end already configured on this experiment; "
                "v1 accepts at most one paired-end step per pipeline. "
                "Build a fresh Experiment if you need a different "
                "layout."
            )

        # Resolve default r2 → r1 BEFORE normalization so the
        # downstream check pins both at the same source shape.
        if r2_length is None:
            r2_length = r1_length

        r1_pairs = _normalize_count(r1_length)
        r2_pairs = _normalize_count(r2_length)
        insert_pairs = _normalize_count(insert_size)

        # r1 / r2 lengths must be strictly positive. `_normalize_count`
        # already rejects negatives but allows 0; the projection
        # kernel needs > 0, so surface the violation at the DSL
        # boundary with a clearer message.
        for value, _w in r1_pairs:
            if value <= 0:
                raise ValueError(
                    f"paired_end r1_length must be positive, got {value}"
                )
        for value, _w in r2_pairs:
            if value <= 0:
                raise ValueError(
                    f"paired_end r2_length must be positive, got {value}"
                )
        # `_normalize_count` already enforces insert_size >= 0;
        # the projection kernel enforces the same.

        # Fixed-value geometry checks. When r1/r2/insert are each
        # single-value distributions we know the exact values the
        # pass will emit, so reject `r1 > insert` / `r2 > insert`
        # here rather than waiting for the engine's per-sample
        # `InvalidDistributionOutput`. Distribution / range cases
        # may sample valid combinations, so we defer those to the
        # engine.
        if len(r1_pairs) == 1 and len(insert_pairs) == 1:
            r1_value = r1_pairs[0][0]
            insert_value = insert_pairs[0][0]
            if r1_value > insert_value:
                raise ValueError(
                    f"paired_end r1_length ({r1_value}) > insert_size "
                    f"({insert_value}); the R1 window would run past "
                    f"the fragment 3' end."
                )
        if len(r2_pairs) == 1 and len(insert_pairs) == 1:
            r2_value = r2_pairs[0][0]
            insert_value = insert_pairs[0][0]
            if r2_value > insert_value:
                raise ValueError(
                    f"paired_end r2_length ({r2_value}) > insert_size "
                    f"({insert_value}); the R2 window would run past "
                    f"the fragment 3' end."
                )

        self._steps.append(
            _PairedEndStep(
                r1_length=r1_pairs,
                r2_length=r2_pairs,
                insert_size=insert_pairs,
            )
        )
        return self

    def _recombine_defaults(self):
        """Lazy-extract the empirical distributions for this
        experiment's DataConfig. Returns ``None`` for raw-RefDataConfig
        experiments. Cached on first call.
        """
        if self._dataconfig is None:
            return None
        from .._dataconfig_extract import extract_recombine_defaults

        return extract_recombine_defaults(self._dataconfig)

    def _build_contracts(self) -> Optional["_engine.ContractSet"]:
        """Synthesize the engine ``ContractSet`` from declared bundles.

        Today only the productive bundle is recognized. Future bundles
        (e.g. ``.in_frame_only()``) would compose here. Returns
        ``None`` when no constraint methods have been called.
        """
        if not self._contracts:
            return None
        # Composition across bundles isn't supported by the engine yet,
        # so for now exactly one bundle is allowed.
        if len(self._contracts) > 1:
            raise NotImplementedError(
                f"composing multiple constraint bundles is not yet supported; "
                f"got {self._contracts!r}"
            )
        bundle = self._contracts[0]
        if bundle == "productive":
            return _engine.productive()
        raise NotImplementedError(
            f"unknown constraint bundle {bundle!r}"
        )

    @staticmethod
    def _default_np_lengths(
        defaults, key: str
    ) -> Tuple[Tuple[int, float], ...]:
        """Pick the NP-length distribution to use when the user
        didn't pass an explicit one: empirical if available, else
        the uniform placeholder.
        """
        if defaults is not None:
            empirical = defaults.get(key)
            if empirical:
                return tuple(empirical)
        return tuple(_DEFAULT_NP_LENGTHS)
