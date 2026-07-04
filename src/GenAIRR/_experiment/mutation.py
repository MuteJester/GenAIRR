from __future__ import annotations

from typing import Dict, Iterable, Optional, Tuple, Union

from .._normalize import _normalize_count
from .._pipeline_ir import (
    _ClonalForkStep,
    _LineageForkStep,
    _MutateStep,
    _RepertoireForkStep,
)
from .._step_validation import (
    _DEFAULT_V_SUBREGION_RATES,
    _validate_segment_rates,
    _validate_v_subregion_rates,
)


class _MutationMixin:
    __slots__ = ()

    def mutate(
        self,
        *,
        model: str = "s5f",
        count: Optional[Union[int, Tuple[int, int], Iterable[Tuple[int, float]]]] = None,
        rate: Optional[float] = None,
        s5f_model: str = "hh_s5f",
        segment_rates: Optional[Dict[str, float]] = None,
        v_subregion_rates: Optional[Dict[str, float]] = None,
    ) -> "Experiment":
        """Append a somatic-hypermutation step.

        ``model`` selects the mutation kernel:
        - ``"s5f"`` (default) — context-dependent SHM via the bundled
          S5F kernel named in ``s5f_model``. Available kernels:
          ``"hh_s5f"``, ``"hh_s5f_60"``, ``"hh_s5f_opposite"``,
          ``"hkl_s5f"``.
        - ``"uniform"`` — position-independent SHM. Each mutated
          position gets a uniformly drawn A/C/G/T replacement.

        **Specify intensity with exactly one of ``rate`` or ``count``.**

        ``rate`` is the per-base mutation rate (e.g. ``0.03`` for 3 %
        SHM, which is roughly memory B-cell SHM). At execute time the
        engine draws ``count ~ Poisson(rate × pool_len)`` against each
        record's current sequence length — so the realized count
        scales with each record's actual length, matching how
        immunologists report SHM in the literature. This is the
        canonical, biology-default form.

        ``count`` is the legacy explicit count distribution, useful
        for benchmark scripts that want a deterministic count
        independent of record length:
        - ``count=15`` — fixed: every simulation gets exactly 15
          mutations.
        - ``count=(5, 25)`` — uniform integer in ``[5, 25]`` (both
          endpoints inclusive).
        - ``count=[(5, 1.0), (10, 2.0), ...]`` — explicit empirical
          ``(count, weight)`` distribution.

        Passing both ``count`` and ``rate`` raises ``ValueError``.
        Passing neither raises ``ValueError``.

        **TCR guard:** somatic hypermutation is a B-cell
        phenomenon — T-cells do not undergo SHM in the periphery.
        Calling ``.mutate()`` on a TCR-configured experiment raises
        ``ValueError`` to prevent silent biological misuse. Use
        ``pcr_amplify`` / ``sequencing_errors`` for sequencing-error
        realism on TCR data instead.
        """
        if self._is_tcr_refdata():
            raise ValueError(
                "mutate(): somatic hypermutation does not occur in TCR "
                "sequences (T-cells lack AID and the SHM machinery). The "
                "configured refdata is a TCR locus. For sequencing-error "
                "realism on TCR data, use pcr_amplify / sequencing_errors "
                "/ polymerase_indels instead."
            )
        model_lc = model.lower()
        if model_lc not in ("uniform", "s5f"):
            raise ValueError(
                f"model must be 'uniform' or 's5f' (got {model!r})"
            )
        if count is not None and rate is not None:
            raise ValueError(
                "mutate(): pass exactly one of `rate` or `count`, not both. "
                "`rate` is the canonical biology default (e.g. rate=0.03 "
                "for 3% SHM); `count` is the explicit per-record count "
                "for benchmark / deterministic-count workflows."
            )
        if count is None and rate is None:
            raise ValueError(
                "mutate(): pass exactly one of `rate` or `count`. "
                "Suggested default: rate=0.03 (~3% SHM, memory B-cell range)."
            )
        if rate is not None:
            if not isinstance(rate, (int, float)) or isinstance(rate, bool):
                raise TypeError(
                    f"rate must be a finite float in [0.0, 1.0], got "
                    f"{type(rate).__name__}"
                )
            if not (0.0 <= float(rate) <= 1.0):
                raise ValueError(
                    f"rate must be in [0.0, 1.0] (got {rate!r}); rate is a "
                    f"per-base mutation probability, not an absolute count."
                )
            seg_rates_tuple = _validate_segment_rates(segment_rates)
            v_sub_rates_tuple = _validate_v_subregion_rates(v_subregion_rates)
            self._check_v_subregion_rates_satisfiable(
                v_subregion_rates, v_sub_rates_tuple
            )
            self._steps.append(
                _MutateStep(
                    model=model_lc,
                    s5f_model_name=s5f_model,
                    rate=float(rate),
                    segment_rates=seg_rates_tuple,
                    v_subregion_rates=v_sub_rates_tuple,
                )
            )
            return self
        pairs = _normalize_count(count)
        seg_rates_tuple = _validate_segment_rates(segment_rates)
        v_sub_rates_tuple = _validate_v_subregion_rates(v_subregion_rates)
        self._check_v_subregion_rates_satisfiable(
            v_subregion_rates, v_sub_rates_tuple
        )
        self._steps.append(
            _MutateStep(
                model=model_lc,
                s5f_model_name=s5f_model,
                count_pairs=pairs,
                segment_rates=seg_rates_tuple,
                v_subregion_rates=v_sub_rates_tuple,
            )
        )
        return self

    def _check_v_subregion_rates_satisfiable(
        self,
        raw_rates: Optional[Dict[str, float]],
        tuple_rates: Tuple[float, float, float, float, float],
    ) -> None:
        """Reject a non-default ``v_subregion_rates`` configuration
        when the bound cartridge has zero annotated V alleles.
        Audit §4: a non-default rate vector against a cartridge
        without subregion annotations is unsatisfiable — no V site
        would ever see a subregion factor, and the user is
        almost certainly building against the wrong cartridge or
        forgot to enable the annotation surface.

        Default rates (omitted kwarg, empty dict, or explicit
        all-ones expansion) skip the check — those are no-ops and
        compose cleanly with any cartridge.
        """
        if raw_rates is None or tuple_rates == _DEFAULT_V_SUBREGION_RATES:
            return
        annotated = 0
        v_total = self._refdata.v_pool_size()
        for v_id in range(v_total):
            if self._refdata.v_allele(v_id).subregions:
                annotated += 1
                # Early exit — we only need at least one annotated.
                return
        # Zero annotated V alleles: the user's rates can never bite.
        raise ValueError(
            "mutate(): v_subregion_rates was supplied but the bound "
            "cartridge carries no V-subregion annotations on any V "
            "allele (annotated_v_count=0 / "
            f"total_v_count={v_total}). The rate vector would be a "
            "deterministic no-op — almost certainly a builder bug. "
            "Either drop v_subregion_rates or use a cartridge with "
            "IMGT-gapped V sequences (the bundled human IGH / IGK / "
            "IGL OGRDB cartridges derive subregions automatically; "
            "see docs/v_region_substructure_audit.md)."
        )

    def _has_clonal_fork(self) -> bool:
        """Whether any clonal fork has already been appended.

        Used by the DSL ordering guards on :meth:`invert_d`,
        :meth:`receptor_revision`, and the clonal fork methods.
        Each fork has a pre-fork parent/founder phase; recombination-time
        mechanisms must be inherited by every descendant, emitted copy,
        or lineage node. Misordered calls used to lower into the wrong
        half and produce records with empty / default fields. The guards
        reject those configurations at the DSL boundary.
        """
        return any(
            isinstance(s, (_ClonalForkStep, _RepertoireForkStep, _LineageForkStep))
            for s in self._steps
        )

    def _is_tcr_refdata(self) -> bool:
        """Detect whether the bound refdata is a TCR locus.

        TCR allele names are prefixed with ``TR`` (TRA, TRB, TRG,
        TRD); BCR with ``IG``. Inspecting the first V allele is
        enough — locus is uniform across the pool. Returns ``False``
        when the V pool is empty (defensive — let downstream errors
        surface that condition).
        """
        if self._refdata.v_pool_size() == 0:
            return False
        first_v_name = self._refdata.v_allele(0).name
        return first_v_name.upper().startswith("TR")
