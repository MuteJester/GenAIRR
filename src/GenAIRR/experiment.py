"""GenAIRR — fluent DSL for receptor-sequence simulation.

The :class:`Experiment` builder lowers a chain of fluent steps into a
Rust ``PassPlan`` (via :mod:`GenAIRR._engine`), compiles that IR into an
owning ``CompiledSimulator``, and runs it to produce a list of
:class:`GenAIRR._engine.Outcome` objects.

Typical usage::

    import GenAIRR as ga
    outcomes = ga.Experiment.on("human_igh").recombine().run(n=100, seed=42)

``Experiment.on`` accepts:

- a config-name string (``"human_igh"``, ``"mouse_tcrb"``, …) — resolved
  through the builtin :mod:`GenAIRR.data` registry,
- a :class:`GenAIRR.DataConfig` object (already-loaded reference data),
- a :class:`GenAIRR._engine.RefDataConfig` (the engine-native form, used
  primarily by tests and advanced callers).

In all three cases the input is normalised to a
:class:`GenAIRR._engine.RefDataConfig` before any pass is appended.
"""
from __future__ import annotations

import warnings
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple, Union

from GenAIRR import _engine  # private Rust extension submodule

from .dataconfig import DataConfig
from ._describe import (
    _describe_experiment_header,
    _describe_step_sequence,
    _format_declared_contracts,
)
from ._normalize import (
    _normalize_count,
    _normalize_lengths,
    _to_immutable_byte_pair_matrix,
    _to_immutable_byte_pairs,
    _to_immutable_pairs,
)
from ._lowering import (
    _extract_invert_d_prob,
    _extract_paired_end_step,
    _extract_receptor_revision_prob,
    _lower_paired_end,
    _lower_recombine,
    lower_step,
)
from ._compiled import (
    CompiledClonalExperiment,
    CompiledExperiment,
    CompiledLineageExperiment,
    CompiledRepertoireExperiment,
)
from ._refdata_resolver import (
    _CONFIG_ALIASES,
    ExperimentInput,
    _coerce_to_refdata_and_dataconfig,
    dataconfig_to_refdata,
)
from ._pipeline_ir import (
    _CORRUPT_KIND_3PRIME_LOSS,
    _CORRUPT_KIND_5PRIME_LOSS,
    _CORRUPT_KIND_CONTAMINANT,
    _CORRUPT_KIND_INDEL,
    _CORRUPT_KIND_NS,
    _CORRUPT_KIND_PCR,
    _CORRUPT_KIND_QUALITY,
    _CORRUPT_KIND_REV_COMP,
    _DEFAULT_NP_LENGTHS,
    _ClonalForkStep,
    _CorruptStep,
    _InvertDStep,
    _LineageForkStep,
    _MutateStep,
    _PairedEndStep,
    _RecombineStep,
    _ReceptorRevisionStep,
    _RepertoireForkStep,
)
from ._step_validation import (
    _DEFAULT_V_SUBREGION_RATES,
    _V_SUBREGION_RATE_ALIASES,
    _V_SUBREGION_RATE_LABELS,
    _descendant_phase_step_classifier,
    _validate_segment_rates,
    _validate_v_subregion_rates,
)

from ._experiment import (
    _ClonalMixin,
    _CompileMixin,
    _ConstraintsMixin,
    _CorruptionMixin,
    _GenotypeAllelesMixin,
    _IntrospectionMixin,
    _MutationMixin,
    _RecombinationMixin,
    _RefdataControlsMixin,
    _RunMixin,
)


# Refdata-resolver helpers (``_CONFIG_ALIASES``, ``_resolve_config_name``,
# ``dataconfig_to_refdata``, ``_coerce_to_refdata_and_dataconfig``,
# ``ExperimentInput``) live in :mod:`._refdata_resolver`. Re-exported
# below so the public surface (``GenAIRR.dataconfig_to_refdata``) and
# legacy internal callers (e.g. ``mcp_server`` reaching for
# ``_CONFIG_ALIASES``) keep working.


class Experiment(
    _CorruptionMixin,
    _ConstraintsMixin,
    _ClonalMixin,
    _MutationMixin,
    _GenotypeAllelesMixin,
    _RecombinationMixin,
    _RefdataControlsMixin,
    _CompileMixin,
    _RunMixin,
    _IntrospectionMixin,
):
    """Fluent builder for a simulation pipeline.

    Build an ``Experiment`` from a config name, a :class:`DataConfig`,
    or a :class:`GenAIRR._engine.RefDataConfig` via :meth:`Experiment.on`,
    chain configuration steps (currently just :meth:`recombine`), then
    call :meth:`run` (or :meth:`compile` for explicit two-stage flow).

    The builder is **stateful but not destructive**: each fluent call
    returns ``self`` after appending a step. The same ``Experiment``
    can be ``compile()``-d and ``run()``-d multiple times with
    different seeds.
    """

    __slots__ = (
        "_refdata",
        "_steps",
        "_dataconfig",
        "_locks",
        "_metadata",
        "_contracts",
        "_allow_curatable_refdata",
        "_genotype",
        "_user_allele_weights_set",
    )

    def __init__(
        self,
        refdata: "_engine.RefDataConfig",
        dataconfig: Optional[DataConfig] = None,
    ) -> None:
        self._refdata = refdata
        self._dataconfig = dataconfig
        self._steps: List[
            Union[
                _RecombineStep,
                _MutateStep,
                _CorruptStep,
                _InvertDStep,
                _ReceptorRevisionStep,
                _PairedEndStep,
            ]
        ] = []
        # Constraint bundles declared via `.productive_only()` etc.
        # Stored as a list so future bundles can compose. Today only
        # the productive bundle is recognized.
        self._contracts: List[str] = []
        # Per-segment allele-lock subsets set by ``.restrict_alleles(...)``. Each
        # entry is ``None`` (no lock — sample uniformly across the pool)
        # or a tuple of allele IDs to sample uniformly across.
        self._locks: Dict[str, Optional[Tuple[int, ...]]] = {
            "V": None,
            "D": None,
            "J": None,
        }
        # sample-level metadata to inject as columns on every
        # AIRR record (e.g. ``sample_id``, ``donor``,
        # ``repertoire_id``, ``cell_id``). Empty by default.
        self._metadata: Dict[str, Any] = {}
        # When True, ``compile()`` / ``run()`` / ``run_records()``
        # runs the refdata gate under the lenient `AllowCuratable`
        # mode. Fatal issues (empty pool, duplicates, invalid byte,
        # anchor out of bounds) still reject; Curatable issues (V
        # anchor not Cys, J anchor unexpected AA, missing anchor)
        # pass. Production users opt in via
        # ``.allow_curatable_refdata()`` when sampling from a real
        # catalogue (bundled mouse_igh / human_tcrb) that includes
        # pseudogene/ORF alleles.
        self._allow_curatable_refdata: bool = False
        # Single-subject diploid genotype attached via ``with_genotype``.
        # ``None`` => the flat (uniform/usage-weighted) allele path runs
        # unchanged. When set, recombination lowers to the phased
        # genotype path (one ``SampleGenotypePass``).
        self._genotype = None
        # True once the user passed an explicit ``*_allele_weights`` to
        # ``recombine`` — distinct from cartridge-usage defaults. Used to
        # enforce mutual exclusion with ``with_genotype``.
        self._user_allele_weights_set: bool = False

    @classmethod
    def on(cls, source: ExperimentInput) -> "Experiment":
        """Start an experiment against the given reference data.

        ``source`` is one of:
        - a config-name string (e.g. ``"human_igh"``),
        - a :class:`GenAIRR.DataConfig` instance,
        - a :class:`GenAIRR._engine.RefDataConfig`.

        When ``source`` is a config name or a ``DataConfig``, the
        underlying empirical distributions (NP lengths, per-gene
        trims) are kept on the experiment so :meth:`recombine` can
        use them as the default sampling distributions. A bare
        ``RefDataConfig`` has no such backing — :meth:`recombine`
        falls through to its uniform ``[0..6]`` placeholder unless
        the caller passes ``np1_lengths`` / ``np2_lengths``
        explicitly.
        """
        refdata, dataconfig = _coerce_to_refdata_and_dataconfig(source)
        return cls(refdata, dataconfig)

    @property
    def chain_type(self) -> str:
        """Chain type of the attached refdata (``"vj"`` or ``"vdj"``)."""
        return self._refdata.chain_type

    @property
    def step_count(self) -> int:
        """Number of fluent steps recorded on this experiment so far."""
        return len(self._steps)

    @property
    def refdata(self) -> "_engine.RefDataConfig":
        """The engine-native ``RefDataConfig`` this experiment is bound to."""
        return self._refdata
