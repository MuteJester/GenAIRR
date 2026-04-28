"""
Compilation — translates Step descriptors into a C-backed simulator.

The compile path: resolve config → serialize to GDC → create CSimulator
→ apply feature flags → return CompiledSimulator handle.

This module is internal. Users interact through :class:`Experiment`.
"""

from __future__ import annotations

from typing import Any, Optional, Union

from .steps import Step, Rearrange

import warnings


# ─── Config name resolution ────────────────────────────────────
# Maps short string aliases (e.g., "human_igh", "mouse_tcrb")
# to builtin DataConfig objects.

_CONFIG_ALIASES = {
    # ── Short aliases (preferred defaults) ──
    # Human (OGRDB preferred for IG, IMGT for TCR)
    "human_igh": "HUMAN_IGH_OGRDB",
    "human_igk": "HUMAN_IGK_OGRDB",
    "human_igl": "HUMAN_IGL_OGRDB",
    "human_tcra": "HUMAN_TCRA_IMGT",
    "human_tcrb": "HUMAN_TCRB_IMGT",
    "human_tcrd": "HUMAN_TCRD_IMGT",
    "human_tcrg": "HUMAN_TCRG_IMGT",
    # Mouse
    "mouse_igh": "MOUSE_IGH_IMGT",
    "mouse_igk": "MOUSE_IGK_IMGT",
    "mouse_igl": "MOUSE_IGL_IMGT",
    "mouse_tcra": "MOUSE_TCRA_IMGT",
    "mouse_tcrb": "MOUSE_TCRB_IMGT",
    "mouse_tcrd": "MOUSE_TCRD_IMGT",
    "mouse_tcrg": "MOUSE_TCRG_IMGT",
    # Rat
    "rat_igh": "RAT_IGH_IMGT",
    "rat_igk": "RAT_IGK_IMGT",
    "rat_igl": "RAT_IGL_IMGT",
    # Rabbit
    "rabbit_igh": "RABBIT_IGH_IMGT",
    "rabbit_igk": "RABBIT_IGK_IMGT",
    "rabbit_igl": "RABBIT_IGL_IMGT",
    "rabbit_tcrb": "RABBIT_TCRB_IMGT",
    # Rhesus macaque
    "rhesus_igh": "RHESUS_IGH_IMGT",
    "rhesus_igk": "RHESUS_IGK_IMGT",
    "rhesus_igl": "RHESUS_IGL_IMGT",
    "rhesus_tcrb": "RHESUS_TCRB_IMGT",
    # Cow
    "cow_igh": "COW_IGH_IMGT",
    "cow_igk": "COW_IGK_IMGT",
    "cow_igl": "COW_IGL_IMGT",
    "cow_tcrb": "COW_TCRB_IMGT",
    # Dog
    "dog_igh": "DOG_IGH_IMGT",
    "dog_igk": "DOG_IGK_IMGT",
    "dog_igl": "DOG_IGL_IMGT",
    "dog_tcrb": "DOG_TCRB_IMGT",
    # Cat
    "cat_igk": "CAT_IGK_IMGT",
    "cat_igl": "CAT_IGL_IMGT",
    "cat_tcrb": "CAT_TCRB_IMGT",
    # Pig
    "pig_igh": "PIG_IGH_IMGT",
    "pig_igk": "PIG_IGK_IMGT",
    "pig_igl": "PIG_IGL_IMGT",
    "pig_tcrb": "PIG_TCRB_IMGT",
}


def _build_full_aliases():
    """Auto-register all builtin config names as aliases."""
    from .data import _CONFIG_NAMES
    for name in _CONFIG_NAMES:
        lower = name.lower()
        if lower not in _CONFIG_ALIASES:
            _CONFIG_ALIASES[lower] = name
        for suffix in ("_imgt", "_ogrdb"):
            if lower.endswith(suffix):
                short = lower[:-len(suffix)]
                if short not in _CONFIG_ALIASES:
                    _CONFIG_ALIASES[short] = name

_build_full_aliases()


def _resolve_config(config: Union[str, Any]) -> Any:
    """Resolve a config name string to a DataConfig object."""
    if isinstance(config, str):
        key = config.lower().replace("-", "_")
        const_name = _CONFIG_ALIASES.get(key)
        if const_name is None:
            from .data import list_configs
            available = list_configs()
            raise ValueError(
                f"Unknown config name {config!r}. "
                f"Available configs ({len(available)} total): "
                f"{', '.join(a.lower() for a in available[:10])}... "
                f"Use GenAIRR.data.list_configs() for full list, "
                f"or pass a DataConfig object directly."
            )
        from . import data
        return getattr(data, const_name)
    return config


# ─── Clause → Step lowering ───────────────────────────────────

def _get_chain_type(config):
    """Return the ChainType of a resolved DataConfig, or None."""
    if hasattr(config, 'metadata') and config.metadata:
        return config.metadata.chain_type
    return None


def _chain_has_d_segment(chain_type) -> bool:
    """Return True if the chain type has a D gene segment."""
    if chain_type is None:
        return True  # unknown → don't warn
    from .dataconfig.enums import ChainType
    return chain_type in (
        ChainType.BCR_HEAVY,
        ChainType.TCR_BETA,
        ChainType.TCR_DELTA,
    )


def _validate_using_alleles(merged_locks: dict, resolved_config) -> None:
    """T2-3: validate allele names in ``using()`` against the config's
    allele pools and raise a single ``ValueError`` with did-you-mean
    suggestions for every bad name.

    Surfaces typos at compile time with a clear message instead of the
    opaque C-level "allele not found in pool" stack trace from
    :meth:`CSimulator.lock_allele`.

    No-op when ``resolved_config`` is ``None`` (e.g. tests that lower
    clauses without a config).
    """
    if resolved_config is None or not merged_locks:
        return

    import difflib

    segment_attr = {"v": "v_alleles", "d": "d_alleles",
                    "j": "j_alleles", "c": "c_alleles"}
    errors: list[str] = []

    for segment, names in merged_locks.items():
        attr = segment_attr.get(segment)
        if attr is None:
            continue  # unknown segment — defer to C-side error
        pool = getattr(resolved_config, attr, None)
        if not pool:
            # Pool is None or empty; the C-side will report a clearer
            # error (e.g. "no D alleles for light chain"). Skip here.
            continue
        valid = {a.name for alleles in pool.values() for a in alleles}
        for name in names:
            if name in valid:
                continue
            suggestions = difflib.get_close_matches(
                name, valid, n=3, cutoff=0.6)
            seg_upper = segment.upper()
            if suggestions:
                errors.append(
                    f"Unknown {seg_upper} allele {name!r}. "
                    f"Did you mean: {', '.join(suggestions)}?")
            else:
                errors.append(
                    f"Unknown {seg_upper} allele {name!r}. "
                    f"No close matches in this config "
                    f"({len(valid)} {seg_upper} alleles available; "
                    f"use GenAIRR.data.<config>.{attr} to inspect).")

    if errors:
        raise ValueError(
            "using() clause references unknown allele(s):\n  - "
            + "\n  - ".join(errors))


def _set_or_warn(current, new_clause, op_name: str):
    """Assign a clause to a slot, warning on overwrite (T2-2).

    The DSL's last-write-wins is preserved for back-compat, but a
    duplicate within the same phase is almost always a typo (e.g.
    ``.mutate(rate(0.01,0.05), rate(0.10,0.20))``) — surface it via a
    ``UserWarning`` so users notice the dropped clause. To silence,
    remove the duplicate or filter ``UserWarning``.
    """
    if current is not None:
        warnings.warn(
            f"{op_name}() specified more than once in the same phase; "
            f"keeping last (last-write-wins). Remove the duplicate clause "
            f"to silence this warning.",
            UserWarning,
            stacklevel=5,
        )
    return new_clause


def _clauses_to_steps(
    recombine_clauses: list,
    mutate_clauses: list,
    prepare_clauses: list,
    sequence_clauses: list,
    observe_clauses: list,
    *,
    resolved_config=None,
) -> list:
    """
    Lower clause objects into a flat list of Step descriptors.

    Applies last-write-wins for unique clause types, segment-level
    merge for ``using()``, and emits compile-time warnings.

    Args:
        recombine_clauses: List of RecombineClause objects.
        mutate_clauses: List of MutateClause objects.
        prepare_clauses: List of PrepareClause objects.
        sequence_clauses: List of SequenceClause objects.
        observe_clauses: List of ObserveClause objects.
        resolved_config: Resolved DataConfig (for chain-type warnings).

    Returns:
        Ordered list of Step descriptors.
    """
    from .ops import (
        UsingClause, DInversionClause, ReceptorRevisionClause,
        ModelClause, RateClause, IsotypeRatesClause, AntigenSelectionClause,
        UMIClause, PCRClause, PrimerMaskClause,
        PairedEndClause, LongReadClause, FivePrimeLossClause,
        ThreePrimeLossClause, QualityProfileClause, ReverseComplementClause,
        ContaminantsClause, IndelsClause, NsClause,
    )
    from .steps import (
        Rearrange, LockAlleles, SimulateDGeneInversion,
        SimulateReceptorRevision, Mutate, SimulateCSR, SelectionPressure,
        PrimerMask, SimulateUMI, PCRAmplification,
        Corrupt5Prime, Corrupt3Prime, CorruptQuality, SimulatePairedEnd,
        SkewBaseComposition, ReverseComplement,
        SpikeContaminants, InsertIndels, InsertNs,
    )

    steps = [Rearrange()]

    chain_type = _get_chain_type(resolved_config) if resolved_config else None

    # ── Recombine ──────────────────────────────────────────────

    # Merge using() clauses at the segment level (last-write-wins per segment)
    merged_locks = {}
    has_d_inversion = None
    has_receptor_revision = None

    for c in recombine_clauses:
        if isinstance(c, UsingClause):
            if c.v is not None:
                merged_locks["v"] = list(c.v)
            if c.d is not None:
                merged_locks["d"] = list(c.d)
            if c.j is not None:
                merged_locks["j"] = list(c.j)
            if c.c is not None:
                merged_locks["c"] = list(c.c)
        elif isinstance(c, DInversionClause):
            has_d_inversion = c
        elif isinstance(c, ReceptorRevisionClause):
            has_receptor_revision = c

    if merged_locks:
        # T2-3: validate allele names before constructing the step so
        # typos surface at compile time with did-you-mean suggestions
        # instead of as an opaque C-side "not found in pool" error.
        _validate_using_alleles(merged_locks, resolved_config)
        steps.append(LockAlleles(locks=merged_locks))

    if has_d_inversion is not None:
        if not _chain_has_d_segment(chain_type):
            warnings.warn(
                "with_d_inversion() has no effect on this chain type "
                "(no D segment). Consider removing it.",
                UserWarning,
                stacklevel=4,
            )
        steps.append(SimulateDGeneInversion(probability=has_d_inversion.prob))

    if has_receptor_revision is not None:
        steps.append(SimulateReceptorRevision(
            probability=has_receptor_revision.prob,
            footprint_min=has_receptor_revision.footprint_min,
            footprint_max=has_receptor_revision.footprint_max,
        ))

    # ── Mutate ─────────────────────────────────────────────────

    if mutate_clauses:
        # Last-write-wins per clause type, with a UserWarning on duplicates
        # (T2-2) so accidental over-specification doesn't drop silently.
        model_clause = None
        rate_clause = None
        isotype_clause = None
        selection_clause = None

        for c in mutate_clauses:
            if isinstance(c, ModelClause):
                model_clause = _set_or_warn(model_clause, c, "model")
            elif isinstance(c, RateClause):
                rate_clause = _set_or_warn(rate_clause, c, "rate")
            elif isinstance(c, IsotypeRatesClause):
                isotype_clause = _set_or_warn(
                    isotype_clause, c, "with_isotype_rates")
            elif isinstance(c, AntigenSelectionClause):
                selection_clause = _set_or_warn(
                    selection_clause, c, "with_antigen_selection")

        model_name = model_clause.name if model_clause is not None else "s5f"
        if rate_clause is not None:
            min_rate, max_rate = rate_clause.min_rate, rate_clause.max_rate
        else:
            min_rate, max_rate = 0.01, 0.05
        has_csr = isotype_clause is not None

        # T2-4: ``model_name`` is now passed through to the C kernel
        # (s5f or uniform) — no silent S5F fallback for "uniform".
        if has_csr:
            steps.append(SimulateCSR(
                min_rate=min_rate, max_rate=max_rate, model=model_name))
        else:
            steps.append(Mutate(
                min_rate=min_rate, max_rate=max_rate, model=model_name))

        if selection_clause is not None:
            steps.append(SelectionPressure(
                strength=selection_clause.strength,
                cdr_r_acceptance=selection_clause.cdr_r_acceptance,
                fwr_r_acceptance=selection_clause.fwr_r_acceptance,
                anchor_r_acceptance=selection_clause.anchor_r_acceptance,
            ))

    # ── Prepare ────────────────────────────────────────────────

    # Last-write-wins per clause type with duplicate-warning (T2-2).
    primer_mask_clause = None
    umi_clause = None
    pcr_clause = None

    for c in prepare_clauses:
        if isinstance(c, PrimerMaskClause):
            primer_mask_clause = _set_or_warn(
                primer_mask_clause, c, "with_primer_mask")
        elif isinstance(c, UMIClause):
            umi_clause = _set_or_warn(umi_clause, c, "with_umi")
        elif isinstance(c, PCRClause):
            pcr_clause = _set_or_warn(pcr_clause, c, "with_pcr")

    if primer_mask_clause is not None:
        steps.append(PrimerMask(mask_length=primer_mask_clause.length))
    if umi_clause is not None:
        steps.append(SimulateUMI(umi_length=umi_clause.length))
    if pcr_clause is not None:
        steps.append(PCRAmplification(
            error_rate=pcr_clause.error_rate,
            n_cycles=pcr_clause.cycles,
        ))

    # ── Sequence ───────────────────────────────────────────────

    # Last-write-wins per clause type with duplicate-warning (T2-2).
    paired_end_clause = None
    long_read_clause = None
    five_prime_clause = None
    three_prime_clause = None
    quality_clause = None
    rc_clause = None

    for c in sequence_clauses:
        if isinstance(c, PairedEndClause):
            paired_end_clause = _set_or_warn(paired_end_clause, c, "paired_end")
        elif isinstance(c, LongReadClause):
            long_read_clause = _set_or_warn(long_read_clause, c, "long_read")
        elif isinstance(c, FivePrimeLossClause):
            five_prime_clause = _set_or_warn(
                five_prime_clause, c, "with_5prime_loss")
        elif isinstance(c, ThreePrimeLossClause):
            three_prime_clause = _set_or_warn(
                three_prime_clause, c, "with_3prime_loss")
        elif isinstance(c, QualityProfileClause):
            quality_clause = _set_or_warn(
                quality_clause, c, "with_quality_profile")
        elif isinstance(c, ReverseComplementClause):
            rc_clause = _set_or_warn(rc_clause, c, "with_reverse_complement")

    # Warn about technology conflict
    if paired_end_clause is not None and long_read_clause is not None:
        warnings.warn(
            "Both paired_end() and long_read() specified — sequencing "
            "technology conflict. Both will be applied.",
            UserWarning,
            stacklevel=4,
        )

    if five_prime_clause is not None:
        steps.append(Corrupt5Prime(
            min_remove=five_prime_clause.min_remove,
            max_remove=five_prime_clause.max_remove,
            min_add=five_prime_clause.min_add,
            max_add=five_prime_clause.max_add,
        ))
    if three_prime_clause is not None:
        steps.append(Corrupt3Prime(
            min_remove=three_prime_clause.min_remove,
            max_remove=three_prime_clause.max_remove,
            min_add=three_prime_clause.min_add,
            max_add=three_prime_clause.max_add,
        ))
    if quality_clause is not None:
        steps.append(CorruptQuality(
            base_error_rate=quality_clause.base,
            peak_error_rate=quality_clause.peak,
        ))
    if paired_end_clause is not None:
        steps.append(SimulatePairedEnd(read_length=paired_end_clause.read_length))
    if long_read_clause is not None:
        steps.append(SkewBaseComposition(
            error_rate=long_read_clause.error_rate,
            min_run_length=long_read_clause.min_run_length,
            insertion_bias=long_read_clause.insertion_bias,
        ))
    if rc_clause is not None:
        steps.append(ReverseComplement(probability=rc_clause.prob))

    # ── Observe ────────────────────────────────────────────────

    # Last-write-wins per clause type with duplicate-warning (T2-2).
    contaminants_clause = None
    indels_clause = None
    ns_clause = None

    for c in observe_clauses:
        if isinstance(c, ContaminantsClause):
            contaminants_clause = _set_or_warn(
                contaminants_clause, c, "with_contaminants")
        elif isinstance(c, IndelsClause):
            indels_clause = _set_or_warn(indels_clause, c, "with_indels")
        elif isinstance(c, NsClause):
            ns_clause = _set_or_warn(ns_clause, c, "with_ns")

    if contaminants_clause is not None:
        steps.append(SpikeContaminants(
            probability=contaminants_clause.rate,
            contaminant_type=contaminants_clause.source,
        ))
    if indels_clause is not None:
        steps.append(InsertIndels(probability=indels_clause.prob))
    if ns_clause is not None:
        steps.append(InsertNs(probability=ns_clause.prob))

    return steps


# ─── CompiledSimulator ────────────────────────────────────────

class CompiledSimulator:
    """
    A compiled simulation backed by the C engine.

    Created by :meth:`Experiment.compile`. Holds a CSimulator handle
    and provides ``simulate()``, ``simulate_to_file()``, ``stream()``,
    and ``set_seed()`` methods.

    **RNG semantics (streaming).** The simulator owns a single PCG32
    stream that persists across ``simulate()`` calls. The seed is
    applied once — on the first simulation call — and successive calls
    advance the same stream:

    .. code-block:: python

        sim = exp.compile(seed=42)
        a = sim.simulate(10)        # records 0..9   from seed 42
        b = sim.simulate(10)        # records 10..19 from seed 42 (NOT 0..9 again)

        # Equivalent two-batch concatenation invariant:
        sim2 = exp.compile(seed=42)
        c = sim2.simulate(20)       # records 0..19 from seed 42
        # → a + b and c contain the same records in the same order.

    To reset the stream mid-life, call :meth:`set_seed`:

    .. code-block:: python

        sim.set_seed(42); a = sim.simulate(10)
        sim.set_seed(42); b = sim.simulate(10)
        # → a and b are byte-identical.

    Auto-seed: ``compile(seed=0)`` (or no ``seed=``) lets the C engine
    pick a seed from wall time XOR the simulator address, so two
    concurrent simulators in the same process get distinct streams.
    """

    def __init__(self, c_sim):
        self._sim = c_sim

    def set_seed(self, seed: int) -> None:
        """
        Reset the RNG stream to a new seed.

        The next call to :meth:`simulate` (or :meth:`simulate_to_file`,
        :meth:`stream`) will re-seed the stream and start fresh. Use
        this to reproduce a deterministic batch after the simulator has
        already been used:

        .. code-block:: python

            sim.set_seed(42); a = sim.simulate(10)
            sim.set_seed(42); b = sim.simulate(10)
            assert a == b

        Args:
            seed: Non-negative 64-bit integer. Pass 0 for the engine's
                "auto-seed from time + simulator address" mode.
        """
        if self._sim is None:
            raise RuntimeError("CompiledSimulator is closed")
        self._sim.set_seed(seed)

    def simulate(self, n: int = 1, *, progress: bool = False) -> Any:
        """
        Run N simulations and return a SimulationResult.

        The RNG stream advances from where the previous ``simulate()``
        call left off. To reset, call :meth:`set_seed` between batches.
        See the class docstring for the full streaming semantics.

        Args:
            n: Number of sequences to generate.
            progress: If True, show a progress bar (requires tqdm).

        Returns:
            SimulationResult with AIRR-format dicts.
        """
        from .result import SimulationResult

        if self._sim is None:
            raise RuntimeError("CompiledSimulator is closed")

        if progress:
            try:
                from tqdm import tqdm
                results = []
                chunk = max(1, n // 100)
                with tqdm(total=n, desc="Simulating") as pbar:
                    remaining = n
                    while remaining > 0:
                        batch = min(chunk, remaining)
                        results.extend(self._sim.simulate(batch))
                        remaining -= batch
                        pbar.update(batch)
            except ImportError:
                results = self._sim.simulate(n)
        else:
            results = self._sim.simulate(n)

        return SimulationResult(results)

    def simulate_to_file(self, n: int, output_path: str) -> int:
        """Write AIRR TSV directly to a file (no Python parsing)."""
        if self._sim is None:
            raise RuntimeError("CompiledSimulator is closed")
        return self._sim.simulate_to_file(n, output_path)

    def stream(self) -> SimulationStream:
        """
        Return a streaming interface for on-demand sequence generation.

        Each call to :meth:`SimulationStream.get` simulates one or more
        sequences directly in memory — no temp files or disk I/O.

        Can be used as a context manager or as an iterator::

            sim = Experiment.on("human_igh").compile(seed=42)

            # Context manager with get()
            with sim.stream() as s:
                seq = s.get()        # single dict
                batch = s.get(10)    # list of 10 dicts

            # Iterator (infinite — break when done)
            for seq in sim.stream():
                process(seq)
                if done:
                    break

        Returns:
            SimulationStream wrapping the C simulator.
        """
        if self._sim is None:
            raise RuntimeError("CompiledSimulator is closed")
        return SimulationStream(self._sim)

    def __iter__(self):
        """Iterate over simulated sequences (infinite stream)."""
        if self._sim is None:
            raise RuntimeError("CompiledSimulator is closed")
        return iter(SimulationStream(self._sim))

    # ── Lifecycle (T2-5) ──────────────────────────────────────────

    def close(self) -> None:
        """Release the underlying C simulator handle.

        Idempotent — safe to call multiple times. After ``close()``,
        :meth:`simulate`, :meth:`stream`, and iteration raise
        ``RuntimeError``. Use this in long-running scripts that build
        many compiled simulators (parameter sweeps, batch jobs) so
        C-side memory is freed deterministically rather than waiting
        on the Python GC.

        Most users should prefer the context-manager form::

            with Experiment.on("human_igh").compile(seed=42) as sim:
                results = sim.simulate(1000)
            # sim is closed here, even if simulate() raised.
        """
        if self._sim is not None:
            self._sim.close()
            self._sim = None

    @property
    def closed(self) -> bool:
        return self._sim is None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False  # don't suppress exceptions

    def __repr__(self):
        return "CompiledSimulator(closed)" if self.closed else "CompiledSimulator()"


class SimulationStream:
    """
    Streaming interface for on-demand sequence generation.

    Generates sequences one at a time directly in memory via the C engine,
    bypassing temp-file I/O. Supports three usage patterns::

        # 1. Context manager with get()
        with sim.stream() as s:
            seq = s.get()        # single AIRR dict
            batch = s.get(10)    # list of 10 AIRR dicts

        # 2. Iterator (infinite)
        for seq in sim.stream():
            ...

        # 3. Direct calls
        s = sim.stream()
        seq = s.get()
    """

    def __init__(self, c_sim):
        self._sim = c_sim

    def get(self, n: int = 1):
        """
        Generate one or more sequences.

        Args:
            n: Number of sequences to generate. Defaults to 1.

        Returns:
            A single AIRR dict if n=1, or a list of AIRR dicts if n>1.
        """
        if n == 1:
            return self._sim.simulate_one()
        return [self._sim.simulate_one() for _ in range(n)]

    def __iter__(self):
        return self

    def __next__(self):
        return self._sim.simulate_one()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def __repr__(self):
        return "SimulationStream()"


# ─── Compile function ─────────────────────────────────────────

def _compile_to_c(
    config: Union[str, Any],
    steps: list[Step],
    *,
    productivity: Optional[Any] = None,
    seed: Optional[int] = None,
) -> CompiledSimulator:
    """
    Compile a list of Steps into a C-backed simulator.

    This is the single compilation path used by Experiment.

    Args:
        config: DataConfig object or string name (e.g. "human_igh").
        steps: Ordered list of Step descriptors.
        productivity: Productivity enum filter mode. None defaults to
            Productivity.PRODUCTIVE_MIXED (no filtering).
        seed: Random seed for reproducibility. None falls back to
            ``GenAIRR.set_seed()``'s global value if one was set;
            otherwise the C engine auto-seeds from time + simulator
            address. An explicit ``seed=0`` is passed through and
            interpreted by the C engine as "auto-seed me".

    Returns:
        CompiledSimulator ready to call simulate().
    """
    from .dataconfig.gdc_io import to_gdc_bytes
    from .dataconfig.enums import ChainType, Productivity
    from ._native import CSimulator
    from .seed import get_seed as _get_global_seed

    if productivity is None:
        productivity = Productivity.PRODUCTIVE_MIXED

    resolved = _resolve_config(config)

    # Determine S5F chain category
    has_mutate = any(
        hasattr(s, '_feature') and s._feature in ('mutate', 'csr')
        for s in steps
    )
    s5f_cat = None
    if has_mutate and hasattr(resolved, 'metadata') and resolved.metadata:
        ct = resolved.metadata.chain_type
        if ct == ChainType.BCR_HEAVY:
            s5f_cat = "heavy"
        elif ct in (ChainType.BCR_LIGHT_KAPPA, ChainType.BCR_LIGHT_LAMBDA):
            s5f_cat = "light"
        else:
            s5f_cat = "heavy"  # TCR defaults to heavy S5F

    # Serialize DataConfig to GDC bytes
    gdc_bytes = to_gdc_bytes(resolved, s5f_chain_category=s5f_cat)

    # Create C simulator from in-memory GDC
    sim = CSimulator(gdc_bytes=gdc_bytes)

    # Apply productivity mode (PRODUCTIVE_MIXED is the implicit default
    # — leaves cfg.features.productivity at PRODUCTIVITY_MIXED == 0).
    if productivity == Productivity.PRODUCTIVE_ONLY:
        sim.set_param("productivity_mode", 1)  # PRODUCTIVITY_PRODUCTIVE_ONLY
    elif productivity == Productivity.NON_PRODUCTIVE_ONLY:
        sim.set_param("productivity_mode", 2)  # PRODUCTIVITY_NON_PRODUCTIVE_ONLY

    # Apply all steps
    for step in steps:
        step.configure(sim)

    # Set seed. When the caller did not pass one, fall back to the
    # global seed set via GenAIRR.set_seed(), if any. If both are
    # None, leave it for the C engine to auto-seed from time +
    # simulator address (see ensure_seeded in genairr_api.c).
    if seed is None:
        seed = _get_global_seed()
    if seed is not None:
        sim.set_seed(seed)

    return CompiledSimulator(sim)
