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
        # Last-write-wins for each clause type
        model_name = "s5f"
        min_rate, max_rate = 0.01, 0.05
        has_csr = False
        selection_clause = None

        for c in mutate_clauses:
            if isinstance(c, ModelClause):
                model_name = c.name
            elif isinstance(c, RateClause):
                min_rate, max_rate = c.min_rate, c.max_rate
            elif isinstance(c, IsotypeRatesClause):
                has_csr = True
            elif isinstance(c, AntigenSelectionClause):
                selection_clause = c

        # Warn about uniform model
        if model_name == "uniform":
            warnings.warn(
                'model("uniform") requested, but the C backend currently uses '
                "S5F motif data from the DataConfig for all mutation. True "
                "position-independent mutation will be implemented in a future "
                "release. Sequences will be generated using S5F mutation patterns.",
                RuntimeWarning,
                stacklevel=4,
            )

        if has_csr:
            steps.append(SimulateCSR(min_rate=min_rate, max_rate=max_rate))
        else:
            steps.append(Mutate(min_rate=min_rate, max_rate=max_rate))

        if selection_clause is not None:
            steps.append(SelectionPressure(
                strength=selection_clause.strength,
                cdr_r_acceptance=selection_clause.cdr_r_acceptance,
                fwr_r_acceptance=selection_clause.fwr_r_acceptance,
            ))

    # ── Prepare ────────────────────────────────────────────────

    # Last-write-wins per clause type
    primer_mask_clause = None
    umi_clause = None
    pcr_clause = None

    for c in prepare_clauses:
        if isinstance(c, PrimerMaskClause):
            primer_mask_clause = c
        elif isinstance(c, UMIClause):
            umi_clause = c
        elif isinstance(c, PCRClause):
            pcr_clause = c

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

    # Last-write-wins per clause type
    paired_end_clause = None
    long_read_clause = None
    five_prime_clause = None
    three_prime_clause = None
    quality_clause = None
    rc_clause = None

    for c in sequence_clauses:
        if isinstance(c, PairedEndClause):
            paired_end_clause = c
        elif isinstance(c, LongReadClause):
            long_read_clause = c
        elif isinstance(c, FivePrimeLossClause):
            five_prime_clause = c
        elif isinstance(c, ThreePrimeLossClause):
            three_prime_clause = c
        elif isinstance(c, QualityProfileClause):
            quality_clause = c
        elif isinstance(c, ReverseComplementClause):
            rc_clause = c

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

    contaminants_clause = None
    indels_clause = None
    ns_clause = None

    for c in observe_clauses:
        if isinstance(c, ContaminantsClause):
            contaminants_clause = c
        elif isinstance(c, IndelsClause):
            indels_clause = c
        elif isinstance(c, NsClause):
            ns_clause = c

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
    and provides ``simulate()`` and ``simulate_to_file()`` methods.
    """

    def __init__(self, c_sim):
        self._sim = c_sim

    def simulate(self, n: int = 1, *, progress: bool = False) -> Any:
        """
        Run N simulations and return a SimulationResult.

        Args:
            n: Number of sequences to generate.
            progress: If True, show a progress bar (requires tqdm).

        Returns:
            SimulationResult with AIRR-format dicts.
        """
        from .result import SimulationResult

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
        return SimulationStream(self._sim)

    def __iter__(self):
        """Iterate over simulated sequences (infinite stream)."""
        return iter(SimulationStream(self._sim))

    def __repr__(self):
        return "CompiledSimulator()"


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
    productive: bool = False,
    seed: Optional[int] = None,
) -> CompiledSimulator:
    """
    Compile a list of Steps into a C-backed simulator.

    This is the single compilation path used by Experiment.

    Args:
        config: DataConfig object or string name (e.g. "human_igh").
        steps: Ordered list of Step descriptors.
        productive: If True, enforce productive sequences (retry loop).
        seed: Random seed for reproducibility. None falls back to
            ``GenAIRR.set_seed()``'s global value if one was set;
            otherwise the C engine auto-seeds from time + simulator
            address. An explicit ``seed=0`` is passed through and
            interpreted by the C engine as "auto-seed me".

    Returns:
        CompiledSimulator ready to call simulate().
    """
    from .dataconfig.gdc_io import to_gdc_bytes
    from .dataconfig.enums import ChainType
    from ._native import CSimulator
    from .seed import get_seed as _get_global_seed

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

    # Apply productive mode
    if productive:
        sim.set_feature("productive", True)

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
