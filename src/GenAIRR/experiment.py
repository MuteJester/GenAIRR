"""
Experiment — fluent DSL for building simulation protocols.

A GenAIRR Experiment reads like a biological protocol specification.
Each phase of the wet-lab process is an explicit section, and clause
objects describe exactly what happens in each phase.

Example::

    from GenAIRR import Experiment
    from GenAIRR.ops import (
        rate, model, with_antigen_selection,
        with_5prime_loss, with_3prime_loss, paired_end,
        with_contaminants, with_indels, with_ns,
    )

    result = (
        Experiment.on("human_igh")
        .mutate(rate(0.02, 0.08), model("s5f"), with_antigen_selection(0.7))
        .sequence(with_5prime_loss(), with_3prime_loss(), paired_end(300))
        .observe(with_contaminants(0.01), with_indels(), with_ns())
        .run(n=1000, seed=42)
    )
"""

from __future__ import annotations

from typing import Any, List, Optional, Union


# ---------------------------------------------------------------------------
# Experiment
# ---------------------------------------------------------------------------

class Experiment:
    """
    Fluent DSL for defining AIRR simulation protocols.

    Build an experiment by chaining phases and clause objects. The experiment
    compiles down to a C-backed simulator for native-speed execution.

    Phases (in biological order):

    - ``.recombine(...)`` — V(D)J recombination configuration
    - ``.mutate(...)`` — Somatic hypermutation
    - ``.prepare(...)`` — Library preparation (primers, UMI, PCR)
    - ``.sequence(...)`` — Sequencing artifacts
    - ``.observe(...)`` — Post-sequencing contamination and errors

    Terminal methods:

    - ``.run()`` — Compile and simulate in one call
    - ``.compile()`` — Get a reusable CompiledSimulator
    """

    def __init__(self, config: Union[str, Any]):
        self._config = config
        self._recombine_clauses: List[Any] = []
        self._mutate_clauses: List[Any] = []
        self._prepare_clauses: List[Any] = []
        self._sequence_clauses: List[Any] = []
        self._observe_clauses: List[Any] = []

    # ── Factory ───────────────────────────────────────────────────

    @classmethod
    def on(cls, config: Union[str, Any]) -> Experiment:
        """
        Start an experiment on a species/chain configuration.

        Args:
            config: DataConfig or string name (e.g., ``"human_igh"``,
                    ``"mouse_igk"``, ``"rabbit_tcrb"``).

        Example::

            exp = Experiment.on("human_igh")
        """
        return cls(config)

    # ── Phase validation ──────────────────────────────────────────

    @staticmethod
    def _validate_phase_clauses(phase_name: str, base_type: type,
                                clauses: tuple) -> None:
        """Check every clause is an instance of the expected base type."""
        for c in clauses:
            if not isinstance(c, base_type):
                raise TypeError(
                    f".{phase_name}() expects {base_type.__name__} objects, "
                    f"got {type(c).__name__}: {c!r}"
                )

    # ==================================================================
    # Phase methods
    # ==================================================================

    def recombine(self, *clauses) -> Experiment:
        """
        V(D)J recombination phase.

        Pass :mod:`GenAIRR.ops` clause objects to configure recombination.

        Example::

            from GenAIRR.ops import using, with_d_inversion

            Experiment.on("human_igh")
            .recombine(using(v="IGHV1-2*01"), with_d_inversion(0.15))
        """
        if clauses:
            from .ops import RecombineClause
            self._validate_phase_clauses("recombine", RecombineClause, clauses)
            self._recombine_clauses.extend(clauses)
        return self

    def mutate(self, *clauses) -> Experiment:
        """
        Somatic hypermutation phase.

        Pass :mod:`GenAIRR.ops` clause objects to configure mutation.

        Example::

            from GenAIRR.ops import rate, model, with_antigen_selection

            Experiment.on("human_igh")
            .mutate(rate(0.02, 0.08), model("s5f"), with_antigen_selection(0.5))
        """
        if clauses:
            from .ops import MutateClause
            self._validate_phase_clauses("mutate", MutateClause, clauses)
            self._mutate_clauses.extend(clauses)
        return self

    def prepare(self, *clauses) -> Experiment:
        """
        Library preparation phase.

        Example::

            from GenAIRR.ops import with_umi, with_pcr, with_primer_mask

            Experiment.on("human_igh")
            .prepare(with_primer_mask(), with_umi(12), with_pcr())
        """
        if clauses:
            from .ops import PrepareClause
            self._validate_phase_clauses("prepare", PrepareClause, clauses)
            self._prepare_clauses.extend(clauses)
        return self

    def sequence(self, *clauses) -> Experiment:
        """
        Sequencing phase.

        Example::

            from GenAIRR.ops import with_5prime_loss, paired_end

            Experiment.on("human_igh")
            .sequence(with_5prime_loss(), paired_end(300))
        """
        if clauses:
            from .ops import SequenceClause
            self._validate_phase_clauses("sequence", SequenceClause, clauses)
            self._sequence_clauses.extend(clauses)
        return self

    def observe(self, *clauses) -> Experiment:
        """
        Post-sequencing observation phase (contamination, indels, Ns).

        Example::

            from GenAIRR.ops import with_contaminants, with_indels, with_ns

            Experiment.on("human_igh")
            .observe(with_contaminants(0.01), with_indels(), with_ns())
        """
        if clauses:
            from .ops import ObserveClause
            self._validate_phase_clauses("observe", ObserveClause, clauses)
            self._observe_clauses.extend(clauses)
        return self

    # ==================================================================
    # Terminal methods
    # ==================================================================

    def _build_steps(self) -> list:
        """Lower clause lists into an ordered list of Step descriptors."""
        from .protocol import _clauses_to_steps, _resolve_config
        resolved = _resolve_config(self._config)
        return _clauses_to_steps(
            self._recombine_clauses,
            self._mutate_clauses,
            self._prepare_clauses,
            self._sequence_clauses,
            self._observe_clauses,
            resolved_config=resolved,
        )

    def compile(self, *, seed: Optional[int] = None,
                productive: bool = False) -> Any:
        """
        Compile this experiment into a C-backed simulator.

        Returns a :class:`CompiledSimulator` that can be called
        repeatedly with ``.simulate(n)``.

        Args:
            seed: Random seed for reproducibility.
            productive: If True, enforce productive rearrangements.

        Example::

            sim = exp.compile(seed=42, productive=True)
            result = sim.simulate(n=1000)
        """
        from .protocol import _compile_to_c
        steps = self._build_steps()
        return _compile_to_c(
            config=self._config,
            steps=steps,
            productive=productive,
            seed=seed,
        )

    def run(self, n: int = 1, *, seed: Optional[int] = None,
            productive: bool = False,
            progress: bool = False) -> Any:
        """
        Compile and run the experiment in one call.

        Args:
            n: Number of sequences to simulate.
            seed: Random seed for reproducibility.
            productive: If True, enforce productive rearrangements.
            progress: If True, show a progress bar (requires tqdm).

        Returns:
            :class:`~GenAIRR.result.SimulationResult`

        Example::

            result = exp.run(n=1000, seed=42, productive=True)
            df = result.to_dataframe()
        """
        sim = self.compile(seed=seed, productive=productive)
        return sim.simulate(n=n, progress=progress)

    # ==================================================================
    # Display
    # ==================================================================

    def _resolve_config_name(self) -> str:
        """Get a display name for the config."""
        if isinstance(self._config, str):
            return self._config
        if hasattr(self._config, 'metadata') and self._config.metadata:
            m = self._config.metadata
            parts = []
            if hasattr(m, 'species') and m.species:
                parts.append(str(m.species.name))
            if hasattr(m, 'chain_type') and m.chain_type:
                parts.append(str(m.chain_type.name))
            if parts:
                return "_".join(parts).lower()
        return repr(self._config)

    def __repr__(self) -> str:
        config_name = self._resolve_config_name()
        w = 52

        lines = []
        lines.append(f"{'':>2}{'':─>{w}}")
        lines.append(f"  Experiment on \"{config_name}\"")
        lines.append(f"{'':>2}{'':─>{w}}")

        phase_defs = [
            ("Recombine", self._recombine_clauses, "V(D)J rearrangement"),
            ("Mutate", self._mutate_clauses, None),
            ("Prepare", self._prepare_clauses, None),
            ("Sequence", self._sequence_clauses, None),
            ("Observe", self._observe_clauses, None),
        ]

        for phase_name, clause_list, default_msg in phase_defs:
            if not clause_list and default_msg is None:
                continue

            lines.append("")
            pad = w - len(phase_name) - 3
            lines.append(f"  ┌─ {phase_name} {'─' * max(0, pad)}┐")

            if clause_list:
                for i, clause in enumerate(clause_list):
                    is_last = (i == len(clause_list) - 1)
                    connector = "└" if is_last else "├"
                    label = clause.__label__()
                    inner_pad = w - len(label) - 5
                    lines.append(
                        f"  │  {connector}── {label}"
                        f"{' ' * max(0, inner_pad)}│"
                    )
            elif default_msg:
                inner_pad = w - len(default_msg) - 5
                lines.append(
                    f"  │  {default_msg}{' ' * max(0, inner_pad)}│"
                )

            lines.append(f"  └{'─' * (w)}┘")

        lines.append("")
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.__repr__()
