from __future__ import annotations

from typing import Any, List

from .._describe import (
    _describe_experiment_header,
    _describe_step_sequence,
    _format_declared_contracts,
)
from .._pipeline_ir import _RecombineStep


class _IntrospectionMixin:
    __slots__ = ()

    def describe(self) -> str:
        """Render a biology-style narrative of this experiment.

        The output is one line per step, prefixed with its position
        in the chain. Locks, allele weights, NP-length distributions,
        SHM kernels, and per-corruption rates are all surfaced. Use
        this to sanity-check what a fluent chain actually encodes —
        if it doesn't read like an immunology protocol, the chain is
        too murky.

        Returns a multi-line string ending without a trailing newline.
        Safe to ``print(exp.describe())``.

        Example::

            >>> print(Experiment.on("human_igh").recombine().mutate(count=(5, 15)).describe())
            Experiment on human_igh (vdj, DataConfig)
              1. V(D)J recombination: sample V/D/J alleles; empirical exonuclease trim (V3', D5', D3', J5'); insert NP1 (0–11 weighted bases) and NP2 (0–11 weighted bases)
              2. Somatic hypermutation (S5F context model, human heavy-chain (HH_S5F)): 5–15 mutations/record
        """
        header = _describe_experiment_header(self._refdata, self._dataconfig)
        if not self._steps and not self._contracts:
            return header + "\n  (no steps appended yet)"

        lines = [header]
        resolved = self._steps_with_locks_resolved()
        body_lines = _describe_step_sequence(resolved, self._refdata.chain_type)
        lines.extend(body_lines)
        contracts_line = _format_declared_contracts(self._contracts)
        if contracts_line:
            lines.append(f"  Constraints: {contracts_line}")
        if self._metadata:
            stamps = ", ".join(f"{k}={v!r}" for k, v in self._metadata.items())
            lines.append(f"  Metadata stamped on every record: {stamps}")
        return "\n".join(lines)

    def _steps_with_locks_resolved(self) -> List[Any]:
        """Return a copy of ``self._steps`` with per-segment allele
        locks (from :meth:`restrict_alleles`) injected into the first
        :class:`_RecombineStep`. Compile-time injection lives in
        :meth:`_build_simulator`; ``describe()`` needs the same
        substitution to render lock info correctly without
        side-effecting the live step list."""
        from dataclasses import replace as _replace

        any_lock = any(self._locks[seg] is not None for seg in ("V", "D", "J"))
        if not any_lock:
            return list(self._steps)
        out: List[Any] = []
        injected = False
        for step in self._steps:
            if not injected and isinstance(step, _RecombineStep):
                step = _replace(
                    step,
                    locks_v=self._locks["V"],
                    locks_d=self._locks["D"],
                    locks_j=self._locks["J"],
                )
                injected = True
            out.append(step)
        return out

    def __repr__(self) -> str:
        return f"<Experiment chain={self.chain_type} steps={self.step_count}>"
