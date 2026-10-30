"""What the engine recorded about each base's origin.

A simulated sequence carries, per base, the germline position it came from --
or nothing, when no germline base produced it. That single fact is enough to
derive the junction's untemplated content from first principles, which is what
makes it an independent check on the recorded np fields.

This module reads a recorded ``Simulation``; it imports nothing from GenAIRR.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BaseOrigin:
    """One base of the final sequence, and where it came from."""

    index: int
    base: int
    germline_position: int | None

    @property
    def untemplated(self) -> bool:
        """True when no germline base produced this one (an N or P addition)."""
        return self.germline_position is None


@dataclass(frozen=True)
class Run:
    """A contiguous stretch of bases sharing an origin."""

    start: int
    bases: bytes

    @property
    def end(self) -> int:
        """One past the last index, so ``[start, end)`` matches region spans."""
        return self.start + len(self.bases)


def base_provenance(simulation) -> tuple[BaseOrigin, ...]:
    """Per-base origin for the whole final sequence, in sequence order."""
    bases = simulation.bases()
    return tuple(
        BaseOrigin(index=i, base=bases[i], germline_position=simulation.germline_position(i))
        for i in range(len(bases))
    )


def untemplated_runs(simulation) -> tuple[Run, ...]:
    """The contiguous stretches of untemplated bases, in sequence order.

    For an ordinary VDJ recombination these are the V-D and D-J junction
    insertions, so the result should equal the recorded ``np1`` and ``np2``.
    Derived here only from which bases lack a germline position.
    """
    runs: list[Run] = []
    current: list[BaseOrigin] = []
    for origin in base_provenance(simulation):
        if origin.untemplated:
            if current and origin.index != current[-1].index + 1:
                runs.append(_as_run(current))
                current = []
            current.append(origin)
        elif current:
            runs.append(_as_run(current))
            current = []
    if current:
        runs.append(_as_run(current))
    return tuple(runs)


def _as_run(origins: list[BaseOrigin]) -> Run:
    return Run(start=origins[0].index, bases=bytes(o.base for o in origins))


def changed_base_count(simulation) -> int:
    """How many positions carry a base different from their germline byte.

    Case-insensitive, because case is a separate channel: a templated base is
    stored lowercase and an SHM write emits uppercase, so a byte comparison
    counts every mutated position as changed even when the drawn base was the
    one already there.

    This is what ``n_mutations`` counts for a single mutation pass, so tests
    whose real subject is something else -- the burden denominator, founder
    evidence, landmark survival -- can derive the expected count instead of
    hardcoding one. They hardcoded the inflated value before the engine
    stopped counting substitutions that changed nothing, which is why they all
    moved at once.

    Not an identity and not a comparison against anything independent, so it
    carries no circularity entry: it reads one recorded fact and counts it.
    """
    bases, germline = simulation.bases(), simulation.germline_bases()
    return sum(1 for i in range(len(bases))
               if chr(bases[i]).upper() != chr(germline[i]).upper())
