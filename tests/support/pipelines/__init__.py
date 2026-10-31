"""Pipelines as data: drawn, applied, and printable as the code that made them.

A finding from a generated pipeline is worth nothing if nobody can reproduce
it, so a ``Pipeline`` is a value object first and an experiment second. It
renders as the exact chain that built it, which is what goes in a bug report
and what becomes a pinned regression test.
"""
from __future__ import annotations

from dataclasses import dataclass

from .registry import LOCI, PASSES, STAGES, Pass, by_name, for_locus

__all__ = ["LOCI", "PASSES", "STAGES", "Pass", "Pipeline", "by_name",
           "for_locus"]


def _render(value) -> str:
    """A literal that evaluates back to ``value``.

    Floats are rounded to six places: a drawn probability prints as
    ``0.8444218515250481`` otherwise, which makes a rendered pipeline unusable
    as something a person reads, and the sixth place changes no behaviour the
    engine has.
    """
    if isinstance(value, float):
        return repr(round(value, 6))
    if isinstance(value, tuple):
        return "(" + ", ".join(_render(item) for item in value) + ")"
    if isinstance(value, list):
        return "[" + ", ".join(_render(item) for item in value) + "]"
    return repr(value)


@dataclass(frozen=True)
class Pipeline:
    """A locus and an ordered list of ``(pass name, kwargs)``."""

    locus: str
    steps: tuple[tuple[str, dict], ...]

    def __post_init__(self) -> None:
        known = by_name()
        for name, _kwargs in self.steps:
            if name not in known:
                raise ValueError(f"{name!r} is not a declared pass")
        stages = [STAGES.index(known[name].stage) for name, _ in self.steps]
        if stages != sorted(stages):
            raise ValueError(
                "steps are out of stage order, so this pipeline does not mean "
                f"what it reads like: {[name for name, _ in self.steps]}")

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(name for name, _kwargs in self.steps)

    def has(self, name: str) -> bool:
        return name in self.names

    def kwargs_for(self, name: str) -> dict:
        for step, kwargs in self.steps:
            if step == name:
                return kwargs
        return {}

    def mechanisms(self) -> frozenset[str]:
        """Every catalog mechanism the passes in this pipeline exercise."""
        known = by_name()
        return frozenset(mechanism for name, _ in self.steps
                         for mechanism in known[name].mechanisms)

    def build(self, *, allow_curatable=None):
        """The Experiment this pipeline describes.

        Imports GenAIRR at call time rather than at module import, because the
        strategies here are also used by tests that only render pipelines.
        """
        from tests.support.builders import experiment as experiment_builder

        experiment = experiment_builder.recombined(
            locus=self.locus, allow_curatable=allow_curatable)
        for name, kwargs in self.steps:
            experiment = getattr(experiment, name)(**kwargs)
        return experiment

    def render(self) -> str:
        """The chain that built this, as code.

        What a failure message prints. Deliberately the builder's own entry
        point rather than ``Experiment.on(...)``, so a reader can paste it into
        this suite and get the same fixture, curatable-reference flag included.
        """
        lines = [f'(experiment_builder.recombined(locus="{self.locus}")']
        for name, kwargs in self.steps:
            arguments = ", ".join(f"{key}={_render(value)}"
                                  for key, value in sorted(kwargs.items()))
            lines.append(f" .{name}({arguments})")
        # Parenthesised, so the render is a single expression a reader can
        # paste and a test can compile. Without the wrapper a multi-line chain
        # is an IndentationError, which the render test caught.
        return "\n".join(lines) + ")"
