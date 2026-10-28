"""Compiled/runtime experiment wrappers, one module per runtime shape."""
from GenAIRR._engine import pass_registry as _pass_registry

from .plain import CompiledExperiment
from .clonal import CompiledClonalExperiment
from .lineage import CompiledLineageExperiment
from .repertoire import CompiledRepertoireExperiment

def pass_registry() -> tuple[str, ...]:
    """Every pass the engine defines, sorted, with no duplicates.

    ``Outcome.pass_names()`` reports only what a given compiled plan
    scheduled; this is the whole set, so a pass nobody schedules is still
    visible to a completeness check.
    """
    return tuple(_pass_registry())


__all__ = [
    "CompiledExperiment",
    "CompiledClonalExperiment",
    "CompiledLineageExperiment",
    "CompiledRepertoireExperiment",
    "pass_registry",
]
