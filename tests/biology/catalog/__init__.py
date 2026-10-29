"""What biology GenAIRR claims to simulate, and which of those claims is verified.

Without this, "which mechanisms have verified truth and which do not" is not a
question anyone can answer: you would have to read the suite and infer it. The
catalog makes the claim explicit and the gap visible.

Coverage is claimed by registration, not by location. A test declares which
mechanism it verifies wherever it lives, so a test under components/ counts the
same as one here, and verification that already exists on the Rust side counts
too. Filing biology tests in a parallel tree would have made the report depend
on a directory rather than on what was actually proved.
"""
from __future__ import annotations

from GenAIRR._compiled import pass_registry

from ._model import Mechanism
from ._registry import CoverageRow, coverage_rows, verifies
from ._rust import (PROJECTION_SOURCES, RUST_SOURCES, RustSource, all_rust_sources,
                    rust_covered_fields, rust_coverage)
from ._stages import (INFRASTRUCTURE, MECHANISMS, SHARED_FIELDS, STAGES,
                      UNCLAIMED_FIELDS)

__all__ = ["Mechanism", "CoverageRow", "RustSource", "STAGES", "INFRASTRUCTURE",
           "PROJECTION_SOURCES", "RUST_SOURCES", "SHARED_FIELDS", "UNCLAIMED_FIELDS",
           "all_mechanisms", "all_rust_sources", "by_id", "classified_passes",
           "coverage_rows", "pass_registry", "rust_coverage",
           "rust_covered_fields", "verifies"]

_BY_STAGE = {stage: i for i, stage in enumerate(STAGES)}


def all_mechanisms() -> tuple[Mechanism, ...]:
    """Every declared mechanism, in developmental order."""
    return tuple(sorted(MECHANISMS, key=lambda m: (_BY_STAGE[m.stage], m.id)))


def by_id(mechanism_id: str) -> Mechanism:
    """The mechanism with this id, or KeyError.

    Raising is deliberate: a typo in a registration must surface, not silently
    count as nothing.
    """
    for mechanism in MECHANISMS:
        if mechanism.id == mechanism_id:
            return mechanism
    raise KeyError(
        f"unknown mechanism {mechanism_id!r}. Known ids: "
        f"{sorted(m.id for m in MECHANISMS)}")


def classified_passes() -> tuple[str, ...]:
    """Every pass the catalog accounts for, as a mechanism or as infrastructure.

    This reads the catalog, not the engine: the returned names come from
    ``MECHANISMS`` and ``INFRASTRUCTURE``, so a pass the engine defines and
    nobody classified does not appear here and cannot be missed by looking at
    the result alone. The registry call is a fail-closed guard and nothing
    more -- a check that cannot read the engine must not report the catalog
    complete.

    What makes an unclassified pass a failure is the pair of tests either side
    of the registry, and neither is optional:

    - ``test_the_registry_names_every_pass_the_engine_defines`` compares the
      registry against every ``impl Pass`` in the engine's own source, so a
      pass added to the crate cannot stay out of the registry.
    - ``test_every_engine_pass_is_classified_exactly_once`` compares the
      registry against this function, so a pass in the registry cannot stay
      out of the catalog.
    """
    pass_registry()                      # fail closed if the engine cannot be read
    claimed = [p for m in MECHANISMS for p in m.engine_passes]
    claimed.extend(INFRASTRUCTURE)
    return tuple(claimed)
