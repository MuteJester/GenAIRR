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
from ._stages import FIELDS_NOT_IN_A_PLAIN_RECORD, INFRASTRUCTURE, MECHANISMS, STAGES

__all__ = ["Mechanism", "CoverageRow", "RustSource", "STAGES", "INFRASTRUCTURE",
           "FIELDS_NOT_IN_A_PLAIN_RECORD", "PROJECTION_SOURCES", "RUST_SOURCES",
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

    Reads the registry so that a pass added to the engine and not classified
    here is a failure rather than an omission nobody sees. Any failure to read
    the registry propagates: a check that cannot see the engine must not report
    the catalog complete.
    """
    pass_registry()                      # fail closed if the engine cannot be read
    claimed = [p for m in MECHANISMS for p in m.engine_passes]
    claimed.extend(INFRASTRUCTURE)
    return tuple(claimed)
