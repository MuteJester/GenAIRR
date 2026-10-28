"""The canonical definition of the reproducibility fixtures.

The generator and the test that checks its output must agree exactly, or a
failure says "this platform diverged" when the truth is "the two scripts
disagree". They therefore share this one definition rather than each spelling
out the configuration.

The previous arrangement is the reason this exists: the fixtures were generated
by a script calling an API that no longer exists (`run_to_file`, `GenAIRR.ops`)
and were never read by any test, so they drifted to 56 columns of a schema the
product no longer produces without anything noticing.
"""
from __future__ import annotations

from pathlib import Path

#: Loci chosen because they diverge for different reasons: heavy chain with a D
#: segment, light chain without one, and a TCR locus whose packaged reference
#: needs curation allowed.
LOCI = ("human_igh", "human_igk", "human_tcrb")
SEED = 42
N = 10

#: A second heavy-chain fixture exercising the features that bypass the ordinary
#: rail: D inversion and receptor revision.
MIXED_LOCUS = "human_igh"
MIXED_N = 20


def golden_name(locus: str, *, n: int = N, seed: int = SEED) -> str:
    return f"seed{seed}_n{n}_{locus}.tsv"


MIXED_NAME = f"seed{SEED}_n{MIXED_N}_{MIXED_LOCUS}_mixed_features.tsv"


def _experiment(locus: str, *, mixed: bool = False):
    from tests.support.builders import experiment as experiment_builder
    exp = experiment_builder.recombined(locus=locus)
    if mixed:
        exp = exp.invert_d(prob=0.5).receptor_revision(prob=0.5)
    return exp


def write_reproducibility_tsv(path: Path, locus: str, *, n: int = N, seed: int = SEED,
                              mixed: bool = False) -> Path:
    """Write the fixture for ``locus`` to ``path`` and return it."""
    exp = _experiment(locus, mixed=mixed)
    exp.compile().run_records(n=n, seed=seed).to_tsv(str(path))
    return path


def every_fixture() -> tuple[tuple[str, str, int, bool], ...]:
    """``(filename, locus, n, mixed)`` for each fixture the suite checks."""
    rows = [(golden_name(locus), locus, N, False) for locus in LOCI]
    rows.append((MIXED_NAME, MIXED_LOCUS, MIXED_N, True))
    return tuple(rows)
