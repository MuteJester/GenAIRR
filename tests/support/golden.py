"""The canonical definition of the reproducibility fixtures.

The generator and the test that checks its output must agree exactly, or a
failure says "this platform diverged" when the truth is "the two scripts
disagree". They therefore share this one definition rather than each spelling
out the configuration.

The previous arrangement is the reason this exists: the fixtures were generated
by a script calling an API that no longer exists (`run_to_file`, `GenAIRR.ops`)
and were never read by any test, so they drifted to 56 columns of a schema the
product no longer produces without anything noticing.

**What the fixtures have to cover.** The cross-platform claim is only made
about code the fixtures run. For a long time every fixture was a bare
`recombined()` pipeline, so the claim covered allele sampling, trimming and NP
generation -- and said nothing about the mutation and corruption passes, which
is where the floating-point draws are. That gap was invisible from inside:
changing the replacement-base model for every substitution pass in the engine
moved no golden byte, because no golden drew a substitution. Each variant below
exists to run code the others do not.
"""
from __future__ import annotations

from pathlib import Path

SEED = 42

#: Loci chosen because they diverge for different reasons: heavy chain with a D
#: segment, light chain without one, and a TCR locus whose packaged reference
#: needs curation allowed.
LOCI = ("human_igh", "human_igk", "human_tcrb")
N = 10

#: Fixture variants: the builder calls layered onto a plain `recombined()`
#: experiment, in order, keyed by the name that appears in the filename. The
#: empty name is the plain pipeline.
#:
#: Counts rather than rates throughout, so every fixture draws the same number
#: of events on every platform regardless of how long the molecule came out.
VARIANTS: dict[str, tuple[tuple[str, dict], ...]] = {
    # Recombination features that bypass the ordinary rail.
    "mixed_features": (("invert_d", {"prob": 0.5}),
                       ("receptor_revision", {"prob": 0.5})),
    # The library-prep stack: the uniform substitution kernel, both error
    # passes that share its code path, and base masking. These are the
    # weighted floating-point draws, which are the ones most likely to differ
    # between platforms and were covered by nothing.
    "library_prep": (("mutate", {"model": "uniform", "count": 6}),
                     ("pcr_amplify", {"count": 4}),
                     ("sequencing_errors", {"count": 4}),
                     ("ambiguous_base_calls", {"count": 2})),
    # The default mutation model, which is what `mutate()` means when nobody
    # names one -- a different kernel from `uniform`, with its own
    # five-mer lookups. Pins the bundled kernel's numbers as well as the code.
    "hypermutation": (("mutate", {"count": 12}),),
}

#: ``(variant, locus, n)`` for every committed fixture.
FIXTURE_SPECS: tuple[tuple[str, str, int], ...] = (
    tuple(("", locus, N) for locus in LOCI)
    + (("mixed_features", "human_igh", 20),
       ("library_prep", "human_igh", 20),
       ("hypermutation", "human_igh", 20))
)


def golden_name(locus: str, *, n: int = N, seed: int = SEED, variant: str = "") -> str:
    suffix = f"_{variant}" if variant else ""
    return f"seed{seed}_n{n}_{locus}{suffix}.tsv"


def _experiment(locus: str, *, variant: str = ""):
    from tests.support.builders import experiment as experiment_builder
    exp = experiment_builder.recombined(locus=locus)
    for step, arguments in VARIANTS.get(variant, ()):
        exp = getattr(exp, step)(**arguments)
    return exp


def write_reproducibility_tsv(path: Path, locus: str, *, n: int = N, seed: int = SEED,
                              variant: str = "") -> Path:
    """Write the fixture for ``locus`` to ``path`` and return it."""
    exp = _experiment(locus, variant=variant)
    exp.compile().run_records(n=n, seed=seed).to_tsv(str(path))
    return path


def every_fixture() -> tuple[tuple[str, str, int, str], ...]:
    """``(filename, locus, n, variant)`` for each fixture the suite checks."""
    return tuple(
        (golden_name(locus, n=n, variant=variant), locus, n, variant)
        for variant, locus, n in FIXTURE_SPECS
    )
