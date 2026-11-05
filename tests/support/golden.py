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

The same gap applied to the clonal path until the ``clonal_*`` variants were
added. Every fixture ran a single molecule, so **no fixture grew a lineage or
drew a clone size**, and that is where the engine's transcendental arithmetic
lives: ``lineage/poisson.rs:26`` (``exp``, every offspring draw),
``lineage/affinity.rs:209,221`` (``exp``, every affinity), and
``lineage/clone_size.rs:33,43,46,52`` (``ln``, ``sqrt``, ``cos``, ``powf``,
``exp``). Rust documents the precision of those as platform-dependent, and a
clone size is ``x.round()`` of a ``powf`` -- so one ulp is a different integer,
not a different last digit. Two determinism-version bumps passed through this
file moving nothing but the version column, which is the symptom: the fixtures
could not see the model they were stamping.
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
    # The clonal paths, which is where the engine's transcendental arithmetic
    # is. Kept small: these are a portability probe, not a biology fixture.
    #
    # A neutral lineage: one `exp` per parent per generation in the Poisson
    # offspring sampler, and the tree topology that follows from it.
    "clonal_lineage": (("clonal_lineage",
                        {"n_clones": 3, "max_generations": 4, "n_max": 50,
                         "n_sample": 5, "rate": 0.01, "lambda_base": 1.8}),),
    # The same with selection on, which adds one `exp` per cell for affinity
    # and makes the genealogy depend on its value.
    "clonal_lineage_selected": (("clonal_lineage",
                                 {"n_clones": 3, "max_generations": 4,
                                  "n_max": 50, "n_sample": 5, "rate": 0.01,
                                  "lambda_base": 1.8, "selection_strength": 0.5,
                                  "beta": 0.001, "mature_substitutions": 5}),),
    # Clone sizes through `powf`. The size is `x.round()` of it, so a one-ulp
    # difference is a different `duplicate_count`, not a rounding artefact.
    "clonal_power_law": (("clonal_repertoire",
                          {"n_clones": 4, "size_distribution": "power_law",
                           "exponent": 2.0, "max_size": 30}),),
    # Clone sizes through Box-Muller and `exp`: `ln`, `sqrt` and `cos` as well,
    # which `powf` alone does not reach.
    "clonal_lognormal": (("clonal_repertoire",
                          {"n_clones": 4, "size_distribution": "lognormal",
                           "mu": 1.0, "sigma": 1.0, "max_size": 30}),),
}

#: The variants above whose `run_records` takes no `n`.
CLONAL = frozenset({"clonal_lineage", "clonal_lineage_selected",
                    "clonal_power_law", "clonal_lognormal"})

#: Variants that fork the pipeline into clones. Their record count follows from
#: the clone structure rather than from a caller's ``n``, so they are named
#: without one and ``run_records`` is called without it.
#: ``(variant, locus, n)`` for every committed fixture. ``n`` is ``None`` for the
#: clonal variants, which do not take one.
FIXTURE_SPECS: tuple[tuple[str, str, int | None], ...] = (
    tuple(("", locus, N) for locus in LOCI)
    + (("mixed_features", "human_igh", 20),
       ("library_prep", "human_igh", 20),
       ("hypermutation", "human_igh", 20),
       ("clonal_lineage", "human_igh", None),
       ("clonal_lineage_selected", "human_igh", None),
       ("clonal_power_law", "human_igh", None),
       ("clonal_lognormal", "human_igh", None))
)


def golden_name(locus: str, *, n: int | None = N, seed: int = SEED,
                variant: str = "") -> str:
    suffix = f"_{variant}" if variant else ""
    size = f"_n{n}" if n is not None else ""
    return f"seed{seed}{size}_{locus}{suffix}.tsv"


def _experiment(locus: str, *, variant: str = ""):
    from tests.support.builders import experiment as experiment_builder
    exp = experiment_builder.recombined(locus=locus)
    for step, arguments in VARIANTS.get(variant, ()):
        exp = getattr(exp, step)(**arguments)
    return exp


def write_reproducibility_tsv(path: Path, locus: str, *, n: int | None = N,
                              seed: int = SEED, variant: str = "") -> Path:
    """Write the fixture for ``locus`` to ``path`` and return it."""
    exp = _experiment(locus, variant=variant)
    compiled = exp.compile()
    result = (compiled.run_records(seed=seed) if variant in CLONAL
              else compiled.run_records(n=n, seed=seed))
    result.to_tsv(str(path))
    return path


def every_fixture() -> tuple[tuple[str, str, int, str], ...]:
    """``(filename, locus, n, variant)`` for each fixture the suite checks."""
    return tuple(
        (golden_name(locus, n=n, variant=variant), locus, n, variant)
        for variant, locus, n in FIXTURE_SPECS
    )
