"""What decides the composition of N-region bases.

Uniform is the intended default, not a gap. GenAIRR provides an empirical model
for human IGH, but a human IGH model would be wrong for other species and
chains, so it is opt-in rather than applied automatically -- a user trains,
fits or designs their own, or stays uniform. `_dataconfig_extract` documents the
cascade as typed -> uniform, and records that auto-lifting the legacy pickle
fields is deliberately deferred because it would silently change output bytes.

So the claims worth pinning are that the default is uniform, and that a
supplied typed model is honoured.
"""
from __future__ import annotations

import collections

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder

pytestmark = pytest.mark.behavior

BASES = "ACGT"


def _np_composition(exp, n=400, seed=13):
    counts = collections.Counter()
    for outcome in exp.run(n=n, seed=seed):
        record = outcome.airr_record(exp.refdata)
        for field in ("np1", "np2"):
            counts.update((record[field] or "").upper())
    total = sum(counts[b] for b in BASES)
    return {b: counts[b] / total for b in BASES}, total


@verifies("junction.n_addition",
          claim="with no NP base model configured, N-region bases are uniform over "
                "A/C/G/T -- the documented typed -> uniform cascade falling through "
                "to uniform, not an unconfigured accident",
          fields=("np1", "np2"))
def test_the_default_n_region_composition_is_uniform():
    composition, total = _np_composition(experiment_builder.recombined())
    assert total > 2000, f"only {total} NP bases sampled; widen the sample"
    for base, share in composition.items():
        assert abs(share - 0.25) < 0.03, (
            f"base {base} is {share:.3f} of N-region content, not ~0.25. "
            f"Full composition: {composition}. The default cascade should land "
            "on UniformBase; a biased default would mean a model was applied "
            "that the user did not ask for.")


@verifies("junction.n_addition",
          claim="a supplied typed NP base model is honoured, so a user can fit a "
                "model for their own species and chain rather than inheriting one "
                "calibrated for human IGH",
          fields=("np1", "np1_length"))
def test_a_supplied_empirical_base_model_is_honoured():
    """The capability that matters. If this did not hold, uniform would not be a
    default but a ceiling -- and an empirical model for human IGH would be the
    only composition anyone could ever simulate.
    """
    import copy

    import GenAIRR as ga
    from GenAIRR.reference_models import NpBaseModelSpec, ReferenceEmpiricalModels

    config = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    existing = getattr(config, "reference_models", None)
    skewed = {"NP1": NpBaseModelSpec(
        kind="empirical_first_base",
        first_base={"A": 0.01, "C": 0.01, "G": 0.97, "T": 0.01},
        transitions=None)}
    config.reference_models = ReferenceEmpiricalModels(
        np_lengths=getattr(existing, "np_lengths", None) or {},
        trims=getattr(existing, "trims", None) or {},
        np_bases=skewed,
    ) if isinstance(existing, ReferenceEmpiricalModels) else ReferenceEmpiricalModels(
        np_bases=skewed)

    exp = ga.Experiment.on(config).recombine()
    composition, total = _np1_composition(exp, n=200, seed=17)
    assert total > 400, f"only {total} NP1 bases sampled; widen the sample"
    assert composition["G"] > 0.80, (
        f"a 97% G first-base model produced G at {composition['G']:.3f}; the "
        f"supplied model was not applied. Composition: {composition}")


def _np1_composition(exp, *, n, seed):
    counts = collections.Counter()
    for outcome in exp.run(n=n, seed=seed):
        counts.update((outcome.airr_record(exp.refdata)["np1"] or "").upper())
    total = sum(counts[b] for b in BASES)
    return {b: counts[b] / total for b in BASES}, total
