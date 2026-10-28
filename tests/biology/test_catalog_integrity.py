"""The catalog must describe the engine that exists, not the one it remembers.

A catalog of biological mechanisms is only worth having if it cannot quietly
fall out of step with the simulation. These tests are what stop it becoming a
comfortable fiction: every pass the engine defines is accounted for exactly
once, every truth field named is a field a record actually carries, and a
missing registry fails loudly instead of skipping.
"""
from __future__ import annotations

import pytest

from GenAIRR._compiled import pass_registry
from tests.biology import catalog

pytestmark = pytest.mark.pin


def test_every_engine_pass_is_classified_exactly_once():
    """Biology or infrastructure, once each.

    Not every pass is biological: the engine documents EchoPass as its
    deterministic transform reference and SampleBasePass as its sampling
    reference. Forcing those into a biological mechanism would make the catalog
    assert something false, so infrastructure is a classification of its own
    rather than a gap.
    """
    claimed = catalog.classified_passes()
    duplicates = sorted({p for p in claimed if list(claimed).count(p) > 1})
    assert duplicates == [], f"passes classified more than once: {duplicates}"
    registry = set(pass_registry())
    assert set(claimed) == registry, (
        f"unclassified passes: {sorted(registry - set(claimed))}\n"
        f"classified but not in the engine: {sorted(set(claimed) - registry)}\n"
        "A new pass must be classified -- as a mechanism, or as infrastructure "
        "with a reason -- before this passes.")


def test_an_unavailable_registry_fails_rather_than_skips(monkeypatch):
    """Fail closed. A check that skips when it cannot see the engine reports the
    catalog complete on no evidence, which is worse than no check."""
    def boom():
        raise RuntimeError("registry unavailable")
    monkeypatch.setattr(catalog, "pass_registry", boom)
    with pytest.raises(RuntimeError):
        catalog.classified_passes()


def test_mechanism_ids_are_unique():
    ids = [m.id for m in catalog.all_mechanisms()]
    assert sorted(ids) == sorted(set(ids))


def test_every_mechanism_id_is_prefixed_with_its_stage():
    wrong = [m.id for m in catalog.all_mechanisms() if not m.id.startswith(f"{m.stage}.")]
    assert wrong == [], f"mechanism ids not prefixed with their stage: {wrong}"


def test_stages_are_exactly_the_declared_developmental_order():
    assert catalog.STAGES == ("germline", "recombination", "junction", "expression",
                              "maturation", "clonality", "assay")
    assert {m.stage for m in catalog.all_mechanisms()} <= set(catalog.STAGES)


def test_every_mechanism_states_a_summary_and_owns_truth_fields():
    thin = [m.id for m in catalog.all_mechanisms() if not m.summary.strip() or not m.truth_fields]
    assert thin == [], f"mechanisms with no summary or no truth fields: {thin}"


def test_every_truth_field_is_a_field_a_record_carries():
    """A typo here would make the catalog claim a mechanism owns a field that
    does not exist, and no test would notice."""
    import GenAIRR as ga
    exp = (ga.Experiment.on("human_igh").recombine().mutate(model="uniform", count=2)
           .polymerase_indels(count=1).sequencing_errors(count=1))
    record = exp.run(n=1, seed=11)[0].airr_record(exp.refdata)
    known = set(record) | catalog.FIELDS_NOT_IN_A_PLAIN_RECORD
    unknown = sorted({f for m in catalog.all_mechanisms() for f in m.truth_fields
                      if f not in known})
    assert unknown == [], f"truth fields no record carries: {unknown}"


def test_infrastructure_passes_each_state_why_they_are_not_biology():
    unexplained = [p for p, reason in catalog.INFRASTRUCTURE.items() if not reason.strip()]
    assert unexplained == [], f"infrastructure passes with no stated reason: {unexplained}"


def test_by_id_finds_a_mechanism_and_rejects_an_unknown_one():
    assert catalog.by_id("junction.n_addition").stage == "junction"
    with pytest.raises(KeyError):
        catalog.by_id("junction.n_additon")
