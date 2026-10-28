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


def test_stages_are_exactly_the_declared_order():
    """Developmental order, then the two that are not developmental stages:
    assay is what the protocol did to the molecule afterwards, and projection is
    the surface a benchmark reads rather than a step anything performs."""
    assert catalog.STAGES == ("germline", "recombination", "junction", "expression",
                              "maturation", "clonality", "assay", "projection")
    assert {m.stage for m in catalog.all_mechanisms()} <= set(catalog.STAGES)


def test_every_mechanism_states_a_summary_and_owns_truth_fields():
    thin = [m.id for m in catalog.all_mechanisms() if not m.summary.strip() or not m.truth_fields]
    assert thin == [], f"mechanisms with no summary or no truth fields: {thin}"


def test_every_truth_field_is_a_field_some_mode_actually_produces():
    """Measured across every run mode, not checked against a hand-written list.

    The earlier version exempted a list of fields "not in a plain record", and
    that list named `generation` -- which exists nowhere in GenAIRR; the field
    is `lineage_generation`. An invented name in the one place the check stops
    looking is the exact failure this test exists to prevent, so the exemption
    is gone and the known set is taken from real records.
    """
    from tests.support.fixtures.record_modes import every_known_field
    known = every_known_field()
    unknown = sorted({f for m in catalog.all_mechanisms() for f in m.truth_fields
                      if f not in known})
    assert unknown == [], f"truth fields no run mode produces: {unknown}"


def test_the_lineage_truth_fields_have_owners():
    """These are what a clonal-inference tool is scored against, and a single
    clonality mechanism owning `clone_id` left all five unowned."""
    from tests.support.fixtures.record_modes import fields_by_mode
    lineage_fields = {f for f in fields_by_mode()["clonal_lineage"] if f.startswith("lineage_")}
    owned = {f for m in catalog.all_mechanisms() for f in m.truth_fields}
    assert lineage_fields <= owned, f"unowned lineage truth: {sorted(lineage_fields - owned)}"


def test_infrastructure_passes_each_state_why_they_are_not_biology():
    unexplained = [p for p, reason in catalog.INFRASTRUCTURE.items() if not reason.strip()]
    assert unexplained == [], f"infrastructure passes with no stated reason: {unexplained}"


def test_by_id_finds_a_mechanism_and_rejects_an_unknown_one():
    assert catalog.by_id("junction.n_addition").stage == "junction"
    with pytest.raises(KeyError):
        catalog.by_id("junction.n_additon")


# ── the projection surface, and truth that is not exported ───────────────────

def test_the_alignment_surface_has_an_owner():
    """These are the fields an aligner benchmark actually compares.

    They were deliberately unowned on the reasoning that no single mechanism
    produces them -- a trim, an N addition, a mutation and an indel all land in
    the same CIGAR. True, but the consequence was that the catalog could read
    fully green while every coordinate a benchmark consumes was wrong. The
    answer is a cross-cutting mechanism with its own claims, not an exemption.
    """
    owned = {f for m in catalog.all_mechanisms() for f in m.truth_fields}
    surface = {
        "v_cigar", "d_cigar", "j_cigar",
        "v_sequence_start", "v_sequence_end", "d_sequence_start", "d_sequence_end",
        "j_sequence_start", "j_sequence_end",
        "v_alignment_start", "v_alignment_end", "d_alignment_start", "d_alignment_end",
        "j_alignment_start", "j_alignment_end",
        "sequence_alignment", "germline_alignment", "germline_alignment_d_mask",
        "sequence_aa",
    }
    assert surface <= owned, f"alignment-surface fields owned by nobody: {sorted(surface - owned)}"


def test_no_exported_field_is_owned_by_nobody():
    """Every field a record carries is either claimed by a mechanism or named
    in the explicit unclaimed record, with a reason. A field nobody owns is a
    field nobody is answerable for."""
    from tests.support.fixtures.record_modes import every_known_field
    owned = {f for m in catalog.all_mechanisms() for f in m.truth_fields}
    unaccounted = sorted(every_known_field() - owned - set(catalog.UNCLAIMED_FIELDS))
    assert unaccounted == [], (
        f"exported fields with no owner and no recorded reason: {unaccounted}")


def test_every_unclaimed_field_states_why():
    assert catalog.UNCLAIMED_FIELDS
    thin = [f for f, reason in catalog.UNCLAIMED_FIELDS.items() if not reason.strip()]
    assert thin == [], f"unclaimed fields with no reason: {thin}"


def test_a_mechanism_whose_truth_is_not_exported_says_so():
    """`germline.genotype` owned only `locus`, and `assay.ambiguous_calls` and
    `assay.duplex_formation` owned only `sequence`, because their real truth is
    not in the record at all. Reporting those as merely untested confuses "no
    one has written the test" with "there is nothing to test against"."""
    gapped = [m for m in catalog.all_mechanisms() if m.truth_gap]
    assert gapped, "no mechanism declares a truth gap, which cannot be right"
    for mechanism in gapped:
        assert mechanism.truth_gap.strip(), mechanism.id


def test_a_truth_gap_reports_as_its_own_status_not_as_unverified():
    rows = {r.mechanism_id: r for r in catalog.coverage_rows({})}
    gapped = [m.id for m in catalog.all_mechanisms() if m.truth_gap]
    for mechanism_id in gapped:
        assert rows[mechanism_id].status == "no-truth-exported", (
            f"{mechanism_id} declares a truth gap but reports "
            f"{rows[mechanism_id].status!r}")


def test_a_truth_gap_stays_visible_even_when_the_mechanism_verifies():
    """A mechanism can be verified on the one borrowed field it owns while the
    truth that actually matters is still not exported. assay.ambiguous_calls is
    exactly that: the claim "an N is not a mutation" can be established against
    the sequence, but no N-count field exists, so "verified" alone would hide
    the gap rather than report it."""
    gapped = {m.id for m in catalog.all_mechanisms() if m.truth_gap}
    rows = {r.mechanism_id: r for r in catalog.coverage_rows(
        {mid: ("x::y",) for mid in gapped})}
    for mechanism_id in gapped:
        assert rows[mechanism_id].truth_gap, (
            f"{mechanism_id} lost its truth gap once a test was registered")


def test_every_claimed_field_is_read_by_the_test_that_claims_it():
    """A registration must claim only fields its test actually inspects.

    This is the failure the whole scheme exists to prevent, and it happened
    anyway: seven registrations claimed eleven fields their bodies never read,
    including end_loss_3_length (zero in every record of its own fixture) and
    the six germline-coordinate record fields, whose values come from the
    column walker rather than the IR track the test inspects. Each produced a
    green or partial row on nothing.

    The decorator is excluded from the searched text, since it names every
    claimed field by construction -- which is why an earlier version of this
    analysis reported no offenders.
    """
    import ast
    import re

    from tests.support import paths

    offenders = {}
    for module in sorted(paths.tests_root().rglob("test_*.py")):
        if "__pycache__" in module.parts:
            continue
        tree = ast.parse(module.read_text(encoding="utf-8"))
        helpers = {node.name: "\n".join(ast.unparse(s) for s in node.body)
                   for node in tree.body
                   if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if not node.name.startswith("test_"):
                continue            # a nested sample used to prove a rejection
            for decorator in node.decorator_list:
                if not (isinstance(decorator, ast.Call)
                        and getattr(decorator.func, "id", "") == "verifies"):
                    continue
                claimed = next((ast.literal_eval(k.value) for k in decorator.keywords
                                if k.arg == "fields"), ())
                body = "\n".join(ast.unparse(s) for s in node.body)
                for name, text in helpers.items():
                    if name != node.name and re.search(rf"\b{name}\b", body):
                        body += "\n" + text
                missing = [f for f in claimed if f not in body]
                if missing:
                    offenders[f"{module.name}::{node.name}"] = missing
    assert offenders == {}, (
        "registrations claiming fields their test never reads:\n"
        + "\n".join(f"  {where}: {fields}" for where, fields in sorted(offenders.items()))
        + "\n\nClaim only what the body inspects. Read fields by name rather than "
          "through an f-string, so this check can see them.")
