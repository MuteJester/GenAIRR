"""The support layer's own contracts.

Each support package ships at least one working member, so later work extends
a pattern instead of inventing one. These tests pin what those members promise.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from tests.support.asserts import records as record_asserts
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import provenance

TESTS_ROOT = Path(__file__).resolve().parents[2]


def _imports(path: Path):
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.ImportFrom) and node.module:
            yield node.module
        elif isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name


def _modules(root: Path):
    return sorted(p for p in root.rglob("*.py") if "__pycache__" not in p.parts)


# ── builders ─────────────────────────────────────────────────────────────────

def test_recombined_builds_a_runnable_experiment():
    outcomes = experiment_builder.recombined().run(n=1, seed=0)
    assert len(outcomes) == 1


def test_recombined_passes_keywords_through_to_recombine():
    exp = experiment_builder.recombined(locus="human_igk")
    assert exp.chain_type == "vj"


# ── oracles ──────────────────────────────────────────────────────────────────

def test_untemplated_runs_recover_the_recorded_np_fields():
    """The oracle's whole claim: np1 and np2 are exactly the bases the engine
    marked as having no germline position.

    Derived from per-base provenance alone -- not by calling the production
    junction builder, and not by re-aligning anything.
    """
    exp = experiment_builder.recombined()
    outcome = exp.run(n=1, seed=0)[0]
    record = outcome.airr_record(exp.refdata)
    runs = provenance.untemplated_runs(outcome.final_simulation())
    assert [r.bases.decode() for r in runs] == [record["np1"], record["np2"]]
    assert [len(r.bases) for r in runs] == [record["np1_length"], record["np2_length"]]


def test_untemplated_runs_are_contiguous_and_ordered():
    exp = experiment_builder.recombined()
    runs = provenance.untemplated_runs(exp.run(n=1, seed=0)[0].final_simulation())
    assert [r.start for r in runs] == sorted(r.start for r in runs)
    for run in runs:
        assert len(run.bases) > 0


def test_templated_bases_carry_a_germline_position():
    exp = experiment_builder.recombined()
    sim = exp.run(n=1, seed=0)[0].final_simulation()
    for entry in provenance.base_provenance(sim):
        assert (entry.germline_position is None) == entry.untemplated


def test_oracles_do_not_import_genairr():
    """An oracle derives truth from the provenance it is handed.

    Importing GenAIRR would let it call the very code it exists to check, and
    a junction oracle that calls the production junction builder proves nothing.
    """
    offenders = {}
    for module in _modules(TESTS_ROOT / "support" / "oracles"):
        bad = sorted(m for m in _imports(module) if m == "GenAIRR" or m.startswith("GenAIRR."))
        if bad:
            offenders[module.name] = bad
    assert offenders == {}, f"oracles importing the code under test: {offenders}"


# ── asserts ──────────────────────────────────────────────────────────────────

def test_field_mismatch_names_the_field_the_claim_and_both_values():
    with pytest.raises(AssertionError) as caught:
        record_asserts.assert_field(
            {"np1": "TCTG"}, "np1", "GGG", because="the junction records untemplated bases")
    message = str(caught.value)
    assert "np1" in message and "TCTG" in message and "GGG" in message
    assert "the junction records untemplated bases" in message


def test_field_match_is_silent():
    record_asserts.assert_field({"np1": "TCTG"}, "np1", "TCTG", because="whatever")


# ── strategies ───────────────────────────────────────────────────────────────

def test_locus_strategy_only_offers_loci_the_suite_can_build():
    from hypothesis import given, settings

    from tests.support.strategies import experiments as experiment_strategies
    seen = set()

    @given(experiment_strategies.loci())
    @settings(max_examples=10, deadline=None)
    def collect(locus):
        seen.add(locus)

    collect()
    assert seen and seen <= set(experiment_strategies.LOCI)


def test_builders_do_not_import_hypothesis():
    """Strategies belong in support/strategies; builders stay plain."""
    offenders = [p.name for p in _modules(TESTS_ROOT / "support" / "builders")
                 if any(m.split(".")[0] == "hypothesis" for m in _imports(p))]
    assert offenders == []


# ── the circularity rule ─────────────────────────────────────────────────────
# The third rule in the oracles package docstring. It exists because the same
# over-claim was made twice: the reconstruction oracle's coordinate
# circularity was found, written down, and then the projection oracle asserted
# the opposite of it. Both times a review caught it after the claim was
# committed.

import re  # noqa: E402


def _reporting_oracles():
    """Oracle modules that report disagreements, so the rule applies to them."""
    import importlib

    root = TESTS_ROOT / "support" / "oracles"
    found = {}
    for module in sorted(root.glob("*.py")):
        if module.name == "__init__.py":
            continue
        source = module.read_text(encoding="utf-8")
        if "class Disagreement" not in source and "class Mismatch" not in source:
            continue
        found[module.stem] = importlib.import_module(
            f"tests.support.oracles.{module.stem}")
    return found


def test_every_reporting_oracle_declares_what_its_identities_do_under_a_both_sides_probe():
    """An identity nobody probed is an identity nobody knows the strength of.

    What this enforces: the mapping exists, covers something, and every note
    carries a measurement rather than a reassurance.

    What it does NOT enforce, stated because this suite's whole subject is
    claims that outrun their checks: the keys are the author's own names for
    the identities, and nothing here ties them to the identity strings the
    module can actually emit. A module could document six identities and emit
    a seventh. That judgement is still the author's -- which is exactly the
    judgement that was wrong twice -- so a reviewer should check the list
    rather than trust this test.
    """
    oracles = _reporting_oracles()
    assert len(oracles) >= 3, (
        f"only found {sorted(oracles)}; the scan is not reading the oracles")
    for name, module in sorted(oracles.items()):
        mapping = getattr(module, "CIRCULARITY", None)
        assert isinstance(mapping, dict) and mapping, (
            f"{name} reports disagreements but declares no CIRCULARITY "
            "mapping. Perturb both sides of each identity, measure what the "
            "oracle says, and record it -- see the oracles package docstring.")
        for identity, note in mapping.items():
            assert identity.strip(), f"{name}: an identity with no name"
            assert len(note.split()) >= 12, (
                f"{name}/{identity!r}: the note is too short to be a "
                f"measurement: {note!r}")
            assert re.search(r"\d", note), (
                f"{name}/{identity!r}: the note carries no number, so no "
                f"probe was measured: {note!r}")
            assert re.search(r"\b(CIRCULAR|NOT CIRCULAR)\b", note), (
                f"{name}/{identity!r}: the note does not say whether the "
                "identity is circular, which is the thing being recorded")


def test_the_circularity_notes_do_not_all_claim_independence():
    """A mapping that says every identity is independent is the shape of the
    defect, not evidence against it.

    Three of the oracles' identities are genuinely circular and say so. If
    that ever became zero, either the oracles were rewritten or somebody
    softened the notes.
    """
    oracles = _reporting_oracles()
    circular = [f"{name}/{identity}"
                for name, module in oracles.items()
                for identity, note in getattr(module, "CIRCULARITY", {}).items()
                if note.lstrip().startswith("CIRCULAR")]
    assert len(circular) >= 3, (
        f"only {len(circular)} identities are recorded as circular: "
        f"{circular}. The coordinate comparisons in reconstruction and "
        "projection are both circular on the coordinate, measured, and the "
        "cigar span identities are true by construction.")
