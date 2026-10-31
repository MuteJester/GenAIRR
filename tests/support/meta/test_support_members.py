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
    """Oracle modules that report disagreements, so the rule applies to them.

    The scan itself lives in the oracles package: the circularity probes need
    the same list, and a test module may not import another test module.
    """
    from tests.support import oracles

    return oracles.reporting_modules()


def test_every_reporting_oracle_declares_what_its_identities_do_under_a_both_sides_probe():
    """An identity nobody probed is an identity nobody knows the strength of.

    What this enforces: the mapping exists, covers something, and every note
    carries a measurement rather than a reassurance.

    Declaring ``IDENTITIES`` is mandatory, and the mapping must cover exactly
    those keys. Two modules were found emitting identities absent from their
    mapping -- projection documented 6 of 13, coordinates 7 of 8 -- which is
    the blind spot the first version of this docstring described as a
    limitation, realized. A third, ``reconstruction``, declared no identities
    at all and said so in a comment, which left two of its comparisons
    unprobed and the mapping unable to notice.

    Four of the oracles additionally validate every reported identity against
    the tuple at construction, which is stronger than this test: it catches an
    identity invented at the call site. ``reconstruction``'s result types carry
    no identity field, so it declares ``REPORTS`` instead and
    ``test_an_oracle_without_an_identity_field_attributes_every_function``
    enforces that.
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
        declared = getattr(module, "IDENTITIES", None)
        assert isinstance(declared, tuple) and declared, (
            f"{name} reports disagreements but declares no IDENTITIES tuple, "
            "so nothing can check that CIRCULARITY covers what it reports. "
            "Name every identity the module can report.")
        assert set(mapping) == set(declared), (
            f"{name}: CIRCULARITY and IDENTITIES disagree.\n"
            f"  emitted but unprobed: {sorted(set(declared) - set(mapping))}\n"
            f"  probed but not emitted: {sorted(set(mapping) - set(declared))}")
        for identity, note in mapping.items():
            assert identity.strip(), f"{name}: an identity with no name"
            assert len(note.split()) >= 12, (
                f"{name}/{identity!r}: the note is too short to be a "
                f"measurement: {note!r}")
            # Any digit satisfies this, including an incidental "3 loci". It
            # catches a note with no measurement at all, which is what five of
            # the first drafts were; it cannot tell a measurement from a
            # mention.
            assert re.search(r"\d", note), (
                f"{name}/{identity!r}: the note carries no number, so no "
                f"probe was measured: {note!r}")
            # "NOT CIRCULAR" already contains "CIRCULAR", so the alternation
            # is for readability rather than matching. The load-bearing check
            # is the startswith() in the test below; a note phrased "This
            # identity is CIRCULAR..." would satisfy this one and evade that,
            # which is a known weakness rather than a guarded one.
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


def _identity_carrying_result(module):
    """The module's result dataclass that carries an ``identity`` field, if any."""
    import dataclasses

    for value in vars(module).values():
        if (dataclasses.is_dataclass(value) and isinstance(value, type)
                and any(f.name == "identity" for f in dataclasses.fields(value))):
            return value
    return None


def test_an_identity_carrying_oracle_refuses_an_identity_it_never_declared():
    """The construction guard, asserted rather than assumed.

    Four oracles validate ``identity`` in ``__post_init__``. That guard is
    what makes their IDENTITIES tuple load-bearing instead of decorative, and
    nothing was checking that it was still there -- the projection oracle
    carried an f-string identity for months, reachable only on a reverse
    complemented record, which is how the gap was found.
    """
    import dataclasses

    probed = []
    for name, module in sorted(_reporting_oracles().items()):
        result = _identity_carrying_result(module)
        if result is None:
            continue
        probed.append(name)
        fields = {f.name: 0 if f.type in ("int", int) else ""
                  for f in dataclasses.fields(result)}
        fields["identity"] = "an identity nobody declared"
        with pytest.raises(ValueError, match="declared identity"):
            result(**fields)
    assert len(probed) >= 3, (
        f"only {probed} were probed; the scan is not finding the result types")


def test_an_oracle_without_an_identity_field_attributes_every_function():
    """Where there is no construction site to guard, every function says what
    it reports.

    ``reconstruction`` returns positions, bare strings and tuples, so no
    ``__post_init__`` can catch an identity invented at the call site. Its
    ``REPORTS`` mapping is the substitute, and this makes it exhaustive: a new
    public comparison fails here until its author attributes it, and the
    attribution has to name a declared -- therefore probed -- identity.

    The previous arrangement was a comment saying coverage could not be
    enforced for this module, under which two of its four comparisons went
    unprobed.
    """
    import inspect

    checked = []
    for name, module in sorted(_reporting_oracles().items()):
        if _identity_carrying_result(module) is not None:
            continue
        checked.append(name)
        reports = getattr(module, "REPORTS", None)
        assert isinstance(reports, dict) and reports, (
            f"{name} has no identity-carrying result type and no REPORTS "
            "mapping, so nothing relates its comparisons to its identities")
        public = {member for member, value in vars(module).items()
                  if inspect.isfunction(value) and not member.startswith("_")
                  and value.__module__ == module.__name__}
        assert set(reports) == public, (
            f"{name}: REPORTS does not account for every public function.\n"
            f"  unattributed: {sorted(public - set(reports))}\n"
            f"  named but absent: {sorted(set(reports) - public)}")
        declared = set(module.IDENTITIES)
        attributed = {identity for value in reports.values() for identity in value}
        assert attributed <= declared, (
            f"{name}: REPORTS names identities that are not declared and so "
            f"are not probed: {sorted(attributed - declared)}")
        assert attributed == declared, (
            f"{name}: declared identities that no function reports: "
            f"{sorted(declared - attributed)}")
    assert checked, (
        "no oracle lacks an identity-carrying result type, so this test is "
        "vacuous -- delete it or restore the case it was written for")
