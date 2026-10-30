"""What a recorded mutation count must mean.

A tool scored against GenAIRR truth can only recover substitutions that are
visible in the bases. If the recorded count includes events that left the
sequence unchanged, every such tool carries a built-in false-negative rate
equal to the discrepancy -- and nothing in the record tells it so.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior

MODELS = ("uniform", "s5f")


@verifies("maturation.hypermutation",
          claim="every counted mutation changed a base, so n_mutations counts "
                "events that actually happened rather than draws that were "
                "made. At one mutation per record it equals the number of "
                "positions differing from their germline byte exactly, for "
                "both models. What it does NOT promise is that a caller can "
                "recover them all: a mutation on an untemplated junction base "
                "has no germline base to differ from, which is intended and "
                "measured in the test below",
          fields=("n_mutations",))
@pytest.mark.parametrize("model", MODELS)
def test_every_counted_mutation_changed_a_base(model):
    """A substitution to the same base is not a mutation under any molecular
    definition: AID deamination repair yields transitions and transversions,
    never a base replaced by itself.

    The engine used to count them. The replacement is drawn from all four
    bases including the original, so roughly a quarter of uniform draws landed
    on the base already present; those writes still happened and were still
    counted. They are now subtracted at the point of the write, which leaves
    the draw, the RNG consumption, the recorded trace and the sequence exactly
    as they were -- only the count moved.

    Compared case-insensitively on purpose: a templated base is stored
    lowercase and an SHM write emits uppercase, so a byte comparison sees a
    change at every mutated position and hides this entirely. The first
    attempt at the fix made exactly that mistake and corrected 2 of 96 no-ops.
    """
    exp = experiment_builder.recombined().mutate(model=model, count=1)
    unchanged = []
    for index, outcome in enumerate(exp.run(n=200, seed=0)):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        bases, germline = simulation.bases(), simulation.germline_bases()
        changed = sum(1 for i in range(len(bases))
                      if chr(bases[i]).upper() != chr(germline[i]).upper())
        if record["n_mutations"] and not changed:
            unchanged.append(index)
        # At one mutation per record there are no repeat hits, so the count
        # and the number of changed positions are the same number.
        assert record["n_mutations"] == changed, (
            f"{model}: record {index} reports {record['n_mutations']} mutations "
            f"with {changed} positions changed")
    assert unchanged == [], (
        f"{model}: {len(unchanged)}/200 records report a mutation that changed "
        "no base at all")


@verifies("maturation.hypermutation",
          claim="a mutation on an untemplated junction base is counted and is "
                "not recoverable from the sequence, which is intended rather "
                "than a defect: there is no germline base at that position to "
                "differ from. Measured at 4% of uniform and 4.5% of s5f "
                "single-mutation records, and it is the whole of the gap "
                "between n_mutations and what a mutation caller can see",
          fields=("n_mutations",))
@pytest.mark.parametrize("model", MODELS)
def test_a_mutation_on_an_untemplated_base_is_counted_and_unrecoverable(model):
    """The part of the gap that is not a defect, pinned so it is not mistaken
    for one again.

    `observable_substitutions` compares only templated positions, because an
    untemplated one has no germline coordinate. A mutation there changes a
    base, so it belongs in n_mutations, and no caller comparing the sequence
    to its germline can recover it. Both statements are true at once.
    """
    exp = experiment_builder.recombined().mutate(model=model, count=1)
    invisible = 0
    checked = 0
    for outcome in exp.run(n=200, seed=0):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        if not record["n_mutations"]:
            continue
        checked += 1
        observable = reconstruction.observable_substitutions(simulation)
        if observable:
            continue
        invisible += 1
        # Every such record must be an untemplated hit, not a no-op: a
        # counted mutation with no change anywhere is the defect the test
        # above forbids.
        bases, germline = simulation.bases(), simulation.germline_bases()
        changed = [i for i in range(len(bases))
                   if chr(bases[i]).upper() != chr(germline[i]).upper()]
        assert changed, f"{model}: a counted mutation changed nothing"
        assert all(simulation.germline_position(i) is None for i in changed), (
            f"{model}: a templated change was invisible to the oracle, which "
            "would mean the oracle is wrong rather than the position untemplated")
    assert checked > 150, f"{model}: only {checked} records carried a mutation"
    assert 0 < invisible < checked * 0.15, (
        f"{model}: {invisible}/{checked} counted mutations were unrecoverable; "
        "expected a few per cent of untemplated hits")


@verifies("maturation.hypermutation",
          claim="the recorded per-segment counts sum to the recorded total, so the "
                "attribution accounts for every event",
          fields=("n_v_mutations", "n_d_mutations", "n_j_mutations", "n_np_mutations"))
@pytest.mark.parametrize("model", MODELS)
def test_per_segment_counts_sum_to_the_total(model):
    exp = experiment_builder.recombined().mutate(model=model, count=4)
    for outcome in exp.run(n=50, seed=2):
        record = outcome.airr_record(exp.refdata)
        parts = sum(record[f] for f in ("n_v_mutations", "n_d_mutations",
                                        "n_j_mutations", "n_np_mutations"))
        assert parts == record["n_mutations"], (
            f"{model}: per-segment counts sum to {parts} but n_mutations is "
            f"{record['n_mutations']}")


@verifies("maturation.hypermutation",
          claim="n_observable_mutations is what a mutation caller can recover: "
                "the templated positions whose base differs from the germline "
                "base recorded there. Checked against an oracle that derives "
                "the same set independently from the per-base provenance, and "
                "bounded above by n_mutations, which additionally counts "
                "mutations on untemplated bases and repeat hits that restored "
                "the germline base",
          fields=("n_observable_mutations",))
@pytest.mark.parametrize("count", [1, 20])
@pytest.mark.parametrize("model", MODELS)
def test_n_observable_mutations_is_what_a_caller_can_recover(model, count):
    """The field exists because `n_mutations` cannot be both things.

    A mutation on an untemplated junction base changes a base, so it belongs
    in `n_mutations`, and has no germline base to differ from, so nothing can
    recover it. Rather than redefine the count, the recoverable number is
    exported alongside it. Scoring a caller against the wrong one of the two
    penalises it for substitutions that left no evidence.
    """
    exp = experiment_builder.recombined().mutate(model=model, count=count)
    totals = {"recorded": 0, "reported": 0, "oracle": 0}
    for outcome in exp.run(n=60, seed=101):
        record = outcome.airr_record(exp.refdata)
        expected = len(reconstruction.observable_substitutions(outcome.final_simulation()))
        assert record["n_observable_mutations"] == expected, (
            f"{model}/{count}: reports {record['n_observable_mutations']} "
            f"recoverable where the oracle derives {expected}")
        assert record["n_observable_mutations"] <= record["n_mutations"], (
            "more recoverable than counted, which cannot happen")
        totals["recorded"] += record["n_mutations"]
        totals["reported"] += record["n_observable_mutations"]
        totals["oracle"] += expected
    assert totals["oracle"] > 0, f"{model}/{count}: nothing was mutated"
    assert totals["reported"] == totals["oracle"]
    # The gap is the point of the field: if it were zero the two counts would
    # be the same number and one of them would be redundant.
    assert totals["reported"] < totals["recorded"], (
        f"{model}/{count}: no gap between {totals['recorded']} counted and "
        f"{totals['reported']} recoverable, so this field establishes nothing")


# Each entry is a pipeline that writes bases for a reason that is not somatic
# hypermutation: PCR misincorporation, a sequencing miscall, an ambiguous base
# call masking a position with N, or wholesale replacement by a background
# molecule. None of them is a mutation a caller should be scored on.
NON_SHM_WRITERS = (
    ("pcr", lambda e: e.pcr_amplify(count=5)),
    ("sequencing", lambda e: e.sequencing_errors(count=5)),
    ("ambiguous", lambda e: e.ambiguous_base_calls(count=4)),
    ("contaminant", lambda e: e.contaminate(prob=1.0)),
)


@verifies("maturation.hypermutation",
          claim="n_observable_mutations is attributable to somatic "
                "hypermutation and to nothing else. A pipeline that writes "
                "bases for a library or sequencing reason -- PCR "
                "misincorporation, a sequencing miscall, an N mask, a "
                "background contaminant -- moves no part of it, so a caller "
                "scored against it is not credited for calling noise",
          fields=("n_observable_mutations",))
@pytest.mark.parametrize("label,corrupt", NON_SHM_WRITERS,
                         ids=[n for n, _ in NON_SHM_WRITERS])
def test_observable_mutations_exclude_every_non_shm_writer(label, corrupt):
    """The field is published as the number to score a mutation caller
    against. Derived from the finished molecule alone it cannot be that: every
    pass that writes a base leaves a position differing from its germline, and
    a count of those differences credits a caller for calling sequencing
    noise as hypermutation. A contaminant record would hand it a target in
    the hundreds against a real burden of five.

    Two assertions, because the two halves fail differently. With no
    mutation pass in the pipeline the only honest answer is zero. With one,
    the recoverable count cannot exceed what was counted.
    """
    unmutated = corrupt(experiment_builder.recombined())
    for index, outcome in enumerate(unmutated.run(n=30, seed=3)):
        record = outcome.airr_record(unmutated.refdata)
        assert record["n_mutations"] == 0, (
            f"{label}: fixture error -- no mutation pass is configured")
        assert record["n_observable_mutations"] == 0, (
            f"{label}: record {index} reports "
            f"{record['n_observable_mutations']} recoverable mutations in a "
            f"pipeline with no mutation pass")

    mutated = corrupt(experiment_builder.recombined().mutate(model="s5f", count=5))
    for index, outcome in enumerate(mutated.run(n=30, seed=3)):
        record = outcome.airr_record(mutated.refdata)
        assert record["n_observable_mutations"] <= record["n_mutations"], (
            f"{label}: record {index} reports "
            f"{record['n_observable_mutations']} recoverable against "
            f"{record['n_mutations']} counted")
