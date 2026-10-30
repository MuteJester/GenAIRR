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
