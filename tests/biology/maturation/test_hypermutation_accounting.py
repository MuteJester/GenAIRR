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


@pytest.mark.xfail(strict=True, reason=(
    "KNOWN DEFECT: the replacement base is drawn from all four bases including "
    "the original, so a substitution that changes nothing is still counted. "
    "Measured on human_igh: uniform 27.3% of count=1 records, s5f 3.3-3.7%; at "
    "rate=0.05 the observable/recorded ratio is 0.713 (uniform) and 0.920 (s5f), "
    "so n_mutations is inflated ~40% under the uniform model. "
    "Not fixed here: UniformBase::sample has no access to the current base, so "
    "excluding it means changing the draw at the call site in "
    "substitute_position_constrained, whose fast path explicitly preserves the "
    "legacy RNG consumption shape and byte-identical output. That is a "
    "reproducibility change for every existing user and needs its own branch."))
@verifies("maturation.hypermutation",
          claim="a recorded substitution changes the base, so n_mutations counts "
                "events a mutation caller can recover from the sequence",
          fields=("n_mutations",))
@pytest.mark.parametrize("model", MODELS)
def test_a_recorded_substitution_changes_the_base(model):
    """A substitution to the same base is not a mutation under any molecular
    definition: AID deamination repair yields transitions and transversions,
    never a base replaced by itself.

    Compared case-insensitively on purpose -- case carries the P/N channel, not
    base identity, and comparing case-sensitively hides this entirely.
    """
    exp = experiment_builder.recombined().mutate(model=model, count=1)
    no_ops = []
    for index, outcome in enumerate(exp.run(n=200, seed=0)):
        record = outcome.airr_record(exp.refdata)
        observable = reconstruction.observable_substitutions(outcome.final_simulation())
        if record["n_mutations"] and not observable:
            no_ops.append(index)
    assert no_ops == [], (
        f"{model}: {len(no_ops)}/200 records report a mutation that changed no base "
        f"({100 * len(no_ops) / 200:.1f}%). A replacement drawn from all four bases "
        "including the original inflates every mutation count.")


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
