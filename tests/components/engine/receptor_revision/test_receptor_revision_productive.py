"""Real-reference regression seeds for productive replacement selection."""

import pytest

import GenAIRR as ga

pytestmark = pytest.mark.behavior


@pytest.mark.parametrize("seed", [4261, 4265, 4267, 4271, 4279])
@pytest.mark.parametrize("strict", [False, True])
def test_revision_samples_a_productive_replacement_before_committing(seed, strict):
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .receptor_revision(prob=1)
        .productive_only()
    )
    result = exp.run_records(n=1, seed=seed, strict=strict)
    rec = result[0]
    assert rec["productive"] is True
    assert rec["receptor_revision_applied"] is True
    assert rec["v_call"] != rec["original_v_call"]
    assert rec["vj_in_frame"] is True and rec["stop_codon"] is False
    assert result.validate_records(exp.refdata).ok
