"""Zero-P VDJ needs complete future support, including the D coding bytes."""
import pytest

import GenAIRR as ga


def ref(inverted=False, good=True):
    cfg = ga.RefDataConfig.vdj()
    cfg.add_v_allele("V*01", "V", b"ATGTGT", anchor=3)
    cfg.add_d_allele("D-bad*01", "D-bad", b"TTA" if inverted else b"TAA")
    if good:
        cfg.add_d_allele("D-good*01", "D-good", b"AAA")
    cfg.add_j_allele("J*01", "J", b"TGGAAA", anchor=0)
    return cfg


def experiment(cfg, original, inverted):
    exp = ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)], np2_lengths=[(0, 1)])
    if original:
        exp.process_coding_ends(reference_ends_are_coding_ends=True)
    else:
        exp.trim(enabled=False)
    if inverted:
        exp.invert_d(prob=1)
    return exp.productive_only()


@pytest.mark.parametrize("original", [False, True])
@pytest.mark.parametrize("inverted", [False, True])
def test_future_stop_in_d_is_removed_before_allele_sampling(original, inverted):
    cfg = ref(inverted)
    compiled = experiment(cfg, original, inverted).compile()
    for strict in (False, True):
        result = compiled.run_records(n=40, seed=0, strict=strict)
        assert all(row["d_call"] == "D-good*01" and row["productive"] is True for row in result)
        assert all("*" not in row["sequence_aa"] for row in result)
        assert result.validate_records(cfg).ok


@pytest.mark.parametrize("original", [False, True])
def test_no_legal_d_completion_is_a_compile_time_error(original):
    with pytest.raises(ValueError, match="no productive completion"):
        experiment(ref(good=False), original, False).compile()


def test_reported_human_zero_n_profile_is_productive_for_every_emission():
    compiled = (ga.Experiment.on("human_igh")
                .recombine(np1_lengths=[(0, 1)], np2_lengths=[(0, 1)])
                .productive_only().compile())
    result = compiled.run_records(n=40, seed=0, strict=True)
    assert all(row["productive"] is True for row in result)
    assert result.validate_records(compiled.refdata).ok
