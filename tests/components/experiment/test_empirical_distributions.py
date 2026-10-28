"""NP-length and trim distributions taken from a DataConfig rather than stated.
"""
from __future__ import annotations
from GenAIRR import _engine as ge
from GenAIRR import CompiledExperiment, Experiment
import GenAIRR as ga  # noqa: E402
from tests.support.fixtures.experiment_dsl import _human_igh_exp

def test_recombine_default_inserts_trim_passes_for_dataconfig():
    # When backed by a DataConfig, recombine() now includes trim
    # passes for V_3, D_5, D_3, J_5 (in that order).
    o = _human_igh_exp().run(n=1, seed=0)[0]
    pn = o.pass_names()
    assert "trim.v_3" in pn
    assert "trim.d_5" in pn
    assert "trim.d_3" in pn
    assert "trim.j_5" in pn


def test_recombine_trim_passes_appear_after_sample_allele_before_assemble():
    o = _human_igh_exp().run(n=1, seed=0)[0]
    pn = o.pass_names()
    # All sample_allele.* must come before all trim.*; all trim.*
    # must come before all assemble.*.
    last_sample = max(i for i, n in enumerate(pn) if n.startswith("sample_allele"))
    first_trim = min(i for i, n in enumerate(pn) if n.startswith("trim."))
    last_trim = max(i for i, n in enumerate(pn) if n.startswith("trim."))
    first_assemble = min(i for i, n in enumerate(pn) if n.startswith("assemble."))
    assert last_sample < first_trim < last_trim < first_assemble


def test_recombine_trim_false_omits_trim_passes():
    o = _human_igh_exp().recombine.__wrapped__ if False else None  # noqa: E501 - placeholder
    # Use a fresh experiment with trim=False explicitly.
    exp = ga.Experiment.on("human_igh").recombine().trim(enabled=False)
    o = exp.run(n=1, seed=0)[0]
    pn = o.pass_names()
    assert not any(n.startswith("trim.") for n in pn)


def test_recombine_with_dataconfig_uses_empirical_np_lengths():
    # Empirical NP1 distribution can produce lengths beyond the
    # uniform [0..6] placeholder. With 50 seeds we should see
    # at least one length > 6 in human_igh (real distribution
    # has a long tail).
    out = _human_igh_exp().run(n=50, seed=0)
    sampled = {len(o.final_simulation().regions()[1]) for o in out}
    assert any(length > 6 for length in sampled), (
        f"expected at least one NP1 > 6 from empirical dist, got {sorted(sampled)}"
    )


def test_recombine_with_raw_refdata_falls_back_to_uniform_np_and_no_trim():
    # Raw RefDataConfig: no DataConfig backing → no trim passes,
    # NP lengths from the uniform [(0..6, 1.0)] placeholder.
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v*01", "v", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j*01", "j", b"TTTAAA", anchor=0)
    o = Experiment.on(cfg).allow_curatable_refdata().recombine().run(n=1, seed=0)[0]
    pn = o.pass_names()
    assert not any(n.startswith("trim.") for n in pn)
    # Pure 5-pass VJ shape, no trims.
    assert pn == [
        "sample_allele.v",
        "sample_allele.j",
        "assemble.v",
        "generate_np.np1",
        "assemble.j",
    ]


def test_recombine_with_raw_refdata_trim_true_is_silent_no_op():
    # trim=True on a raw RefDataConfig is silent — there's no data
    # to source distributions from, so the step compiles to the
    # untrimmed plan rather than raising.
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v*01", "v", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j*01", "j", b"TTTAAA", anchor=0)
    o = Experiment.on(cfg).allow_curatable_refdata().recombine().trim(enabled=True).run(n=1, seed=0)[0]
    assert not any(n.startswith("trim.") for n in o.pass_names())


def test_recombine_explicit_np1_lengths_overrides_empirical_default():
    # User-supplied np1_lengths wins over the empirical default
    # (which would otherwise draw from human_igh's empirical NP1
    # distribution).
    #
    # An earlier version of this test chained `_human_igh_exp().recombine(...)`
    # to verify the second recombine's NP1 distribution overrode the first.
    # That pattern is invalid under the compile-time trim-after-assembly
    # check (the second recombine's trim passes land after the first
    # recombine's assemble passes). The real intent — "explicit np1_lengths
    # makes the recorded NP1 length match the supplied distribution" — is
    # independent of chaining, so we use a single recombine with the
    # explicit override.
    out = (
        Experiment.on("human_igh")
        .recombine(np1_lengths=[(0, 1.0)])
        .run(n=10, seed=0)
    )
    for o in out:
        np1_lens = [
            r.value for r in o.trace().prefix_query("np.np1.length")
        ]
        assert np1_lens, "expected at least one np.np1.length trace entry"
        assert np1_lens[-1] == 0


def test_recombine_d_trims_are_capped_per_segment():
    # Min D allele length in human_igh is ~10-12bp. With margin=3
    # and two_sided cap, the per-side cap is (min_d - 3) // 2 ≈ 3-4.
    # Across 100 seeds no D_5 or D_3 trim should exceed this.
    out = _human_igh_exp().run(n=100, seed=0)
    d5_max = max(o.trace().find("trim.d_5").value for o in out)
    d3_max = max(o.trace().find("trim.d_3").value for o in out)
    # Conservative upper bound: should never approach min_d (~10).
    assert d5_max <= 7, f"D_5 trim {d5_max} too aggressive"
    assert d3_max <= 7, f"D_3 trim {d3_max} too aggressive"


def test_recombine_does_not_panic_across_many_seeds():
    # The whole point of trim caps: no AssembleSegmentPass panic.
    out = _human_igh_exp().run(n=200, seed=0)
    assert len(out) == 200


def test_recombine_with_empirical_data_composes_with_productive_contract():
    out = _human_igh_exp().productive_only().run(n=10, seed=0)
    assert len(out) == 10


def test_extract_recombine_defaults_returns_distributions_for_human_igh():
    # Direct unit test on the helper.
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    cfg = ga.HUMAN_IGH_OGRDB
    defaults = extract_recombine_defaults(cfg)
    assert defaults["np1"] is not None and len(defaults["np1"]) > 0
    assert defaults["np2"] is not None and len(defaults["np2"]) > 0
    assert defaults["trim_v_3"] is not None
    assert defaults["trim_d_5"] is not None
    assert defaults["trim_d_3"] is not None
    assert defaults["trim_j_5"] is not None


def test_extract_recombine_defaults_returns_none_for_vj_d_trims():
    # IGK is a VJ chain — d_alleles is empty, so D trim distributions
    # come back as None.
    from GenAIRR._dataconfig_extract import extract_recombine_defaults

    cfg = ga.HUMAN_IGK_OGRDB
    defaults = extract_recombine_defaults(cfg)
    assert defaults["trim_d_5"] is None
    assert defaults["trim_d_3"] is None


def test_recombine_vj_uses_two_trim_passes():
    # Light chain: V_3 + J_5, no D trims.
    o = ga.Experiment.on("human_igk").recombine().run(n=1, seed=0)[0]
    pn = o.pass_names()
    trims = [n for n in pn if n.startswith("trim.")]
    assert trims == ["trim.v_3", "trim.j_5"]


def test_recombine_step_count_after_trim_does_not_change():
    # `step_count` is the number of FLUENT calls, not the number of
    # passes. recombine() with or without trims is one step.
    a = ga.Experiment.on("human_igh").recombine()
    b = ga.Experiment.on("human_igh").recombine().trim(enabled=False)
    assert a.step_count == 1
    assert b.step_count == 1


def test_empirical_recombine_is_deterministic_under_same_seed():
    a = _human_igh_exp().run(n=5, seed=0xCAFE)
    b = _human_igh_exp().run(n=5, seed=0xCAFE)
    for x, y in zip(a, b):
        assert x.final_simulation().bases() == y.final_simulation().bases()
