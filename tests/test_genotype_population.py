"""Population-prior sampling: Genotype.sample(...)."""
import math

import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def test_sample_is_fully_specified_and_deterministic():
    cfg = _cfg()
    g1 = Genotype.sample(cfg, seed=5)
    g2 = Genotype.sample(cfg, seed=5)
    assert g1.to_table() == g2.to_table()  # deterministic
    for gene in cfg.v_alleles:
        assert g1.is_specified("V", gene)
    for gene in cfg.j_alleles:
        assert g1.is_specified("J", gene)
    for gene in cfg.d_alleles:
        assert g1.is_specified("D", gene)
    res = ga.Experiment.on(cfg).with_genotype(g1).recombine().run_records(n=20, seed=1)
    assert len(res) == 20


def test_different_seeds_differ():
    cfg = _cfg()
    assert Genotype.sample(cfg, seed=1).to_table() != Genotype.sample(cfg, seed=2).to_table()


def test_segments_must_cover_required():
    cfg = _cfg()  # VDJ
    with pytest.raises(ValueError, match="required segment"):
        Genotype.sample(cfg, seed=0, segments_to_sample=("V",))  # omits D, J
    with pytest.raises(ValueError, match="no .* segment|required segment|unknown segment"):
        Genotype.sample(cfg, seed=0, segments_to_sample=("V", "D", "J", "C"))


def test_frequencies_bias_homozygosity():
    cfg = _cfg()
    vg = next(g for g, al in cfg.v_alleles.items() if len(al) >= 2)
    a0, a1 = (a.name for a in cfg.v_alleles[vg][:2])
    freqs = {"V": {vg: {a0: 100.0, a1: 1.0}}}  # heavily favour a0
    homo_a0 = sum(
        1
        for s in range(60)
        if Genotype.sample(cfg, seed=s, allele_frequencies=freqs).carried_alleles("V", vg) == {a0}
    )
    assert homo_a0 > 40  # dominant allele -> usually homozygous-common


def test_zero_weight_excludes_allele():
    cfg = _cfg()
    vg = next(g for g, al in cfg.v_alleles.items() if len(al) >= 2)
    a0, a1 = (a.name for a in cfg.v_alleles[vg][:2])
    freqs = {"V": {vg: {a0: 1.0, a1: 0.0}}}  # a1 excluded
    for s in range(40):
        assert a1 not in Genotype.sample(cfg, seed=s, allele_frequencies=freqs).carried_alleles("V", vg)


def test_frequency_validation():
    cfg = _cfg()
    vg = next(iter(cfg.v_alleles))
    a0 = cfg.v_alleles[vg][0].name
    with pytest.raises(ValueError, match="not a known"):
        Genotype.sample(cfg, seed=0, allele_frequencies={"V": {vg: {"NOPE*9": 1.0}}})
    with pytest.raises(ValueError, match="must be > 0"):
        Genotype.sample(cfg, seed=0, allele_frequencies={"V": {vg: {a0: 0.0}}})
    with pytest.raises(ValueError, match="finite and >= 0"):
        Genotype.sample(cfg, seed=0, allele_frequencies={"V": {vg: {a0: -1.0}}})


def test_flat_gene_shape_accepted():
    cfg = _cfg()
    vg = next(iter(cfg.v_alleles))
    a0 = cfg.v_alleles[vg][0].name
    g = Genotype.sample(cfg, seed=0, allele_frequencies={vg: {a0: 1.0}})  # flat, unambiguous
    assert g.carried_alleles("V", vg) <= {a0}


def test_usage_as_prior_requires_typed_usage():
    cfg = _cfg()  # no typed allele_usage
    with pytest.raises(ValueError, match="allele_usage|usage_as_prior"):
        Genotype.sample(cfg, seed=0, allele_frequencies="usage_as_prior")


def test_deletion_prob_per_gene():
    cfg = _cfg()
    drop = list(cfg.v_alleles)[0]
    other = list(cfg.v_alleles)[1]
    g = Genotype.sample(cfg, seed=0, haplotype_deletion_prob={"V": {drop: 1.0}})
    assert g.carried_alleles("V", drop) == set()  # fully deleted
    assert g.carried_alleles("V", other)          # others still carried


def test_deletion_validation():
    cfg = _cfg()
    with pytest.raises(ValueError, match=r"in \[0, 1\]"):
        Genotype.sample(cfg, seed=0, haplotype_deletion_prob=1.5)
    with pytest.raises(ValueError, match="unknown"):
        Genotype.sample(cfg, seed=0, haplotype_deletion_prob={"V": {"NOPE": 0.5}})


def test_ensure_viable_raises_when_all_deleted():
    cfg = _cfg()
    with pytest.raises(ValueError, match="viable|deletion"):
        Genotype.sample(cfg, seed=0, haplotype_deletion_prob=1.0, max_resamples=10)
    g = Genotype.sample(cfg, seed=0, haplotype_deletion_prob=1.0, ensure_viable=False)
    assert all(g.carried_alleles("V", gene) == set() for gene in cfg.v_alleles)


def test_max_resamples_validation():
    cfg = _cfg()
    with pytest.raises(ValueError, match="max_resamples"):
        Genotype.sample(cfg, seed=0, max_resamples=0)


def test_nested_freq_validation_no_silent_ignore():
    cfg = _cfg()
    vg = next(iter(cfg.v_alleles))
    with pytest.raises(ValueError, match="not in the cartridge"):
        Genotype.sample(cfg, seed=0, allele_frequencies={"V": {"NOPE": {"NOPE*01": 1.0}}})
    with pytest.raises(ValueError, match="must be a mapping"):
        Genotype.sample(cfg, seed=0, allele_frequencies={"V": []})
    with pytest.raises(ValueError, match="non-empty mapping"):
        Genotype.sample(cfg, seed=0, allele_frequencies={"V": {vg: None}})
    with pytest.raises(ValueError, match="non-empty mapping"):
        Genotype.sample(cfg, seed=0, allele_frequencies={"V": {vg: {}}})
    with pytest.raises(ValueError, match="must be a mapping"):
        Genotype.sample(cfg, seed=0, allele_frequencies=True)


def test_deletion_dict_non_mapping_raises():
    cfg = _cfg()
    with pytest.raises(ValueError, match="must be a mapping"):
        Genotype.sample(cfg, seed=0, haplotype_deletion_prob={"V": []})


def test_segments_dedup_and_order_invariant():
    cfg = _cfg()
    with pytest.raises(ValueError, match="duplicate segment"):
        Genotype.sample(cfg, seed=0, segments_to_sample=("V", "V", "D", "J"))
    # order-invariant: same seed, reordered segments -> identical genotype
    a = Genotype.sample(cfg, seed=3, segments_to_sample=("V", "D", "J"))
    b = Genotype.sample(cfg, seed=3, segments_to_sample=("J", "D", "V"))
    assert a.to_table() == b.to_table()


def test_zero_weight_chromosome_never_expressed():
    cfg = _cfg()
    # 100% express chromosome 1; chromosome 0 must never appear in records.
    g = Genotype.sample(cfg, seed=4, chromosome_weights=(0.0, 1.0))
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(
        n=60, seed=1, expose_provenance=True
    )
    assert all(r["haplotype"] == 1 for r in res)


def test_zero_weight_chromosome_not_counted_viable():
    cfg = _cfg()
    # Delete heavily on chromosome 1 only, but give all weight to chromosome 1:
    # the only complete haplotype (chrom 0) has zero weight -> not viable -> raises.
    with pytest.raises(ValueError, match="positive-weight|viable"):
        Genotype.sample(
            cfg,
            seed=0,
            haplotype_deletion_prob={"J": {gene: 1.0 for gene in cfg.j_alleles}},
            chromosome_weights=(1.0, 0.0),
            max_resamples=5,
        )


def test_chromosome_weights_validation_in_sample():
    cfg = _cfg()
    with pytest.raises(ValueError, match="finite"):
        Genotype.sample(cfg, seed=0, chromosome_weights=(float("nan"), 1.0))
    with pytest.raises(ValueError, match="non-negative"):
        Genotype.sample(cfg, seed=0, chromosome_weights=(-1.0, 1.0))


def test_conditioned_vs_unconditioned_viability():
    cfg = _cfg()
    # ensure_viable=False can yield an infeasible (no complete haplotype) draw at
    # high deletion; ensure_viable=True never does.
    g_uncond = Genotype.sample(cfg, seed=0, haplotype_deletion_prob=1.0, ensure_viable=False)
    assert not g_uncond._is_viable(cfg)  # all deleted -> not viable
    with pytest.raises(ValueError):
        Genotype.sample(cfg, seed=0, haplotype_deletion_prob=1.0, max_resamples=5)


def test_end_to_end_truth_calls_are_carried():
    cfg = _cfg()
    g = Genotype.sample(cfg, seed=11, haplotype_deletion_prob=0.1)
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(
        n=200, seed=2, expose_provenance=True
    )
    for r in res:
        for seg, col in (("V", "truth_v_call"), ("D", "truth_d_call"), ("J", "truth_j_call")):
            call = r[col]
            if not call:
                continue
            gene = call.split("*")[0]
            assert call in g.carried_alleles(seg, gene), (seg, call)
