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
