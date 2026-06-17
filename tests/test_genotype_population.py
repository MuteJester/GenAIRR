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
