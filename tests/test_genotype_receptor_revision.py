"""Genotype-aware (same-haplotype) receptor revision."""
import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def _two_allele_v_gene(cfg):
    return next(g for g, al in cfg.v_alleles.items() if len(al) >= 2)


def test_receptor_revision_same_haplotype_must_be_bool():
    cfg = _cfg()
    with pytest.raises(ValueError, match="same_haplotype"):
        ga.Experiment.on(cfg).recombine().receptor_revision(prob=0.5, same_haplotype="yes")


def test_genotype_plus_receptor_revision_compiles_and_runs():
    cfg = _cfg()
    g = Genotype.sample(cfg, seed=0, subject_id="A")
    res = (ga.Experiment.on(cfg).with_genotype(g)
           .recombine().receptor_revision(prob=1.0)
           .run_records(n=10, seed=1, expose_provenance=True))
    assert len(res) == 10


def test_genotype_plus_clonal_fork_still_rejected():
    cfg = _cfg()
    g = Genotype.sample(cfg, seed=0, subject_id="A")
    with pytest.raises(ValueError, match="clonal|expand_clones|clonal_lineage"):
        (ga.Experiment.on(cfg).with_genotype(g)
         .recombine().clonal_lineage(n_clones=2).compile())


def test_run_cohort_plus_receptor_revision_runs():
    cfg = _cfg()
    gs = [Genotype.sample(cfg, seed=s, subject_id=f"D{s}") for s in range(2)]
    c = (ga.Experiment.on(cfg).recombine().receptor_revision(prob=1.0)
         .run_cohort(gs, n_per_subject=5, seed=0))
    assert len(c) == 10
