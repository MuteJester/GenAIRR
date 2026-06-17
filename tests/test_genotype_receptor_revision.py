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
