"""Cartridge genotype plane: PopulationGenotypeModel + Genotype.sample consumption."""
import pickle

import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype
from GenAIRR.genotype_priors import PopulationGenotypeModel, PopulationNovelAllele


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def _vg_two_alleles(cfg):
    """A V gene with >= 2 alleles, plus its first two allele names."""
    vg = next(g for g, al in cfg.v_alleles.items() if len(al) >= 2)
    a0, a1 = (a.name for a in cfg.v_alleles[vg][:2])
    return vg, a0, a1


def test_model_requires_identity():
    m = PopulationGenotypeModel(model_id="", source="x")
    with pytest.raises(ValueError, match="model_id"):
        m.validate()
    m = PopulationGenotypeModel(model_id="x", source="")
    with pytest.raises(ValueError, match="source"):
        m.validate()


def test_model_validates_shapes():
    cfg = _cfg()
    vg, a0, a1 = _vg_two_alleles(cfg)
    ok = PopulationGenotypeModel(
        model_id="toy", source="unit-test",
        allele_frequencies={"V": {vg: {a0: 2.0, a1: 1.0}}},
        haplotype_deletion_prob={"V": {vg: 0.1}},
        chromosome_weights=(0.5, 0.5),
    )
    ok.validate(chain_type="vdj")  # no raise

    with pytest.raises(ValueError, match="finite"):
        PopulationGenotypeModel(model_id="t", source="s",
            allele_frequencies={"V": {vg: {a0: float("nan")}}}).validate()
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        PopulationGenotypeModel(model_id="t", source="s",
            haplotype_deletion_prob={"V": {vg: 1.5}}).validate()
    with pytest.raises(ValueError, match="non-negative"):
        PopulationGenotypeModel(model_id="t", source="s",
            chromosome_weights=(-1.0, 1.0)).validate()
    with pytest.raises(ValueError, match="at least one"):
        PopulationGenotypeModel(model_id="t", source="s",
            allele_frequencies={"V": {vg: {a0: 0.0, a1: 0.0}}}).validate()


def test_novel_shape_validation():
    base = PopulationNovelAllele(name="IGHV1-2*99", segment="V",
        base_allele="IGHV1-2*02", sequence="ACGT", frequency=1.0)
    PopulationGenotypeModel(model_id="t", source="s", novel_alleles=[base]).validate()

    with pytest.raises(ValueError, match="must be > 0"):
        PopulationGenotypeModel(model_id="t", source="s", novel_alleles=[
            PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHV1-2*02", sequence="ACGT", frequency=0.0)]).validate()
    with pytest.raises(ValueError, match="A/C/G/T|DNA"):
        PopulationGenotypeModel(model_id="t", source="s", novel_alleles=[
            PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHV1-2*02", sequence="ACGX", frequency=1.0)]).validate()
    with pytest.raises(ValueError, match="duplicate|unique"):
        PopulationGenotypeModel(model_id="t", source="s", novel_alleles=[
            PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHV1-2*02", sequence="ACGT", frequency=1.0),
            PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHV1-2*02", sequence="ACGA", frequency=1.0)]).validate()


def test_d_on_vj_rejected():
    m = PopulationGenotypeModel(model_id="t", source="s",
        allele_frequencies={"D": {"IGHD1-1": {"IGHD1-1*01": 1.0}}})
    with pytest.raises(ValueError, match="VJ|D-segment|D segment"):
        m.validate(chain_type="vj")
