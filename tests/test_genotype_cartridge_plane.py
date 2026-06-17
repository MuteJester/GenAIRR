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


def test_d_on_vj_rejected_via_chaintype_object():
    class _FakeChainType:  # ChainType-like: exposes has_d
        has_d = False
    m = PopulationGenotypeModel(model_id="t", source="s",
        haplotype_deletion_prob={"D": {"IGHD1-1": 0.5}})
    with pytest.raises(ValueError, match="VJ|D-segment|D segment"):
        m.validate(chain_type=_FakeChainType())


def test_gene_keys_must_be_strings():
    with pytest.raises(ValueError, match="gene"):
        PopulationGenotypeModel(model_id="t", source="s",
            allele_frequencies={"V": {123: {"IGHV1-2*02": 1.0}}}).validate()
    with pytest.raises(ValueError, match="gene"):
        PopulationGenotypeModel(model_id="t", source="s",
            haplotype_deletion_prob={"V": {123: 0.1}}).validate()


def test_allow_nonfunctional_must_be_bool():
    with pytest.raises(ValueError, match="allow_nonfunctional"):
        PopulationGenotypeModel(model_id="t", source="s", novel_alleles=[
            PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHV1-2*02", sequence="ACGT", frequency=1.0,
                allow_nonfunctional="yes")]).validate()


def test_content_checksum_is_canonical():
    vg = "IGHV1-2"
    m1 = PopulationGenotypeModel(model_id="m", source="s",
        allele_frequencies={"V": {vg: {"IGHV1-2*02": 2, "IGHV1-2*04": 1}}},
        novel_alleles=[PopulationNovelAllele(name="IGHV1-2*99", segment="V",
            base_allele="IGHV1-2*02", sequence="acgt", frequency=1.0)])
    # same content, different dict insertion order + int-vs-float + DNA case
    m2 = PopulationGenotypeModel(model_id="m", source="s",
        allele_frequencies={"V": {vg: {"IGHV1-2*04": 1.0, "IGHV1-2*02": 2.0}}},
        novel_alleles=[PopulationNovelAllele(name="IGHV1-2*99", segment="V",
            base_allele="IGHV1-2*02", sequence="ACGT", frequency=1.0)])
    assert m1.content_checksum() == m2.content_checksum()

    m3 = PopulationGenotypeModel(model_id="m", source="s",
        allele_frequencies={"V": {vg: {"IGHV1-2*02": 3.0, "IGHV1-2*04": 1.0}}})
    assert m3.content_checksum() != m1.content_checksum()


def test_genotype_priors_field_default_and_checksum_invariant():
    cfg = _cfg()
    # default-new and bundled both have no plane
    assert getattr(cfg, "genotype_priors", "MISSING") is None
    base_checksum = cfg.compute_checksum()

    import copy
    cfg2 = copy.deepcopy(cfg)
    cfg2.genotype_priors = None  # explicitly None must not change the checksum
    assert cfg2.compute_checksum() == base_checksum

    cfg3 = copy.deepcopy(cfg)
    cfg3.genotype_priors = PopulationGenotypeModel(model_id="m", source="s")
    assert cfg3.compute_checksum() != base_checksum  # a real plane is cartridge identity


def test_genotype_priors_pickle_round_trip():
    import copy
    cfg = copy.deepcopy(_cfg())
    m = PopulationGenotypeModel(model_id="m", source="s",
        haplotype_deletion_prob={"V": {next(iter(cfg.v_alleles)): 0.2}})
    cfg.genotype_priors = m
    back = pickle.loads(pickle.dumps(cfg, protocol=4))
    assert back.genotype_priors.model_id == "m"
    assert back.genotype_priors.content_checksum() == m.content_checksum()


def test_manifest_genotype_priors_block():
    import copy
    cfg = copy.deepcopy(_cfg())
    block = cfg.cartridge_manifest()["models"]["genotype_priors"]
    assert block["available"] is False

    vg = next(iter(cfg.v_alleles))
    cfg.genotype_priors = PopulationGenotypeModel(
        model_id="m1", source="VDJbase-toy", version="1",
        allele_frequencies={"V": {vg: {cfg.v_alleles[vg][0].name: 1.0}}},
        haplotype_deletion_prob={"V": {vg: 0.1}},
        novel_alleles=[PopulationNovelAllele(name="IGHV1-2*99", segment="V",
            base_allele="IGHV1-2*02", sequence="ACGT", frequency=1.0)],
    )
    block = cfg.cartridge_manifest()["models"]["genotype_priors"]
    assert block["available"] is True
    assert block["model_id"] == "m1"
    assert block["source"] == "VDJbase-toy"
    assert block["model_checksum"] == cfg.genotype_priors.content_checksum()
    assert block["freq_gene_counts"]["V"] == 1
    assert block["deletion_gene_counts"]["V"] == 1
    assert block["novel_allele_count"] == 1
    assert block["chromosome_weights"] == [0.5, 0.5]
    assert block["source_field"] == "DataConfig.genotype_priors"
