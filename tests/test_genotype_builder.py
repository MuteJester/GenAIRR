"""Builder-level tests for GenAIRR.genotype.Genotype (PR1)."""
import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def test_homozygous_then_subject_builds_diploid_genotype():
    cfg = _cfg()
    v_gene = next(iter(cfg.v_alleles))
    a1 = cfg.v_alleles[v_gene][0].name
    g = Genotype.from_dataconfig(cfg).homozygous(v_gene, a1).with_subject("S1")
    assert g.subject_id == "S1"
    assert g.carried_alleles("V", v_gene) == {a1}


def test_unknown_allele_name_raises():
    cfg = _cfg()
    v_gene = next(iter(cfg.v_alleles))
    with pytest.raises(ValueError, match="not a known"):
        Genotype.from_dataconfig(cfg).homozygous(v_gene, "IGHV-NOPE*99")


def test_unknown_gene_raises():
    cfg = _cfg()
    with pytest.raises(ValueError, match="not a known"):
        Genotype.from_dataconfig(cfg).homozygous("NOSUCHGENE", "x*01")


def test_strict_genotype_reports_unspecified_gene():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg)
    assert g.is_specified("V", next(iter(cfg.v_alleles))) is False


def test_complete_from_reference_specifies_every_gene():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference()
    assert all(g.is_specified("V", gene) for gene in cfg.v_alleles)
    assert all(g.is_specified("J", gene) for gene in cfg.j_alleles)


def test_heterozygous_to_table_reports_zygosity():
    cfg = _cfg()
    v_gene = next(iter(cfg.v_alleles))
    names = [a.name for a in cfg.v_alleles[v_gene]]
    if len(names) < 2:
        pytest.skip("need >=2 alleles for heterozygous test")
    g = Genotype.from_dataconfig(cfg).heterozygous(v_gene, names[0], names[1])
    rows = [r for r in g.to_table() if r["gene"] == v_gene]
    assert rows and rows[0]["zygosity"] == "heterozygous"
    assert g.carried_alleles("V", v_gene) == {names[0], names[1]}


def test_delete_gene_one_haplotype_keeps_the_other():
    cfg = _cfg()
    v_gene = next(iter(cfg.v_alleles))
    a1 = cfg.v_alleles[v_gene][0].name
    g = Genotype.from_dataconfig(cfg).homozygous(v_gene, a1).delete_gene(v_gene, haplotype=1)
    # haplotype 0 still carries a1; haplotype 1 deleted
    rows = [r for r in g.to_table() if r["gene"] == v_gene]
    assert rows[0]["haplotype_0"] == [a1]
    assert rows[0]["haplotype_1"] == []


def test_permissive_is_flagged():
    cfg = _cfg()
    g = Genotype.permissive(cfg)
    assert g.is_permissive is True
    assert Genotype.from_dataconfig(cfg).is_permissive is False
