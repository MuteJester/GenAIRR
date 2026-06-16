"""Novel / private allele support on genotypes (PR additions)."""
import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def test_add_novel_allele_synthesizes_from_base_and_mutations():
    cfg = _cfg()
    base = cfg.v_alleles["IGHVF1-G1"][0]
    g = Genotype.from_dataconfig(cfg).add_novel_allele(
        "IGHVF1-G1*i01", base="IGHVF1-G1*01", mutations=[(120, "T"), (130, "A")]
    )
    assert g.has_novel()
    assert "IGHVF1-G1*i01" in g.novel_allele_names()
    nv = g._novel["IGHVF1-G1*i01"]["allele"]
    assert nv.ungapped_seq[120] == "T" and nv.ungapped_seq[130] == "A"
    assert len(nv.ungapped_seq) == len(base.ungapped_seq)  # point mutations: no length change
    assert nv.gene == "IGHVF1-G1" and nv.anchor == base.anchor  # metadata inherited


def test_novel_allele_validation():
    cfg = _cfg()
    G = lambda: Genotype.from_dataconfig(cfg)
    with pytest.raises(ValueError, match="base allele"):
        G().add_novel_allele("X*i01", base="NOPE*01", mutations=[(1, "T")])
    with pytest.raises(ValueError, match="collides"):
        G().add_novel_allele("IGHVF1-G1*01", base="IGHVF1-G1*01", mutations=[(1, "T")])
    with pytest.raises(ValueError, match="out of range"):
        G().add_novel_allele("X*i01", base="IGHVF1-G1*01", mutations=[(99999, "T")])
    with pytest.raises(ValueError, match="exactly one"):
        G().add_novel_allele("X*i01", base="IGHVF1-G1*01")  # neither mutations nor sequence
    with pytest.raises(ValueError, match="identical"):
        # mutate to the same base it already is
        b = cfg.v_alleles["IGHVF1-G1"][0].ungapped_seq.upper()
        G().add_novel_allele("X*i01", base="IGHVF1-G1*01", mutations=[(0, b[0])])


def test_novel_allele_placed_and_sampled_appears_in_airr():
    cfg = _cfg()
    gene = "IGHVF1-G1"
    base_seq = cfg.v_alleles[gene][0].ungapped_seq.upper()
    pos = 120
    new = "T" if base_seq[pos] != "T" else "A"
    g = (
        Genotype.from_dataconfig(cfg)
        .add_novel_allele(f"{gene}*i01", base=f"{gene}*01", mutations=[(pos, new)])
        .complete_from_reference()
        .homozygous(gene, f"{gene}*i01")   # carry ONLY the novel allele for this gene
        .with_subject("DONOR_N")
    )
    res = (
        ga.Experiment.on(cfg)
        .with_genotype(g)
        .recombine()
        .run_records(n=300, seed=3, expose_provenance=True)
    )
    # every read assigned to this gene must be the novel allele (truth)
    truth_for_gene = {
        r["truth_v_call"] for r in res if r["truth_v_call"].startswith(gene + "*")
    }
    assert truth_for_gene == {f"{gene}*i01"}, truth_for_gene
    # and the engine actually produced the mutated base in those reads
    novel_reads = [r for r in res if r["truth_v_call"] == f"{gene}*i01"]
    assert novel_reads, "expected some reads from the novel allele"
    assert any(r["sequence"][pos].upper() == new for r in novel_reads)


def test_genotype_without_novel_is_unaffected():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("S1")
    assert g.has_novel() is False
    # still runs (uses base refdata)
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(n=10, seed=1)
    assert len(res) == 10
