"""Novel / private allele support on genotypes."""
import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


# Guaranteed-safe SNPs for IGHVF1-G1*01: wobble of non-T-starting codons in
# the framework — cannot create a stop, anchor untouched.
_SAFE = [(38, "C"), (41, "A")]


def test_add_novel_allele_synthesizes_from_base_and_mutations():
    cfg = _cfg()
    base = cfg.v_alleles["IGHVF1-G1"][0]
    g = Genotype.from_dataconfig(cfg).add_novel_allele(
        "IGHVF1-G1*i01", base="IGHVF1-G1*01", mutations=_SAFE
    )
    assert g.has_novel()
    assert "IGHVF1-G1*i01" in g.novel_allele_names()
    nv = g._novel["IGHVF1-G1*i01"]["allele"]
    assert nv.ungapped_seq[38] == "C" and nv.ungapped_seq[41] == "A"
    assert len(nv.ungapped_seq) == len(base.ungapped_seq)  # substitution-only
    assert nv.gene == "IGHVF1-G1" and nv.anchor == base.anchor  # gene/anchor inherited
    # gapped sequence projected consistently (ungapped derived from it)
    assert nv.gapped_seq.replace(".", "") == nv.ungapped_seq
    assert g._novel["IGHVF1-G1*i01"]["functional"] is True


def test_novel_name_gene_must_match_base_gene():
    cfg = _cfg()
    with pytest.raises(ValueError, match="implies gene"):
        Genotype.from_dataconfig(cfg).add_novel_allele(
            "IGHVF1-G2*i01", base="IGHVF1-G1*01", mutations=_SAFE  # name gene != base gene
        )


def test_novel_allele_basic_validation():
    cfg = _cfg()
    G = lambda: Genotype.from_dataconfig(cfg)
    with pytest.raises(ValueError, match="base allele"):
        G().add_novel_allele("NOPE*i01", base="NOPE*01", mutations=[(1, "T")])
    with pytest.raises(ValueError, match="collides with a catalogue"):
        G().add_novel_allele("IGHVF1-G1*01", base="IGHVF1-G1*01", mutations=_SAFE)
    with pytest.raises(ValueError, match="out of range"):
        G().add_novel_allele("IGHVF1-G1*i01", base="IGHVF1-G1*01", mutations=[(99999, "T")])
    with pytest.raises(ValueError, match="exactly one"):
        G().add_novel_allele("IGHVF1-G1*i01", base="IGHVF1-G1*01")
    with pytest.raises(ValueError, match="must be an int"):
        G().add_novel_allele("IGHVF1-G1*i01", base="IGHVF1-G1*01", mutations=[(1.5, "A")])


def test_cross_segment_name_collision_rejected():
    # A novel can't be named for a different gene/segment than its base, so
    # naming a V-derived novel after a real J allele is rejected outright —
    # the truth table can never mislabel the real J row as novel.
    cfg = _cfg()
    j_name = next(iter(cfg.j_alleles[next(iter(cfg.j_alleles))])).name
    with pytest.raises(ValueError):
        Genotype.from_dataconfig(cfg).add_novel_allele(
            j_name, base="IGHVF1-G1*01", mutations=_SAFE
        )


def test_nonfunctional_novel_rejected_by_default_and_allowed_explicitly():
    cfg = _cfg()
    base = cfg.v_alleles["IGHVF1-G1"][0]
    # Force a stop codon at framework codon 13 (positions 39,40,41 -> TAA).
    stop = [(39, "T"), (40, "A"), (41, "A")]
    with pytest.raises(ValueError, match="non-functional.*stop codon"):
        Genotype.from_dataconfig(cfg).add_novel_allele(
            "IGHVF1-G1*i01", base="IGHVF1-G1*01", mutations=stop
        )
    # Explicit override keeps it, marked non-functional.
    g = Genotype.from_dataconfig(cfg).add_novel_allele(
        "IGHVF1-G1*i01", base="IGHVF1-G1*01", mutations=stop, allow_nonfunctional=True
    )
    assert g._novel["IGHVF1-G1*i01"]["functional"] is False


def test_broken_anchor_codon_rejected():
    cfg = _cfg()
    base = cfg.v_alleles["IGHVF1-G1"][0]
    a = base.anchor
    # rewrite the conserved Cys anchor codon to GGG (Gly)
    with pytest.raises(ValueError, match="conserved anchor codon"):
        Genotype.from_dataconfig(cfg).add_novel_allele(
            "IGHVF1-G1*i01", base="IGHVF1-G1*01",
            mutations=[(a, "G"), (a + 1, "G"), (a + 2, "G")],
        )


def test_novel_allele_placed_and_sampled_appears_in_airr():
    cfg = _cfg()
    gene = "IGHVF1-G1"
    g = (
        Genotype.from_dataconfig(cfg)
        .add_novel_allele(f"{gene}*i01", base=f"{gene}*01", mutations=_SAFE)
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
    truth_for_gene = {
        r["truth_v_call"] for r in res if r["truth_v_call"].startswith(gene + "*")
    }
    assert truth_for_gene == {f"{gene}*i01"}, truth_for_gene
    novel_reads = [r for r in res if r["truth_v_call"] == f"{gene}*i01"]
    assert novel_reads
    assert any(r["sequence"][38].upper() == "C" for r in novel_reads)


def test_to_table_and_tsv_expose_novel(tmp_path):
    cfg = _cfg()
    gene = "IGHVF1-G1"
    g = (
        Genotype.from_dataconfig(cfg)
        .add_novel_allele(f"{gene}*i01", base=f"{gene}*01", mutations=_SAFE)
        .heterozygous(gene, f"{gene}*01", f"{gene}*i01")
    )
    row = next(r for r in g.to_table() if r["gene"] == gene)
    assert row["novel"] == [f"{gene}*i01"]
    p = tmp_path / "truth.tsv"
    g.to_tsv(str(p))
    header = p.read_text().splitlines()[0].split("\t")
    assert "novel" in header


def test_genotype_without_novel_is_unaffected():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("S1")
    assert g.has_novel() is False
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(n=10, seed=1)
    assert len(res) == 10
