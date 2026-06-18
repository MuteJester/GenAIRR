"""DSL-level guards for Experiment.with_genotype (PR1)."""
import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def _full_genotype():
    return Genotype.from_dataconfig(_cfg()).complete_from_reference().with_subject("S1")


def test_with_genotype_then_restrict_alleles_raises():
    g = _full_genotype()
    v_gene = next(iter(_cfg().v_alleles))
    a1 = _cfg().v_alleles[v_gene][0].name
    with pytest.raises(ValueError, match="mutually exclusive"):
        ga.Experiment.on(_cfg()).with_genotype(g).restrict_alleles(v=a1)


def test_restrict_alleles_then_with_genotype_raises():
    g = _full_genotype()
    v_gene = next(iter(_cfg().v_alleles))
    a1 = _cfg().v_alleles[v_gene][0].name
    with pytest.raises(ValueError, match="mutually exclusive"):
        ga.Experiment.on(_cfg()).restrict_alleles(v=a1).with_genotype(g)


def test_recombine_weights_then_with_genotype_raises():
    g = _full_genotype()
    v_gene = next(iter(_cfg().v_alleles))
    a1 = _cfg().v_alleles[v_gene][0].name
    with pytest.raises(ValueError, match="mutually exclusive"):
        ga.Experiment.on(_cfg()).recombine(v_allele_weights={a1: 2.0}).with_genotype(g)


def test_with_genotype_then_recombine_weights_raises():
    g = _full_genotype()
    v_gene = next(iter(_cfg().v_alleles))
    a1 = _cfg().v_alleles[v_gene][0].name
    with pytest.raises(ValueError, match="mutually exclusive"):
        ga.Experiment.on(_cfg()).with_genotype(g).recombine(v_allele_weights={a1: 2.0})


def test_receptor_revision_with_genotype_compiles_and_runs():
    # Genotype-aware receptor revision is now supported: the replacement V is
    # restricted to carried alleles on the drawn chromosome (no longer rejected).
    g = _full_genotype()
    exp = ga.Experiment.on(_cfg()).with_genotype(g).recombine().receptor_revision(prob=0.5)
    exp.compile()  # no raise
    res = exp.run_records(n=10, seed=1, expose_provenance=True)
    assert len(res) == 10


def test_genotype_with_clonal_fork_raises_at_compile():
    g = _full_genotype()
    exp = (
        ga.Experiment.on(_cfg())
        .with_genotype(g)
        .recombine()
        .clonal_lineage(n_clones=2)
    )
    with pytest.raises(ValueError, match="not supported together with"):
        exp.compile()


def test_cartridge_hash_mismatch_raises():
    # Genotype built on IGH, attached to a TCRB experiment → mismatch.
    g = Genotype.from_dataconfig(_cfg()).complete_from_reference()
    other = gdata.HUMAN_TCRB_IMGT
    with pytest.raises(ValueError, match="different cartridge|content hash"):
        ga.Experiment.on(other).with_genotype(g)
