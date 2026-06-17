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


def _carried_on_hap(g, seg, hap):
    """All allele names carried on chromosome `hap` for a segment (union over
    genes), read from the genotype's per-haplotype slots."""
    out = set()
    for _gene, haps in g._slots[seg].items():
        out.update(a for (a, _c, _w) in haps[hap])
    return out


def test_same_haplotype_revision_truth_is_carried_changed_and_on_drawn_chromosome():
    cfg = _cfg()
    g = Genotype.sample(cfg, seed=11, subject_id="A")
    res = (ga.Experiment.on(cfg).with_genotype(g)
           .recombine().receptor_revision(prob=1.0)
           .run_records(n=200, seed=3, expose_provenance=True))
    applied = [r for r in res if r.get("receptor_revision_applied")]
    assert applied, "prob=1 with eligible alternates must revise some records"
    for r in applied:
        tv, ov = r["truth_v_call"], r["original_v_call"]
        assert tv and ov and tv != ov                 # revision changed the V
        # revised V is carried on the SAME chromosome the rearrangement drew from
        assert tv in _carried_on_hap(g, "V", r["haplotype"])


def test_same_haplotype_false_admits_other_chromosome():
    cfg = _cfg()
    g = Genotype.sample(cfg, seed=12, subject_id="A")
    res = (ga.Experiment.on(cfg).with_genotype(g)
           .recombine().receptor_revision(prob=1.0, same_haplotype=False)
           .run_records(n=200, seed=4, expose_provenance=True))
    applied = [r for r in res if r.get("receptor_revision_applied")]
    assert applied
    for r in applied:
        # under both-haplotype mode the revised V is carried somewhere in the
        # diploid genotype (union), and still differs from the original
        gene = r["truth_v_call"].split("*")[0]
        assert r["truth_v_call"] in g.carried_alleles("V", gene)
        assert r["truth_v_call"] != r["original_v_call"]


def test_novel_v_can_be_revision_target():
    cfg = _cfg()
    vg = _two_allele_v_gene(cfg)
    base = cfg.v_alleles[vg][0].name
    bs = next(a for a in cfg.v_alleles[vg] if a.name == base).ungapped_seq.upper()
    novel_seq = None
    for pos in range(len(bs)):
        for nt in "ACGT":
            if nt == bs[pos]:
                continue
            cand = bs[:pos] + nt + bs[pos + 1:]
            try:
                Genotype.from_dataconfig(cfg).add_novel_allele(
                    f"{vg}*97", base=base, sequence=cand, segment="V")
                novel_seq = cand
                break
            except ValueError:
                continue
        if novel_seq:
            break
    novel = f"{vg}*97"
    g = (Genotype.from_dataconfig(cfg)
         .add_novel_allele(novel, base=base, sequence=novel_seq, segment="V")
         .heterozygous(vg, base, novel).complete_from_reference().with_subject("A"))
    res = (ga.Experiment.on(cfg).with_genotype(g)
           .recombine().receptor_revision(prob=1.0)
           .run_records(n=200, seed=5, expose_provenance=True))
    # the carried novel can appear as a revised truth V
    assert any(r.get("truth_v_call") == novel for r in res)
