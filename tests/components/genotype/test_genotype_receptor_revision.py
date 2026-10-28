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


def _equal_len_two_allele_v_gene(cfg):
    for gene, al in cfg.v_alleles.items():
        if len(al) >= 2 and len(al[0].ungapped_seq) == len(al[1].ungapped_seq):
            return gene, al[0].name, al[1].name
    raise AssertionError("no V gene with two equal-length alleles")


def test_same_haplotype_false_admits_other_chromosome():
    # Heterozygous V gene (a0 on hap0, a1 on hap1, equal length); force the
    # rearrangement to always draw hap0 via chromosome_weights=(1, 0). a1 is
    # carried ONLY on hap1, so it can never be the original V — its appearance as
    # a revised truth V proves same_haplotype=False pulled across chromosomes.
    cfg = _cfg()
    vg, a0, a1 = _equal_len_two_allele_v_gene(cfg)
    g = (Genotype.from_dataconfig(cfg)
         .heterozygous(vg, a0, a1, segment="V")
         .complete_from_reference()
         .chromosome_weights(1.0, 0.0)
         .with_subject("A"))
    res = (ga.Experiment.on(cfg).with_genotype(g)
           .recombine().receptor_revision(prob=1.0, same_haplotype=False)
           .run_records(n=400, seed=4, expose_provenance=True))
    applied = [r for r in res if r.get("receptor_revision_applied")]
    assert applied
    # every rearrangement drew hap0, so any revision target carried only on hap1
    # (here a1) is a genuine cross-chromosome replacement.
    assert all(r["haplotype"] == 0 for r in res)
    assert any(r["truth_v_call"] == a1 for r in applied), (
        "same_haplotype=False must admit the opposite chromosome's allele")
    # sanity: a1 is indeed hap1-exclusive in this genotype
    assert a1 in _carried_on_hap(g, "V", 1) and a1 not in _carried_on_hap(g, "V", 0)


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
