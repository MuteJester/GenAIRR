"""End-to-end phased-genotype recombination tests (PR1)."""
import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def test_phased_recombine_only_emits_carried_allele_for_overridden_gene():
    cfg = _cfg()
    # Pick a V gene with >=2 alleles and carry ONLY its second allele.
    v_gene = next(g for g, al in cfg.v_alleles.items() if len(al) >= 2)
    names = [a.name for a in cfg.v_alleles[v_gene]]
    carried = names[1]
    g = (
        Genotype.from_dataconfig(cfg)
        .complete_from_reference("homozygous_first_reference")
        .homozygous(v_gene, carried)
        .with_subject("S1")
    )
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(
        n=300, seed=1, expose_provenance=True
    )
    # Use the ground-truth call (truth_v_call) — the evidence-based v_call
    # can be ambiguous (comma-joined) between similar alleles.
    seen_for_gene = {
        r["truth_v_call"] for r in res if r["truth_v_call"].startswith(v_gene + "*")
    }
    # Only the carried allele of that gene may appear (never names[0]).
    assert seen_for_gene <= {carried}, seen_for_gene


def test_records_carry_subject_and_haplotype_and_result_exposes_genotype():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("S1")
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(n=20, seed=3)
    assert all(r["subject_id"] == "S1" for r in res)
    assert all(r["haplotype"] in (0, 1) for r in res)
    assert res.genotypes is not None
    assert res.genotypes[0].subject_id == "S1"


def test_no_genotype_result_has_no_genotypes_and_no_haplotype_field():
    cfg = _cfg()
    res = ga.Experiment.on(cfg).recombine().run_records(n=10, seed=3)
    assert res.genotypes is None
    assert "subject_id" not in res[0]
    assert "haplotype" not in res[0]


def test_phased_run_is_deterministic_under_same_seed():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("S1")
    exp = ga.Experiment.on(cfg).with_genotype(g).recombine()
    a = exp.run_records(n=40, seed=77)
    b = exp.run_records(n=40, seed=77)
    assert [r["v_call"] for r in a] == [r["v_call"] for r in b]
    assert [r["j_call"] for r in a] == [r["j_call"] for r in b]


def test_deleted_gene_is_never_sampled():
    cfg = _cfg()
    drop = list(cfg.v_alleles)[1]
    g = (
        Genotype.from_dataconfig(cfg)
        .complete_from_reference()
        .delete_gene(drop, haplotype="both", segment="V")
        .with_subject("S1")
    )
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(
        n=300, seed=5, expose_provenance=True
    )
    assert all(not r["truth_v_call"].startswith(drop + "*") for r in res)


def test_heterozygous_expression_is_roughly_balanced():
    cfg = _cfg()
    v_gene = next(g for g, al in cfg.v_alleles.items() if len(al) >= 2)
    names = [a.name for a in cfg.v_alleles[v_gene]]
    a0, a1 = names[0], names[1]
    g = Genotype.from_dataconfig(cfg).complete_from_reference()
    for other in cfg.v_alleles:
        if other != v_gene:
            g.delete_gene(other, haplotype="both", segment="V")
    g.heterozygous(v_gene, a0, a1).with_subject("S1")
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(
        n=400, seed=9, expose_provenance=True
    )
    calls = [r["truth_v_call"] for r in res]
    assert set(calls) <= {a0, a1}, set(calls)
    frac0 = calls.count(a0) / len(calls)
    assert 0.35 < frac0 < 0.65, frac0


def test_vdj_phasing_links_v_and_j_to_one_chromosome():
    """The core phasing guarantee: with hap0 = {V:gv0, J:gj0} and
    hap1 = {V:gv1, J:gj1}, every rearrangement's V and J must come from
    the SAME chromosome — never a cross pairing."""
    cfg = _cfg()
    vgenes, jgenes = list(cfg.v_alleles), list(cfg.j_alleles)
    gv0, gv1 = vgenes[0], vgenes[1]
    gj0, gj1 = jgenes[0], jgenes[1]
    g = Genotype.from_dataconfig(cfg).complete_from_reference()
    # Keep only gv0/gv1 (V) and gj0/gj1 (J); delete everything else.
    for gene in cfg.v_alleles:
        if gene not in (gv0, gv1):
            g.delete_gene(gene, "both", segment="V")
    for gene in cfg.j_alleles:
        if gene not in (gj0, gj1):
            g.delete_gene(gene, "both", segment="J")
    # Phase: hap0 carries gv0 + gj0; hap1 carries gv1 + gj1.
    g.delete_gene(gv0, haplotype=1, segment="V")
    g.delete_gene(gv1, haplotype=0, segment="V")
    g.delete_gene(gj0, haplotype=1, segment="J")
    g.delete_gene(gj1, haplotype=0, segment="J")
    g.with_subject("S1")
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(
        n=200, seed=4, expose_provenance=True
    )
    for r in res:
        v_is_0 = r["truth_v_call"].startswith(gv0 + "*")
        j_is_0 = r["truth_j_call"].startswith(gj0 + "*")
        assert v_is_0 == j_is_0, (r["truth_v_call"], r["truth_j_call"])
        assert (r["haplotype"] == 0) == v_is_0


def test_one_dead_haplotype_uses_the_live_one_under_productive_only():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference()
    # Delete every J gene on haplotype 1 → only haplotype 0 can produce
    # a rearrangement; productive_only must still succeed via haplotype 0.
    for j_gene in cfg.j_alleles:
        g.delete_gene(j_gene, haplotype=1, segment="J")
    g.with_subject("S1")
    res = (
        ga.Experiment.on(cfg)
        .productive_only()
        .with_genotype(g)
        .recombine()
        .run_records(n=40, seed=2)
    )
    assert len(res) == 40
    assert all(r["haplotype"] == 0 for r in res)
