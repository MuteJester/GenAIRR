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
        .complete_from_reference("homozygous_common")
        .homozygous(v_gene, carried)
        .with_subject("S1")
    )
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(n=300, seed=1)
    seen_for_gene = {
        r["v_call"] for r in res if r["v_call"].startswith(v_gene + "*")
    }
    # Only the carried allele of that gene may appear (never names[0]).
    assert seen_for_gene <= {carried}, seen_for_gene


def test_phased_run_is_deterministic_under_same_seed():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("S1")
    exp = ga.Experiment.on(cfg).with_genotype(g).recombine()
    a = exp.run_records(n=40, seed=77)
    b = exp.run_records(n=40, seed=77)
    assert [r["v_call"] for r in a] == [r["v_call"] for r in b]
    assert [r["j_call"] for r in a] == [r["j_call"] for r in b]
