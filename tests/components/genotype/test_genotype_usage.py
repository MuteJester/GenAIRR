"""Genotype gene-usage wiring is active when the cartridge authors a
typed ``reference_models.allele_usage`` plane (review #3).

This pins that the Rust gene-level usage aggregation is *reachable* and
*effective* — not just code-complete. Bundled configs that don't author
allele_usage fall back to uniform-over-present-genes (× copy dosage),
which is documented, not tested here.
"""
import dataclasses

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype
from GenAIRR.reference_models import AlleleUsageSpec, ReferenceEmpiricalModels


def _cfg_with_v_usage(target_allele_name, weight=1000.0):
    base = gdata.HUMAN_IGH_OGRDB
    rm = ReferenceEmpiricalModels(allele_usage=AlleleUsageSpec(v={target_allele_name: weight}))
    return dataclasses.replace(base, reference_models=rm)


def test_typed_allele_usage_biases_genotype_gene_choice():
    base = gdata.HUMAN_IGH_OGRDB
    # Heavily weight the first allele of one V gene.
    target_gene = list(base.v_alleles)[10]
    target_allele = base.v_alleles[target_gene][0].name
    cfg = _cfg_with_v_usage(target_allele, weight=1000.0)

    g = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("S1")
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(
        n=300, seed=1, expose_provenance=True
    )
    frac_target = sum(
        1 for r in res if r["truth_v_call"].startswith(target_gene + "*")
    ) / len(res)
    # With ~52 V genes uniform would give ~0.02; a 1000x usage weight on
    # this gene must dominate.
    assert frac_target > 0.5, frac_target
