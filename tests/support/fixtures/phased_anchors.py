"""A cartridge whose V and J alleles sit at known reading-frame phases.

Used where behavior depends on the phase an anchor lands in rather than on
the allele itself.
"""
from __future__ import annotations

import GenAIRR as ga
from GenAIRR.alleles import JAllele, VAllele
from GenAIRR.genotype import Genotype
from tests.support.cases.productive import synthetic_cartridge

def _phased_fixture():
    base = synthetic_cartridge({"V_3": [(1, 1)]})
    v0 = VAllele("IGKV-PHASE0*01", "TGT", 3, anchor_override=0)
    v2 = VAllele("IGKV-PHASE2*01", "TGTCC", 5, anchor_override=0)
    j0 = JAllele("IGKJ-PHASE0*01", "TTTGGG", 6, anchor_override=0)
    j2 = JAllele("IGKJ-PHASE2*01", "CCTTTGGG", 8, anchor_override=2)
    cfg = ga.DataConfig(
        name="productive-p-phased-fixture", metadata=base.metadata,
        v_alleles={a.gene: [a] for a in (v0, v2)},
        j_alleles={a.gene: [a] for a in (j0, j2)},
        reference_models=base.reference_models,
    )
    genotype = (
        Genotype.from_dataconfig(cfg).complete_from_reference()
        .delete_gene(v2.gene, haplotype=0, segment="V")
        .delete_gene(j2.gene, haplotype=0, segment="J")
        .delete_gene(j0.gene, haplotype=1, segment="J")
    )
    # P=1,N=0 allows V0+J2 and V2+J0. Hap0={V0,J0} is
    # impossible; hap1={V0,V2,J2} permits V0 only. Marginal allele
    # support across both chromosomes incorrectly accepts both hap0
    # and hap1's V2, then fails after earlier choices are committed.
    return cfg, genotype, v0.name, j2.name
