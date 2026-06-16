//! `SampleGenotypePass` — phased, genotype-aware V(D)J allele sampling.
//!
//! Replaces the three flat `SampleAllelePass` passes when a genotype is
//! attached. In ONE pass it:
//!   1. draws the chromosome (haplotype) once, restricted to **viable**
//!      haplotypes (those carrying a feasible allele for every required
//!      segment under the active contracts + feasibility), and
//!   2. for each of V/(D)/J, samples a gene present on that chromosome
//!      (usage-weighted) then the allele within the gene slot
//!      (single-copy slots are deterministic), assigning the slot.
//!
//! The chromosome is a local variable — there is no cross-pass state to
//! share (each pass gets a fresh per-pass trace), which is exactly why
//! haplotype + gene + allele all live in one pass. The canonical
//! `sample_allele.{seg}` choice + the slot assignment are emitted just
//! like the flat path, so every downstream pass (assemble/trim/AIRR/
//! replay) is unchanged. Gene + within-slot choices are recorded as
//! extra addresses for provenance and replay.
use std::sync::Arc;

use crate::address::{self, ChoiceAddress, ChoiceAddressPattern};
use crate::assignment::AlleleInstance;
use crate::contract::ChoiceContext;
use crate::dist::FilteredSampleError;
use crate::genotype::Genotype;
use crate::ir::{Segment, Simulation, SimulationBuilder};
use crate::pass::{AlleleIdSupport, Pass, PassCompileFact, PassContext, PassEffect, PassError};
use crate::refdata::{AlleleId, GeneId};
use crate::rng::Rng;
use crate::trace::ChoiceValue;

pub struct SampleGenotypePass {
    genotype: Arc<Genotype>,
    d_required: bool,
    // Per-segment gene-usage weights (empty => uniform over present genes).
    usage_v: Vec<(GeneId, f64)>,
    usage_d: Vec<(GeneId, f64)>,
    usage_j: Vec<(GeneId, f64)>,
}

impl SampleGenotypePass {
    pub fn new(
        genotype: Arc<Genotype>,
        d_required: bool,
        usage_v: Vec<(GeneId, f64)>,
        usage_d: Vec<(GeneId, f64)>,
        usage_j: Vec<(GeneId, f64)>,
    ) -> Self {
        Self {
            genotype,
            d_required,
            usage_v,
            usage_d,
            usage_j,
        }
    }

    fn segments(&self) -> Vec<Segment> {
        if self.d_required {
            vec![Segment::V, Segment::D, Segment::J]
        } else {
            vec![Segment::V, Segment::J]
        }
    }

    fn usage(&self, seg: Segment) -> &[(GeneId, f64)] {
        match seg {
            Segment::V => &self.usage_v,
            Segment::D => &self.usage_d,
            Segment::J => &self.usage_j,
            _ => &[],
        }
    }

    fn usage_of(&self, seg: Segment, g: GeneId) -> f64 {
        let table = self.usage(seg);
        if table.is_empty() {
            1.0
        } else {
            table
                .iter()
                .find(|(gg, _)| *gg == g)
                .map(|(_, w)| *w)
                .unwrap_or(1.0)
        }
    }

    fn vseg(seg: Segment) -> address::VdjSegment {
        seg.try_into().expect("V/D/J segment")
    }

    fn allele_feasible(
        &self,
        seg: Segment,
        id: AlleleId,
        sim: &Simulation,
        ctx: &PassContext,
    ) -> bool {
        let choice = ChoiceValue::AlleleId(id.index());
        let vseg = Self::vseg(seg);
        let contract_ok = ctx.contracts.map_or(true, |k| {
            k.admits_typed(
                sim,
                ctx.refdata,
                ChoiceContext::none().with_address(ChoiceAddress::SampleAllele(vseg)),
                &choice,
            )
            .is_ok()
        });
        let feasible_ok = ctx.feasibility.map_or(true, |f| {
            f.admits(
                ctx.pass_index,
                sim,
                ctx.refdata,
                address::sample_allele_vdj(seg),
                &choice,
            )
        });
        contract_ok && feasible_ok
    }

    /// Carried alleles on chromosome `c` for `seg` that pass contracts +
    /// feasibility.
    fn feasible_alleles(
        &self,
        c: usize,
        seg: Segment,
        sim: &Simulation,
        ctx: &PassContext,
    ) -> Vec<AlleleId> {
        self.genotype
            .haplotype(c)
            .carried_alleles(seg)
            .into_iter()
            .filter(|id| self.allele_feasible(seg, *id, sim, ctx))
            .collect()
    }

    fn is_viable(&self, c: usize, sim: &Simulation, ctx: &PassContext) -> bool {
        self.segments()
            .iter()
            .all(|seg| !self.feasible_alleles(c, *seg, sim, ctx).is_empty())
    }

    fn viable_set(&self, sim: &Simulation, ctx: &PassContext) -> Vec<usize> {
        (0..2).filter(|&c| self.is_viable(c, sim, ctx)).collect()
    }

    fn draw_haplotype(&self, viable: &[usize], rng: &mut Rng) -> usize {
        let w = self.genotype.chromosome_weights();
        let total: f64 = viable.iter().map(|&c| w[c] as f64).sum();
        if total <= 0.0 {
            return *viable.first().expect("viable non-empty");
        }
        let mut x = rng.next_f64() * total;
        for &c in viable {
            x -= w[c] as f64;
            if x < 0.0 {
                return c;
            }
        }
        *viable.last().expect("viable non-empty")
    }

    fn weighted_pick<T: Copy>(items: &[(T, f64)], rng: &mut Rng) -> T {
        let total: f64 = items.iter().map(|(_, w)| *w).sum();
        if total <= 0.0 {
            return items.first().expect("non-empty").0;
        }
        let mut x = rng.next_f64() * total;
        for (t, w) in items {
            x -= *w;
            if x < 0.0 {
                return *t;
            }
        }
        items.last().expect("non-empty").0
    }

    fn infeasible_error(&self) -> PassError {
        PassError::constraint_sampling(
            self.name(),
            "sample_haplotype",
            FilteredSampleError::EmptyAdmissibleSupport,
        )
    }

    fn commit(&self, seg: Segment, sim: Simulation, id: AlleleId, ctx: &mut PassContext) -> Simulation {
        let mut b = SimulationBuilder::from_simulation(sim);
        if ctx.event_log_sink.is_some() {
            b.attach_event_log_observer();
        }
        b.assign_allele(seg, AlleleInstance::new(id));
        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(b.seal_event_log_observer());
        }
        b.seal()
    }

    /// Union of carried alleles across both haplotypes for a segment —
    /// the support advertised to the feasibility/schedule analyzer.
    fn union_support(&self, seg: Segment) -> Vec<(AlleleId, f64)> {
        let mut ids: Vec<AlleleId> = Vec::new();
        for c in 0..2 {
            for id in self.genotype.haplotype(c).carried_alleles(seg) {
                if !ids.contains(&id) {
                    ids.push(id);
                }
            }
        }
        ids.into_iter().map(|id| (id, 1.0)).collect()
    }

    /// Live (fresh-RNG) sampling of one segment within chromosome `c`.
    fn sample_segment_live(
        &self,
        seg: Segment,
        c: usize,
        sim: Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        let hap = self.genotype.haplotype(c);
        let genes: Vec<(GeneId, f64)> = hap
            .present_genes(seg)
            .filter(|g| {
                hap.slot(seg, *g)
                    .iter()
                    .any(|cp| self.allele_feasible(seg, cp.allele, &sim, ctx))
            })
            .map(|g| (g, self.usage_of(seg, g)))
            .collect();
        if genes.is_empty() {
            return Err(PassError::constraint_sampling(
                self.name(),
                address::sample_allele_vdj(seg),
                FilteredSampleError::EmptyAdmissibleSupport,
            ));
        }
        let gene = Self::weighted_pick(&genes, ctx.rng);

        let slot: Vec<(AlleleId, f64)> = hap
            .slot(seg, gene)
            .iter()
            .filter(|cp| self.allele_feasible(seg, cp.allele, &sim, ctx))
            .map(|cp| (cp.allele, cp.weight as f64 * cp.copies as f64))
            .collect();
        let vseg = Self::vseg(seg);
        let id = if slot.len() == 1 {
            slot[0].0
        } else {
            let chosen = Self::weighted_pick(&slot, ctx.rng);
            ctx.trace.record_choice(
                ChoiceAddress::SampleAlleleInSlot(vseg),
                ChoiceValue::AlleleId(chosen.index()),
            );
            chosen
        };
        ctx.trace
            .record_choice(ChoiceAddress::SampleGene(vseg), ChoiceValue::GeneId(gene.index()));
        ctx.trace
            .record_choice(ChoiceAddress::SampleAllele(vseg), ChoiceValue::AlleleId(id.index()));
        Ok(self.commit(seg, sim, id, ctx))
    }

    /// Replay (trace-injected) sampling of one segment within `c`.
    fn sample_segment_replay(
        &self,
        seg: Segment,
        c: usize,
        sim: Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        let vseg = Self::vseg(seg);
        let gene_idx = ctx
            .replay_cursor
            .as_deref_mut()
            .expect("replay cursor present")
            .expect_gene_id(ChoiceAddress::SampleGene(vseg))
            .map_err(|r| PassError::replay(self.name(), r))?;
        let gene = GeneId::new(gene_idx);
        let slot = self.genotype.haplotype(c).slot(seg, gene);
        let id = if slot.len() == 1 {
            slot[0].allele
        } else {
            let a = ctx
                .replay_cursor
                .as_deref_mut()
                .expect("replay cursor present")
                .expect_allele_id(ChoiceAddress::SampleAlleleInSlot(vseg))
                .map_err(|r| PassError::replay(self.name(), r))?;
            ctx.trace
                .record_choice(ChoiceAddress::SampleAlleleInSlot(vseg), ChoiceValue::AlleleId(a));
            AlleleId::new(a)
        };
        ctx.trace
            .record_choice(ChoiceAddress::SampleGene(vseg), ChoiceValue::GeneId(gene_idx));
        ctx.trace
            .record_choice(ChoiceAddress::SampleAllele(vseg), ChoiceValue::AlleleId(id.index()));
        Ok(self.commit(seg, sim, id, ctx))
    }
}

impl Pass for SampleGenotypePass {
    fn name(&self) -> &str {
        "sample_genotype"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_checked(sim, ctx)
            .expect("SampleGenotypePass permissive execution must not error")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        // Decide the chromosome (replay consumes; live draws among viable).
        let c = if ctx.replay_cursor.is_some() {
            let recorded = ctx
                .replay_cursor
                .as_deref_mut()
                .expect("replay cursor present")
                .expect_haplotype(ChoiceAddress::SampleHaplotype)
                .map_err(|r| PassError::replay(self.name(), r))?;
            let viable = self.viable_set(sim, ctx);
            if !viable.contains(&(recorded as usize)) {
                return Err(self.infeasible_error());
            }
            ctx.trace
                .record_choice(ChoiceAddress::SampleHaplotype, ChoiceValue::Haplotype(recorded));
            recorded as usize
        } else {
            let viable = self.viable_set(sim, ctx);
            if viable.is_empty() {
                return Err(self.infeasible_error());
            }
            let c = self.draw_haplotype(&viable, ctx.rng);
            ctx.trace
                .record_choice(ChoiceAddress::SampleHaplotype, ChoiceValue::Haplotype(c as u8));
            c
        };

        let mut current = sim.clone();
        let replaying = ctx.replay_cursor.is_some();
        for seg in self.segments() {
            current = if replaying {
                self.sample_segment_replay(seg, c, current, ctx)?
            } else {
                self.sample_segment_live(seg, c, current, ctx)?
            };
        }
        Ok(current)
    }

    fn declared_choice_patterns(&self) -> Vec<ChoiceAddressPattern> {
        let mut patterns = vec![ChoiceAddressPattern::SampleHaplotype];
        for seg in self.segments() {
            let vseg = Self::vseg(seg);
            patterns.push(ChoiceAddressPattern::SampleGene(vseg));
            patterns.push(ChoiceAddressPattern::SampleAlleleInSlot(vseg));
            patterns.push(ChoiceAddressPattern::SampleAllele(vseg));
        }
        patterns
    }

    fn effects(&self) -> Vec<PassEffect> {
        self.segments()
            .into_iter()
            .map(PassEffect::AssignAllele)
            .collect()
    }

    fn compile_facts(&self) -> Vec<PassCompileFact> {
        self.segments()
            .into_iter()
            .map(|seg| PassCompileFact::AlleleSampleSupport {
                segment: seg,
                support: AlleleIdSupport::from_weighted_pairs(Some(self.union_support(seg))),
            })
            .collect()
    }
}

#[cfg(test)]
pub(crate) mod test_support {
    use super::*;
    use crate::genotype::{GeneCopy, Haplotype};

    pub fn copy(id: u32) -> GeneCopy {
        GeneCopy {
            allele: AlleleId::new(id),
            copies: 1,
            weight: 1.0,
        }
    }

    /// hap0 carries V (gene0 -> allele 0) + J (gene0 -> allele 100);
    /// hap1 carries V but no J → only hap0 viable.
    pub fn geno_chrom1_deletes_j() -> Genotype {
        let mut h0 = Haplotype::new();
        h0.set(Segment::V, GeneId::new(0), vec![copy(0)]);
        h0.set(Segment::J, GeneId::new(0), vec![copy(100)]);
        let mut h1 = Haplotype::new();
        h1.set(Segment::V, GeneId::new(0), vec![copy(1)]);
        Genotype::new([h0, h1], [0.5, 0.5], Some("S1".into()), "sha256:test".into())
    }

    /// Neither haplotype carries a J gene → no viable haplotype.
    pub fn geno_both_delete_j() -> Genotype {
        let mut h0 = Haplotype::new();
        h0.set(Segment::V, GeneId::new(0), vec![copy(0)]);
        let mut h1 = Haplotype::new();
        h1.set(Segment::V, GeneId::new(0), vec![copy(1)]);
        Genotype::new([h0, h1], [0.5, 0.5], Some("S1".into()), "sha256:test".into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;

    #[test]
    fn draws_only_the_viable_haplotype_and_assigns_carried_alleles() {
        let g = Arc::new(test_support::geno_chrom1_deletes_j());
        let pass = SampleGenotypePass::new(g, false, vec![], vec![], vec![]);
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));
        for seed in 0..30u64 {
            let outcome = PassRuntime::execute(&plan, Simulation::new(), seed);
            match outcome.trace.find("sample_haplotype").unwrap().value {
                ChoiceValue::Haplotype(c) => assert_eq!(c, 0, "seed {seed}"),
                _ => panic!("wrong variant"),
            }
            let sim = outcome.final_simulation();
            assert_eq!(
                sim.assignments.get(Segment::V).unwrap().allele_id,
                AlleleId::new(0)
            );
            assert_eq!(
                sim.assignments.get(Segment::J).unwrap().allele_id,
                AlleleId::new(100)
            );
            // single-copy slots => no within-slot record
            assert!(outcome.trace.find("sample_allele_in_slot.v").is_none());
        }
    }

    #[test]
    fn records_canonical_sample_allele_and_gene_addresses() {
        let g = Arc::new(test_support::geno_chrom1_deletes_j());
        let pass = SampleGenotypePass::new(g, false, vec![], vec![], vec![]);
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);
        match outcome.trace.find("sample_allele.v").unwrap().value {
            ChoiceValue::AlleleId(id) => assert_eq!(id, 0),
            _ => panic!("expected AlleleId at sample_allele.v"),
        }
        match outcome.trace.find("sample_gene.v").unwrap().value {
            ChoiceValue::GeneId(g) => assert_eq!(g, 0),
            _ => panic!("expected GeneId at sample_gene.v"),
        }
    }

    #[test]
    fn errors_when_no_haplotype_viable() {
        let g = Arc::new(test_support::geno_both_delete_j());
        let pass = SampleGenotypePass::new(g, false, vec![], vec![], vec![]);
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));
        let result =
            PassRuntime::execute_strict_with_context(&plan, Simulation::new(), 0, None, None);
        assert!(result.is_err(), "expected genotype-infeasibility error");
    }
}
