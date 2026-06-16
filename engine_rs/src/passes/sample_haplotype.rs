//! `SampleHaplotypePass` — phased-genotype chromosome draw.
//!
//! Runs once per rearrangement, before the per-segment gene/allele
//! passes. It picks the chromosome (0 or 1) the rearrangement will draw
//! its V/D/J from, restricted to **viable** haplotypes — those that
//! carry at least one feasible allele for every required segment under
//! the active contracts + feasibility. The chosen chromosome is recorded
//! to the trace at `sample_haplotype`; the gene/allele passes read it
//! back from there (there is no per-simulation scratch state).
use std::sync::Arc;

use crate::address::{self, ChoiceAddress};
use crate::contract::ChoiceContext;
use crate::dist::FilteredSampleError;
use crate::genotype::Genotype;
use crate::ir::{Segment, Simulation};
use crate::pass::{Pass, PassContext, PassError};
use crate::refdata::AlleleId;
use crate::rng::Rng;
use crate::trace::ChoiceValue;

pub struct SampleHaplotypePass {
    genotype: Arc<Genotype>,
    d_required: bool,
}

impl SampleHaplotypePass {
    pub fn new(genotype: Arc<Genotype>, d_required: bool) -> Self {
        Self {
            genotype,
            d_required,
        }
    }

    fn choice_address(&self) -> ChoiceAddress {
        ChoiceAddress::SampleHaplotype
    }

    /// Carried alleles on chromosome `c` for `seg` that pass the active
    /// contracts + feasibility. With neither active, all carried alleles
    /// are admissible.
    fn feasible_alleles(
        &self,
        c: usize,
        seg: Segment,
        sim: &Simulation,
        ctx: &PassContext,
    ) -> Vec<AlleleId> {
        let carried = self.genotype.haplotype(c).carried_alleles(seg);
        let addr = address::sample_allele_vdj(seg);
        let vseg: address::VdjSegment = seg.try_into().expect("V/D/J segment");
        carried
            .into_iter()
            .filter(|id| {
                let choice = ChoiceValue::AlleleId(id.index());
                let contract_ok = ctx.contracts.map_or(true, |k| {
                    k.admits_typed(
                        sim,
                        ctx.refdata,
                        ChoiceContext::none()
                            .with_address(ChoiceAddress::SampleAllele(vseg)),
                        &choice,
                    )
                    .is_ok()
                });
                let feasible_ok = ctx.feasibility.map_or(true, |feas| {
                    feas.admits(ctx.pass_index, sim, ctx.refdata, addr, &choice)
                });
                contract_ok && feasible_ok
            })
            .collect()
    }

    fn is_viable(&self, c: usize, sim: &Simulation, ctx: &PassContext) -> bool {
        if self.feasible_alleles(c, Segment::V, sim, ctx).is_empty() {
            return false;
        }
        if self.feasible_alleles(c, Segment::J, sim, ctx).is_empty() {
            return false;
        }
        if self.d_required && self.feasible_alleles(c, Segment::D, sim, ctx).is_empty() {
            return false;
        }
        true
    }

    fn viable_set(&self, sim: &Simulation, ctx: &PassContext) -> Vec<usize> {
        (0..2).filter(|&c| self.is_viable(c, sim, ctx)).collect()
    }

    fn draw_from(&self, viable: &[usize], rng: &mut Rng) -> usize {
        let weights = self.genotype.chromosome_weights();
        let total: f64 = viable.iter().map(|&c| weights[c] as f64).sum();
        if total <= 0.0 {
            return *viable.first().expect("viable non-empty");
        }
        let mut x = rng.next_f64() * total;
        for &c in viable {
            x -= weights[c] as f64;
            if x < 0.0 {
                return c;
            }
        }
        *viable.last().expect("viable non-empty")
    }

    fn infeasible_error(&self) -> PassError {
        // Surfaced as a constraint-sampling error at `sample_haplotype`;
        // the message names the empty admissible support. The subject id
        // is included in the address-bearing diagnostics downstream.
        PassError::constraint_sampling(
            "sample_haplotype",
            "sample_haplotype",
            FilteredSampleError::EmptyAdmissibleSupport,
        )
    }
}

impl Pass for SampleHaplotypePass {
    fn name(&self) -> &str {
        "sample_haplotype"
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_checked(sim, ctx)
            .expect("SampleHaplotypePass permissive execution must not error")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        // Replay: consume + revalidate the recorded chromosome.
        if ctx.replay_cursor.is_some() {
            let c = ctx
                .replay_cursor
                .as_deref_mut()
                .expect("replay cursor present")
                .expect_haplotype(self.choice_address())
                .map_err(|r| PassError::replay(self.name(), r))?;
            let viable = self.viable_set(sim, ctx);
            if !viable.contains(&(c as usize)) {
                return Err(self.infeasible_error());
            }
            ctx.trace
                .record_choice(self.choice_address(), ChoiceValue::Haplotype(c));
            return Ok(sim.clone());
        }

        let viable = self.viable_set(sim, ctx);
        if viable.is_empty() {
            return Err(self.infeasible_error());
        }
        let c = self.draw_from(&viable, ctx.rng);
        ctx.trace
            .record_choice(self.choice_address(), ChoiceValue::Haplotype(c as u8));
        Ok(sim.clone())
    }

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![address::ChoiceAddressPattern::SampleHaplotype]
    }
}

#[cfg(test)]
pub(crate) mod test_support {
    use super::*;
    use crate::genotype::{GeneCopy, Haplotype};
    use crate::refdata::GeneId;

    fn copy(id: u32) -> GeneCopy {
        GeneCopy {
            allele: AlleleId::new(id),
            copies: 1,
            weight: 1.0,
        }
    }

    /// hap0 carries V gene0 + J gene0; hap1 carries V gene0 but DELETES
    /// J → only hap0 is viable (J-less hap1 can't make a rearrangement).
    pub fn geno_chrom1_deletes_j() -> Genotype {
        let mut h0 = Haplotype::new();
        h0.set(Segment::V, GeneId::new(0), vec![copy(0)]);
        h0.set(Segment::J, GeneId::new(0), vec![copy(2)]);
        let mut h1 = Haplotype::new();
        h1.set(Segment::V, GeneId::new(0), vec![copy(1)]);
        // J intentionally absent on h1.
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
    use crate::ir::Simulation;
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;

    #[test]
    fn draws_only_the_viable_haplotype_when_one_is_dead() {
        let geno = test_support::geno_chrom1_deletes_j();
        let pass = SampleHaplotypePass::new(Arc::new(geno), false);
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));
        for seed in 0..50u64 {
            let outcome = PassRuntime::execute(&plan, Simulation::new(), seed);
            match outcome.trace.find("sample_haplotype").unwrap().value {
                ChoiceValue::Haplotype(c) => assert_eq!(c, 0, "seed {seed}"),
                _ => panic!("wrong variant at sample_haplotype"),
            }
        }
    }

    #[test]
    fn errors_when_no_haplotype_viable() {
        let geno = test_support::geno_both_delete_j();
        let pass = SampleHaplotypePass::new(Arc::new(geno), false);
        let mut plan = PassPlan::new();
        plan.push(Box::new(pass));
        let result =
            PassRuntime::execute_strict_with_context(&plan, Simulation::new(), 0, None, None);
        assert!(result.is_err(), "expected genotype-infeasibility error");
    }
}
