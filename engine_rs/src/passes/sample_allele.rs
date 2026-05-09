//! `SampleAllelePass` — recombination-stage allele sampling (C.5).

use crate::assignment::AlleleInstance;
use crate::dist::{sample_filtered_result, Distribution, FilteredSampleError};
use crate::ir::{Segment, Simulation};
use crate::pass::{AlleleIdSupport, Pass, PassCompileFact, PassContext, PassEffect, PassError};
use crate::refdata::AlleleId;
use crate::trace::ChoiceValue;

/// Sample one allele from a pool distribution and assign it to the
/// simulation's V / D / J slot.
///
/// The pass is parameterized by:
/// - **Segment** — must be `V`, `D`, or `J` (NP segments and `C`
///   are not allowed; assemble passes don't run for NP, and
///   constant-region sampling is a future concern).
/// - **Distribution** — any `Box<dyn Distribution<Output = AlleleId>>`,
///   most commonly an `AllelePoolDist` (C.3) constructed against
///   the segment's pool in the active `RefDataConfig`. Construction
///   discipline guarantees every sampled `AlleleId` is in-bounds
///   for that pool.
///
/// On execute the pass:
/// 1. Draws one `AlleleId` from the distribution via `ctx.rng`.
/// 2. Records `ChoiceValue::AlleleId(id.index())` to the trace at
///    address `"sample_allele.{segment}"`.
/// 3. Constructs a fresh `AlleleInstance` (zero trims) and assigns
///    it to the simulation's slot for `segment`.
///
/// The pass does *not* read allele bases — that work belongs to
/// the assembly pass (C.8). Sampling here only chooses the id.
pub struct SampleAllelePass {
    segment: Segment,
    distribution: Box<dyn Distribution<Output = AlleleId>>,
}

impl SampleAllelePass {
    /// Construct a sampling pass for the given segment.
    ///
    /// Panics if `segment` is anything other than V, D, or J.
    /// The constructor catches the misuse at plan-build time so
    /// the pass-execute path is panic-free for valid plans.
    pub fn new(segment: Segment, distribution: Box<dyn Distribution<Output = AlleleId>>) -> Self {
        match segment {
            Segment::V | Segment::D | Segment::J => {}
            _ => panic!(
                "SampleAllelePass: segment must be V, D, or J — got {:?}",
                segment
            ),
        }
        Self {
            segment,
            distribution,
        }
    }

    /// Inspect the configured segment.
    pub fn segment(&self) -> Segment {
        self.segment
    }

    /// The hierarchical-string address (D3) at which this pass
    /// records its choice. Same string as `name()` since the pass
    /// makes exactly one choice per execution.
    fn address(&self) -> &'static str {
        match self.segment {
            Segment::V => "sample_allele.v",
            Segment::D => "sample_allele.d",
            Segment::J => "sample_allele.j",
            // Unreachable due to constructor validation; if this
            // fires it's a code defect, not a user error.
            _ => unreachable!("SampleAllelePass with non-V/D/J segment"),
        }
    }

    fn sample_allele(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<AlleleId, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;
        let feasibility = ctx.feasibility;
        let pass_index = ctx.pass_index;

        if contracts.is_some() || feasibility.is_some() {
            match sample_filtered_result(ctx.rng, self.distribution.as_ref(), |candidate| {
                let choice = ChoiceValue::AlleleId(candidate.index());
                let contract_ok = contracts.map_or(true, |contracts| {
                    contracts
                        .admits(sim, refdata, self.address(), &choice)
                        .is_ok()
                });
                let feasible_ok = feasibility.map_or(true, |feasibility| {
                    feasibility.admits(pass_index, sim, refdata, self.address(), &choice)
                });
                contract_ok && feasible_ok
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if strict => {
                    return Err(self.constraint_sampling_error(reason));
                }
                Err(_) => {}
            }
        }

        Ok(self.distribution.sample(ctx.rng))
    }

    fn constraint_sampling_error(&self, reason: FilteredSampleError) -> PassError {
        PassError::constraint_sampling(self.address(), self.address(), reason)
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        let id = self.sample_allele(sim, ctx, strict)?;
        ctx.trace
            .record(self.address(), ChoiceValue::AlleleId(id.index()));
        Ok(sim.with_allele_assigned(self.segment, AlleleInstance::new(id)))
    }
}

impl Pass for SampleAllelePass {
    fn name(&self) -> &str {
        self.address()
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("SampleAllelePass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choices(&self) -> Vec<String> {
        vec![self.address().to_string()]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::AssignAllele(self.segment)]
    }

    fn compile_facts(&self) -> Vec<PassCompileFact> {
        vec![PassCompileFact::AlleleSampleSupport {
            segment: self.segment,
            support: AlleleIdSupport::from_weighted_pairs(self.distribution.support()),
        }]
    }
}

#[cfg(test)]
pub(super) mod test_support {
    use super::*;
    use crate::refdata::{Allele, AllelePool};

    /// Build an allele pool of `n` named alleles for testing.
    pub fn make_test_pool(n: usize, segment: Segment) -> AllelePool {
        let mut p = AllelePool::new();
        for i in 0..n {
            let _ = p.push(Allele {
                name: format!("test_allele_{}*01", i),
                gene: format!("test_allele_{}", i),
                seq: vec![b'A'; 30],
                segment,
                anchor: Some(10),
            });
        }
        p
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::make_test_pool;
    use super::*;
    use crate::dist::AllelePoolDist;
    use crate::pass::{PassPlan, PassRuntime};

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn sample_allele_pass_rejects_np1() {
        let pool = make_test_pool(1, Segment::V);
        let _ = SampleAllelePass::new(Segment::Np1, Box::new(AllelePoolDist::uniform(&pool)));
    }

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn sample_allele_pass_rejects_np2() {
        let pool = make_test_pool(1, Segment::V);
        let _ = SampleAllelePass::new(Segment::Np2, Box::new(AllelePoolDist::uniform(&pool)));
    }

    #[test]
    fn sample_allele_pass_assigns_to_correct_slot_for_v() {
        let pool = make_test_pool(1, Segment::V);
        let dist = Box::new(AllelePoolDist::uniform(&pool));
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(Segment::V, dist)));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 42);

        let final_sim = outcome.final_simulation();
        // Single-allele dist always returns AlleleId(0).
        assert_eq!(
            final_sim.assignments.v.unwrap().allele_id,
            crate::refdata::AlleleId::new(0)
        );
        assert!(final_sim.assignments.d.is_none());
        assert!(final_sim.assignments.j.is_none());
        // Default trims.
        assert_eq!(final_sim.assignments.v.unwrap().trim_5, 0);
        assert_eq!(final_sim.assignments.v.unwrap().trim_3, 0);
    }

    #[test]
    fn sample_allele_pass_records_to_trace_at_segment_address() {
        let v_pool = make_test_pool(1, Segment::V);
        let d_pool = make_test_pool(1, Segment::D);
        let j_pool = make_test_pool(1, Segment::J);

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::uniform(&d_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&j_pool)),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        // Three choices, one per pass, at the canonical addresses.
        assert_eq!(outcome.trace.len(), 3);
        for addr in ["sample_allele.v", "sample_allele.d", "sample_allele.j"] {
            let rec = outcome
                .trace
                .find(addr)
                .unwrap_or_else(|| panic!("missing {}", addr));
            match rec.value {
                ChoiceValue::AlleleId(id) => {
                    assert_eq!(id, 0); // single-allele pool
                }
                _ => panic!("wrong variant at {}", addr),
            }
        }
    }

    #[test]
    fn sample_allele_pass_is_deterministic_under_same_seed() {
        let pool = make_test_pool(10, Segment::V);
        let mut plan_a = PassPlan::new();
        plan_a.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&pool)),
        )));
        let mut plan_b = PassPlan::new();
        plan_b.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&pool)),
        )));

        let oa = PassRuntime::execute(&plan_a, Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan_b, Simulation::new(), 0xc0ff_ee);

        assert_eq!(
            oa.final_simulation().assignments.v.unwrap().allele_id,
            ob.final_simulation().assignments.v.unwrap().allele_id
        );
        assert_eq!(oa.trace.choices()[0].value, ob.trace.choices()[0].value);
    }

    #[test]
    fn sample_allele_pass_full_recombination_chain_for_vdj() {
        // V + D + J sampling for a heavy chain: all three
        // assignments should be populated after the plan runs.
        let v_pool = make_test_pool(5, Segment::V);
        let d_pool = make_test_pool(3, Segment::D);
        let j_pool = make_test_pool(2, Segment::J);

        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::D,
            Box::new(AllelePoolDist::uniform(&d_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&j_pool)),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 42);

        let sim = outcome.final_simulation();
        assert!(sim.assignments.v.is_some());
        assert!(sim.assignments.d.is_some());
        assert!(sim.assignments.j.is_some());
        assert!(sim.assignments.c.is_none());

        // Sampled ids are in-bounds for their pools (D-binding from C.3).
        assert!(sim.assignments.v.unwrap().allele_id.as_usize() < 5);
        assert!(sim.assignments.d.unwrap().allele_id.as_usize() < 3);
        assert!(sim.assignments.j.unwrap().allele_id.as_usize() < 2);
    }

    #[test]
    fn sample_allele_pass_declared_choices_returns_address() {
        let pool = make_test_pool(1, Segment::V);
        let pass_v = SampleAllelePass::new(Segment::V, Box::new(AllelePoolDist::uniform(&pool)));
        assert_eq!(
            pass_v.declared_choices(),
            vec!["sample_allele.v".to_string()]
        );

        let d_pool = make_test_pool(1, Segment::D);
        let pass_d = SampleAllelePass::new(Segment::D, Box::new(AllelePoolDist::uniform(&d_pool)));
        assert_eq!(
            pass_d.declared_choices(),
            vec!["sample_allele.d".to_string()]
        );
    }

    #[test]
    fn sample_allele_pass_segment_accessor() {
        let pool = make_test_pool(1, Segment::J);
        let pass = SampleAllelePass::new(Segment::J, Box::new(AllelePoolDist::uniform(&pool)));
        assert_eq!(pass.segment(), Segment::J);
        assert_eq!(pass.name(), "sample_allele.j");
    }
}
