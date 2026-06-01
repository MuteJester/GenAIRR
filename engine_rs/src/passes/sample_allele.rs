//! `SampleAllelePass` — recombination-stage allele sampling (C.5).

use crate::address;
use crate::assignment::AlleleInstance;
use crate::contract::ChoiceContext;
use crate::dist::{sample_filtered_result, Distribution, FilteredSampleError};
use crate::ir::{Segment, Simulation, SimulationBuilder};
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
        address::sample_allele_vdj(self.segment)
    }

    fn typed_segment(&self) -> address::VdjSegment {
        match self.segment {
            Segment::V => address::VdjSegment::V,
            Segment::D => address::VdjSegment::D,
            Segment::J => address::VdjSegment::J,
            Segment::Np1 | Segment::Np2 => unreachable!("constructor rejects NP segments"),
        }
    }

    fn choice_address(&self) -> address::ChoiceAddress {
        address::ChoiceAddress::SampleAllele(self.typed_segment())
    }

    /// Constraint-aware allele draw.
    ///
    /// **v3.0 documented exception to the global
    /// constrain-before-propose invariant.** The other constrained
    /// samplers — TrimPass, GenerateNP (length and base),
    /// ContaminantPass, the indel tuple sampler, the substitution
    /// helpers — all skip the slot on permissive empty-support
    /// rather than falling back to an unconstrained draw. Allele
    /// sampling cannot do the same: every downstream pass
    /// (assembly, trim, NP generation) requires an `AlleleInstance`
    /// in the segment slot, so a "skip" semantics would panic the
    /// pipeline at the next assemble. The architecturally honest
    /// alternatives are:
    /// 1. Strict-error in *both* modes (forces upstream
    ///    re-planning).
    /// 2. Sentinel `AlleleId(0)` (semantically wrong — the chosen
    ///    allele is arbitrary, not "the first one").
    /// 3. Keep the unconstrained fallback in permissive mode,
    ///    document it as the exception.
    ///
    /// Option 3 is the lowest-blast-radius choice and what's
    /// implemented here. Strict mode honors the rule (returns
    /// `ConstraintSampling`); permissive mode falls back to the
    /// raw natural draw. A future v3.x phase can revisit this if
    /// allele feasibility filtering becomes more load-bearing.
    /// Contracts may narrow allele IDs (for example
    /// `AnchorPreserved` rejects candidates whose anchors cannot be
    /// retained), and the runtime feasibility context may narrow
    /// them further.
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
                        .admits_typed(
                            sim,
                            refdata,
                            ChoiceContext::none().with_address(self.choice_address()),
                            &choice,
                        )
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
                // v3.0 documented exception: see method-level doc
                // above. Allele assignment is load-bearing for
                // downstream passes; a permissive skip would panic
                // the pipeline.
                Err(_) => {}
            }
        }

        Ok(self.distribution.sample(ctx.rng))
    }

    fn constraint_sampling_error(&self, reason: FilteredSampleError) -> PassError {
        PassError::constraint_sampling(self.address(), self.address(), reason)
    }

    /// Apply the chosen allele assignment via a `SimulationBuilder`
    /// so the [`SimulationEvent::AssignmentChanged`] event flows
    /// through every attached sink. The persistent mutation
    /// itself is still `Simulation::with_allele_assigned`, called
    /// inside [`SimulationBuilder::assign_allele`]; this wrap is
    /// purely about routing the consequence onto the typed event
    /// channel.
    ///
    /// When `ctx.event_log_sink` is supplied (compiled execution
    /// path), the captured event is forwarded into it so the
    /// pass's `EventRecord` carries the consequence on
    /// `simulation_events`.
    fn commit_assignment(
        &self,
        sim: Simulation,
        id: AlleleId,
        ctx: &mut PassContext,
    ) -> Simulation {
        let mut builder = SimulationBuilder::from_simulation(sim);
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }
        builder.assign_allele(self.segment, AlleleInstance::new(id));
        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(builder.seal_event_log_observer());
        }
        builder.seal()
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        _strict: bool,
    ) -> Result<Simulation, PassError> {
        // Trace-injected replay (Option B): consume the recorded
        // allele id from the cursor, then run the same admissibility
        // chain a fresh draw would have to pass. The replay contract
        // is **trace supplies proposals; engine validation decides
        // whether they apply** — a recorded id that no longer
        // belongs in refdata, in the distribution's support, or
        // under the active contracts / feasibility is a structured
        // error, not silently-applied corruption.
        if ctx.replay_cursor.is_some() {
            let id_index = ctx
                .replay_cursor
                .as_deref_mut()
                .expect("replay_cursor presence checked above")
                .expect_allele_id(self.choice_address())
                .map_err(|reason| PassError::replay(self.name(), reason))?;
            self.validate_replay_candidate(sim, ctx, id_index)?;
            ctx.trace
                .record_choice(self.choice_address(), ChoiceValue::AlleleId(id_index));
            return Ok(self.commit_assignment(sim.clone(), AlleleId::new(id_index), ctx));
        }

        let id = self.sample_allele(sim, ctx, _strict)?;
        ctx.trace
            .record_choice(self.choice_address(), ChoiceValue::AlleleId(id.index()));
        Ok(self.commit_assignment(sim.clone(), id, ctx))
    }

    /// Run the full admissibility chain against a replay candidate:
    ///
    /// 1. The allele id must resolve in `ctx.refdata` for this
    ///    segment (when refdata is supplied).
    /// 2. The id must lie in the natural distribution's enumerable
    ///    support — a recorded value outside the support implies
    ///    either a refdata swap or a tampered trace.
    /// 3. Any active contracts must admit the candidate via
    ///    `admits_typed` (same surface as the fresh-sampling path's
    ///    contract filter).
    /// 4. The runtime feasibility tracker must admit the candidate
    ///    (same call as the fresh-sampling path).
    ///
    /// Mismatches surface as the structured `PassError` variant
    /// that best names *what* failed:
    /// - unknown allele id → `MissingAllele`
    /// - out-of-support / feasibility-rejected →
    ///   `InvalidDistributionOutput` with a specific `reason`
    ///   string (so replay diagnostics can distinguish them from
    ///   live-sampling failures by tag, not by error variant)
    /// - contract rejection → `ContractViolation` wrapping the
    ///   single violation `admits_typed` short-circuited on
    ///
    /// When `support()` returns `None` (distribution too large to
    /// enumerate), step 2 is skipped — there's no usable oracle.
    /// Today's `AllelePoolDist` always enumerates, so this branch
    /// is reachable only for future distribution kinds.
    fn validate_replay_candidate(
        &self,
        sim: &Simulation,
        ctx: &PassContext,
        id_index: u32,
    ) -> Result<(), PassError> {
        let allele_id = AlleleId::new(id_index);
        let choice = ChoiceValue::AlleleId(id_index);

        if let Some(refdata) = ctx.refdata {
            if refdata.get(self.segment, allele_id).is_none() {
                return Err(PassError::missing_allele(
                    self.name(),
                    self.segment,
                    id_index,
                ));
            }
        }

        if let Some(support) = self.distribution.support() {
            let in_support = support
                .iter()
                .any(|(id, weight)| id.index() == id_index && *weight > 0.0);
            if !in_support {
                return Err(PassError::invalid_distribution_output(
                    self.name(),
                    self.address(),
                    id_index as i64,
                    "allele_id_not_in_distribution_support",
                ));
            }
        }

        if let Some(contracts) = ctx.contracts {
            contracts
                .admits_typed(
                    sim,
                    ctx.refdata,
                    ChoiceContext::none().with_address(self.choice_address()),
                    &choice,
                )
                .map_err(|violation| {
                    PassError::contract_violation(self.name(), vec![violation])
                })?;
        }

        if let Some(feasibility) = ctx.feasibility {
            if !feasibility.admits(
                ctx.pass_index,
                sim,
                ctx.refdata,
                self.address(),
                &choice,
            ) {
                return Err(PassError::invalid_distribution_output(
                    self.name(),
                    self.address(),
                    id_index as i64,
                    "feasibility_rejected",
                ));
            }
        }

        Ok(())
    }
}

impl Pass for SampleAllelePass {
    fn name(&self) -> &str {
        self.address()
    }

    fn parameter_signature(&self) -> String {
        // The sampling distribution is opaque (typically the pool's
        // uniform `AllelePoolDist`). We deliberately do NOT fold
        // the distribution's `support()` into the signature: the
        // pool identity is already covered by
        // `refdata_content_hash`, and emitting the entire allele-id
        // ↔ weight table per V/D/J pass would inflate signatures
        // without adding information beyond what the content hash
        // already gates. The pass's segment is encoded in
        // `name()` (`sample_allele.v` / `.d` / `.j`), so no
        // additional parameter signature is needed.
        String::new()
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

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![address::ChoiceAddressPattern::SampleAllele(
            self.typed_segment(),
        )]
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
                functional_status: None,
                subregions: Vec::new(),
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
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;

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
            final_sim
                .assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .allele_id,
            crate::refdata::AlleleId::new(0)
        );
        assert!(final_sim.assignments.get(Segment::D).is_none());
        assert!(final_sim.assignments.get(Segment::J).is_none());
        // Default trims.
        assert_eq!(
            final_sim
                .assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .trim_5,
            0
        );
        assert_eq!(
            final_sim
                .assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .trim_3,
            0
        );
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
            oa.final_simulation()
                .assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .allele_id,
            ob.final_simulation()
                .assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .allele_id
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
        assert!(sim.assignments.has(Segment::V));
        assert!(sim.assignments.has(Segment::D));
        assert!(sim.assignments.has(Segment::J));
        assert_eq!(sim.assignments.iter().count(), 3);

        // Sampled ids are in-bounds for their pools (D-binding from C.3).
        assert!(
            sim.assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .allele_id
                .as_usize()
                < 5
        );
        assert!(
            sim.assignments
                .get(Segment::D)
                .copied()
                .unwrap()
                .allele_id
                .as_usize()
                < 3
        );
        assert!(
            sim.assignments
                .get(Segment::J)
                .copied()
                .unwrap()
                .allele_id
                .as_usize()
                < 2
        );
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

    // ──────────────────────────────────────────────────────────────
    // Replay validation parity
    //
    // These tests pin the replay contract: **trace supplies the
    // proposal; engine validation decides whether it applies**. A
    // recorded `AlleleId` must pass refdata + support + contract +
    // feasibility checks before `sim.with_allele_assigned` runs.
    // ──────────────────────────────────────────────────────────────

    use crate::address::ChoiceAddress;
    use crate::contract::{ChoiceContext, Contract, ContractSet, ContractViolation};
    use crate::pass::PassContext;
    use crate::refdata::{ChainType, RefDataConfig};
    use crate::replay::TraceCursor;
    use crate::rng::Rng;
    use crate::trace::{ChoiceRecord, Trace};

    /// Drive the pass through the replay path with the given
    /// pre-recorded `AlleleId` choice, optional contracts, and the
    /// supplied refdata. Returns the result + a snapshot of any
    /// trace records the pass wrote.
    fn run_replay(
        pass: &SampleAllelePass,
        sim: Simulation,
        recorded_id: u32,
        refdata: Option<&RefDataConfig>,
        contracts: Option<&ContractSet>,
    ) -> (Result<Simulation, PassError>, Trace) {
        let addr = pass.choice_address();
        let records = vec![ChoiceRecord::new(
            addr.to_string(),
            ChoiceValue::AlleleId(recorded_id),
        )];
        let mut cursor = TraceCursor::from_owned(records);
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let result;
        {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata,
                contracts,
                feasibility: None,
                reference_index: None,
                replay_cursor: Some(&mut cursor),
                event_log_sink: None,
            };
            result = pass.execute_with_sampling_mode(&sim, &mut ctx, true);
        }
        (result, trace)
    }

    fn refdata_with_v_pool(n_alleles: usize) -> RefDataConfig {
        let mut rd = RefDataConfig::empty(ChainType::Vdj);
        for i in 0..n_alleles {
            let _ = rd.v_pool.push(crate::refdata::Allele {
                name: format!("test_allele_{}*01", i),
                gene: format!("test_allele_{}", i),
                seq: vec![b'A'; 30],
                segment: Segment::V,
                anchor: Some(10),
                functional_status: None,
                subregions: Vec::new(),
            });
        }
        rd
    }

    #[test]
    fn replay_happy_path_assigns_recorded_id_when_all_checks_pass() {
        // Pool of 3 V alleles; replay recorded id 1; matching
        // refdata. Should land cleanly with the trace re-emitted.
        let pool = make_test_pool(3, Segment::V);
        let pass =
            SampleAllelePass::new(Segment::V, Box::new(AllelePoolDist::uniform(&pool)));
        let rd = refdata_with_v_pool(3);

        let (result, trace) =
            run_replay(&pass, Simulation::new(), 1, Some(&rd), None);
        let sim = result.expect("happy-path replay must succeed");

        assert_eq!(
            sim.assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .allele_id
                .index(),
            1
        );
        // Trace re-emitted with the same id.
        let rec = trace.find("sample_allele.v").unwrap();
        assert_eq!(rec.value, ChoiceValue::AlleleId(1));
    }

    #[test]
    fn replay_rejects_id_missing_from_refdata() {
        // Recorded id 5; refdata only has 3 alleles → MissingAllele.
        let pool = make_test_pool(3, Segment::V);
        let pass =
            SampleAllelePass::new(Segment::V, Box::new(AllelePoolDist::uniform(&pool)));
        let rd = refdata_with_v_pool(3);

        let (result, _) = run_replay(&pass, Simulation::new(), 5, Some(&rd), None);
        match result.unwrap_err() {
            PassError::MissingAllele {
                segment,
                allele_id,
                pass_name,
            } => {
                assert_eq!(segment, Segment::V);
                assert_eq!(allele_id, 5);
                assert_eq!(pass_name, "sample_allele.v");
            }
            other => panic!("expected MissingAllele, got {other:?}"),
        }
    }

    #[test]
    fn replay_rejects_id_outside_distribution_support() {
        // Distribution support enumerates ids 0..3. Replay recorded
        // id 7. Even if refdata happens to contain a 7th allele,
        // the distribution can't have drawn it — the recorded value
        // didn't come from this plan's natural draw.
        let pool = make_test_pool(3, Segment::V);
        let pass =
            SampleAllelePass::new(Segment::V, Box::new(AllelePoolDist::uniform(&pool)));
        let rd = refdata_with_v_pool(10); // wider refdata so MissingAllele isn't first

        let (result, _) = run_replay(&pass, Simulation::new(), 7, Some(&rd), None);
        match result.unwrap_err() {
            PassError::InvalidDistributionOutput {
                pass_name,
                address,
                value,
                reason,
            } => {
                assert_eq!(pass_name, "sample_allele.v");
                assert_eq!(address, "sample_allele.v");
                assert_eq!(value, 7);
                assert_eq!(reason, "allele_id_not_in_distribution_support");
            }
            other => panic!("expected InvalidDistributionOutput, got {other:?}"),
        }
    }

    #[test]
    fn replay_rejects_id_when_contract_short_circuits() {
        // Stub contract that rejects id 1 at `sample_allele.v`.
        // Replay records id 1 → ContractViolation surfaced (rather
        // than silent application).
        struct RejectAlleleId1;
        impl Contract for RejectAlleleId1 {
            fn name(&self) -> &str {
                "reject_allele_id_1"
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), ContractViolation> {
                Ok(())
            }
            fn admits_typed(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
                _context: ChoiceContext<'_>,
                candidate: &ChoiceValue,
            ) -> Result<(), ContractViolation> {
                if let ChoiceValue::AlleleId(1) = candidate {
                    return Err(ContractViolation::new(self.name(), "id 1 rejected"));
                }
                Ok(())
            }
        }

        let pool = make_test_pool(3, Segment::V);
        let pass =
            SampleAllelePass::new(Segment::V, Box::new(AllelePoolDist::uniform(&pool)));
        let rd = refdata_with_v_pool(3);
        let contracts = ContractSet::new().with(Box::new(RejectAlleleId1));

        let (result, _) =
            run_replay(&pass, Simulation::new(), 1, Some(&rd), Some(&contracts));
        match result.unwrap_err() {
            PassError::ContractViolation {
                pass_name,
                violations,
            } => {
                assert_eq!(pass_name, "sample_allele.v");
                assert_eq!(violations.len(), 1);
                assert_eq!(violations[0].contract_name, "reject_allele_id_1");
            }
            other => panic!("expected ContractViolation, got {other:?}"),
        }
        // Sanity: a valid id under the same contract still lands.
        let (ok_result, _) =
            run_replay(&pass, Simulation::new(), 0, Some(&rd), Some(&contracts));
        let _ = ok_result.expect("id 0 should pass the contract and apply");
        // And we know we used the choice_address — match the
        // generate_np pattern: the address is the canonical one.
        let _ = ChoiceAddress::SampleAllele(address::VdjSegment::V);
    }
}
