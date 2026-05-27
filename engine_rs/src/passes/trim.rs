//! `TrimPass` — recombination-stage trim sampling (C.6).

use crate::address;
use crate::assignment::TrimEnd;
use crate::contract::{LengthSupport, TrimTarget};
use crate::dist::{sample_filtered_with_policy, Distribution, EmptySupport};
use crate::ir::{Segment, Simulation, SimulationBuilder};
use crate::pass::{
    IntegerSupport, Pass, PassCompileFact, PassContext, PassEffect, PassError, PassRequirement,
};
use crate::trace::ChoiceValue;

/// v3.0 empty-support policy for trim sampling. A zero-length trim
/// is the only architectural no-op for a length sampler — same
/// shape as [`GenerateNPPass`](crate::passes::GenerateNPPass)'s
/// length policy.
///
/// **Permissive degradation note:** the per-candidate predicate
/// in [`TrimPass::sample_trim`] combines the contract bundle's
/// `admissible_trim_lengths` support with the runtime
/// `FeasibilityContext::admits` filter. When their intersection
/// with the natural distribution is empty, the sentinel `0` is
/// written *without re-checking* feasibility — applying the
/// sentinel skips the predicate entirely. The built-in
/// feasibility filter always admits length 0 (zero trim
/// preserves all downstream geometry), so this is safe under
/// the canonical pipeline. A custom feasibility filter that
/// specifically rejects length 0 would see its rejection
/// silently bypassed under permissive empty-support; that
/// scenario should be detected by upstream plan validation
/// (a feasibility filter incompatible with the sentinel is a
/// plan-construction issue, not a runtime fallback bug).
const TRIM_LENGTH_EMPTY_SUPPORT: EmptySupport<i64> = EmptySupport::Sentinel(0);

/// Sample a trim amount from a distribution and apply it to the
/// assigned allele on the given segment / end.
///
/// The pass is parameterized by:
/// - **Segment** — must be V, D, or J (NP and C are rejected at
///   construction).
/// - **End** — `TrimEnd::Five` or `TrimEnd::Three`. All six
///   `(segment, end)` combinations are syntactically valid; whether
///   biology uses them is up to the plan author.
/// - **Distribution** — any `Box<dyn Distribution<Output = i64>>`,
///   typically an `EmpiricalLengthDist` constructed at plan-build
///   time with the empirical trim distribution for that segment.
///
/// On execute the pass:
/// 1. Draws one `i64` trim value from the distribution via
///    `ctx.rng`.
/// 2. Validates `0 <= value <= u16::MAX` — distributions that
///    return negative or oversized trims are caller bugs and
///    panic loudly with the offending value.
/// 3. Records `ChoiceValue::Int(value)` to the trace at address
///    `"trim.{segment}_{end}"` (e.g. `"trim.v_3"`, `"trim.j_5"`).
/// 4. Returns `sim.with_trim(segment, end, value as u16)`.
///
/// **Pre-condition:** an allele must already be assigned to
/// `segment` (i.e., the matching `SampleAllelePass` ran earlier in
/// the plan). `Simulation::with_trim` panics if the slot is empty,
/// and so does this pass by extension.
pub struct TrimPass {
    segment: Segment,
    end: TrimEnd,
    distribution: Box<dyn Distribution<Output = i64>>,
}

impl TrimPass {
    /// Construct a trim pass.
    ///
    /// Panics if `segment` is anything other than V, D, or J.
    pub fn new(
        segment: Segment,
        end: TrimEnd,
        distribution: Box<dyn Distribution<Output = i64>>,
    ) -> Self {
        match segment {
            Segment::V | Segment::D | Segment::J => {}
            _ => panic!("TrimPass: segment must be V, D, or J — got {:?}", segment),
        }
        Self {
            segment,
            end,
            distribution,
        }
    }

    pub fn segment(&self) -> Segment {
        self.segment
    }

    pub fn end(&self) -> TrimEnd {
        self.end
    }

    /// The hierarchical-string address (D3) at which this pass
    /// records its choice. Same string as `name()`.
    fn address(&self) -> &'static str {
        address::trim_vdj(self.segment, self.end)
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
        address::ChoiceAddress::Trim {
            segment: self.typed_segment(),
            end: self.end,
        }
    }

    fn sample_trim(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<i64, PassError> {
        let refdata = ctx.refdata;
        let contracts = ctx.contracts;
        let feasibility = ctx.feasibility;
        let pass_index = ctx.pass_index;

        if contracts.is_some() || feasibility.is_some() {
            // v3.0 constrain-before-propose: ask the bundle once
            // for the typed [`LengthSupport`] over `[0, u16::MAX]`,
            // then filter the natural distribution to that support.
            // Anchor-preservation gives a closed-form
            // `Full(anchor)` / `Full(allele_len − anchor − 3)`
            // upper bound without per-candidate trait dispatch;
            // contracts that don't override the hook fall through
            // to the trait default `Full(requested_max)`, so the
            // composed support is the tightest opinion in the
            // bundle.
            let trim_support = contracts
                .map(|c| {
                    c.admissible_trim_lengths(
                        sim,
                        refdata,
                        TrimTarget {
                            segment: self.segment,
                            end: self.end,
                        },
                        u16::MAX as u32,
                    )
                })
                .unwrap_or(LengthSupport::Full(u16::MAX as u32));

            let pass_name = self.name().to_string();
            let address = self.address();
            let outcome = sample_filtered_with_policy(
                ctx.rng,
                self.distribution.as_ref(),
                |candidate| {
                    // Candidates outside `[0, u16::MAX]` are
                    // invalid; surface as a downstream
                    // `InvalidDistributionOutput`. Drop them from
                    // the filtered support here.
                    if *candidate < 0 || *candidate > u16::MAX as i64 {
                        return false;
                    }
                    let length = *candidate as u32;
                    let in_support = match &trim_support {
                        LengthSupport::Full(max) => length <= *max,
                        LengthSupport::Subset(set) => set.binary_search(&length).is_ok(),
                        LengthSupport::Empty => false,
                    };
                    if !in_support {
                        return false;
                    }
                    feasibility.map_or(true, |feasibility| {
                        feasibility.admits(
                            pass_index,
                            sim,
                            refdata,
                            address,
                            &ChoiceValue::Int(*candidate),
                        )
                    })
                },
                strict,
                &pass_name,
                address,
                TRIM_LENGTH_EMPTY_SUPPORT,
            )?;
            return Ok(outcome.expect("TRIM_LENGTH_EMPTY_SUPPORT is Sentinel(0), never Skip"));
        }

        Ok(self.distribution.sample(ctx.rng))
    }

    fn execute_with_validation(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        if strict && sim.assignments.get(self.segment).is_none() {
            return Err(PassError::missing_assignment(self.name(), self.segment));
        }

        // Trace-injected replay: consume the recorded trim length
        // from the cursor instead of drawing from the distribution.
        // The downstream validation + `sim.with_trim(...)` path runs
        // unchanged — bad recorded values still surface as
        // `InvalidDistributionOutput` rather than corrupting the IR.
        let value = if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            cursor
                .expect_int(self.choice_address())
                .map_err(|reason| PassError::replay(self.name(), reason))?
        } else {
            self.sample_trim(sim, ctx, strict)?
        };
        if strict && value < 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                self.address(),
                value,
                "negative_trim",
            ));
        }
        if strict && value > u16::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                self.address(),
                value,
                "trim_exceeds_u16",
            ));
        }

        assert!(
            value >= 0,
            "TrimPass({}): distribution returned negative trim {}",
            self.address(),
            value
        );
        assert!(
            value <= u16::MAX as i64,
            "TrimPass({}): distribution returned trim {} > u16::MAX",
            self.address(),
            value
        );

        ctx.trace
            .record_choice(self.choice_address(), ChoiceValue::Int(value));
        // Route the trim update through the builder so the
        // `TrimChanged` event flows onto every attached sink. The
        // persistent mutation is still `Simulation::with_trim`
        // (called inside `SimulationBuilder::update_trim`).
        // Forward the captured event into `ctx.event_log_sink` when
        // the caller supplied one — the per-pass `EventRecord`
        // then carries the consequence on `simulation_events`.
        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }
        builder.update_trim(self.segment, self.end, value as u16);
        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(builder.seal_event_log_observer());
        }
        Ok(builder.seal())
    }
}

impl Pass for TrimPass {
    fn name(&self) -> &str {
        self.address()
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_validation(sim, ctx, false)
            .expect("TrimPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_validation(sim, ctx, true)
    }

    fn declared_choice_patterns(&self) -> Vec<crate::address::ChoiceAddressPattern> {
        vec![crate::address::ChoiceAddressPattern::Trim {
            segment: self.typed_segment(),
            end: self.end,
        }]
    }

    fn requirements(&self) -> Vec<PassRequirement> {
        vec![PassRequirement::AlleleAssignment(self.segment)]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::TrimAllele(self.segment)]
    }

    fn compile_facts(&self) -> Vec<PassCompileFact> {
        vec![PassCompileFact::TrimSupport {
            segment: self.segment,
            end: self.end,
            support: IntegerSupport::from_weighted_pairs(self.distribution.support()),
        }]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::{AllelePoolDist, EmpiricalLengthDist, FilteredSampleError};
    use crate::pass::testing::PassRuntime;
    use crate::pass::{PassError, PassPlan};
    use crate::passes::sample_allele::test_support::make_test_pool;
    use crate::passes::SampleAllelePass;

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn trim_pass_rejects_np1() {
        let _ = TrimPass::new(
            Segment::Np1,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        );
    }

    #[test]
    fn trim_pass_addresses_cover_all_six_combinations() {
        let dist = || Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)]));

        let cases = [
            (Segment::V, TrimEnd::Five, "trim.v_5"),
            (Segment::V, TrimEnd::Three, "trim.v_3"),
            (Segment::D, TrimEnd::Five, "trim.d_5"),
            (Segment::D, TrimEnd::Three, "trim.d_3"),
            (Segment::J, TrimEnd::Five, "trim.j_5"),
            (Segment::J, TrimEnd::Three, "trim.j_3"),
        ];
        for (seg, end, expected_addr) in cases {
            let pass = TrimPass::new(seg, end, dist());
            assert_eq!(pass.name(), expected_addr);
            assert_eq!(pass.declared_choices(), vec![expected_addr.to_string()]);
        }
    }

    #[test]
    fn trim_pass_applies_trim_and_records_choice() {
        // V allele assigned, then trim_3 sampled.
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            // Single-value dist: always returns 4.
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        // Trim recorded at canonical address.
        let rec = outcome.trace.find("trim.v_3").expect("trim recorded");
        assert_eq!(rec.value, ChoiceValue::Int(4));

        // Trim applied to the V allele instance.
        let v = outcome
            .final_simulation()
            .assignments
            .get(Segment::V)
            .copied()
            .unwrap();
        assert_eq!(v.trim_3, 4);
        assert_eq!(v.trim_5, 0);
    }

    #[test]
    #[should_panic(expected = "no instance assigned to segment")]
    fn trim_pass_panics_when_no_allele_assigned() {
        // The trim pass relies on `Simulation::with_trim`, which
        // panics when no instance is assigned to the segment. This
        // test pins the propagation.
        let mut plan = PassPlan::new();
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        )));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 0);
    }

    #[test]
    fn trim_pass_strict_errors_when_no_allele_assigned() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        )));

        let err = PassRuntime::execute_strict_with_context(&plan, Simulation::new(), 0, None, None)
            .unwrap_err();

        assert_eq!(err.pass_name(), "trim.v_3");
        assert!(matches!(
            err,
            PassError::MissingAssignment {
                segment: Segment::V,
                ..
            }
        ));
    }

    #[test]
    #[should_panic(expected = "negative trim")]
    fn trim_pass_panics_on_negative_distribution_output() {
        // Runtime defensive check. UniformInt::new(-5, -1) always
        // returns a negative trim, which TrimPass must reject loudly.
        use crate::dist::UniformInt;
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(UniformInt::new(-5, -1)),
        )));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 1);
    }

    #[test]
    fn trim_pass_strict_errors_on_negative_distribution_output() {
        use crate::dist::UniformInt;
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(UniformInt::new(-5, -1)),
        )));

        let err = PassRuntime::execute_strict_with_context(&plan, Simulation::new(), 1, None, None)
            .unwrap_err();

        assert_eq!(err.pass_name(), "trim.v_3");
        assert_eq!(err.address(), "trim.v_3");
        assert!(matches!(
            err,
            PassError::InvalidDistributionOutput { value, .. } if value < 0
        ));
    }

    #[test]
    fn trim_pass_full_vdj_chain_records_six_choices() {
        // Heavy chain shape: sample V/D/J, then trim each side.
        let v_pool = make_test_pool(1, Segment::V);
        let d_pool = make_test_pool(1, Segment::D);
        let j_pool = make_test_pool(1, Segment::J);
        let dist = || Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)]));

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
        plan.push(Box::new(TrimPass::new(Segment::V, TrimEnd::Three, dist())));
        plan.push(Box::new(TrimPass::new(Segment::D, TrimEnd::Five, dist())));
        plan.push(Box::new(TrimPass::new(Segment::D, TrimEnd::Three, dist())));
        plan.push(Box::new(TrimPass::new(Segment::J, TrimEnd::Five, dist())));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 42);

        // 3 sampling + 4 trim = 7 trace entries.
        assert_eq!(outcome.trace.len(), 7);

        // All four trim addresses present with value 2.
        for addr in ["trim.v_3", "trim.d_5", "trim.d_3", "trim.j_5"] {
            assert_eq!(outcome.trace.find(addr).unwrap().value, ChoiceValue::Int(2));
        }

        // Trims applied to the assignments.
        let sim = outcome.final_simulation();
        assert_eq!(sim.assignments.get(Segment::V).copied().unwrap().trim_3, 2);
        assert_eq!(sim.assignments.get(Segment::D).copied().unwrap().trim_5, 2);
        assert_eq!(sim.assignments.get(Segment::D).copied().unwrap().trim_3, 2);
        assert_eq!(sim.assignments.get(Segment::J).copied().unwrap().trim_5, 2);
    }

    #[test]
    fn trim_pass_is_deterministic_under_same_seed() {
        let v_pool = make_test_pool(1, Segment::V);
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(SampleAllelePass::new(
                Segment::V,
                Box::new(AllelePoolDist::uniform(&v_pool)),
            )));
            p.push(Box::new(TrimPass::new(
                Segment::V,
                TrimEnd::Three,
                Box::new(EmpiricalLengthDist::from_pairs(vec![
                    (0, 1.0),
                    (1, 2.0),
                    (2, 3.0),
                    (5, 1.0),
                ])),
            )));
            p
        };

        let oa = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);
        assert_eq!(
            oa.final_simulation()
                .assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .trim_3,
            ob.final_simulation()
                .assignments
                .get(Segment::V)
                .copied()
                .unwrap()
                .trim_3
        );
        assert_eq!(oa.trace.choices(), ob.trace.choices());
    }

    // v3.0 Phase E: trim lengths sampled from natural ∩ admissible.
    // Pin that the bundle's `admissible_trim_lengths` is honored:
    // under productive() the V's 5' trim must stay ≤ anchor offset,
    // even when the natural distribution allows larger values.
    mod productive_trim {
        use super::*;
        use crate::assignment::AlleleInstance;
        use crate::contract::productive;
        use crate::refdata::{Allele, AlleleId, AllelePool, ChainType, RefDataConfig};

        fn vj_cfg_with_v_anchor(v_anchor: u16) -> RefDataConfig {
            let mut cfg = RefDataConfig::empty(ChainType::Vj);
            // V allele = 9 bytes, anchor at `v_anchor`.
            let mut seq = Vec::new();
            for _ in 0..v_anchor as usize {
                seq.push(b'A');
            }
            seq.extend_from_slice(b"TGT"); // anchor codon Cys
            for _ in (v_anchor as usize + 3)..9 {
                seq.push(b'C');
            }
            assert_eq!(seq.len(), 9);
            let _ = cfg.v_pool.push(Allele {
                name: "v_trim*01".into(),
                gene: "v_trim".into(),
                seq,
                segment: Segment::V,
                anchor: Some(v_anchor),
            });
            let _ = cfg.j_pool.push(Allele {
                name: "j_trim*01".into(),
                gene: "j_trim".into(),
                seq: b"TGGAAA".to_vec(),
                segment: Segment::J,
                anchor: Some(0),
            });
            cfg
        }

        fn make_v_assigned_sim() -> (RefDataConfig, Simulation) {
            let cfg = vj_cfg_with_v_anchor(4);
            // V anchor at 4 → max admissible 5' trim = 4.
            let sim = Simulation::new()
                .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
                .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
            (cfg, sim)
        }

        #[test]
        fn five_prime_trim_under_productive_never_exceeds_anchor_offset() {
            // Natural distribution puts mass on lengths 0..=8.
            // Under productive() the support must clamp to 0..=4
            // (anchor at offset 4). Verified empirically across
            // many seeds.
            let (cfg, sim) = make_v_assigned_sim();
            let contracts = productive();
            for seed in 0..256u64 {
                let mut plan = PassPlan::new();
                plan.push(Box::new(TrimPass::new(
                    Segment::V,
                    TrimEnd::Five,
                    Box::new(EmpiricalLengthDist::from_pairs(vec![
                        (0, 1.0),
                        (1, 1.0),
                        (2, 1.0),
                        (3, 1.0),
                        (4, 1.0),
                        (5, 1.0),
                        (6, 1.0),
                        (7, 1.0),
                        (8, 1.0),
                    ])),
                )));
                let outcome = PassRuntime::execute_with_context(
                    &plan,
                    sim.clone(),
                    seed,
                    Some(&cfg),
                    Some(&contracts),
                );
                let trim = outcome
                    .final_simulation()
                    .assignments
                    .get(Segment::V)
                    .copied()
                    .unwrap()
                    .trim_5;
                assert!(
                    trim <= 4,
                    "seed {} produced trim_5 = {} > anchor offset 4",
                    seed,
                    trim
                );
                assert!(contracts
                    .verify(outcome.final_simulation(), Some(&cfg))
                    .is_ok());
            }
        }

        #[test]
        fn five_prime_trim_strict_errors_when_support_disjoint_from_natural() {
            // Natural distribution puts mass ONLY on out-of-support
            // values (5..=8 under V anchor at 4). Strict mode must
            // surface `EmptyAdmissibleSupport`; the contract's
            // typed support has no overlap with the natural support.
            let (cfg, sim) = make_v_assigned_sim();
            let contracts = productive();
            let mut plan = PassPlan::new();
            plan.push(Box::new(TrimPass::new(
                Segment::V,
                TrimEnd::Five,
                Box::new(EmpiricalLengthDist::from_pairs(vec![
                    (5, 1.0),
                    (6, 1.0),
                    (7, 1.0),
                    (8, 1.0),
                ])),
            )));
            let err = PassRuntime::execute_strict_with_context(
                &plan,
                sim,
                0,
                Some(&cfg),
                Some(&contracts),
            )
            .unwrap_err();
            assert_eq!(err.pass_name(), "trim.v_5");
            assert_eq!(err.address(), "trim.v_5");
            assert_eq!(
                err.constraint_reason(),
                Some(FilteredSampleError::EmptyAdmissibleSupport)
            );
        }

        #[test]
        fn five_prime_trim_permissive_noops_when_support_disjoint_from_natural() {
            // Same disjoint support as the strict test above, but
            // permissive mode must not fall back to an unconstrained
            // trim. The pass records/apply a zero-length no-op.
            let (cfg, sim) = make_v_assigned_sim();
            let contracts = productive();
            let mut plan = PassPlan::new();
            plan.push(Box::new(TrimPass::new(
                Segment::V,
                TrimEnd::Five,
                Box::new(EmpiricalLengthDist::from_pairs(vec![
                    (5, 1.0),
                    (6, 1.0),
                    (7, 1.0),
                    (8, 1.0),
                ])),
            )));

            let outcome =
                PassRuntime::execute_with_context(&plan, sim, 0, Some(&cfg), Some(&contracts));

            assert_eq!(
                outcome.trace.find("trim.v_5").unwrap().value,
                ChoiceValue::Int(0)
            );
            assert_eq!(
                outcome
                    .final_simulation()
                    .assignments
                    .get(Segment::V)
                    .copied()
                    .unwrap()
                    .trim_5,
                0
            );
            assert!(contracts
                .verify(outcome.final_simulation(), Some(&cfg))
                .is_ok());
        }

        #[test]
        fn no_contracts_path_unaffected_by_admissible_trim_lengths() {
            // Without contracts, the constrain-before-propose
            // helper is bypassed entirely — the raw natural
            // distribution governs the draw.
            let cfg = vj_cfg_with_v_anchor(4);
            let _ = cfg;
            let sim = Simulation::new()
                .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));
            let mut plan = PassPlan::new();
            plan.push(Box::new(TrimPass::new(
                Segment::V,
                TrimEnd::Five,
                Box::new(EmpiricalLengthDist::from_pairs(vec![(7, 1.0)])),
            )));
            let outcome = PassRuntime::execute(&plan, sim, 0);
            assert_eq!(
                outcome
                    .final_simulation()
                    .assignments
                    .get(Segment::V)
                    .copied()
                    .unwrap()
                    .trim_5,
                7,
            );
        }

        // Suppress unused-import warning on `AllelePool` — the
        // `make_test_pool` helper isn't used here but the test
        // module re-uses the parent's imports.
        #[allow(dead_code)]
        const _: Option<AllelePool> = None;
    }
}
