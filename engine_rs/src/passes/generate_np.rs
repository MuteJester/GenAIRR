//! `GenerateNPPass` — TdT-like N-nucleotide region generation (C.7).

use crate::address;
use crate::dist::Distribution;
use crate::ir::{Segment, Simulation};
use crate::pass::{IntegerSupport, Pass, PassCompileFact, PassContext, PassEffect, PassError};

// `execute_with_sampling_mode` and the constraint-aware sampling
// helpers live in submodules so the codon-rail + admit-mask observer
// wiring (/3) and the `sample_base_with_admit_mask` bypass
// stay focused. Both submodules `impl GenerateNPPass` here.
mod execution;
mod sampling;

/// Generate one NP (non-template) region — a stretch of bases
/// interpolated between adjacent V/D/J segments by the TdT enzyme
/// during recombination.
///
/// Single pass produces single IR revision (D2 boundary), but
/// internally records *each* atomic random choice to the trace at
/// its own indexed address. Length, then N bases. The trace stays
/// faithful at the choice level even though the IR revision count
/// stays at the biological-event level.
///
/// Parameterized by:
/// - **np_segment** — `Segment::Np1` or `Segment::Np2`. Other
///   segments rejected at construction.
/// - **length_dist** — `Box<dyn Distribution<Output = i64>>`.
///   Empirical NP length distribution. Negative or oversized
///   lengths are caller bugs and panic at execute time.
/// - **base_dist** — `Box<dyn Distribution<Output = u8>>`.
///   Per-base distribution (typically `UniformBase` when no
///   empirical TdT model is configured). When contracts are active
///   and the distribution exposes finite support, each base draw
///   is filtered through `ContractSet::admits` before being
///   recorded.
///
/// On execute:
/// 1. Sample length L from `length_dist`. Validate `0 <= L <= u32::MAX`.
/// 2. Record `ChoiceValue::Int(L)` at `"np.{np1|np2}.length"`.
/// 3. For each of L positions: sample a base, record at
///    `"np.{np1|np2}.bases[i]"`, push a synthetic nucleotide
///    (`Nucleotide::synthetic(base, np_segment, N_NUC)`) onto the
///    pool.
/// 4. Construct `Region(np_segment, [start, start+L))` with
///    `frame_phase = 0` (cross-region frame chaining is C.8's
///    responsibility). No codon-rail data is stored on the region;
///    callers that need translation call
///    [`crate::ir::compute_codon_rail`] on demand.
/// 5. Append the region to the sequence.
pub struct GenerateNPPass {
    np_segment: Segment,
    length_dist: Box<dyn Distribution<Output = i64>>,
    base_dist: Box<dyn Distribution<Output = u8>>,
}

impl GenerateNPPass {
    /// Construct an NP-generation pass.
    ///
    /// Panics if `np_segment` is not `Segment::Np1` or
    /// `Segment::Np2`.
    pub fn new(
        np_segment: Segment,
        length_dist: Box<dyn Distribution<Output = i64>>,
        base_dist: Box<dyn Distribution<Output = u8>>,
    ) -> Self {
        match np_segment {
            Segment::Np1 | Segment::Np2 => {}
            _ => panic!(
                "GenerateNPPass: np_segment must be Np1 or Np2 — got {:?}",
                np_segment
            ),
        }
        Self {
            np_segment,
            length_dist,
            base_dist,
        }
    }

    pub fn np_segment(&self) -> Segment {
        self.np_segment
    }

    pub(super) fn length_address(&self) -> &'static str {
        address::np_length_region(self.np_segment)
    }

    pub(super) fn typed_np_segment(&self) -> address::NpSegment {
        match self.np_segment {
            Segment::Np1 => address::NpSegment::Np1,
            Segment::Np2 => address::NpSegment::Np2,
            Segment::V | Segment::D | Segment::J => {
                unreachable!("constructor rejects V/D/J segments")
            }
        }
    }

    pub(super) fn length_choice_address(&self) -> address::ChoiceAddress {
        address::ChoiceAddress::NpLength(self.typed_np_segment())
    }

    pub(super) fn pass_name(&self) -> &'static str {
        address::generate_np_region(self.np_segment)
    }
}

impl Pass for GenerateNPPass {
    fn name(&self) -> &str {
        self.pass_name()
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("GenerateNPPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        let segment = self.typed_np_segment();
        vec![
            address::ChoiceAddressPattern::NpLength(segment),
            address::ChoiceAddressPattern::NpBase(segment),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::AppendRegion(self.np_segment)]
    }

    fn compile_facts(&self) -> Vec<PassCompileFact> {
        vec![PassCompileFact::NpLengthSupport {
            segment: self.np_segment,
            support: IntegerSupport::from_weighted_pairs(self.length_dist.support()),
        }]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::{EmpiricalLengthDist, UniformBase};
    use crate::ir::{flag, NucHandle};
    use crate::pass::testing::PassRuntime;
    use crate::pass::{PassError, PassPlan};
    use crate::passes::EchoPass;
    use crate::trace::ChoiceValue;

    #[test]
    #[should_panic(expected = "np_segment must be Np1 or Np2")]
    fn generate_np_pass_rejects_v_segment() {
        let _ = GenerateNPPass::new(
            Segment::V,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
    }

    #[test]
    #[should_panic(expected = "np_segment must be Np1 or Np2")]
    fn generate_np_pass_rejects_d_segment() {
        let _ = GenerateNPPass::new(
            Segment::D,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
    }

    #[test]
    fn generate_np_pass_zero_length_creates_empty_region() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 0);

        let sim = outcome.final_simulation();
        // Pool is unchanged.
        assert_eq!(sim.pool.len(), 0);
        // One region was added (NP1, empty).
        assert_eq!(sim.sequence.region_count(), 1);
        let r = &sim.sequence.regions[0];
        assert_eq!(r.segment, Segment::Np1);
        assert!(r.is_empty());
        // Trace recorded length, no bases.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(
            outcome.trace.find("np.np1.length").unwrap().value,
            ChoiceValue::Int(0)
        );
    }

    #[test]
    fn generate_np_pass_pushes_n_bases_with_correct_metadata() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 1234);

        let sim = outcome.final_simulation();
        assert_eq!(sim.pool.len(), 5);
        assert_eq!(sim.sequence.region_count(), 1);

        // Region covers all 5 nucleotides.
        let r = &sim.sequence.regions[0];
        assert_eq!(r.segment, Segment::Np1);
        assert_eq!(r.start, NucHandle::new(0));
        assert_eq!(r.end, NucHandle::new(5));

        // Each pushed nucleotide is a synthetic NP1 N-nuc.
        for i in 0..5 {
            let n = sim.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.segment, Segment::Np1);
            assert!(n.germline_pos.is_none());
            assert!(n.flags.contains(flag::N_NUC));
            assert!(matches!(n.base, b'A' | b'C' | b'G' | b'T'));
        }
    }

    #[test]
    fn generate_np_pass_records_length_and_indexed_bases() {
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np2,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 7);

        // Length entry + 3 base entries = 4 trace records.
        assert_eq!(outcome.trace.len(), 4);
        assert_eq!(
            outcome.trace.find("np.np2.length").unwrap().value,
            ChoiceValue::Int(3)
        );
        for i in 0..3 {
            let addr = format!("np.np2.bases[{}]", i);
            let rec = outcome
                .trace
                .find(&addr)
                .unwrap_or_else(|| panic!("missing {}", addr));
            assert!(matches!(rec.value, ChoiceValue::Base(_)));
        }
    }

    #[test]
    fn generate_np_pass_trace_value_matches_pushed_base_at_each_position() {
        // Faithfulness: trace[i] equals pool[i].base for each NP base.
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(8, 1.0)])),
            Box::new(UniformBase),
        )));
        let outcome = PassRuntime::execute(&plan, Simulation::new(), 99);
        let sim = outcome.final_simulation();

        for i in 0..8 {
            let addr = format!("np.np1.bases[{}]", i);
            let recorded = match outcome.trace.find(&addr).unwrap().value {
                ChoiceValue::Base(b) => b,
                _ => panic!("wrong variant"),
            };
            let pushed = sim.pool.get(NucHandle::new(i)).unwrap().base;
            assert_eq!(
                recorded, pushed,
                "trace at {} = {} but pool[{}] = {}",
                addr, recorded as char, i, pushed as char
            );
        }
    }

    #[test]
    fn generate_np_pass_is_deterministic_under_same_seed() {
        let plan = || {
            let mut p = PassPlan::new();
            p.push(Box::new(GenerateNPPass::new(
                Segment::Np1,
                Box::new(EmpiricalLengthDist::from_pairs(vec![
                    (1, 1.0),
                    (3, 2.0),
                    (5, 1.0),
                ])),
                Box::new(UniformBase),
            )));
            p
        };
        let oa = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);
        let ob = PassRuntime::execute(&plan(), Simulation::new(), 0xc0ff_ee);

        assert_eq!(
            oa.final_simulation().pool.len(),
            ob.final_simulation().pool.len()
        );
        assert_eq!(oa.trace.choices(), ob.trace.choices());
    }

    #[test]
    fn generate_np_pass_declared_choices_lists_length_and_bases_pattern() {
        let pass = GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
            Box::new(UniformBase),
        );
        let declared = pass.declared_choices();
        assert_eq!(declared.len(), 2);
        assert!(declared.contains(&"np.np1.length".to_string()));
        assert!(declared.contains(&"np.np1.bases[0..n]".to_string()));
    }

    #[test]
    fn generate_np_pass_appends_region_after_existing_pool_state() {
        // Pre-populate the pool with some nucleotides, then run NP
        // generation. The new region should start at the current
        // pool boundary, not at handle 0.
        let mut plan = PassPlan::new();
        // Three Echo passes push 3 nucleotides as germline.
        plan.push(Box::new(EchoPass::new(b'A', 0, Segment::V)));
        plan.push(Box::new(EchoPass::new(b'C', 1, Segment::V)));
        plan.push(Box::new(EchoPass::new(b'G', 2, Segment::V)));
        // Then NP1 of length 4.
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
            Box::new(UniformBase),
        )));

        let outcome = PassRuntime::execute(&plan, Simulation::new(), 5);
        let sim = outcome.final_simulation();

        // Total pool: 3 V + 4 NP1 = 7.
        assert_eq!(sim.pool.len(), 7);
        // One region (NP1) — Echo doesn't create regions.
        assert_eq!(sim.sequence.region_count(), 1);
        let r = &sim.sequence.regions[0];
        assert_eq!(r.start, NucHandle::new(3));
        assert_eq!(r.end, NucHandle::new(7));
        assert_eq!(r.segment, Segment::Np1);
    }

    #[test]
    #[should_panic(expected = "negative value")]
    fn generate_np_pass_panics_on_negative_length() {
        use crate::dist::UniformInt;
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(UniformInt::new(-3, -1)),
            Box::new(UniformBase),
        )));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 0);
    }

    #[test]
    fn generate_np_pass_strict_errors_on_negative_length() {
        use crate::dist::UniformInt;
        let mut plan = PassPlan::new();
        plan.push(Box::new(GenerateNPPass::new(
            Segment::Np1,
            Box::new(UniformInt::new(-3, -1)),
            Box::new(UniformBase),
        )));

        let err = PassRuntime::execute_strict_with_context(&plan, Simulation::new(), 0, None, None)
            .unwrap_err();

        assert_eq!(err.pass_name(), "generate_np.np1");
        assert_eq!(err.address(), "np.np1.length");
        assert!(matches!(
            err,
            PassError::InvalidDistributionOutput { value, .. } if value < 0
        ));
    }

    // ── Trace-injected replay (Tier 3 NP base) ──────────────────

    use crate::address::{ChoiceAddress, NpSegment};
    use crate::assignment::AlleleInstance;
    use crate::contract::{productive, ProductiveJunctionFrame};
    use crate::contract::ContractSet;
    use crate::ir::{Nucleotide, Region, Segment};
    use crate::pass::PassContext;
    use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};
    use crate::replay::{ReplayError, TraceCursor};
    use crate::rng::Rng;
    use crate::trace::{ChoiceRecord, Trace};

    fn rec(addr: ChoiceAddress, v: ChoiceValue) -> ChoiceRecord {
        ChoiceRecord::new(addr.to_string(), v)
    }

    fn make_np_pass(np_length: i64) -> GenerateNPPass {
        GenerateNPPass::new(
            Segment::Np1,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(np_length, 1.0)])),
            Box::new(UniformBase),
        )
    }

    fn run_np_replay(
        pass: &GenerateNPPass,
        sim: Simulation,
        records: Vec<ChoiceRecord>,
        contracts: Option<&ContractSet>,
        refdata: Option<&RefDataConfig>,
    ) -> (Result<Simulation, PassError>, Trace, u64) {
        let mut cursor = TraceCursor::from_owned(records);
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let result = {
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
            pass.execute_checked(&sim, &mut ctx)
        };
        let words = rng.words_consumed();
        (result, trace, words)
    }

    // ── 1. Valid canonical recorded NP base applies with zero RNG ──

    #[test]
    fn np_base_replay_consumes_canonical_recorded_base_with_zero_rng() {
        let pass = make_np_pass(3);
        let bases = b"ACG";
        let mut records = vec![rec(
            ChoiceAddress::NpLength(NpSegment::Np1),
            ChoiceValue::Int(3),
        )];
        for (i, b) in bases.iter().enumerate() {
            records.push(rec(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np1,
                    index: i as u32,
                },
                ChoiceValue::Base(*b),
            ));
        }

        let (result, trace, rng_words) =
            run_np_replay(&pass, Simulation::new(), records, None, None);
        let next = result.unwrap();

        // Three NP nucleotides pushed with the recorded bytes.
        assert_eq!(next.pool.len(), 3);
        for (i, b) in bases.iter().enumerate() {
            assert_eq!(next.pool.get(NucHandle::new(i as u32)).unwrap().base, *b);
            let addr = format!("np.np1.bases[{}]", i);
            assert_eq!(trace.find(&addr).unwrap().value, ChoiceValue::Base(*b));
        }
        assert_eq!(rng_words, 0);
    }

    // ── 2. Wrong address / wrong kind → PassError::Replay ──

    #[test]
    fn np_base_replay_wrong_kind_surfaces_replay_error() {
        let pass = make_np_pass(1);
        let records = vec![
            rec(ChoiceAddress::NpLength(NpSegment::Np1), ChoiceValue::Int(1)),
            // Wrong kind: Int instead of Base.
            rec(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np1,
                    index: 0,
                },
                ChoiceValue::Int(0),
            ),
        ];
        let (result, _, _) = run_np_replay(&pass, Simulation::new(), records, None, None);
        match result.unwrap_err() {
            PassError::Replay { reason, .. } => {
                assert!(matches!(reason, ReplayError::ValueKindMismatch { .. }));
            }
            other => panic!("expected Replay::ValueKindMismatch, got {other:?}"),
        }
    }

    #[test]
    fn np_base_replay_wrong_address_surfaces_replay_error() {
        let pass = make_np_pass(1);
        let records = vec![
            rec(ChoiceAddress::NpLength(NpSegment::Np1), ChoiceValue::Int(1)),
            // Wrong address — second NP record points at index 1
            // when the pass expects index 0.
            rec(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np1,
                    index: 1,
                },
                ChoiceValue::Base(b'A'),
            ),
        ];
        let (result, _, _) = run_np_replay(&pass, Simulation::new(), records, None, None);
        match result.unwrap_err() {
            PassError::Replay { reason, .. } => {
                assert!(matches!(reason, ReplayError::AddressMismatch { .. }));
            }
            other => panic!("expected Replay::AddressMismatch, got {other:?}"),
        }
    }

    // ── 3. Contract-rejected canonical base errors instead of force-apply ──

    /// Build a VJ fixture set up exactly so the NP1 codon-frame
    /// filter under `ProductiveJunctionFrame` rejects a specific
    /// candidate. We pre-assemble V[0..3] = "TGT" (Cys) and assign
    /// J's anchor at index 0 → frame contribution from V[anchor..end]
    /// is 3 bp, J[0..anchor] is 0 bp, so total framing requires
    /// `NP1 % 3 == 0`. With NP1 length = 1, a recorded NP base of
    /// any value violates frame… but `ProductiveJunctionFrame`
    /// filters NP **length**, not NP base, so this test instead
    /// uses `NoStopCodonInJunction`: a recorded NP base that
    /// completes a stop codon must be rejected.
    fn np_stop_filter_fixture() -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        // V anchor codon TGT (Cys) at pool [0..3). Anchor offset 0.
        let _ = cfg.v_pool.push(Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: b"TGT".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        // J anchor codon TGG (Trp). Anchor at pool offset 0.
        // With NP1 length = 1 inserted between V and J, framing
        // breaks anyway — we only need J to define the junction.
        let _ = cfg.j_pool.push(Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TGG".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        // Push V[0..3] = TGT.
        for (i, b) in b"TGT".iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
                *b,
                i as u16,
                Segment::V,
            ));
            sim = next;
        }
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(3),
        ));
        sim = sim
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
        (cfg, sim)
    }

    #[test]
    fn np_base_replay_contract_rejects_inadmissible_base() {
        // Under productive(), ProductiveJunctionFrame restricts the
        // NP **length** to 3 (so the junction length stays
        // divisible by 3). With length=3 the per-NP-base filter
        // from NoStopCodonInJunction can still reject specific
        // bases that complete a stop codon. We construct: V anchor
        // codon TGT (Cys), then NP1[0..3] under draw. The junction
        // is V_anchor[0..3]=TGT + NP1[0..3] + J_anchor=TGG (3 bp).
        // Total junction length 9 — divisible by 3.
        //
        // The junction codon at NP1[0..3] starts at NP1 offset 0
        // (junction codon position 1 = NP1 bases 0..3 form codon 2
        // of the junction). The first NP1 base sits at junction-
        // codon offset 0 (start of codon 2). For a stop codon to
        // be formed at codon 2, NP1[0..3] needs to be TAA / TAG /
        // TGA.
        //
        // Simpler: use a fixture that makes NP1[0]=T force a stop
        // when NP1[1..3] are committed later. The cleanest direct
        // test of the contract reject is to record an NP base the
        // contract bundle would reject. Use TAA: NP1 length 3 and
        // recorded base = T at index 0.
        //
        // But contracts only reject AT THE FINAL CODON, which means
        // the per-NP-base check rejects only when the partial codon
        // already completed. With NP1 length 3 starting at pool
        // index 3 (after V[0..3]=TGT), the bases form codon
        // [V[anchor]=T, NP1[0]=?, NP1[1]=?]. Wait — V anchor codon
        // is positions 0..3, then NP1[0..3], then J. Junction codons:
        //   codon 0 = V[0..3] = TGT
        //   codon 1 = NP1[0..3]
        //   codon 2 = J[0..3] = TGG
        // So NP1[0..3] is one full junction codon. If NP1 ends up
        // TAA, the contract bundle rejects.
        //
        // For NP base replay to test the contract rejection, we
        // need a fixture where the contract evaluates a SPECIFIC
        // recorded NP base. Easiest: use a custom contract that
        // rejects NP1[0] = T.
        struct RejectNpFirstBaseT;
        impl crate::contract::Contract for RejectNpFirstBaseT {
            fn name(&self) -> &str {
                "reject_np_first_base_t"
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), crate::contract::ContractViolation> {
                Ok(())
            }
            fn admits_typed(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
                context: crate::contract::ChoiceContext<'_>,
                candidate: &ChoiceValue,
            ) -> Result<(), crate::contract::ContractViolation> {
                if let Some(ChoiceAddress::NpBase { index: 0, .. }) = context.address {
                    if let ChoiceValue::Base(b'T') = candidate {
                        return Err(crate::contract::ContractViolation::new(
                            self.name(),
                            "rejected",
                        ));
                    }
                }
                Ok(())
            }
        }

        let pass = make_np_pass(1);
        let contracts = ContractSet::new().with(Box::new(RejectNpFirstBaseT));
        let records = vec![
            rec(ChoiceAddress::NpLength(NpSegment::Np1), ChoiceValue::Int(1)),
            // T at index 0 — the custom contract rejects.
            rec(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np1,
                    index: 0,
                },
                ChoiceValue::Base(b'T'),
            ),
        ];
        let (result, _, _) =
            run_np_replay(&pass, Simulation::new(), records, Some(&contracts), None);
        match result.unwrap_err() {
            PassError::ConstraintSampling { address, .. } => {
                assert_eq!(address, "np.np1.bases[0]");
            }
            other => panic!("expected ConstraintSampling, got {other:?}"),
        }
    }

    // ── 4. N sentinel accepted on permissive empty-support ──

    #[test]
    fn np_base_replay_accepts_n_sentinel_when_contract_rejects_all_bases() {
        // Build a custom contract that rejects EVERY canonical NP1
        // base. That's exactly the empty-support condition the
        // fresh-RNG slow path emits the `N` sentinel for. Replay
        // should accept the recorded `N` only under this condition.
        struct RejectEveryNpBase;
        impl crate::contract::Contract for RejectEveryNpBase {
            fn name(&self) -> &str {
                "reject_every_np_base"
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), crate::contract::ContractViolation> {
                Ok(())
            }
            fn admits_typed(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
                context: crate::contract::ChoiceContext<'_>,
                _candidate: &ChoiceValue,
            ) -> Result<(), crate::contract::ContractViolation> {
                if let Some(ChoiceAddress::NpBase { .. }) = context.address {
                    return Err(crate::contract::ContractViolation::new(
                        self.name(),
                        "rejected",
                    ));
                }
                Ok(())
            }
        }

        let pass = make_np_pass(1);
        let contracts = ContractSet::new().with(Box::new(RejectEveryNpBase));
        let records = vec![
            rec(ChoiceAddress::NpLength(NpSegment::Np1), ChoiceValue::Int(1)),
            // Recorded as the policy sentinel `N` — valid because
            // contracts reject every canonical candidate.
            rec(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np1,
                    index: 0,
                },
                ChoiceValue::Base(b'N'),
            ),
        ];
        let (result, trace, _) =
            run_np_replay(&pass, Simulation::new(), records, Some(&contracts), None);
        let next = result.unwrap();
        assert_eq!(next.pool.len(), 1);
        assert_eq!(next.pool.get(NucHandle::new(0)).unwrap().base, b'N');
        assert_eq!(
            trace.find("np.np1.bases[0]").unwrap().value,
            ChoiceValue::Base(b'N'),
        );
    }

    // ── 5. N rejected when support is non-empty ──

    #[test]
    fn np_base_replay_rejects_n_when_canonical_support_is_non_empty() {
        // No contracts → base_dist.support() = [A, C, G, T] all
        // available. The fresh sampler would NEVER have emitted N
        // here. Replay must reject.
        let pass = make_np_pass(1);
        let records = vec![
            rec(ChoiceAddress::NpLength(NpSegment::Np1), ChoiceValue::Int(1)),
            rec(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np1,
                    index: 0,
                },
                ChoiceValue::Base(b'N'),
            ),
        ];
        let (result, _, _) = run_np_replay(&pass, Simulation::new(), records, None, None);
        match result.unwrap_err() {
            PassError::ConstraintSampling { address, .. } => {
                assert_eq!(address, "np.np1.bases[0]");
            }
            other => panic!("expected ConstraintSampling, got {other:?}"),
        }
    }

    #[test]
    fn np_base_replay_rejects_n_when_contract_admits_some_base() {
        // Contracts are active but still admit some canonical base.
        // The sampler would have picked one of the admitted; the
        // policy sentinel branch never fires. Recording N is a
        // mismatch with the live sampler's behavior.
        let pass = make_np_pass(1);
        let contracts = ContractSet::new().with(Box::new(ProductiveJunctionFrame::new()));
        let records = vec![
            rec(ChoiceAddress::NpLength(NpSegment::Np1), ChoiceValue::Int(1)),
            rec(
                ChoiceAddress::NpBase {
                    segment: NpSegment::Np1,
                    index: 0,
                },
                ChoiceValue::Base(b'N'),
            ),
        ];
        let (cfg, sim) = np_stop_filter_fixture();
        let (result, _, _) =
            run_np_replay(&pass, sim, records, Some(&contracts), Some(&cfg));
        match result.unwrap_err() {
            PassError::ConstraintSampling { address, .. } => {
                assert_eq!(address, "np.np1.bases[0]");
            }
            other => panic!("expected ConstraintSampling, got {other:?}"),
        }
        // Keep `productive` helper load-bearing (re-exported via tests).
        let _ = productive;
    }

    // ── Pass-level event-log emission ────────────────────────────

    #[test]
    fn generate_np_pass_emits_region_added_with_expected_np_region() {
        // GenerateNPPass with length=3, starting from an empty
        // simulation, must emit exactly one `RegionAdded` event
        // carrying an Np1 region from pool index 0..3.
        let pass = make_np_pass(3);
        let sim = Simulation::new();
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let mut captured: Vec<crate::ir::SimulationEvent> = Vec::new();
        let result = {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: None,
                contracts: None,
                feasibility: None,
                reference_index: None,
                replay_cursor: None,
                event_log_sink: Some(&mut captured),
            };
            pass.execute_checked(&sim, &mut ctx)
        };
        let sealed = result.expect("permissive NP execute should succeed");

        // The region the pass added.
        assert_eq!(sealed.sequence.region_count(), 1);
        let expected_region = sealed.sequence.regions[0].clone();
        assert_eq!(expected_region.segment, Segment::Np1);
        assert_eq!(expected_region.start, NucHandle::new(0));
        assert_eq!(expected_region.end, NucHandle::new(3));

        // Exactly one `RegionAdded` event with byte-identical
        // payload. (The base-push builder also emits one
        // `BasePushed` per NP base — 3 here — so the captured
        // stream contains those plus the region event.)
        let region_events: Vec<_> = captured
            .iter()
            .filter(|e| matches!(e, crate::ir::SimulationEvent::RegionAdded { .. }))
            .collect();
        assert_eq!(region_events.len(), 1);
        assert_eq!(
            *region_events[0],
            crate::ir::SimulationEvent::RegionAdded {
                region: expected_region
            }
        );
        let pushed = captured
            .iter()
            .filter(|e| matches!(e, crate::ir::SimulationEvent::BasePushed { .. }))
            .count();
        assert_eq!(pushed, 3, "3 NP bases pushed");
    }
}
