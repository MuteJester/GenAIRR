//! `AssembleSegmentPass` — copy a germline allele slice into the pool (C.8).

use crate::address;
#[cfg(test)]
use crate::ir::{GermlinePos, NucHandle};
use crate::ir::{Segment, Simulation};
use crate::pass::{Pass, PassContext, PassEffect, PassError, PassRequirement};

mod execution;

/// Assemble one germline segment (V, D, or J) from its assigned
/// `AlleleInstance` into the simulation pool.
///
/// Reads the allele's reference bases from
/// `ctx.refdata.expect(...)`, slices off `trim_5` from the start
/// and `trim_3` from the end, and pushes the retained bases as
/// germline-derived nucleotides into the pool. Then constructs a
/// `Region` for the segment with `frame_phase` chained from the
/// cumulative length of all prior regions. No codon-rail data is
/// stored on the region — callers that need translation call
/// [`crate::ir::compute_codon_rail`] against the post-assembly pool.
///
/// **Determinism:** the pass makes no random choices.
/// `declared_choices()` returns the empty vector — there is
/// nothing to record to the trace.
///
/// **Pre-conditions:**
/// - An allele must be assigned to `segment` (panics otherwise —
///   propagated from `Simulation::assignments.get`).
/// - `ctx.refdata` must be `Some(...)` (panics otherwise —
///   `PassContext::refdata` is `Option` so this pass can be
///   constructed in any plan, but it requires the data at
///   execute time).
/// - `trim_5 + trim_3 <= allele.len()` (panics otherwise — caller
///   bug, asserted defensively).
///
/// **Frame phase:** computed as `(cumulative prior region length)
/// % 3`. This anchors the codon frame at the start of the *first*
/// region. Real biology anchors at the V Cys; that anchor-aware
/// offset is a future refinement. The raw cumulative-length-mod-3
/// is correct as long as the user isn't reading frame-relative
/// properties yet.
pub struct AssembleSegmentPass {
    segment: Segment,
}

impl AssembleSegmentPass {
    /// Construct an assembly pass for the given segment.
    /// Panics if `segment` is not V, D, or J.
    pub fn new(segment: Segment) -> Self {
        match segment {
            Segment::V | Segment::D | Segment::J => {}
            _ => panic!(
                "AssembleSegmentPass: segment must be V, D, or J — got {:?}",
                segment
            ),
        }
        Self { segment }
    }

    pub fn segment(&self) -> Segment {
        self.segment
    }

    // `execute_with_validation` lives in the [`execution`] submodule;
    // see `execution.rs` for the body.
}

impl Pass for AssembleSegmentPass {
    fn name(&self) -> &str {
        address::assemble_vdj(self.segment)
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_validation(sim, ctx, false)
            .expect("AssembleSegmentPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_validation(sim, ctx, true)
    }

    fn declared_choices(&self) -> Vec<String> {
        // Deterministic — no choices.
        Vec::new()
    }

    fn requirements(&self) -> Vec<PassRequirement> {
        vec![
            PassRequirement::RefData,
            PassRequirement::AlleleAssignment(self.segment),
        ]
    }

    fn effects(&self) -> Vec<PassEffect> {
        vec![PassEffect::AssembleSegment(self.segment)]
    }
}

#[cfg(test)]
pub(super) mod test_support {
    use crate::refdata::{Allele, ChainType, RefDataConfig};

    use crate::ir::Segment;

    /// Build a tiny reference config with one V, D, and J allele
    /// each of known length / sequence for assembly tests.
    pub fn make_test_refdata() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_test*01".to_string(),
            gene: "v_test".to_string(),
            seq: b"AAACCCGGG".to_vec(), // 9 bases
            segment: Segment::V,
            anchor: Some(6), // C of CCC at position 6
            functional_status: None,
            subregions: Vec::new(),
        });
        let _ = cfg.d_pool.push(Allele {
            name: "d_test*01".to_string(),
            gene: "d_test".to_string(),
            seq: b"TTTTTT".to_vec(), // 6 bases
            segment: Segment::D,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_test*01".to_string(),
            gene: "j_test".to_string(),
            seq: b"GGGCCC".to_vec(), // 6 bases
            segment: Segment::J,
            anchor: Some(0),
            functional_status: None,
            subregions: Vec::new(),
        });
        cfg
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::make_test_refdata;
    use super::*;
    use crate::assignment::TrimEnd;
    use crate::dist::{AllelePoolDist, EmpiricalLengthDist};
    use crate::pass::testing::PassRuntime;
    use crate::pass::{PassError, PassPlan};
    use crate::passes::sample_allele::test_support::make_test_pool;
    use crate::passes::{SampleAllelePass, TrimPass};

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn assemble_segment_pass_rejects_np1() {
        let _ = AssembleSegmentPass::new(Segment::Np1);
    }

    #[test]
    #[should_panic(expected = "PassContext.refdata is None")]
    fn assemble_segment_pass_panics_without_refdata() {
        // Use plain `execute` (refdata = None) — should panic.
        let mut plan = PassPlan::new();
        let v_pool = make_test_pool(1, Segment::V);
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        let _ = PassRuntime::execute(&plan, Simulation::new(), 0);
    }

    #[test]
    fn assemble_segment_strict_errors_without_refdata() {
        let mut plan = PassPlan::new();
        let v_pool = make_test_pool(1, Segment::V);
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let err = PassRuntime::execute_strict_with_context(&plan, Simulation::new(), 0, None, None)
            .unwrap_err();

        assert_eq!(err.pass_name(), "assemble.v");
        assert!(matches!(err, PassError::MissingRefData { .. }));
    }

    #[test]
    #[should_panic(expected = "no allele assigned")]
    fn assemble_segment_pass_panics_without_allele_assignment() {
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        let _ = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
    }

    #[test]
    fn assemble_segment_strict_errors_without_allele_assignment() {
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let err = PassRuntime::execute_strict_with_context(
            &plan,
            Simulation::new(),
            0,
            Some(&refdata),
            None,
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "assemble.v");
        assert!(matches!(
            err,
            PassError::MissingAssignment {
                segment: Segment::V,
                ..
            }
        ));
    }

    #[test]
    #[should_panic(expected = "exceeds allele length")]
    fn assemble_segment_pass_panics_when_trim_exceeds_length() {
        let refdata = make_test_refdata();
        // V is 9 bases. Trim 5 + 3 = 8, leaves 1 base — fine.
        // Trim 5' = 5, 3' = 5 → sum 10 > 9 → panic.
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        // Note: SampleAllelePass uses make_test_pool, but AssembleSegmentPass
        // reads from refdata. They MUST be the same allele for the
        // test to succeed; we rely on test_refdata having allele 0
        // be the V whose length is 9.
        let _ = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
    }

    #[test]
    fn assemble_segment_strict_errors_when_trim_exceeds_length() {
        let refdata = make_test_refdata();
        let v_pool = make_test_pool(1, Segment::V);
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let err = PassRuntime::execute_strict_with_context(
            &plan,
            Simulation::new(),
            0,
            Some(&refdata),
            None,
        )
        .unwrap_err();

        assert_eq!(err.pass_name(), "assemble.v");
        assert!(matches!(err, PassError::InvalidPlanState { .. }));
    }

    #[test]
    fn assemble_segment_pass_v_no_trim_pushes_full_allele() {
        let refdata = make_test_refdata();
        // Use the V pool from refdata directly so SampleAllele sees
        // the same alleles as Assemble. We sample id 0 (the only V).
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        // V allele was 9 bases; no trim → 9 nucleotides in pool.
        assert_eq!(sim.pool.len(), 9);
        // One region (V) of length 9.
        assert_eq!(sim.sequence.region_count(), 1);
        let r = &sim.sequence.regions[0];
        assert_eq!(r.segment, Segment::V);
        assert_eq!(r.len(), 9);
        assert_eq!(r.frame_phase, 0); // first region

        // Bases match the allele sequence exactly, with germline_pos
        // matching position-in-allele.
        let expected = b"AAACCCGGG";
        for i in 0..9 {
            let n = sim.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.base, expected[i as usize]);
            assert_eq!(n.germline, expected[i as usize]);
            assert_eq!(n.germline_pos, GermlinePos::pos(i as u16));
            assert_eq!(n.segment, Segment::V);
        }
    }

    #[test]
    fn assemble_segment_pass_v_with_trim_pushes_post_trim_slice() {
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(2, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        // Allele = "AAACCCGGG" (9 bases); trim_5=2, trim_3=3.
        // Slice = "ACCCG" (positions 2..6, length 4 — wait let me recompute).
        // Actually 9 - 2 - 3 = 4 → slice is positions [2, 6) = "ACCC".
        assert_eq!(sim.pool.len(), 4);
        let expected = b"ACCC";
        for i in 0..4 {
            let n = sim.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.base, expected[i as usize]);
            assert_eq!(n.germline_pos, GermlinePos::pos((i + 2) as u16)); // shifted by trim_5
        }
    }

    #[test]
    fn assemble_segment_pass_chains_frame_phase_across_regions() {
        // V (9 bases, frame_phase=0) + assemble J (6 bases,
        // frame_phase = 9 % 3 = 0). So both regions have phase 0.
        // Test that the chain is correct after assembly.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        assert_eq!(sim.pool.len(), 15); // 9 V + 6 J
        assert_eq!(sim.sequence.region_count(), 2);
        assert_eq!(sim.sequence.regions[0].segment, Segment::V);
        assert_eq!(sim.sequence.regions[0].frame_phase, 0);
        assert_eq!(sim.sequence.regions[1].segment, Segment::J);
        // J starts after 9 V bases; 9 % 3 = 0 → frame_phase 0.
        assert_eq!(sim.sequence.regions[1].frame_phase, 0);
    }

    #[test]
    fn assemble_segment_pass_chains_frame_phase_with_offset() {
        // Build a plan where V has 2 bases (after heavy trim), so
        // J's frame_phase = 2 % 3 = 2. Tests the non-zero phase chain.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(SampleAllelePass::new(
            Segment::J,
            Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
        )));
        // V is 9 bases. Trim 5'=4, 3'=3 → slice has 9-4-3 = 2 bases.
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Five,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
        )));
        plan.push(Box::new(TrimPass::new(
            Segment::V,
            TrimEnd::Three,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();

        assert_eq!(sim.sequence.regions[0].len(), 2);
        assert_eq!(sim.sequence.regions[0].frame_phase, 0);
        // J's first base sits at position 2 of the cumulative chain;
        // 2 % 3 = 2 → frame_phase 2.
        assert_eq!(sim.sequence.regions[1].frame_phase, 2);
    }

    #[test]
    fn assemble_segment_pass_does_not_record_to_trace() {
        // Determinism contract for transform passes: no trace
        // entries.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 7, &refdata);

        // Trace contains only the SampleAllele entry, none from Assemble.
        assert_eq!(outcome.trace.len(), 1);
        assert_eq!(outcome.trace.choices()[0].address, "sample_allele.v");
    }

    #[test]
    fn assemble_segment_pass_declared_choices_is_empty() {
        let pass = AssembleSegmentPass::new(Segment::V);
        assert!(pass.declared_choices().is_empty());
    }

    #[test]
    fn assemble_segment_pass_codon_rail_computed_on_assembled_region() {
        // V allele "AAACCCGGG" (9 bases), no trim, frame_phase 0
        // → codons AAA, CCC, GGG → K, P, G.
        let refdata = make_test_refdata();
        let mut plan = PassPlan::new();
        plan.push(Box::new(SampleAllelePass::new(
            Segment::V,
            Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
        )));
        plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));

        let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);
        let sim = outcome.final_simulation();
        // Codon rail is computed on demand at consumer call sites.
        let rail = crate::ir::compute_codon_rail(&sim.sequence.regions[0], &sim.pool);
        assert_eq!(rail.amino_acids, b"KPG");
        assert!(rail.stop_codon_positions.is_empty());
    }

    // ── Pass-level event-log emission ────────────────────────────

    #[test]
    fn assemble_segment_pass_emits_region_added_with_expected_v_region() {
        // AssembleSegmentPass for a 9-base V allele with zero
        // trims must emit exactly one `RegionAdded` event whose
        // region matches the V region the pass actually installs.
        use crate::assignment::AlleleInstance;
        use crate::pass::PassContext;
        use crate::refdata::AlleleId;
        use crate::rng::Rng;
        use crate::trace::Trace;

        let refdata = make_test_refdata();
        // Pre-assign V allele 0 (zero trims via `AlleleInstance::new`).
        let sim = Simulation::new().with_allele_assigned(
            Segment::V,
            AlleleInstance::new(AlleleId::new(0)),
        );

        let pass = AssembleSegmentPass::new(Segment::V);
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let mut captured: Vec<crate::ir::SimulationEvent> = Vec::new();
        let result = {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: Some(&refdata),
                contracts: None,
                feasibility: None,
                reference_index: None,
                replay_cursor: None,
                event_log_sink: Some(&mut captured),
            };
            pass.execute_checked(&sim, &mut ctx)
        };
        let sealed = result.expect("assembly should succeed with V allele assigned");

        // The pass should have appended exactly one region.
        assert_eq!(sealed.sequence.region_count(), 1);
        let expected_region = sealed.sequence.regions[0].clone();
        assert_eq!(expected_region.segment, Segment::V);
        assert_eq!(expected_region.start, NucHandle::new(0));
        assert_eq!(expected_region.end, NucHandle::new(9));

        // Exactly one `RegionAdded` event with byte-identical
        // payload. (The base-push builder also emits one
        // `BasePushed` per germline base — 9 here — so the
        // captured stream contains those plus the region event.)
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
        // Sanity: the BasePushed events are also in the stream.
        let pushed = captured
            .iter()
            .filter(|e| matches!(e, crate::ir::SimulationEvent::BasePushed { .. }))
            .count();
        assert_eq!(pushed, 9, "9 germline bases pushed during V assembly");
    }

    // ── Slice B: D-inversion assembly emission ───────────────────

    /// Tiny VDJ cartridge with a caller-chosen D allele sequence.
    /// V and J share a 9-base / 6-base default so the fixtures stay
    /// human-readable; tests of D orientation only consult the D
    /// pool, so V / J shape doesn't matter beyond having any allele
    /// available for the chain.
    fn make_refdata_with_d_seq(d_seq: &[u8]) -> crate::refdata::RefDataConfig {
        use crate::refdata::{Allele, ChainType, RefDataConfig};
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        cfg.rules.v_anchor.required = false;
        cfg.rules.j_anchor.required = false;
        let _ = cfg.v_pool.push(Allele {
            name: "v*01".into(),
            gene: "v".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        });
        let _ = cfg.d_pool.push(Allele {
            name: "d*01".into(),
            gene: "d".into(),
            seq: d_seq.to_vec(),
            segment: Segment::D,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"GGGCCC".to_vec(),
            segment: Segment::J,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        });
        cfg
    }

    /// Drive `AssembleSegmentPass(segment)` against a pre-staged
    /// `Simulation` that already carries the allele assignment we
    /// want to exercise. No SampleAllelePass; no reference walker.
    /// Returns the post-pass simulation.
    fn run_assemble_against_sim(
        sim: Simulation,
        segment: Segment,
        refdata: &crate::refdata::RefDataConfig,
    ) -> Simulation {
        use crate::pass::PassContext;
        use crate::rng::Rng;
        use crate::trace::Trace;

        let pass = AssembleSegmentPass::new(segment);
        let mut trace = Trace::new();
        let mut rng = Rng::new(0);
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: Some(refdata),
            contracts: None,
            feasibility: None,
            reference_index: None,
            replay_cursor: None,
            event_log_sink: None,
        };
        pass.execute(&sim, &mut ctx)
    }

    #[test]
    fn forward_d_emits_unchanged_slice_with_ascending_germline_pos() {
        // Pin the byte-identical guarantee for the forward path:
        // every existing caller commits Forward, so emitted bytes
        // must equal the retained slice with ascending germline_pos
        // and zero flags.
        use crate::assignment::AlleleInstance;
        use crate::refdata::AlleleId;

        let refdata = make_refdata_with_d_seq(b"ACGTTA");
        let sim = Simulation::new()
            .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)));
        let after = run_assemble_against_sim(sim, Segment::D, &refdata);

        assert_eq!(after.pool.len(), 6);
        let expected = b"ACGTTA";
        for i in 0..6 {
            let n = after.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.base, expected[i as usize]);
            assert_eq!(n.germline, expected[i as usize]);
            assert_eq!(n.germline_pos, GermlinePos::pos(i as u16));
            assert_eq!(n.segment, Segment::D);
            assert!(!n.flags.contains(crate::ir::flag::INVERTED));
        }
    }

    #[test]
    fn reverse_complement_d_emits_rc_bytes_with_descending_germline_pos() {
        // D = `ACGTTA` (6 bases), no trims, orientation
        // ReverseComplement. Expected emitted bytes are the
        // reverse-complement walk: position 5=A→T, 4=T→A, 3=T→A,
        // 2=G→C, 1=C→G, 0=A→T → `TAACGT`. germline_pos descends
        // 5,4,3,2,1,0. Every base carries NucFlags::INVERTED.
        use crate::assignment::{AlleleInstance, SegmentOrientation};
        use crate::refdata::AlleleId;

        let refdata = make_refdata_with_d_seq(b"ACGTTA");
        let sim = Simulation::new()
            .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_orientation(Segment::D, SegmentOrientation::ReverseComplement);
        let after = run_assemble_against_sim(sim, Segment::D, &refdata);

        assert_eq!(after.pool.len(), 6);
        let expected_bases = b"TAACGT";
        let expected_germline_pos = [5u16, 4, 3, 2, 1, 0];
        for i in 0..6 {
            let n = after.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(
                n.base, expected_bases[i as usize],
                "byte {i} should be the reverse-complement of the source allele",
            );
            // `germline` mirrors the emitted base — the pre-mutation
            // reference at this pool slot, which under inversion is
            // already the complemented byte.
            assert_eq!(n.germline, expected_bases[i as usize]);
            assert_eq!(
                n.germline_pos,
                GermlinePos::pos(expected_germline_pos[i as usize]),
                "germline_pos at pool index {i} must point to the source allele position",
            );
            assert_eq!(n.segment, Segment::D);
            assert!(
                n.flags.contains(crate::ir::flag::INVERTED),
                "every byte from an inverted D assembly must carry NucFlags::INVERTED",
            );
        }
    }

    #[test]
    fn reverse_complement_d_respects_trims() {
        // D = `ACGTTAGC` (8 bases), trim_5 = 1, trim_3 = 2.
        // Retained allele slice = positions [1, 6) = `CGTTA`.
        // Forward emission would push `CGTTA` with germline_pos
        // 1,2,3,4,5. Inverted emission walks the retained slice
        // from position 5 down to 1, complementing each byte:
        //   pos 5 = A → T,  pos 4 = T → A,  pos 3 = T → A,
        //   pos 2 = G → C,  pos 1 = C → G
        // → emitted bases `TAACG`, germline_pos descending 5,4,3,2,1.
        use crate::assignment::{AlleleInstance, SegmentOrientation};
        use crate::refdata::AlleleId;

        let refdata = make_refdata_with_d_seq(b"ACGTTAGC");
        let inst = AlleleInstance::new(AlleleId::new(0))
            .with_trim_5(1)
            .with_trim_3(2);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::D, inst)
            .with_allele_orientation(Segment::D, SegmentOrientation::ReverseComplement);
        let after = run_assemble_against_sim(sim, Segment::D, &refdata);

        assert_eq!(after.pool.len(), 5, "post-trim slice is 8 - 1 - 2 = 5 bases");
        let expected_bases = b"TAACG";
        let expected_germline_pos = [5u16, 4, 3, 2, 1];
        for i in 0..5 {
            let n = after.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.base, expected_bases[i as usize]);
            assert_eq!(
                n.germline_pos,
                GermlinePos::pos(expected_germline_pos[i as usize]),
                "inverted-D walk must descend over the retained slice's \
                 original-allele coordinates only — never out of the slice",
            );
            assert!(n.flags.contains(crate::ir::flag::INVERTED));
        }
        // And the per-region attrs match the forward case
        // structurally — same length, same segment, same start
        // position in the pool.
        let region = &after.sequence.regions[0];
        assert_eq!(region.segment, Segment::D);
        assert_eq!(region.start, NucHandle::new(0));
        assert_eq!(region.end, NucHandle::new(5));
    }

    #[test]
    fn v_orientation_reverse_complement_is_ignored_by_assembly() {
        // V (and J — see neighbouring test) are deliberately
        // scope-narrowed out of Slice B / v1 D inversion. Even when
        // a caller manually flips V's orientation via the IR
        // primitive, assembly must continue to emit the forward
        // slice. This guards against accidental broadening before
        // an explicit design pass adds V/J inversion.
        use crate::assignment::{AlleleInstance, SegmentOrientation};
        use crate::refdata::AlleleId;

        // Build a tiny cfg with a non-palindromic V so forward vs
        // RC is observable byte-wise.
        let refdata = {
            use crate::refdata::{Allele, ChainType, RefDataConfig};
            let mut cfg = RefDataConfig::empty(ChainType::Vdj);
            cfg.rules.v_anchor.required = false;
            cfg.rules.j_anchor.required = false;
            let _ = cfg.v_pool.push(Allele {
                name: "v*01".into(),
                gene: "v".into(),
                seq: b"ACGTAA".to_vec(),
                segment: Segment::V,
                anchor: None,
                functional_status: None,
                subregions: Vec::new(),
            });
            let _ = cfg.d_pool.push(Allele {
                name: "d*01".into(),
                gene: "d".into(),
                seq: b"GGG".to_vec(),
                segment: Segment::D,
                anchor: None,
                functional_status: None,
                subregions: Vec::new(),
            });
            let _ = cfg.j_pool.push(Allele {
                name: "j*01".into(),
                gene: "j".into(),
                seq: b"AAA".to_vec(),
                segment: Segment::J,
                anchor: None,
                functional_status: None,
                subregions: Vec::new(),
            });
            cfg
        };

        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_orientation(Segment::V, SegmentOrientation::ReverseComplement);
        let after = run_assemble_against_sim(sim, Segment::V, &refdata);

        // Forward emission expected — RC orientation on V is
        // intentionally a no-op.
        assert_eq!(after.pool.len(), 6);
        let expected = b"ACGTAA";
        for i in 0..6 {
            let n = after.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.base, expected[i as usize]);
            assert_eq!(n.germline_pos, GermlinePos::pos(i as u16));
            assert!(
                !n.flags.contains(crate::ir::flag::INVERTED),
                "V emission must NEVER carry INVERTED in v1 — D inversion only",
            );
        }
    }

    #[test]
    fn j_orientation_reverse_complement_is_ignored_by_assembly() {
        // Mirror of the V test. Same scope rationale: J inversion
        // is out of scope for v1; the IR setter is defensive but
        // assembly must not act on it for J.
        use crate::assignment::{AlleleInstance, SegmentOrientation};
        use crate::refdata::AlleleId;

        let refdata = {
            use crate::refdata::{Allele, ChainType, RefDataConfig};
            let mut cfg = RefDataConfig::empty(ChainType::Vdj);
            cfg.rules.v_anchor.required = false;
            cfg.rules.j_anchor.required = false;
            let _ = cfg.v_pool.push(Allele {
                name: "v*01".into(),
                gene: "v".into(),
                seq: b"AAA".to_vec(),
                segment: Segment::V,
                anchor: None,
                functional_status: None,
                subregions: Vec::new(),
            });
            let _ = cfg.d_pool.push(Allele {
                name: "d*01".into(),
                gene: "d".into(),
                seq: b"GGG".to_vec(),
                segment: Segment::D,
                anchor: None,
                functional_status: None,
                subregions: Vec::new(),
            });
            let _ = cfg.j_pool.push(Allele {
                name: "j*01".into(),
                gene: "j".into(),
                seq: b"ACGTAA".to_vec(),
                segment: Segment::J,
                anchor: None,
                functional_status: None,
                subregions: Vec::new(),
            });
            cfg
        };

        let sim = Simulation::new()
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_orientation(Segment::J, SegmentOrientation::ReverseComplement);
        let after = run_assemble_against_sim(sim, Segment::J, &refdata);

        assert_eq!(after.pool.len(), 6);
        let expected = b"ACGTAA";
        for i in 0..6 {
            let n = after.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(n.base, expected[i as usize]);
            assert_eq!(n.germline_pos, GermlinePos::pos(i as u16));
            assert!(!n.flags.contains(crate::ir::flag::INVERTED));
        }
    }

    #[test]
    fn forward_d_path_is_byte_identical_to_pre_slice_b_emission() {
        // Lock-in test for the byte-identical Slice B guarantee:
        // an `AlleleInstance::new(...)` (default Forward) must
        // produce the same emitted pool — base / germline /
        // germline_pos / segment / flags — as the legacy single-
        // path emission did before Slice B introduced the branch.
        //
        // The comparison fixture is just the existing
        // `make_test_refdata` D allele (`TTTTTT`, 6 bases,
        // anchor=None), driven through the trim path with the
        // same trims any production pipeline would use.
        use crate::assignment::{AlleleInstance, SegmentOrientation};
        use crate::refdata::AlleleId;

        let refdata = make_test_refdata();
        // Forward (default) path — control.
        let default_inst = AlleleInstance::new(AlleleId::new(0));
        let sim_forward = Simulation::new()
            .with_allele_assigned(Segment::D, default_inst);
        let after_forward = run_assemble_against_sim(sim_forward, Segment::D, &refdata);

        // Forward with explicit Forward orientation — must match
        // the implicit-Forward control byte-for-byte.
        let explicit_inst = AlleleInstance::new(AlleleId::new(0))
            .with_orientation(SegmentOrientation::Forward);
        let sim_explicit = Simulation::new()
            .with_allele_assigned(Segment::D, explicit_inst);
        let after_explicit = run_assemble_against_sim(sim_explicit, Segment::D, &refdata);

        assert_eq!(after_forward.pool.len(), after_explicit.pool.len());
        for i in 0..after_forward.pool.len() as u32 {
            let a = after_forward.pool.get(NucHandle::new(i)).unwrap();
            let b = after_explicit.pool.get(NucHandle::new(i)).unwrap();
            assert_eq!(a.base, b.base);
            assert_eq!(a.germline, b.germline);
            assert_eq!(a.germline_pos, b.germline_pos);
            assert_eq!(a.segment, b.segment);
            assert_eq!(a.flags, b.flags);
        }
    }
}
