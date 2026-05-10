//! End-to-end heavy-chain recombination integration test.
//!
//! Composes the full recombination pass family into a single VDJ
//! recombination plan and verifies the result has the expected
//! structure: V → NP1 → D → NP2 → J regions in order, correct
//! pool length, correct germline_pos provenance on V/D/J bases,
//! NP bases tagged as synthetic N-nucs, frame phases chained
//! correctly across regions, and the `AnchorPreserved` contract
//! satisfied on the assembled IR.
//!
//! Lives in `engine_rs/tests/` so it exercises the public
//! crate API the way an external consumer would.

use genairr_engine::assignment::TrimEnd;
use genairr_engine::contract::{AnchorPreserved, Contract};
use genairr_engine::dist::{AllelePoolDist, EmpiricalLengthDist, UniformBase};
use genairr_engine::ir::{flag, NucHandle, Nucleotide, Segment, Simulation};
use genairr_engine::pass::{PassPlan, PassRuntime};
use genairr_engine::passes::{AssembleSegmentPass, GenerateNPPass, SampleAllelePass, TrimPass};
use genairr_engine::refdata::{Allele, ChainType, RefDataConfig};
use genairr_engine::trace::ChoiceValue;

/// Build a tiny synthetic heavy-chain reference config.
///
/// Sequences are deliberately short and have lengths that are
/// multiples of 3 so the frame-phase math at region boundaries is
/// trivial to reason about by hand.
///
/// - V: "AAACCCGGG" (9 bp, anchor at position 6)
/// - D: "TTTTTT" (6 bp, anchorless)
/// - J: "GGGCCC" (6 bp, anchor at position 0)
fn make_synthetic_vdj_refdata() -> RefDataConfig {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(Allele {
        name: "v_test*01".to_string(),
        gene: "v_test".to_string(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
    });
    let _ = cfg.d_pool.push(Allele {
        name: "d_test*01".to_string(),
        gene: "d_test".to_string(),
        seq: b"TTTTTT".to_vec(),
        segment: Segment::D,
        anchor: None,
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j_test*01".to_string(),
        gene: "j_test".to_string(),
        seq: b"GGGCCC".to_vec(),
        segment: Segment::J,
        anchor: Some(0),
    });
    cfg
}

/// Build a complete heavy-chain recombination plan against the
/// synthetic refdata. No trims, NP1/NP2 = 3 bases each.
///
/// Pass order is the canonical biological flow:
///   1. SampleV / SampleD / SampleJ
///   2. (no trim — defaults to 0)
///   3. AssembleV → GenerateNP1 → AssembleD → GenerateNP2 → AssembleJ
fn build_basic_vdj_plan(refdata: &RefDataConfig) -> PassPlan {
    let mut plan = PassPlan::new();

    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::uniform(&refdata.d_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
    )));

    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np2,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    plan
}

#[test]
fn end_to_end_vdj_produces_expected_region_structure() {
    let refdata = make_synthetic_vdj_refdata();
    let plan = build_basic_vdj_plan(&refdata);
    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &refdata);

    let sim = outcome.final_simulation();

    // 8 passes → 9 IR revisions (initial + 8). Note: SampleAllele,
    // AssembleSegment, and GenerateNP each produce one revision —
    // 3 + 5 = 8 passes, 9 revisions.
    assert_eq!(outcome.revisions.len(), 9);

    // 5 regions: V → NP1 → D → NP2 → J.
    assert_eq!(sim.sequence.region_count(), 5);
    let segments: Vec<Segment> = sim.sequence.regions.iter().map(|r| r.segment).collect();
    assert_eq!(
        segments,
        vec![
            Segment::V,
            Segment::Np1,
            Segment::D,
            Segment::Np2,
            Segment::J
        ]
    );

    // Pool: 9 (V) + 3 (NP1) + 6 (D) + 3 (NP2) + 6 (J) = 27 bases.
    assert_eq!(sim.pool.len(), 27);
}

#[test]
fn end_to_end_vdj_region_ranges_chain_correctly() {
    let refdata = make_synthetic_vdj_refdata();
    let plan = build_basic_vdj_plan(&refdata);
    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &refdata);
    let sim = outcome.final_simulation();

    // V: [0, 9).
    assert_eq!(sim.sequence.regions[0].start, NucHandle::new(0));
    assert_eq!(sim.sequence.regions[0].end, NucHandle::new(9));
    // NP1: [9, 12).
    assert_eq!(sim.sequence.regions[1].start, NucHandle::new(9));
    assert_eq!(sim.sequence.regions[1].end, NucHandle::new(12));
    // D: [12, 18).
    assert_eq!(sim.sequence.regions[2].start, NucHandle::new(12));
    assert_eq!(sim.sequence.regions[2].end, NucHandle::new(18));
    // NP2: [18, 21).
    assert_eq!(sim.sequence.regions[3].start, NucHandle::new(18));
    assert_eq!(sim.sequence.regions[3].end, NucHandle::new(21));
    // J: [21, 27).
    assert_eq!(sim.sequence.regions[4].start, NucHandle::new(21));
    assert_eq!(sim.sequence.regions[4].end, NucHandle::new(27));
}

#[test]
fn end_to_end_vdj_frame_phases_are_zero_for_all_multiple_of_three_lengths() {
    // Every region's start position is a multiple of 3 (because all
    // segment lengths in this test are multiples of 3), so every
    // region's frame_phase should be 0.
    let refdata = make_synthetic_vdj_refdata();
    let plan = build_basic_vdj_plan(&refdata);
    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &refdata);
    let sim = outcome.final_simulation();

    for (i, r) in sim.sequence.regions.iter().enumerate() {
        assert_eq!(
            r.frame_phase, 0,
            "region[{}] ({:?}) frame_phase should be 0 in this test fixture",
            i, r.segment
        );
    }
}

#[test]
fn end_to_end_vdj_germline_provenance_is_preserved_for_alleles() {
    let refdata = make_synthetic_vdj_refdata();
    let plan = build_basic_vdj_plan(&refdata);
    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &refdata);
    let sim = outcome.final_simulation();

    // V positions [0, 9) have germline_pos = original allele
    // position (also 0..9 since trim_5 == 0).
    for i in 0..9 {
        let n = sim.pool.get(NucHandle::new(i)).unwrap();
        assert_eq!(n.segment, Segment::V);
        assert_eq!(n.germline_pos, i as u16);
        assert_eq!(n.base, b"AAACCCGGG"[i as usize]);
    }

    // NP1 positions [9, 12): synthetic, no germline provenance.
    for i in 9..12 {
        let n = sim.pool.get(NucHandle::new(i)).unwrap();
        assert_eq!(n.segment, Segment::Np1);
        assert_eq!(n.germline_pos, Nucleotide::NO_GERMLINE_POS);
        assert!(n.flags.contains(flag::N_NUC));
        assert!(matches!(n.base, b'A' | b'C' | b'G' | b'T'));
    }

    // D positions [12, 18): germline_pos = position-in-D-allele
    // (0..6 since trim_5 == 0).
    for (offset, pool_idx) in (12..18).enumerate() {
        let n = sim.pool.get(NucHandle::new(pool_idx)).unwrap();
        assert_eq!(n.segment, Segment::D);
        assert_eq!(n.germline_pos, offset as u16);
        assert_eq!(n.base, b'T');
    }

    // NP2 positions [18, 21): synthetic.
    for i in 18..21 {
        let n = sim.pool.get(NucHandle::new(i)).unwrap();
        assert_eq!(n.segment, Segment::Np2);
        assert_eq!(n.germline_pos, Nucleotide::NO_GERMLINE_POS);
        assert!(n.flags.contains(flag::N_NUC));
    }

    // J positions [21, 27): germline_pos 0..6 (no trim).
    for (offset, pool_idx) in (21..27).enumerate() {
        let n = sim.pool.get(NucHandle::new(pool_idx)).unwrap();
        assert_eq!(n.segment, Segment::J);
        assert_eq!(n.germline_pos, offset as u16);
        assert_eq!(n.base, b"GGGCCC"[offset]);
    }
}

#[test]
fn end_to_end_vdj_anchors_are_preserved_under_no_trim() {
    let refdata = make_synthetic_vdj_refdata();
    let plan = build_basic_vdj_plan(&refdata);
    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &refdata);
    let sim = outcome.final_simulation();

    let v_contract = AnchorPreserved::new(Segment::V);
    let j_contract = AnchorPreserved::new(Segment::J);

    assert!(v_contract.verify(sim, Some(&refdata)).is_ok());
    assert!(j_contract.verify(sim, Some(&refdata)).is_ok());
}

#[test]
fn end_to_end_vdj_anchor_violation_detected_when_v_trim_3_eats_anchor() {
    // V allele has anchor at position 6 in a 9-base allele. trim_3 = 4
    // means retained_end = 5, which is less than anchor + 3 = 9 → violation.
    let refdata = make_synthetic_vdj_refdata();
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(4, 1.0)])),
    )));
    // Note: no Assemble pass — the contract works on the
    // assignments alone.
    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 0, &refdata);

    let v_contract = AnchorPreserved::new(Segment::V);
    let result = v_contract.verify(outcome.final_simulation(), Some(&refdata));
    assert!(result.is_err());
    let err = result.unwrap_err();
    assert_eq!(err.contract_name, "anchor_preserved.v");
}

#[test]
fn end_to_end_vdj_trace_contains_expected_addresses() {
    let refdata = make_synthetic_vdj_refdata();
    let plan = build_basic_vdj_plan(&refdata);
    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &refdata);

    // Sampling: 3 entries (sample_allele.v/d/j).
    assert!(outcome.trace.find("sample_allele.v").is_some());
    assert!(outcome.trace.find("sample_allele.d").is_some());
    assert!(outcome.trace.find("sample_allele.j").is_some());

    // NP generation: 2 length entries + 6 base entries (3 each
    // for NP1 and NP2).
    assert!(outcome.trace.find("np.np1.length").is_some());
    assert!(outcome.trace.find("np.np2.length").is_some());
    for i in 0..3 {
        assert!(
            outcome
                .trace
                .find(&format!("np.np1.bases[{}]", i))
                .is_some(),
            "missing np.np1.bases[{}]",
            i
        );
        assert!(
            outcome
                .trace
                .find(&format!("np.np2.bases[{}]", i))
                .is_some(),
            "missing np.np2.bases[{}]",
            i
        );
    }

    // No trim entries (no TrimPass in plan).
    assert!(outcome.trace.find("trim.v_3").is_none());
    assert!(outcome.trace.find("trim.j_5").is_none());

    // No assembly entries (assembly is a transform pass).
    assert!(outcome.trace.find("assemble.v").is_none());

    // Total: 3 + 2 + 6 = 11 trace entries.
    assert_eq!(outcome.trace.len(), 11);
}

#[test]
fn end_to_end_vdj_codon_rail_correct_for_v_region() {
    // V region with frame_phase 0 + bases "AAACCCGGG" → codons
    // AAA, CCC, GGG → K, P, G.
    let refdata = make_synthetic_vdj_refdata();
    let plan = build_basic_vdj_plan(&refdata);
    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &refdata);
    let sim = outcome.final_simulation();
    let v = &sim.sequence.regions[0];
    assert_eq!(v.segment, Segment::V);
    assert_eq!(v.amino_acids, b"KPG");
    assert!(v.stop_codon_positions.is_empty());
}

#[test]
fn end_to_end_vdj_is_deterministic_under_same_seed() {
    let refdata = make_synthetic_vdj_refdata();

    let oa = PassRuntime::execute_with_refdata(
        &build_basic_vdj_plan(&refdata),
        Simulation::new(),
        0xc0ff_ee,
        &refdata,
    );
    let ob = PassRuntime::execute_with_refdata(
        &build_basic_vdj_plan(&refdata),
        Simulation::new(),
        0xc0ff_ee,
        &refdata,
    );

    // Trace identical.
    assert_eq!(oa.trace.choices(), ob.trace.choices());

    // Final pool identical, byte for byte.
    assert_eq!(
        oa.final_simulation().pool.len(),
        ob.final_simulation().pool.len()
    );
    for i in 0..oa.final_simulation().pool.len() {
        let h = NucHandle::new(i as u32);
        let na = oa.final_simulation().pool.get(h).unwrap();
        let nb = ob.final_simulation().pool.get(h).unwrap();
        assert_eq!(na.base, nb.base);
        assert_eq!(na.segment, nb.segment);
        assert_eq!(na.germline_pos, nb.germline_pos);
        assert_eq!(na.flags, nb.flags);
    }

    // Region structure identical.
    assert_eq!(oa.final_simulation().sequence.region_count(), 5);
    assert_eq!(
        oa.final_simulation().sequence.region_count(),
        ob.final_simulation().sequence.region_count()
    );
    for (ra, rb) in oa
        .final_simulation()
        .sequence
        .regions
        .iter()
        .zip(ob.final_simulation().sequence.regions.iter())
    {
        assert_eq!(ra.segment, rb.segment);
        assert_eq!(ra.start, rb.start);
        assert_eq!(ra.end, rb.end);
        assert_eq!(ra.frame_phase, rb.frame_phase);
        assert_eq!(ra.amino_acids, rb.amino_acids);
    }
}

#[test]
fn end_to_end_vdj_with_trims_assembles_post_trim_slices_and_chains_phases() {
    // Heavy chain with non-zero trims that introduce a non-zero
    // frame-phase chain. Trim V trim_3 = 1, J trim_5 = 1. NP1 = 3, NP2 = 3.
    //
    // V slice: 9 - 0 - 1 = 8 bases (AAACCCGG)
    // D slice: 6 bases (TTTTTT)
    // J slice: 6 - 1 - 0 = 5 bases (GGCCC)
    //
    // Pool: 8 + 3 + 6 + 3 + 5 = 25.
    // Frame phases:
    //   V starts at 0 → phase 0
    //   NP1 starts at 8 → 8 % 3 = 2
    //   D starts at 11 → 11 % 3 = 2
    //   NP2 starts at 17 → 17 % 3 = 2
    //   J starts at 20 → 20 % 3 = 2
    let refdata = make_synthetic_vdj_refdata();
    let mut plan = PassPlan::new();
    plan.push(Box::new(SampleAllelePass::new(
        Segment::V,
        Box::new(AllelePoolDist::uniform(&refdata.v_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::D,
        Box::new(AllelePoolDist::uniform(&refdata.d_pool)),
    )));
    plan.push(Box::new(SampleAllelePass::new(
        Segment::J,
        Box::new(AllelePoolDist::uniform(&refdata.j_pool)),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::V,
        TrimEnd::Three,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
    )));
    plan.push(Box::new(TrimPass::new(
        Segment::J,
        TrimEnd::Five,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::V)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np1,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::D)));
    plan.push(Box::new(GenerateNPPass::new(
        Segment::Np2,
        Box::new(EmpiricalLengthDist::from_pairs(vec![(3, 1.0)])),
        Box::new(UniformBase),
    )));
    plan.push(Box::new(AssembleSegmentPass::new(Segment::J)));

    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 7, &refdata);
    let sim = outcome.final_simulation();

    assert_eq!(sim.pool.len(), 25);
    assert_eq!(sim.sequence.region_count(), 5);
    assert_eq!(sim.sequence.regions[0].len(), 8);
    assert_eq!(sim.sequence.regions[1].len(), 3);
    assert_eq!(sim.sequence.regions[2].len(), 6);
    assert_eq!(sim.sequence.regions[3].len(), 3);
    assert_eq!(sim.sequence.regions[4].len(), 5);

    // Frame phases: V=0, then chain (V len 8) → 2, 2, 2, 2.
    assert_eq!(sim.sequence.regions[0].frame_phase, 0);
    assert_eq!(sim.sequence.regions[1].frame_phase, 2);
    assert_eq!(sim.sequence.regions[2].frame_phase, 2);
    assert_eq!(sim.sequence.regions[3].frame_phase, 2);
    assert_eq!(sim.sequence.regions[4].frame_phase, 2);

    // V anchor still preserved (anchor at 6, retained range
    // [0, 8) → 6 + 3 = 9 > 8 → BUT the anchor codon ends at 9
    // which is past the trim cut, so this should *fail*).
    let v_contract = AnchorPreserved::new(Segment::V);
    assert!(v_contract.verify(sim, Some(&refdata)).is_err());

    // J anchor: anchor at 0 with trim_5 = 1 → anchor 5'-trimmed
    // → fails.
    let j_contract = AnchorPreserved::new(Segment::J);
    assert!(j_contract.verify(sim, Some(&refdata)).is_err());

    // The bases of the V slice are verified by germline_pos.
    // Position 0..8 of V are AAACCCGG (positions 0..8 of the original
    // 9-base allele).
    let v_expected = b"AAACCCGG";
    for i in 0..8 {
        let n = sim.pool.get(NucHandle::new(i)).unwrap();
        assert_eq!(n.base, v_expected[i as usize]);
        assert_eq!(n.germline_pos, i as u16);
    }
}

#[test]
fn end_to_end_vdj_trace_value_addresses_match_sample_allele_ids() {
    // The AlleleId recorded in the trace must equal the AlleleId
    // ultimately found in the simulation's assignments — sanity
    // check on the addressed-trace faithfulness for sampling.
    let refdata = make_synthetic_vdj_refdata();
    let plan = build_basic_vdj_plan(&refdata);
    let outcome = PassRuntime::execute_with_refdata(&plan, Simulation::new(), 42, &refdata);
    let sim = outcome.final_simulation();

    for (segment, addr) in [
        (Segment::V, "sample_allele.v"),
        (Segment::D, "sample_allele.d"),
        (Segment::J, "sample_allele.j"),
    ] {
        let recorded = match outcome.trace.find(addr).unwrap().value {
            ChoiceValue::AlleleId(id) => id,
            _ => panic!("wrong variant at {}", addr),
        };
        let assigned = sim.assignments.get(segment).unwrap().allele_id.index();
        assert_eq!(
            recorded, assigned,
            "trace at {} = {} but assignment.allele_id = {}",
            addr, recorded, assigned
        );
    }
}
