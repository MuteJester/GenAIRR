use super::super::{build_airr_record, AirrRecord};
use super::*;
use crate::live_call::{
    AlleleBitSet, EvidenceScore, HypothesisFlags, PlacementHypothesis, SegmentCalls,
    SegmentLiveCall,
};

#[test]
fn airr_call_projection_falls_back_to_origin_assignment_without_live_call() {
    let (cfg, sim, _v0, _v1) = call_projection_fixture();
    let outcome = outcome_from_sim(sim);
    let rec = build_airr_record(&outcome, &cfg, "fallback");

    assert_eq!(rec.v_call, "IGHV1-1*01");
    assert_eq!(rec.locus, "IGH");
}

#[test]
fn airr_call_projection_prefers_live_call_allele_set() {
    let (cfg, sim, v0, v1) = call_projection_fixture();
    let hypothesis = PlacementHypothesis::new(
        Segment::V,
        0,
        6,
        0,
        6,
        AlleleBitSet::from_ids(cfg.v_pool.len(), [v0, v1]),
        EvidenceScore::exact(6, 0),
        HypothesisFlags::EMPTY,
    );
    let call = SegmentLiveCall::from_hypotheses(Segment::V, cfg.v_pool.len(), vec![hypothesis], 1);
    let live = SegmentCalls::empty().with_segment_call(call);
    let outcome = outcome_from_sim(sim.with_segment_calls(live));

    let rec = build_airr_record(&outcome, &cfg, "live");

    assert_eq!(rec.v_call, "IGHV1-1*01,IGHV1-1*02");
    assert_eq!(rec.locus, "IGH");
}

#[test]
fn airr_call_projection_falls_back_to_truth_when_live_call_is_unsupported() {
    let (cfg, sim, _v0, _v1) = call_projection_fixture();
    let call = SegmentLiveCall::from_hypotheses(Segment::V, cfg.v_pool.len(), Vec::new(), 1);
    let live = SegmentCalls::empty().with_segment_call(call);
    let outcome = outcome_from_sim(sim.with_segment_calls(live));

    let rec = build_airr_record(&outcome, &cfg, "unsupported");

    assert_eq!(rec.v_call, "IGHV1-1*01");
}

// ──────────────────────────────────────────────────────────────────
// End-loss interaction with the fallback germline-coord path.
//
// Investigation of audit §6.3 surfaced that the trim-only fallback
// already matches the walker's AIRR semantic: the walker encodes
// a 5'-end-loss as `D` (deletion) ops at the start of V's CIGAR,
// keeping `v_germline_start = v_trim_5` regardless of how many
// pool bytes the end-loss physically removed. The 3' J side
// mirrors. So the fallback's pre-existing math (also keying off
// `v_trim_5` / `j_trim_3`) was already coherent.
//
// The tests below pin that parity invariant — for the same trace,
// walker-path and fallback-path germline coords match.
// ──────────────────────────────────────────────────────────────────

use crate::address::ChoiceAddress;
use crate::ir::flag;
use crate::trace::{ChoiceValue, Trace};

/// Build a synthetic fixture that forces the V-side fallback:
/// assigns a V allele but DOES NOT add a V region to the sim.
/// Result: the walker never iterates V → `ref_ranges[V]` stays
/// `None` → the trim+end-loss fallback fires. A small non-V
/// region keeps the sequence non-empty so the AIRR builder's
/// pre-recombine early return doesn't short-circuit.
///
/// This is the canonical way to exercise the fallback branch in
/// isolation. In production with real assembly the V region is
/// always present and the walker handles end-loss directly; the
/// fallback is the defensive last-mile for partial-state outcomes.
fn fallback_forcing_v_fixture() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let v0 = cfg.v_pool.push(Allele {
        name: "v_fallback*01".into(),
        gene: "v_fallback".into(),
        seq: b"AAACCCGGG".to_vec(), // 9 bases
        segment: Segment::V,
        anchor: Some(6),
        functional_status: None,
        subregions: Vec::new(),
    });

    let mut sim = Simulation::new();
    // Push a few bytes carrying Np1 segment tag so the sequence is
    // non-empty but no V region exists. Walker iterates regions —
    // with no V region, ref_ranges[V] stays None.
    for _ in 0..3 {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
            b'A',
            Segment::Np1,
            flag::N_NUC,
        ));
        sim = next;
    }
    sim = sim
        .with_region_added(Region::new(
            Segment::Np1,
            NucHandle::new(0),
            NucHandle::new(3),
        ))
        .with_allele_assigned(Segment::V, AlleleInstance::new(v0));

    (cfg, sim)
}

/// J-side mirror of [`fallback_forcing_v_fixture`]: assigns a J
/// allele but adds only a non-J region, so `ref_ranges[J]` stays
/// `None` and the fallback fires.
fn fallback_forcing_j_fixture() -> (RefDataConfig, Simulation) {
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let j0 = cfg.j_pool.push(Allele {
        name: "j_fallback*01".into(),
        gene: "j_fallback".into(),
        seq: b"TTTAAA".to_vec(), // 6 bases
        segment: Segment::J,
        anchor: Some(0),
        functional_status: None,
        subregions: Vec::new(),
    });

    let mut sim = Simulation::new();
    for _ in 0..3 {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
            b'A',
            Segment::Np1,
            flag::N_NUC,
        ));
        sim = next;
    }
    sim = sim
        .with_region_added(Region::new(
            Segment::Np1,
            NucHandle::new(0),
            NucHandle::new(3),
        ))
        .with_allele_assigned(Segment::J, AlleleInstance::new(j0));

    (cfg, sim)
}

fn make_trace_with_end_loss_5(amount: i64) -> Trace {
    use crate::address::PrimeEnd;
    let mut trace = Trace::new();
    trace.record_choice(
        ChoiceAddress::CorruptEndLoss(PrimeEnd::Five),
        ChoiceValue::Int(amount),
    );
    trace
}

fn make_trace_with_end_loss_3(amount: i64) -> Trace {
    use crate::address::PrimeEnd;
    let mut trace = Trace::new();
    trace.record_choice(
        ChoiceAddress::CorruptEndLoss(PrimeEnd::Three),
        ChoiceValue::Int(amount),
    );
    trace
}

#[test]
fn fallback_v_germline_coords_unaffected_by_end_loss_match_walker() {
    // The walker encodes 5'-end-loss as leading `D` ops in V's
    // CIGAR and keeps `v_germline_start = v_trim_5`. The fallback
    // — which fires when the walker can't compute ref_ranges (no
    // V region in the sim, allele assigned for projection only) —
    // matches that convention by keying off `v_trim_5` /
    // `v_trim_3` and ignoring `end_loss_5_length` for coord math.
    // The provenance is carried by `end_loss_5_length` itself.
    let (cfg, sim) = fallback_forcing_v_fixture();
    let outcome = Outcome {
        revisions: vec![sim],
        pass_names: Vec::new(),
        trace: make_trace_with_end_loss_5(2),
        events: Vec::new(),
    };

    let rec = build_airr_record(&outcome, &cfg, "fallback_end_loss");

    // Recombination trim is independent provenance.
    assert_eq!(rec.v_trim_5, 0);
    assert_eq!(rec.v_trim_3, 0);
    // End-loss field surfaces from the trace (audit §6.1 fix).
    assert_eq!(rec.end_loss_5_length, 2);
    // Germline coords pin the AIRR convention: 0..allele_len,
    // unchanged by end-loss. The walker would emit '2D' at CIGAR
    // start; the fallback would emit the same coords with no
    // CIGAR ops (no walker columns), matching the trim-only
    // semantic.
    assert_eq!(rec.v_germline_start, Some(0));
    assert_eq!(rec.v_germline_end, Some(9));
}

#[test]
fn fallback_j_germline_coords_unaffected_by_end_loss_match_walker() {
    // J-side mirror. 3'-end-loss does NOT shift `j_germline_end`
    // in either the walker path (trailing 'D' ops) or the
    // fallback (trim-only). End-loss provenance lives in
    // `end_loss_3_length`.
    let (cfg, sim) = fallback_forcing_j_fixture();
    let outcome = Outcome {
        revisions: vec![sim],
        pass_names: Vec::new(),
        trace: make_trace_with_end_loss_3(2),
        events: Vec::new(),
    };

    let rec = build_airr_record(&outcome, &cfg, "fallback_j_endloss");

    assert_eq!(rec.j_trim_5, 0);
    assert_eq!(rec.j_trim_3, 0);
    assert_eq!(rec.end_loss_3_length, 2);
    assert_eq!(rec.j_germline_start, Some(0));
    // Full allele coverage (6 bases) — end-loss doesn't shift the
    // germline_end coord, only the AIRR projection's CIGAR.
    assert_eq!(rec.j_germline_end, Some(6));
}

#[test]
fn walker_path_keeps_v_germline_start_at_zero_under_5prime_end_loss() {
    // Sanity check on the WALKER path: a real post-end-loss pool
    // (germline_pos[2..9] after deleting V[0..2]) still produces
    // `v_germline_start = 0` because the walker fills the
    // deletion gap with 'D' ops. This is the AIRR convention the
    // fallback above mirrors.
    let mut cfg = RefDataConfig::empty(ChainType::Vj);
    let v0 = cfg.v_pool.push(Allele {
        name: "v_walker*01".into(),
        gene: "v_walker".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
        functional_status: None,
        subregions: Vec::new(),
    });
    let mut sim = Simulation::new();
    for (i, &b) in b"AAACCCGGG".iter().enumerate() {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    sim = sim
        .with_region_added(Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9)))
        .with_allele_assigned(Segment::V, AlleleInstance::new(v0))
        .with_indel_deleted(0)
        .with_indel_deleted(0);

    let outcome = Outcome {
        revisions: vec![sim],
        pass_names: Vec::new(),
        trace: make_trace_with_end_loss_5(2),
        events: Vec::new(),
    };

    let rec = build_airr_record(&outcome, &cfg, "walker_end_loss");

    // End-loss provenance still surfaces.
    assert_eq!(rec.end_loss_5_length, 2);
    // But germline coords describe the full allele coverage, not
    // the post-end-loss read coverage. The walker emits '2D' at
    // CIGAR start.
    assert_eq!(
        rec.v_germline_start,
        Some(0),
        "walker should keep v_germline_start at 0 (end-loss encoded \
         as leading 'D' ops, not a coord shift)"
    );
    assert_eq!(rec.v_germline_end, Some(9));
}

#[test]
fn fallback_and_walker_paths_agree_on_germline_coords_under_end_loss() {
    // Parity invariant: same end-loss recorded in the trace,
    // walker path and fallback path produce the same
    // `v_germline_start/end`. This is the AIRR-consistency check
    // — the audit's §6.3 concern was that they might diverge; in
    // practice the trim-only fallback math already matches the
    // walker's deletion-op encoding without modification.
    let mut cfg_walker = RefDataConfig::empty(ChainType::Vj);
    let v0 = cfg_walker.v_pool.push(Allele {
        name: "v_parity*01".into(),
        gene: "v_parity".into(),
        seq: b"AAACCCGGG".to_vec(),
        segment: Segment::V,
        anchor: Some(6),
        functional_status: None,
        subregions: Vec::new(),
    });
    let mut sim_walker = Simulation::new();
    for (i, &b) in b"AAACCCGGG".iter().enumerate() {
        let (next, _) =
            sim_walker.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim_walker = next;
    }
    sim_walker = sim_walker
        .with_region_added(Region::new(Segment::V, NucHandle::new(0), NucHandle::new(9)))
        .with_allele_assigned(Segment::V, AlleleInstance::new(v0))
        .with_indel_deleted(0)
        .with_indel_deleted(0);
    let walker_outcome = Outcome {
        revisions: vec![sim_walker],
        pass_names: Vec::new(),
        trace: make_trace_with_end_loss_5(2),
        events: Vec::new(),
    };

    let (cfg_fb, sim_fb) = fallback_forcing_v_fixture();
    let fallback_outcome = Outcome {
        revisions: vec![sim_fb],
        pass_names: Vec::new(),
        trace: make_trace_with_end_loss_5(2),
        events: Vec::new(),
    };

    let walker_rec = build_airr_record(&walker_outcome, &cfg_walker, "walker");
    let fallback_rec = build_airr_record(&fallback_outcome, &cfg_fb, "fallback");

    assert_eq!(
        fallback_rec.v_germline_start, walker_rec.v_germline_start,
        "fallback and walker paths must agree on v_germline_start"
    );
    assert_eq!(
        fallback_rec.v_germline_end, walker_rec.v_germline_end,
        "fallback and walker paths must agree on v_germline_end"
    );
    assert_eq!(fallback_rec.end_loss_5_length, walker_rec.end_loss_5_length);
}

// ──────────────────────────────────────────────────────────────────
// Slice E: AirrRecord.d_inverted is sourced from the final IR
// ──────────────────────────────────────────────────────────────────

/// Tiny VDJ fixture with a D allele assigned in default Forward
/// orientation. Used by the `d_inverted` tests; the test then
/// optionally flips the orientation via the persistent setter
/// before building the AIRR record.
fn vdj_with_d_assigned_fixture() -> (RefDataConfig, Simulation, AlleleId) {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    // Anchorless V/D/J to keep the fixture minimal — the AIRR
    // builder doesn't engage the anchor rule here.
    cfg.rules.v_anchor.required = false;
    cfg.rules.j_anchor.required = false;
    let _ = cfg.v_pool.push(Allele {
        name: "v*01".into(),
        gene: "v".into(),
        seq: b"AAACCC".to_vec(),
        segment: Segment::V,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let d_id = cfg.d_pool.push(Allele {
        name: "d*01".into(),
        gene: "d".into(),
        seq: b"ACGTTA".to_vec(),
        segment: Segment::D,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j*01".into(),
        gene: "j".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });

    // Minimal sim: pool empty, assignment present. The builder reads
    // `sim.assignments.get(Segment::D)` for `d_inverted` — pool /
    // region structure doesn't matter for the field's provenance.
    let sim = Simulation::new()
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::D, AlleleInstance::new(d_id))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
    (cfg, sim, d_id)
}

#[test]
fn d_inverted_is_false_for_default_forward_d_assignment() {
    let (cfg, sim, _) = vdj_with_d_assigned_fixture();
    let outcome = outcome_from_sim(sim);
    let rec = build_airr_record(&outcome, &cfg, "fwd-d");
    assert!(
        !rec.d_inverted,
        "default AlleleInstance::new(D) defaults to Forward; \
         AirrRecord.d_inverted must reflect that as false",
    );
}

#[test]
fn d_inverted_is_true_when_d_assignment_orientation_is_reverse() {
    use crate::assignment::SegmentOrientation;
    let (cfg, sim, _) = vdj_with_d_assigned_fixture();
    let sim = sim.with_allele_orientation(Segment::D, SegmentOrientation::ReverseComplement);
    let outcome = outcome_from_sim(sim);
    let rec = build_airr_record(&outcome, &cfg, "rc-d");
    assert!(
        rec.d_inverted,
        "with D orientation flipped to ReverseComplement, \
         AirrRecord.d_inverted must read true straight from the IR",
    );
}

#[test]
fn d_inverted_only_reads_d_orientation_not_v_or_j() {
    // Defensive: a future Slice that extends inversion to V/J must
    // make a conscious choice about projecting it. Slice E reads
    // ONLY the D-segment orientation. Setting V or J to RC must
    // NOT flip `d_inverted`.
    use crate::assignment::SegmentOrientation;
    let (cfg, sim, _) = vdj_with_d_assigned_fixture();
    let sim = sim
        .with_allele_orientation(Segment::V, SegmentOrientation::ReverseComplement)
        .with_allele_orientation(Segment::J, SegmentOrientation::ReverseComplement);
    let outcome = outcome_from_sim(sim);
    let rec = build_airr_record(&outcome, &cfg, "vj-rc-but-d-fwd");
    assert!(
        !rec.d_inverted,
        "d_inverted must remain false when only V or J orientation is reversed",
    );
}

#[test]
fn d_inverted_is_false_for_vj_chain_with_no_d_assignment() {
    // VJ chains have no D pool. The builder's `unwrap_or(false)`
    // handles the absent D assignment without panicking.
    let (cfg, sim) = anchor_record_fixture();
    let outcome = outcome_from_sim(sim);
    let rec = build_airr_record(&outcome, &cfg, "vj-no-d");
    assert!(
        !rec.d_inverted,
        "VJ chains lack a D assignment; d_inverted must default to false",
    );
}

#[test]
fn validator_flags_d_inverted_mismatch_when_record_disagrees_with_sim() {
    // The validator is a defence-in-depth check against forks of
    // the AIRR builder that forget to populate `d_inverted` or
    // post-processing that silently flips the bit. Forge a
    // mismatch by manually editing the field after the canonical
    // builder ran and confirm the issue surfaces.
    use crate::airr_record::{validate_airr_record, RecordValidationIssue};
    use crate::assignment::SegmentOrientation;
    let (cfg, sim, _) = vdj_with_d_assigned_fixture();
    let sim = sim.with_allele_orientation(Segment::D, SegmentOrientation::ReverseComplement);
    let outcome = outcome_from_sim(sim);
    let mut rec = build_airr_record(&outcome, &cfg, "tampered");
    // Builder set this to true; tamper to simulate a post-build edit.
    assert!(rec.d_inverted);
    rec.d_inverted = false;
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::DInvertedMismatch {
                reported: false,
                expected: true,
            }
        )),
        "validator must flag the tampered d_inverted; got {issues:?}",
    );
}

#[test]
fn validator_passes_clean_d_inverted_record() {
    // The clean record built by `build_airr_record` against the
    // same `Outcome` must never trip `DInvertedMismatch` — proves
    // the validator's expected-value calculation agrees with the
    // builder's source. Pin this in BOTH the forward and inverted
    // branches.
    use crate::airr_record::{validate_airr_record, RecordValidationIssue};
    use crate::assignment::SegmentOrientation;

    for orientation in [SegmentOrientation::Forward, SegmentOrientation::ReverseComplement] {
        let (cfg, sim, _) = vdj_with_d_assigned_fixture();
        let sim = sim.with_allele_orientation(Segment::D, orientation);
        let outcome = outcome_from_sim(sim);
        let rec = build_airr_record(&outcome, &cfg, "clean");
        let issues = validate_airr_record(&rec, &outcome, &cfg);
        assert!(
            !issues.iter().any(|i| matches!(
                i,
                RecordValidationIssue::DInvertedMismatch { .. }
            )),
            "validator must accept a clean record for orientation \
             {orientation:?}; got {issues:?}",
        );
    }
}

// ──────────────────────────────────────────────────────────────────
// Slice E: AirrRecord.receptor_revision_applied + original_v_call
// ──────────────────────────────────────────────────────────────────

/// Helper: build a VDJ fixture plus an Outcome that lets the
/// receptor-revision tests shape the exact projection state the
/// builder reads.
///
/// Receptor-revision provenance moved from the trace to the IR
/// (Bug D fix): `applied=true` is encoded as
/// `sim.assignments.v.receptor_revision_original_id = Some(...)`,
/// and `original_v_call` projects from that id.
/// `sample_allele_v_original_id` is the **pre-revision** V allele
/// id to install in the IR slot — distinct from the *current*
/// `assignments.v.allele_id`, which is whatever
/// `vdj_with_d_assigned_fixture` set up.
///
/// For non-revision cases (`applied=None` or `Some(false)`), the
/// IR slot is left untouched. The trace records that used to be
/// the source of truth are gone — the builder no longer reads
/// them — but we still leave the helper signature symmetrical
/// to make the callsite intent clear.
fn receptor_revision_outcome_with_trace(
    applied: Option<bool>,
    sample_allele_v_original_id: Option<u32>,
) -> (RefDataConfig, Outcome) {
    let (cfg, mut sim, _) = vdj_with_d_assigned_fixture();
    if matches!(applied, Some(true)) {
        let original_id = crate::refdata::AlleleId::new(
            sample_allele_v_original_id.unwrap_or(0),
        );
        let new_v = sim
            .assignments
            .get(Segment::V)
            .copied()
            .expect("vdj fixture must have V assigned")
            .with_receptor_revision_original_id(original_id);
        sim.assignments = sim.assignments.with_assigned(Segment::V, new_v);
    }
    let outcome = Outcome {
        revisions: vec![sim],
        pass_names: Vec::new(),
        trace: Trace::new(),
        events: Vec::new(),
    };
    (cfg, outcome)
}

#[test]
fn receptor_revision_defaults_false_and_empty_when_no_trace_record() {
    // No receptor-revision step ran; the trace lacks the address.
    // The builder must default applied=false and original_v_call="".
    let (cfg, outcome) = receptor_revision_outcome_with_trace(None, None);
    let rec = build_airr_record(&outcome, &cfg, "no-revision");
    assert!(!rec.receptor_revision_applied);
    assert_eq!(rec.original_v_call, "");
}

#[test]
fn receptor_revision_applied_false_keeps_original_v_call_empty() {
    // The pass ran but rolled applied=false. The original_v_call
    // must still be empty — the empty sentinel reads as "no
    // revision happened" regardless of whether `sample_allele.v`
    // is in the trace.
    let (cfg, outcome) = receptor_revision_outcome_with_trace(Some(false), Some(0));
    let rec = build_airr_record(&outcome, &cfg, "applied-false");
    assert!(!rec.receptor_revision_applied);
    assert_eq!(rec.original_v_call, "");
}

#[test]
fn receptor_revision_applied_true_populates_original_v_call_from_trace() {
    // The pass fired. `original_v_call` resolves from the trace's
    // `sample_allele.v` record (the V the recombine pass chose
    // before revision overwrote it).
    let (cfg, outcome) = receptor_revision_outcome_with_trace(Some(true), Some(0));
    let rec = build_airr_record(&outcome, &cfg, "applied-true");
    assert!(rec.receptor_revision_applied);
    assert_eq!(rec.original_v_call, "v*01");
}

#[test]
fn receptor_revision_applied_true_but_original_id_unresolvable_leaves_empty() {
    // Edge: IR provenance carries an `original_id` that does not
    // resolve in the current refdata (e.g. a refdata swap between
    // record and replay). The builder returns "" for
    // `original_v_call` and lets the validator surface
    // `OriginalVCallMismatch` if the consumer cares. Same shape
    // as the trace-era "no sample_allele.v record" case but
    // mapped to the IR-sourced model.
    let (cfg, mut sim, _) = vdj_with_d_assigned_fixture();
    let bogus_id = crate::refdata::AlleleId::new(9999); // out-of-range
    let new_v = sim
        .assignments
        .get(Segment::V)
        .copied()
        .expect("vdj fixture must have V assigned")
        .with_receptor_revision_original_id(bogus_id);
    sim.assignments = sim.assignments.with_assigned(Segment::V, new_v);
    let outcome = Outcome {
        revisions: vec![sim],
        pass_names: Vec::new(),
        trace: Trace::new(),
        events: Vec::new(),
    };
    let rec = build_airr_record(&outcome, &cfg, "applied-bogus-original");
    assert!(rec.receptor_revision_applied);
    assert_eq!(rec.original_v_call, "");
}

#[test]
fn validator_flags_receptor_revision_applied_mismatch() {
    use crate::airr_record::{validate_airr_record, RecordValidationIssue};
    let (cfg, outcome) = receptor_revision_outcome_with_trace(Some(true), Some(0));
    let mut rec = build_airr_record(&outcome, &cfg, "tampered");
    // Tamper the applied bool to the wrong value.
    assert!(rec.receptor_revision_applied);
    rec.receptor_revision_applied = false;
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReceptorRevisionAppliedMismatch {
                reported: false,
                expected: true,
            }
        )),
        "validator should flag a tampered receptor_revision_applied; got {issues:?}",
    );
}

#[test]
fn validator_flags_original_v_call_mismatch() {
    use crate::airr_record::{validate_airr_record, RecordValidationIssue};
    let (cfg, outcome) = receptor_revision_outcome_with_trace(Some(true), Some(0));
    let mut rec = build_airr_record(&outcome, &cfg, "tampered");
    assert_eq!(rec.original_v_call, "v*01");
    rec.original_v_call = "bogus*99".into();
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::OriginalVCallMismatch { expected, .. } if expected == "v*01"
        )),
        "validator should flag a tampered original_v_call; got {issues:?}",
    );
}

#[test]
fn validator_accepts_clean_no_revision_record() {
    use crate::airr_record::{validate_airr_record, RecordValidationIssue};
    let (cfg, outcome) = receptor_revision_outcome_with_trace(None, None);
    let rec = build_airr_record(&outcome, &cfg, "no-revision");
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        !issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReceptorRevisionAppliedMismatch { .. }
                | RecordValidationIssue::OriginalVCallMismatch { .. }
        )),
        "validator must accept a clean no-revision record; got {issues:?}",
    );
}

#[test]
fn validator_accepts_clean_revised_record() {
    use crate::airr_record::{validate_airr_record, RecordValidationIssue};
    let (cfg, outcome) = receptor_revision_outcome_with_trace(Some(true), Some(0));
    let rec = build_airr_record(&outcome, &cfg, "clean-revised");
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        !issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReceptorRevisionAppliedMismatch { .. }
                | RecordValidationIssue::OriginalVCallMismatch { .. }
        )),
        "validator must accept a clean revised record; got {issues:?}",
    );
}

// ──────────────────────────────────────────────────────────────────
// Slice A: paired-end schema + no-layout-default invariant
// ──────────────────────────────────────────────────────────────────

#[test]
fn paired_end_fields_default_to_empty_under_baseline_build() {
    // The builder must not populate any paired-end field on a
    // record that never went through a projection layer. Pinned
    // for every existing fixture so a future Slice A regression
    // (e.g. accidentally defaulting `read_layout = "single_end"`)
    // surfaces here.
    let (cfg, sim, _) = vdj_with_d_assigned_fixture();
    let outcome = outcome_from_sim(sim);
    let rec = build_airr_record(&outcome, &cfg, "baseline-paired-end");
    assert_eq!(rec.read_layout, "");
    assert_eq!(rec.r1_sequence, "");
    assert_eq!(rec.r2_sequence, "");
    assert_eq!(rec.r1_start, None);
    assert_eq!(rec.r1_end, None);
    assert_eq!(rec.r2_start, None);
    assert_eq!(rec.r2_end, None);
    assert_eq!(rec.insert_size, 0);
}

#[test]
fn validator_accepts_baseline_record_with_default_paired_end_fields() {
    use crate::airr_record::{validate_airr_record, RecordValidationIssue};
    let (cfg, sim, _) = vdj_with_d_assigned_fixture();
    let outcome = outcome_from_sim(sim);
    let rec = build_airr_record(&outcome, &cfg, "baseline-paired-end");
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        !issues
            .iter()
            .any(|i| matches!(i, RecordValidationIssue::PairedEndFieldWithoutLayout { .. })),
        "validator must accept a baseline record with default paired-end fields; \
         got {issues:?}",
    );
}

#[test]
fn validator_flags_each_paired_end_field_populated_without_layout() {
    // Tamper each of the seven non-layout fields one at a time;
    // every tamper must surface the matching
    // `PairedEndFieldWithoutLayout` variant. Pins the
    // per-field arm of the C6 check.
    use crate::airr_record::{
        validate_airr_record, PairedEndField, RecordValidationIssue,
    };
    let (cfg, sim, _) = vdj_with_d_assigned_fixture();
    let outcome = outcome_from_sim(sim);

    // Sanity: read_layout stays empty for every case below.
    let cases: Vec<(PairedEndField, Box<dyn Fn(&mut AirrRecord)>)> = vec![
        (
            PairedEndField::R1Sequence,
            Box::new(|r| r.r1_sequence = "ACGT".into()),
        ),
        (
            PairedEndField::R2Sequence,
            Box::new(|r| r.r2_sequence = "ACGT".into()),
        ),
        (
            PairedEndField::R1Start,
            Box::new(|r| r.r1_start = Some(0)),
        ),
        (
            PairedEndField::R1End,
            Box::new(|r| r.r1_end = Some(4)),
        ),
        (
            PairedEndField::R2Start,
            Box::new(|r| r.r2_start = Some(5)),
        ),
        (
            PairedEndField::R2End,
            Box::new(|r| r.r2_end = Some(9)),
        ),
        (
            PairedEndField::InsertSize,
            Box::new(|r| r.insert_size = 9),
        ),
    ];

    for (field, tamper) in cases {
        let mut rec = build_airr_record(&outcome, &cfg, "tampered");
        assert_eq!(rec.read_layout, "", "fixture must keep read_layout empty");
        tamper(&mut rec);
        let issues = validate_airr_record(&rec, &outcome, &cfg);
        assert!(
            issues.iter().any(|i| matches!(
                i,
                RecordValidationIssue::PairedEndFieldWithoutLayout { field: f } if *f == field
            )),
            "validator missed PairedEndFieldWithoutLayout for {field:?}; got {issues:?}",
        );
    }
}

#[test]
fn validator_no_longer_defers_geometry_after_slice_b_landed() {
    // Slice A originally deferred all geometry checks when
    // `read_layout` was set; Slice B replaces that deference with
    // active geometry checks. This regression test pins the
    // post-Slice-B contract: a `paired_end` layout with garbage
    // geometry surfaces at least one of the four reserved
    // variants. The Slice A "no-op when layout set" assertion is
    // intentionally obsolete now — pin the inverse so a future
    // refactor that accidentally restores the deference fails
    // here first.
    use crate::airr_record::{validate_airr_record, RecordValidationIssue};
    let (cfg, sim, _) = vdj_with_d_assigned_fixture();
    let outcome = outcome_from_sim(sim);
    let mut rec = build_airr_record(&outcome, &cfg, "post-slice-b");
    rec.read_layout = "paired_end".into();
    rec.r1_sequence = "garbage".into();
    rec.insert_size = -42;

    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReadWindowOutOfBounds { .. }
                | RecordValidationIssue::ReadSequenceMismatch { .. }
                | RecordValidationIssue::ReadInsertSizeMismatch { .. }
        )),
        "Slice B must surface geometry mismatches for `paired_end` records; \
         got {issues:?}",
    );
}

#[test]
fn paired_end_field_str_round_trip_pins_serialization_names() {
    // The PyO3 dict serializer uses `PairedEndField::as_str()` as
    // the `reported` and `details.field` payload. Renames break
    // downstream MCP / dashboard consumers; pin the eight strings
    // here so a refactor surfaces.
    use crate::airr_record::PairedEndField;
    assert_eq!(PairedEndField::ReadLayout.as_str(), "read_layout");
    assert_eq!(PairedEndField::R1Sequence.as_str(), "r1_sequence");
    assert_eq!(PairedEndField::R2Sequence.as_str(), "r2_sequence");
    assert_eq!(PairedEndField::R1Start.as_str(), "r1_start");
    assert_eq!(PairedEndField::R1End.as_str(), "r1_end");
    assert_eq!(PairedEndField::R2Start.as_str(), "r2_start");
    assert_eq!(PairedEndField::R2End.as_str(), "r2_end");
    assert_eq!(PairedEndField::InsertSize.as_str(), "insert_size");
}

// ──────────────────────────────────────────────────────────────────
// Slice B: projection kernel + per-layout geometry validator
// ──────────────────────────────────────────────────────────────────
//
// Tests live behind these helpers so the projection kernel itself
// can stay `pub(crate)` — no PyO3 surface yet, no public way to
// create paired-end records (per the Slice B scope). Synthetic
// records are constructed directly from the no-layout baseline and
// then patched with the kernel's output (or hand-tampered for the
// validator-rejection tests).

fn synthetic_record_with_sequence(s: &str) -> AirrRecord {
    let (cfg, sim, _) = vdj_with_d_assigned_fixture();
    let outcome = outcome_from_sim(sim);
    let mut rec = build_airr_record(&outcome, &cfg, "synthetic");
    rec.sequence = s.to_string();
    rec.sequence_length = s.len() as i64;
    rec
}

fn validate_synthetic(rec: &AirrRecord) -> Vec<crate::airr_record::RecordValidationIssue> {
    // Use the bundled VDJ fixture as the refdata/outcome backdrop;
    // the paired-end checks only read fields off the record so the
    // backdrop's content doesn't affect the result.
    let (cfg, sim, _) = vdj_with_d_assigned_fixture();
    let outcome = outcome_from_sim(sim);
    crate::airr_record::validate_airr_record(rec, &outcome, &cfg)
}

fn apply_projection(rec: &mut AirrRecord, p: &crate::airr_record::sequence::PairedEndProjection) {
    rec.read_layout = p.read_layout.clone();
    rec.r1_sequence = p.r1_sequence.clone();
    rec.r2_sequence = p.r2_sequence.clone();
    rec.r1_start = p.r1_start;
    rec.r1_end = p.r1_end;
    rec.r2_start = p.r2_start;
    rec.r2_end = p.r2_end;
    rec.insert_size = p.insert_size;
}

#[test]
fn projection_kernel_produces_documented_windows_for_simple_sequence() {
    // sequence = "AAAACCGGTT" (len 10), r1=3, r2=3, insert=8
    //   r1_start=0, r1_end=3  -> r1_sequence = "AAA"
    //   r2_start=5, r2_end=8  -> r2_inner    = "CGG"
    //                             revcomp     = "CCG"
    //   insert_size = 8
    use crate::airr_record::sequence::project_paired_end_layout;
    let rec = synthetic_record_with_sequence("AAAACCGGTT");
    let p = project_paired_end_layout(&rec, 3, 3, 8).expect("valid inputs");
    assert_eq!(p.read_layout, "paired_end");
    assert_eq!(p.r1_sequence, "AAA");
    assert_eq!(p.r2_sequence, "CCG");
    assert_eq!(p.r1_start, Some(0));
    assert_eq!(p.r1_end, Some(3));
    assert_eq!(p.r2_start, Some(5));
    assert_eq!(p.r2_end, Some(8));
    assert_eq!(p.insert_size, 8);
}

#[test]
fn projection_r2_is_reverse_complemented_not_just_reversed() {
    // sequence chosen so the r2 inner ("CTAG") differs under
    // reverse-only vs. reverse-complement:
    //   reverse("CTAG") = "GATC"
    //   revcomp("CTAG") = "CTAG"  (palindromic by RC)
    // Pick a non-palindromic inner instead so the two transforms
    // diverge unambiguously.
    use crate::airr_record::sequence::project_paired_end_layout;
    // sequence = "TTTTACCGGT" (len 10), r1=2, r2=4, insert=10
    //   r2_inner = sequence[6:10] = "CGGT"
    //   reverse("CGGT") = "TGGC"   (not what we want)
    //   revcomp("CGGT") = "ACCG"   (what r2_sequence must be)
    let rec = synthetic_record_with_sequence("TTTTACCGGT");
    let p = project_paired_end_layout(&rec, 2, 4, 10).expect("valid inputs");
    assert_eq!(p.r2_sequence, "ACCG");
    assert_ne!(p.r2_sequence, "TGGC", "R2 must be reverse-complement, not reverse-only");
}

#[test]
fn projection_only_reads_record_sequence_post_end_loss_or_revcomp() {
    // Pin §6.1: the kernel's coordinate space is `rec.sequence` —
    // whatever the pipeline finalised before paired-end runs.
    // Simulate "end-loss already happened" by handing the kernel
    // a short rec.sequence; the projection must not reach into
    // any longer or different upstream buffer.
    use crate::airr_record::sequence::project_paired_end_layout;
    let rec = synthetic_record_with_sequence("AAACCC");
    let p = project_paired_end_layout(&rec, 2, 2, 6).expect("valid inputs");
    assert_eq!(p.r1_sequence, "AA");
    assert_eq!(p.r2_sequence, "GGG"[..2].to_string());
    // r2_inner = "CC" -> revcomp = "GG"
    assert_eq!(p.r2_sequence, "GG");
    // Pool length / pre-end-loss bytes are NOT visible to the
    // kernel — only `rec.sequence` is. Pin this contract by
    // verifying that the projection's windows index strictly
    // into the short string.
    assert!(p.r2_end.unwrap() as i64 <= rec.sequence_length);
}

#[test]
fn projection_rejects_non_positive_r1_length() {
    use crate::airr_record::sequence::{project_paired_end_layout, PairedEndProjectionError};
    use crate::airr_record::PairedEndRead;
    let rec = synthetic_record_with_sequence("AAAACCGGTT");
    let err = project_paired_end_layout(&rec, 0, 3, 8).unwrap_err();
    assert_eq!(
        err,
        PairedEndProjectionError::NonPositiveLength { side: PairedEndRead::R1, length: 0 }
    );
    let err = project_paired_end_layout(&rec, -1, 3, 8).unwrap_err();
    assert!(matches!(
        err,
        PairedEndProjectionError::NonPositiveLength { side: PairedEndRead::R1, length: -1 }
    ));
}

#[test]
fn projection_rejects_non_positive_r2_length() {
    use crate::airr_record::sequence::{project_paired_end_layout, PairedEndProjectionError};
    use crate::airr_record::PairedEndRead;
    let rec = synthetic_record_with_sequence("AAAACCGGTT");
    let err = project_paired_end_layout(&rec, 3, 0, 8).unwrap_err();
    assert_eq!(
        err,
        PairedEndProjectionError::NonPositiveLength { side: PairedEndRead::R2, length: 0 }
    );
}

#[test]
fn projection_rejects_insert_size_outside_sequence_bounds() {
    use crate::airr_record::sequence::{project_paired_end_layout, PairedEndProjectionError};
    let rec = synthetic_record_with_sequence("AAAACCGGTT"); // len 10
    // Past sequence end.
    let err = project_paired_end_layout(&rec, 3, 3, 11).unwrap_err();
    assert_eq!(
        err,
        PairedEndProjectionError::InsertSizeOutOfBounds {
            insert_size: 11,
            sequence_length: 10
        }
    );
    // Negative.
    let err = project_paired_end_layout(&rec, 3, 3, -1).unwrap_err();
    assert!(matches!(
        err,
        PairedEndProjectionError::InsertSizeOutOfBounds {
            insert_size: -1,
            ..
        }
    ));
}

#[test]
fn projection_rejects_read_length_exceeding_insert_size() {
    use crate::airr_record::sequence::{project_paired_end_layout, PairedEndProjectionError};
    use crate::airr_record::PairedEndRead;
    let rec = synthetic_record_with_sequence("AAAACCGGTT");
    // r1_length > insert_size
    let err = project_paired_end_layout(&rec, 9, 3, 8).unwrap_err();
    assert!(matches!(
        err,
        PairedEndProjectionError::ReadExceedsInsert {
            side: PairedEndRead::R1,
            length: 9,
            insert_size: 8
        }
    ));
    // r2_length > insert_size
    let err = project_paired_end_layout(&rec, 3, 9, 8).unwrap_err();
    assert!(matches!(
        err,
        PairedEndProjectionError::ReadExceedsInsert {
            side: PairedEndRead::R2,
            length: 9,
            insert_size: 8
        }
    ));
}

#[test]
fn projection_allows_r1_plus_r2_overlapping_within_insert() {
    // Audit §6.2: r1+r2 > insert_size is allowed (reads overlap
    // in the middle). The kernel succeeds; the validator sees two
    // overlapping windows and accepts.
    use crate::airr_record::sequence::project_paired_end_layout;
    let rec = synthetic_record_with_sequence("AAAACCGGTT");
    let p = project_paired_end_layout(&rec, 6, 6, 10).expect("overlap allowed");
    assert_eq!(p.r1_start, Some(0));
    assert_eq!(p.r1_end, Some(6));
    assert_eq!(p.r2_start, Some(4));
    assert_eq!(p.r2_end, Some(10));
}

#[test]
fn validator_accepts_record_populated_via_projection_kernel() {
    use crate::airr_record::sequence::project_paired_end_layout;
    let mut rec = synthetic_record_with_sequence("AAAACCGGTT");
    let p = project_paired_end_layout(&rec, 3, 3, 8).unwrap();
    apply_projection(&mut rec, &p);
    let issues = validate_synthetic(&rec);
    let paired_end_issues: Vec<_> = issues
        .iter()
        .filter(|i| {
            matches!(
                i,
                crate::airr_record::RecordValidationIssue::PairedEndFieldWithoutLayout { .. }
                    | crate::airr_record::RecordValidationIssue::ReadWindowOutOfBounds { .. }
                    | crate::airr_record::RecordValidationIssue::ReadSequenceMismatch { .. }
                    | crate::airr_record::RecordValidationIssue::ReadInsertSizeMismatch { .. }
                    | crate::airr_record::RecordValidationIssue::ReadLayoutMismatch { .. }
            )
        })
        .collect();
    assert!(
        paired_end_issues.is_empty(),
        "validator must accept a clean kernel-derived record; got {paired_end_issues:?}",
    );
}

#[test]
fn validator_flags_read_window_out_of_bounds_for_r1() {
    use crate::airr_record::{PairedEndRead, RecordValidationIssue};
    use crate::airr_record::sequence::project_paired_end_layout;
    let mut rec = synthetic_record_with_sequence("AAAACCGGTT");
    let p = project_paired_end_layout(&rec, 3, 3, 8).unwrap();
    apply_projection(&mut rec, &p);
    // Tamper R1 window to overshoot the sequence.
    rec.r1_end = Some(99);
    let issues = validate_synthetic(&rec);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReadWindowOutOfBounds { side: PairedEndRead::R1, .. }
        )),
        "expected ReadWindowOutOfBounds for R1; got {issues:?}",
    );
}

#[test]
fn validator_flags_read_window_out_of_bounds_when_coord_missing() {
    // Missing coord (None) while read_layout = "paired_end" must
    // surface as out-of-bounds with the sentinel -1 in the
    // structured payload.
    use crate::airr_record::{PairedEndRead, RecordValidationIssue};
    use crate::airr_record::sequence::project_paired_end_layout;
    let mut rec = synthetic_record_with_sequence("AAAACCGGTT");
    let p = project_paired_end_layout(&rec, 3, 3, 8).unwrap();
    apply_projection(&mut rec, &p);
    rec.r2_start = None;
    let issues = validate_synthetic(&rec);
    let oob = issues.iter().find(|i| matches!(
        i,
        RecordValidationIssue::ReadWindowOutOfBounds { side: PairedEndRead::R2, .. }
    ));
    let issue = oob.expect("expected ReadWindowOutOfBounds for R2 missing-coord");
    match issue {
        RecordValidationIssue::ReadWindowOutOfBounds { start, end, .. } => {
            assert_eq!(*start, -1, "missing start must surface as sentinel -1");
            assert_eq!(*end, 8, "end stays at the populated value");
        }
        _ => unreachable!(),
    }
}

#[test]
fn validator_flags_read_sequence_mismatch_for_r1() {
    use crate::airr_record::{PairedEndRead, RecordValidationIssue};
    use crate::airr_record::sequence::project_paired_end_layout;
    let mut rec = synthetic_record_with_sequence("AAAACCGGTT");
    let p = project_paired_end_layout(&rec, 3, 3, 8).unwrap();
    apply_projection(&mut rec, &p);
    // Bytes diverge from sequence[r1_start..r1_end] = "AAA".
    rec.r1_sequence = "GGG".to_string();
    let issues = validate_synthetic(&rec);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReadSequenceMismatch {
                side: PairedEndRead::R1, expected, ..
            } if expected == "AAA"
        )),
        "expected ReadSequenceMismatch for R1 with expected=\"AAA\"; got {issues:?}",
    );
}

#[test]
fn validator_flags_read_sequence_mismatch_for_r2() {
    use crate::airr_record::{PairedEndRead, RecordValidationIssue};
    use crate::airr_record::sequence::project_paired_end_layout;
    let mut rec = synthetic_record_with_sequence("AAAACCGGTT");
    let p = project_paired_end_layout(&rec, 3, 3, 8).unwrap();
    apply_projection(&mut rec, &p);
    // Pretend the user supplied reverse-only instead of
    // reverse-complement for R2. The validator must catch it.
    let r2_inner = &rec.sequence[5..8]; // "CGG"
    let reversed_only: String = r2_inner.chars().rev().collect(); // "GGC"
    rec.r2_sequence = reversed_only;
    let issues = validate_synthetic(&rec);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReadSequenceMismatch {
                side: PairedEndRead::R2, ..
            }
        )),
        "expected ReadSequenceMismatch for R2 (reverse-only vs. revcomp); got {issues:?}",
    );
}

#[test]
fn validator_flags_read_insert_size_mismatch() {
    use crate::airr_record::RecordValidationIssue;
    use crate::airr_record::sequence::project_paired_end_layout;
    let mut rec = synthetic_record_with_sequence("AAAACCGGTT");
    let p = project_paired_end_layout(&rec, 3, 3, 8).unwrap();
    apply_projection(&mut rec, &p);
    // insert_size must equal r2_end (=8) per §8. Bump it.
    rec.insert_size = 7;
    let issues = validate_synthetic(&rec);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReadInsertSizeMismatch { reported: 7, expected: 8 }
        )),
        "expected ReadInsertSizeMismatch reported=7 expected=8; got {issues:?}",
    );
}

#[test]
fn validator_flags_read_layout_mismatch_for_unknown_value() {
    use crate::airr_record::RecordValidationIssue;
    let mut rec = synthetic_record_with_sequence("AAAACCGGTT");
    rec.read_layout = "pair_end".to_string(); // typo
    let issues = validate_synthetic(&rec);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReadLayoutMismatch { reported, .. } if reported == "pair_end"
        )),
        "expected ReadLayoutMismatch for unknown layout; got {issues:?}",
    );
}

#[test]
fn validator_accepts_single_end_layout_as_reserved_no_op() {
    // "single_end" is documented as reserved (§2.2); Slice B
    // doesn't validate geometry for it. Pin the reserved status
    // so a future implementer doesn't accidentally land
    // single-end checks before the dedicated audit slice.
    use crate::airr_record::RecordValidationIssue;
    let mut rec = synthetic_record_with_sequence("AAAACCGGTT");
    rec.read_layout = "single_end".to_string();
    let issues = validate_synthetic(&rec);
    assert!(
        !issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReadLayoutMismatch { .. }
                | RecordValidationIssue::PairedEndFieldWithoutLayout { .. }
                | RecordValidationIssue::ReadWindowOutOfBounds { .. }
                | RecordValidationIssue::ReadSequenceMismatch { .. }
                | RecordValidationIssue::ReadInsertSizeMismatch { .. }
        )),
        "single_end layout must be a no-op in Slice B; got {issues:?}",
    );
}

#[test]
fn validator_accepts_baseline_no_layout_unchanged_from_slice_a() {
    // Regression pin: the Slice B dispatch must not change the
    // Slice A no-layout invariant. A clean baseline record stays
    // green.
    use crate::airr_record::RecordValidationIssue;
    let rec = synthetic_record_with_sequence("AAAACCGGTT");
    assert_eq!(rec.read_layout, "");
    let issues = validate_synthetic(&rec);
    assert!(
        !issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::PairedEndFieldWithoutLayout { .. }
                | RecordValidationIssue::ReadWindowOutOfBounds { .. }
                | RecordValidationIssue::ReadSequenceMismatch { .. }
                | RecordValidationIssue::ReadInsertSizeMismatch { .. }
                | RecordValidationIssue::ReadLayoutMismatch { .. }
        )),
        "Slice B regressed the Slice A baseline-clean invariant; got {issues:?}",
    );
}

// ──────────────────────────────────────────────────────────────────
// Slice C: trace-driven builder integration
//
// End-to-end pinning: `PairedEndSamplingPass` records three Ints on
// the trace; `build_airr_record` reads those records back and
// applies the projection kernel. Replay reproduces the eight AIRR
// fields bit-for-bit. Tests run through the real `PassRuntime`
// harness rather than the synthetic-record helpers above so the
// trace→builder roundtrip is exercised at the same depth a real
// Slice D DSL caller would.
// ──────────────────────────────────────────────────────────────────

fn vdj_recombine_sim() -> (RefDataConfig, Simulation) {
    // Tiny VDJ fixture with a fully assembled molecule so the AIRR
    // builder lands a non-empty `rec.sequence` (long enough that
    // 5/5/12 paired-end geometry fits).
    use crate::ir::{NucHandle, Nucleotide, Region, Segment};
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    cfg.rules.v_anchor.required = false;
    cfg.rules.j_anchor.required = false;
    let _ = cfg.v_pool.push(Allele {
        name: "v*01".into(),
        gene: "v".into(),
        seq: b"AAAAAA".to_vec(),
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
        seq: b"CCCCCC".to_vec(),
        segment: Segment::J,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });

    let mut sim = Simulation::new();
    for (i, b) in b"AAAAAAGGGCCCCCC".iter().enumerate() {
        let seg = if i < 6 {
            Segment::V
        } else if i < 9 {
            Segment::D
        } else {
            Segment::J
        };
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(*b, i as u16, seg));
        sim = next;
    }
    sim = sim
        .with_region_added(Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6)))
        .with_region_added(Region::new(Segment::D, NucHandle::new(6), NucHandle::new(9)))
        .with_region_added(Region::new(Segment::J, NucHandle::new(9), NucHandle::new(15)))
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

    (cfg, sim)
}

fn run_with_paired_end_pass(
    r1: i64,
    r2: i64,
    insert: i64,
) -> (RefDataConfig, crate::pass::Outcome) {
    use crate::dist::EmpiricalLengthDist;
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;
    use crate::passes::paired_end::{PairedEndLayoutSpec, PairedEndSamplingPass};

    let (cfg, sim) = vdj_recombine_sim();
    let mut plan = PassPlan::new();
    plan.push(Box::new(PairedEndSamplingPass::new(
        PairedEndLayoutSpec::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(r1, 1.0)])),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(r2, 1.0)])),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(insert, 1.0)])),
        ),
    )));
    let outcome = PassRuntime::execute(&plan, sim, 0);
    (cfg, outcome)
}

#[test]
fn slice_c_fresh_run_populates_airr_record_via_projection_kernel() {
    // r1=5, r2=5, insert=12, sequence length 15. Kernel returns
    // r1=sequence[0..5], r2=revcomp(sequence[7..12]).
    let (cfg, outcome) = run_with_paired_end_pass(5, 5, 12);
    let rec = build_airr_record(&outcome, &cfg, "fresh");
    assert_eq!(rec.read_layout, "paired_end");
    assert_eq!(rec.r1_start, Some(0));
    assert_eq!(rec.r1_end, Some(5));
    assert_eq!(rec.r2_start, Some(7));
    assert_eq!(rec.r2_end, Some(12));
    assert_eq!(rec.insert_size, 12);
    // sequence = "AAAAAAGGGCCCCCC"
    //   r1 = sequence[0..5]  = "AAAAA"
    //   r2_inner = sequence[7..12] = "GGCCC"
    //   revcomp("GGCCC") = "GGGCC"
    assert_eq!(rec.r1_sequence, "AAAAA");
    assert_eq!(rec.r2_sequence, "GGGCC");
}

#[test]
fn slice_c_validator_accepts_trace_populated_record() {
    use crate::airr_record::{validate_airr_record, RecordValidationIssue};
    let (cfg, outcome) = run_with_paired_end_pass(5, 5, 12);
    let rec = build_airr_record(&outcome, &cfg, "trace-populated");
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        !issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::PairedEndFieldWithoutLayout { .. }
                | RecordValidationIssue::ReadWindowOutOfBounds { .. }
                | RecordValidationIssue::ReadSequenceMismatch { .. }
                | RecordValidationIssue::ReadInsertSizeMismatch { .. }
                | RecordValidationIssue::ReadLayoutMismatch { .. }
        )),
        "validator must accept a clean trace-populated paired-end record; got {issues:?}",
    );
}

#[test]
fn slice_c_builder_is_pure_function_of_trace_so_replay_reproduces_fields() {
    // Pin the replay-bit-for-bit invariant by demonstrating the
    // composition: (a) the pass-level replay test
    // (`replay_consumes_recorded_values_in_order` in
    // `passes::paired_end::tests`) proves the three Ints round-
    // trip; (b) this test proves the builder is a pure function
    // of the trace — same trace + same `Simulation` ⇒ identical
    // eight AIRR fields. Composition gives bit-for-bit replay
    // through the full builder path.
    use crate::address::ChoiceAddress;
    use crate::pass::Outcome;
    use crate::trace::{ChoiceValue, Trace};

    let (cfg, base_sim) = vdj_recombine_sim();
    let mut trace_a = Trace::new();
    trace_a.record_choice(ChoiceAddress::PairedEndR1Length, ChoiceValue::Int(5));
    trace_a.record_choice(ChoiceAddress::PairedEndR2Length, ChoiceValue::Int(5));
    trace_a.record_choice(ChoiceAddress::PairedEndInsertSize, ChoiceValue::Int(12));
    let outcome_a = Outcome {
        revisions: vec![base_sim.clone()],
        pass_names: Vec::new(),
        trace: trace_a.clone(),
        events: Vec::new(),
    };
    let outcome_b = Outcome {
        revisions: vec![base_sim],
        pass_names: Vec::new(),
        trace: trace_a, // identical trace → identical record
        events: Vec::new(),
    };
    let rec_a = build_airr_record(&outcome_a, &cfg, "a");
    let rec_b = build_airr_record(&outcome_b, &cfg, "b");

    assert_eq!(rec_a.read_layout, rec_b.read_layout);
    assert_eq!(rec_a.r1_sequence, rec_b.r1_sequence);
    assert_eq!(rec_a.r2_sequence, rec_b.r2_sequence);
    assert_eq!(rec_a.r1_start, rec_b.r1_start);
    assert_eq!(rec_a.r1_end, rec_b.r1_end);
    assert_eq!(rec_a.r2_start, rec_b.r2_start);
    assert_eq!(rec_a.r2_end, rec_b.r2_end);
    assert_eq!(rec_a.insert_size, rec_b.insert_size);
    // Sanity: the projection did actually fire (otherwise the
    // equality would be trivially true via defaults).
    assert_eq!(rec_a.read_layout, "paired_end");
    assert_eq!(rec_a.r1_sequence, "AAAAA");
}

#[test]
fn slice_c_builder_tags_layout_when_kernel_rejects_insert_past_sequence() {
    // Insert size 30 against a 15-byte sequence — kernel
    // returns InsertSizeOutOfBounds. Builder tags
    // read_layout="paired_end" but leaves coords None; validator
    // surfaces ReadWindowOutOfBounds for both sides.
    use crate::airr_record::{validate_airr_record, PairedEndRead, RecordValidationIssue};
    use crate::dist::EmpiricalLengthDist;
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;
    use crate::passes::paired_end::{PairedEndLayoutSpec, PairedEndSamplingPass};

    // The PairedEndSamplingPass's relationship check requires
    // r1_length <= insert_size, so we use r1=r2=15 with
    // insert_size=15 to clear the pass-level check; then we
    // bypass the kernel's seq-length check by hand-crafting a
    // shorter sequence below.
    //
    // Instead: have the sampling pass record an insert larger
    // than the sequence_length. The sampling pass doesn't know
    // the seq length, so this is allowed at the pass layer.
    // The kernel catches it at AIRR build.
    let (cfg, sim) = vdj_recombine_sim();
    let mut plan = PassPlan::new();
    plan.push(Box::new(PairedEndSamplingPass::new(
        PairedEndLayoutSpec::new(
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(5, 1.0)])),
            Box::new(EmpiricalLengthDist::from_pairs(vec![(99, 1.0)])),
        ),
    )));
    let outcome = PassRuntime::execute(&plan, sim, 0);
    let rec = build_airr_record(&outcome, &cfg, "kernel-rejected");
    // Layout was set but coords stayed None — kernel rejected.
    assert_eq!(rec.read_layout, "paired_end");
    assert_eq!(rec.r1_start, None);
    assert_eq!(rec.r2_start, None);

    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::ReadWindowOutOfBounds {
                side: PairedEndRead::R1, ..
            }
        )),
        "kernel-rejected record must surface ReadWindowOutOfBounds for R1; got {issues:?}",
    );
}

#[test]
fn slice_c_baseline_run_without_pass_keeps_airr_fields_at_defaults() {
    use crate::airr_record::{validate_airr_record, RecordValidationIssue};
    use crate::pass::testing::PassRuntime;
    use crate::pass::PassPlan;
    // Empty plan; no PairedEndSamplingPass.
    let (cfg, sim) = vdj_recombine_sim();
    let plan = PassPlan::new();
    let outcome = PassRuntime::execute(&plan, sim, 0);
    let rec = build_airr_record(&outcome, &cfg, "baseline");
    assert_eq!(rec.read_layout, "");
    assert_eq!(rec.r1_sequence, "");
    assert_eq!(rec.r2_sequence, "");
    assert_eq!(rec.r1_start, None);
    assert_eq!(rec.r1_end, None);
    assert_eq!(rec.r2_start, None);
    assert_eq!(rec.r2_end, None);
    assert_eq!(rec.insert_size, 0);
    // And the trace carries no paired_end.* records.
    for choice in outcome.trace.choices() {
        assert!(
            !choice.address.starts_with("paired_end."),
            "baseline outcome emitted paired_end.* address {:?}",
            choice.address,
        );
    }
    // Validator runs clean (regression check matching the Slice
    // A behaviour now end-to-end through the projection).
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        !issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::PairedEndFieldWithoutLayout { .. }
                | RecordValidationIssue::ReadWindowOutOfBounds { .. }
                | RecordValidationIssue::ReadSequenceMismatch { .. }
                | RecordValidationIssue::ReadInsertSizeMismatch { .. }
                | RecordValidationIssue::ReadLayoutMismatch { .. }
        )),
        "baseline outcome surfaced paired-end issues: {issues:?}",
    );
}


// ──────────────────────────────────────────────────────────────────
// V-subregion mutation counters — validator mismatch detection
// (Slice — V-Subregion Mutation Counters)
// ──────────────────────────────────────────────────────────────────

/// Build a `(RefDataConfig, Outcome)` where the V allele has a
/// six-base sequence with two equal-sized IMGT subregion intervals
/// (FWR1 = [0, 3), CDR1 = [3, 6)). The outcome carries a single
/// `mutate.s5f` `EventRecord` containing one `BaseChanged` event
/// in V at `germline_pos = 4` — i.e. inside CDR1.
fn v_subregion_counter_outcome_fixture_with_pos(
    germline_pos: Option<u16>,
) -> (RefDataConfig, Outcome) {
    use crate::event::{EventRecord, StateSummary, TraceSpan};
    use crate::ir::SimulationEvent;
    use crate::refdata::{VSubregion, VSubregionLabel};

    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    cfg.rules.v_anchor.required = false;
    cfg.rules.j_anchor.required = false;
    let _ = cfg.v_pool.push(Allele {
        name: "v*01".into(),
        gene: "v".into(),
        seq: b"AAACCC".to_vec(),
        segment: Segment::V,
        anchor: None,
        functional_status: None,
        subregions: vec![
            VSubregion {
                label: VSubregionLabel::Fwr1,
                start: 0,
                end: 3,
            },
            VSubregion {
                label: VSubregionLabel::Cdr1,
                start: 3,
                end: 6,
            },
        ],
    });
    let _ = cfg.d_pool.push(Allele {
        name: "d*01".into(),
        gene: "d".into(),
        seq: b"ACGTTA".to_vec(),
        segment: Segment::D,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j*01".into(),
        gene: "j".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });

    // Push a six-base V pool so the AIRR builder doesn't early-return
    // on `rec.sequence.is_empty()`. The bases mirror the V allele's
    // `seq` ("AAACCC") with germline_pos = 0..6 covering the two
    // subregion intervals (FWR1 = [0, 3), CDR1 = [3, 6)).
    let mut sim = Simulation::new();
    for (i, b) in b"AAACCC".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(
            crate::ir::Nucleotide::germline(*b, i as u16, Segment::V),
        );
        sim = next;
    }
    let v_region = crate::ir::Region::new(
        Segment::V,
        crate::ir::NucHandle::new(0),
        crate::ir::NucHandle::new(6),
    );
    let sim = sim
        .with_region_added(v_region)
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

    let state_pre = StateSummary::from_simulation(&sim);
    let state_post = state_pre.clone();
    let event_record = EventRecord::pass_committed(
        0,
        crate::address::MUTATE_S5F,
        Vec::new(),
        TraceSpan::new(0, 0),
        state_pre,
        state_post,
        vec![SimulationEvent::BaseChanged {
            handle: crate::ir::NucHandle::new(0),
            old_base: b'A',
            new_base: b'C',
            segment: Segment::V,
            germline_pos,
        }],
    );

    let outcome = Outcome {
        revisions: vec![sim],
        pass_names: vec![crate::address::MUTATE_S5F.to_string()],
        trace: Trace::new(),
        events: vec![event_record],
    };
    (cfg, outcome)
}

fn v_subregion_counter_outcome_fixture() -> (RefDataConfig, Outcome) {
    // Default: germline_pos = 4 → inside CDR1 = [3, 6).
    v_subregion_counter_outcome_fixture_with_pos(Some(4))
}

#[test]
fn v_subregion_counters_aggregate_to_correct_bucket() {
    use crate::airr_record::build_airr_record;
    let (cfg, outcome) = v_subregion_counter_outcome_fixture();
    let rec = build_airr_record(&outcome, &cfg, "single-cdr1");
    // The single V `BaseChanged` event at germline_pos=4 must
    // route to CDR1 (interval [3, 6)).
    assert_eq!(rec.n_v_mutations, 1);
    assert_eq!(rec.n_fwr1_mutations, 0);
    assert_eq!(rec.n_cdr1_mutations, 1);
    assert_eq!(rec.n_fwr2_mutations, 0);
    assert_eq!(rec.n_cdr2_mutations, 0);
    assert_eq!(rec.n_fwr3_mutations, 0);
    assert_eq!(rec.n_v_unannotated_mutations, 0);
}

#[test]
fn v_subregion_counters_route_to_unannotated_when_germline_pos_is_none() {
    use crate::airr_record::build_airr_record;
    let (cfg, outcome) = v_subregion_counter_outcome_fixture_with_pos(None);
    let rec = build_airr_record(&outcome, &cfg, "indel-base");
    assert_eq!(rec.n_v_mutations, 1);
    assert_eq!(rec.n_v_unannotated_mutations, 1);
    assert_eq!(rec.n_cdr1_mutations, 0);
}

#[test]
fn validator_flags_tampered_n_cdr1_mutations() {
    use crate::airr_record::{build_airr_record, validate_airr_record, RecordValidationIssue};
    let (cfg, outcome) = v_subregion_counter_outcome_fixture();
    let mut rec = build_airr_record(&outcome, &cfg, "tampered-cdr1");
    // Engine projection assigned the V SHM event to CDR1.
    // Tamper: shove it into FWR2 instead.
    assert_eq!(rec.n_cdr1_mutations, 1);
    rec.n_cdr1_mutations = 0;
    rec.n_fwr2_mutations = 1;
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::NCdr1MutationsMismatch {
                reported: 0,
                event_count: 1,
            }
        )),
        "validator should flag tampered n_cdr1_mutations; got {issues:?}",
    );
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::NFwr2MutationsMismatch {
                reported: 1,
                event_count: 0,
            }
        )),
        "validator should flag the symmetric FWR2 mismatch; got {issues:?}",
    );
}

#[test]
fn validator_flags_v_subregion_partition_sum_mismatch() {
    use crate::airr_record::{build_airr_record, validate_airr_record, RecordValidationIssue};
    let (cfg, outcome) = v_subregion_counter_outcome_fixture();
    let mut rec = build_airr_record(&outcome, &cfg, "tampered-partition");
    // Bump a subregion counter without touching n_v_mutations
    // — the partition invariant breaks.
    rec.n_fwr3_mutations += 7;
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::VSubregionMutationCountSumMismatch {
                reported_v_total,
                sum_of_subregion_buckets,
            } if *reported_v_total == 1 && *sum_of_subregion_buckets == 8
        )),
        "validator should flag V-subregion partition mismatch; got {issues:?}",
    );
}

fn p_addition_outcome_fixture() -> (RefDataConfig, Outcome) {
    use crate::address::{PEnd, P_ADDITION_V_3};
    use crate::event::{EventRecord, StateSummary, TraceSpan};
    use crate::ir::{Region, SimulationEvent};

    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    cfg.rules.v_anchor.required = false;
    cfg.rules.j_anchor.required = false;
    let _ = cfg.v_pool.push(Allele {
        name: "v*01".into(),
        gene: "v".into(),
        seq: b"AAACCC".to_vec(),
        segment: Segment::V,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.d_pool.push(Allele {
        name: "d*01".into(),
        gene: "d".into(),
        seq: b"ACGTTA".to_vec(),
        segment: Segment::D,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });
    let _ = cfg.j_pool.push(Allele {
        name: "j*01".into(),
        gene: "j".into(),
        seq: b"TTTAAA".to_vec(),
        segment: Segment::J,
        anchor: None,
        functional_status: None,
        subregions: Vec::new(),
    });

    let mut sim = Simulation::new();
    for (i, b) in b"AAACCC".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(
            crate::ir::Nucleotide::germline(*b, i as u16, Segment::V),
        );
        sim = next;
    }
    let v_region = Region::new(
        Segment::V,
        crate::ir::NucHandle::new(0),
        crate::ir::NucHandle::new(6),
    );
    let sim = sim
        .with_region_added(v_region)
        .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)))
        .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

    // Emit a P-region event of length 3 at V_3 — the validator's
    // recompute should sum to 3.
    let state_pre = StateSummary::from_simulation(&sim);
    let state_post = state_pre.clone();
    let p_v_3_region = Region::new(
        Segment::V,
        crate::ir::NucHandle::new(6),
        crate::ir::NucHandle::new(9),
    );
    let event_record = EventRecord::pass_committed(
        0,
        P_ADDITION_V_3,
        Vec::new(),
        TraceSpan::new(0, 0),
        state_pre,
        state_post,
        vec![SimulationEvent::PRegionAdded {
            end: PEnd::V3,
            region: p_v_3_region,
        }],
    );
    let outcome = Outcome {
        revisions: vec![sim],
        pass_names: vec![P_ADDITION_V_3.to_string()],
        trace: Trace::new(),
        events: vec![event_record],
    };
    (cfg, outcome)
}

#[test]
fn validator_flags_p_length_mismatch_on_tampered_field() {
    use crate::address::PEnd;
    use crate::airr_record::{build_airr_record, validate_airr_record, RecordValidationIssue};
    let (cfg, outcome) = p_addition_outcome_fixture();
    let mut rec = build_airr_record(&outcome, &cfg, "tampered-p-v-3");
    assert_eq!(rec.p_v_3_length, 3, "fixture builder should compute 3");
    // Tamper: claim 99 P_V3 bytes but the ledger says 3.
    rec.p_v_3_length = 99;
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::PLengthMismatch {
                end: PEnd::V3,
                reported: 99,
                event_count: 3,
            }
        )),
        "validator should flag PLengthMismatch on tampered p_v_3_length; got {issues:?}",
    );
}

#[test]
fn validator_accepts_clean_p_length_projection() {
    use crate::airr_record::{build_airr_record, validate_airr_record, RecordValidationIssue};
    let (cfg, outcome) = p_addition_outcome_fixture();
    let rec = build_airr_record(&outcome, &cfg, "clean-p");
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        !issues
            .iter()
            .any(|i| matches!(i, RecordValidationIssue::PLengthMismatch { .. })),
        "clean P-addition projection must validate; got {issues:?}",
    );
}

#[test]
fn validator_accepts_clean_v_subregion_projection() {
    use crate::airr_record::{build_airr_record, validate_airr_record, RecordValidationIssue};
    let (cfg, outcome) = v_subregion_counter_outcome_fixture();
    let rec = build_airr_record(&outcome, &cfg, "clean");
    let issues = validate_airr_record(&rec, &outcome, &cfg);
    assert!(
        !issues.iter().any(|i| matches!(
            i,
            RecordValidationIssue::NFwr1MutationsMismatch { .. }
                | RecordValidationIssue::NCdr1MutationsMismatch { .. }
                | RecordValidationIssue::NFwr2MutationsMismatch { .. }
                | RecordValidationIssue::NCdr2MutationsMismatch { .. }
                | RecordValidationIssue::NFwr3MutationsMismatch { .. }
                | RecordValidationIssue::NVUnannotatedMutationsMismatch { .. }
                | RecordValidationIssue::VSubregionMutationCountSumMismatch { .. }
        )),
        "engine-projected V-subregion partition must validate clean; got {issues:?}",
    );
}
