use super::super::build_airr_record;
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
