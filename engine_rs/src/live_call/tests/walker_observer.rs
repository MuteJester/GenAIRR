//! Streamed-vs-from-scratch equivalence tests for the live-call walker.
//!
//! These tests anchor the walker refactor invariant: the
//! `WalkerObserverState` accumulating per-allele scores inline with
//! `SimulationBuilder::push_nucleotide` must produce the **exact same**
//! `SegmentLiveCall` (byte-for-byte equality on every public field)
//! that the from-scratch `call_from_region` walker produces on the
//! fully-assembled region. The from-scratch walker is the property
//! oracle and lives unchanged at [`super::super::walker`].
//!
//! Each test:
//!
//! 1. Builds a hand-crafted refdata config + assembled region (the same
//!    fixtures the existing `caller.rs` tests use).
//! 2. Drives the observer through the same byte sequence via a
//!    `SimulationBuilder` and seals it (running the post-seal extensions).
//! 3. Calls `assembled_segment_live_call` for the same fixture and
//!    compares the two `SegmentLiveCall`s on every field.
//!
//! The equality check covers `seq_start`, `seq_end`, `ref_start`,
//! `ref_end`, `allele_call`, `confidence`, `boundary_summary`,
//! `hypotheses`, and `evidence_version` — exactly the surface the
//! user flagged as the load-bearing ambiguity-tracking contract.

use super::super::walker_observer::{SealedWalkerState, WalkerObserverState};
use super::super::{
    assembled_segment_live_call, ReferenceMatchIndex, SegmentLiveCall, SegmentRefIndex,
};
use super::{allele, simulation_with_region};
use crate::ir::{flag, NucFlags, NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::refdata::{ChainType, RefDataConfig};

/// Drive a `WalkerObserverState` through every base of `sim`'s latest
/// region for `segment`, then seal + run extensions, returning the
/// final `SegmentLiveCall`.
///
/// Mirrors the assembly-pass code path that constructs the observer
/// call. Importantly, it reads from a *fully built* simulation (the
/// fixture builders in `tests/mod.rs` use the persistent
/// `with_nucleotide_pushed` API), which means the extension walks have
/// access to neighbour NP regions just like they would at assembly
/// time when run from `execute_streaming`.
fn observer_call_for(
    sim: &Simulation,
    segment_index: &SegmentRefIndex,
    segment: Segment,
    evidence_version: u64,
) -> SegmentLiveCall {
    let region = sim
        .sequence
        .regions
        .iter()
        .rev()
        .find(|r| r.segment == segment)
        .expect("test fixture must have a region for the requested segment");
    let seq_start = region.start.index();
    let seq_end = region.end.index();

    let mut observer = WalkerObserverState::new(segment_index, seq_start);
    for seq_pos in seq_start..seq_end {
        let nuc = *sim
            .pool
            .get(NucHandle::new(seq_pos))
            .expect("region range must be in pool");
        observer.on_base_pushed(&nuc);
    }
    let sealed = observer.seal(seq_end);
    sealed.finalize_with_extensions(sim, segment_index, evidence_version, seq_start, seq_end)
}

fn assert_calls_equal(observer: &SegmentLiveCall, oracle: &SegmentLiveCall) {
    assert_eq!(observer.segment, oracle.segment, "segment mismatch");
    assert_eq!(
        observer.evidence_version, oracle.evidence_version,
        "evidence_version mismatch"
    );
    assert_eq!(
        observer.confidence, oracle.confidence,
        "confidence mismatch"
    );
    assert_eq!(
        observer.allele_call, oracle.allele_call,
        "allele_call (tie-set bitset) mismatch"
    );
    assert_eq!(
        observer.boundary_summary, oracle.boundary_summary,
        "boundary_summary mismatch"
    );
    assert_eq!(
        observer.hypotheses, oracle.hypotheses,
        "hypotheses vector mismatch (includes seq_start/seq_end/ref_start/ref_end/score/flags)"
    );
}

#[test]
fn observer_matches_oracle_trim_induced_ambiguity() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(allele(Segment::V, "v1*01", b"GGAAACCC"));
    let _ = cfg.v_pool.push(allele(Segment::V, "v2*01", b"TTAAACCC"));
    let _ = cfg.v_pool.push(allele(Segment::V, "v3*01", b"GGTATCCC"));
    let index = ReferenceMatchIndex::build(&cfg);
    let sim = simulation_with_region(Segment::V, b"AAACCC", 2);

    let observer = observer_call_for(&sim, &index.v, Segment::V, 1);
    let oracle = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("oracle must produce a call");
    assert_calls_equal(&observer, &oracle);
}

#[test]
fn observer_matches_oracle_mutated_base_disambiguates() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(allele(Segment::V, "v1*01", b"AACT"));
    let _ = cfg.v_pool.push(allele(Segment::V, "v1*02", b"AAGT"));
    let index = ReferenceMatchIndex::build(&cfg);
    let sim = simulation_with_region(Segment::V, b"AACT", 0);
    let sim = sim.with_base_changed(NucHandle::new(2), b'G');

    let observer = observer_call_for(&sim, &index.v, Segment::V, 1);
    let oracle = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("oracle must produce a call");
    assert_calls_equal(&observer, &oracle);
}

#[test]
fn observer_matches_oracle_truth_allele_retained_under_single_mutation() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(allele(Segment::V, "v*01", b"ACGTAC"));
    let _ = cfg.v_pool.push(allele(Segment::V, "v*02", b"ACGTAG"));
    let _ = cfg.v_pool.push(allele(Segment::V, "v*03", b"AGGTAC"));
    let index = ReferenceMatchIndex::build(&cfg);

    let sim = simulation_with_region(Segment::V, b"ACGTAC", 0);
    let sim = sim.with_base_changed(NucHandle::new(2), b'T');

    let observer = observer_call_for(&sim, &index.v, Segment::V, 1);
    let oracle = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("oracle must produce a call");
    assert_calls_equal(&observer, &oracle);
}

#[test]
fn observer_matches_oracle_flip_toward_other_allele() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(allele(Segment::V, "v*01", b"ACGTAC"));
    let _ = cfg.v_pool.push(allele(Segment::V, "v*02", b"AGGTGA"));
    let index = ReferenceMatchIndex::build(&cfg);

    let sim = simulation_with_region(Segment::V, b"ACGTAC", 0);
    let sim = sim.with_base_changed(NucHandle::new(1), b'G');
    let sim = sim.with_base_changed(NucHandle::new(4), b'G');
    let sim = sim.with_base_changed(NucHandle::new(5), b'A');

    let observer = observer_call_for(&sim, &index.v, Segment::V, 1);
    let oracle = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("oracle must produce a call");
    assert_calls_equal(&observer, &oracle);
}

#[test]
fn observer_matches_oracle_full_pool_when_zero_score() {
    // Every base flipped to something inconsistent with every allele.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(allele(Segment::V, "v*01", b"AACT"));
    let _ = cfg.v_pool.push(allele(Segment::V, "v*02", b"GGCA"));
    let index = ReferenceMatchIndex::build(&cfg);

    let sim = simulation_with_region(Segment::V, b"AACT", 0);
    let sim = sim.with_base_changed(NucHandle::new(0), b'T');
    let sim = sim.with_base_changed(NucHandle::new(1), b'T');
    let sim = sim.with_base_changed(NucHandle::new(2), b'A');
    let sim = sim.with_base_changed(NucHandle::new(3), b'G');

    let observer = observer_call_for(&sim, &index.v, Segment::V, 1);
    let oracle = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("oracle must produce a call");
    assert_calls_equal(&observer, &oracle);
}

#[test]
fn observer_matches_oracle_n_wildcard_path() {
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(allele(Segment::V, "v1*01", b"AACT"));
    let _ = cfg.v_pool.push(allele(Segment::V, "v1*02", b"AAGT"));
    let index = ReferenceMatchIndex::build(&cfg);
    let sim = simulation_with_region(Segment::V, b"AACT", 0);
    let sim = sim.with_base_changed(NucHandle::new(2), b'N');

    let observer = observer_call_for(&sim, &index.v, Segment::V, 1);
    let oracle = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("oracle must produce a call");
    assert_calls_equal(&observer, &oracle);
}

#[test]
fn observer_matches_oracle_trim_ambiguity_narrows_via_np_extension() {
    // V trimmed to a prefix that's ambiguous over three alleles; an NP
    // region whose first bases extend one allele's suffix breaks the tie.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(allele(Segment::V, "v*01", b"AACCAA"));
    let _ = cfg.v_pool.push(allele(Segment::V, "v*02", b"AACCTG"));
    let _ = cfg.v_pool.push(allele(Segment::V, "v*03", b"AACCGG"));
    let index = ReferenceMatchIndex::build(&cfg);

    let mut sim = Simulation::new();
    for (i, &b) in b"AACC".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(4),
    ));
    for &b in &[b'T', b'G'] {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::synthetic(b, Segment::Np1, NucFlags::empty()));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::Np1,
        NucHandle::new(4),
        NucHandle::new(6),
    ));

    let observer = observer_call_for(&sim, &index.v, Segment::V, 1);
    let oracle = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("oracle must produce a call");
    assert_calls_equal(&observer, &oracle);
}

#[test]
fn observer_matches_oracle_cross_segment_nucleotide_produces_unsupported() {
    // Manually build a region where one of the nucleotides claims the
    // wrong segment — the walker treats this as malformed IR.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(allele(Segment::V, "v*01", b"AACT"));
    let index = ReferenceMatchIndex::build(&cfg);

    let mut sim = Simulation::new();
    let bases = [
        Nucleotide::germline(b'A', 0, Segment::V),
        Nucleotide::germline(b'A', 1, Segment::V),
        // Cross-segment nucleotide inside the V region -> malformed IR.
        Nucleotide::synthetic(b'C', Segment::Np1, NucFlags::empty()),
        Nucleotide::germline(b'T', 3, Segment::V),
    ];
    for n in bases {
        let (next, _) = sim.with_nucleotide_pushed(n);
        sim = next;
    }
    let sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(4),
    ));

    let observer = observer_call_for(&sim, &index.v, Segment::V, 1);
    let oracle = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("oracle must produce a call");
    assert_calls_equal(&observer, &oracle);
    assert_eq!(observer.hypotheses, Vec::new(), "unsupported call has no hypotheses");
}

#[test]
fn observer_matches_oracle_indel_insert_skipped_silently() {
    // An indel-inserted nucleotide (GermlinePos::NONE) inside the V
    // region must be skipped by the scoring loop but still
    // contribute to the seq range. The walker tolerates this; the
    // observer must too.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(allele(Segment::V, "v*01", b"AACT"));
    let index = ReferenceMatchIndex::build(&cfg);

    let mut sim = Simulation::new();
    let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b'A', 0, Segment::V));
    sim = next;
    let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b'A', 1, Segment::V));
    sim = next;
    // Indel-inserted nucleotide carrying segment::V but no germline_pos.
    let (next, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
        b'G',
        Segment::V,
        flag::INDEL_INSERTED,
    ));
    sim = next;
    let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b'C', 2, Segment::V));
    sim = next;
    let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b'T', 3, Segment::V));
    sim = next;
    let sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(5),
    ));

    let observer = observer_call_for(&sim, &index.v, Segment::V, 1);
    let oracle = assembled_segment_live_call(&sim, &index, Segment::V, 1)
        .expect("oracle must produce a call");
    assert_calls_equal(&observer, &oracle);
}

#[test]
fn observer_matches_oracle_d_segment_left_extension_into_np1() {
    // Build V + NP1 + D and verify the observer's D call (with the
    // left extension into NP1) matches the from-scratch oracle.
    let mut cfg = RefDataConfig::empty(ChainType::Vdj);
    let _ = cfg.v_pool.push(allele(Segment::V, "v*01", b"AACT"));
    let _ = cfg.d_pool.push(allele(Segment::D, "d*01", b"GGCCAA"));
    let _ = cfg.d_pool.push(allele(Segment::D, "d*02", b"GGCCTT"));
    let index = ReferenceMatchIndex::build(&cfg);

    let mut sim = Simulation::new();
    for (i, &b) in b"AACT".iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(4),
    ));
    for &b in &[b'G', b'G'] {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::synthetic(b, Segment::Np1, NucFlags::empty()));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::Np1,
        NucHandle::new(4),
        NucHandle::new(6),
    ));
    // D allele "d*01" trimmed_5=2 → bases CCAA at ref pos 2,3,4,5.
    for (i, &b) in b"CCAA".iter().enumerate() {
        let (next, _) =
            sim.with_nucleotide_pushed(Nucleotide::germline(b, (i + 2) as u16, Segment::D));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::D,
        NucHandle::new(6),
        NucHandle::new(10),
    ));

    let observer = observer_call_for(&sim, &index.d, Segment::D, 1);
    let oracle = assembled_segment_live_call(&sim, &index, Segment::D, 1)
        .expect("oracle must produce a call");
    assert_calls_equal(&observer, &oracle);
}
