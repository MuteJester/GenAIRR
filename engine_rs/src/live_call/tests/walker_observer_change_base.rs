//! property tests for `WalkerObserverState::on_base_changed`.
//!
//! Drives the streaming walker observer through a `push_nucleotide`
//! → `change_base` sequence on a `SimulationBuilder`, then verifies
//! the sealed `SegmentLiveCall` equals the from-scratch oracle
//! `call_from_region` running on the mutated simulation.

use super::allele;
use crate::ir::{Nucleotide, Region, Segment, Simulation, SimulationBuilder};
use crate::live_call::reference_index::{SegmentRefIndex, DEFAULT_REFERENCE_KMER_LEN};
use crate::live_call::walker_observer::SealedWalkerState;
use crate::refdata::AllelePool;

fn build_v_index() -> SegmentRefIndex {
    // Three V alleles, length 4, differing at one or two ref positions.
    let mut pool = AllelePool::new();
    let _ = pool.push(allele(Segment::V, "v1*01", b"AACT")); // id 0
    let _ = pool.push(allele(Segment::V, "v1*02", b"AACG")); // id 1, differs from v1*01 at pos 3
    let _ = pool.push(allele(Segment::V, "v2*01", b"AGCT")); // id 2, differs at pos 1
    SegmentRefIndex::build(Segment::V, &pool, DEFAULT_REFERENCE_KMER_LEN)
}

/// Drive a walker observer through `bases` then apply `(handle, new_base)`
/// edits, returning the sealed scores + match counts.
fn drive(
    index: &SegmentRefIndex,
    bases: &[u8],
    edits: &[(u32, u8)],
) -> (Vec<u32>, u32, u32) {
    let mut builder = SimulationBuilder::from_simulation(Simulation::new());
    builder.attach_walker_observer(index, 0);
    for (i, &b) in bases.iter().enumerate() {
        builder.push_nucleotide(Nucleotide::germline(b, i as u16, Segment::V));
    }
    for &(handle, new_base) in edits {
        builder.change_base(crate::ir::NucHandle::new(handle), new_base);
    }
    let sealed = builder.seal_walker_observer(bases.len() as u32);
    match sealed {
        SealedWalkerState::Resolved(r) => (
            r.scores,
            r.informative_matches,
            r.wildcard_matches,
        ),
        SealedWalkerState::Unsupported { .. } | SealedWalkerState::Unresolved { .. } => {
            panic!("test expected Resolved state")
        }
    }
}

/// Build a fresh observer driven over `bases` with no edits — the
/// score state any post-edit observer should converge to.
fn fresh_after(index: &SegmentRefIndex, bases: &[u8]) -> (Vec<u32>, u32, u32) {
    drive(index, bases, &[])
}

#[test]
fn change_base_no_op_when_old_equals_new_base() {
    let idx = build_v_index();
    // AACT — matches v1*01 fully (4 informative matches), v1*02 at 0..3 (3 matches), v2*01 at 0,2,3 (3 matches).
    let bases = b"AACT";
    let (scores_no_edit, inform_no, wild_no) = fresh_after(&idx, bases);
    let (scores_edit, inform_e, wild_e) = drive(&idx, bases, &[(2, b'C')]); // no-op: pos 2 already 'C'
    assert_eq!(scores_no_edit, scores_edit);
    assert_eq!(inform_no, inform_e);
    assert_eq!(wild_no, wild_e);
}

#[test]
fn change_base_at_distinguishing_position_matches_fresh_drive() {
    let idx = build_v_index();
    // Start with AACT; change pos 3 from T → G.
    // Equivalent to driving fresh over AACG.
    let start = b"AACT";
    let target = b"AACG";
    let (scores_target, inform_target, wild_target) = fresh_after(&idx, target);
    let (scores_edit, inform_edit, wild_edit) = drive(&idx, start, &[(3, b'G')]);
    assert_eq!(scores_target, scores_edit, "scores diverge after edit");
    assert_eq!(inform_target, inform_edit);
    assert_eq!(wild_target, wild_edit);
}

#[test]
fn multiple_edits_compose_to_fresh_drive() {
    let idx = build_v_index();
    // Start with AACT; change pos 1 (A→G) AND pos 3 (T→T no-op) AND pos 3 (T→G).
    // Final sequence should be AGCG.
    let start = b"AACT";
    let target = b"AGCG";
    let (scores_target, inform_target, wild_target) = fresh_after(&idx, target);
    let (scores_edit, inform_edit, wild_edit) = drive(
        &idx,
        start,
        &[(1, b'G'), (3, b'T'), (3, b'G')],
    );
    assert_eq!(scores_target, scores_edit);
    assert_eq!(inform_target, inform_edit);
    assert_eq!(wild_target, wild_edit);
}

#[test]
fn change_base_to_n_treats_as_wildcard_evidence() {
    let idx = build_v_index();
    // Start with AACT; change pos 1 (A → N). N matches all alleles via
    // the wildcard codepath (all four base bitsets unioned).
    let start = b"AACT";
    let target = b"ANCT";
    let (scores_target, inform_target, wild_target) = fresh_after(&idx, target);
    let (scores_edit, inform_edit, wild_edit) = drive(&idx, start, &[(1, b'N')]);
    assert_eq!(scores_target, scores_edit);
    assert_eq!(inform_target, inform_edit);
    assert_eq!(wild_target, wild_edit);
}

#[test]
fn change_base_against_callfromregion_oracle() {
    use crate::ir::NucHandle;
    use crate::live_call::with_assembled_segment_live_call;
    use crate::live_call::ReferenceMatchIndex;

    // Build a real ReferenceMatchIndex via the empty-config helper
    // + manual allele pushes.
    let mut refdata = crate::refdata::RefDataConfig::empty(crate::refdata::ChainType::Vdj);
    let _ = refdata.v_pool.push(allele(Segment::V, "v1*01", b"AACT"));
    let _ = refdata.v_pool.push(allele(Segment::V, "v1*02", b"AACG"));
    let _ = refdata.v_pool.push(allele(Segment::V, "v2*01", b"AGCT"));
    let ref_index = ReferenceMatchIndex::build(&refdata);

    // Build the same sim through the builder + observer path with
    // an edit, then compare its `SegmentLiveCall` against the
    // from-scratch oracle on the equivalent mutated Simulation.
    let mut builder = SimulationBuilder::from_simulation(Simulation::new());
    let seg_index = ref_index.get(Segment::V).expect("V index exists");
    builder.attach_walker_observer(seg_index, 0);
    for (i, &b) in b"AACT".iter().enumerate() {
        builder.push_nucleotide(Nucleotide::germline(b, i as u16, Segment::V));
    }
    // Mutate pos 3 from T to G.
    builder.change_base(NucHandle::new(3), b'G');
    let sealed = builder.seal_walker_observer(4);
    let observer_call = sealed.finalize_with_extensions(
        builder.peek(),
        seg_index,
        /* evidence_version */ 1,
        /* seq_start */ 0,
        /* seq_end */ 4,
    );

    // Oracle: build a fresh sim with the mutated bytes and call
    // `with_assembled_segment_live_call`, which under the
    // observer-aware fast path WOULD reuse a staged call — but
    // we're feeding a sim without a staged call so it routes
    // through the slow `call_from_region` oracle.
    let mut oracle_sim = Simulation::new();
    for (i, &b) in b"AACG".iter().enumerate() {
        let (next, _) = oracle_sim.with_nucleotide_pushed(Nucleotide::germline(
            b,
            i as u16,
            Segment::V,
        ));
        oracle_sim = next;
    }
    oracle_sim = oracle_sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(4),
    ));
    let oracle_sim = with_assembled_segment_live_call(&oracle_sim, &ref_index, Segment::V);
    let oracle_call = oracle_sim.live_calls.as_ref().unwrap().v.as_ref().unwrap();

    // The two paths must agree on the tie-set, boundaries, and
    // evidence-score shape.
    assert_eq!(
        observer_call.allele_call.to_ids(),
        oracle_call.allele_call.to_ids(),
        "observer-via-change_base tie-set diverges from oracle"
    );
    assert_eq!(
        observer_call.boundary_summary, oracle_call.boundary_summary,
        "observer-via-change_base boundary diverges from oracle"
    );
    assert_eq!(observer_call.confidence, oracle_call.confidence);
}
