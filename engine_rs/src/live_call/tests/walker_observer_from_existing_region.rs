//! property tests for
//! `WalkerObserverState::from_existing_region`.
//!
//! Rebuild an observer from an already-assembled `Region` + the
//! current pool, then apply further edits via `on_base_changed`,
//! and verify the result matches a fresh observer driven over the
//! equivalent byte stream + edits — and, when finalized through
//! `with_assembled_segment_live_call`, matches the from-scratch
//! `call_from_region` oracle.

use super::allele;
use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation, SimulationBuilder};
use crate::live_call::reference_index::{SegmentRefIndex, DEFAULT_REFERENCE_KMER_LEN};
use crate::live_call::walker_observer::{SealedWalkerState, WalkerObserverState};
use crate::refdata::AllelePool;

fn build_v_index() -> SegmentRefIndex {
    let mut pool = AllelePool::new();
    let _ = pool.push(allele(Segment::V, "v1*01", b"AACT"));
    let _ = pool.push(allele(Segment::V, "v1*02", b"AACG"));
    let _ = pool.push(allele(Segment::V, "v2*01", b"AGCT"));
    SegmentRefIndex::build(Segment::V, &pool, DEFAULT_REFERENCE_KMER_LEN)
}

/// Build a small assembled sim by pushing bases through a fresh
/// observer (via builder.push_nucleotide), seal it, and return both
/// the sealed observer state and the sim itself.
fn assemble_via_observer(
    idx: &SegmentRefIndex,
    bases: &[u8],
) -> (Simulation, Region, SealedWalkerState) {
    let mut builder = SimulationBuilder::from_simulation(Simulation::new());
    builder.attach_walker_observer(idx, 0);
    for (i, &b) in bases.iter().enumerate() {
        builder.push_nucleotide(Nucleotide::germline(b, i as u16, Segment::V));
    }
    let sealed = builder.seal_walker_observer(bases.len() as u32);
    let sim = builder.seal();
    let region = Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(bases.len() as u32),
    );
    let sim_with_region = sim.with_region_added(region.clone());
    (sim_with_region, region, sealed)
}

fn drain_resolved_state(sealed: SealedWalkerState) -> (Vec<u32>, u32, u32) {
    match sealed {
        SealedWalkerState::Resolved(r) => {
            (r.scores, r.informative_matches, r.wildcard_matches)
        }
        _ => panic!("expected Resolved"),
    }
}

#[test]
fn rebuild_matches_fresh_observer_at_seal() {
    let idx = build_v_index();
    let bases = b"AACT";
    let (sim, region, fresh_sealed) = assemble_via_observer(&idx, bases);

    // Rebuild a new observer from the existing region and seal it
    // immediately. Should produce identical state to the fresh one.
    let rebuilt = WalkerObserverState::from_existing_region(&idx, &sim, &region);
    let rebuilt_sealed = rebuilt.seal(bases.len() as u32);

    let (s_fresh, i_fresh, w_fresh) = drain_resolved_state(fresh_sealed);
    let (s_rebuilt, i_rebuilt, w_rebuilt) = drain_resolved_state(rebuilt_sealed);
    assert_eq!(s_fresh, s_rebuilt, "scores diverge");
    assert_eq!(i_fresh, i_rebuilt, "informative_matches diverge");
    assert_eq!(w_fresh, w_rebuilt, "wildcard_matches diverge");
}

#[test]
fn rebuild_supports_subsequent_on_base_changed() {
    let idx = build_v_index();
    let bases = b"AACT";
    let (sim, region, _) = assemble_via_observer(&idx, bases);

    // Rebuild and apply an edit: pos 3 T → G.
    let mut rebuilt = WalkerObserverState::from_existing_region(&idx, &sim, &region);
    use crate::ir::builder::IrEventObserver;
    let old = *sim.pool.get(NucHandle::new(3)).unwrap();
    rebuilt.on_base_changed(NucHandle::new(3), &old, b'G');
    let rebuilt_sealed = rebuilt.seal(bases.len() as u32);

    // Compare to a fresh observer driven over the mutated byte stream.
    let (_, _, fresh_sealed) = assemble_via_observer(&idx, b"AACG");
    let (s_fresh, i_fresh, w_fresh) = drain_resolved_state(fresh_sealed);
    let (s_rebuilt, i_rebuilt, w_rebuilt) = drain_resolved_state(rebuilt_sealed);
    assert_eq!(s_fresh, s_rebuilt);
    assert_eq!(i_fresh, i_rebuilt);
    assert_eq!(w_fresh, w_rebuilt);
}

#[test]
fn rebuild_plus_extensions_matches_call_from_region_oracle() {
    use crate::live_call::with_assembled_segment_live_call;
    use crate::live_call::ReferenceMatchIndex;
    use crate::refdata::{ChainType, RefDataConfig};

    let mut refdata = RefDataConfig::empty(ChainType::Vdj);
    let _ = refdata.v_pool.push(allele(Segment::V, "v1*01", b"AACT"));
    let _ = refdata.v_pool.push(allele(Segment::V, "v1*02", b"AACG"));
    let _ = refdata.v_pool.push(allele(Segment::V, "v2*01", b"AGCT"));
    let ref_index = ReferenceMatchIndex::build(&refdata);
    let seg_index = ref_index.get(Segment::V).unwrap();

    // Assemble a V region via the persistent IR, then rebuild a
    // walker observer from it and apply an edit. Compare the
    // observer-produced live call to what `with_assembled_segment_live_call`
    // would compute on the mutated sim from scratch.
    let bases = b"AACT";
    let mut sim = Simulation::new();
    for (i, &b) in bases.iter().enumerate() {
        let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
            b,
            i as u16,
            Segment::V,
        ));
        sim = next;
    }
    sim = sim.with_region_added(Region::new(
        Segment::V,
        NucHandle::new(0),
        NucHandle::new(bases.len() as u32),
    ));

    // Rebuild observer from the assembled sim. Edit pos 3 T → G.
    let region_v = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V)
        .unwrap()
        .clone();
    let mut rebuilt = WalkerObserverState::from_existing_region(seg_index, &sim, &region_v);
    use crate::ir::builder::IrEventObserver;
    let old = *sim.pool.get(NucHandle::new(3)).unwrap();
    rebuilt.on_base_changed(NucHandle::new(3), &old, b'G');
    let rebuilt_sealed = rebuilt.seal(bases.len() as u32);
    // Apply the same edit to the sim so `finalize_with_extensions`
    // sees the mutated pool.
    let mutated_sim = sim.with_base_changed(NucHandle::new(3), b'G');
    let observer_call = rebuilt_sealed.finalize_with_extensions(
        &mutated_sim,
        seg_index,
        /* evidence_version */ 1,
        0,
        bases.len() as u32,
    );

    // Oracle: with_assembled_segment_live_call on the mutated sim.
    let oracle_sim =
        with_assembled_segment_live_call(&mutated_sim, &ref_index, Segment::V);
    let oracle_call = oracle_sim.segment_calls.get(Segment::V).unwrap();

    assert_eq!(
        observer_call.allele_call.to_ids(),
        oracle_call.allele_call.to_ids(),
        "rebuild-via-from_existing_region tie-set diverges from oracle"
    );
    assert_eq!(observer_call.boundary_summary, oracle_call.boundary_summary);
    assert_eq!(observer_call.confidence, oracle_call.confidence);
}

#[test]
fn rebuild_empty_region_yields_unresolved() {
    let idx = build_v_index();
    let sim = Simulation::new();
    let empty_region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(0));
    let rebuilt = WalkerObserverState::from_existing_region(&idx, &sim, &empty_region);
    match rebuilt.seal(0) {
        SealedWalkerState::Unresolved { .. } => {}
        other => panic!("expected Unresolved for empty region, got {:?}", std::mem::discriminant(&other)),
    }
}
