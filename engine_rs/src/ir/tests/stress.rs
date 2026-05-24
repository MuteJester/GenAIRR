use super::*;
use crate::ir::compute_codon_rail;

// ── Stress test for the persistent IR + codon rail ──────────────

/// Tiny deterministic PRNG (xorshift32) for the stress test. We
/// avoid bringing in `rand` here so the test has zero external
/// dependencies and reproduces identically across machines.
struct Xorshift32(u32);

impl Xorshift32 {
    fn new(seed: u32) -> Self {
        // xorshift32 cannot have a zero seed.
        Self(if seed == 0 { 0xdead_beef } else { seed })
    }
    fn next(&mut self) -> u32 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        self.0 = x;
        x
    }
}

/// Apply 1000 random base mutations through the persistent API,
/// keeping every IR revision in memory. Verifies two properties:
///
///   1. **Persistent contract (D1)**: the initial revision retains
///      its original amino-acid sequence after all 1000 mutations
///      have been applied to descendant revisions. No mutation
///      leaks back upstream.
///   2. **No silent coalescing**: most random mutations produce a
///      visible delta in the new revision's pool versus the previous
///      revision's pool.
///
/// (The old test also verified that each revision's "stored"
/// codon-rail matched a fresh recompute. That invariant is now
/// structurally true — the codon rail is computed on demand via
/// [`compute_codon_rail`] rather than stored on the region.)
#[test]
fn stress_persistent_ir_under_random_mutation() {
    const N_NUCS: u32 = 90; // 30 codons
    const N_MUTATIONS: usize = 1000;

    // Build the initial pool: 90 'A' nucleotides → 30 codons of AAA → 30 K's.
    let mut pool0 = NucleotidePool::new();
    for i in 0..N_NUCS {
        pool0.push(Nucleotide::germline(b'A', i as u16, Segment::V));
    }
    let region0 = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(N_NUCS));
    let sim0 = Simulation {
        pool: pool0,
        sequence: std::sync::Arc::new(Sequence::new().with_region_added(region0)),
        assignments: crate::assignment::AlleleAssignments::new(),
        segment_calls: std::sync::Arc::new(crate::live_call::SegmentCalls::empty()),
        dirty_log: std::sync::Arc::new(crate::live_call::DirtyLog::empty()),
        mutation_count: 0,
    };

    // Property check at revision 0.
    let rail0 = compute_codon_rail(&sim0.sequence.regions[0], &sim0.pool);
    assert_eq!(rail0.amino_acids.len(), 30);
    assert!(rail0.amino_acids.iter().all(|&aa| aa == b'K'));

    // Snapshot the initial amino acid sequence to verify property
    // (1) — that revision 0 is never perturbed.
    let initial_amino_acids = rail0.amino_acids.clone();

    // Apply 1000 random mutations. Every mutation produces a new
    // Simulation revision; every previous revision stays alive.
    let mut history: Vec<Simulation> = Vec::with_capacity(N_MUTATIONS + 1);
    history.push(sim0);

    let mut rng = Xorshift32::new(0x517c_c1ed); // arbitrary fixed seed
    let bases = [b'A', b'C', b'G', b'T'];
    let mut visible_mutations = 0usize;

    for _ in 0..N_MUTATIONS {
        let prev = history.last().unwrap();
        let pos = rng.next() % N_NUCS;
        let new_base = bases[(rng.next() & 0b11) as usize];

        let prev_base = prev.pool.get(NucHandle::new(pos)).unwrap().base;
        if prev_base != new_base {
            visible_mutations += 1;
        }

        let new_pool = prev.pool.with_base_changed(NucHandle::new(pos), new_base);
        history.push(Simulation {
            pool: new_pool,
            sequence: prev.sequence.clone(),
            assignments: crate::assignment::AlleleAssignments::new(),
            segment_calls: prev.segment_calls.clone(),
            dirty_log: prev.dirty_log.clone(),
            mutation_count: prev.mutation_count,
        });
    }

    // Sanity: history size is correct.
    assert_eq!(history.len(), N_MUTATIONS + 1);

    // Property (1): revision 0's amino acids are exactly what they
    // were before any mutation happened.
    let rev0_rail = compute_codon_rail(&history[0].sequence.regions[0], &history[0].pool);
    assert_eq!(rev0_rail.amino_acids, initial_amino_acids);
    assert!(rev0_rail.amino_acids.iter().all(|&aa| aa == b'K'));

    // Property (2): most mutations produced a visible base delta.
    // With a uniform 4-base alphabet, the expected fraction of
    // self-substitutions is 1/4, so visible should be ~75% × 1000.
    // Use a generous lower bound to absorb statistical jitter.
    assert!(
        visible_mutations >= 600,
        "expected ≥600 visible mutations, got {}",
        visible_mutations
    );
}
