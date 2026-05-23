use std::sync::Arc;

use super::{NucHandle, Nucleotide, NucleotidePool, Region, Segment, Sequence};

/// Default initial capacity for `Simulation::new()` pools.
const DEFAULT_POOL_CAPACITY: usize = 500;

/// The root of one simulation. Owns the nucleotide pool, assembled
/// sequence, allele assignments, and optional live-call evidence.
#[derive(Clone, Debug, Default)]
pub struct Simulation {
    /// Arena of all nucleotides.
    pub pool: NucleotidePool,

    /// The assembled sequence (initially empty).
    pub sequence: Arc<Sequence>,

    /// V/D/J/C allele instances sampled for this simulation.
    pub assignments: crate::assignment::AlleleAssignments,

    /// Dynamic V/D/J call evidence over the current observed sequence.
    pub live_calls: Option<Arc<crate::live_call::LiveCallState>>,
}

impl Simulation {
    /// Empty simulation.
    pub fn new() -> Self {
        Self::with_capacity(DEFAULT_POOL_CAPACITY)
    }

    /// Empty simulation with `nuc_capacity` reserved up front.
    pub fn with_capacity(nuc_capacity: usize) -> Self {
        Self {
            pool: NucleotidePool::with_capacity(nuc_capacity),
            sequence: Arc::new(Sequence::new()),
            assignments: crate::assignment::AlleleAssignments::new(),
            live_calls: None,
        }
    }

    /// Return a new simulation with dynamic live-call evidence set.
    pub fn with_live_calls(&self, live_calls: crate::live_call::LiveCallState) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments,
            live_calls: Some(Arc::new(live_calls)),
        }
    }

    #[must_use = "with_nucleotide_pushed returns (Simulation, NucHandle); \
                  both must be used or destructured. Drop the handle \
                  explicitly with `let (next, _) = ...` if it isn't needed."]
    pub fn with_nucleotide_pushed(&self, n: Nucleotide) -> (Self, NucHandle) {
        let (pool, h) = self.pool.with_pushed(n);
        (
            Self {
                pool,
                sequence: self.sequence.clone(),
                assignments: self.assignments,
                live_calls: self.live_calls.clone(),
            },
            h,
        )
    }

    #[must_use]
    pub fn with_nucleotides_extended<I: IntoIterator<Item = Nucleotide>>(
        &self,
        iter: I,
    ) -> (Self, std::ops::Range<u32>) {
        let (pool, range) = self.pool.with_pushed_many(iter);
        (
            Self {
                pool,
                sequence: self.sequence.clone(),
                assignments: self.assignments,
                live_calls: self.live_calls.clone(),
            },
            range,
        )
    }

    /// Return a new simulation with the nucleotide at `handle` replaced.
    pub fn with_nucleotide_changed(&self, handle: NucHandle, new_n: Nucleotide) -> Self {
        let new_pool = self.pool.with_nucleotide_changed(handle, new_n);
        let new_sequence = refresh_regions_covering(&self.sequence, &new_pool, handle);
        Self {
            pool: new_pool,
            sequence: Arc::new(new_sequence),
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }

    /// Return a new simulation with only the base at `handle` changed.
    pub fn with_base_changed(&self, handle: NucHandle, new_base: u8) -> Self {
        let new_pool = self.pool.with_base_changed(handle, new_base);
        let new_sequence = refresh_regions_covering(&self.sequence, &new_pool, handle);
        Self {
            pool: new_pool,
            sequence: Arc::new(new_sequence),
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }

    /// Return a new simulation with `region` added to the sequence.
    pub fn with_region_added(&self, region: Region) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: Arc::new(self.sequence.with_region_added(region)),
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }

    /// Phase 7: replace the region whose `(segment, start)` matches
    /// `(replacement.segment, replacement.start)` with `replacement`.
    ///
    /// Used by post-assembly mutation passes (S5F, …) to write
    /// observer-produced codon-rail data back into the existing
    /// region without going through `with_region_added` (which
    /// would create a duplicate). If no matching region is found
    /// the returned simulation is unchanged — same behaviour as
    /// the pre-Phase-7 path where the rail update was deferred to
    /// the post-pass refresh.
    pub fn with_region_replaced_for_segment(&self, replacement: Region) -> Self {
        let idx = match self.sequence.regions.iter().position(|r| {
            r.segment == replacement.segment && r.start == replacement.start
        }) {
            Some(idx) => idx,
            None => return self.clone(),
        };
        Self {
            pool: self.pool.clone(),
            sequence: Arc::new(self.sequence.with_region_replaced(idx, replacement)),
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }

    /// Return a new simulation with `instance` assigned to `segment`.
    pub fn with_allele_assigned(
        &self,
        segment: Segment,
        instance: crate::assignment::AlleleInstance,
    ) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments.with_assigned(segment, instance),
            live_calls: self.live_calls.clone(),
        }
    }

    /// Return a new simulation with the trim at `segment` updated.
    pub fn with_trim(&self, segment: Segment, end: crate::assignment::TrimEnd, value: u16) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments.with_trim(segment, end, value),
            live_calls: self.live_calls.clone(),
        }
    }

    /// Return a new simulation with `n` inserted at pool position `at`.
    pub fn with_indel_inserted(&self, at: u32, n: Nucleotide) -> Self {
        let new_pool = self.pool.with_inserted(at, n);
        let adjusted = self.sequence.with_indel_adjusted(at, 1);
        let new_regions: Vec<Region> = adjusted
            .regions
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let original = &self.sequence.regions[i];
                if r.start != original.start
                    || r.end != original.end
                    || r.frame_phase != original.frame_phase
                {
                    r.with_codon_rail_recomputed(&new_pool)
                } else {
                    r.clone()
                }
            })
            .collect();
        Self {
            pool: new_pool,
            sequence: Arc::new(Sequence {
                regions: new_regions,
            }),
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }

    /// Return a new simulation with the nucleotide at pool position `at` removed.
    pub fn with_indel_deleted(&self, at: u32) -> Self {
        let new_pool = self.pool.with_deleted(at);
        let adjusted = self.sequence.with_indel_adjusted(at, -1);
        let new_regions: Vec<Region> = adjusted
            .regions
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let original = &self.sequence.regions[i];
                if r.start != original.start
                    || r.end != original.end
                    || r.frame_phase != original.frame_phase
                {
                    r.with_codon_rail_recomputed(&new_pool)
                } else {
                    r.clone()
                }
            })
            .collect();
        Self {
            pool: new_pool,
            sequence: Arc::new(Sequence {
                regions: new_regions,
            }),
            assignments: self.assignments,
            live_calls: self.live_calls.clone(),
        }
    }
}

fn refresh_regions_covering(seq: &Sequence, pool: &NucleotidePool, handle: NucHandle) -> Sequence {
    let h = handle.index();
    let new_regions: Vec<Region> = seq
        .regions
        .iter()
        .map(|r| {
            if h >= r.start.index() && h < r.end.index() {
                r.with_codon_rail_recomputed(pool)
            } else {
                r.clone()
            }
        })
        .collect();
    Sequence {
        regions: new_regions,
    }
}
