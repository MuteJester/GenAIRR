use std::sync::Arc;

use super::{NucHandle, Nucleotide, NucleotidePool, Region, Segment, Sequence};

/// Default initial capacity for `Simulation::new()` pools.
const DEFAULT_POOL_CAPACITY: usize = 500;

/// The root of one simulation. Owns the nucleotide pool, assembled
/// sequence, allele assignments, and three derived-state sidecars
/// (live calls, dirty-window log, mutation counter).
///
/// The sidecars used to live in a single `Option<Arc<LiveCallState>>`
/// grab bag whose six+ defensive `.as_ref().map(...).unwrap_or_default()`
/// readers were a textbook anti-pattern. They're now three focused
/// fields — see `segment_calls`, `dirty_log`, `mutation_count` for
/// the writer-reader contract of each.
#[derive(Clone, Debug, Default)]
pub struct Simulation {
    /// Arena of all nucleotides.
    pub pool: NucleotidePool,

    /// The assembled sequence (initially empty).
    pub sequence: Arc<Sequence>,

    /// V/D/J/C allele instances sampled for this simulation.
    pub assignments: crate::assignment::AlleleAssignments,

    /// V/D/J live-call evidence and its monotonic version chain.
    /// Writers: `with_assembled_segment_live_call`,
    /// `SimulationBuilder::seal_with_committed_live_calls`. Readers:
    /// AIRR projection, the live-call walker fast-path, downstream
    /// metadata-truth tests.
    pub segment_calls: Arc<crate::live_call::SegmentCalls>,

    /// Per-pass dirty-window message log. Writer:
    /// `DirtySignalObserver` (drained at seal time). Reader:
    /// `LiveCallRefreshHook` (then cleared via [`Self::with_dirty_log`]).
    pub dirty_log: Arc<crate::live_call::DirtyLog>,

    /// Per-record running count of base mutations applied by S5F /
    /// uniform mutation passes. Read by AIRR projection's
    /// `n_mutations` field.
    pub mutation_count: u32,
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
            segment_calls: Arc::new(crate::live_call::SegmentCalls::empty()),
            dirty_log: Arc::new(crate::live_call::DirtyLog::empty()),
            mutation_count: 0,
        }
    }

    /// Return a new simulation with the given V/D/J live-call state.
    pub fn with_segment_calls(&self, calls: crate::live_call::SegmentCalls) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments,
            segment_calls: Arc::new(calls),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
        }
    }

    /// Return a new simulation with the given dirty-window log.
    pub fn with_dirty_log(&self, log: crate::live_call::DirtyLog) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: Arc::new(log),
            mutation_count: self.mutation_count,
        }
    }

    /// Return a new simulation with the given mutation count.
    pub fn with_mutation_count(&self, count: u32) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: count,
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
                segment_calls: self.segment_calls.clone(),
                dirty_log: self.dirty_log.clone(),
                mutation_count: self.mutation_count,
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
                segment_calls: self.segment_calls.clone(),
                dirty_log: self.dirty_log.clone(),
                mutation_count: self.mutation_count,
            },
            range,
        )
    }

    /// Return a new simulation with the nucleotide at `handle` replaced.
    pub fn with_nucleotide_changed(&self, handle: NucHandle, new_n: Nucleotide) -> Self {
        Self {
            pool: self.pool.with_nucleotide_changed(handle, new_n),
            sequence: self.sequence.clone(),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
        }
    }

    /// Return a new simulation with only the base at `handle` changed.
    pub fn with_base_changed(&self, handle: NucHandle, new_base: u8) -> Self {
        Self {
            pool: self.pool.with_base_changed(handle, new_base),
            sequence: self.sequence.clone(),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
        }
    }

    /// Return a new simulation with `region` added to the sequence.
    pub fn with_region_added(&self, region: Region) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: Arc::new(self.sequence.with_region_added(region)),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
        }
    }

    /// Replace the region whose `(segment, start)` matches
    /// `(replacement.segment, replacement.start)` with `replacement`.
    ///
    /// If no matching region is found the returned simulation is
    /// unchanged.
    pub fn with_region_replaced_for_segment(&self, replacement: Region) -> Self {
        let idx = match self
            .sequence
            .regions
            .iter()
            .position(|r| r.segment == replacement.segment && r.start == replacement.start)
        {
            Some(idx) => idx,
            None => return self.clone(),
        };
        Self {
            pool: self.pool.clone(),
            sequence: Arc::new(self.sequence.with_region_replaced(idx, replacement)),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
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
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
        }
    }

    /// Return a new simulation with the trim at `segment` updated.
    pub fn with_trim(&self, segment: Segment, end: crate::assignment::TrimEnd, value: u16) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments.with_trim(segment, end, value),
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
        }
    }

    /// Return a new simulation with `n` inserted at pool position
    /// `at`. The persistent layer adjusts every region's `(start,
    /// end)` for the +1 shift; codon-rail metadata is no longer
    /// stored on `Region`, so callers that need it compute on
    /// demand via [`crate::ir::compute_codon_rail`].
    pub fn with_indel_inserted(&self, at: u32, n: Nucleotide) -> Self {
        let new_pool = self.pool.with_inserted(at, n);
        let adjusted = self.sequence.with_indel_adjusted(at, 1);
        Self {
            pool: new_pool,
            sequence: Arc::new(adjusted),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
        }
    }

    /// Return a new simulation with the nucleotide at pool position
    /// `at` removed. Same range-shift semantics as
    /// [`Self::with_indel_inserted`].
    pub fn with_indel_deleted(&self, at: u32) -> Self {
        let new_pool = self.pool.with_deleted(at);
        let adjusted = self.sequence.with_indel_adjusted(at, -1);
        Self {
            pool: new_pool,
            sequence: Arc::new(adjusted),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
        }
    }
}
