use std::sync::Arc;

use super::{NucHandle, Nucleotide, NucleotidePool, PoolRange, Region, Segment, Sequence};

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

    /// Donor genotype's immutable union of carried alleles across both
    /// chromosomes. None leaves inference reference-wide.
    pub live_call_support: Option<Arc<crate::live_call::LiveCallSupport>>,

    /// Optional retained-outcome choices and persistent P-base ancestry.
    /// Ordinary simulations pay no allocation for this opt-in model.
    pub coding_end_provenance: Option<Arc<super::CodingEndProvenance>>,
    /// Optional immutable repair founder and persistent descendant observation map.
    pub duplex_provenance: Option<Arc<super::DuplexProvenance>>,

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
    /// Cumulative realized SHM event burden per input-molecule nucleotide:
    /// sum of `applied_count / input_pool_len` over mutation stages.
    /// Observation loss and structural edits preserve earlier contributions.
    /// Lineage projection resets this with its net-difference counter.
    pub mutation_burden: f64,
}

impl Simulation {
    pub fn with_duplex_provenance(
        &self,
        provenance: Arc<super::DuplexProvenance>,
    ) -> Result<Self, String> {
        if self.duplex_provenance.is_some() {
            return Err("duplex founder ancestry is already installed".into());
        }
        provenance.validate_initial_pool(self.pool.as_slice())?;
        let mut next = self.clone();
        next.duplex_provenance = Some(provenance);
        Ok(next)
    }

    pub fn coding_end_choice(&self, segment: Segment) -> Option<&super::CodingEndChoice> {
        self.coding_end_provenance
            .as_ref()?
            .choices
            .iter()
            .rev()
            .find(|choice| choice.segment == segment)
    }

    pub fn with_coding_end_choice(
        &self,
        segment: Segment,
        allele_id: crate::refdata::AlleleId,
        five: i64,
        three: i64,
        outcome_index: usize,
    ) -> Self {
        let mut next = self.clone();
        let provenance = Arc::make_mut(
            next.coding_end_provenance
                .get_or_insert_with(|| Arc::new(super::CodingEndProvenance::default())),
        );
        provenance.choices.push(super::CodingEndChoice {
            choice_id: provenance.choices.len(),
            segment,
            allele_id,
            five,
            three,
            outcome_index,
        });
        next
    }

    pub fn live_call_alleles(&self, segment: Segment) -> Option<&[crate::refdata::AlleleId]> {
        self.live_call_support
            .as_ref()
            .and_then(|support| support.alleles(segment))
    }

    pub(crate) fn with_live_call_support(
        &self,
        support: Arc<crate::live_call::LiveCallSupport>,
    ) -> Self {
        let mut next = self.clone();
        if self.live_call_support.as_deref() != Some(support.as_ref()) {
            next.segment_calls = Arc::new(crate::live_call::SegmentCalls::empty());
        }
        next.live_call_support = Some(support);
        next
    }

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
            mutation_burden: 0.0,
            live_call_support: None,
            coding_end_provenance: None,
            duplex_provenance: None,
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
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            coding_end_provenance: self.coding_end_provenance.clone(),
            duplex_provenance: self.duplex_provenance.clone(),
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
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            coding_end_provenance: self.coding_end_provenance.clone(),
            duplex_provenance: self.duplex_provenance.clone(),
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
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            coding_end_provenance: self.coding_end_provenance.clone(),
            duplex_provenance: self.duplex_provenance.clone(),
        }
    }

    /// Replace the SHM burden when recording a mutation stage or synthesizing
    /// lineage net-difference provenance. This is not an observed mismatch rate.
    pub(crate) fn with_mutation_burden(&self, burden: f64) -> Self {
        let mut next = self.clone();
        next.mutation_burden = burden;
        next
    }

    #[must_use = "with_nucleotide_pushed returns (Simulation, NucHandle); \
                  both must be used or destructured. Drop the handle \
                  explicitly with `let (next, _) = ...` if it isn't needed."]
    pub fn with_nucleotide_pushed(&self, n: Nucleotide) -> (Self, NucHandle) {
        let (pool, h) = self.pool.with_pushed(n);
        let duplex_provenance = self
            .duplex_provenance
            .as_ref()
            .map(|p| Arc::new(p.with_inserted(h.index(), n)));
        (
            Self {
                pool,
                sequence: self.sequence.clone(),
                assignments: self.assignments,
                segment_calls: self.segment_calls.clone(),
                dirty_log: self.dirty_log.clone(),
                mutation_count: self.mutation_count,
                mutation_burden: self.mutation_burden,
                live_call_support: self.live_call_support.clone(),
                coding_end_provenance: self.coding_end_provenance.clone(),
                duplex_provenance,
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
        let duplex_provenance = self.duplex_provenance.as_ref().map(|p| {
            Arc::new(p.with_replacement(
                range.start,
                range.start,
                &pool.as_slice()[range.start as usize..range.end as usize],
                super::DuplexRemovalCause::Replacement,
            ))
        });
        (
            Self {
                pool,
                sequence: self.sequence.clone(),
                assignments: self.assignments,
                segment_calls: self.segment_calls.clone(),
                dirty_log: self.dirty_log.clone(),
                mutation_count: self.mutation_count,
                mutation_burden: self.mutation_burden,
                live_call_support: self.live_call_support.clone(),
                coding_end_provenance: self.coding_end_provenance.clone(),
                duplex_provenance,
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
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            coding_end_provenance: self.coding_end_provenance.clone(),
            duplex_provenance: self.duplex_provenance.as_ref().map(|p| {
                Arc::new(
                    p.with_substitution(
                        handle.index(),
                        *self
                            .pool
                            .get(handle)
                            .expect("nucleotide handle outside pool"),
                        new_n.base,
                    ),
                )
            }),
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
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            coding_end_provenance: self.coding_end_provenance.clone(),
            duplex_provenance: self.duplex_provenance.as_ref().map(|p| {
                Arc::new(
                    p.with_substitution(
                        handle.index(),
                        *self
                            .pool
                            .get(handle)
                            .expect("nucleotide handle outside pool"),
                        new_base,
                    ),
                )
            }),
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
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            coding_end_provenance: self.coding_end_provenance.clone(),
            duplex_provenance: self.duplex_provenance.clone(),
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
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            coding_end_provenance: self.coding_end_provenance.clone(),
            duplex_provenance: self.duplex_provenance.clone(),
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
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            coding_end_provenance: self.coding_end_provenance.clone(),
            duplex_provenance: self.duplex_provenance.clone(),
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
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            coding_end_provenance: self.coding_end_provenance.clone(),
            duplex_provenance: self.duplex_provenance.clone(),
        }
    }

    /// Return a new simulation with the orientation at `segment`
    /// updated. Mirror of [`Self::with_trim`]; persistent IR contract.
    ///
    /// Slice A wiring — `AssembleSegmentPass` does not yet read the
    /// orientation, so calling this setter is observationally a
    /// no-op on the assembled output. The setter exists so future
    /// slices (B onwards) have a stable surface to commit D
    /// inversion decisions through.
    pub fn with_allele_orientation(
        &self,
        segment: Segment,
        orientation: crate::assignment::SegmentOrientation,
    ) -> Self {
        Self {
            pool: self.pool.clone(),
            sequence: self.sequence.clone(),
            assignments: self.assignments.with_orientation(segment, orientation),
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            coding_end_provenance: self.coding_end_provenance.clone(),
            duplex_provenance: self.duplex_provenance.clone(),
        }
    }

    /// Return a new simulation with `n` inserted at pool position
    /// `at`. The persistent layer adjusts every region's `(start,
    /// end)` for the +1 shift; codon-rail metadata is no longer
    /// stored on `Region`, so callers that need it compute on
    /// demand via [`crate::ir::compute_codon_rail`].
    pub fn with_indel_inserted(&self, at: u32, n: Nucleotide) -> Self {
        let new_pool = self.pool.with_inserted(at, n);
        let adjusted = if self.duplex_provenance.is_some() {
            duplex_sequence_after_insertion(&self.sequence, at, n.segment)
        } else {
            self.sequence.with_indel_adjusted(at, 1)
        };
        Self {
            pool: new_pool,
            sequence: Arc::new(adjusted),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            duplex_provenance: self
                .duplex_provenance
                .as_ref()
                .map(|p| Arc::new(p.with_inserted(at, n))),
            coding_end_provenance: self
                .coding_end_provenance
                .as_ref()
                .map(|p| Arc::new(p.shifted_for_replacement(at, at, 1))),
        }
    }

    /// Return a new simulation with the **unique** region for
    /// `segment` structurally replaced by `replacement` bytes.
    ///
    /// The pool span the old region occupied is excised; `replacement`
    /// is spliced in at the same start position. Downstream regions
    /// shift by `replacement.len() - old_region.len()` via the
    /// [`PoolRange::after_segment_replacement`] kernel.
    ///
    /// Returns `(new_simulation, old_region, new_region)`. The
    /// caller (`SimulationBuilder::replace_segment`) uses the
    /// returned regions to compute `bytes_delta` and broadcast the
    /// matching `SimulationEvent::SegmentReplaced`.
    ///
    /// **Pre-conditions:** the sequence holds **exactly one**
    /// region for `segment`. Panics otherwise — receptor revision
    /// Slice A is defined for the heavy-chain V-only case where
    /// this invariant is structural; future slices (clonal
    /// expansion, repeat alleles) will revisit the API rather than
    /// loosen the gate here.
    ///
    /// The codon-rail recompute happens via
    /// `Sequence::with_frame_phases_recomputed` inside
    /// `with_segment_replaced`; consumers that need the translated
    /// rail compute on demand via [`crate::ir::compute_codon_rail`].
    #[must_use = "with_segment_replaced returns (Simulation, old Region, new Region); \
                  all three must be used or destructured."]
    pub fn with_segment_replaced(
        &self,
        segment: Segment,
        replacement: Vec<Nucleotide>,
    ) -> (Self, Region, Region) {
        let mut matches = self
            .sequence
            .regions
            .iter()
            .filter(|r| r.segment == segment);
        let old_region = matches
            .next()
            .cloned()
            .unwrap_or_else(|| {
                panic!(
                    "Simulation::with_segment_replaced: no region for {segment:?} in sequence with {} region(s)",
                    self.sequence.regions.len(),
                )
            });
        assert!(
            matches.next().is_none(),
            "Simulation::with_segment_replaced: multiple regions for {segment:?}; receptor revision is defined for single-region segments",
        );

        let old_range = PoolRange::new(old_region.start.index(), old_region.end.index());
        let new_len = replacement.len() as u32;
        let new_pool = self.pool.with_range_replaced(old_range, &replacement);

        let new_region = Region::new(
            segment,
            NucHandle::new(old_range.start),
            NucHandle::new(old_range.start + new_len),
        );

        let new_sequence =
            self.sequence
                .with_segment_replaced(segment, old_range, new_region.clone());
        // The new_region returned to the caller must reflect the
        // frame-phase recompute that `with_segment_replaced` runs;
        // fetch the post-recompute region back out so the event
        // payload matches what's installed in the sequence.
        let installed_new_region = new_sequence
            .regions
            .iter()
            .find(|r| r.segment == segment)
            .cloned()
            .expect("Sequence::with_segment_replaced preserves the segment's region");

        let new_sim = Self {
            pool: new_pool,
            sequence: Arc::new(new_sequence),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            duplex_provenance: self.duplex_provenance.as_ref().map(|p| {
                Arc::new(p.with_replacement(
                    old_range.start,
                    old_range.end,
                    &replacement,
                    super::DuplexRemovalCause::Replacement,
                ))
            }),
            coding_end_provenance: self.coding_end_provenance.as_ref().map(|p| {
                Arc::new(p.shifted_for_replacement(old_range.start, old_range.end, new_len))
            }),
        };

        (new_sim, old_region, installed_new_region)
    }

    /// Return a new simulation with the nucleotide at pool position
    /// `at` removed. Same range-shift semantics as
    /// [`Self::with_indel_inserted`].
    pub fn with_indel_deleted(&self, at: u32) -> Self {
        self.with_query_deleted(at, super::DuplexRemovalCause::StructuralDeletion)
    }

    /// Explicit missing read coverage, distinct from a structural deletion.
    pub fn with_observation_deleted(&self, at: u32) -> Self {
        self.with_query_deleted(at, super::DuplexRemovalCause::ObservationLoss)
    }

    fn with_query_deleted(&self, at: u32, cause: super::DuplexRemovalCause) -> Self {
        let new_pool = self.pool.with_deleted(at);
        let adjusted = self.sequence.with_indel_adjusted(at, -1);
        Self {
            pool: new_pool,
            sequence: Arc::new(adjusted),
            assignments: self.assignments,
            segment_calls: self.segment_calls.clone(),
            dirty_log: self.dirty_log.clone(),
            mutation_count: self.mutation_count,
            mutation_burden: self.mutation_burden,
            live_call_support: self.live_call_support.clone(),
            duplex_provenance: self.duplex_provenance.as_ref().map(|p| {
                Arc::new(
                    p.with_removed(
                        at,
                        *self
                            .pool
                            .get(NucHandle::new(at))
                            .expect("deletion outside pool"),
                        cause,
                    ),
                )
            }),
            coding_end_provenance: self
                .coding_end_provenance
                .as_ref()
                .map(|p| Arc::new(p.shifted_for_replacement(at, at + 1, 0))),
        }
    }
}

/// A duplex observation has an explicit partition, including zero-width
/// sources. Keep its declared insertion owner even at the outer query edge
/// or when every source has become empty. Ordinary IR range semantics stay
/// unchanged.
fn duplex_sequence_after_insertion(sequence: &Sequence, at: u32, segment: Segment) -> Sequence {
    let Some(owner) = sequence.regions.iter().position(|region| {
        region.segment == segment && region.start.index() <= at && at <= region.end.index()
    }) else {
        // Invalid low-level owners remain detectable by the layout validator.
        return sequence.with_indel_adjusted(at, 1);
    };
    Sequence {
        regions: sequence
            .regions
            .iter()
            .enumerate()
            .map(|(index, region)| {
                let start = region.start.index();
                let end = region.end.index();
                let (start, end) = if index == owner {
                    (start, end + 1)
                } else if at < start || (at == start && index > owner) {
                    (start + 1, end + 1)
                } else {
                    (start, end)
                };
                Region {
                    start: NucHandle::new(start),
                    end: NucHandle::new(end),
                    ..region.clone()
                }
            })
            .collect(),
    }
    .with_frame_phases_recomputed()
}
