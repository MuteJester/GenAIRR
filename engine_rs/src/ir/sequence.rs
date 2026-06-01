use super::{NucHandle, PoolRange, Region, Segment};

/// The assembled sequence — the root structural entity at the
/// biological product level.
#[derive(Clone, Debug, Default)]
pub struct Sequence {
    /// Regions in biological assembly order.
    pub regions: Vec<Region>,
}

impl Sequence {
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of regions currently in the sequence.
    pub fn region_count(&self) -> usize {
        self.regions.len()
    }

    /// Return a new sequence with `region` appended.
    pub fn with_region_added(&self, region: Region) -> Self {
        let mut next = self.clone();
        next.regions.push(region);
        next
    }

    /// Return a new sequence with the region at `idx` replaced by `region`.
    pub fn with_region_replaced(&self, idx: usize, region: Region) -> Self {
        let mut next = self.clone();
        next.regions[idx] = region;
        next
    }

    /// Return a new sequence whose regions have frame phases recomputed.
    pub fn with_frame_phases_recomputed(&self) -> Self {
        let mut cumulative_len = 0u64;
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                let next = r.with_frame_phase((cumulative_len % 3) as u8);
                cumulative_len = cumulative_len.saturating_add(next.len() as u64);
                next
            })
            .collect();
        Self { regions }
    }

    /// Return a new sequence with every region's range adjusted for
    /// a single-nucleotide indel at `pos`. `delta` must be `+1`
    /// (insertion) or `-1` (deletion); the per-region shift rule
    /// lives in [`PoolRange::after_insertion`] / [`PoolRange::after_deletion`].
    pub fn with_indel_adjusted(&self, pos: u32, delta: i32) -> Self {
        debug_assert!(
            delta == 1 || delta == -1,
            "with_indel_adjusted only supports ±1 deltas, got {delta}"
        );
        let new_regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                let range = PoolRange::new(r.start.index(), r.end.index());
                let shifted = if delta == 1 {
                    range.after_insertion(pos)
                } else {
                    range.after_deletion(pos)
                };
                Region {
                    start: NucHandle::new(shifted.start),
                    end: NucHandle::new(shifted.end),
                    ..r.clone()
                }
            })
            .collect();
        Sequence {
            regions: new_regions,
        }
        .with_frame_phases_recomputed()
    }

    /// Return a new sequence with the single region for `segment`
    /// replaced by `new_region`, and every other region's range
    /// adjusted via [`PoolRange::after_segment_replacement`] for
    /// the swap of `replaced` for the bytes now occupying
    /// `[replaced.start, replaced.start + new_region.len())`.
    ///
    /// Panics if `segment` has zero or multiple regions; receptor
    /// revision Slice A is defined for single-region segments and
    /// the [`crate::ir::Simulation::with_segment_replaced`] gate
    /// enforces this before dispatching here.
    pub fn with_segment_replaced(
        &self,
        segment: Segment,
        replaced: PoolRange,
        new_region: Region,
    ) -> Self {
        let matches = self
            .regions
            .iter()
            .filter(|r| r.segment == segment)
            .count();
        assert!(
            matches == 1,
            "Sequence::with_segment_replaced: expected exactly 1 region for {segment:?}, found {matches}",
        );
        let new_len = new_region.len();
        let new_regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                if r.segment == segment {
                    new_region.clone()
                } else {
                    let range = PoolRange::new(r.start.index(), r.end.index());
                    let shifted = range.after_segment_replacement(replaced, new_len);
                    Region {
                        start: NucHandle::new(shifted.start),
                        end: NucHandle::new(shifted.end),
                        ..r.clone()
                    }
                }
            })
            .collect();
        Sequence {
            regions: new_regions,
        }
        .with_frame_phases_recomputed()
    }
}
