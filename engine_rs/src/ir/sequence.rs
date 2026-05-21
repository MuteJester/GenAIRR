use super::{NucHandle, Region};

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

    /// Return a new sequence with every region's range adjusted for an indel.
    pub fn with_indel_adjusted(&self, pos: u32, delta: i32) -> Self {
        let new_regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                let r_start = r.start.index();
                let r_end = r.end.index();

                if r_end <= pos {
                    r.clone()
                } else if r_start > pos {
                    Region {
                        start: NucHandle::new(shift_pos(r_start, delta)),
                        end: NucHandle::new(shift_pos(r_end, delta)),
                        ..r.clone()
                    }
                } else {
                    Region {
                        end: NucHandle::new(shift_pos(r_end, delta)),
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

fn shift_pos(pos: u32, delta: i32) -> u32 {
    if delta >= 0 {
        pos.saturating_add(delta as u32)
    } else {
        pos.saturating_sub((-delta) as u32)
    }
}
