//! Validate the execution IR partition without inferring source ownership from
//! DNA. Empty source spans remain structural boundaries after loss or deletion.
use crate::ir::{Segment, Simulation};

mod history;

pub(crate) fn validate_layout(sim: &Simulation) -> Result<(), String> {
    let rank = |segment| match segment {
        Segment::V => 0,
        Segment::Np1 => 1,
        Segment::D => 2,
        Segment::Np2 => 3,
        Segment::J => 4,
    };
    let mut seen = [false; 5];
    let mut previous = None;
    let mut cursor = 0usize;
    for region in &sim.sequence.regions {
        let position = rank(region.segment);
        if seen[position] || previous.is_some_and(|p| p >= position) {
            return Err("duplex structural regions are duplicated or out of source order".into());
        }
        seen[position] = true;
        previous = Some(position);
        let (start, end) = (region.start.index() as usize, region.end.index() as usize);
        if start != cursor || start > end || end > sim.pool.len() {
            return Err("duplex structural regions do not partition the current query".into());
        }
        if sim.pool.as_slice()[start..end]
            .iter()
            .any(|n| n.segment != region.segment)
        {
            return Err("duplex nucleotide placement disagrees with its structural region".into());
        }
        cursor = end;
    }
    if !seen[0] || !seen[2] || !seen[4] {
        return Err("duplex structural source region is missing".into());
    }
    if cursor != sim.pool.len() {
        return Err("duplex structural regions do not cover the current query".into());
    }
    history::validate_history(sim)
}
