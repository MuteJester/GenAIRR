use crate::ir::{translate_codon, NucHandle, Segment, Simulation, AMINO_STOP};
use crate::refdata::{Allele, RefDataConfig};

use super::model::{JunctionSlot, SlotSource, StaticViolation};

pub(super) fn np_segment_label(segment: Segment) -> &'static str {
    match segment {
        Segment::Np1 => "np1",
        Segment::Np2 => "np2",
        _ => "np?",
    }
}

pub(super) fn retained_bounds(allele_len: u32, trim_5: u16, trim_3: u16) -> Option<(u32, u32)> {
    let trim_5 = trim_5 as u32;
    let trim_3 = trim_3 as u32;
    if trim_5 + trim_3 > allele_len {
        return None;
    }
    Some((trim_5, allele_len - trim_3))
}

/// Append the NP region's slots into `slots`. Returns `true` when
/// the NP region is fully represented (either already in the pool
/// or this is the active draw segment with a known total length).
/// Returns `false` when the slow path would also terminate the
/// junction-buffer construction here.
pub(super) fn append_np_slots(
    slots: &mut Vec<JunctionSlot>,
    sim: &Simulation,
    np_segment: Segment,
    cand_segment: Segment,
    cand_total_len: u32,
) -> bool {
    if let Some(region) = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == np_segment)
    {
        let region_start = region.start.index();
        let region_len = region.end.index().saturating_sub(region_start);
        for i in 0..region_len {
            slots.push(JunctionSlot {
                source: SlotSource::Np {
                    segment: np_segment,
                    np_index: i,
                    pool_pos: region_start + i,
                },
            });
        }
        return true;
    }

    if cand_segment == np_segment {
        let base = sim.pool.len() as u32;
        for i in 0..cand_total_len {
            slots.push(JunctionSlot {
                source: SlotSource::Np {
                    segment: np_segment,
                    np_index: i,
                    pool_pos: base + i,
                },
            });
        }
        true
    } else {
        false
    }
}

pub(super) fn append_d_body(
    slots: &mut Vec<JunctionSlot>,
    fixed_bytes: &mut Vec<u8>,
    sim: &Simulation,
    refdata: &RefDataConfig,
) -> Option<()> {
    let d_inst = sim.assignments.get(Segment::D).copied()?;
    let d_allele = refdata.get(Segment::D, d_inst.allele_id)?;
    let (d_retained_start, d_retained_end) =
        retained_bounds(d_allele.len(), d_inst.trim_5, d_inst.trim_3)?;
    if d_retained_end <= d_retained_start {
        return Some(());
    }

    if let Some(region) = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::D)
    {
        let pool_len = d_retained_end.saturating_sub(d_retained_start);
        let pool_end = region
            .start
            .index()
            .saturating_add(pool_len)
            .min(region.end.index());
        for pos in region.start.index()..pool_end {
            let byte = sim
                .pool
                .get(NucHandle::new(pos))
                .map(|n| n.base)
                .unwrap_or(b'N');
            let fixed_index = fixed_bytes.len() as u32;
            fixed_bytes.push(byte);
            slots.push(JunctionSlot {
                source: SlotSource::Fixed { fixed_index },
            });
        }
    } else {
        for pos in d_retained_start..d_retained_end {
            let byte = d_allele.seq.get(pos as usize).copied().unwrap_or(b'N');
            let fixed_index = fixed_bytes.len() as u32;
            fixed_bytes.push(byte);
            slots.push(JunctionSlot {
                source: SlotSource::Fixed { fixed_index },
            });
        }
    }

    Some(())
}

pub(super) fn append_j_head_from_refdata(
    slots: &mut Vec<JunctionSlot>,
    fixed_bytes: &mut Vec<u8>,
    j_allele: &Allele,
    j_anchor: u32,
    j_junction_end: u32,
    j_retained_start: u32,
) {
    if j_anchor >= j_junction_end {
        return;
    }

    // Match the slow path's append_fixed_segment for J:
    // allele.seq[j_trim_5..min(j_anchor+3, j_retained_end)].
    for pos in j_retained_start..j_junction_end {
        let byte = j_allele.seq.get(pos as usize).copied().unwrap_or(b'N');
        let fixed_index = fixed_bytes.len() as u32;
        fixed_bytes.push(byte);
        slots.push(JunctionSlot {
            source: SlotSource::Fixed { fixed_index },
        });
    }
}

pub(super) fn precompute_static_violation(
    slots: &[JunctionSlot],
    fixed_bytes: &[u8],
) -> Option<StaticViolation> {
    let n = slots.len();
    let mut off = 0usize;
    while off + 3 <= n {
        let s0 = slots[off].source;
        let s1 = slots[off + 1].source;
        let s2 = slots[off + 2].source;
        if let (
            SlotSource::Fixed { fixed_index: i0 },
            SlotSource::Fixed { fixed_index: i1 },
            SlotSource::Fixed { fixed_index: i2 },
        ) = (s0, s1, s2)
        {
            let b0 = fixed_bytes[i0 as usize];
            let b1 = fixed_bytes[i1 as usize];
            let b2 = fixed_bytes[i2 as usize];
            if translate_codon(b0, b1, b2) == AMINO_STOP {
                return Some(StaticViolation {
                    offset: off as u32,
                    codon: [b0, b1, b2],
                });
            }
        }
        off += 3;
    }
    None
}
