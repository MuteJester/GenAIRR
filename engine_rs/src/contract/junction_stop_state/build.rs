use crate::ir::{NucHandle, Segment, Simulation};
use crate::refdata::RefDataConfig;

use super::layout::{
    append_d_body, append_j_head_from_refdata, append_np_slots, precompute_static_violation,
    retained_bounds,
};
use super::model::{JunctionSlot, SlotSource};
use super::JunctionStopState;

impl JunctionStopState {
    /// Build state for an NP draw loop.
    ///
    /// `np_segment` is the NP region currently being drawn (`Np1`
    /// or `Np2`). The state's junction layout truncates at the
    /// boundary the slow path would: an NP1 draw doesn't include
    /// `J-head` (NP2 isn't drawn yet, codons past the candidate
    /// boundary are mostly incomplete), but an NP2 draw includes
    /// everything because NP1 is already committed.
    ///
    /// `np_total_len` is the total planned NP draw count for the
    /// pass — equivalent to `ChoiceContext::draw_count` for the
    /// individual candidate calls.
    ///
    /// Returns `None` when the slow path's `hypothetical_junction_bases`
    /// would also return `None` (missing V/J assignment, anchor,
    /// retained bounds violation). The pass falls back to the slow
    /// path in that case.
    pub fn build(
        sim: &Simulation,
        refdata: &RefDataConfig,
        np_segment: Segment,
        np_total_len: u32,
    ) -> Option<Self> {
        let v_inst = sim.assignments.get(Segment::V).copied()?;
        let j_inst = sim.assignments.get(Segment::J).copied()?;
        let v_allele = refdata.get(Segment::V, v_inst.allele_id)?;
        let j_allele = refdata.get(Segment::J, j_inst.allele_id)?;
        let v_anchor = v_allele.anchor? as u32;
        let j_anchor = j_allele.anchor? as u32;

        let (v_retained_start, v_retained_end) =
            retained_bounds(v_allele.len(), v_inst.trim_5, v_inst.trim_3)?;
        let (j_retained_start, j_retained_end) =
            retained_bounds(j_allele.len(), j_inst.trim_5, j_inst.trim_3)?;
        if v_anchor < v_retained_start || v_anchor >= v_retained_end {
            return None;
        }
        if j_anchor < j_retained_start {
            return None;
        }

        let v_region = sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)?;
        let v_anchor_pool = v_region.start.index() + (v_anchor - v_retained_start);
        let v_end_pool = v_region.end.index();
        if v_anchor_pool > v_end_pool {
            return None;
        }

        let has_d = sim.assignments.has(Segment::D);
        let mut slots: Vec<JunctionSlot> = Vec::new();
        let mut fixed_bytes: Vec<u8> = Vec::new();

        for pos in v_anchor_pool..v_end_pool {
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

        let np1_handled = append_np_slots(&mut slots, sim, Segment::Np1, np_segment, np_total_len);
        if !np1_handled {
            let static_violation = precompute_static_violation(&slots, &fixed_bytes);
            return Some(Self {
                slots,
                fixed_bytes,
                static_violation,
            });
        }

        if has_d {
            append_d_body(&mut slots, &mut fixed_bytes, sim, refdata)?;

            let np2_handled =
                append_np_slots(&mut slots, sim, Segment::Np2, np_segment, np_total_len);
            if !np2_handled {
                let static_violation = precompute_static_violation(&slots, &fixed_bytes);
                return Some(Self {
                    slots,
                    fixed_bytes,
                    static_violation,
                });
            }
        }

        let j_junction_end = j_anchor.saturating_add(3).min(j_retained_end);
        append_j_head_from_refdata(
            &mut slots,
            &mut fixed_bytes,
            j_allele,
            j_anchor,
            j_junction_end,
            j_retained_start,
        );

        let static_violation = precompute_static_violation(&slots, &fixed_bytes);

        Some(Self {
            slots,
            fixed_bytes,
            static_violation,
        })
    }
}
