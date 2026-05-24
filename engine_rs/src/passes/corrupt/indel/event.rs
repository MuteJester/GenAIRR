use crate::address;
use crate::ir::{flag, NucHandle, Nucleotide, Segment, Simulation};
use crate::pass::{PassContext, PassError};
use crate::passes::mutation_transaction::MutationTransaction;
use crate::trace::ChoiceValue;

use super::IndelPass;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum IndelEvent {
    Insertion { site: u32, base: u8 },
    Deletion { site: Option<u32> },
}

impl IndelPass {
    /// Best-effort segment provenance for an inserted observation base.
    ///
    /// Insertions are synthetic, but downstream metadata still needs the
    /// nucleotide's segment to agree with the sequence context it lands in.
    /// This mirrors `Sequence::with_indel_adjusted`: if the insertion point
    /// is inside a region, or exactly at the start of a following region,
    /// the inserted nucleotide belongs to that region. Outside all regions,
    /// fall back to the nearest nucleotide's segment.
    pub(super) fn insertion_segment(sim: &Simulation, at: u32) -> Segment {
        for region in &sim.sequence.regions {
            let start = region.start.index();
            let end = region.end.index();
            if start <= at && at < end {
                return region.segment;
            }
        }

        if at < sim.pool.len() as u32 {
            return sim
                .pool
                .get(NucHandle::new(at))
                .map(|n| n.segment)
                .unwrap_or(Segment::V);
        }

        if at > 0 {
            return sim
                .pool
                .get(NucHandle::new(at - 1))
                .map(|n| n.segment)
                .unwrap_or(Segment::V);
        }

        Segment::V
    }

    pub(super) fn record_event(ctx: &mut PassContext, index: u32, event: IndelEvent) {
        match event {
            IndelEvent::Insertion { site, base } => {
                ctx.trace
                    .record(address::corrupt_indel_kind(index), ChoiceValue::Bool(true));
                ctx.trace.record(
                    address::corrupt_indel_site(index),
                    ChoiceValue::Int(site as i64),
                );
                ctx.trace
                    .record(address::corrupt_indel_base(index), ChoiceValue::Base(base));
            }
            IndelEvent::Deletion { site } => {
                ctx.trace
                    .record(address::corrupt_indel_kind(index), ChoiceValue::Bool(false));
                ctx.trace.record(
                    address::corrupt_indel_site(index),
                    ChoiceValue::Int(site.map(i64::from).unwrap_or(-1)),
                );
            }
        }
    }

    /// Apply an indel event through the active `MutationTransaction`,
    /// emitting `on_indel_inserted` / `on_indel_deleted` events to
    /// attached observers.
    ///
    /// Default observer impls are no-ops, so observer state becomes
    /// stale after the call — the post-pass
    /// `PassEffect::StructuralIndel` refresh handles re-derivation
    /// from scratch.
    pub(super) fn apply_event_via_tx(
        tx: &mut MutationTransaction<'_, '_>,
        event: IndelEvent,
    ) -> Result<(), PassError> {
        match event {
            IndelEvent::Insertion { site, base } => {
                let segment = Self::insertion_segment(tx.peek(), site);
                let new_nuc = Nucleotide::synthetic(base, segment, flag::INDEL_INSERTED);
                tx.insert_base(site, new_nuc)
            }
            IndelEvent::Deletion { site: Some(site) } => tx.delete_base(site).map(|_| ()),
            IndelEvent::Deletion { site: None } => Ok(()),
        }
    }
}
