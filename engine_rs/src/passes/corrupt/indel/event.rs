use crate::address;
use crate::ir::{flag, NucHandle, Nucleotide, Segment, Simulation};
use crate::pass::{PassContext, PassError};
use crate::passes::mutation_transaction::MutationTransaction;
use crate::trace::ChoiceValue;

use super::IndelPass;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum IndelEvent {
    Insertion {
        site: u32,
        base: u8,
    },
    Deletion {
        site: Option<u32>,
    },
    /// Permissive-mode reduce-and-skip: the count slot is consumed
    /// but no structural change is committed. Used when an active
    /// contract bundle rejects every (kind, site, base) candidate
    /// and the runtime is in permissive mode. Honors the contract
    /// rather than falling through to unconstrained sampling.
    NoOp,
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
                ctx.trace.record_choice(
                    address::ChoiceAddress::CorruptIndelKind(index),
                    ChoiceValue::Bool(true),
                );
                ctx.trace.record_choice(
                    address::ChoiceAddress::CorruptIndelSite(index),
                    ChoiceValue::Int(site as i64),
                );
                ctx.trace.record_choice(
                    address::ChoiceAddress::CorruptIndelBase(index),
                    ChoiceValue::Base(base),
                );
            }
            IndelEvent::Deletion { site } => {
                ctx.trace.record_choice(
                    address::ChoiceAddress::CorruptIndelKind(index),
                    ChoiceValue::Bool(false),
                );
                ctx.trace.record_choice(
                    address::ChoiceAddress::CorruptIndelSite(index),
                    ChoiceValue::Int(site.map(i64::from).unwrap_or(-1)),
                );
            }
            IndelEvent::NoOp => {
                // Permissive-mode reduce-and-skip. Record the site
                // as -1 (matching the "no-site" sentinel used by an
                // empty-pool deletion) so trace consumers can see
                // the slot was consumed without a structural change.
                ctx.trace.record_choice(
                    address::ChoiceAddress::CorruptIndelKind(index),
                    ChoiceValue::Bool(false),
                );
                ctx.trace.record_choice(
                    address::ChoiceAddress::CorruptIndelSite(index),
                    ChoiceValue::Int(-1),
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
            // Permissive-mode reduce-and-skip: the contract bundle
            // rejected every (kind, site, base) candidate for this
            // slot. Honor the contract by committing no structural
            // change; the trace records the slot via record_event so
            // consumers can audit the reduction.
            IndelEvent::NoOp => Ok(()),
        }
    }
}
