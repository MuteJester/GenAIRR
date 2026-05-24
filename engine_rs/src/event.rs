//! Committed simulation event ledger.
//!
//! Events are not a pub/sub bus and they are not emitted directly by
//! individual passes. The compiled simulator creates one event only
//! after a pass has completed, strict fences have passed, and the
//! trace delta/state revision are ready to commit atomically.

use crate::ir::Simulation;
use crate::pass::PassEffect;

/// Stable class of a committed event.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum EventKind {
    /// One pass committed a state transition.
    PassCommitted,
}

/// Half-open trace range covered by one committed event.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct TraceSpan {
    pub start: usize,
    pub end: usize,
}

impl TraceSpan {
    pub fn new(start: usize, end: usize) -> Self {
        assert!(
            start <= end,
            "TraceSpan start must be <= end (got {}..{})",
            start,
            end
        );
        Self { start, end }
    }

    pub fn len(&self) -> usize {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }
}

/// Lightweight state summary recorded before and after a committed pass.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct StateSummary {
    pub pool_len: usize,
    pub region_count: usize,
    pub assigned_allele_count: usize,
}

impl StateSummary {
    pub fn from_simulation(sim: &Simulation) -> Self {
        let assigned_allele_count = sim.assignments.iter().count();

        Self {
            pool_len: sim.pool.len(),
            region_count: sim.sequence.region_count(),
            assigned_allele_count,
        }
    }
}

/// One durable event in the committed run ledger.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EventRecord {
    pub pass_index: usize,
    pub pass_name: String,
    pub kind: EventKind,
    pub effects: Vec<PassEffect>,
    pub trace_span: TraceSpan,
    pub pre: StateSummary,
    pub post: StateSummary,
}

impl EventRecord {
    pub fn pass_committed(
        pass_index: usize,
        pass_name: impl Into<String>,
        effects: Vec<PassEffect>,
        trace_span: TraceSpan,
        pre: StateSummary,
        post: StateSummary,
    ) -> Self {
        Self {
            pass_index,
            pass_name: pass_name.into(),
            kind: EventKind::PassCommitted,
            effects,
            trace_span,
            pre,
            post,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{NucFlags, Nucleotide, Segment};

    #[test]
    fn trace_span_len_uses_half_open_interval() {
        let span = TraceSpan::new(3, 8);
        assert_eq!(span.len(), 5);
        assert!(!span.is_empty());

        let empty = TraceSpan::new(4, 4);
        assert_eq!(empty.len(), 0);
        assert!(empty.is_empty());
    }

    #[test]
    #[should_panic(expected = "TraceSpan start must be <= end")]
    fn trace_span_rejects_inverted_range() {
        let _ = TraceSpan::new(9, 1);
    }

    #[test]
    fn state_summary_reflects_lightweight_simulation_shape() {
        let sim = Simulation::new();
        assert_eq!(
            StateSummary::from_simulation(&sim),
            StateSummary {
                pool_len: 0,
                region_count: 0,
                assigned_allele_count: 0,
            }
        );

        let (sim, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
            b'A',
            Segment::Np1,
            NucFlags::empty(),
        ));
        assert_eq!(StateSummary::from_simulation(&sim).pool_len, 1);
    }
}
