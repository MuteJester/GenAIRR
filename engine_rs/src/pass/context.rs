use crate::contract::ContractSet;
use crate::feasibility::FeasibilityContext;
use crate::ir::SimulationEvent;
use crate::live_call::ReferenceMatchIndex;
use crate::refdata::RefDataConfig;
use crate::replay::TraceCursor;
use crate::rng::Rng;
use crate::trace::Trace;

/// Context handed to every `Pass::execute` invocation.
#[non_exhaustive]
pub struct PassContext<'a> {
    pub trace: &'a mut Trace,
    pub rng: &'a mut Rng,
    /// Zero-based index of the pass currently executing inside the plan.
    pub pass_index: usize,
    pub refdata: Option<&'a RefDataConfig>,
    /// Active contract set. `None` means no contracts are active.
    pub contracts: Option<&'a ContractSet>,
    /// Compile-derived downstream feasibility domains.
    pub feasibility: Option<&'a FeasibilityContext>,
    /// Compile-derived per-segment reference index (k-mer + base-pos
    /// bitsets) used by the live-call walker. `Some` only when the
    /// pass is running through the compiled execution path; the
    /// test-only `PassRuntime::execute_with_refdata` leaves this
    /// `None`, in which case observer-aware passes fall back to their
    /// batch-push code paths.
    pub reference_index: Option<&'a ReferenceMatchIndex>,
    /// Active trace-injected replay cursor. `Some(...)` puts the
    /// runtime into **consume mode** (Option B): each sampling site
    /// must consume the next pre-recorded `ChoiceRecord` from the
    /// cursor instead of drawing from `rng`. `None` is the normal
    /// fresh-RNG path. The cursor is mutated as each consumption
    /// advances its position, hence the `&mut`.
    pub replay_cursor: Option<&'a mut TraceCursor>,
    /// Caller-supplied capture buffer for [`SimulationEvent`]s fired
    /// from a pass's internal [`crate::ir::SimulationBuilder`]s. When
    /// `Some(&mut buf)`, migrated passes attach an
    /// [`crate::ir::event_log_observer::EventLogObserver`] to their
    /// internal builder(s), drain it after seal, and append the
    /// captured events to `buf` in emission order. Used by
    /// pass-level event-log tests; production callers pass `None`
    /// and the per-pass observer attachment becomes a zero-cost
    /// no-op.
    pub event_log_sink: Option<&'a mut Vec<SimulationEvent>>,
}
