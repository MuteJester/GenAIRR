use crate::contract::ContractSet;
use crate::feasibility::FeasibilityContext;
use crate::live_call::ReferenceMatchIndex;
use crate::refdata::RefDataConfig;
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
}
