//! `EffectHook` — runtime-managed reactions to a just-committed
//! pass.
//!
//! In Bevy, when a system mutates a `Component`, the schedule
//! invokes the next system that declared `Changed<Component>` — the
//! cross-system coupling is data, not code. The engine equivalent:
//! after each pass commits, the runtime iterates registered hooks
//! and gives each one both the pass's static
//! [`PassCompileEffect`] declarations **and** the runtime
//! [`crate::ir::SimulationEvent`] stream it emitted. Hooks decide
//! for themselves what to react to.
//!
//! Today's only consumer ([`crate::live_call::LiveCallRefreshHook`])
//! reads the event stream — **not** the compile-effect list — to
//! decide which V/D/J segments to refresh. That inversion is the
//! load-bearing architectural commitment: runtime derived-state
//! refresh follows what *actually happened* (the event stream),
//! not what a pass *claimed it would do* (the compile-effect
//! declaration).
//!
//! Adding a future derived-state consumer — e.g. a junction cache,
//! a per-allele evidence summariser — is one new file implementing
//! `EffectHook`, plus one `compile()`-site registration. No edits to
//! the executor. New hooks should prefer reading
//! `simulation_events`; the `compile_effects` parameter is retained
//! on the trait signature for hooks that legitimately need the
//! declarative surface (e.g. a schedule-aware diagnostic).

use crate::ir::{Simulation, SimulationEvent};
use crate::live_call::ReferenceMatchIndex;

use super::PassCompileEffect;

/// Borrowed context handed to every [`EffectHook::apply`] call.
///
/// Extensible by adding fields as new hooks need more inputs.
/// Today's only consumer (the live-call refresh) needs the
/// reference match index; future hooks may need refdata, contracts,
/// or the trace.
#[derive(Copy, Clone)]
pub struct HookContext<'a> {
    /// Compile-time reference-data index used by the live-call
    /// walker. `None` for test paths that compile without refdata.
    pub reference_index: Option<&'a ReferenceMatchIndex>,
}

/// A runtime-registered reaction to a just-committed pass.
///
/// One hook instance lives in the `CompiledSimulator`. After each
/// pass commits, the executor calls `apply` on every hook with the
/// pass's freshly-emitted [`PassEffect`]s **and** its
/// [`SimulationEvent`] stream. Hooks return the (possibly updated)
/// simulation; multiple hooks chain in registration order.
///
/// ## `compile_effects` vs. `events` — two complementary surfaces
///
/// - `compile_effects` are **static, compile-time declarations**
///   from [`crate::pass::Pass::effects`] (see
///   [`PassCompileEffect`]). They flow into the schedule analyzer
///   for ordering / dependency reasoning. Hooks that only need a
///   coarse "this category of pass happened" cue may still read
///   them — but they describe *intent*, not what happened.
/// - `events` is the **runtime consequence stream** the pass
///   actually emitted via its internal builders. This is the
///   source of truth for derived-state refresh: a hook reacting
///   to `BaseChanged` events knows exactly which sites were
///   edited, while a hook reading `PassCompileEffect::EditBases`
///   only knows a pass *intended* to edit.
///
/// New hooks should prefer `events`. The `compile_effects`
/// parameter stays on the signature so future schedule-only hooks
/// (e.g. a "count how many passes declared `TrimAllele`"
/// diagnostic) can still see the declarative surface.
pub trait EffectHook: Send + Sync {
    /// Stable, human-readable identifier for the hook. Used in
    /// diagnostics and tests.
    fn name(&self) -> &str;

    /// React to the just-committed pass.
    ///
    /// `compile_effects` are the pass's static declarations (kept
    /// for scheduling-aware hooks). `events` is the ordered stream
    /// of consequences the pass actually emitted via its internal
    /// builders — the source of truth for derived-state refresh.
    fn apply(
        &self,
        sim: Simulation,
        compile_effects: &[PassCompileEffect],
        events: &[SimulationEvent],
        ctx: HookContext,
    ) -> Simulation;
}
