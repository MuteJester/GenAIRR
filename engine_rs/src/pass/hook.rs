//! `EffectHook` — runtime-managed reactions to pass effects.
//!
//! In Bevy, when a system mutates a `Component`, the schedule
//! invokes the next system that declared `Changed<Component>` — the
//! cross-system coupling is data, not code. The engine equivalent:
//! after each pass commits, the runtime iterates registered hooks
//! and gives each one the freshly-emitted [`PassEffect`] list. Hooks
//! decide for themselves whether to react and how.
//!
//! Stage 2 of the dep-graph scheduler migration moves the previously
//! hand-coded `apply_live_call_updates` `match` (in `compiled/execute.rs`)
//! into a single [`crate::live_call::LiveCallRefreshHook`]. The
//! cross-segment biology rules ("assembling D refreshes V", "indel
//! refreshes V/D/J unconditionally") now live next to the rest of
//! the live-call module, where they're discoverable by the people
//! editing that module.
//!
//! Adding a future derived-state consumer — e.g. a junction cache,
//! a per-allele evidence summariser — is one new file implementing
//! `EffectHook`, plus one `compile()`-site registration. No edits to
//! the executor.

use crate::ir::Simulation;
use crate::live_call::ReferenceMatchIndex;

use super::PassEffect;

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

/// A runtime-registered reaction to `PassEffect`s.
///
/// One hook instance lives in the `CompiledSimulator`. After each
/// pass commits, the executor calls `apply` on every hook with the
/// freshly-emitted effects. Hooks return the (possibly updated)
/// simulation; multiple hooks chain in registration order.
///
/// Hooks decide internally whether to react based on the `effects`
/// slice. Returning the input simulation unchanged is a no-op.
pub trait EffectHook: Send + Sync {
    /// Stable, human-readable identifier for the hook. Used in
    /// diagnostics and tests.
    fn name(&self) -> &str;

    /// React to the effects emitted by the just-committed pass.
    ///
    /// Hooks may inspect `effects` to decide whether to do work
    /// (matching the structure of the previous `apply_live_call_updates`
    /// `match`). They may also inspect the input simulation's state
    /// (e.g. `sim.live_calls.dirty_windows`) for finer-grained
    /// dispatch.
    fn apply(&self, sim: Simulation, effects: &[PassEffect], ctx: HookContext) -> Simulation;
}
