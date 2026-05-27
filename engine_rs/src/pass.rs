//! Pass manager — orchestration types for executing a fixed sequence
//! of passes, recording choices, and preserving IR history.
//!
//! The biology lives in concrete pass implementations under
//! `crate::passes`. This module owns the shared execution contract:
//! metadata, context, errors, the `Pass` trait, ordered plans,
//! outcomes, and the direct runtime.

mod context;
mod error;
mod hook;
mod metadata;
mod outcome;
mod schedule;
mod support;
mod traits;

pub use context::PassContext;
pub use error::PassError;
pub use hook::{EffectHook, HookContext};
pub use metadata::{PassCompileEffect, PassEffect, PassRequirement};
pub use outcome::Outcome;
pub use schedule::{NodeId, Schedule, ScheduleError};
pub use support::{AlleleIdSupport, IntegerSupport, PassCompileFact};
pub use traits::Pass;

/// Test-only execution helpers. Hidden from rustdoc and intentionally
/// nested under `testing::` so the import path itself signals
/// non-production usage. A guard test in `testing::lockdown` walks
/// `src/` and asserts no production source file (anything outside a
/// `#[cfg(test)]` scope) names `pass::testing::PassRuntime`.
///
/// Production callers must go through
/// [`crate::compiled::CompiledSimulator`] (or `OwnedCompiledSimulator`),
/// which threads the full schedule topo-sort, `analyze_plan`
/// validation, `ReferenceMatchIndex`, and effect-hook pipeline that
/// `PassRuntime` deliberately omits.
#[doc(hidden)]
pub mod testing;

/// Back-compat alias for the dependency-graph schedule. The plan/list
/// terminology is being phased out in favour of [`Schedule`] — both
/// now name the same type.
pub type PassPlan = Schedule;

#[cfg(test)]
mod tests;
