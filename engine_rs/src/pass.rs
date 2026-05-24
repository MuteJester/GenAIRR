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
mod runtime;
mod schedule;
mod support;
mod traits;

pub use context::PassContext;
pub use error::PassError;
pub use hook::{EffectHook, HookContext};
pub use metadata::{PassEffect, PassRequirement};
pub use outcome::Outcome;
pub use runtime::PassRuntime;
pub use schedule::{NodeId, Schedule, ScheduleError};
pub use support::{AlleleIdSupport, IntegerSupport, PassCompileFact};
pub use traits::Pass;

/// Back-compat alias for the dependency-graph schedule. The plan/list
/// terminology is being phased out in favour of [`Schedule`] — both
/// now name the same type.
pub type PassPlan = Schedule;

#[cfg(test)]
mod tests;
