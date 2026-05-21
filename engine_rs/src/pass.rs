//! Pass manager — orchestration types for executing a fixed sequence
//! of passes, recording choices, and preserving IR history.
//!
//! The biology lives in concrete pass implementations under
//! `crate::passes`. This module owns the shared execution contract:
//! metadata, context, errors, the `Pass` trait, ordered plans,
//! outcomes, and the direct runtime.

mod context;
mod error;
mod metadata;
mod outcome;
mod plan;
mod runtime;
mod support;
mod traits;

pub use context::PassContext;
pub use error::PassError;
pub use metadata::{PassEffect, PassRequirement};
pub use outcome::Outcome;
pub use plan::PassPlan;
pub use runtime::PassRuntime;
pub use support::{AlleleIdSupport, IntegerSupport, PassCompileFact};
pub use traits::Pass;

#[cfg(test)]
mod tests;
