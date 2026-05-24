//! Test-only execution helpers. See the module-level doc on
//! [`crate::pass::testing`] for the rationale and production-path
//! pointers.

mod runtime;

pub use runtime::PassRuntime;

#[cfg(test)]
mod lockdown;
