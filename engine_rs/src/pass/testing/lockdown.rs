//! Guard: no production source names `pass::testing::PassRuntime`.
//!
//! The reviewer flagged `PassRuntime` as a foot-gun — a runtime
//! entrypoint that bypasses `analyze_plan`, schedule topo-sort,
//! effect hooks, and `ReferenceMatchIndex`. Today every caller is
//! inside a `#[cfg(test)]` scope or under `engine_rs/tests/`, but
//! visibility can't enforce that on its own (`pub` is required so
//! integration tests can import the symbol).
//!
//! The source walk, the test-scope detection and the test-file
//! conventions live in [`crate::source_guard`], shared with the
//! raw-seed guard in [`crate::rng::discipline`]. This file is the
//! policy: what is looked for, and which files may name it anyway.

use crate::source_guard::mentions_in_production;

const REFERENCE: &str = "pass::testing::PassRuntime";

/// Files allowed to name the path unconditionally: the module's own
/// files plus the public re-export site.
const MODULE_INTERNAL: &[&str] = &[
    "src/pass.rs",
    "src/pass/testing/mod.rs",
    "src/pass/testing/runtime.rs",
    "src/pass/testing/lockdown.rs",
];

#[test]
fn no_production_source_references_pass_testing_passruntime() {
    let violations: Vec<String> = mentions_in_production(REFERENCE)
        .into_iter()
        .filter(|m| !MODULE_INTERNAL.contains(&m.file.as_str()))
        .map(|m| m.to_string())
        .collect();

    assert!(
        violations.is_empty(),
        "Production source files referencing `pass::testing::PassRuntime` outside a \
         `#[cfg(test)]` scope. Use `CompiledSimulator` in production code. Offenders:\n{}",
        violations.join("\n")
    );
}
