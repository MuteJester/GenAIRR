//! Compile-time-ish guard: assert that no production source file
//! names `pass::testing::PassRuntime`.
//!
//! The reviewer flagged `PassRuntime` as a foot-gun — a runtime
//! entrypoint that bypasses `analyze_plan`, schedule topo-sort,
//! effect hooks, and `ReferenceMatchIndex`. Today every caller is
//! inside a `#[cfg(test)]` scope or under `engine_rs/tests/`, but
//! visibility can't enforce that on its own (`pub` is required so
//! integration tests can import the symbol).
//!
//! This test walks `src/` at test time, finds every reference to
//! `pass::testing::PassRuntime`, and asserts each one sits inside a
//! `#[cfg(test)] mod ... { ... }` scope via brace tracking. Files
//! entirely under a `tests/` directory or named `tests.rs` get a
//! pass — they're test-only by convention. The module's own files
//! (the runtime impl, this guard, the re-export site) are
//! whitelisted explicitly.

use std::fs;
use std::path::{Path, PathBuf};

const REFERENCE: &str = "pass::testing::PassRuntime";

/// Files allowed to name the path unconditionally: the module's own
/// files plus the public re-export site.
fn is_module_internal(path: &Path) -> bool {
    let rel = path.strip_prefix(env!("CARGO_MANIFEST_DIR")).unwrap_or(path);
    let rel_str = rel.to_string_lossy().replace('\\', "/");
    matches!(
        rel_str.as_str(),
        "src/pass.rs"
            | "src/pass/testing/mod.rs"
            | "src/pass/testing/runtime.rs"
            | "src/pass/testing/lockdown.rs"
    )
}

fn is_under_tests_dir(path: &Path) -> bool {
    path.components()
        .any(|c| c.as_os_str().to_string_lossy() == "tests")
}

fn is_dedicated_test_file(path: &Path) -> bool {
    if is_under_tests_dir(path) {
        return true;
    }
    path.file_name()
        .and_then(|s| s.to_str())
        .map(|n| n == "tests.rs")
        .unwrap_or(false)
}

/// Per-line check: walk the file forward tracking brace depth. Each
/// time a `#[cfg(test)]` attribute precedes a `mod ... {` opener,
/// remember the depth that just opened — that scope and everything
/// inside it is "test scope". When the depth drops back below that
/// mark, the scope closes. At the target line, return whether any
/// active test scope is on the stack.
///
/// This is heuristic — it doesn't handle braces in string literals
/// or character literals or comments — but is correct for the pass
/// files in this repo, which use `mod tests { ... }` blocks without
/// such pathologies. Good enough for a guard test.
fn line_is_in_cfg_test_mod(body: &str, target_line_idx: usize) -> bool {
    let lines: Vec<&str> = body.lines().collect();
    let mut depth: i32 = 0;
    let mut test_scope_depths: Vec<i32> = Vec::new();
    let mut pending_cfg_test = false;

    for (line_no, line) in lines.iter().enumerate() {
        if line_no == target_line_idx {
            return !test_scope_depths.is_empty();
        }

        let trimmed = line.trim();
        if trimmed.starts_with("//") {
            continue;
        }

        if trimmed.starts_with("#[cfg(test)]") {
            pending_cfg_test = true;
            continue;
        }

        let is_mod_decl = trimmed.starts_with("mod ")
            || trimmed.starts_with("pub mod ")
            || trimmed.starts_with("pub(crate) mod ")
            || trimmed.starts_with("pub(super) mod ");
        let mut opened_test_scope_this_line = pending_cfg_test && is_mod_decl;

        for ch in line.chars() {
            match ch {
                '{' => {
                    depth += 1;
                    if opened_test_scope_this_line {
                        test_scope_depths.push(depth);
                        opened_test_scope_this_line = false;
                    }
                }
                '}' => {
                    if let Some(&d) = test_scope_depths.last() {
                        if d == depth {
                            test_scope_depths.pop();
                        }
                    }
                    depth -= 1;
                }
                _ => {}
            }
        }

        // Clear pending cfg(test) on any non-attribute, non-blank,
        // non-comment line — attributes can stack but only attach to
        // the next item.
        if !trimmed.starts_with("#[") {
            pending_cfg_test = false;
        }
    }

    false
}

fn collect_rs_files(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).expect("read_dir under src/") {
        let entry = entry.expect("dir entry");
        let path = entry.path();
        if path.is_dir() {
            collect_rs_files(&path, out);
        } else if path.extension().and_then(|s| s.to_str()) == Some("rs") {
            out.push(path);
        }
    }
}

#[test]
fn no_production_source_references_pass_testing_passruntime() {
    let src_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let mut files = Vec::new();
    collect_rs_files(&src_dir, &mut files);

    let mut violations = Vec::new();
    for path in files {
        if is_module_internal(&path) {
            continue;
        }
        let body = match fs::read_to_string(&path) {
            Ok(s) => s,
            Err(_) => continue,
        };
        if !body.contains(REFERENCE) {
            continue;
        }
        // Dedicated test files (under tests/ or named tests.rs) are
        // gated by their parent's `#[cfg(test)] mod ...;` declaration.
        if is_dedicated_test_file(&path) {
            continue;
        }
        for (line_no, line) in body.lines().enumerate() {
            if !line.contains(REFERENCE) {
                continue;
            }
            let trimmed = line.trim_start();
            // Allow doc-comment mentions ("see PassRuntime").
            if trimmed.starts_with("//") {
                continue;
            }
            if line_is_in_cfg_test_mod(&body, line_no) {
                continue;
            }
            violations.push(format!(
                "{}:{} → {}",
                path.display(),
                line_no + 1,
                line.trim()
            ));
        }
    }

    assert!(
        violations.is_empty(),
        "Production source files referencing `pass::testing::PassRuntime` outside a \
         `#[cfg(test)]` scope. Use `CompiledSimulator` in production code. Offenders:\n{}",
        violations.join("\n")
    );
}

#[cfg(test)]
mod self_tests {
    use super::line_is_in_cfg_test_mod;

    #[test]
    fn detects_reference_inside_cfg_test_mod() {
        let body = "\
fn prod() {}

#[cfg(test)]
mod tests {
    fn uses_passruntime() {
        let _ = pass::testing::PassRuntime;
    }
}
";
        // Line 5 (0-indexed) = "        let _ = pass::testing::PassRuntime;"
        let lines: Vec<&str> = body.lines().collect();
        let target = lines
            .iter()
            .position(|l| l.contains("pass::testing::PassRuntime"))
            .unwrap();
        assert!(line_is_in_cfg_test_mod(body, target));
    }

    #[test]
    fn detects_reference_outside_cfg_test_mod() {
        let body = "\
fn prod() {
    let _ = pass::testing::PassRuntime;
}

#[cfg(test)]
mod tests {
    fn other() {}
}
";
        let lines: Vec<&str> = body.lines().collect();
        let target = lines
            .iter()
            .position(|l| l.contains("pass::testing::PassRuntime"))
            .unwrap();
        assert!(!line_is_in_cfg_test_mod(body, target));
    }

    #[test]
    fn handles_plain_mod_without_cfg_test_attribute() {
        let body = "\
mod tests {
    let _ = pass::testing::PassRuntime;
}
";
        let lines: Vec<&str> = body.lines().collect();
        let target = lines
            .iter()
            .position(|l| l.contains("pass::testing::PassRuntime"))
            .unwrap();
        assert!(!line_is_in_cfg_test_mod(body, target));
    }
}
