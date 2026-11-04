//! Shared machinery for guards that read the crate's own source.
//!
//! Some invariants cannot be stated in the type system or in
//! visibility. "No production code constructs a generator from a raw
//! seed" is one: `Rng::new` must stay `pub` so integration tests can
//! build a generator, which means `pub` cannot also mean "tests
//! only". A test that walks `src/` and fails on an unsanctioned call
//! site can state it.
//!
//! Three guards in this crate work that way — the `PassRuntime`
//! lockdown in [`crate::pass::testing`], the two builder-mutator
//! lockdowns in [`crate::passes::mutation_transaction`], and the
//! raw-seed allowlist in [`crate::rng::discipline`] — and each needs
//! the same three things: the list of source files, a way to tell a
//! dedicated test file from a production one, and a way to tell
//! whether a line sits inside a `#[cfg(test)]` scope. They each used
//! to carry their own copy, which meant a fix to one brace tracker
//! left the others wrong. There is one implementation here, and the
//! guards hold only their policy: what to look for, and who may name
//! it anyway.
//!
//! The layering is deliberate: [`mentions_in_body`] is pure and takes
//! source as a string, so a guard's scanner can be tested against
//! synthetic code rather than against whatever the repository
//! happens to contain today. [`mentions`] walks the tree and calls
//! it.
//!
//! Compiled under `cfg(test)` only; nothing here ships.

use std::fs;
use std::path::{Path, PathBuf};

/// Dedicated test files, by crate convention: anything under a
/// `tests/` directory, `tests.rs`, and `test_support.rs`. Each is
/// reachable only through a `#[cfg(test)] mod ...;` declaration in
/// its parent, so visibility already guarantees they are test-only.
const TEST_ONLY_FILE_NAMES: &[&str] = &["tests.rs", "test_support.rs"];

/// One line of source naming something a guard looks for.
#[derive(Debug, Clone)]
pub struct Mention {
    /// Path relative to the crate root, forward-slashed on every
    /// platform so an allowlist can be written once and read on any.
    pub file: String,
    /// 1-based, so it reads like an editor's gutter.
    pub line: usize,
    /// The line itself, trimmed. The original text, not the
    /// comment-stripped copy the match was made against.
    pub text: String,
    /// Whether this line is test-only — inside a `#[cfg(test)]`
    /// scope, gated directly by a `#[cfg(test)]` attribute, or in a
    /// dedicated test file. A guard that bans something outright
    /// ignores this; a production-path guard filters on it.
    pub in_test_scope: bool,
}

impl std::fmt::Display for Mention {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{} → {}", self.file, self.line, self.text)
    }
}

/// Every line in `src/` naming any of `needles`, each tagged with
/// whether it is test-only.
///
/// A guard filters the result against its own policy; what counts as
/// sanctioned is the guard's business, not this function's.
pub fn mentions(needles: &[&str]) -> Vec<Mention> {
    let mut found = Vec::new();
    for path in source_files() {
        let body = match fs::read_to_string(&path) {
            Ok(s) => s,
            Err(_) => continue,
        };
        if !needles.iter().any(|n| body.contains(n)) {
            continue;
        }
        let file = crate_relative(&path);
        let mut in_file = mentions_in_body(&file, &body, needles);
        if is_dedicated_test_file(&path) {
            // The file is reached only through a `#[cfg(test)] mod`
            // declaration, so nothing in it is production code
            // whatever the braces inside it say.
            for mention in &mut in_file {
                mention.in_test_scope = true;
            }
        }
        found.append(&mut in_file);
    }
    found
}

/// [`mentions`] for a single needle, keeping only the production
/// lines — the common case.
pub fn mentions_in_production(needle: &str) -> Vec<Mention> {
    mentions(&[needle])
        .into_iter()
        .filter(|m| !m.in_test_scope)
        .collect()
}

/// Every line of `body` naming any of `needles`, tagged with whether
/// it is inside a `#[cfg(test)]` scope.
///
/// Matching happens against the line with any trailing `//` comment
/// removed, so prose that names the thing is not read as a use of
/// it. `file` is used only to label the result.
///
/// Test scope is tracked by walking forward and counting braces.
/// Each time a `#[cfg(test)]` attribute precedes a brace-opening
/// item, the depth that just opened is remembered; when the depth
/// drops back below the mark, the scope closes. An attribute that
/// gates a non-braced item (`#[cfg(test)] use ...;`) covers that one
/// line and does not make the rest of the file test-only — the case
/// that is easiest to get wrong, and the reason this is in one place.
///
/// Heuristic: braces inside string literals, character literals and
/// block comments are counted. That is correct for this crate, which
/// writes `mod tests { ... }` without such pathologies, and a guard
/// confused by a future pathology fails loudly rather than quietly
/// passing a violation.
pub fn mentions_in_body(file: &str, body: &str, needles: &[&str]) -> Vec<Mention> {
    let mut found = Vec::new();
    let mut depth: i32 = 0;
    let mut test_scope_depths: Vec<i32> = Vec::new();
    let mut pending_cfg_test = false;

    for (index, line) in body.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.starts_with("//") {
            continue;
        }

        if trimmed.starts_with("#[cfg(test)]") {
            pending_cfg_test = true;
            continue;
        }

        // Snapshot before this line's own braces move the depth: a
        // line that opens a test scope is itself test scope, and a
        // line gated directly by a pending attribute is too.
        let in_test_scope = !test_scope_depths.is_empty() || pending_cfg_test;

        let is_mod_decl = trimmed.starts_with("mod ")
            || trimmed.starts_with("pub mod ")
            || trimmed.starts_with("pub(crate) mod ")
            || trimmed.starts_with("pub(super) mod ");
        let mut opens_test_scope = pending_cfg_test && is_mod_decl;

        let code = strip_inline_comment(line);
        for ch in code.chars() {
            match ch {
                '{' => {
                    depth += 1;
                    if opens_test_scope {
                        test_scope_depths.push(depth);
                        opens_test_scope = false;
                    }
                }
                '}' => {
                    if test_scope_depths.last() == Some(&depth) {
                        test_scope_depths.pop();
                    }
                    depth -= 1;
                }
                _ => {}
            }
        }

        // An attribute may stack with other attributes but attaches
        // only to the next item, so any other content clears it.
        if !trimmed.starts_with("#[") {
            pending_cfg_test = false;
        }

        if needles.iter().any(|n| code.contains(n)) {
            found.push(Mention {
                file: file.to_string(),
                line: index + 1,
                text: trimmed.to_string(),
                in_test_scope,
            });
        }
    }

    found
}

/// Every `.rs` file under the crate's `src/`, in no particular order.
pub fn source_files() -> Vec<PathBuf> {
    let mut out = Vec::new();
    collect_rs_files(&Path::new(env!("CARGO_MANIFEST_DIR")).join("src"), &mut out);
    out
}

/// `src/rng/substream.rs` — forward slashes on every platform, so a
/// guard's allowlist is a list of literals rather than a path dance.
pub fn crate_relative(path: &Path) -> String {
    path.strip_prefix(env!("CARGO_MANIFEST_DIR"))
        .unwrap_or(path)
        .to_string_lossy()
        .replace('\\', "/")
}

/// Whether `path` is a dedicated test file by crate convention. See
/// [`TEST_ONLY_FILE_NAMES`].
pub fn is_dedicated_test_file(path: &Path) -> bool {
    if path
        .components()
        .any(|c| c.as_os_str().to_string_lossy() == "tests")
    {
        return true;
    }
    path.file_name()
        .and_then(|s| s.to_str())
        .map(|n| TEST_ONLY_FILE_NAMES.contains(&n))
        .unwrap_or(false)
}

/// The part of `line` before any `//`. Does not handle `//` inside a
/// string literal; a guard's needle would have to appear in the same
/// line's prose for that to matter.
fn strip_inline_comment(line: &str) -> &str {
    match line.find("//") {
        Some(idx) => &line[..idx],
        None => line,
    }
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

#[cfg(test)]
mod tests {
    use super::{crate_relative, is_dedicated_test_file, mentions_in_body, source_files};
    use std::path::Path;

    const NEEDLE: &str = "pass::testing::PassRuntime";

    /// Whether the one line naming [`NEEDLE`] is test scope.
    fn needle_is_test_scope(body: &str) -> bool {
        let found = mentions_in_body("synthetic.rs", body, &[NEEDLE]);
        assert_eq!(found.len(), 1, "fixture should name the needle once");
        found[0].in_test_scope
    }

    #[test]
    fn a_line_inside_a_cfg_test_mod_is_test_scope() {
        assert!(needle_is_test_scope(
            "\
fn prod() {}

#[cfg(test)]
mod tests {
    fn uses_passruntime() {
        let _ = pass::testing::PassRuntime;
    }
}
"
        ));
    }

    #[test]
    fn a_line_before_the_cfg_test_mod_is_production_scope() {
        assert!(!needle_is_test_scope(
            "\
fn prod() {
    let _ = pass::testing::PassRuntime;
}

#[cfg(test)]
mod tests {
    fn other() {}
}
"
        ));
    }

    #[test]
    fn a_mod_named_tests_without_the_attribute_is_production_scope() {
        assert!(!needle_is_test_scope(
            "\
mod tests {
    let _ = pass::testing::PassRuntime;
}
"
        ));
    }

    #[test]
    fn a_scope_that_has_closed_again_is_not_still_open() {
        assert!(!needle_is_test_scope(
            "\
#[cfg(test)]
mod tests {
    fn inside() {}
}

fn after() {
    let _ = pass::testing::PassRuntime;
}
"
        ));
    }

    #[test]
    fn a_cfg_test_gated_use_does_not_make_the_rest_of_the_file_test_only() {
        // The sticky-state bug: `#[cfg(test)] use ...;` opens no
        // brace, so production code after it must still be checked.
        assert!(!needle_is_test_scope(
            "\
#[cfg(test)]
use crate::ir::NucHandle;

fn production() {
    let _ = pass::testing::PassRuntime;
}
"
        ));
    }

    #[test]
    fn a_line_gated_directly_by_the_attribute_is_test_scope() {
        assert!(needle_is_test_scope(
            "\
#[cfg(test)]
use pass::testing::PassRuntime;

fn production() {}
"
        ));
    }

    #[test]
    fn prose_naming_the_needle_is_not_a_use_of_it() {
        assert!(mentions_in_body(
            "synthetic.rs",
            "\
//! See pass::testing::PassRuntime for why.
/// Prefer CompiledSimulator over pass::testing::PassRuntime.
fn prod() {
    let _ = other(); // not pass::testing::PassRuntime
}
",
            &[NEEDLE]
        )
        .is_empty());
    }

    #[test]
    fn several_needles_are_found_in_one_walk_of_the_body() {
        let found = mentions_in_body(
            "synthetic.rs",
            "\
fn prod(sim: &Simulation) -> Simulation {
    sim.with_allele_assigned(Segment::V, instance)
        .with_trim(Segment::V, TrimEnd::Five, 2)
}
",
            &[".with_allele_assigned(", ".with_trim("],
        );
        assert_eq!(found.len(), 2);
        assert!(found.iter().all(|m| !m.in_test_scope));
        assert_eq!(found[0].line, 2);
        assert_eq!(found[1].line, 3);
    }

    #[test]
    fn the_walk_finds_this_crates_own_modules() {
        let names: Vec<String> = source_files().iter().map(|p| crate_relative(p)).collect();
        assert!(
            names.iter().any(|n| n == "src/rng.rs"),
            "the walk should reach src/rng.rs; found {} files",
            names.len()
        );
        assert!(
            names.iter().any(|n| n == "src/rng/substream.rs"),
            "the walk should recurse into subdirectories"
        );
    }

    #[test]
    fn dedicated_test_files_are_recognised_by_every_convention() {
        assert!(is_dedicated_test_file(Path::new(
            "src/dist/tests/uniform.rs"
        )));
        assert!(is_dedicated_test_file(Path::new(
            "src/passes/invert_d/tests.rs"
        )));
        assert!(is_dedicated_test_file(Path::new(
            "src/passes/test_support.rs"
        )));
        assert!(!is_dedicated_test_file(Path::new("src/passes/invert_d.rs")));
    }
}
