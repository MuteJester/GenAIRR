//! `MutationTransaction` — scoped mutation API for the 8 mutation
//! passes (uniform / S5F / PCR / quality / ncorrupt / contaminant /
//! indel / end_loss).
//!
//! Pre-refactor every mutation pass duplicated the same ~30-line
//! protocol: clone the simulation, open a `SimulationBuilder`,
//! `attach_standard_mutation_observers`, loop over base edits routing
//! each through `builder.change_base`, branch on
//! `ctx.reference_index.is_some()` to pick between
//! `seal_with_committed_live_calls` vs `seal`, then conditionally
//! apply `with_mutation_count`. Five implicit invariants per pass.
//!
//! `MutationTransaction` owns that protocol once. Passes call
//! `MutationTransaction::open(sim, ctx, ...)`, submit typed mutation
//! commands (`substitute_base`, `insert_base`, `delete_base`), and
//! `tx.commit()` returns the sealed `Simulation`. Out-of-range
//! handles return `Result<_, PassError>` instead of panicking
//! through `builder.change_base()`'s `.expect()`. Contract filtering
//! for substitution flows through the TX so contract violations
//! become structured `PassError`s in strict mode without each pass
//! having to call `sample_targeted_base` manually.
//!
//! ## Module layout
//!
//! - [`mod.rs`](self): struct definition, lifecycle (`open` /
//!   `commit`), focused accessors (`peek` / `rng` / `trace` /
//!   `ctx` / `split_borrows`), mutation-count sidecar
//!   (`add_to_mutation_count`), module docs, and the
//!   [`boundary_lockdown`] test.
//! - [`substitution`]: per-site base substitution API —
//!   `substitute_base`, `substitute_position_constrained`,
//!   `force_substitute_base`, `substitute_base_admitting`, plus
//!   the private `sample_with_filter` helper and the
//!   `SUBSTITUTE_BASE_EMPTY_SUPPORT` policy const.
//! - [`indel`]: structural-edit API — `insert_base`,
//!   `insert_base_admitting`, `delete_base`,
//!   `delete_base_admitting`.

mod indel;
mod substitution;

use crate::ir::{Simulation, SimulationBuilder};
use crate::pass::{PassContext, PassError};

/// Scoped mutation session for one pass execution.
///
/// Open with [`Self::open`], submit mutations via
/// [`Self::substitute_base`] / [`Self::insert_base`] /
/// [`Self::delete_base`], optionally call [`Self::add_to_mutation_count`],
/// and finalize with [`Self::commit`].
///
/// The TX borrows `ctx` mutably for its lifetime so the pass can't
/// touch trace / RNG directly while a TX is open. Use [`Self::rng`] /
/// [`Self::trace`] / [`Self::ctx`] to reach into the borrowed context
/// when the pass needs to sample a site or record a trace value
/// before calling a mutation method.
// Some helper variants are intentionally available before every pass
// needs them. The lockdown test at the bottom of this file enforces
// that production pass code mutates through this transaction layer.
#[allow(dead_code)]
pub(crate) struct MutationTransaction<'a, 'idx> {
    pub(super) builder: SimulationBuilder<'idx>,
    pub(super) ctx: &'a mut PassContext<'idx>,
    pub(super) pass_name: String,
    /// Counter for mutations that count toward AIRR `n_mutations`.
    /// Bumped by [`Self::add_to_mutation_count`]; written to
    /// `sealed.mutation_count` at commit.
    pub(super) mutations_applied: u32,
    /// Controls whether contract-filter exhaustion returns
    /// `PassError` or the pass-specific permissive sentinel
    /// (`Ok(false)` / no-op). Out-of-range handles are always
    /// errors; accepting them would corrupt the pool.
    pub(super) strict: bool,
}

#[allow(dead_code)]
impl<'a, 'idx> MutationTransaction<'a, 'idx> {
    /// Open a transaction for a mutation pass. Clones `sim` into a
    /// `SimulationBuilder` and attaches the standard mutation
    /// observers (walker for live-call updates + dirty-signal for
    /// refresh narrowing) when a reference index is available on
    /// the context.
    pub(crate) fn open(
        sim: &Simulation,
        ctx: &'a mut PassContext<'idx>,
        pass_name: &str,
        strict: bool,
    ) -> Self {
        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        builder.attach_standard_mutation_observers(ctx.reference_index);
        // Attach an event-log observer to the in-flight builder
        // when the caller supplied `ctx.event_log_sink`. The TX
        // drains this on commit so every `BaseChanged` /
        // `IndelInserted` / `IndelDeleted` fired by per-site
        // substitutions and indels flows through to the
        // caller-supplied buffer (and, in the compiled execution
        // path, ends up on the pass's `EventRecord`).
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }
        Self {
            builder,
            ctx,
            pass_name: pass_name.to_string(),
            mutations_applied: 0,
            strict,
        }
    }

    /// Read-only snapshot of the in-progress simulation. Sees every
    /// mutation committed via the TX so far.
    #[inline]
    pub(crate) fn peek(&self) -> &Simulation {
        self.builder.peek()
    }

    /// RNG borrow for site sampling. The TX holds `ctx` mutably for
    /// its lifetime; this is the accessor for the inner RNG.
    #[inline]
    pub(crate) fn rng(&mut self) -> &mut crate::rng::Rng {
        self.ctx.rng
    }

    /// Trace-injected replay cursor borrow. `Some(_)` when the run is
    /// in consume-trace mode; `None` for fresh-RNG runs. Passes that
    /// own their site sampling (e.g. N-corrupt, contaminant) branch
    /// on this to consume the recorded site/base from the cursor
    /// before calling the relevant substitution helper. Substitution
    /// helpers that own both site AND base sampling (notably
    /// [`Self::substitute_position_constrained`]) consult the cursor
    /// internally, so passes that call them don't need to touch this.
    #[inline]
    pub(crate) fn replay_cursor(&mut self) -> Option<&mut crate::replay::TraceCursor> {
        self.ctx.replay_cursor.as_deref_mut()
    }

    /// Trace borrow for recording site choices etc. before submitting
    /// the matching mutation command.
    #[inline]
    pub(crate) fn trace(&mut self) -> &mut crate::trace::Trace {
        self.ctx.trace
    }

    /// Full PassContext borrow. Use sparingly — prefer the focused
    /// accessors above. Provided for passes that need refdata or
    /// contracts access (e.g. for compile-fact registration during
    /// execution).
    #[inline]
    #[allow(dead_code)]
    pub(crate) fn ctx(&mut self) -> &mut PassContext<'idx> {
        self.ctx
    }

    /// Split-borrow accessor for passes whose iteration interleaves
    /// reads from the in-progress simulation, reads from contracts
    /// / refdata, and writes to RNG / trace — typically S5F-style
    /// per-step sampling. Returns disjoint borrows so the borrow
    /// checker accepts both being live at once.
    #[inline]
    pub(crate) fn split_borrows<'b>(&'b mut self) -> (&'b Simulation, &'b mut PassContext<'idx>) {
        (self.builder.peek(), self.ctx)
    }

    /// Record that `n` mutations applied during this pass should
    /// count toward AIRR `n_mutations`. Used by SHM passes (S5F,
    /// uniform). PCR / quality / N-corrupt / contaminant do NOT
    /// call this — they're sequencing artifacts, not biological
    /// mutations.
    pub(crate) fn add_to_mutation_count(&mut self, n: u32) {
        self.mutations_applied = self.mutations_applied.saturating_add(n);
    }

    /// Commit the transaction. Seals the builder, applying the
    /// reference-index-aware path when `ctx.reference_index.is_some()`
    /// so observer-staged live calls absorb cleanly. Applies the
    /// mutation count sidecar update if any mutations were tagged
    /// via [`Self::add_to_mutation_count`]. Returns the finished
    /// `Simulation`.
    ///
    /// Returns `Result<Simulation, _>` rather than the multi-pass
    /// `Outcome` aggregate because Pass::execute returns Simulation
    /// — the runtime loop accumulates per-pass Outcomes externally
    /// via PassContext / EventRecord, and a TX-produced single-pass
    /// "Outcome" would just be a sealed sim wrapped in mostly-empty
    /// metadata. If a future refactor evolves Pass::execute to
    /// return a richer per-pass result, this method's return type
    /// can grow without changing call sites that use `?` propagation.
    pub(crate) fn commit(self) -> Result<Simulation, PassError> {
        let MutationTransaction {
            mut builder,
            ctx,
            mutations_applied,
            ..
        } = self;

        // Drain the in-flight event log (if attached at open time)
        // *before* sealing, since `seal*` consume the builder. The
        // captured events are the per-site substitution / indel
        // consequences fired during the TX lifetime.
        if ctx.event_log_sink.is_some() {
            let events = builder.seal_event_log_observer();
            if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
                sink.extend(events);
            }
        }

        let mut sealed = if let Some(ref_index) = ctx.reference_index {
            builder.seal_with_committed_live_calls(ref_index)
        } else {
            builder.seal()
        };

        if mutations_applied > 0 {
            // Route the mutation-count bump through a fresh
            // `SimulationBuilder` so the `MutationCountChanged`
            // event flows to any attached sink. Zero-delta commits
            // are skipped entirely (no event, no builder wrap) —
            // mutation count is a derived metadata sidecar and an
            // unchanged value is not a state-change consequence.
            let new_count = sealed
                .mutation_count
                .saturating_add(mutations_applied);
            let mut count_builder = SimulationBuilder::from_simulation(sealed);
            if ctx.event_log_sink.is_some() {
                count_builder.attach_event_log_observer();
            }
            count_builder.set_mutation_count(new_count);
            if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
                sink.extend(count_builder.seal_event_log_observer());
            }
            sealed = count_builder.seal();
        }

        Ok(sealed)
    }
}

#[cfg(test)]
mod boundary_lockdown {
    //! Guard test that enforces the pass-author convention encoded
    //! in [`crate::ir::builder::SimulationBuilder`]: no file under
    //! `src/passes/` other than the `mutation_transaction/` module
    //! tree may call `.change_base(`, `.insert_indel(`, or
    //! `.delete_indel(` directly. New passes must route mutations
    //! through [`MutationTransaction`](super::MutationTransaction).
    //!
    //! Visibility lockdown via `pub(super)` isn't an option because
    //! `MutationTransaction` itself lives in `src/passes/`, not in
    //! `src/ir/`. So we enforce the rule with a grep-style assertion
    //! that runs in CI.
    //!
    //! Walks `src/passes/` from the manifest dir, skipping the
    //! whole `mutation_transaction/` subtree (the only legitimate
    //! callers), and asserts the forbidden call patterns are
    //! absent. Doc-comment references are excluded by matching
    //! only lines that contain `.` immediately followed by the
    //! method name and `(`.
    use std::fs;
    use std::path::Path;

    const FORBIDDEN: &[&str] = &[".change_base(", ".insert_indel(", ".delete_indel("];

    /// Returns `true` if `path` is anywhere under
    /// `src/passes/mutation_transaction/` (or is the legacy
    /// `mutation_transaction.rs` flat file — guard for compatibility
    /// during in-progress splits).
    fn is_inside_mutation_transaction(path: &Path) -> bool {
        let mut components = path.components();
        let mut prev: Option<std::path::Component<'_>> = None;
        while let Some(component) = components.next() {
            if let std::path::Component::Normal(name) = component {
                let s = name.to_string_lossy();
                if s == "mutation_transaction" || s == "mutation_transaction.rs" {
                    // Confirm the previous component (if any) was
                    // a directory named `passes` — defends against
                    // a hypothetical homonym elsewhere.
                    if let Some(std::path::Component::Normal(prev_name)) = prev {
                        if prev_name.to_string_lossy() == "passes" {
                            return true;
                        }
                    }
                }
            }
            prev = Some(component);
        }
        false
    }

    fn visit(dir: &Path, violations: &mut Vec<String>) {
        for entry in fs::read_dir(dir).expect("read_dir under src/passes") {
            let entry = entry.expect("dir entry");
            let path = entry.path();
            if is_inside_mutation_transaction(&path) {
                continue;
            }
            if path.is_dir() {
                visit(&path, violations);
                continue;
            }
            if path.extension().and_then(|s| s.to_str()) != Some("rs") {
                continue;
            }
            let body = match fs::read_to_string(&path) {
                Ok(s) => s,
                Err(_) => continue,
            };
            for (line_no, line) in body.lines().enumerate() {
                let trimmed = line.trim_start();
                // Skip doc-comment references like `builder.change_base()`.
                if trimmed.starts_with("//") || trimmed.starts_with("///") {
                    continue;
                }
                for pattern in FORBIDDEN {
                    if line.contains(pattern) {
                        violations.push(format!(
                            "{}:{} → {}",
                            path.display(),
                            line_no + 1,
                            line.trim()
                        ));
                    }
                }
            }
        }
    }

    #[test]
    fn no_pass_calls_low_level_builder_mutators_directly() {
        let passes_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/passes");
        let mut violations = Vec::new();
        visit(&passes_dir, &mut violations);
        assert!(
            violations.is_empty(),
            "Pass files calling SimulationBuilder mutation methods directly — \
             use MutationTransaction instead. Offenders:\n{}",
            violations.join("\n")
        );
    }

    // ──────────────────────────────────────────────────────────────
    // Production-path lockdown for persistent `Simulation::with_*`
    // mutators that now have event-emitting `SimulationBuilder`
    // equivalents.
    //
    // The previous slice landed builder methods that emit
    // `SimulationEvent`s (`assign_allele`, `update_trim`,
    // `add_region`, `replace_region`, `set_mutation_count`). The
    // runtime live-call refresh trusts those events — not
    // `PassCompileEffect` declarations — so a production pass that
    // calls `sim.with_allele_assigned(...)` directly silently
    // skips the refresh.
    //
    // The full-stack policy-conformance test catches this at run
    // time. This lockdown catches it at *code shape* before a pass
    // even reaches the audit plan.
    //
    // Unlike the low-level lockdown above, the forbidden methods
    // here ARE legitimately used in test scaffolding (every
    // `#[cfg(test)] mod tests` block builds fixtures via
    // `Simulation::with_*`). So the scanner skips lines inside
    // `#[cfg(test)]`-gated regions, identified via a simple
    // brace-depth tracker keyed off the attribute.
    // ──────────────────────────────────────────────────────────────

    /// Persistent `Simulation::with_*` mutators that have event-
    /// emitting `SimulationBuilder` equivalents. Production code in
    /// `src/passes/` (outside `mutation_transaction/`) must route
    /// through the builder; the test scaffolds are still free to
    /// build initial fixtures via these.
    const FORBIDDEN_PERSISTENT_MUTATORS: &[&str] = &[
        ".with_allele_assigned(",
        ".with_trim(",
        ".with_region_added(",
        ".with_region_replaced_for_segment(",
        ".with_mutation_count(",
    ];

    /// Strip an inline `//` line comment off a source line. Returns
    /// the prefix before `//`. Doesn't handle `//` inside string
    /// literals — that's a rare-enough case in pass code that the
    /// false-positive risk is negligible, and the scanner's
    /// downstream pattern check would still need a matching
    /// `.with_*(` literal in the source line to surface a false
    /// violation.
    fn strip_inline_comment(line: &str) -> &str {
        match line.find("//") {
            Some(idx) => &line[..idx],
            None => line,
        }
    }

    /// Scan `body` (the file source) for [`FORBIDDEN_PERSISTENT_MUTATORS`]
    /// outside `#[cfg(test)]`-gated regions, appending one violation
    /// string per offender to `violations`.
    ///
    /// The cfg(test) tracker is a tiny state machine: it counts
    /// `{` / `}` per line (after stripping inline comments) and
    /// when it encounters a `#[cfg(test)]` attribute followed by a
    /// brace-opening item, it pushes the pre-open brace depth onto
    /// a stack. Lines inside any active test region skip
    /// enforcement. The region pops when the brace depth returns
    /// to the pushed value.
    pub(super) fn scan_for_persistent_mutator_violations(
        path: &Path,
        body: &str,
        violations: &mut Vec<String>,
    ) {
        let mut brace_depth: i32 = 0;
        let mut test_stack: Vec<i32> = Vec::new();
        let mut pending_test_attr = false;

        for (line_no, line) in body.lines().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }

            let is_comment = trimmed.starts_with("//");

            // cfg(test) attribute? Defers to the next non-blank
            // code line, which may or may not open a brace.
            if !is_comment && trimmed.starts_with("#[cfg(test)]") {
                pending_test_attr = true;
                continue;
            }

            // Snapshot whether this line is inside a test region
            // BEFORE updating state, so the brace counter on this
            // line itself doesn't mis-attribute.
            let line_in_test = !test_stack.is_empty() || pending_test_attr;

            // Brace-count off a comment-stripped copy.
            let scan = strip_inline_comment(line);
            let opens = scan.matches('{').count() as i32;
            let closes = scan.matches('}').count() as i32;

            // If a cfg(test) attribute was pending and this line
            // opens a brace, the attribute gates the new braced
            // region. Record the pre-open depth so we know when to
            // pop. If no brace opens, the attribute gated a non-
            // braced item (`use`, `mod foo;`) and is simply
            // consumed.
            if pending_test_attr && opens > 0 {
                test_stack.push(brace_depth);
            }
            pending_test_attr = false;

            brace_depth += opens;
            brace_depth -= closes;

            // Pop any test region(s) we just exited (brace depth
            // dropped back to the depth at which the region
            // opened).
            while test_stack.last() == Some(&brace_depth) {
                test_stack.pop();
            }

            if is_comment {
                continue;
            }
            if line_in_test {
                continue;
            }

            for pattern in FORBIDDEN_PERSISTENT_MUTATORS {
                if line.contains(pattern) {
                    violations.push(format!(
                        "{}:{} → {}",
                        path.display(),
                        line_no + 1,
                        line.trim()
                    ));
                }
            }
        }
    }

    /// Returns `true` if `path` is in a directory or file the codebase
    /// reserves for test code reached only via `#[cfg(test)] mod ...;`
    /// declarations from production files. The persistent-mutator
    /// lockdown is a *production-path* guard, so these paths are
    /// allowed to call `Simulation::with_*` directly for fixture
    /// setup.
    ///
    /// Skip rules:
    /// - any path component named `tests` (covers `passes/foo/tests/`,
    ///   `passes/foo/tests/sub/`, etc.)
    /// - filenames `tests.rs` and `test_support.rs` (the canonical
    ///   convention for in-file test fixtures across this crate)
    fn is_test_only_path(path: &Path) -> bool {
        if path
            .components()
            .any(|c| matches!(c, std::path::Component::Normal(n) if n == "tests"))
        {
            return true;
        }
        if let Some(name) = path.file_name().and_then(|s| s.to_str()) {
            if name == "tests.rs" || name == "test_support.rs" {
                return true;
            }
        }
        false
    }

    fn visit_persistent(dir: &Path, violations: &mut Vec<String>) {
        for entry in fs::read_dir(dir).expect("read_dir under src/passes") {
            let entry = entry.expect("dir entry");
            let path = entry.path();
            if is_inside_mutation_transaction(&path) {
                continue;
            }
            if is_test_only_path(&path) {
                continue;
            }
            if path.is_dir() {
                visit_persistent(&path, violations);
                continue;
            }
            if path.extension().and_then(|s| s.to_str()) != Some("rs") {
                continue;
            }
            let body = match fs::read_to_string(&path) {
                Ok(s) => s,
                Err(_) => continue,
            };
            scan_for_persistent_mutator_violations(&path, &body, violations);
        }
    }

    #[test]
    fn no_pass_calls_persistent_with_star_mutators_in_production_code() {
        let passes_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/passes");
        let mut violations = Vec::new();
        visit_persistent(&passes_dir, &mut violations);
        assert!(
            violations.is_empty(),
            "Production pass code calling persistent `Simulation::with_*` mutators \
             directly — route through `SimulationBuilder` (e.g. \
             `builder.assign_allele(...)`, `builder.update_trim(...)`, \
             `builder.add_region(...)`, `builder.replace_region(...)`, \
             `builder.set_mutation_count(...)`) so the corresponding \
             `SimulationEvent` reaches the runtime derived-state refresh.\n\
             Offenders:\n{}",
            violations.join("\n")
        );
    }

    // ──────────────────────────────────────────────────────────────
    // Self-test for the scanner itself: it must catch a forbidden
    // call in production-position and not flag the same call inside
    // a `#[cfg(test)]` region. Run against synthetic source strings
    // so this test doesn't depend on the live repo state.
    // ──────────────────────────────────────────────────────────────

    #[test]
    fn scanner_flags_persistent_mutator_in_production_position() {
        let synthetic = r#"
fn execute(&self, sim: &Simulation) -> Simulation {
    sim.with_allele_assigned(Segment::V, instance)
}
"#;
        let mut violations = Vec::new();
        scan_for_persistent_mutator_violations(
            Path::new("synthetic.rs"),
            synthetic,
            &mut violations,
        );
        assert_eq!(violations.len(), 1);
        assert!(violations[0].contains("with_allele_assigned"));
    }

    #[test]
    fn scanner_allows_persistent_mutator_inside_cfg_test_block() {
        let synthetic = r#"
fn execute(&self, sim: &Simulation) -> Simulation {
    builder.assign_allele(Segment::V, instance);
    builder.seal()
}

#[cfg(test)]
mod tests {
    fn fixture() -> Simulation {
        let sim = Simulation::new();
        // Test scaffold — building initial state, not a runtime
        // mutation, so the lockdown should NOT fire here.
        sim.with_allele_assigned(Segment::V, instance)
            .with_trim(Segment::V, TrimEnd::Five, 2)
            .with_region_added(region)
            .with_mutation_count(3)
    }
}
"#;
        let mut violations = Vec::new();
        scan_for_persistent_mutator_violations(
            Path::new("synthetic.rs"),
            synthetic,
            &mut violations,
        );
        assert!(
            violations.is_empty(),
            "Scanner must not flag persistent mutators inside #[cfg(test)] regions. Got: {:?}",
            violations
        );
    }

    #[test]
    fn scanner_allows_cfg_test_gated_use_statements() {
        // `#[cfg(test)] use ...;` does not open a brace and must
        // not put the scanner into a sticky test-region state.
        // Subsequent production code in the same file MUST still
        // be checked.
        let synthetic = r#"
#[cfg(test)]
use crate::ir::NucHandle;

fn production(sim: &Simulation) -> Simulation {
    sim.with_region_added(region)
}
"#;
        let mut violations = Vec::new();
        scan_for_persistent_mutator_violations(
            Path::new("synthetic.rs"),
            synthetic,
            &mut violations,
        );
        assert_eq!(
            violations.len(),
            1,
            "Scanner must not stay in test-region state after a \
             cfg(test)-gated `use`. Got: {:?}",
            violations
        );
        assert!(violations[0].contains("with_region_added"));
    }
}

#[cfg(test)]
mod commit_event_tests {
    //! Pinning the commit-path event contract:
    //!
    //! - `mutations_applied > 0` → one `MutationCountChanged`
    //!   event per commit, with `old` taken from the sim's prior
    //!   `mutation_count` (so multi-pass accumulation is reflected
    //!   in the payload), `new == old + mutations_applied`, and
    //!   `delta == mutations_applied` cast to `i32`.
    //! - `mutations_applied == 0` → no event at all (the commit
    //!   path skips both the builder wrap and the
    //!   `with_mutation_count` call).
    //! - When `ctx.event_log_sink` is `None`, the commit path
    //!   still bumps the count but the broadcast is a no-op.

    use super::*;
    use crate::address::ChoiceAddress;
    use crate::ir::{NucHandle, Nucleotide, Segment, Simulation, SimulationEvent};
    use crate::rng::Rng;
    use crate::trace::Trace;

    fn sim_with_bases(bases: &[u8]) -> Simulation {
        let mut sim = Simulation::new();
        for (i, &b) in bases.iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
                b,
                i as u16,
                Segment::V,
            ));
            sim = next;
        }
        sim
    }

    /// Open a transaction against `sim` with `event_log_sink`
    /// attached, call the supplied `apply` closure to drive
    /// mutations, then commit and return `(sealed, captured_events)`.
    fn run_with_capture<F>(
        sim: Simulation,
        apply: F,
    ) -> (Simulation, Vec<SimulationEvent>)
    where
        F: for<'a, 'idx> FnOnce(&mut MutationTransaction<'a, 'idx>),
    {
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let mut captured: Vec<SimulationEvent> = Vec::new();
        let sealed;
        {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: None,
                contracts: None,
                feasibility: None,
                reference_index: None,
                replay_cursor: None,
                event_log_sink: Some(&mut captured),
            };
            let mut tx = MutationTransaction::open(&sim, &mut ctx, "test.commit", false);
            apply(&mut tx);
            sealed = tx.commit().expect("commit should succeed");
        }
        (sealed, captured)
    }

    #[test]
    fn commit_with_one_applied_substitution_emits_mutation_count_changed() {
        // Pre-state mutation_count = 0; one substitution tagged
        // via add_to_mutation_count(1). Commit should fire exactly
        // one `MutationCountChanged { old: 0, new: 1, delta: 1 }`.
        let sim = sim_with_bases(b"AAAA");
        let (sealed, events) = run_with_capture(sim, |tx| {
            // The substitution itself produces a `BaseChanged`
            // event (from the previous slice). We assert on the
            // mutation-count event specifically below; both
            // events should appear in the captured stream.
            tx.substitute_base_admitting(NucHandle::new(1), b'G', ChoiceAddress::MutateUniformBase(0)).unwrap();
            tx.add_to_mutation_count(1);
        });

        assert_eq!(sealed.mutation_count, 1);

        // Find the MutationCountChanged event among the captured
        // stream and assert its payload exactly.
        let mc = events
            .iter()
            .find(|e| matches!(e, SimulationEvent::MutationCountChanged { .. }))
            .expect("commit must emit MutationCountChanged when mutations_applied > 0");
        assert_eq!(
            *mc,
            SimulationEvent::MutationCountChanged {
                old: 0,
                new: 1,
                delta: 1,
            }
        );
        // Exactly one such event — the commit path is the only
        // emitter of this variant today.
        let count_events: usize = events
            .iter()
            .filter(|e| matches!(e, SimulationEvent::MutationCountChanged { .. }))
            .count();
        assert_eq!(count_events, 1);
    }

    #[test]
    fn commit_with_zero_applied_mutations_emits_no_mutation_count_event() {
        // No `add_to_mutation_count` call → commit skips the
        // entire mutation-count code path. No event should appear.
        let sim = sim_with_bases(b"AAAA");
        let (sealed, events) = run_with_capture(sim, |_tx| {
            // Intentionally no mutations and no add_to_mutation_count.
        });

        assert_eq!(sealed.mutation_count, 0);
        assert!(
            !events
                .iter()
                .any(|e| matches!(e, SimulationEvent::MutationCountChanged { .. })),
            "zero-delta commit must not emit MutationCountChanged; got events: {events:?}"
        );
    }

    #[test]
    fn commit_accumulates_mutation_count_across_passes() {
        // Sim already has mutation_count = 2 (from a hypothetical
        // earlier pass). This commit adds 3 more. The event must
        // carry the absolute old/new values plus the +3 delta,
        // not start from zero.
        let sim = sim_with_bases(b"AAAA").with_mutation_count(2);
        let (sealed, events) = run_with_capture(sim, |tx| {
            tx.substitute_base_admitting(NucHandle::new(0), b'C', ChoiceAddress::MutateUniformBase(0)).unwrap();
            tx.substitute_base_admitting(NucHandle::new(1), b'C', ChoiceAddress::MutateUniformBase(1)).unwrap();
            tx.substitute_base_admitting(NucHandle::new(2), b'C', ChoiceAddress::MutateUniformBase(2)).unwrap();
            tx.add_to_mutation_count(3);
        });

        assert_eq!(sealed.mutation_count, 5);

        let mc = events
            .iter()
            .find(|e| matches!(e, SimulationEvent::MutationCountChanged { .. }))
            .expect("commit must emit MutationCountChanged");
        assert_eq!(
            *mc,
            SimulationEvent::MutationCountChanged {
                old: 2,
                new: 5,
                delta: 3,
            }
        );
    }

    #[test]
    fn commit_without_event_log_sink_still_bumps_count() {
        // Production-path parity: when no sink is attached, the
        // count must still update (broadcast is a no-op). This is
        // the path every mutation pass runs in production today.
        let sim = sim_with_bases(b"AAAA");
        let mut trace = Trace::new();
        let mut rng = Rng::new(0);
        let sealed;
        {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata: None,
                contracts: None,
                feasibility: None,
                reference_index: None,
                replay_cursor: None,
                event_log_sink: None,
            };
            let mut tx = MutationTransaction::open(&sim, &mut ctx, "test.commit", false);
            tx.substitute_base_admitting(NucHandle::new(0), b'T', ChoiceAddress::MutateUniformBase(0)).unwrap();
            tx.add_to_mutation_count(1);
            sealed = tx.commit().expect("commit succeeds");
        }
        assert_eq!(sealed.mutation_count, 1);
    }
}
