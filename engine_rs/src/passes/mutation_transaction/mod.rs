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
pub(crate) mod substitution;

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
    /// Substitutions whose drawn base equalled the base already present.
    /// Subtracted in `commit`: a replacement that changes nothing is not a
    /// mutation, and counting it inflated `n_mutations` by roughly 26 points
    /// under the uniform model.
    pub(super) no_op_substitutions: u32,
    /// Input molecule length before this pass can change observation coverage.
    pub(super) input_pool_len: usize,
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
        builder.set_duplex_edit_context(Some(pass_name));
        Self {
            builder,
            ctx,
            pass_name: pass_name.to_string(),
            mutations_applied: 0,
            no_op_substitutions: 0,
            input_pool_len: sim.pool.len(),
            strict,
        }
    }

    /// Preserve sampled attempts separately from applied substitution events.
    pub(crate) fn record_duplex_attempt_stage(
        &mut self,
        stage: crate::ir::DuplexAttemptStage,
        count: u32,
    ) {
        self.builder.record_duplex_attempt_stage(stage, count);
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
    /// count toward AIRR `n_mutations`. Substitutions that drew the
    /// base already present are subtracted in `commit`, so a pass
    /// may pass its full written count here. Used by SHM passes (S5F,
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
            no_op_substitutions,
            input_pool_len,
            ..
        } = self;

        // Realized-and-visible semantics: a substitution that drew the base
        // already present wrote a byte but changed nothing, so it is not a
        // mutation. The draw itself is untouched -- same RNG consumption,
        // same recorded trace, same sequence -- which is why this does not
        // change any simulated output, only the count of it.
        let mutations_applied = mutations_applied.saturating_sub(no_op_substitutions);

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
            let mut count_builder = SimulationBuilder::from_simulation(sealed);
            if ctx.event_log_sink.is_some() {
                count_builder.attach_event_log_observer();
            }
            count_builder.add_shm_mutations(mutations_applied, input_pool_len);
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
    //! Two guards over the pass-author convention encoded in
    //! [`crate::ir::builder::SimulationBuilder`]. Both are policy
    //! only: the source walk, the test-scope tracking and the
    //! test-file conventions live in [`crate::source_guard`], shared
    //! with the other guards in this crate.
    //!
    //! Visibility lockdown via `pub(super)` isn't an option for
    //! either, because `MutationTransaction` itself lives in
    //! `src/passes/`, not in `src/ir/`.
    //!
    //! **The low-level guard**: no file under `src/passes/` other
    //! than the `mutation_transaction/` tree may call
    //! `.change_base(`, `.insert_indel(` or `.delete_indel(` at all.
    //! New passes route mutations through
    //! [`MutationTransaction`](super::MutationTransaction). Test
    //! scaffolds are not exempt — there is no legitimate reason to
    //! reach these from a pass file.
    //!
    //! **The persistent-mutator guard**: no *production* code under
    //! `src/passes/` may call the `Simulation::with_*` mutators that
    //! now have event-emitting builder equivalents. The runtime
    //! live-call refresh trusts those events, not
    //! `PassCompileEffect` declarations, so a pass that calls
    //! `sim.with_allele_assigned(...)` directly silently skips the
    //! refresh. Unlike the low-level guard this one does exempt test
    //! scope: every `#[cfg(test)] mod tests` block builds fixtures
    //! through `Simulation::with_*`, which is what those methods are
    //! for. The full-stack policy-conformance test catches the same
    //! mistake at run time; this catches it at code shape, before a
    //! pass reaches the audit plan.

    use crate::source_guard::{mentions, Mention};

    /// Low-level builder mutators, forbidden outright in pass files.
    const FORBIDDEN: &[&str] = &[".change_base(", ".insert_indel(", ".delete_indel("];

    /// Persistent `Simulation::with_*` mutators that have event-
    /// emitting `SimulationBuilder` equivalents. Production code in
    /// `src/passes/` must route through the builder; test scaffolds
    /// are still free to build initial fixtures with these.
    const FORBIDDEN_PERSISTENT_MUTATORS: &[&str] = &[
        ".with_allele_assigned(",
        ".with_trim(",
        ".with_region_added(",
        ".with_region_replaced_for_segment(",
        ".with_mutation_count(",
    ];

    /// A pass file this module's guards apply to: under
    /// `src/passes/`, outside the `mutation_transaction/` tree that
    /// legitimately owns these calls.
    fn is_guarded_pass_file(mention: &Mention) -> bool {
        mention.file.starts_with("src/passes/")
            && !mention.file.starts_with("src/passes/mutation_transaction")
    }

    #[test]
    fn no_pass_calls_low_level_builder_mutators_directly() {
        let violations: Vec<String> = mentions(FORBIDDEN)
            .into_iter()
            .filter(is_guarded_pass_file)
            .map(|m| m.to_string())
            .collect();
        assert!(
            violations.is_empty(),
            "Pass files calling SimulationBuilder mutation methods directly — \
             use MutationTransaction instead. Offenders:\n{}",
            violations.join("\n")
        );
    }

    #[test]
    fn no_pass_calls_persistent_with_star_mutators_in_production_code() {
        let violations: Vec<String> = mentions(FORBIDDEN_PERSISTENT_MUTATORS)
            .into_iter()
            .filter(|m| !m.in_test_scope && is_guarded_pass_file(m))
            .map(|m| m.to_string())
            .collect();
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

    /// The guards above are policy over [`crate::source_guard`],
    /// whose scanner is tested against synthetic source in its own
    /// module. What is left to pin here is that each policy reads
    /// the `in_test_scope` flag the way its doc comment says: the
    /// low-level guard ignores it, the persistent-mutator guard
    /// honours it.
    #[test]
    fn the_two_policies_disagree_about_test_scope_by_design() {
        use crate::source_guard::mentions_in_body;

        let body = "\
#[cfg(test)]
mod tests {
    fn fixture() -> Simulation {
        builder.change_base(handle, base);
        sim.with_allele_assigned(Segment::V, instance)
    }
}
";
        let low_level = mentions_in_body("src/passes/synthetic.rs", body, FORBIDDEN);
        assert_eq!(low_level.len(), 1, "the fixture calls one forbidden method");
        assert!(
            low_level.iter().all(is_guarded_pass_file),
            "a pass file outside mutation_transaction/ is guarded"
        );

        let persistent = mentions_in_body(
            "src/passes/synthetic.rs",
            body,
            FORBIDDEN_PERSISTENT_MUTATORS,
        );
        assert_eq!(persistent.len(), 1);
        assert!(
            persistent.iter().all(|m| m.in_test_scope),
            "the same fixture is test scope, so the persistent guard passes it"
        );
    }

    #[test]
    fn the_mutation_transaction_tree_is_the_one_exemption() {
        let owner = Mention {
            file: "src/passes/mutation_transaction/substitution.rs".to_string(),
            line: 1,
            text: String::new(),
            in_test_scope: false,
        };
        let other_pass = Mention {
            file: "src/passes/mutate/uniform.rs".to_string(),
            line: 1,
            text: String::new(),
            in_test_scope: false,
        };
        let outside_passes = Mention {
            file: "src/ir/builder.rs".to_string(),
            line: 1,
            text: String::new(),
            in_test_scope: false,
        };
        assert!(!is_guarded_pass_file(&owner));
        assert!(is_guarded_pass_file(&other_pass));
        assert!(!is_guarded_pass_file(&outside_passes));
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
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        sim
    }

    /// Open a transaction against `sim` with `event_log_sink`
    /// attached, call the supplied `apply` closure to drive
    /// mutations, then commit and return `(sealed, captured_events)`.
    fn run_with_capture<F>(sim: Simulation, apply: F) -> (Simulation, Vec<SimulationEvent>)
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
            tx.substitute_base_admitting(
                NucHandle::new(1),
                b'G',
                ChoiceAddress::MutateUniformBase(0),
            )
            .unwrap();
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
            tx.substitute_base_admitting(
                NucHandle::new(0),
                b'C',
                ChoiceAddress::MutateUniformBase(0),
            )
            .unwrap();
            tx.substitute_base_admitting(
                NucHandle::new(1),
                b'C',
                ChoiceAddress::MutateUniformBase(1),
            )
            .unwrap();
            tx.substitute_base_admitting(
                NucHandle::new(2),
                b'C',
                ChoiceAddress::MutateUniformBase(2),
            )
            .unwrap();
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
            tx.substitute_base_admitting(
                NucHandle::new(0),
                b'T',
                ChoiceAddress::MutateUniformBase(0),
            )
            .unwrap();
            tx.add_to_mutation_count(1);
            sealed = tx.commit().expect("commit succeeds");
        }
        assert_eq!(sealed.mutation_count, 1);
    }
}
