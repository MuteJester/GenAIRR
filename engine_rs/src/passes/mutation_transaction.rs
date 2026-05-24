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

use crate::contract::ChoiceContext;
use crate::dist::{sample_filtered_result, Distribution};
use crate::ir::Simulation;
use crate::ir::{NucHandle, Nucleotide, SimulationBuilder};
use crate::pass::{PassContext, PassError};
use crate::trace::ChoiceValue;

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
// Phase 1: scaffolding only; no pass uses this yet. Phase 2
// migrations will remove the allow(dead_code).
#[allow(dead_code)]
pub(crate) struct MutationTransaction<'a, 'idx> {
    builder: SimulationBuilder<'idx>,
    ctx: &'a mut PassContext<'idx>,
    pass_name: String,
    /// Counter for mutations that count toward AIRR `n_mutations`.
    /// Bumped by [`Self::add_to_mutation_count`]; written to
    /// `sealed.mutation_count` at commit.
    mutations_applied: u32,
    /// Controls whether contract-filter exhaustion / out-of-range
    /// handles return `PassError` or fall back to permissive
    /// behavior (silent fallback to unconstrained sample for
    /// substitution; the apply still happens for invalid handles
    /// — no, that would be unsafe. Permissive mode still rejects
    /// invalid handles since allowing them would corrupt the pool).
    strict: bool,
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
    pub(crate) fn split_borrows<'b>(
        &'b mut self,
    ) -> (&'b Simulation, &'b mut PassContext<'idx>) {
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

    /// Substitute the base at `site` with a value drawn from
    /// `base_dist`, contract-filtered if a `ContractSet` is active
    /// on the context.
    ///
    /// Behavior:
    /// - Contract-filtered sampling: if `ctx.contracts.is_some()`,
    ///   the candidate base is drawn from `base_dist` with
    ///   `ChoiceContext::targeted_base_substitution(draw_index,
    ///   draw_count, site)` and any candidate the active contracts
    ///   reject is filtered out.
    /// - Strict mode + filter exhausted: returns `PassError::ConstraintSampling`.
    /// - Permissive mode + filter exhausted: falls back to an
    ///   unconstrained `base_dist.sample(rng)` draw (matching the
    ///   pre-TX `sample_targeted_base` behavior).
    /// - Optional `value_transform` is applied to the drawn base
    ///   before commitment (e.g. `lowercase_base` for quality
    ///   errors).
    /// - Records the chosen base to the trace at `address`.
    /// - Applies the substitution via the inner builder (which fires
    ///   `on_base_changed` to observers).
    ///
    /// Returns the chosen base on success. Returns `Err` for
    /// strict-mode contract failures or for an out-of-range
    /// `site` (which the old `builder.change_base().expect()` would
    /// have panicked on).
    pub(crate) fn substitute_base(
        &mut self,
        site: NucHandle,
        base_dist: &dyn Distribution<Output = u8>,
        address: &str,
        draw_index: u32,
        draw_count: u32,
        value_transform: Option<fn(u8) -> u8>,
    ) -> Result<u8, PassError> {
        // Validate site is in pool range BEFORE sampling (so we
        // don't waste an RNG draw on an invalid handle).
        let pool_len = self.builder.peek().pool.len() as u32;
        if site.index() >= pool_len {
            return Err(PassError::invalid_plan_state(
                self.pass_name.to_string(),
                format!(
                    "substitute_base: site handle {} out of pool range [0, {})",
                    site.index(),
                    pool_len
                ),
            ));
        }

        let transform = value_transform.unwrap_or(identity_base);
        let base = self.sample_with_filter(base_dist, address, draw_index, draw_count, site, transform)?;
        let base = transform(base);

        // Record the chosen base to the trace.
        self.ctx
            .trace
            .record(address, ChoiceValue::Base(base));

        // Apply through the builder. Safe — handle range validated
        // above; `change_base` will succeed.
        self.builder.change_base(site, base);

        Ok(base)
    }

    /// Substitute the base at `site` with a known value, without
    /// running it through any distribution or contract filter.
    /// Used by passes whose mutation is unconditional (e.g.
    /// N-corruption always writes `b'N'`). Does NOT record to the
    /// trace — caller decides what to record. Validates the handle
    /// range.
    pub(crate) fn substitute_base_fixed(
        &mut self,
        site: NucHandle,
        new_base: u8,
    ) -> Result<(), PassError> {
        let pool_len = self.builder.peek().pool.len() as u32;
        if site.index() >= pool_len {
            return Err(PassError::invalid_plan_state(
                self.pass_name.clone(),
                format!(
                    "substitute_base_fixed: site handle {} out of pool range [0, {})",
                    site.index(),
                    pool_len
                ),
            ));
        }
        self.builder.change_base(site, new_base);
        Ok(())
    }

    /// Insert a nucleotide at pool position `at`. Position must be
    /// in `[0, pool_len]` (insertion at the end is allowed). Returns
    /// `Err` for out-of-range positions.
    pub(crate) fn insert_base(
        &mut self,
        at: u32,
        nucleotide: Nucleotide,
    ) -> Result<(), PassError> {
        let pool_len = self.builder.peek().pool.len() as u32;
        if at > pool_len {
            return Err(PassError::invalid_plan_state(
                self.pass_name.to_string(),
                format!(
                    "insert_base: position {} out of pool range [0, {}]",
                    at, pool_len
                ),
            ));
        }
        self.builder.insert_indel(at, nucleotide);
        Ok(())
    }

    /// Delete the nucleotide at pool position `at`. Returns the
    /// removed nucleotide, or `Err` if `at >= pool_len` or the pool
    /// is empty.
    pub(crate) fn delete_base(&mut self, at: u32) -> Result<Nucleotide, PassError> {
        let pool_len = self.builder.peek().pool.len() as u32;
        if at >= pool_len {
            return Err(PassError::invalid_plan_state(
                self.pass_name.to_string(),
                format!(
                    "delete_base: position {} out of pool range [0, {})",
                    at, pool_len
                ),
            ));
        }
        self.builder.delete_indel(at).ok_or_else(|| {
            PassError::invalid_plan_state(
                self.pass_name.to_string(),
                format!("delete_base: pool returned None for position {}", at),
            )
        })
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
            builder,
            ctx,
            mutations_applied,
            ..
        } = self;

        let mut sealed = if let Some(ref_index) = ctx.reference_index {
            builder.seal_with_committed_live_calls(ref_index)
        } else {
            builder.seal()
        };

        if mutations_applied > 0 {
            sealed = sealed.with_mutation_count(
                sealed.mutation_count.saturating_add(mutations_applied),
            );
        }

        Ok(sealed)
    }

    /// Sample a candidate base from `base_dist`, contract-filtered
    /// against the active contract set when present. Returns the
    /// pre-transform base; caller applies the transform.
    fn sample_with_filter(
        &mut self,
        base_dist: &dyn Distribution<Output = u8>,
        addr: &str,
        draw_index: u32,
        draw_count: u32,
        site: NucHandle,
        transform: fn(u8) -> u8,
    ) -> Result<u8, PassError> {
        let contracts = self.ctx.contracts;
        let refdata = self.ctx.refdata;
        if let Some(contracts) = contracts {
            let context =
                ChoiceContext::targeted_base_substitution(draw_index, draw_count, site);
            match sample_filtered_result(self.ctx.rng, base_dist, |candidate: &u8| {
                let transformed = transform(*candidate);
                contracts
                    .admits_with_context(
                        self.builder.peek(),
                        refdata,
                        addr,
                        &ChoiceValue::Base(transformed),
                        context,
                    )
                    .is_ok()
            }) {
                Ok(value) => return Ok(value),
                Err(reason) if self.strict => {
                    return Err(PassError::constraint_sampling(
                        self.pass_name.to_string(),
                        addr.to_string(),
                        reason,
                    ));
                }
                Err(_) => {
                    // Permissive: fall through to unconstrained sample.
                }
            }
        }
        Ok(base_dist.sample(self.ctx.rng))
    }
}

#[inline]
fn identity_base(base: u8) -> u8 {
    base
}

#[cfg(test)]
mod boundary_lockdown {
    //! Guard test that enforces the pass-author convention encoded
    //! in [`crate::ir::builder::SimulationBuilder`]: no file under
    //! `src/passes/` other than this one may call `.change_base(`,
    //! `.insert_indel(`, or `.delete_indel(` directly. New passes
    //! must route mutations through [`MutationTransaction`].
    //!
    //! Visibility lockdown via `pub(super)` isn't an option because
    //! `MutationTransaction` itself lives in `src/passes/`, not in
    //! `src/ir/`. So we enforce the rule with a grep-style assertion
    //! that runs in CI.
    //!
    //! Walks `src/passes/` from the manifest dir, ignoring this file
    //! (the sole legitimate caller), and asserts the forbidden call
    //! patterns are absent. Doc-comment references are excluded by
    //! matching only lines that contain `.` immediately followed by
    //! the method name and `(`.
    use std::fs;
    use std::path::Path;

    const FORBIDDEN: &[&str] = &[".change_base(", ".insert_indel(", ".delete_indel("];

    fn visit(dir: &Path, violations: &mut Vec<String>) {
        for entry in fs::read_dir(dir).expect("read_dir under src/passes") {
            let entry = entry.expect("dir entry");
            let path = entry.path();
            if path.is_dir() {
                visit(&path, violations);
                continue;
            }
            if path.extension().and_then(|s| s.to_str()) != Some("rs") {
                continue;
            }
            if path.file_name().and_then(|s| s.to_str())
                == Some("mutation_transaction.rs")
            {
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
}

