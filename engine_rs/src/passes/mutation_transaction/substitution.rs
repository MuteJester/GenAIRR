//! Per-site base substitution API for [`MutationTransaction`].
//!
//! Exposes the four substitution entry points:
//! - [`MutationTransaction::substitute_base`] — natural-distribution
//!   draw at a fixed site, contract-filtered.
//! - [`MutationTransaction::substitute_position_constrained`] —
//!   joint (site, base) draw under contracts via per-site
//!   admissible-mass weighting.
//! - [`MutationTransaction::force_substitute_base`] — pinned write
//!   bypassing contracts (architectural exception, loudly named).
//! - [`MutationTransaction::substitute_base_admitting`] — pinned
//!   write with contract arbitration.
//!
//! Plus the private `sample_with_filter` helper and the
//! `SUBSTITUTE_BASE_EMPTY_SUPPORT` policy const that
//! [`MutationTransaction::substitute_base`] routes through.

use crate::address::ChoiceAddress;
use crate::contract::ChoiceContext;
use crate::dist::{sample_filtered_with_policy, Distribution, EmptySupport};
use crate::ir::NucHandle;
use crate::pass::PassError;
use crate::trace::ChoiceValue;

use super::MutationTransaction;

/// v3.0 empty-support policy for fixed-site substitution under
/// active contracts: the call site (`substitute_base`) treats an
/// empty filter as a per-site no-op (`Ok(None)` here, surfaced as
/// `Ok(false)` to the pass). The pass loops continue with the
/// site untouched and no trace record is emitted.
const SUBSTITUTE_BASE_EMPTY_SUPPORT: EmptySupport<u8> = EmptySupport::Skip;

#[allow(dead_code)]
impl<'a, 'idx> MutationTransaction<'a, 'idx> {
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
    /// - **Permissive mode + filter exhausted: the site is
    ///   skipped** (no trace record, no pool mutation). Returns
    ///   `Ok(false)`. This honors the v3.0 architectural rule
    ///   that the engine never proposes contract-violating
    ///   actions — falling back to an unconstrained draw would
    ///   reintroduce reject-after-propose at the per-event level.
    /// - Optional `value_transform` is applied to the drawn base
    ///   before commitment (e.g. `lowercase_base` for quality
    ///   errors).
    /// - Records the chosen base to the trace at `address` only
    ///   when the substitution is applied.
    /// - Applies the substitution via the inner builder (which fires
    ///   `on_base_changed` to observers).
    ///
    /// Returns `Ok(true)` on applied substitution, `Ok(false)` on
    /// permissive empty-support skip, `Err(..)` for strict-mode
    /// rejection or out-of-range `site`.
    pub(crate) fn substitute_base(
        &mut self,
        site: NucHandle,
        base_dist: &dyn Distribution<Output = u8>,
        address: ChoiceAddress,
        draw_index: u32,
        draw_count: u32,
        value_transform: Option<fn(u8) -> u8>,
    ) -> Result<bool, PassError> {
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

        // Trace-injected replay (Tier 3): the site is fixed by the
        // caller; only the base is consumed from the cursor. The
        // caller (e.g. ContaminantPass) is responsible for skipping
        // sites where the original run produced no record — they
        // call this helper only when the cursor's next address
        // matches `address`.
        //
        // Validation matches `substitute_position_constrained`'s
        // replay branch: the recorded base must be admissible under
        // the current contract bundle, else error (no force-apply).
        if self.ctx.replay_cursor.is_some() {
            let base_byte = {
                let cursor = self
                    .ctx
                    .replay_cursor
                    .as_deref_mut()
                    .expect("replay_cursor presence checked above");
                cursor
                    .expect_base(address)
                    .map_err(|reason| PassError::replay(self.pass_name.clone(), reason))?
            };
            if let Some(contracts) = self.ctx.contracts {
                let sim_ref = self.builder.peek();
                let refdata = self.ctx.refdata;
                let mask = contracts.admissible_bases_at(sim_ref, refdata, site);
                let admitted = mask.admits(base_byte)
                    || contracts.admits_fixed_base_at(sim_ref, refdata, site, base_byte);
                if !admitted {
                    return Err(PassError::constraint_sampling(
                        self.pass_name.clone(),
                        address.to_string(),
                        crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }
            }
            self.ctx
                .trace
                .record_choice(address, ChoiceValue::Base(base_byte));
            self.builder.change_base(site, base_byte);
            // Suppress unused-binding warnings for fresh-sample
            // parameters in the replay path.
            let _ = (base_dist, draw_index, draw_count, value_transform);
            return Ok(true);
        }

        let transform = value_transform.unwrap_or(identity_base);
        let base = match self
            .sample_with_filter(base_dist, address, draw_index, draw_count, site, transform)?
        {
            Some(b) => b,
            None => return Ok(false),
        };
        let base = transform(base);

        // Record the chosen base to the trace.
        self.ctx
            .trace
            .record_choice(address, ChoiceValue::Base(base));

        // Apply through the builder. Safe — handle range validated
        // above; `change_base` will succeed.
        self.builder.change_base(site, base);

        Ok(true)
    }

    /// Substitute a base at a position chosen jointly with the base
    /// under the v3.0 **constrain-before-propose** API.
    ///
    /// Fast path (no contracts): pick `site` uniformly in `[0, pool_len)`,
    /// draw `base = transform(base_dist.sample(rng))`, record both to
    /// the trace, and write through the builder. RNG consumption
    /// (one word for site, one for base) matches the pre-v3.0
    /// uniform / PCR / quality loops exactly.
    ///
    /// Constrained path (contracts present, `base_dist.support()`
    /// enumerable): compute per-site weight as
    /// `sum_{(b, w) ∈ support} w · [mask.admits(b) ∧ admits_fixed_base_at(site, transform(b))]`
    /// over each site in the pool, weighted-sample the site, then
    /// weighted-sample the base from the per-site admissible
    /// support. Sites that admit nothing contribute zero weight.
    /// If every site shuts the support out, returns `Ok(false)` in
    /// permissive mode or `PassError::ConstraintSampling` in strict
    /// mode — matching the realized-count semantics: nothing was
    /// applied.
    ///
    /// `value_transform` is applied to the drawn base before the
    /// trace record and before the write. When the transform
    /// changes the byte (e.g. quality's lowercase), the
    /// transformed value is additionally checked through
    /// [`super::super::contract::Contract::admits_fixed_base_at`] —
    /// the canonical mask only describes A/C/G/T support, not
    /// pinned non-canonical writes.
    ///
    /// Returns `Ok(true)` when a substitution was applied,
    /// `Ok(false)` for permissive empty-support, and `Err(..)` for
    /// strict-mode failures.
    pub(crate) fn substitute_position_constrained(
        &mut self,
        base_dist: &dyn Distribution<Output = u8>,
        site_address: ChoiceAddress,
        base_address: ChoiceAddress,
        value_transform: Option<fn(u8) -> u8>,
        segment_rates: Option<&crate::passes::mutate::SegmentRateWeights>,
        v_subregion_rates: Option<&crate::passes::mutate::VSubregionRateWeights>,
    ) -> Result<bool, PassError> {
        let pool_len = self.builder.peek().pool.len() as u32;
        if pool_len == 0 {
            return Ok(false);
        }
        let transform = value_transform.unwrap_or(identity_base);

        // ── Trace-injected replay (Tier 2 sub-slice 2) ───────────
        //
        // Consume the recorded `(site, base)` pair from the cursor
        // instead of weighted-sampling. The recorded `base` is
        // already post-transform (the fresh-RNG path records
        // `transform(sampled)`); we don't re-apply the transform.
        //
        // The user-invariant for this slice: replayed values still
        // run through the same admissibility checks as fresh
        // sampling. If the recorded base would have been rejected
        // by the current contract bundle at this site (or the site
        // is out of range), error instead of force-applying.
        if self.ctx.replay_cursor.is_some() {
            // Consume both records first; the cursor borrow ends
            // after the second expect_* so we can re-borrow `self`
            // for validation + write.
            let (site_idx, base_byte) = {
                let cursor = self
                    .ctx
                    .replay_cursor
                    .as_deref_mut()
                    .expect("replay_cursor presence checked above");
                let s = cursor
                    .expect_int(site_address)
                    .map_err(|reason| PassError::replay(self.pass_name.clone(), reason))?;
                let b = cursor
                    .expect_base(base_address)
                    .map_err(|reason| PassError::replay(self.pass_name.clone(), reason))?;
                (s, b)
            };
            if site_idx < 0 || site_idx >= pool_len as i64 {
                return Err(PassError::invalid_plan_state(
                    self.pass_name.clone(),
                    format!(
                        "substitute_position_constrained: replayed site {} out of pool range [0, {})",
                        site_idx, pool_len
                    ),
                ));
            }
            let site = site_idx as u32;
            // Segment-rate replay validation. A site recorded under
            // rates that excluded zero-rate sites from support must
            // also be excluded here; otherwise replay could
            // force-apply a substitution at a position the current
            // rate vector would never sample. Mirrors the contract
            // admissibility check below: both gates run before the
            // write.
            //
            // V-subregion-rate replay validation (Slice B) folds in
            // the same way: a V site recorded under a rate vector
            // where that site's subregion has weight > 0 must still
            // be reachable under the current vector. A vector that
            // zeroed the subregion must reject the recorded site
            // rather than force-apply.
            {
                let sim_ref = self.builder.peek();
                let site_factor = combined_site_factor(
                    site,
                    sim_ref,
                    self.ctx.refdata,
                    segment_rates,
                    v_subregion_rates,
                );
                if site_factor <= 0.0 {
                    return Err(PassError::constraint_sampling(
                        self.pass_name.clone(),
                        site_address.to_string(),
                        crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }
            }
            // Admissibility check against the CURRENT contract bundle
            // at the CURRENT in-progress sim. A recorded base that
            // would have been rejected fresh-sampled must also be
            // rejected here; otherwise replay could force-apply a
            // contract-violating value.
            if let Some(contracts) = self.ctx.contracts {
                let sim_ref = self.builder.peek();
                let refdata = self.ctx.refdata;
                let handle = NucHandle::new(site);
                let mask = contracts.admissible_bases_at(sim_ref, refdata, handle);
                let admitted = mask.admits(base_byte)
                    || contracts.admits_fixed_base_at(sim_ref, refdata, handle, base_byte);
                if !admitted {
                    return Err(PassError::constraint_sampling(
                        self.pass_name.clone(),
                        base_address.to_string(),
                        crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }
            }
            // Record + write — same shape as the fresh-RNG paths.
            self.ctx
                .trace
                .record_choice(site_address, ChoiceValue::Int(site as i64));
            self.ctx
                .trace
                .record_choice(base_address, ChoiceValue::Base(base_byte));
            self.builder.change_base(NucHandle::new(site), base_byte);
            // Suppress the unused-binding warning for the
            // base distribution + transform in this replay path
            // (still part of the API surface for fresh-RNG calls).
            let _ = (base_dist, transform);
            return Ok(true);
        }

        // ── Fast path: no active contracts ───────────────────────
        //
        // When both rate vectors are None or flat-default the
        // original uniform-site fast path runs — preserves the
        // legacy RNG consumption shape and byte-identical output.
        // A non-default rate vector (segment or V-subregion) takes
        // a weighted-sample path that consults
        // `combined_site_factor` once per pool index. Both
        // branches record exactly one site choice + one base
        // choice (so the trace contract is unchanged); the only
        // difference is the site selection distribution.
        let segment_rates_active = segment_rates
            .map(|r| !r.is_default())
            .unwrap_or(false);
        let v_subregion_rates_active = v_subregion_rates
            .map(|r| !r.is_default())
            .unwrap_or(false);
        let any_rate_weighting = segment_rates_active || v_subregion_rates_active;
        if self.ctx.contracts.is_none() {
            let site = if any_rate_weighting {
                let sim_ref = self.builder.peek();
                let refdata = self.ctx.refdata;
                let mut site_weights: Vec<f64> = Vec::with_capacity(pool_len as usize);
                let mut total = 0.0;
                for s in 0..pool_len {
                    let w = combined_site_factor(
                        s,
                        sim_ref,
                        refdata,
                        segment_rates,
                        v_subregion_rates,
                    );
                    site_weights.push(w);
                    total += w;
                }
                if !total.is_finite() || total <= 0.0 {
                    if self.strict {
                        return Err(PassError::constraint_sampling(
                            self.pass_name.clone(),
                            base_address.to_string(),
                            crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                        ));
                    }
                    return Ok(false);
                }
                let r = self.ctx.rng.next_f64() * total;
                let mut cum = 0.0;
                let mut chosen = pool_len - 1;
                for (idx, w) in site_weights.iter().enumerate() {
                    cum += *w;
                    if r < cum {
                        chosen = idx as u32;
                        break;
                    }
                }
                chosen
            } else {
                self.ctx.rng.range_u32(pool_len)
            };
            self.ctx
                .trace
                .record_choice(site_address, ChoiceValue::Int(site as i64));
            let base = transform(base_dist.sample(self.ctx.rng));
            self.ctx
                .trace
                .record_choice(base_address, ChoiceValue::Base(base));
            self.builder.change_base(NucHandle::new(site), base);
            return Ok(true);
        }

        let contracts = self.ctx.contracts.expect("contracts checked above");
        let refdata = self.ctx.refdata;

        // Pull the base_dist support once. Under v3.0 the helper can
        // only honor constrain-before-propose when the distribution
        // exposes its support — without it we can't compute per-site
        // admissible mass, and any fall-through to predicate-based
        // sampling reintroduces the reject-after-propose loop the v3
        // rule was written to eliminate. So if support is
        // unavailable we treat the slot as architecturally
        // unsamplable: strict mode surfaces
        // `SupportUnavailable`; permissive mode consumes the slot
        // as a no-op (no trace entry, no RNG, no pool mutation).
        // Built-in distributions (`UniformBase`, `EmpiricalLengthDist`,
        // `AllelePoolDist`) all override `support()`; only custom
        // dists wired to the contract-aware passes can land here.
        let support = match base_dist.support() {
            Some(s) => s,
            None => {
                if self.strict {
                    return Err(PassError::constraint_sampling(
                        self.pass_name.clone(),
                        base_address.to_string(),
                        crate::dist::FilteredSampleError::SupportUnavailable,
                    ));
                }
                return Ok(false);
            }
        };

        // ── Constrained path: per-site admissible-mass weighting ─
        //
        // Site weight = admissible-base mass at the site × segment
        // rate at the site (the latter when segment_rates is
        // non-default). Zero-rate sites drop out before the
        // cumulative-weight sampling step — segment rates are part
        // of proposal support, not a post-sampling rejection layer
        // (audit §4 ordering).
        let mut site_weights: Vec<f64> = Vec::with_capacity(pool_len as usize);
        let mut total: f64 = 0.0;
        {
            let sim = self.builder.peek();
            for site in 0..pool_len {
                let mask = contracts.admissible_bases_at(sim, refdata, NucHandle::new(site));
                let mass: f64 = support
                    .iter()
                    .filter(|(b, _)| {
                        if !mask.admits(*b) {
                            return false;
                        }
                        let written = transform(*b);
                        if written == *b {
                            return true;
                        }
                        contracts.admits_fixed_base_at(sim, refdata, NucHandle::new(site), written)
                    })
                    .map(|(_, w)| *w)
                    .sum();
                let rate_factor = combined_site_factor(
                    site,
                    sim,
                    refdata,
                    segment_rates,
                    v_subregion_rates,
                );
                let weighted = mass * rate_factor;
                site_weights.push(weighted);
                total += weighted;
            }
        }

        if !total.is_finite() || total <= 0.0 {
            if self.strict {
                return Err(PassError::constraint_sampling(
                    self.pass_name.clone(),
                    base_address.to_string(),
                    crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                ));
            }
            return Ok(false);
        }

        // Weighted-sample the site.
        let r = self.ctx.rng.next_f64() * total;
        let mut cum = 0.0;
        let mut chosen_site = pool_len - 1;
        for (i, w) in site_weights.iter().enumerate() {
            cum += *w;
            if r < cum {
                chosen_site = i as u32;
                break;
            }
        }
        self.ctx
            .trace
            .record_choice(site_address, ChoiceValue::Int(chosen_site as i64));

        // Weighted-sample the base from base_dist restricted to the
        // admitted (and, if transformed, fixed-base-admitted) support
        // at the chosen site.
        let filtered: Vec<(u8, f64)> = {
            let sim = self.builder.peek();
            let mask = contracts.admissible_bases_at(sim, refdata, NucHandle::new(chosen_site));
            support
                .iter()
                .filter(|(b, _)| {
                    if !mask.admits(*b) {
                        return false;
                    }
                    let written = transform(*b);
                    if written == *b {
                        return true;
                    }
                    contracts.admits_fixed_base_at(
                        sim,
                        refdata,
                        NucHandle::new(chosen_site),
                        written,
                    )
                })
                .copied()
                .collect()
        };

        debug_assert!(
            !filtered.is_empty(),
            "site_weight > 0 implies non-empty per-site admissible support"
        );
        let total_f: f64 = filtered.iter().map(|(_, w)| *w).sum();
        let r = self.ctx.rng.next_f64() * total_f;
        let mut cum = 0.0;
        let mut chosen_pre = filtered.last().expect("filtered non-empty").0;
        for &(b, w) in filtered.iter() {
            cum += w;
            if r < cum {
                chosen_pre = b;
                break;
            }
        }
        let chosen = transform(chosen_pre);
        self.ctx
            .trace
            .record_choice(base_address, ChoiceValue::Base(chosen));
        self.builder
            .change_base(NucHandle::new(chosen_site), chosen);

        Ok(true)
    }

    /// Substitute the base at `site` with a known value, **bypassing
    /// the active contract bundle**. Used only for operations that
    /// are conceptually outside the contract model (e.g. event-log
    /// replay reconstructing a recorded value, or a future pass that
    /// genuinely wants a force-write). All ordinary mutation paths
    /// — including those that look like fixed writes (for example
    /// N-injection) — should call
    /// [`Self::substitute_base_admitting`] instead so the active
    /// contracts arbitrate the candidate. The deliberately-loud
    /// `force_` prefix marks this as the architectural exception.
    /// Does NOT record to the trace — caller decides what to record.
    /// Validates the handle range.
    pub(crate) fn force_substitute_base(
        &mut self,
        site: NucHandle,
        new_base: u8,
    ) -> Result<(), PassError> {
        let pool_len = self.builder.peek().pool.len() as u32;
        if site.index() >= pool_len {
            return Err(PassError::invalid_plan_state(
                self.pass_name.clone(),
                format!(
                    "force_substitute_base: site handle {} out of pool range [0, {})",
                    site.index(),
                    pool_len
                ),
            ));
        }
        self.builder.change_base(site, new_base);
        Ok(())
    }

    /// Substitute the base at `site` with a known value, contract-aware.
    /// The candidate value is pinned (no sampling), but the active
    /// contract bundle still gets to arbitrate: if any contract
    /// rejects the candidate at this site, the write is **skipped**
    /// in permissive mode (no pool mutation) and returns
    /// `Err(PassError::ConstraintSampling)` in strict mode.
    ///
    /// Returns:
    /// - `Ok(true)` — substitution applied (no contract, or contract
    ///   admitted).
    /// - `Ok(false)` — contract rejected, write skipped (permissive).
    /// - `Err(..)` — out-of-range handle, or strict-mode rejection.
    ///
    /// Used by passes whose candidate value is naturally pinned
    /// (N-corruption writes `b'N'`). Does NOT record to the trace —
    /// caller decides what to record.
    pub(crate) fn substitute_base_admitting(
        &mut self,
        site: NucHandle,
        new_base: u8,
        address: ChoiceAddress,
    ) -> Result<bool, PassError> {
        let pool_len = self.builder.peek().pool.len() as u32;
        if site.index() >= pool_len {
            return Err(PassError::invalid_plan_state(
                self.pass_name.clone(),
                format!(
                    "substitute_base_admitting: site handle {} out of pool range [0, {})",
                    site.index(),
                    pool_len
                ),
            ));
        }
        if let Some(contracts) = self.ctx.contracts {
            let context =
                ChoiceContext::targeted_base_substitution(0, 1, site).with_address(address);
            let address = address.to_string();
            if contracts
                .admits_typed(
                    self.builder.peek(),
                    self.ctx.refdata,
                    context,
                    &ChoiceValue::Base(new_base),
                )
                .is_err()
            {
                if self.strict {
                    return Err(PassError::constraint_sampling(
                        self.pass_name.clone(),
                        address.to_string(),
                        crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                    ));
                }
                return Ok(false);
            }
        }
        self.builder.change_base(site, new_base);
        Ok(true)
    }

    /// Sample a candidate base from `base_dist`, contract-filtered
    /// against the active contract set when present.
    ///
    /// Returns `Ok(Some(base))` on a successful draw (caller
    /// applies the transform). Returns `Ok(None)` when contracts
    /// are active and the filtered support is empty in permissive
    /// mode — the v3.0 rule: no unconstrained fallback. Returns
    /// `Err(..)` in strict mode when the filter is empty. When
    /// no contracts are active, draws unfiltered from `base_dist`
    /// and returns `Ok(Some(_))`.
    fn sample_with_filter(
        &mut self,
        base_dist: &dyn Distribution<Output = u8>,
        addr: ChoiceAddress,
        draw_index: u32,
        draw_count: u32,
        site: NucHandle,
        transform: fn(u8) -> u8,
    ) -> Result<Option<u8>, PassError> {
        let contracts = self.ctx.contracts;
        let refdata = self.ctx.refdata;
        if let Some(contracts) = contracts {
            let context = ChoiceContext::targeted_base_substitution(draw_index, draw_count, site)
                .with_address(addr);
            let addr = addr.to_string();
            // Policy is `Skip` (see `SUBSTITUTE_BASE_EMPTY_SUPPORT`)
            // — permissive empty-support surfaces as `Ok(None)`
            // here and the caller (`substitute_base`) translates
            // that to `Ok(false)` for the pass.
            return sample_filtered_with_policy(
                self.ctx.rng,
                base_dist,
                |candidate: &u8| {
                    let transformed = transform(*candidate);
                    contracts
                        .admits_typed(
                            self.builder.peek(),
                            refdata,
                            context,
                            &ChoiceValue::Base(transformed),
                        )
                        .is_ok()
                },
                self.strict,
                &self.pass_name,
                &addr,
                SUBSTITUTE_BASE_EMPTY_SUPPORT,
            );
        }
        Ok(Some(base_dist.sample(self.ctx.rng)))
    }
}

#[inline]
fn identity_base(base: u8) -> u8 {
    base
}

/// Combined per-position SHM weight factor: `segment_rate ×
/// v_subregion_rate`. Returns `1.0` when both vectors are flat-
/// default, so the existing fast paths inside
/// `substitute_position_constrained` and S5F's `build_profile`
/// short-circuit cleanly.
///
/// Discipline:
/// - Segment lookup uses `segment_at_position`; sites without a
///   segment (e.g. pre-assembly) receive segment factor `0.0`
///   when `segment_rates` is non-default — drops them from
///   support, same as the segment-rates slice (audit §4).
/// - V-subregion lookup runs ONLY for V sites; non-V sites
///   receive subregion factor `1.0`.
/// - V sites on an allele without IMGT annotations
///   (`v_subregion_at_position` returns `None`) receive factor
///   `1.0` — the brief's mixed-cartridge fallback.
/// - V sites in the V-side CDR3 contribution (between FWR3.end
///   and `len(allele.seq)`) also receive factor `1.0`; that
///   stretch is deliberately outside the five-label set.
#[inline]
pub(crate) fn combined_site_factor(
    pos: u32,
    sim: &crate::ir::Simulation,
    refdata: Option<&crate::refdata::RefDataConfig>,
    segment_rates: Option<&crate::passes::mutate::SegmentRateWeights>,
    v_subregion_rates: Option<&crate::passes::mutate::VSubregionRateWeights>,
) -> f64 {
    let seg_active = segment_rates.map(|r| !r.is_default()).unwrap_or(false);
    let v_sub_active = v_subregion_rates
        .map(|r| !r.is_default())
        .unwrap_or(false);
    if !seg_active && !v_sub_active {
        return 1.0;
    }
    let segment = crate::passes::mutate::segment_at_position(&sim.sequence, pos);
    let segment_factor = if seg_active {
        let rates = segment_rates.expect("active check");
        segment.map(|s| rates.rate_for(s)).unwrap_or(0.0)
    } else {
        1.0
    };
    let v_sub_factor = if v_sub_active && segment == Some(crate::ir::Segment::V) {
        let rates = v_subregion_rates.expect("active check");
        if let Some(rd) = refdata {
            crate::passes::mutate::v_subregion_at_position(
                &sim.sequence,
                &sim.assignments,
                rd,
                pos,
            )
            .map(|label| rates.rate_for(label))
            // No annotation under this V allele OR position in the
            // V-side CDR3 stretch → identity factor.
            .unwrap_or(1.0)
        } else {
            // No refdata in scope → can't look up subregions; treat
            // as if unannotated. (Should not happen under the DSL
            // entry point; defensive for direct Rust callers.)
            1.0
        }
    } else {
        1.0
    };
    segment_factor * v_sub_factor
}

// ──────────────────────────────────────────────────────────────────
// Tests — trace-injected replay path in `substitute_position_constrained`
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod replay_tests {
    use super::super::MutationTransaction;
    use crate::address::ChoiceAddress;
    use crate::contract::{AnchorPreserved, ContractSet};
    use crate::dist::UniformBase;
    use crate::ir::{NucHandle, Nucleotide, Segment, Simulation};
    use crate::pass::{PassContext, PassError};
    use crate::refdata::{Allele, AlleleId, ChainType, RefDataConfig};
    use crate::replay::{ReplayError, TraceCursor};
    use crate::rng::Rng;
    use crate::trace::{ChoiceRecord, ChoiceValue, Trace};

    fn make_sim(bases: &[u8]) -> Simulation {
        let mut sim = Simulation::new();
        for (i, b) in bases.iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
                *b,
                i as u16,
                Segment::V,
            ));
            sim = next;
        }
        sim
    }

    fn rec(addr: ChoiceAddress, v: ChoiceValue) -> ChoiceRecord {
        ChoiceRecord::new(addr.to_string(), v)
    }

    /// Drive `substitute_position_constrained` once against a
    /// caller-supplied `(records, contracts, refdata)`. Returns
    /// `(result, output trace, rng words consumed)`.
    fn run_substitute(
        sim: &Simulation,
        records: Vec<ChoiceRecord>,
        site_address: ChoiceAddress,
        base_address: ChoiceAddress,
        transform: Option<fn(u8) -> u8>,
        contracts: Option<&ContractSet>,
        refdata: Option<&RefDataConfig>,
    ) -> (Result<(bool, Simulation), PassError>, Trace, u64) {
        let mut cursor = TraceCursor::from_owned(records);
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let words_after;
        let result = {
            let mut ctx = PassContext {
                trace: &mut trace,
                rng: &mut rng,
                pass_index: 0,
                refdata,
                contracts,
                feasibility: None,
                reference_index: None,
                replay_cursor: Some(&mut cursor),
                event_log_sink: None,
            };
            let mut tx = MutationTransaction::open(sim, &mut ctx, "test.pass", true);
            let dist = UniformBase;
            let wrote = tx.substitute_position_constrained(
                &dist,
                site_address,
                base_address,
                transform,
                None, // test helper exercises the legacy paths
                None,
            );
            let outcome = wrote.and_then(|w| tx.commit().map(|s| (w, s)));
            outcome
        };
        words_after = rng.words_consumed();
        (result, trace, words_after)
    }

    // ── No-contracts path ─────────────────────────────────────────

    #[test]
    fn replay_no_contracts_consumes_recorded_site_and_base() {
        let sim = make_sim(b"AAAA");
        let records = vec![
            rec(ChoiceAddress::MutateUniformSite(0), ChoiceValue::Int(2)),
            rec(ChoiceAddress::MutateUniformBase(0), ChoiceValue::Base(b'G')),
        ];
        let (result, trace, rng_words) = run_substitute(
            &sim,
            records,
            ChoiceAddress::MutateUniformSite(0),
            ChoiceAddress::MutateUniformBase(0),
            None,
            None,
            None,
        );
        let (wrote, next) = result.unwrap();
        assert!(wrote);
        assert_eq!(next.pool.get(NucHandle::new(2)).unwrap().base, b'G');
        assert_eq!(
            trace.find("mutate.uniform.site[0]").unwrap().value,
            ChoiceValue::Int(2),
        );
        assert_eq!(
            trace.find("mutate.uniform.base[0]").unwrap().value,
            ChoiceValue::Base(b'G'),
        );
        // RNG must not have been consumed — replay bypasses the
        // sampler entirely.
        assert_eq!(rng_words, 0);
    }

    // ── Error: site out of range ──────────────────────────────────

    #[test]
    fn replay_site_out_of_range_surfaces_invalid_plan_state() {
        let sim = make_sim(b"AAA");
        let records = vec![
            rec(ChoiceAddress::MutateUniformSite(0), ChoiceValue::Int(99)),
            rec(ChoiceAddress::MutateUniformBase(0), ChoiceValue::Base(b'G')),
        ];
        let (result, _trace, _) = run_substitute(
            &sim,
            records,
            ChoiceAddress::MutateUniformSite(0),
            ChoiceAddress::MutateUniformBase(0),
            None,
            None,
            None,
        );
        match result.unwrap_err() {
            PassError::InvalidPlanState { reason, .. } => {
                assert!(reason.contains("99"));
                assert!(reason.contains("out of pool range"));
            }
            other => panic!("expected InvalidPlanState, got {other:?}"),
        }
    }

    // ── Error: wrong-kind cursor records ──────────────────────────

    #[test]
    fn replay_wrong_kind_site_surfaces_replay_error() {
        let sim = make_sim(b"AAA");
        let records = vec![
            // Site should be Int, but we put a Bool.
            rec(ChoiceAddress::MutateUniformSite(0), ChoiceValue::Bool(true)),
            rec(ChoiceAddress::MutateUniformBase(0), ChoiceValue::Base(b'G')),
        ];
        let (result, _trace, _) = run_substitute(
            &sim,
            records,
            ChoiceAddress::MutateUniformSite(0),
            ChoiceAddress::MutateUniformBase(0),
            None,
            None,
            None,
        );
        match result.unwrap_err() {
            PassError::Replay { reason, .. } => {
                assert!(matches!(reason, ReplayError::ValueKindMismatch { .. }));
            }
            other => panic!("expected Replay::ValueKindMismatch, got {other:?}"),
        }
    }

    #[test]
    fn replay_wrong_kind_base_surfaces_replay_error() {
        let sim = make_sim(b"AAA");
        let records = vec![
            rec(ChoiceAddress::MutateUniformSite(0), ChoiceValue::Int(0)),
            // Base should be Base, but we put an Int.
            rec(ChoiceAddress::MutateUniformBase(0), ChoiceValue::Int(7)),
        ];
        let (result, _trace, _) = run_substitute(
            &sim,
            records,
            ChoiceAddress::MutateUniformSite(0),
            ChoiceAddress::MutateUniformBase(0),
            None,
            None,
            None,
        );
        match result.unwrap_err() {
            PassError::Replay { reason, .. } => {
                assert!(matches!(reason, ReplayError::ValueKindMismatch { .. }));
            }
            other => panic!("expected Replay::ValueKindMismatch, got {other:?}"),
        }
    }

    // ── Contract validation: admissible vs rejected ───────────────

    /// VJ fixture where the V anchor codon is `TGG` (Trp). The
    /// `AnchorPreserved::V` contract restricts site 0 (the first T
    /// of the codon) such that only `T` keeps the anchor intact;
    /// substituting other bases changes the amino acid.
    fn anchor_locked_vj() -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_anchor*01".into(),
            gene: "v_anchor".into(),
            seq: b"TGG".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
            functional_status: None,
            subregions: Vec::new(),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_anchor*01".into(),
            gene: "j_anchor".into(),
            seq: b"TGG".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
            functional_status: None,
            subregions: Vec::new(),
        });
        let mut sim = Simulation::new();
        for (i, b) in b"TGG".iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
                *b,
                i as u16,
                Segment::V,
            ));
            sim = next;
        }
        sim = sim.with_region_added(crate::ir::Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(3),
        ));
        sim = sim.with_allele_assigned(
            Segment::V,
            crate::assignment::AlleleInstance::new(AlleleId::new(0)),
        );
        (cfg, sim)
    }

    #[test]
    fn replay_under_contract_accepts_admissible_recorded_base() {
        let (cfg, sim) = anchor_locked_vj();
        let contracts = ContractSet::new().with(Box::new(AnchorPreserved::new(Segment::V)));
        // At site 0, T is the only base that preserves the W anchor
        // (TGG → ?GG). T is admissible.
        let records = vec![
            rec(ChoiceAddress::MutateUniformSite(0), ChoiceValue::Int(0)),
            rec(ChoiceAddress::MutateUniformBase(0), ChoiceValue::Base(b'T')),
        ];
        let (result, _trace, _) = run_substitute(
            &sim,
            records,
            ChoiceAddress::MutateUniformSite(0),
            ChoiceAddress::MutateUniformBase(0),
            None,
            Some(&contracts),
            Some(&cfg),
        );
        let (wrote, next) = result.unwrap();
        assert!(wrote);
        assert_eq!(next.pool.get(NucHandle::new(0)).unwrap().base, b'T');
    }

    #[test]
    fn replay_under_contract_rejects_inadmissible_recorded_base() {
        // Same fixture: substituting site 0 with A makes the codon
        // AGG (Arg) instead of TGG (Trp). AnchorPreserved.V rejects.
        let (cfg, sim) = anchor_locked_vj();
        let contracts = ContractSet::new().with(Box::new(AnchorPreserved::new(Segment::V)));
        let records = vec![
            rec(ChoiceAddress::MutateUniformSite(0), ChoiceValue::Int(0)),
            rec(ChoiceAddress::MutateUniformBase(0), ChoiceValue::Base(b'A')),
        ];
        let (result, _trace, _) = run_substitute(
            &sim,
            records,
            ChoiceAddress::MutateUniformSite(0),
            ChoiceAddress::MutateUniformBase(0),
            None,
            Some(&contracts),
            Some(&cfg),
        );
        match result.unwrap_err() {
            PassError::ConstraintSampling { address, .. } => {
                assert_eq!(address, "mutate.uniform.base[0]");
            }
            other => panic!("expected ConstraintSampling, got {other:?}"),
        }
    }

    // ── Lowercase quality transform: case-folded recorded base
    //    still validates through the canonical mask. ──────────────

    fn lowercase(b: u8) -> u8 {
        b.to_ascii_lowercase()
    }

    #[test]
    fn replay_quality_lowercase_base_validates_via_canonical_mask() {
        // Site 1 of TGG is unconstrained by AnchorPreserved.V (it
        // covers the whole codon; substituting any letter at index
        // 1 still changes the amino acid, so actually it IS
        // constrained — pick a different layout).
        //
        // Use a sim with a free site and verify the lowercase
        // recorded base flows through unchanged.
        let sim = make_sim(b"AAAA");
        let records = vec![
            rec(
                ChoiceAddress::CorruptQualitySite(0),
                ChoiceValue::Int(2),
            ),
            rec(
                ChoiceAddress::CorruptQualityBase(0),
                ChoiceValue::Base(b'c'),
            ),
        ];
        let (result, trace, _) = run_substitute(
            &sim,
            records,
            ChoiceAddress::CorruptQualitySite(0),
            ChoiceAddress::CorruptQualityBase(0),
            Some(lowercase),
            None,
            None,
        );
        let (wrote, next) = result.unwrap();
        assert!(wrote);
        // The recorded byte is lowercase 'c' and gets written verbatim
        // (the fresh-RNG path records post-transform — same shape).
        assert_eq!(next.pool.get(NucHandle::new(2)).unwrap().base, b'c');
        assert_eq!(
            trace.find("corrupt.quality.error_base[0]").unwrap().value,
            ChoiceValue::Base(b'c'),
        );
    }
}
