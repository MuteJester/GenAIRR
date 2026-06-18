//! `ReceptorRevisionPass` — post-recombination V-segment replacement
//! (receptor revision Slice C).
//!
//! Models the biological process by which an already-assembled
//! receptor undergoes a second VDJ-recombination-like event that
//! swaps the V segment for a different germline allele while leaving
//! the D/J/NP regions intact. v1 is intentionally narrow:
//!
//! - **VDJ chains only.** The pass declares
//!   `PassRequirement::AlleleAssignment(V)`; downstream sanity is the
//!   caller's job (Slice D's DSL gate).
//! - **V replacement only.** D / J / C / NP replacement are not in
//!   scope; the same persistent `replace_segment` primitive would
//!   cover them, but Slice C ships the narrow case.
//! - **Same-length retained slice.** The candidate V allele must
//!   satisfy `allele.len() >= old_v_region.len()`; the 3' trim is
//!   derived so that `allele.len() - trim_3 == old_v_region.len()`
//!   (5' trim stays at 0). This means downstream pool positions don't
//!   shift, and the AIRR validator's region-length invariants stay
//!   intact without a Slice E migration.
//!
//! ## Pipeline position
//!
//! The pass runs **after** initial V recombination (sample+assemble)
//! and **before** SHM / corruption. The schedule analyser orders it
//! after `AssembleSegmentPass(V)` via the
//! `PassRequirement::AlleleAssignment(V)` chain; placing it before
//! mutation passes is the caller's responsibility (Slice D wires
//! this in the DSL).
//!
//! ## Trace addresses (Slice C)
//!
//! - `receptor_revision.applied` — `Bool`. Always recorded.
//! - `receptor_revision.v_allele` — `AlleleId`. Only when applied.
//! - `receptor_revision.v_trim_3` — `Int`. Only when applied.
//!
//! ## Event stream (when applied)
//!
//! In order, through the same `SimulationBuilder`:
//!
//! 1. `AssignmentChanged { V, old: Some(prev), new: new_inst }`
//! 2. `TrimChanged { V, Three, old: Some(prev.trim_3), new: derived_trim_3 }`
//! 3. `SegmentReplaced { V, old_region, new_region, bytes_delta: 0 }`
//!
//! The pass is wired through the `receptor_revision()` DSL method and
//! the AIRR record exposes `receptor_revision_applied` + `original_v_call`
//! (pre-revision V); `v_call` is the post-revision V. A genotype-aware
//! variant (see `GenotypeVConstraint`) restricts the replacement V to the
//! carried alleles on the drawn rearrangement chromosome.

use crate::address;
use crate::assignment::{AlleleInstance, TrimEnd};
use crate::contract::ChoiceContext;
use crate::dist::Distribution;
use crate::ir::{Nucleotide, Segment, Simulation, SimulationBuilder};
use crate::pass::{Pass, PassContext, PassError, PassRequirement};
use crate::refdata::{Allele, AlleleId, RefDataConfig};
use crate::trace::ChoiceValue;

/// Per-haplotype carried-V candidates for genotype-aware receptor revision.
/// `per_hap[c]` is the list of `(allele, mass)` (mass = weight * copies)
/// carried on chromosome `c`. `same_haplotype` restricts the draw to the
/// rearrangement chromosome; otherwise both haplotypes are aggregated.
pub struct GenotypeVConstraint {
    pub per_hap: [Vec<(AlleleId, f64)>; 2],
    pub same_haplotype: bool,
}

pub struct ReceptorRevisionPass {
    prob: f64,
    v_distribution: Box<dyn Distribution<Output = AlleleId>>,
    genotype_v: Option<GenotypeVConstraint>,
}

impl ReceptorRevisionPass {
    /// Construct a receptor-revision pass.
    ///
    /// `prob` is the per-simulation probability of applying the
    /// revision. `v_distribution` produces replacement V allele
    /// candidates; the pass filters to those whose retained slice
    /// matches the old V region's length.
    ///
    /// Panics if `prob` is not finite or falls outside `[0.0, 1.0]`.
    pub fn new(prob: f64, v_distribution: Box<dyn Distribution<Output = AlleleId>>) -> Self {
        assert!(
            prob.is_finite() && (0.0..=1.0).contains(&prob),
            "ReceptorRevisionPass: prob must be in [0.0, 1.0], got {}",
            prob,
        );
        Self {
            prob,
            v_distribution,
            genotype_v: None,
        }
    }

    /// Read-only accessor for tests / report builders.
    pub fn prob(&self) -> f64 {
        self.prob
    }

    /// Attach a genotype-V constraint, switching the pass to genotype-aware
    /// candidate selection (carried alleles on the drawn chromosome, current
    /// V excluded). Builder — returns `self`.
    #[must_use]
    pub fn with_genotype_constraint(mut self, constraint: GenotypeVConstraint) -> Self {
        self.genotype_v = Some(constraint);
        self
    }

    /// Aggregated, current-excluded, same-length-eligible candidate pool for
    /// genotype-aware revision. Aggregates mass by `AlleleId` (summing across
    /// haplotypes when `!same_haplotype`, and collapsing duplicates within a
    /// haplotype), drops the currently-assigned `current_v`, and keeps only
    /// alleles whose length admits the retained slice.
    fn genotype_eligible(
        &self,
        constraint: &GenotypeVConstraint,
        c: usize,
        current_v: AlleleId,
        refdata: &RefDataConfig,
        old_v_len: u32,
    ) -> Vec<(AlleleId, f64)> {
        use std::collections::BTreeMap;
        let mut agg: BTreeMap<u32, f64> = BTreeMap::new();
        let haps: &[Vec<(AlleleId, f64)>] = if constraint.same_haplotype {
            std::slice::from_ref(&constraint.per_hap[c])
        } else {
            &constraint.per_hap
        };
        for list in haps {
            for (id, mass) in list {
                if *mass > 0.0 {
                    *agg.entry(id.index()).or_insert(0.0) += *mass;
                }
            }
        }
        agg.into_iter()
            .filter(|(idx, _)| *idx != current_v.index())
            .filter(|(idx, _)| {
                refdata
                    .get(Segment::V, AlleleId::new(*idx))
                    .map(|a| a.len() >= old_v_len)
                    .unwrap_or(false)
            })
            .map(|(idx, mass)| (AlleleId::new(idx), mass))
            .collect()
    }

    fn execute_genotype_aware(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
        constraint: &GenotypeVConstraint,
    ) -> Result<Simulation, PassError> {
        if !sim.assignments.has(Segment::V) {
            return Err(PassError::missing_assignment(self.name(), Segment::V));
        }
        let refdata = ctx
            .refdata
            .ok_or_else(|| PassError::missing_refdata(self.name()))?;
        let old_v_len = self.old_v_region_len(sim)?;
        let current = sim.assignments.get(Segment::V).expect("V assigned");
        let current_v = current.allele_id;
        let c = current.haplotype.ok_or_else(|| {
            PassError::invalid_plan_state(
                self.name(),
                "genotype-aware receptor revision requires a haplotype-stamped V \
                 assignment (SampleGenotypePass must run first)",
            )
        })? as usize;

        let replaying = ctx.replay_cursor.is_some();
        let coin = if replaying {
            ctx.replay_cursor
                .as_deref_mut()
                .expect("cursor")
                .expect_bool(address::ChoiceAddress::ReceptorRevisionApplied)
                .map_err(|r| PassError::replay(self.name(), r))?
        } else {
            ctx.rng.next_f64() < self.prob
        };

        if replaying {
            ctx.trace.record_choice(
                address::ChoiceAddress::ReceptorRevisionApplied,
                ChoiceValue::Bool(coin),
            );
            if !coin {
                return Ok(sim.clone());
            }
            let (new_id, derived_trim_3) = self.consume_replay_choices_genotype(
                constraint, c, current_v, refdata, old_v_len, ctx,
            )?;
            return self.finish_apply(sim, refdata, new_id, derived_trim_3, c, strict, ctx);
        }

        // Fresh: applied = coin AND an eligible alternate exists.
        let eligible = self.genotype_eligible(constraint, c, current_v, refdata, old_v_len);
        if coin && eligible.is_empty() {
            if strict {
                return Err(PassError::constraint_sampling(
                    self.name(),
                    address::ChoiceAddress::ReceptorRevisionVAllele.to_string(),
                    crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                ));
            }
            ctx.trace.record_choice(
                address::ChoiceAddress::ReceptorRevisionApplied,
                ChoiceValue::Bool(false),
            );
            return Ok(sim.clone());
        }
        let applied = coin && !eligible.is_empty();
        ctx.trace.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(applied),
        );
        if !applied {
            return Ok(sim.clone());
        }
        let total: f64 = eligible.iter().map(|(_, m)| m).sum();
        let mut u = ctx.rng.next_f64() * total;
        let mut new_id = eligible.last().expect("eligible non-empty").0;
        for (id, m) in &eligible {
            u -= *m;
            if u <= 0.0 {
                new_id = *id;
                break;
            }
        }
        let new_allele = refdata
            .get(Segment::V, new_id)
            .ok_or_else(|| PassError::missing_allele(self.name(), Segment::V, new_id.index()))?;
        let derived_trim_3 = new_allele.len() - old_v_len;
        self.finish_apply(sim, refdata, new_id, derived_trim_3, c, strict, ctx)
    }

    /// Record the typed choices + commit the replacement (shared by the fresh
    /// and replay genotype paths). Stamps the original rearrangement haplotype
    /// `c` onto the replacement instance.
    fn finish_apply(
        &self,
        sim: &Simulation,
        refdata: &RefDataConfig,
        new_id: AlleleId,
        derived_trim_3: u32,
        c: usize,
        strict: bool,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        if derived_trim_3 > u16::MAX as u32 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::ReceptorRevisionVTrim3.to_string(),
                derived_trim_3 as i64,
                "trim_exceeds_u16",
            ));
        }
        let new_allele = refdata
            .get(Segment::V, new_id)
            .ok_or_else(|| PassError::missing_allele(self.name(), Segment::V, new_id.index()))?
            .clone();
        ctx.trace.record_choice(
            address::ChoiceAddress::ReceptorRevisionVAllele,
            ChoiceValue::AlleleId(new_id.index()),
        );
        ctx.trace.record_choice(
            address::ChoiceAddress::ReceptorRevisionVTrim3,
            ChoiceValue::Int(derived_trim_3 as i64),
        );
        let old_v_len = self.old_v_region_len(sim)?;
        // Strict post-event contract arbitration — parity with the non-genotype
        // path: build a hypothetical post-replacement IR (no trace/sink side
        // effects) and ask the active contracts whether the receptor-revised
        // state is admissible (e.g. productive_only). Permissive skips this.
        if strict && ctx.contracts.is_some() {
            let mut hypothetical_ctx = PassContext {
                trace: ctx.trace,
                rng: ctx.rng,
                pass_index: ctx.pass_index,
                refdata: ctx.refdata,
                contracts: ctx.contracts,
                feasibility: ctx.feasibility,
                reference_index: ctx.reference_index,
                replay_cursor: None,
                event_log_sink: None,
            };
            let hypothetical = self.commit_replacement(
                sim.clone(),
                &new_allele,
                new_id,
                derived_trim_3 as u16,
                old_v_len,
                Some(c as u8),
                &mut hypothetical_ctx,
            );
            self.validate_post_event_contracts(sim, &hypothetical, ctx)?;
        }
        Ok(self.commit_replacement(
            sim.clone(),
            &new_allele,
            new_id,
            derived_trim_3 as u16,
            old_v_len,
            Some(c as u8),
            ctx,
        ))
    }

    /// Replay validation for genotype-aware revision: the recorded
    /// `(v_allele, v_trim_3)` must be carried by the allowed pool (per
    /// `same_haplotype`), differ from the current V, resolve in refdata, and
    /// retain exactly `old_v_len` bytes. "Trace proposes, engine validates."
    fn consume_replay_choices_genotype(
        &self,
        constraint: &GenotypeVConstraint,
        c: usize,
        current_v: AlleleId,
        refdata: &RefDataConfig,
        old_v_len: u32,
        ctx: &mut PassContext,
    ) -> Result<(AlleleId, u32), PassError> {
        let id_index = ctx
            .replay_cursor
            .as_deref_mut()
            .expect("cursor")
            .expect_allele_id(address::ChoiceAddress::ReceptorRevisionVAllele)
            .map_err(|r| PassError::replay(self.name(), r))?;
        let trim_i64 = ctx
            .replay_cursor
            .as_deref_mut()
            .expect("cursor")
            .expect_int(address::ChoiceAddress::ReceptorRevisionVTrim3)
            .map_err(|r| PassError::replay(self.name(), r))?;
        let id = AlleleId::new(id_index);

        // Resolve the allele FIRST so an unresolvable (out-of-range) recorded id
        // surfaces as `missing_allele` — matching the non-genotype replay path's
        // diagnostic contract — before the eligible/current/length checks.
        let allele = refdata
            .get(Segment::V, id)
            .ok_or_else(|| PassError::missing_allele(self.name(), Segment::V, id_index))?;

        let eligible = self.genotype_eligible(constraint, c, current_v, refdata, old_v_len);
        if !eligible.iter().any(|(eid, _)| eid.index() == id_index) {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::ReceptorRevisionVAllele.to_string(),
                id_index as i64,
                if id_index == current_v.index() {
                    "receptor_revision_allele_equals_current"
                } else {
                    "receptor_revision_allele_not_carried_on_haplotype"
                },
            ));
        }
        if trim_i64 < 0 || trim_i64 > u16::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::ReceptorRevisionVTrim3.to_string(),
                trim_i64,
                "trim_out_of_range",
            ));
        }
        let trim_3 = trim_i64 as u32;
        let retained = (allele.len() as i64).checked_sub(trim_3 as i64).unwrap_or(-1);
        if retained != old_v_len as i64 {
            return Err(PassError::invalid_plan_state(
                self.name(),
                format!(
                    "replay length mismatch: allele {} length {} trim_3 {} retains {}, expected {}",
                    allele.name,
                    allele.len(),
                    trim_3,
                    retained,
                    old_v_len
                ),
            ));
        }
        Ok((id, trim_3))
    }

    /// Resolve the single V region in `sim`. Receptor revision v1
    /// requires exactly one — zero means assembly never ran, more
    /// than one violates the one-region-per-segment invariant the
    /// audit pinned (`MultipleRegionsForSegment`).
    fn old_v_region_len(&self, sim: &Simulation) -> Result<u32, PassError> {
        let mut iter = sim.sequence.regions.iter().filter(|r| r.segment == Segment::V);
        let first = iter.next().ok_or_else(|| {
            PassError::invalid_plan_state(
                self.name(),
                "no V region in sequence — AssembleSegmentPass(V) must run first",
            )
        })?;
        if iter.next().is_some() {
            return Err(PassError::invalid_plan_state(
                self.name(),
                "multiple V regions in sequence — receptor revision is defined for one-region V only",
            ));
        }
        Ok(first.len())
    }

    /// Sample a replacement V allele constrained to same-retained-
    /// length candidates. Strict mode errors on empty support;
    /// permissive mode falls back to the raw distribution draw
    /// (mirrors `SampleAllelePass`'s v3.0 documented exception for
    /// allele sampling — a "skip" semantics would leave the V
    /// segment in an inconsistent state).
    fn sample_replacement_allele(
        &self,
        refdata: &RefDataConfig,
        old_v_len: u32,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<AlleleId, PassError> {
        let support = self.v_distribution.support().ok_or_else(|| {
            PassError::invalid_plan_state(
                self.name(),
                "V distribution does not enumerate support; \
                 receptor revision requires an enumerable distribution",
            )
        })?;

        let eligible: Vec<(AlleleId, f64)> = support
            .into_iter()
            .filter(|(id, weight)| {
                if *weight <= 0.0 {
                    return false;
                }
                refdata
                    .get(Segment::V, *id)
                    .map(|allele| allele.len() >= old_v_len)
                    .unwrap_or(false)
            })
            .collect();

        if eligible.is_empty() {
            if strict {
                return Err(PassError::constraint_sampling(
                    self.name(),
                    address::ChoiceAddress::ReceptorRevisionVAllele.to_string(),
                    crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                ));
            }
            // Permissive fall-through: take the raw natural draw
            // even if it's length-ineligible. Caller will surface
            // the mismatch downstream when the post-state contract
            // check (or in v1, the slice() overflow) fires.
            return Ok(self.v_distribution.sample(ctx.rng));
        }

        let total_weight: f64 = eligible.iter().map(|(_, w)| w).sum();
        let mut u = ctx.rng.next_f64() * total_weight;
        for (id, weight) in &eligible {
            u -= weight;
            if u <= 0.0 {
                return Ok(*id);
            }
        }
        Ok(eligible.last().expect("eligible non-empty").0)
    }

    /// Apply the replacement to the simulation through a single
    /// builder so the three consequence events (`AssignmentChanged`,
    /// `TrimChanged`, `SegmentReplaced`) fan out to every attached
    /// sink in canonical order.
    ///
    /// Returns the post-replacement `Simulation`. The caller is
    /// responsible for forwarding the captured event log to
    /// `ctx.event_log_sink` if one is attached.
    fn commit_replacement(
        &self,
        sim: Simulation,
        new_allele: &Allele,
        new_id: AlleleId,
        derived_trim_3: u16,
        old_v_len: u32,
        haplotype: Option<u8>,
        ctx: &mut PassContext,
    ) -> Simulation {
        // Capture the pre-revision V identity **before** the builder
        // mutates the slot. This is the IR-side source of truth for
        // the AIRR `original_v_call` projection — read by
        // `airr_record::builder` from the post-replacement
        // `assignments.v.receptor_revision_original_id`. If the
        // current V already carries a preserved original id (a
        // hypothetical chained revision), keep the *first* original
        // rather than overwriting with the intermediate identity.
        // Today the DSL rejects chained revisions, but the property
        // stays stable under any future relaxation.
        let preserved_original_id = sim
            .assignments
            .get(Segment::V)
            .map(|inst| inst.receptor_revision_original_id.unwrap_or(inst.allele_id))
            .expect(
                "ReceptorRevisionPass::commit_replacement called \
                 without a V assignment; \
                 execute_with_sampling_mode must guarantee this",
            );

        let mut builder = SimulationBuilder::from_simulation(sim);
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }

        // 1. Reassign V to the new allele, threading the pre-
        //    revision id into the persistent provenance slot. Trim
        //    resets to 0 here; `update_trim` below installs the
        //    derived value but preserves every other field on the
        //    instance (including the provenance) per the
        //    `AlleleInstance::with_trim_3` contract.
        // Preserve the rearrangement-chromosome provenance: receptor revision
        // does not change which chromosome the receptor was rearranged on, even
        // when a cross-haplotype (same_haplotype=false) replacement is drawn.
        let mut new_inst =
            AlleleInstance::new(new_id).with_receptor_revision_original_id(preserved_original_id);
        if let Some(h) = haplotype {
            new_inst = new_inst.with_haplotype(h);
        }
        builder.assign_allele(Segment::V, new_inst);
        // 2. Commit the derived 3' trim.
        builder.update_trim(Segment::V, TrimEnd::Three, derived_trim_3);
        // 3. Build the replacement nucleotides — the retained slice
        //    is `allele.seq[0..old_v_len]`. 5' trim stays at 0 in
        //    receptor-revision v1; the same-length constraint
        //    guarantees the slice fits.
        let replacement: Vec<Nucleotide> = (0..old_v_len)
            .map(|i| {
                Nucleotide::germline(new_allele.seq[i as usize], i as u16, Segment::V)
            })
            .collect();
        let _ = builder.replace_segment(Segment::V, replacement);

        if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
            sink.extend(builder.seal_event_log_observer());
        }
        builder.seal()
    }

    /// Strict post-event contract arbitration. Builds the post-
    /// replacement IR, asks the contract bundle whether it admits
    /// the receptor-revised state, and surfaces a structured
    /// `ContractViolation` on rejection. Permissive mode skips this
    /// check — receptor revision is a rare event, and the
    /// post-event filter would just bypass it silently anyway.
    fn validate_post_event_contracts(
        &self,
        pre_sim: &Simulation,
        post_sim: &Simulation,
        ctx: &PassContext,
    ) -> Result<(), PassError> {
        let Some(contracts) = ctx.contracts else {
            return Ok(());
        };
        let context = ChoiceContext::none()
            .with_address(address::ChoiceAddress::ReceptorRevisionApplied);
        contracts
            .admits_post_event(pre_sim, post_sim, ctx.refdata, context)
            .map_err(|violation| PassError::contract_violation(self.name(), vec![violation]))
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // Genotype-aware mode: restrict the replacement V to carried alleles on
        // the drawn rearrangement chromosome (same-haplotype), excluding the
        // current V. The non-genotype path below is unchanged (byte-identical).
        if let Some(constraint) = self.genotype_v.as_ref() {
            return self.execute_genotype_aware(sim, ctx, strict, constraint);
        }
        // Pre-conditions: V must be assigned and refdata must be
        // available. The schedule analyser already orders us after
        // `AssembleSegmentPass(V)`, but a hand-built plan can skip
        // that — surface a structured error rather than panicking
        // inside `refdata.get`.
        if !sim.assignments.has(Segment::V) {
            return Err(PassError::missing_assignment(self.name(), Segment::V));
        }
        let refdata = ctx
            .refdata
            .ok_or_else(|| PassError::missing_refdata(self.name()))?;

        let old_v_len = self.old_v_region_len(sim)?;

        // Replay path consumes the recorded Bool first; the
        // applied/not-applied branch then mirrors the fresh path.
        let replaying = ctx.replay_cursor.is_some();
        let applied = if replaying {
            ctx.replay_cursor
                .as_deref_mut()
                .expect("replay_cursor checked above")
                .expect_bool(address::ChoiceAddress::ReceptorRevisionApplied)
                .map_err(|reason| PassError::replay(self.name(), reason))?
        } else {
            ctx.rng.next_f64() < self.prob
        };
        ctx.trace.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(applied),
        );

        if !applied {
            return Ok(sim.clone());
        }

        // Choose the replacement allele + derive the 3' trim. The
        // replay path additionally validates that the recorded
        // values still satisfy the same-length retained constraint
        // — a refdata swap that shortens the recorded allele below
        // `old_v_len` surfaces as `InvalidDistributionOutput`
        // rather than silently corrupting the V slice.
        let (new_id, derived_trim_3) = if replaying {
            self.consume_replay_choices(refdata, old_v_len, ctx)?
        } else {
            let id = self.sample_replacement_allele(refdata, old_v_len, ctx, strict)?;
            let allele = refdata
                .get(Segment::V, id)
                .ok_or_else(|| PassError::missing_allele(self.name(), Segment::V, id.index()))?;
            // Same-length constraint: permissive fall-through above
            // can hand us a too-short allele. Strict mode would
            // have errored earlier; permissive surfaces as
            // `invalid_plan_state` so the caller sees the
            // contradiction rather than a panic on the trim
            // subtraction below.
            if allele.len() < old_v_len {
                return Err(PassError::invalid_plan_state(
                    self.name(),
                    format!(
                        "sampled V allele length {} < old V region length {}; \
                         no admissible candidate under same-length constraint",
                        allele.len(),
                        old_v_len
                    ),
                ));
            }
            let trim = allele.len() - old_v_len;
            (id, trim)
        };

        // Validate the derived 3' trim fits in `u16` (the persistent
        // trim slot's width). Cannot overflow in practice — alleles
        // longer than `u16::MAX` would have failed validation upstream
        // — but the check is cheap and gives the strict caller a
        // structured error if a future refdata source breaks the
        // invariant.
        if derived_trim_3 > u16::MAX as u32 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::ReceptorRevisionVTrim3.to_string(),
                derived_trim_3 as i64,
                "trim_exceeds_u16",
            ));
        }

        let new_allele = refdata
            .get(Segment::V, new_id)
            .ok_or_else(|| PassError::missing_allele(self.name(), Segment::V, new_id.index()))?
            .clone();

        // Record the typed choices before the IR mutation so the
        // trace order matches the event order downstream consumers
        // (e.g. MCP audit) will see.
        ctx.trace.record_choice(
            address::ChoiceAddress::ReceptorRevisionVAllele,
            ChoiceValue::AlleleId(new_id.index()),
        );
        ctx.trace.record_choice(
            address::ChoiceAddress::ReceptorRevisionVTrim3,
            ChoiceValue::Int(derived_trim_3 as i64),
        );

        // For strict mode, build a hypothetical post-sim and ask
        // the active contract bundle whether the receptor-revised
        // state is admissible. The check runs against the same
        // `commit_replacement` path the real apply uses below, so
        // we don't drift between "what we validated" and "what we
        // committed".
        //
        // The hypothetical commit forwards no event-log sink — the
        // real commit below is what fires events. This keeps the
        // contract check side-effect-free on the trace / sink
        // surfaces.
        if strict && ctx.contracts.is_some() {
            let mut hypothetical_ctx = PassContext {
                trace: ctx.trace,
                rng: ctx.rng,
                pass_index: ctx.pass_index,
                refdata: ctx.refdata,
                contracts: ctx.contracts,
                feasibility: ctx.feasibility,
                reference_index: ctx.reference_index,
                replay_cursor: None,
                event_log_sink: None,
            };
            let hypothetical = self.commit_replacement(
                sim.clone(),
                &new_allele,
                new_id,
                derived_trim_3 as u16,
                old_v_len,
                None,
                &mut hypothetical_ctx,
            );
            self.validate_post_event_contracts(sim, &hypothetical, ctx)?;
        }

        Ok(self.commit_replacement(
            sim.clone(),
            &new_allele,
            new_id,
            derived_trim_3 as u16,
            old_v_len,
            None,
            ctx,
        ))
    }

    /// Consume the (allele_id, trim_3) records the trace carries
    /// when `applied=true`, and validate them against the current
    /// IR state. The validation chain matches the fresh path's
    /// guarantees: allele resolves in refdata, retained length
    /// matches the old V region. A mismatch surfaces as the
    /// structured `PassError` variant that best names *what*
    /// failed.
    fn consume_replay_choices(
        &self,
        refdata: &RefDataConfig,
        old_v_len: u32,
        ctx: &mut PassContext,
    ) -> Result<(AlleleId, u32), PassError> {
        let cursor = ctx
            .replay_cursor
            .as_deref_mut()
            .expect("consume_replay_choices invoked outside replay");
        let id_index = cursor
            .expect_allele_id(address::ChoiceAddress::ReceptorRevisionVAllele)
            .map_err(|reason| PassError::replay(self.name(), reason))?;
        let trim_3_i64 = cursor
            .expect_int(address::ChoiceAddress::ReceptorRevisionVTrim3)
            .map_err(|reason| PassError::replay(self.name(), reason))?;

        let id = AlleleId::new(id_index);
        let allele = refdata
            .get(Segment::V, id)
            .ok_or_else(|| PassError::missing_allele(self.name(), Segment::V, id_index))?;

        if trim_3_i64 < 0 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::ReceptorRevisionVTrim3.to_string(),
                trim_3_i64,
                "negative_trim",
            ));
        }
        if trim_3_i64 > u16::MAX as i64 {
            return Err(PassError::invalid_distribution_output(
                self.name(),
                address::ChoiceAddress::ReceptorRevisionVTrim3.to_string(),
                trim_3_i64,
                "trim_exceeds_u16",
            ));
        }
        let trim_3 = trim_3_i64 as u32;

        // Same-length retained-slice invariant. A recorded
        // `(allele, trim_3)` whose retained length disagrees with
        // the current V region must error rather than silently
        // produce a mismatched-length replacement.
        let retained = (allele.len() as i64).checked_sub(trim_3 as i64).unwrap_or(-1);
        if retained != old_v_len as i64 {
            return Err(PassError::invalid_plan_state(
                self.name(),
                format!(
                    "replay length mismatch: allele {} of length {} with trim_3 {} retains \
                     {} bytes, expected {} (current V region length)",
                    allele.name,
                    allele.len(),
                    trim_3,
                    retained,
                    old_v_len
                ),
            ));
        }

        Ok((id, trim_3))
    }
}

impl Pass for ReceptorRevisionPass {
    fn name(&self) -> &str {
        address::RECEPTOR_REVISION
    }

    fn parameter_signature(&self) -> String {
        // The V replacement distribution is over the same V pool
        // covered by `refdata_content_hash`; skip it (same
        // rationale as `SampleAllelePass`) and pin only `prob`.
        let base = crate::passes::paramsig::fmt_prob("prob", self.prob);
        match &self.genotype_v {
            None => base,
            Some(g) => {
                // Genotype-aware mode pins the constraint state directly so
                // two genotypes with different carried-V sets (or different
                // same_haplotype) produce distinct plan signatures (replay-
                // cache correctness), rather than relying on the genotype
                // pass's own signature.
                use std::fmt::Write;
                let mut s = base;
                let _ = write!(s, "|geno_rr|same={}", g.same_haplotype);
                for (c, list) in g.per_hap.iter().enumerate() {
                    let mut rows: Vec<(u32, u64)> =
                        list.iter().map(|(id, m)| (id.index(), m.to_bits())).collect();
                    rows.sort_unstable();
                    let _ = write!(s, "|h{}={:?}", c, rows);
                }
                s
            }
        }
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("ReceptorRevisionPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
        vec![
            address::ChoiceAddressPattern::ReceptorRevisionApplied,
            address::ChoiceAddressPattern::ReceptorRevisionVAllele,
            address::ChoiceAddressPattern::ReceptorRevisionVTrim3,
        ]
    }

    fn requirements(&self) -> Vec<PassRequirement> {
        // The engine-enforced requirement is `AlleleAssignment(V)` (+ RefData to
        // resolve the replacement allele's bytes); the schedule analyser derives
        // the dependency from this declaration. The stronger "V region already
        // assembled" guarantee is a lowering/schedule-position contract (Python
        // inlines the pass after `push_assemble("J")`), backstopped at runtime by
        // `old_v_region_len()` which errors if no V region is present.
        vec![
            PassRequirement::RefData,
            PassRequirement::AlleleAssignment(Segment::V),
        ]
    }

    fn effects(&self) -> Vec<crate::pass::PassEffect> {
        // Empty by design — receptor revision is a runtime-event-
        // driven pass. The `LiveCallRefreshHook` reacts to the
        // `SegmentReplaced` event (Slice B) without needing a
        // matching compile-effect. Declaring `AssignAllele(V)` /
        // `TrimAllele(V)` / `AppendRegion(V)` here would mislead
        // the schedule analyser into treating receptor revision
        // as initial recombination.
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::contract::ContractSet;
    use crate::dist::AllelePoolDist;
    use crate::ir::{NucHandle, Region, Segment, SimulationEvent};
    use crate::pass::testing::PassRuntime;
    use crate::pass::{PassError, PassPlan};
    use crate::refdata::{Allele, AllelePool, ChainType, RefDataConfig};
    use crate::replay::TraceCursor;
    use crate::rng::Rng;
    use crate::trace::Trace;

    fn allele(name: &str, seq: &[u8]) -> Allele {
        Allele {
            name: name.to_string(),
            gene: name.split('*').next().unwrap_or(name).to_string(),
            seq: seq.to_vec(),
            segment: Segment::V,
            anchor: None,
            functional_status: None,
            subregions: Vec::new(),
        }
    }

    /// Reference data with two V alleles whose retained 6-byte
    /// prefix (under 0 trim) differs. The original assembled V is
    /// `AAAAAA` matching allele V0 exactly; replacement candidate
    /// V1 carries `GGGGGGCC` — length 8, retains 6 bytes at
    /// `trim_3 = 2`.
    fn two_v_refdata() -> (RefDataConfig, AlleleId, AlleleId) {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let v0 = cfg.v_pool.push(allele("V1*01", b"AAAAAA"));
        let v1 = cfg.v_pool.push(allele("V1*02", b"GGGGGGCC"));
        (cfg, v0, v1)
    }

    /// Build a sim with the V slot already assigned to V0 and a
    /// V region of length 6 (`AAAAAA`). Doubles for tests as the
    /// post-recombination starting state.
    fn sim_v_assembled(v0: AlleleId) -> Simulation {
        let mut sim = Simulation::new();
        for (i, &b) in b"AAAAAA".iter().enumerate() {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::germline(
                b,
                i as u16,
                Segment::V,
            ));
            sim = next;
        }
        sim.with_allele_assigned(Segment::V, AlleleInstance::new(v0))
            .with_region_added(Region::new(
                Segment::V,
                NucHandle::new(0),
                NucHandle::new(6),
            ))
    }

    fn v_only_pool(cfg: &RefDataConfig) -> AllelePool {
        cfg.v_pool.clone()
    }

    fn run_pass(
        prob: f64,
        seed: u64,
        cfg: RefDataConfig,
        sim: Simulation,
    ) -> (Trace, Simulation) {
        let mut plan = PassPlan::new();
        let pool = v_only_pool(&cfg);
        plan.push(Box::new(ReceptorRevisionPass::new(
            prob,
            Box::new(AllelePoolDist::uniform(&pool)),
        )));
        let outcome = PassRuntime::execute_with_refdata(&plan, sim, seed, &cfg);
        let final_sim = outcome.final_simulation().clone();
        (outcome.trace, final_sim)
    }

    fn run_with_ctx(
        pass: &ReceptorRevisionPass,
        cfg: &RefDataConfig,
        contracts: Option<&ContractSet>,
        initial: Simulation,
        cursor: Option<&mut TraceCursor>,
        event_sink: Option<&mut Vec<SimulationEvent>>,
    ) -> Result<(Trace, Simulation), PassError> {
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: Some(cfg),
            contracts,
            feasibility: None,
            reference_index: None,
            replay_cursor: cursor,
            event_log_sink: event_sink,
        };
        let next = pass.execute_checked(&initial, &mut ctx)?;
        Ok((trace, next))
    }

    // ── Construction guards ──────────────────────────────────────

    #[test]
    #[should_panic(expected = "prob must be in [0.0, 1.0]")]
    fn receptor_revision_rejects_out_of_range_prob() {
        let (cfg, _, _) = two_v_refdata();
        let _ = ReceptorRevisionPass::new(
            1.5,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        );
    }

    #[test]
    #[should_panic(expected = "prob must be in [0.0, 1.0]")]
    fn receptor_revision_rejects_negative_prob() {
        let (cfg, _, _) = two_v_refdata();
        let _ = ReceptorRevisionPass::new(
            -0.5,
            Box::new(AllelePoolDist::uniform(&cfg.v_pool)),
        );
    }

    // ── genotype-aware candidate selection ──────────────────────

    fn constraint(per_hap: [Vec<(AlleleId, f64)>; 2], same: bool) -> GenotypeVConstraint {
        GenotypeVConstraint { per_hap, same_haplotype: same }
    }

    fn geno_sim(v0: AlleleId) -> Simulation {
        sim_v_assembled(v0)
            .with_allele_assigned(Segment::V, AlleleInstance::new(v0).with_haplotype(0))
    }

    fn geno_pass(cfg: &RefDataConfig, per_hap: [Vec<(AlleleId, f64)>; 2], same: bool) -> ReceptorRevisionPass {
        ReceptorRevisionPass::new(1.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)))
            .with_genotype_constraint(constraint(per_hap, same))
    }

    #[test]
    fn genotype_same_haplotype_excludes_current_and_restricts_to_chromosome() {
        let (cfg, v0, v1) = two_v_refdata();
        // hap0 carries {V0, V1}; hap1 carries {V0}. Current is V0 on hap0.
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        let (trace, after) = run_with_ctx(&pass, &cfg, None, geno_sim(v0), None, None).unwrap();
        assert_eq!(trace.find("receptor_revision.applied").unwrap().value, ChoiceValue::Bool(true));
        // exclude-current: the only eligible alternate on hap0 is V1
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v1);
        assert_eq!(after.assignments.get(Segment::V).unwrap().receptor_revision_original_id, Some(v0));
        assert_eq!(after.assignments.get(Segment::V).unwrap().haplotype, Some(0));
    }

    fn run_permissive(pass: &ReceptorRevisionPass, cfg: &RefDataConfig, initial: Simulation) -> (Trace, Simulation) {
        let mut trace = Trace::new();
        let mut rng = Rng::new(0xc0ff_ee);
        let mut ctx = PassContext {
            trace: &mut trace,
            rng: &mut rng,
            pass_index: 0,
            refdata: Some(cfg),
            contracts: None,
            feasibility: None,
            reference_index: None,
            replay_cursor: None,
            event_log_sink: None,
        };
        let next = pass.execute(&initial, &mut ctx);
        (trace, next)
    }

    #[test]
    fn genotype_no_eligible_alternate_permissive_applied_false() {
        let (cfg, v0, _v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0)], vec![(v0, 1.0)]], true);
        let (trace, after) = run_permissive(&pass, &cfg, geno_sim(v0));
        assert_eq!(trace.find("receptor_revision.applied").unwrap().value, ChoiceValue::Bool(false));
        assert!(trace.find("receptor_revision.v_allele").is_none());
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v0);
    }

    #[test]
    fn genotype_no_eligible_alternate_strict_errors() {
        let (cfg, v0, _v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0)], vec![(v0, 1.0)]], true);
        let err = run_with_ctx(&pass, &cfg, None, geno_sim(v0), None, None).unwrap_err();
        assert!(matches!(err, PassError::ConstraintSampling { .. }), "got {err:?}");
    }

    #[test]
    fn genotype_both_haplotypes_admits_other_chromosome_allele() {
        let (cfg, v0, v1) = two_v_refdata();
        // V1 carried only on hap1; same_haplotype=false aggregates both.
        let pass = geno_pass(&cfg, [vec![(v0, 1.0)], vec![(v1, 1.0)]], false);
        let (_t, after) = run_with_ctx(&pass, &cfg, None, geno_sim(v0), None, None).unwrap();
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v1);
        // haplotype provenance preserved as the original rearrangement chromosome (0)
        assert_eq!(after.assignments.get(Segment::V).unwrap().haplotype, Some(0));
    }

    #[test]
    fn genotype_missing_haplotype_stamp_errors() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        // sim_v_assembled assigns V0 WITHOUT a haplotype stamp.
        let err = run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), None, None).unwrap_err();
        assert!(matches!(err, PassError::InvalidPlanState { .. }), "got {err:?}");
    }

    #[test]
    fn genotype_strict_post_event_contract_rejection_errors() {
        use crate::contract::{Contract, ContractViolation};

        // A contract whose verify() always rejects the post-event state.
        struct RejectAll;
        impl Contract for RejectAll {
            fn name(&self) -> &str {
                "reject_all_test"
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), ContractViolation> {
                Err(ContractViolation::new(self.name(), "rejected by test"))
            }
        }

        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        let contracts = ContractSet::new().with(Box::new(RejectAll));
        // Strict mode (execute_checked) with a rejecting contract must surface a
        // contract violation — parity with the non-genotype path.
        let err =
            run_with_ctx(&pass, &cfg, Some(&contracts), geno_sim(v0), None, None).unwrap_err();
        assert!(matches!(err, PassError::ContractViolation { .. }), "got {err:?}");
    }

    #[test]
    fn genotype_replay_strict_post_event_contract_rejection_errors() {
        use crate::contract::{Contract, ContractViolation};

        struct RejectAll;
        impl Contract for RejectAll {
            fn name(&self) -> &str {
                "reject_all_test"
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), ContractViolation> {
                Err(ContractViolation::new(self.name(), "rejected by test"))
            }
        }

        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        let contracts = ContractSet::new().with(Box::new(RejectAll));
        // A valid applied=true replay (v1, trim 2) must STILL be rejected by the
        // post-event contract in strict mode — replay parity with the fresh path.
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(v1.index()), Some(2)));
        let err = run_with_ctx(&pass, &cfg, Some(&contracts), geno_sim(v0), Some(&mut cursor), None)
            .unwrap_err();
        assert!(matches!(err, PassError::ContractViolation { .. }), "got {err:?}");
    }

    // ── genotype-aware replay validation ────────────────────────

    fn replay_records(applied: bool, allele: Option<u32>, trim: Option<i64>) -> Vec<crate::trace::ChoiceRecord> {
        let mut t = Trace::new();
        t.record_choice(address::ChoiceAddress::ReceptorRevisionApplied, ChoiceValue::Bool(applied));
        if let Some(a) = allele {
            t.record_choice(address::ChoiceAddress::ReceptorRevisionVAllele, ChoiceValue::AlleleId(a));
        }
        if let Some(tr) = trim {
            t.record_choice(address::ChoiceAddress::ReceptorRevisionVTrim3, ChoiceValue::Int(tr));
        }
        t.choices().to_vec()
    }

    #[test]
    fn replay_genotype_valid_replacement_reproduces() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(v1.index()), Some(2)));
        let (_t, after) = run_with_ctx(&pass, &cfg, None, geno_sim(v0), Some(&mut cursor), None).unwrap();
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v1);
        assert!(cursor.is_drained());
    }

    #[test]
    fn replay_genotype_allele_not_carried_on_haplotype_errors() {
        let (cfg, v0, v1) = two_v_refdata();
        // hap0 carries only V0; recorded V1 is not carried on the drawn chromosome.
        let pass = geno_pass(&cfg, [vec![(v0, 1.0)], vec![(v0, 1.0), (v1, 1.0)]], true);
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(v1.index()), Some(2)));
        let err = run_with_ctx(&pass, &cfg, None, geno_sim(v0), Some(&mut cursor), None).unwrap_err();
        assert!(matches!(err, PassError::InvalidDistributionOutput { .. }), "got {err:?}");
    }

    #[test]
    fn replay_genotype_unresolvable_allele_id_reports_missing_allele() {
        let (cfg, v0, _v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0)], vec![(v0, 1.0)]], true);
        // out-of-range recorded v_allele -> missing_allele (not "not_carried"),
        // matching the non-genotype replay diagnostic contract.
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(999), Some(0)));
        let err = run_with_ctx(&pass, &cfg, None, geno_sim(v0), Some(&mut cursor), None).unwrap_err();
        match err {
            PassError::MissingAllele { allele_id, .. } => assert_eq!(allele_id, 999),
            other => panic!("expected MissingAllele, got {other:?}"),
        }
    }

    #[test]
    fn replay_genotype_equals_current_allele_errors() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(v0.index()), Some(0)));
        let err = run_with_ctx(&pass, &cfg, None, geno_sim(v0), Some(&mut cursor), None).unwrap_err();
        assert!(matches!(err, PassError::InvalidDistributionOutput { .. }), "got {err:?}");
    }

    #[test]
    fn replay_genotype_trim_length_mismatch_errors() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = geno_pass(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true);
        // V1 len 8, old 6 -> trim must be 2; record 3 (retains 5) => mismatch.
        let mut cursor = TraceCursor::from_owned(replay_records(true, Some(v1.index()), Some(3)));
        let err = run_with_ctx(&pass, &cfg, None, geno_sim(v0), Some(&mut cursor), None).unwrap_err();
        assert!(matches!(err, PassError::InvalidPlanState { .. }), "got {err:?}");
    }

    #[test]
    fn genotype_signature_differs_by_candidate_set_and_same_haplotype() {
        let (cfg, v0, v1) = two_v_refdata();
        let no_geno = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)))
            .parameter_signature();
        let g_true = geno_pass_prob(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], true).parameter_signature();
        let g_false = geno_pass_prob(&cfg, [vec![(v0, 1.0), (v1, 1.0)], vec![(v0, 1.0)]], false).parameter_signature();
        let g_other = geno_pass_prob(&cfg, [vec![(v0, 2.0)], vec![(v1, 1.0)]], true).parameter_signature();
        assert_ne!(g_true, no_geno);
        assert_ne!(g_true, g_false);
        assert_ne!(g_true, g_other);
    }

    fn geno_pass_prob(cfg: &RefDataConfig, per_hap: [Vec<(AlleleId, f64)>; 2], same: bool) -> ReceptorRevisionPass {
        ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)))
            .with_genotype_constraint(constraint(per_hap, same))
    }

    // ── prob=0: no replacement ──────────────────────────────────

    #[test]
    fn prob_zero_records_applied_false_and_no_mutation_events() {
        let (cfg, v0, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        let mut events: Vec<SimulationEvent> = Vec::new();
        let (trace, after) = run_with_ctx(
            &pass,
            &cfg,
            None,
            sim_v_assembled(v0),
            None,
            Some(&mut events),
        )
        .unwrap();

        // Applied=false recorded.
        let rec = trace
            .find("receptor_revision.applied")
            .expect("applied Bool must be recorded even when prob = 0");
        assert_eq!(rec.value, ChoiceValue::Bool(false));
        // No allele/trim records.
        assert!(trace.find("receptor_revision.v_allele").is_none());
        assert!(trace.find("receptor_revision.v_trim_3").is_none());

        // No state-changing events.
        assert!(events.iter().all(|e| !matches!(
            e,
            SimulationEvent::AssignmentChanged { .. }
                | SimulationEvent::TrimChanged { .. }
                | SimulationEvent::SegmentReplaced { .. }
        )));

        // V assignment + pool bytes unchanged.
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v0);
        let bases: Vec<u8> = after.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"AAAAAA");
    }

    // ── prob=1: replacement ─────────────────────────────────────

    #[test]
    fn prob_one_records_three_choices_and_emits_three_events() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(1.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        let mut events: Vec<SimulationEvent> = Vec::new();
        let (trace, after) = run_with_ctx(
            &pass,
            &cfg,
            None,
            sim_v_assembled(v0),
            None,
            Some(&mut events),
        )
        .unwrap();

        // applied = true
        assert_eq!(
            trace.find("receptor_revision.applied").unwrap().value,
            ChoiceValue::Bool(true),
        );
        // The only length-eligible candidate at old_v_len=6 is V1
        // (length 8). V0 is length 6 too — same length is eligible
        // (trim_3=0), so RNG decides. Either way, the recorded
        // pair must satisfy `len - trim_3 == 6`.
        let recorded_id = match trace.find("receptor_revision.v_allele").unwrap().value {
            ChoiceValue::AlleleId(id) => id,
            _ => panic!("expected AlleleId record"),
        };
        let recorded_trim = match trace.find("receptor_revision.v_trim_3").unwrap().value {
            ChoiceValue::Int(t) => t,
            _ => panic!("expected Int record"),
        };
        let recorded_allele = cfg
            .get(Segment::V, AlleleId::new(recorded_id))
            .expect("recorded allele must resolve");
        assert_eq!(
            recorded_allele.len() as i64 - recorded_trim,
            6,
            "retained length must equal old V region length"
        );
        let _unused = v1;

        // Exactly one of each of the three replacement events.
        let assignments = events
            .iter()
            .filter(|e| matches!(e, SimulationEvent::AssignmentChanged { segment: Segment::V, .. }))
            .count();
        let trims = events
            .iter()
            .filter(|e| matches!(e, SimulationEvent::TrimChanged { segment: Segment::V, end: TrimEnd::Three, .. }))
            .count();
        let replaces = events
            .iter()
            .filter(|e| matches!(e, SimulationEvent::SegmentReplaced { segment: Segment::V, .. }))
            .count();
        assert_eq!(assignments, 1);
        assert_eq!(trims, 1);
        assert_eq!(replaces, 1);

        // The committed V assignment matches the recorded id.
        assert_eq!(
            after.assignments.get(Segment::V).unwrap().allele_id.index(),
            recorded_id,
        );
        // Pool length unchanged (same-length constraint).
        assert_eq!(after.pool.len(), 6);
        // The replacement bytes equal the recorded allele's
        // 6-byte prefix.
        let bases: Vec<u8> = after.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(bases, recorded_allele.seq[..6].to_vec());
    }

    // ── Replay ──────────────────────────────────────────────────

    #[test]
    fn replay_applied_true_reproduces_replacement_without_rng() {
        let (cfg, v0, v1) = two_v_refdata();
        // prob=0.0 would normally fire applied=false; the trace
        // overrides via cursor. Pins "trace is the source of truth".
        let pass = ReceptorRevisionPass::new(0.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(true),
        );
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVAllele,
            ChoiceValue::AlleleId(v1.index()),
        );
        // V1 length 8, old V length 6 → trim_3 = 2.
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVTrim3,
            ChoiceValue::Int(2),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let (trace, after) =
            run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), Some(&mut cursor), None).unwrap();

        assert_eq!(
            trace.find("receptor_revision.applied").unwrap().value,
            ChoiceValue::Bool(true),
        );
        assert_eq!(
            after.assignments.get(Segment::V).unwrap().allele_id,
            v1,
        );
        assert_eq!(
            after.assignments.get(Segment::V).unwrap().trim_3,
            2,
        );
        // V1's retained 6-byte prefix is "GGGGGG".
        let bases: Vec<u8> = after.pool.as_slice().iter().map(|n| n.base).collect();
        assert_eq!(&bases, b"GGGGGG");
        assert!(cursor.is_drained());
    }

    #[test]
    fn replay_applied_false_consumes_only_bool_record() {
        let (cfg, v0, _) = two_v_refdata();
        // prob=1.0 would normally fire applied=true; trace says
        // false → no allele/trim consumption.
        let pass = ReceptorRevisionPass::new(1.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(false),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let (trace, after) =
            run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), Some(&mut cursor), None).unwrap();

        assert_eq!(
            trace.find("receptor_revision.applied").unwrap().value,
            ChoiceValue::Bool(false),
        );
        assert!(trace.find("receptor_revision.v_allele").is_none());
        // Unchanged.
        assert_eq!(after.assignments.get(Segment::V).unwrap().allele_id, v0);
        assert!(cursor.is_drained());
    }

    #[test]
    fn replay_missing_allele_record_after_applied_true_errors() {
        let (cfg, v0, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(true),
        );
        // Missing allele + trim records.
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let err = run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), Some(&mut cursor), None)
            .unwrap_err();
        match err {
            PassError::Replay { pass_name, reason } => {
                assert_eq!(pass_name, "receptor_revision");
                let msg = format!("{reason}");
                assert!(
                    msg.contains("receptor_revision.v_allele") || msg.contains("exhausted"),
                    "expected replay error to mention v_allele or exhaustion, got: {msg}",
                );
            }
            other => panic!("expected PassError::Replay, got {other:?}"),
        }
    }

    #[test]
    fn replay_unknown_allele_id_errors() {
        let (cfg, v0, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(true),
        );
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVAllele,
            ChoiceValue::AlleleId(999),
        );
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVTrim3,
            ChoiceValue::Int(0),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let err = run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), Some(&mut cursor), None)
            .unwrap_err();
        match err {
            PassError::MissingAllele { pass_name, segment, allele_id } => {
                assert_eq!(pass_name, "receptor_revision");
                assert_eq!(segment, Segment::V);
                assert_eq!(allele_id, 999);
            }
            other => panic!("expected PassError::MissingAllele, got {other:?}"),
        }
    }

    #[test]
    fn replay_trim_causing_length_mismatch_errors() {
        let (cfg, v0, v1) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.0, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut input = Trace::new();
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionApplied,
            ChoiceValue::Bool(true),
        );
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVAllele,
            ChoiceValue::AlleleId(v1.index()),
        );
        // V1 has length 8, old V region is 6. A trim_3 of 3 would
        // retain 5 bytes — mismatch.
        input.record_choice(
            address::ChoiceAddress::ReceptorRevisionVTrim3,
            ChoiceValue::Int(3),
        );
        let mut cursor = TraceCursor::from_owned(input.choices().to_vec());

        let err = run_with_ctx(&pass, &cfg, None, sim_v_assembled(v0), Some(&mut cursor), None)
            .unwrap_err();
        match err {
            PassError::InvalidPlanState { pass_name, reason } => {
                assert_eq!(pass_name, "receptor_revision");
                assert!(reason.contains("length mismatch"), "got: {reason}");
            }
            other => panic!("expected PassError::InvalidPlanState, got {other:?}"),
        }
    }

    // ── Plan-state guards ───────────────────────────────────────

    #[test]
    fn missing_v_assignment_errors_in_checked_path() {
        let (cfg, _, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        let err = run_with_ctx(&pass, &cfg, None, Simulation::new(), None, None).unwrap_err();
        match err {
            PassError::MissingAssignment { pass_name, segment } => {
                assert_eq!(pass_name, "receptor_revision");
                assert_eq!(segment, Segment::V);
            }
            other => panic!("expected MissingAssignment, got {other:?}"),
        }
    }

    #[test]
    fn missing_refdata_errors_in_checked_path() {
        // Build a sim with V assigned but execute without refdata.
        let (cfg, v0, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));

        let mut trace = Trace::new();
        let mut rng = Rng::new(0);
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
        let err = pass
            .execute_checked(&sim_v_assembled(v0), &mut ctx)
            .unwrap_err();
        match err {
            PassError::MissingRefData { pass_name } => {
                assert_eq!(pass_name, "receptor_revision");
            }
            other => panic!("expected MissingRefData, got {other:?}"),
        }
    }

    #[test]
    fn no_v_region_errors_in_checked_path() {
        let (cfg, v0, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        // Assignment without a region — assembly never ran.
        let sim = Simulation::new().with_allele_assigned(Segment::V, AlleleInstance::new(v0));
        let err = run_with_ctx(&pass, &cfg, None, sim, None, None).unwrap_err();
        match err {
            PassError::InvalidPlanState { pass_name, reason } => {
                assert_eq!(pass_name, "receptor_revision");
                assert!(reason.contains("no V region"), "got: {reason}");
            }
            other => panic!("expected InvalidPlanState, got {other:?}"),
        }
    }

    // ── Pass metadata ───────────────────────────────────────────

    #[test]
    fn declares_three_choice_patterns_in_order() {
        let (cfg, _, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        assert_eq!(
            pass.declared_choice_patterns(),
            vec![
                address::ChoiceAddressPattern::ReceptorRevisionApplied,
                address::ChoiceAddressPattern::ReceptorRevisionVAllele,
                address::ChoiceAddressPattern::ReceptorRevisionVTrim3,
            ],
        );
    }

    #[test]
    fn declares_refdata_and_v_assignment_requirements() {
        let (cfg, _, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        assert_eq!(
            pass.requirements(),
            vec![
                PassRequirement::RefData,
                PassRequirement::AlleleAssignment(Segment::V),
            ],
        );
    }

    #[test]
    fn declares_no_compile_effects() {
        // Pinning the empty declaration — receptor revision is
        // event-driven, the schedule analyser must not treat it
        // as initial recombination.
        let (cfg, _, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        assert!(pass.effects().is_empty());
    }

    #[test]
    fn pass_name_is_stable() {
        // The frozen pass name flows through `pass_plan_signature`
        // — a rename here breaks every existing trace.
        let (cfg, _, _) = two_v_refdata();
        let pass = ReceptorRevisionPass::new(0.5, Box::new(AllelePoolDist::uniform(&cfg.v_pool)));
        assert_eq!(pass.name(), "receptor_revision");
    }

    // ── Full plan round-trip (V assembled → replacement) ────────

    #[test]
    fn full_runtime_round_trip_with_prob_one() {
        let (cfg, v0, _) = two_v_refdata();
        let (trace, after) = run_pass(1.0, 0, cfg.clone(), sim_v_assembled(v0));

        // Three records present.
        assert!(trace.find("receptor_revision.applied").is_some());
        assert!(trace.find("receptor_revision.v_allele").is_some());
        assert!(trace.find("receptor_revision.v_trim_3").is_some());

        // V region length unchanged (same-length constraint).
        let v_region = after
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)
            .expect("V region must remain after replacement");
        assert_eq!(v_region.len(), 6);
        assert_eq!(after.pool.len(), 6);
    }
}
