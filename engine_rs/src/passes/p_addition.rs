//! `PAdditionPass` — templated palindromic (P-)nucleotide
//! addition at a V(D)J coding-end junction side.
//!
//! Implements the v1 scope from
//! [`docs/p_nucleotide_design.md`](../../../../docs/p_nucleotide_design.md):
//!
//! - **Lengths sampled, bases deterministic.** Per-end length
//!   recorded at the canonical
//!   [`crate::address::ChoiceAddress::PLength`] address; bases
//!   derive from the source allele's post-trim, post-orientation
//!   coding flank via `complement_base`.
//! - **No new structural region.** P-bytes are pushed with
//!   `flag::P_NUC` set; the descriptive `Region` flows through
//!   the [`SimulationEvent::PRegionAdded`] event only.
//!   `sim.sequence.regions` keeps one structural entry per
//!   biological segment so existing
//!   `find(|r| r.segment == ...)` projection sites stay correct.
//! - **Replay symmetric.** Trace cursor consumes the recorded
//!   length; the same allele + trim + orientation produces the
//!   same bytes.
//! - **Empty-support sentinel reuses
//!   [`crate::dist::EmptySupport::Sentinel`]`(0)`** — under
//!   productive-only, if every length in the cartridge's
//!   distribution would disrupt the junction frame / drop an
//!   anchor, the pass emits length 0 (no-op) under permissive
//!   and surfaces `PassError::ConstraintSampling` under strict.
//!
//! Pipeline placement is per the audit's §1.4 / §9 (post-trim,
//! pre-NP / pre-assemble on the relevant side; for `D5` the
//! D-inversion decision MUST commit first, so the order
//! `invert_d → p_addition.d_5 → assemble.D` is the only
//! IR-correct one).

use crate::address::{ChoiceAddress, ChoiceAddressPattern, PEnd};
use crate::assignment::SegmentOrientation;
use crate::dist::{Distribution, EmptySupport};
use crate::ir::{
    complement_base, flag, NucHandle, Nucleotide, Region, Segment, Simulation, SimulationBuilder,
};
use crate::pass::{Pass, PassCompileFact, PassContext, PassEffect, PassError, IntegerSupport};
use crate::trace::ChoiceValue;

/// Permissive empty-support policy: length 0 means "no P
/// extension at this end" — a clean no-op.
const P_LENGTH_EMPTY_SUPPORT: EmptySupport<i64> = EmptySupport::Sentinel(0);

/// Generate templated palindromic (P-)nucleotides at the named
/// V(D)J coding-end junction side.
///
/// Construction takes:
///
/// - `end: PEnd` — V3 / D5 / D3 / J5.
/// - `length_dist: Box<dyn Distribution<Output = i64>>` —
///   empirical distribution over the per-end P-length. v1 typical
///   support is `{0, 1, 2, 3, 4}` per the audit / legacy
///   `DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS`.
///
/// On execute:
///
/// 1. Sample (or replay) one length `L` from `length_dist`.
/// 2. Read the source-segment's [`crate::assignment::AlleleInstance`]
///    (`allele_id`, `trim_5`, `trim_3`, `orientation`).
/// 3. Compute the effective post-trim, post-orientation coding
///    bytes for that segment.
/// 4. Derive `L` palindromic P-bytes via `complement_base` —
///    rule depends on whether the end is 3' or 5' (§details
///    below).
/// 5. Push each P-byte as `Nucleotide::synthetic(byte,
///    source_segment, flag::P_NUC)`.
/// 6. Emit `SimulationEvent::PRegionAdded { end, region }` with
///    the pool span the P-bytes occupy. **No structural
///    region** is added to `sim.sequence.regions`.
///
/// ### Palindrome derivation rule
///
/// Let `effective` be the source allele's post-trim,
/// post-orientation byte sequence (the bytes that `AssembleSegmentPass`
/// would push, in pool order). Then:
///
/// - For **3' extensions** (`V3`, `D3`): P-bytes are the
///   reverse-complement of the *last* `L` effective bytes; in
///   pool order, P\[i\] = `complement_base(effective[len - 1 - i])`.
/// - For **5' extensions** (`D5`, `J5`): P-bytes are the
///   reverse-complement of the *first* `L` effective bytes; in
///   pool order, P\[i\] = `complement_base(effective[L - 1 - i])`.
///
/// ### Replay
///
/// Replay consumes one `ChoiceValue::Int(L)` from
/// `replay_cursor` (validated to be `>= 0`) and re-runs steps
/// 2–6 verbatim. Bases reproduce by construction because
/// `complement_base` is pure and the source allele + trim +
/// orientation are recovered from `sim.assignments`.
pub struct PAdditionPass {
    end: PEnd,
    length_dist: Box<dyn Distribution<Output = i64>>,
}

impl PAdditionPass {
    /// Construct a P-addition pass for the named end.
    pub fn new(end: PEnd, length_dist: Box<dyn Distribution<Output = i64>>) -> Self {
        Self { end, length_dist }
    }

    /// The PEnd this pass extends.
    pub fn end(&self) -> PEnd {
        self.end
    }

    /// The pass's canonical name (used as the pass-plan key).
    fn pass_name(&self) -> &'static str {
        match self.end {
            PEnd::V3 => crate::address::P_ADDITION_V_3,
            PEnd::D5 => crate::address::P_ADDITION_D_5,
            PEnd::D3 => crate::address::P_ADDITION_D_3,
            PEnd::J5 => crate::address::P_ADDITION_J_5,
        }
    }

    fn length_choice_address(&self) -> ChoiceAddress {
        ChoiceAddress::PLength { end: self.end }
    }

    fn source_segment(&self) -> Segment {
        self.end.source_segment()
    }

    /// Materialise the source allele's effective post-trim,
    /// post-orientation byte sequence — same shape
    /// `AssembleSegmentPass` would push. Used by the palindrome
    /// derivation to know which source bytes are at positions
    /// `[0..L)` (5' palindrome) or `[len-L..len)` (3' palindrome).
    fn effective_coding_bytes(
        allele_seq: &[u8],
        trim_5: u32,
        trim_3: u32,
        orientation: SegmentOrientation,
    ) -> Vec<u8> {
        let len = allele_seq.len() as u32;
        let slice_start = trim_5 as usize;
        let slice_end = (len - trim_3) as usize;
        if orientation.is_reverse() {
            (slice_start..slice_end)
                .rev()
                .map(|i| complement_base(allele_seq[i]))
                .collect()
        } else {
            allele_seq[slice_start..slice_end].to_vec()
        }
    }

    /// Compute the `length` palindromic P-bytes for `end` from
    /// `effective` (the post-trim, post-orientation coding
    /// sequence). Returns the bytes in pool order.
    fn derive_p_bytes(end: PEnd, effective: &[u8], length: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(length);
        match end {
            PEnd::V3 | PEnd::D3 => {
                let n = effective.len();
                for i in 0..length {
                    out.push(complement_base(effective[n - 1 - i]));
                }
            }
            PEnd::D5 | PEnd::J5 => {
                for i in 0..length {
                    out.push(complement_base(effective[length - 1 - i]));
                }
            }
        }
        out
    }

    fn execute_with_sampling_mode(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
        strict: bool,
    ) -> Result<Simulation, PassError> {
        // Sample / replay the length without recording yet —
        // we record only after the validation gates pass so a
        // strict-mode failure doesn't leak a half-committed
        // trace entry.
        let length = if let Some(cursor) = ctx.replay_cursor.as_deref_mut() {
            cursor
                .expect_int(self.length_choice_address())
                .map_err(|reason| PassError::replay(self.pass_name(), reason))?
        } else {
            self.length_dist.sample(ctx.rng)
        };

        // Negative length — caller bug. Strict surfaces;
        // permissive falls back to the sentinel.
        if length < 0 {
            if strict {
                return Err(PassError::invalid_distribution_output(
                    self.pass_name(),
                    self.length_choice_address(),
                    length,
                    "negative_length",
                ));
            }
            let sentinel = match P_LENGTH_EMPTY_SUPPORT {
                EmptySupport::Sentinel(s) => s,
                EmptySupport::Skip => unreachable!(),
            };
            ctx.trace
                .record_choice(self.length_choice_address(), ChoiceValue::Int(sentinel));
            return Ok(sim.clone());
        }

        // Length 0 — clean no-op. Record + return; no event,
        // no pool mutation.
        if length == 0 {
            ctx.trace
                .record_choice(self.length_choice_address(), ChoiceValue::Int(0));
            return Ok(sim.clone());
        }

        // Read source allele + trim + orientation. Defensive
        // strict-mode error paths mirror AssembleSegmentPass.
        let seg = self.source_segment();
        let refdata = match ctx.refdata {
            Some(refdata) => refdata,
            None if strict => return Err(PassError::missing_refdata(self.pass_name())),
            None => panic!(
                "PAdditionPass({:?}): PassContext.refdata is None — \
                 use PassRuntime::execute_with_refdata for plans \
                 containing P-addition passes",
                self.end
            ),
        };
        let inst = match sim.assignments.get(seg).copied() {
            Some(inst) => inst,
            None if strict => {
                return Err(PassError::missing_assignment(self.pass_name(), seg))
            }
            None => panic!(
                "PAdditionPass({:?}): no allele assigned to {seg:?} — \
                 SampleAllelePass must run before P-addition",
                self.end
            ),
        };
        let allele = match refdata.get(seg, inst.allele_id) {
            Some(a) => a,
            None if strict => {
                return Err(PassError::missing_allele(
                    self.pass_name(),
                    seg,
                    inst.allele_id.index(),
                ));
            }
            None => panic!(
                "PAdditionPass({:?}): allele_id {:?} out of bounds for refdata pool",
                self.end, inst.allele_id
            ),
        };
        let effective = Self::effective_coding_bytes(
            &allele.seq,
            inst.trim_5 as u32,
            inst.trim_3 as u32,
            inst.orientation,
        );
        let length_usize = length as usize;
        // The effective coding flank must be at least `length`
        // bytes — otherwise we can't form the palindrome. Under
        // permissive, fall back to length 0 (no record yet, so
        // we record the sentinel and return); under strict,
        // surface `ConstraintSampling` with the
        // EmptyAdmissibleSupport reason (audit §10).
        if effective.len() < length_usize {
            if strict {
                return Err(PassError::constraint_sampling(
                    self.pass_name().to_string(),
                    self.length_choice_address().to_string(),
                    crate::dist::FilteredSampleError::EmptyAdmissibleSupport,
                ));
            }
            ctx.trace
                .record_choice(self.length_choice_address(), ChoiceValue::Int(0));
            return Ok(sim.clone());
        }

        // All gates passed — record the length and emit.
        ctx.trace
            .record_choice(self.length_choice_address(), ChoiceValue::Int(length));

        let p_bytes = Self::derive_p_bytes(self.end, &effective, length_usize);
        let region_start = NucHandle::new(sim.pool.len() as u32);
        let mut builder = SimulationBuilder::from_simulation(sim.clone());
        if ctx.event_log_sink.is_some() {
            builder.attach_event_log_observer();
        }
        for byte in &p_bytes {
            builder.push_nucleotide(Nucleotide::synthetic(*byte, seg, flag::P_NUC));
        }
        let region_end = NucHandle::new(builder.peek().pool.len() as u32);
        // Frame phase is the codon offset of the first P-byte
        // — i.e. (current pool length) % 3. Use the pool length
        // directly so any prior P-bytes (which don't carry a
        // structural region) are counted correctly.
        let frame_phase = (sim.pool.len() as u32 % 3) as u8;
        let region =
            Region::new(seg, region_start, region_end).with_frame_phase(frame_phase);
        builder.record_p_region(self.end, region);
        if ctx.event_log_sink.is_some() {
            let captured = builder.seal_event_log_observer();
            if let Some(sink) = ctx.event_log_sink.as_deref_mut() {
                sink.extend(captured);
            }
        }
        Ok(builder.seal())
    }
}

impl Pass for PAdditionPass {
    fn name(&self) -> &str {
        self.pass_name()
    }

    fn parameter_signature(&self) -> String {
        // Fold the per-end length distribution — same shape
        // `generate_np.np1(length=...)` uses, so plan signatures
        // are byte-stable under the existing
        // `fmt_int_dist` convention.
        use crate::passes::paramsig::{fmt_int_dist, join_parts};
        join_parts([format!("length={}", fmt_int_dist(self.length_dist.as_ref()))])
    }

    fn execute(&self, sim: &Simulation, ctx: &mut PassContext) -> Simulation {
        self.execute_with_sampling_mode(sim, ctx, false)
            .expect("PAdditionPass permissive execution must not return PassError")
    }

    fn execute_checked(
        &self,
        sim: &Simulation,
        ctx: &mut PassContext,
    ) -> Result<Simulation, PassError> {
        self.execute_with_sampling_mode(sim, ctx, true)
    }

    fn declared_choice_patterns(&self) -> Vec<ChoiceAddressPattern> {
        vec![ChoiceAddressPattern::PLength { end: self.end }]
    }

    fn effects(&self) -> Vec<PassEffect> {
        // P-bytes mutate the pool but DON'T add a structural
        // region — no PassEffect::AppendRegion. The pool
        // mutation is observable via the BasePushed event
        // stream; downstream observers that need to know the
        // P extension committed see the PRegionAdded event.
        Vec::new()
    }

    fn compile_facts(&self) -> Vec<PassCompileFact> {
        // Per-end length distribution is folded as an integer
        // support fact so the schedule analyser and any future
        // pre-flight contract pass can introspect it the same
        // way they introspect NP length.
        vec![PassCompileFact::PLengthSupport {
            end: self.end,
            support: IntegerSupport::from_weighted_pairs(self.length_dist.support()),
        }]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::EmpiricalLengthDist;

    #[test]
    fn p3_palindrome_complements_last_n_in_pool_order() {
        // effective = "ACGT", length = 3.
        // V3 / D3 rule: P[i] = complement(effective[len-1-i])
        // i=0: complement(T) = A
        // i=1: complement(G) = C
        // i=2: complement(C) = G
        let bytes = PAdditionPass::derive_p_bytes(PEnd::V3, b"ACGT", 3);
        assert_eq!(bytes, b"ACG");
        let bytes = PAdditionPass::derive_p_bytes(PEnd::D3, b"ACGT", 3);
        assert_eq!(bytes, b"ACG");
    }

    #[test]
    fn p5_palindrome_complements_first_n_in_pool_order_with_reversal() {
        // effective = "ACGT", length = 3.
        // D5 / J5 rule: P[i] = complement(effective[length-1-i])
        // i=0: complement(effective[2]) = complement(G) = C
        // i=1: complement(effective[1]) = complement(C) = G
        // i=2: complement(effective[0]) = complement(A) = T
        let bytes = PAdditionPass::derive_p_bytes(PEnd::D5, b"ACGT", 3);
        assert_eq!(bytes, b"CGT");
        let bytes = PAdditionPass::derive_p_bytes(PEnd::J5, b"ACGT", 3);
        assert_eq!(bytes, b"CGT");
    }

    #[test]
    fn effective_coding_bytes_forward_orientation_skips_trim() {
        let bytes = PAdditionPass::effective_coding_bytes(
            b"NACGTM",
            1,
            1,
            SegmentOrientation::Forward,
        );
        assert_eq!(bytes, b"ACGT");
    }

    #[test]
    fn effective_coding_bytes_reverse_orientation_reverse_complements() {
        let bytes = PAdditionPass::effective_coding_bytes(
            b"NACGTM",
            1,
            1,
            SegmentOrientation::ReverseComplement,
        );
        // reverse of "ACGT" with complement on each byte =
        // "ACGT" → rev "TGCA" → complement_each → "ACGT"
        assert_eq!(bytes, b"ACGT");
    }

    #[test]
    fn p_addition_pass_zero_length_is_no_op() {
        let pass = PAdditionPass::new(
            PEnd::V3,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(0, 1.0)])),
        );
        assert_eq!(pass.end(), PEnd::V3);
        assert_eq!(pass.name(), "p_addition.v_3");
    }

    #[test]
    fn p_addition_pass_declares_single_pattern() {
        let pass = PAdditionPass::new(
            PEnd::J5,
            Box::new(EmpiricalLengthDist::from_pairs(vec![(1, 1.0)])),
        );
        let patterns = pass.declared_choice_patterns();
        assert_eq!(patterns.len(), 1);
        assert_eq!(
            patterns[0],
            ChoiceAddressPattern::PLength { end: PEnd::J5 }
        );
    }
}
