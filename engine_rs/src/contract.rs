//! Contracts — composable predicates over the simulation IR.
//!
//! ## What a contract is
//!
//! Per design doc §6 + D6 + D7, a contract is a first-class
//! predicate that asserts something about the simulation. Contracts
//! have two architectural modes:
//!
//! 1. **Verify** — given a simulation IR (and optionally the
//!    reference data), decide whether the contract holds. Returns
//!    `Ok(())` if satisfied, `Err(ContractViolation)` with a
//!    structured reason if not. Used at any point where invariants
//!    need to be checked: build-time validation (D7 Phase 1),
//!    debug introspection, post-pipeline assertions.
//!
//! 2. **Filter** — given a candidate sampling action and the
//!    current state, decide whether the action is admissible
//!    *before* sampling. This is what makes constraint-aware
//!    sampling work (D6 — `respect=[productive()]`). Filter mode
//!    arrives in Phase D when contracts get wired into the
//!    `Distribution` trait surface.
//!
//! ## Phase C.9 scope
//!
//! Verify mode only. Filter mode is deliberately deferred to Phase
//! D where it has somewhere to plug into. The `Contract` trait
//! defined here is forwards-compatible: filter methods can be
//! added as defaulted trait methods later without breaking
//! existing contract implementations.
//!
//! ## First concrete contract — `AnchorPreserved`
//!
//! Verifies that an assigned allele's anchor codon (the conserved
//! Cys for V or W/F for J) remains within the post-trim retained
//! slice. Trims that erase the anchor produce a `ContractViolation`.

use crate::ir::{translate_codon, NucHandle, Segment, Simulation, AMINO_STOP};
use crate::junction::compute_junction;
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

// ──────────────────────────────────────────────────────────────────
// ContractViolation — structured diagnostic for a failed verify
// ──────────────────────────────────────────────────────────────────

/// One reason a contract verification failed.
///
/// Carries a stable contract identifier (`contract_name`) for
/// programmatic dispatch and a human-readable `reason` string for
/// diagnostics. Phase D's strict-mode failure shape (D7) builds
/// on this.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ContractViolation {
    pub contract_name: String,
    pub reason: String,
}

impl ContractViolation {
    pub fn new(contract_name: impl Into<String>, reason: impl Into<String>) -> Self {
        Self {
            contract_name: contract_name.into(),
            reason: reason.into(),
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// ChoiceContext — optional execution context for candidate filtering
// ──────────────────────────────────────────────────────────────────

/// Extra context for a candidate choice being filtered by contracts.
///
/// Plain addressed choices only carry `"np.np1.bases[3]"` plus the
/// candidate value. Some contracts need bounded local execution
/// context to remain precise:
/// - NP-base filters need the current draw index and total planned
///   draw count to distinguish known future fixed bases from still
///   random future bases.
/// - Site-based transforms need the target nucleotide handle so a
///   contract can evaluate the exact post-candidate local state.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct ChoiceContext {
    pub draw_index: Option<u32>,
    pub draw_count: Option<u32>,
    pub target: Option<NucHandle>,
}

impl ChoiceContext {
    pub const fn none() -> Self {
        Self {
            draw_index: None,
            draw_count: None,
            target: None,
        }
    }

    pub const fn indexed(draw_index: u32, draw_count: u32) -> Self {
        Self {
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: None,
        }
    }

    pub const fn indexed_target(draw_index: u32, draw_count: u32, target: NucHandle) -> Self {
        Self {
            draw_index: Some(draw_index),
            draw_count: Some(draw_count),
            target: Some(target),
        }
    }
}

impl Default for ChoiceContext {
    fn default() -> Self {
        Self::none()
    }
}

// ──────────────────────────────────────────────────────────────────
// Contract trait
// ──────────────────────────────────────────────────────────────────

/// A predicate over the simulation IR.
///
/// Phase C.9 surfaces only `verify`. Phase D will add filter
/// methods (and possibly an `upstream_bound` for backward
/// constraint propagation) as defaulted trait methods so existing
/// contract implementations continue to compile.
pub trait Contract {
    /// Stable, human-readable identifier for this contract.
    /// Convention: hierarchical-string addresses matching D3 (e.g.,
    /// `"anchor_preserved.v"`, `"productive_junction_frame"`).
    fn name(&self) -> &str;

    /// Verify mode: does this contract hold for `sim`?
    ///
    /// `refdata` is optional because some contracts don't need it
    /// (e.g., a contract that only inspects the assembled bases).
    /// Contracts that *do* need it should treat `None` as
    /// "insufficient data — skip" by returning `Ok(())`, not by
    /// panicking. This keeps verification safe to call in any
    /// context.
    fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation>;

    /// Filter mode: would `candidate` at `address` keep this
    /// contract satisfiable given the current `sim` state?
    ///
    /// Used by constraint-aware sampling (D.6 onward): before a
    /// sampling pass commits to a draw, it asks every active
    /// contract `admits(sim, refdata, addr, candidate)` and only
    /// accepts candidates that all contracts admit.
    ///
    /// **Default behaviour**: `Ok(())` — "always allow." Contracts
    /// that can usefully prune candidates at sampling time
    /// override this method. Contracts that can only check after
    /// a transform applies (e.g., `NoStopCodonInJunction` looking
    /// at codons that don't exist yet) keep the default; their
    /// constraints get enforced via `verify` post-hoc, or by
    /// future contract-aware base-sampling passes (Phase E).
    ///
    /// **Returning `Err` is not a fatal failure** — it's "this
    /// specific candidate is inadmissible." The caller (the
    /// sampling pass via the `ContractSet` from D.5) treats it as
    /// a filter signal and tries another candidate.
    ///
    /// `address` follows the D3 hierarchical-string convention:
    /// `"trim.v_3"`, `"np.np1.length"`, `"sample_allele.v"`, etc.
    /// Contracts dispatch on prefix to handle the addresses they
    /// care about and ignore the rest.
    fn admits(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        // Defaulted: ignore inputs, always allow.
        let _ = (sim, refdata, address, candidate);
        Ok(())
    }

    /// Context-aware filter mode. Defaults to the simpler `admits`
    /// implementation so existing contracts only override this when
    /// they need execution-local metadata such as indexed draw count.
    fn admits_with_context(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        let _ = context;
        self.admits(sim, refdata, address, candidate)
    }
}

// ──────────────────────────────────────────────────────────────────
// AnchorPreserved — anchor codon must remain in the retained slice
// ──────────────────────────────────────────────────────────────────

/// Verifies that the assigned allele's anchor codon (3 bases
/// starting at `allele.anchor`) sits entirely within the post-trim
/// retained slice for a given segment.
///
/// Vacuously satisfied (returns `Ok`) when:
/// - no allele is assigned to `segment` yet,
/// - no `RefDataConfig` is provided to look the allele up in,
/// - the allele has no anchor (`Allele::anchor == None` —
///   pseudogenes / partial alleles).
///
/// Returns a `ContractViolation` when:
/// - `trim_5 > allele.anchor` (anchor 5'-trimmed away), or
/// - `allele.anchor + 3 > allele.len() - trim_3` (anchor
///   3'-trimmed away).
pub struct AnchorPreserved {
    segment: Segment,
}

impl AnchorPreserved {
    /// Construct the contract for a given segment.
    /// Panics if `segment` is not V, D, or J.
    pub fn new(segment: Segment) -> Self {
        match segment {
            Segment::V | Segment::D | Segment::J => {}
            _ => panic!(
                "AnchorPreserved: segment must be V, D, or J — got {:?}",
                segment
            ),
        }
        Self { segment }
    }

    pub fn segment(&self) -> Segment {
        self.segment
    }

    fn name_for(segment: Segment) -> &'static str {
        match segment {
            Segment::V => "anchor_preserved.v",
            Segment::D => "anchor_preserved.d",
            Segment::J => "anchor_preserved.j",
            _ => unreachable!("AnchorPreserved with non-V/D/J segment"),
        }
    }
}

impl Contract for AnchorPreserved {
    fn name(&self) -> &str {
        Self::name_for(self.segment)
    }

    fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation> {
        // No allele yet → vacuously satisfied.
        let inst = match sim.assignments.get(self.segment) {
            None => return Ok(()),
            Some(i) => i,
        };

        // No refdata → can't verify; treat as satisfied.
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };

        // Allele not in pool → can't verify; treat as satisfied.
        let allele = match refdata.get(self.segment, inst.allele_id) {
            None => return Ok(()),
            Some(a) => a,
        };

        // Anchorless allele → can't violate.
        let anchor = match allele.anchor {
            None => return Ok(()),
            Some(a) => a as u32,
        };

        let trim_5 = inst.trim_5 as u32;
        let trim_3 = inst.trim_3 as u32;
        let allele_len = allele.len();
        let retained_end = allele_len.saturating_sub(trim_3);

        // Anchor codon spans [anchor, anchor + 3). Must be fully
        // within the retained slice [trim_5, retained_end).
        if anchor < trim_5 {
            return Err(ContractViolation::new(
                self.name(),
                format!(
                    "anchor at allele position {} but trim_5 = {} \
                     (anchor 5'-trimmed away from {} '{}')",
                    anchor, trim_5, allele.name, allele.gene
                ),
            ));
        }
        if anchor + 3 > retained_end {
            return Err(ContractViolation::new(
                self.name(),
                format!(
                    "anchor codon at allele positions [{}, {}) but retained \
                     slice ends at {} (anchor 3'-trimmed from {} '{}')",
                    anchor,
                    anchor + 3,
                    retained_end,
                    allele.name,
                    allele.gene
                ),
            ));
        }

        Ok(())
    }
}

// ──────────────────────────────────────────────────────────────────
// ProductiveJunctionFrame — junction length must be divisible by 3
// ──────────────────────────────────────────────────────────────────

/// Verifies that the junction (V Cys → J W/F + 3) has a length
/// divisible by 3 — i.e., the codon frame closes cleanly across
/// the V/NP/D/NP/J junction.
///
/// Vacuously satisfied (returns `Ok`) when:
/// - no `RefDataConfig` is provided to look anchors up in,
/// - the junction is undefined (no V or J assignment, anchorless
///   allele, V or J region not yet assembled, etc. — see
///   `compute_junction`).
///
/// Returns a `ContractViolation` when the junction is materialized
/// and its length is not divisible by 3.
///
/// **Note:** this contract checks frame *only*. Stop codons inside
/// an in-frame junction are the `NoStopCodonInJunction` contract's
/// concern (D.3). The `productive()` bundle (D.5) composes both.
pub struct ProductiveJunctionFrame;

impl ProductiveJunctionFrame {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ProductiveJunctionFrame {
    fn default() -> Self {
        Self::new()
    }
}

impl Contract for ProductiveJunctionFrame {
    fn name(&self) -> &str {
        "productive_junction_frame"
    }

    fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation> {
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };

        let junction = match compute_junction(sim, refdata) {
            None => return Ok(()),
            Some(j) => j,
        };

        if junction.is_in_frame() {
            return Ok(());
        }

        Err(ContractViolation::new(
            self.name(),
            format!(
                "junction length {} is not divisible by 3 (start={}, end={})",
                junction.length,
                junction.start.index(),
                junction.end.index()
            ),
        ))
    }

    fn admits(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        // Filter NP length samples to values that produce an in-frame
        // junction. The architectural pivot from D.4 → D.6: instead
        // of verifying after-the-fact, prune the candidate distribution
        // before sampling.
        //
        // We can compute the hypothetical junction length when:
        // - the NP being sampled is the LAST event before J assembly,
        // - we have V assembled (so V's pool position is known),
        // - we have V and J anchors known.
        //
        // For VJ chains: NP1 is the last event before J → filter np.np1.length.
        // For VDJ chains: NP2 is the last event before J → filter np.np2.length.
        // VDJ NP1 has D + NP2 + J between it and the junction end — NP2 will
        // compensate, so we don't filter NP1 in the VDJ case.
        //
        // **Pre-condition (asserted in debug builds):** the J region
        // must NOT be assembled yet. The hypothetical-J-start math
        // below assumes `sim.pool.len() + length` is where J will
        // start; that's only true if J hasn't already been pushed
        // into the pool. Phase A-D plans honor this by always
        // running NP generation before J assembly. A future plan
        // (e.g., a Phase E pass that pushes contaminant or adapter
        // bases into the pool before J assembly) that violates this
        // invariant would silently produce wrong frame predictions.
        // The debug assertion catches that drift in dev builds; in
        // release builds the contract may produce an over-eager
        // rejection (false-negative admits) but never a hidden
        // soundness violation.
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };

        let length = match candidate {
            ChoiceValue::Int(n) if *n >= 0 => *n as u32,
            _ => return Ok(()),
        };

        let is_vj = sim.assignments.d.is_none();
        let applicable = match address {
            "np.np1.length" => is_vj,
            "np.np2.length" => !is_vj,
            _ => false,
        };
        if !applicable {
            return Ok(());
        }

        // Need V/J alleles + anchors + V region in pool.
        let v_inst = match sim.assignments.v {
            None => return Ok(()),
            Some(v) => v,
        };
        let j_inst = match sim.assignments.j {
            None => return Ok(()),
            Some(j) => j,
        };
        let v_allele = match refdata.get(Segment::V, v_inst.allele_id) {
            None => return Ok(()),
            Some(a) => a,
        };
        let j_allele = match refdata.get(Segment::J, j_inst.allele_id) {
            None => return Ok(()),
            Some(a) => a,
        };
        let v_anchor = match v_allele.anchor {
            None => return Ok(()),
            Some(a) => a as u32,
        };
        let j_anchor = match j_allele.anchor {
            None => return Ok(()),
            Some(a) => a as u32,
        };

        let v_trim_5 = v_inst.trim_5 as u32;
        let j_trim_5 = j_inst.trim_5 as u32;
        if v_trim_5 > v_anchor || j_trim_5 > j_anchor {
            return Ok(());
        }

        let v_region = match sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)
        {
            None => return Ok(()),
            Some(r) => r,
        };
        let v_anchor_pool = v_region.start.index() + (v_anchor - v_trim_5);

        // J must not yet be assembled — see method-level doc above.
        debug_assert!(
            sim.sequence.regions.iter().all(|r| r.segment != Segment::J),
            "ProductiveJunctionFrame::admits at {}: J region is already \
             assembled in sim.sequence.regions; the hypothetical-J-start \
             math below assumes J has not been added to the pool yet. \
             This indicates a pass-ordering bug — NP generation must run \
             before J assembly.",
            address
        );

        // Hypothetical J start: pool grows by `length` more bases
        // when the NP gets generated, then J starts immediately
        // after.
        let hypothetical_j_start = sim.pool.len() as u32 + length;
        let hypothetical_j_anchor_pool = hypothetical_j_start + (j_anchor - j_trim_5);
        // Defensive: junction must be a positive window.
        if hypothetical_j_anchor_pool + 3 <= v_anchor_pool {
            return Ok(());
        }
        let junction_length = (hypothetical_j_anchor_pool + 3) - v_anchor_pool;

        if junction_length % 3 == 0 {
            Ok(())
        } else {
            Err(ContractViolation::new(
                self.name(),
                format!(
                    "NP length {} at {} would produce out-of-frame \
                     junction (hypothetical length {})",
                    length, address, junction_length
                ),
            ))
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// NoStopCodonInJunction — junction codons must translate to amino acids
// ──────────────────────────────────────────────────────────────────

/// Verifies that no codon inside the junction translates to a stop
/// (TAA, TAG, TGA). Walks codons from `junction.start` in steps of
/// 3 and checks each translation.
///
/// Vacuously satisfied (returns `Ok`) when:
/// - no `RefDataConfig` is provided,
/// - the junction is undefined (`compute_junction` returns `None`),
/// - the junction is out of frame — that's a different contract's
///   concern (`ProductiveJunctionFrame`); we pass on rather than
///   double-reporting.
///
/// Returns a `ContractViolation` with the first stop's position
/// and codon string when a stop is detected.
///
/// **Note:** the contract reads from `sim.pool`, not from the
/// allele reference. Mutations / NP base substitutions are
/// reflected — this is the *current* pool's stop-codon status, not
/// the germline's.
pub struct NoStopCodonInJunction;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum NpFilterCandidate {
    Length {
        segment: Segment,
        length: u32,
    },
    Base {
        segment: Segment,
        index: u32,
        total_len: Option<u32>,
        base: u8,
    },
}

impl NoStopCodonInJunction {
    pub fn new() -> Self {
        Self
    }

    fn parse_np_base_address(address: &str) -> Option<(Segment, u32)> {
        let (segment, prefix) = if address.starts_with("np.np1.bases[") {
            (Segment::Np1, "np.np1.bases[")
        } else if address.starts_with("np.np2.bases[") {
            (Segment::Np2, "np.np2.bases[")
        } else {
            return None;
        };
        let rest = address.strip_prefix(prefix)?;
        let idx = rest.strip_suffix(']')?.parse::<u32>().ok()?;
        Some((segment, idx))
    }

    fn parse_np_length_address(address: &str) -> Option<Segment> {
        match address {
            "np.np1.length" => Some(Segment::Np1),
            "np.np2.length" => Some(Segment::Np2),
            _ => None,
        }
    }

    fn parse_targeted_substitution_base_address(address: &str) -> Option<u32> {
        for prefix in [
            "mutate.uniform.base[",
            "mutate.s5f.base[",
            "corrupt.pcr.error_base[",
            "corrupt.contaminant.bases[",
        ] {
            if let Some(rest) = address.strip_prefix(prefix) {
                return rest.strip_suffix(']')?.parse::<u32>().ok();
            }
        }
        None
    }

    fn parse_filter_candidate(
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Option<NpFilterCandidate> {
        if let Some(segment) = Self::parse_np_length_address(address) {
            let length = match candidate {
                ChoiceValue::Int(n) if *n >= 0 => *n as u32,
                _ => return None,
            };
            return Some(NpFilterCandidate::Length { segment, length });
        }

        let (segment, index_from_address) = Self::parse_np_base_address(address)?;
        let base = match candidate {
            ChoiceValue::Base(b) => *b,
            _ => return None,
        };
        let index = context.draw_index.unwrap_or(index_from_address);
        Some(NpFilterCandidate::Base {
            segment,
            index,
            total_len: context.draw_count,
            base,
        })
    }

    fn append_pool_range(out: &mut Vec<Option<u8>>, sim: &Simulation, start: u32, end: u32) {
        for pos in start..end {
            out.push(sim.pool.get(NucHandle::new(pos)).map(|n| n.base));
        }
    }

    fn retained_bounds(allele_len: u32, trim_5: u16, trim_3: u16) -> Option<(u32, u32)> {
        let trim_5 = trim_5 as u32;
        let trim_3 = trim_3 as u32;
        if trim_5 + trim_3 > allele_len {
            return None;
        }
        Some((trim_5, allele_len - trim_3))
    }

    fn append_ref_slice(out: &mut Vec<Option<u8>>, seq: &[u8], start: u32, end: u32) {
        let end = end.min(seq.len() as u32);
        for pos in start..end {
            out.push(seq.get(pos as usize).copied());
        }
    }

    fn append_fixed_segment(
        out: &mut Vec<Option<u8>>,
        sim: &Simulation,
        refdata: &RefDataConfig,
        segment: Segment,
        allele_end_exclusive: Option<u32>,
    ) -> Option<()> {
        let inst = sim.assignments.get(segment).copied()?;
        let allele = refdata.get(segment, inst.allele_id)?;
        let (retained_start, retained_end) =
            Self::retained_bounds(allele.len(), inst.trim_5, inst.trim_3)?;
        let allele_end = allele_end_exclusive
            .unwrap_or(retained_end)
            .min(retained_end);
        if allele_end <= retained_start {
            return Some(());
        }

        if let Some(region) = sim.sequence.regions.iter().find(|r| r.segment == segment) {
            let pool_len = allele_end.saturating_sub(retained_start);
            let pool_end = region
                .start
                .index()
                .saturating_add(pool_len)
                .min(region.end.index());
            Self::append_pool_range(out, sim, region.start.index(), pool_end);
        } else {
            Self::append_ref_slice(out, &allele.seq, retained_start, allele_end);
        }

        Some(())
    }

    fn append_np_segment(
        out: &mut Vec<Option<u8>>,
        sim: &Simulation,
        segment: Segment,
        np_candidate: NpFilterCandidate,
    ) -> Option<bool> {
        if let Some(region) = sim.sequence.regions.iter().find(|r| r.segment == segment) {
            Self::append_pool_range(out, sim, region.start.index(), region.end.index());
            return Some(true);
        }

        match np_candidate {
            NpFilterCandidate::Length {
                segment: candidate_segment,
                length,
            } if candidate_segment == segment => {
                out.extend((0..length).map(|_| None));
                Some(true)
            }
            NpFilterCandidate::Base {
                segment: candidate_segment,
                index,
                total_len,
                base,
            } if candidate_segment == segment => {
                let count_to_append = total_len.unwrap_or(index.saturating_add(1));
                if count_to_append <= index {
                    return Some(false);
                }
                let current_np_start = (sim.pool.len() as u32).checked_sub(index)?;
                for i in 0..count_to_append {
                    if i < index {
                        let pos = current_np_start + i;
                        out.push(sim.pool.get(NucHandle::new(pos)).map(|n| n.base));
                    } else if i == index {
                        out.push(Some(base));
                    } else {
                        out.push(None);
                    }
                }
                Some(total_len.is_some())
            }
            _ => Some(false),
        }
    }

    fn hypothetical_junction_bases(
        sim: &Simulation,
        refdata: &RefDataConfig,
        np_candidate: NpFilterCandidate,
    ) -> Option<Vec<Option<u8>>> {
        let v_inst = sim.assignments.v?;
        let j_inst = sim.assignments.j?;
        let v_allele = refdata.get(Segment::V, v_inst.allele_id)?;
        let j_allele = refdata.get(Segment::J, j_inst.allele_id)?;
        let v_anchor = v_allele.anchor? as u32;
        let j_anchor = j_allele.anchor? as u32;

        let (v_retained_start, v_retained_end) =
            Self::retained_bounds(v_allele.len(), v_inst.trim_5, v_inst.trim_3)?;
        let (j_retained_start, j_retained_end) =
            Self::retained_bounds(j_allele.len(), j_inst.trim_5, j_inst.trim_3)?;
        if v_anchor < v_retained_start || v_anchor >= v_retained_end {
            return None;
        }
        if j_anchor < j_retained_start {
            return None;
        }

        let v_region = sim
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)?;
        let v_anchor_pool = v_region.start.index() + (v_anchor - v_retained_start);

        let mut bases = Vec::new();
        Self::append_pool_range(&mut bases, sim, v_anchor_pool, v_region.end.index());

        if sim.assignments.d.is_some() {
            if !Self::append_np_segment(&mut bases, sim, Segment::Np1, np_candidate)? {
                return Some(bases);
            }
            Self::append_fixed_segment(&mut bases, sim, refdata, Segment::D, None)?;
            if !Self::append_np_segment(&mut bases, sim, Segment::Np2, np_candidate)? {
                return Some(bases);
            }
        } else if !Self::append_np_segment(&mut bases, sim, Segment::Np1, np_candidate)? {
            return Some(bases);
        }

        let j_junction_end = j_anchor.saturating_add(3).min(j_retained_end);
        Self::append_fixed_segment(&mut bases, sim, refdata, Segment::J, Some(j_junction_end))?;

        Some(bases)
    }

    fn reject_known_stop(
        &self,
        address: &str,
        bases: &[Option<u8>],
    ) -> Result<(), ContractViolation> {
        let mut offset = 0usize;
        while offset + 3 <= bases.len() {
            if let (Some(b1), Some(b2), Some(b3)) =
                (bases[offset], bases[offset + 1], bases[offset + 2])
            {
                if translate_codon(b1, b2, b3) == AMINO_STOP {
                    return Err(ContractViolation::new(
                        self.name(),
                        format!(
                            "candidate at {} would force stop codon '{}{}{}' at junction offset {}",
                            address, b1 as char, b2 as char, b3 as char, offset
                        ),
                    ));
                }
            }
            offset += 3;
        }
        Ok(())
    }

    fn admits_np_candidate(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        let np_candidate = match Self::parse_filter_candidate(address, candidate, context) {
            None => return Ok(()),
            Some(c) => c,
        };
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };
        let bases = match Self::hypothetical_junction_bases(sim, refdata, np_candidate) {
            None => return Ok(()),
            Some(b) => b,
        };
        self.reject_known_stop(address, &bases)
    }

    fn pool_base_with_candidate(
        sim: &Simulation,
        handle: NucHandle,
        target: NucHandle,
        candidate_base: u8,
    ) -> Option<u8> {
        if handle == target {
            Some(candidate_base)
        } else {
            sim.pool.get(handle).map(|n| n.base)
        }
    }

    fn admits_targeted_substitution_candidate(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        if Self::parse_targeted_substitution_base_address(address).is_none() {
            return Ok(());
        }

        let candidate_base = match candidate {
            ChoiceValue::Base(b) => *b,
            _ => return Ok(()),
        };
        let target = match context.target {
            Some(h) => h,
            None => return Ok(()),
        };
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };
        let junction = match compute_junction(sim, refdata) {
            None => return Ok(()),
            Some(j) => j,
        };

        // Frame violations are owned by ProductiveJunctionFrame. This
        // contract only reasons over codon-shaped triples.
        if !junction.is_in_frame() {
            return Ok(());
        }

        let start = junction.start.index();
        let end = junction.end.index();
        let target_idx = target.index();
        if target_idx < start || target_idx >= end {
            return Ok(());
        }

        let mut pos = start;
        while pos + 3 <= end {
            let b1 =
                Self::pool_base_with_candidate(sim, NucHandle::new(pos), target, candidate_base);
            let b2 = Self::pool_base_with_candidate(
                sim,
                NucHandle::new(pos + 1),
                target,
                candidate_base,
            );
            let b3 = Self::pool_base_with_candidate(
                sim,
                NucHandle::new(pos + 2),
                target,
                candidate_base,
            );

            if let (Some(b1), Some(b2), Some(b3)) = (b1, b2, b3) {
                if translate_codon(b1, b2, b3) == AMINO_STOP {
                    return Err(ContractViolation::new(
                        self.name(),
                        format!(
                            "candidate at {} for target {} would leave stop codon '{}{}{}' at junction position {}",
                            address,
                            target_idx,
                            b1 as char,
                            b2 as char,
                            b3 as char,
                            pos
                        ),
                    ));
                }
            }
            pos += 3;
        }

        Ok(())
    }
}

impl Default for NoStopCodonInJunction {
    fn default() -> Self {
        Self::new()
    }
}

impl Contract for NoStopCodonInJunction {
    fn name(&self) -> &str {
        "no_stop_codon_in_junction"
    }

    fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), ContractViolation> {
        let refdata = match refdata {
            None => return Ok(()),
            Some(r) => r,
        };

        let junction = match compute_junction(sim, refdata) {
            None => return Ok(()),
            Some(j) => j,
        };

        // Out of frame → not our concern. ProductiveJunctionFrame
        // is the contract that owns the frame check; here we'd just
        // be walking codon-shaped triples that aren't real codons.
        if !junction.is_in_frame() {
            return Ok(());
        }

        let start = junction.start.index();
        let end = junction.end.index();

        let mut pos = start;
        while pos + 3 <= end {
            // Defensive lookups — pool pointers from the junction
            // window should always be valid for a well-formed
            // simulation, but we don't want to panic on a malformed
            // input. Treat missing handles as ambiguous.
            let b1 = sim
                .pool
                .get(NucHandle::new(pos))
                .map(|n| n.base)
                .unwrap_or(b'N');
            let b2 = sim
                .pool
                .get(NucHandle::new(pos + 1))
                .map(|n| n.base)
                .unwrap_or(b'N');
            let b3 = sim
                .pool
                .get(NucHandle::new(pos + 2))
                .map(|n| n.base)
                .unwrap_or(b'N');

            let aa = translate_codon(b1, b2, b3);
            if aa == AMINO_STOP {
                return Err(ContractViolation::new(
                    self.name(),
                    format!(
                        "stop codon at junction position {}: '{}{}{}'",
                        pos, b1 as char, b2 as char, b3 as char
                    ),
                ));
            }
            pos += 3;
        }

        Ok(())
    }

    fn admits(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        self.admits_np_candidate(sim, refdata, address, candidate, ChoiceContext::none())
    }

    fn admits_with_context(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        self.admits_np_candidate(sim, refdata, address, candidate, context)?;
        self.admits_targeted_substitution_candidate(sim, refdata, address, candidate, context)
    }
}

// ──────────────────────────────────────────────────────────────────
// ContractSet — composition of multiple contracts (D.5)
// ──────────────────────────────────────────────────────────────────

/// A set of contracts that all must hold for the simulation.
///
/// Composition is *intersection*: `verify` succeeds only if every
/// contained contract verifies; `admits` accepts a candidate only
/// if every contained contract admits it. This mirrors the
/// architectural commitment from D6 — "Multiple contract specs
/// intersect — all must hold throughout the simulation."
///
/// **`verify` collects all violations** (not first-only) so
/// diagnostic output shows the full picture of what failed.
///
/// **`admits` short-circuits** on the first inadmissible contract
/// because sampling speed matters and one violation is enough to
/// reject the candidate.
pub struct ContractSet {
    contracts: Vec<Box<dyn Contract>>,
}

impl ContractSet {
    /// Empty set — admits everything, verifies anything.
    pub fn new() -> Self {
        Self {
            contracts: Vec::new(),
        }
    }

    /// Builder-style append. Returns `self` for chaining.
    #[must_use]
    pub fn with(mut self, contract: Box<dyn Contract>) -> Self {
        self.contracts.push(contract);
        self
    }

    /// In-place append. Returns `&mut self` for builder-mut chains.
    pub fn add(&mut self, contract: Box<dyn Contract>) -> &mut Self {
        self.contracts.push(contract);
        self
    }

    /// Number of contained contracts.
    pub fn len(&self) -> usize {
        self.contracts.len()
    }

    /// Whether the set contains zero contracts.
    pub fn is_empty(&self) -> bool {
        self.contracts.is_empty()
    }

    /// Iterate over the contained contracts.
    pub fn iter(&self) -> impl Iterator<Item = &dyn Contract> {
        self.contracts.iter().map(|c| c.as_ref())
    }

    /// Verify every contract. Returns `Ok(())` only if all pass.
    /// On failure returns the *complete* list of violations, so
    /// callers see every contract that failed in one pass.
    pub fn verify(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
    ) -> Result<(), Vec<ContractViolation>> {
        let mut violations = Vec::new();
        for c in &self.contracts {
            if let Err(v) = c.verify(sim, refdata) {
                violations.push(v);
            }
        }
        if violations.is_empty() {
            Ok(())
        } else {
            Err(violations)
        }
    }

    /// Test whether `candidate` at `address` is admissible by every
    /// contract. Short-circuits on the first violator (sampling
    /// hot path). Returns the violator's `ContractViolation`.
    pub fn admits(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
    ) -> Result<(), ContractViolation> {
        for c in &self.contracts {
            c.admits(sim, refdata, address, candidate)?;
        }
        Ok(())
    }

    /// Context-aware variant of [`ContractSet::admits`].
    pub fn admits_with_context(
        &self,
        sim: &Simulation,
        refdata: Option<&RefDataConfig>,
        address: &str,
        candidate: &ChoiceValue,
        context: ChoiceContext,
    ) -> Result<(), ContractViolation> {
        for c in &self.contracts {
            c.admits_with_context(sim, refdata, address, candidate, context)?;
        }
        Ok(())
    }
}

impl Default for ContractSet {
    fn default() -> Self {
        Self::new()
    }
}

// ──────────────────────────────────────────────────────────────────
// productive() — the canonical productive-sequence bundle
// ──────────────────────────────────────────────────────────────────

/// The canonical productive-sequence contract bundle:
///
/// 1. `ProductiveJunctionFrame` — junction length divisible by 3
/// 2. `NoStopCodonInJunction` — no stops inside the junction
/// 3. `AnchorPreserved::V` — V Cys codon retained after V trim
/// 4. `AnchorPreserved::J` — J W/F codon retained after J trim
///
/// All four must hold for a sequence to be considered productive
/// in the standard biological sense. The DSL `respect=[productive()]`
/// (Phase F) compiles to this bundle internally.
pub fn productive() -> ContractSet {
    ContractSet::new()
        .with(Box::new(ProductiveJunctionFrame::new()))
        .with(Box::new(NoStopCodonInJunction::new()))
        .with(Box::new(AnchorPreserved::new(Segment::V)))
        .with(Box::new(AnchorPreserved::new(Segment::J)))
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::ir::Segment;
    use crate::refdata::{Allele, AlleleId, AllelePool, ChainType, RefDataConfig};

    fn make_v_anchor_at(len: u32, anchor: Option<u16>) -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_test*01".to_string(),
            gene: "v_test".to_string(),
            seq: vec![b'A'; len as usize],
            segment: Segment::V,
            anchor,
        });
        cfg
    }

    #[test]
    #[should_panic(expected = "segment must be V, D, or J")]
    fn anchor_preserved_rejects_np_segment() {
        let _ = AnchorPreserved::new(Segment::Np1);
    }

    #[test]
    fn anchor_preserved_name_per_segment() {
        assert_eq!(
            AnchorPreserved::new(Segment::V).name(),
            "anchor_preserved.v"
        );
        assert_eq!(
            AnchorPreserved::new(Segment::D).name(),
            "anchor_preserved.d"
        );
        assert_eq!(
            AnchorPreserved::new(Segment::J).name(),
            "anchor_preserved.j"
        );
    }

    #[test]
    fn anchor_preserved_no_allele_assigned_passes_vacuously() {
        let cfg = make_v_anchor_at(10, Some(5));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new();

        // No assignment → vacuously satisfied.
        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_no_refdata_passes_vacuously() {
        let cfg = make_v_anchor_at(10, Some(5));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        // Refdata None → can't verify; treat as satisfied.
        assert!(contract.verify(&sim, None).is_ok());
        // Sanity: WITH refdata it would also pass (no trims).
        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_anchorless_allele_passes_vacuously() {
        let cfg = make_v_anchor_at(10, None); // anchorless
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_zero_trims_and_anchor_in_range_passes() {
        // 10-base allele, anchor at 3. No trims. Anchor codon spans
        // [3, 6) which fits entirely in [0, 10).
        let cfg = make_v_anchor_at(10, Some(3));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_violates_when_5prime_trim_eats_anchor() {
        // Anchor at 3. trim_5 = 4 → anchor 5'-trimmed away.
        use crate::assignment::TrimEnd;

        let cfg = make_v_anchor_at(10, Some(3));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_trim(Segment::V, TrimEnd::Five, 4);

        let err = contract.verify(&sim, Some(&cfg)).unwrap_err();
        assert_eq!(err.contract_name, "anchor_preserved.v");
        assert!(
            err.reason.contains("trim_5"),
            "reason should mention trim_5, got: {}",
            err.reason
        );
    }

    #[test]
    fn anchor_preserved_violates_when_3prime_trim_eats_anchor() {
        // 10-base allele, anchor at 6. Anchor codon spans [6, 9).
        // trim_3 = 2 → retained_end = 8, so anchor codon's last
        // position (8) is at the boundary — anchor + 3 = 9 > 8 → fail.
        use crate::assignment::TrimEnd;

        let cfg = make_v_anchor_at(10, Some(6));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_trim(Segment::V, TrimEnd::Three, 2);

        let err = contract.verify(&sim, Some(&cfg)).unwrap_err();
        assert_eq!(err.contract_name, "anchor_preserved.v");
        assert!(
            err.reason.contains("retained slice ends at"),
            "reason should describe retained slice, got: {}",
            err.reason
        );
    }

    #[test]
    fn anchor_preserved_anchor_at_5prime_boundary_passes() {
        // anchor at 0 → trim_5 = 0 leaves anchor at the very start.
        let cfg = make_v_anchor_at(10, Some(0));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_anchor_at_3prime_boundary_passes() {
        // 10-base allele, anchor at 7. Anchor codon spans [7, 10).
        // trim_3 = 0 → retained_end = 10. Boundary case — passes.
        let cfg = make_v_anchor_at(10, Some(7));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn anchor_preserved_one_more_trim_just_past_boundary_violates() {
        // Same as above but trim_3 = 1 → retained_end = 9, anchor
        // codon end (10) > 9 → fail.
        use crate::assignment::TrimEnd;

        let cfg = make_v_anchor_at(10, Some(7));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_trim(Segment::V, TrimEnd::Three, 1);

        assert!(contract.verify(&sim, Some(&cfg)).is_err());
    }

    #[test]
    fn anchor_preserved_works_for_j_segment_too() {
        // Build a J pool with an anchor at position 0 (typical J
        // anchor is near the 5' end). trim_5 = 0 → anchor preserved.
        // trim_5 = 1 → anchor 5'-trimmed → fail.
        use crate::assignment::TrimEnd;

        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.j_pool.push(Allele {
            name: "j_test*01".to_string(),
            gene: "j_test".to_string(),
            seq: vec![b'A'; 12],
            segment: Segment::J,
            anchor: Some(0),
        });
        let contract = AnchorPreserved::new(Segment::J);

        let ok = Simulation::new()
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
        assert!(contract.verify(&ok, Some(&cfg)).is_ok());

        let bad = ok.with_trim(Segment::J, TrimEnd::Five, 1);
        assert!(contract.verify(&bad, Some(&cfg)).is_err());
    }

    #[test]
    fn contract_violation_construction_round_trip() {
        let v = ContractViolation::new("test.contract", "something failed");
        assert_eq!(v.contract_name, "test.contract");
        assert_eq!(v.reason, "something failed");
    }

    #[test]
    fn contract_works_through_box_dyn() {
        // Trait-object usage is dyn-compatible.
        let cfg = make_v_anchor_at(10, Some(3));
        let contract: Box<dyn Contract> = Box::new(AnchorPreserved::new(Segment::V));
        let sim = Simulation::new()
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)));

        assert!(contract.verify(&sim, Some(&cfg)).is_ok());
        assert_eq!(contract.name(), "anchor_preserved.v");
    }

    // ── ProductiveJunctionFrame tests (D.2) ───────────────────────

    use crate::ir::{NucHandle, Nucleotide, Region};

    /// Build a complete V+J refdata with controllable anchors.
    /// V: "AAACCCGGG" (9bp, anchor at v_anchor)
    /// J: "TTTAAA" (6bp, anchor at j_anchor)
    fn make_vj_for_frame_test(v_anchor: Option<u16>, j_anchor: Option<u16>) -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(crate::refdata::ChainType::Vdj);
        let _ = cfg.v_pool.push(crate::refdata::Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: b"AAACCCGGG".to_vec(),
            segment: Segment::V,
            anchor: v_anchor,
        });
        let _ = cfg.j_pool.push(crate::refdata::Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: j_anchor,
        });
        cfg
    }

    fn make_assembled_sim(
        v_pool_start: u32,
        v_len: u32,
        j_pool_start: u32,
        j_len: u32,
        v_inst: AlleleInstance,
        j_inst: AlleleInstance,
    ) -> Simulation {
        let mut sim = Simulation::new();
        for i in 0..v_pool_start {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'X', i as u16, Segment::V));
            sim = next;
        }
        for i in 0..v_len {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'A', i as u16, Segment::V));
            sim = next;
        }
        for _ in v_pool_start + v_len..j_pool_start {
            let (next, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
                b'a',
                Segment::Np1,
                crate::ir::flag::N_NUC,
            ));
            sim = next;
        }
        for i in 0..j_len {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b'T', i as u16, Segment::J));
            sim = next;
        }
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(v_pool_start),
            NucHandle::new(v_pool_start + v_len),
        ));
        sim = sim.with_region_added(Region::new(
            Segment::J,
            NucHandle::new(j_pool_start),
            NucHandle::new(j_pool_start + j_len),
        ));
        sim.with_allele_assigned(Segment::V, v_inst)
            .with_allele_assigned(Segment::J, j_inst)
    }

    #[test]
    fn productive_junction_frame_name_is_canonical() {
        let c = ProductiveJunctionFrame::new();
        assert_eq!(c.name(), "productive_junction_frame");
    }

    #[test]
    fn productive_junction_frame_no_refdata_passes_vacuously() {
        let c = ProductiveJunctionFrame::new();
        let sim = Simulation::new();
        assert!(c.verify(&sim, None).is_ok());
    }

    #[test]
    fn productive_junction_frame_no_v_assignment_passes_vacuously() {
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let c = ProductiveJunctionFrame::new();
        let sim = Simulation::new()
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn productive_junction_frame_anchorless_v_passes_vacuously() {
        // Junction undefined → contract has nothing to check.
        let cfg = make_vj_for_frame_test(None, Some(0));
        let c = ProductiveJunctionFrame::new();
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 9, 6, v_inst, j_inst);
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn productive_junction_frame_in_frame_junction_passes() {
        // V anchor 6 in 9bp V → V anchor pool = 6.
        // J anchor 0 in 6bp J → J anchor pool = 9.
        // Junction = [6, 12). Length = 6 (divisible by 3).
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let c = ProductiveJunctionFrame::new();
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 9, 6, v_inst, j_inst);
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn productive_junction_frame_out_of_frame_junction_violates() {
        // V anchor 6, J anchor 2, J trim_5=1 → J anchor in trimmed J = 1.
        // J at pool [9, 14), J anchor pool = 9 + 1 = 10. End = 13.
        // V anchor pool = 6. Junction = [6, 13). Length 7 → not divisible by 3.
        let cfg = make_vj_for_frame_test(Some(6), Some(2));
        let c = ProductiveJunctionFrame::new();
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0)).with_trim_5(1);
        let sim = make_assembled_sim(0, 9, 9, 5, v_inst, j_inst);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert_eq!(err.contract_name, "productive_junction_frame");
        assert!(
            err.reason.contains("not divisible by 3"),
            "reason should mention frame: {}",
            err.reason
        );
        assert!(
            err.reason.contains("7"),
            "reason should mention length 7: {}",
            err.reason
        );
    }

    #[test]
    fn productive_junction_frame_in_frame_with_np_padding_passes() {
        // V at [0, 9). NP1 padding 3bp at [9, 12). J at [12, 18).
        // V anchor pool = 6. J anchor pool = 12. Junction = [6, 15). Length 9.
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let c = ProductiveJunctionFrame::new();
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 12, 6, v_inst, j_inst);
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn productive_junction_frame_works_through_box_dyn() {
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 9, 6, v_inst, j_inst);

        let c: Box<dyn Contract> = Box::new(ProductiveJunctionFrame::new());
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
        assert_eq!(c.name(), "productive_junction_frame");
    }

    // ── NoStopCodonInJunction tests (D.3) ──────────────────────────

    /// Build a refdata with V allele "AAACCCxxx" (anchor at 6, so
    /// junction starts at AAACCCxxx[6..9] = "xxx") and J allele
    /// "yyyAAA" (anchor at 0, so junction ends at yyyAAA[0..3] +
    /// 3 = first 3 chars). The exact bases at the V anchor codon
    /// and J anchor codon control whether the junction has stops.
    fn make_vj_with_anchor_codons(
        v_anchor_codon: &[u8; 3],
        j_anchor_codon: &[u8; 3],
    ) -> RefDataConfig {
        let mut v_seq = b"AAACCC".to_vec();
        v_seq.extend_from_slice(v_anchor_codon);
        // V is now 9 bases, anchor at 6.

        let mut j_seq = j_anchor_codon.to_vec();
        j_seq.extend_from_slice(b"AAA");
        // J is now 6 bases, anchor at 0.

        let mut cfg = RefDataConfig::empty(crate::refdata::ChainType::Vdj);
        let _ = cfg.v_pool.push(crate::refdata::Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: v_seq,
            segment: Segment::V,
            anchor: Some(6),
        });
        let _ = cfg.j_pool.push(crate::refdata::Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: j_seq,
            segment: Segment::J,
            anchor: Some(0),
        });
        cfg
    }

    /// Build a simulation that copies the V/J sequences from `cfg`
    /// into the pool with the given pool placement, no NP padding.
    /// V anchor codon ends at V[8], J anchor codon at J[0..3], so
    /// junction is `[V_anchor_pool, J_anchor_pool + 3) = [6, 12)`,
    /// which spans the V anchor codon (3 bases) + 3 bases of J.
    /// Length 6 → in frame.
    fn make_assembled_sim_from_refdata(cfg: &RefDataConfig) -> Simulation {
        let v_allele = cfg.v_pool.get(AlleleId::new(0)).unwrap();
        let j_allele = cfg.j_pool.get(AlleleId::new(0)).unwrap();

        let mut sim = Simulation::new();
        for (i, &b) in v_allele.seq.iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        for (i, &b) in j_allele.seq.iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::J));
            sim = next;
        }
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(v_allele.seq.len() as u32),
        ));
        sim = sim.with_region_added(Region::new(
            Segment::J,
            NucHandle::new(v_allele.seq.len() as u32),
            NucHandle::new((v_allele.seq.len() + j_allele.seq.len()) as u32),
        ));
        sim.with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)))
    }

    #[test]
    fn no_stop_codon_name_is_canonical() {
        assert_eq!(
            NoStopCodonInJunction::new().name(),
            "no_stop_codon_in_junction"
        );
    }

    #[test]
    fn no_stop_codon_no_refdata_passes_vacuously() {
        let c = NoStopCodonInJunction::new();
        assert!(c.verify(&Simulation::new(), None).is_ok());
    }

    #[test]
    fn no_stop_codon_undefined_junction_passes_vacuously() {
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        let c = NoStopCodonInJunction::new();
        // No assignments at all → junction undefined.
        assert!(c.verify(&Simulation::new(), Some(&cfg)).is_ok());
    }

    #[test]
    fn no_stop_codon_in_frame_no_stops_passes() {
        // V anchor codon GGG (Gly), J anchor codon TTT (Phe).
        // Junction codons: GGG, TTT → no stops.
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn no_stop_codon_admits_rejects_uniform_mutation_base_that_creates_stop() {
        // Junction codons start as TAC, GGG. Mutating the third base
        // of TAC to A would create TAA.
        let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "mutate.uniform.base[0]",
                &ChoiceValue::Base(b'A'),
                ChoiceContext::indexed_target(0, 1, NucHandle::new(8)),
            )
            .unwrap_err();

        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);
        assert!(c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "mutate.uniform.base[0]",
                &ChoiceValue::Base(b'C'),
                ChoiceContext::indexed_target(0, 1, NucHandle::new(8)),
            )
            .is_ok());
    }

    #[test]
    fn no_stop_codon_admits_rejects_s5f_mutation_base_that_creates_stop() {
        // Same target-site contract path as uniform mutation, but
        // dispatched through the biologically important S5F address.
        let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "mutate.s5f.base[0]",
                &ChoiceValue::Base(b'A'),
                ChoiceContext::indexed_target(0, 1, NucHandle::new(8)),
            )
            .unwrap_err();

        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_admits_rejects_pcr_error_base_that_creates_stop() {
        // Observation-stage substitutions use the same target-site
        // contract path as biological mutations.
        let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "corrupt.pcr.error_base[0]",
                &ChoiceValue::Base(b'A'),
                ChoiceContext::indexed_target(0, 1, NucHandle::new(8)),
            )
            .unwrap_err();

        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_admits_rejects_contaminant_base_that_creates_stop() {
        // Contaminant replacement is also a target-site substitution.
        let cfg = make_vj_with_anchor_codons(b"AAA", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "corrupt.contaminant.bases[6]",
                &ChoiceValue::Base(b'T'),
                ChoiceContext::indexed_target(6, 12, NucHandle::new(6)),
            )
            .unwrap_err();

        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_v_anchor_taa_violates() {
        // V anchor codon TAA (stop).
        let cfg = make_vj_with_anchor_codons(b"TAA", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("stop"), "{}", err.reason);
        assert!(err.reason.contains("TAA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_v_anchor_tag_violates() {
        let cfg = make_vj_with_anchor_codons(b"TAG", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert!(err.reason.contains("TAG"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_v_anchor_tga_violates() {
        let cfg = make_vj_with_anchor_codons(b"TGA", b"GGG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert!(err.reason.contains("TGA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_j_anchor_stop_violates() {
        // V anchor GGG (fine), J anchor TAA (stop).
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TAA");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert!(err.reason.contains("TAA"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_first_stop_reported() {
        // V anchor TAA (first), J anchor TAG (second). Should report V's TAA.
        let cfg = make_vj_with_anchor_codons(b"TAA", b"TAG");
        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);

        let err = c.verify(&sim, Some(&cfg)).unwrap_err();
        assert!(err.reason.contains("TAA"), "{}", err.reason);
        // Position 6 is the V anchor pool position.
        assert!(err.reason.contains("position 6"), "{}", err.reason);
    }

    #[test]
    fn no_stop_codon_out_of_frame_junction_passes_vacuously() {
        // J anchor at 2 → J anchor in pool = 9 + 2 = 11.
        // V anchor at 6 → V anchor pool = 6. Junction = [6, 14). Length 8.
        // 8 % 3 != 0 → out of frame → contract returns Ok (deferred to
        // ProductiveJunctionFrame).
        let mut cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        // Move J anchor to 2 to break frame.
        let j_allele = cfg.j_pool.get(AlleleId::new(0)).unwrap().clone();
        cfg.j_pool = AllelePool::from_vec(vec![crate::refdata::Allele {
            anchor: Some(2),
            ..j_allele
        }]);

        let c = NoStopCodonInJunction::new();
        let sim = make_assembled_sim_from_refdata(&cfg);
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn no_stop_codon_works_through_box_dyn() {
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        let sim = make_assembled_sim_from_refdata(&cfg);
        let c: Box<dyn Contract> = Box::new(NoStopCodonInJunction::new());
        assert!(c.verify(&sim, Some(&cfg)).is_ok());
        assert_eq!(c.name(), "no_stop_codon_in_junction");
    }

    fn make_partial_np_stop_filter_case() -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(crate::refdata::ChainType::Vj);
        let _ = cfg.v_pool.push(crate::refdata::Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: b"GGGTA".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = cfg.j_pool.push(crate::refdata::Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        for (i, &b) in b"GGGTA".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(5),
        ));
        sim = sim
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

        (cfg, sim)
    }

    fn make_partial_np1_future_d_stop_case(v_seq: &[u8]) -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(crate::refdata::ChainType::Vdj);
        let _ = cfg.v_pool.push(crate::refdata::Allele {
            name: "v_test*01".into(),
            gene: "v_test".into(),
            seq: v_seq.to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = cfg.d_pool.push(crate::refdata::Allele {
            name: "d_test*01".into(),
            gene: "d_test".into(),
            seq: b"AAC".to_vec(),
            segment: Segment::D,
            anchor: None,
        });
        let _ = cfg.j_pool.push(crate::refdata::Allele {
            name: "j_test*01".into(),
            gene: "j_test".into(),
            seq: b"TTTAAA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        for (i, &b) in v_seq.iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        sim = sim.with_region_added(Region::new(
            Segment::V,
            NucHandle::new(0),
            NucHandle::new(v_seq.len() as u32),
        ));
        sim = sim
            .with_allele_assigned(Segment::V, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::D, AlleleInstance::new(AlleleId::new(0)))
            .with_allele_assigned(Segment::J, AlleleInstance::new(AlleleId::new(0)));

        (cfg, sim)
    }

    #[test]
    fn no_stop_codon_admits_rejects_np_base_that_completes_stop() {
        let (cfg, sim) = make_partial_np_stop_filter_case();
        let c = NoStopCodonInJunction::new();

        let err = c
            .admits(
                &sim,
                Some(&cfg),
                "np.np1.bases[0]",
                &ChoiceValue::Base(b'A'),
            )
            .unwrap_err();
        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);

        assert!(c
            .admits(
                &sim,
                Some(&cfg),
                "np.np1.bases[0]",
                &ChoiceValue::Base(b'C')
            )
            .is_ok());
    }

    #[test]
    fn no_stop_codon_admits_ignores_np_base_until_codon_complete() {
        let (cfg, sim) = make_partial_np_stop_filter_case();
        let c = NoStopCodonInJunction::new();
        let (sim, _) = sim.with_nucleotide_pushed(Nucleotide::synthetic(
            b'T',
            Segment::Np1,
            crate::ir::flag::N_NUC,
        ));

        // Candidate position is offset 6 from the junction start, so it
        // starts a new codon rather than completing one. The future third
        // base will be filtered when that draw is reached.
        assert!(c
            .admits(
                &sim,
                Some(&cfg),
                "np.np1.bases[1]",
                &ChoiceValue::Base(b'A')
            )
            .is_ok());
    }

    #[test]
    fn no_stop_codon_admits_rejects_np1_base_that_forces_future_d_stop() {
        let (cfg, sim) = make_partial_np1_future_d_stop_case(b"GGG");
        let c = NoStopCodonInJunction::new();

        let err = c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "np.np1.bases[0]",
                &ChoiceValue::Base(b'T'),
                ChoiceContext::indexed(0, 1),
            )
            .unwrap_err();
        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);

        assert!(c
            .admits_with_context(
                &sim,
                Some(&cfg),
                "np.np1.bases[0]",
                &ChoiceValue::Base(b'C'),
                ChoiceContext::indexed(0, 1),
            )
            .is_ok());
    }

    #[test]
    fn no_stop_codon_admits_rejects_zero_np1_length_that_forces_future_d_stop() {
        let (cfg, sim) = make_partial_np1_future_d_stop_case(b"GGGT");
        let c = NoStopCodonInJunction::new();

        let err = c
            .admits(&sim, Some(&cfg), "np.np1.length", &ChoiceValue::Int(0))
            .unwrap_err();
        assert_eq!(err.contract_name, "no_stop_codon_in_junction");
        assert!(err.reason.contains("TAA"), "{}", err.reason);

        assert!(c
            .admits(&sim, Some(&cfg), "np.np1.length", &ChoiceValue::Int(1),)
            .is_ok());
    }

    // ── Contract::admits default behaviour (D.4) ──────────────────

    #[test]
    fn default_admits_returns_ok_for_anchor_preserved() {
        let cfg = make_v_anchor_at(10, Some(3));
        let contract = AnchorPreserved::new(Segment::V);
        let sim = Simulation::new();

        // Default admits ignores the address + candidate and returns Ok.
        assert!(contract
            .admits(&sim, Some(&cfg), "trim.v_3", &ChoiceValue::Int(5))
            .is_ok());
        assert!(contract
            .admits(&sim, None, "np.np1.length", &ChoiceValue::Int(3))
            .is_ok());
        assert!(contract
            .admits(
                &sim,
                Some(&cfg),
                "sample_allele.v",
                &ChoiceValue::AlleleId(0)
            )
            .is_ok());
    }

    #[test]
    fn default_admits_returns_ok_for_productive_junction_frame() {
        let contract = ProductiveJunctionFrame::new();
        let sim = Simulation::new();

        // Same default behaviour for the new contract — D.6 will
        // override `admits` on this contract specifically; for now
        // it inherits the trait default.
        assert!(contract
            .admits(&sim, None, "np.np1.length", &ChoiceValue::Int(7))
            .is_ok());
    }

    #[test]
    fn no_stop_codon_admits_passes_vacuously_without_refdata() {
        let contract = NoStopCodonInJunction::new();
        let sim = Simulation::new();

        assert!(contract
            .admits(&sim, None, "np.np1.bases[0]", &ChoiceValue::Base(b'T'))
            .is_ok());
    }

    #[test]
    fn admits_works_through_box_dyn() {
        let contract: Box<dyn Contract> = Box::new(AnchorPreserved::new(Segment::V));
        let sim = Simulation::new();
        assert!(contract
            .admits(&sim, None, "trim.v_3", &ChoiceValue::Int(0))
            .is_ok());
    }

    // ── ContractSet tests (D.5) ────────────────────────────────────

    #[test]
    fn contract_set_empty_admits_everything() {
        let s = ContractSet::new();
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);

        let sim = Simulation::new();
        assert!(s.verify(&sim, None).is_ok());
        assert!(s
            .admits(&sim, None, "any.address", &ChoiceValue::Int(42))
            .is_ok());
    }

    #[test]
    fn contract_set_with_chains_contracts() {
        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(ProductiveJunctionFrame::new()));
        assert_eq!(s.len(), 2);
        assert!(!s.is_empty());

        let names: Vec<&str> = s.iter().map(|c| c.name()).collect();
        assert_eq!(
            names,
            vec!["anchor_preserved.v", "productive_junction_frame"]
        );
    }

    #[test]
    fn contract_set_add_in_place() {
        let mut s = ContractSet::new();
        s.add(Box::new(AnchorPreserved::new(Segment::V)));
        s.add(Box::new(AnchorPreserved::new(Segment::J)));
        assert_eq!(s.len(), 2);
    }

    #[test]
    fn contract_set_verify_succeeds_when_all_pass() {
        let cfg = make_vj_for_frame_test(Some(6), Some(0));
        let v_inst = AlleleInstance::new(AlleleId::new(0));
        let j_inst = AlleleInstance::new(AlleleId::new(0));
        let sim = make_assembled_sim(0, 9, 9, 6, v_inst, j_inst);

        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(ProductiveJunctionFrame::new()));
        assert!(s.verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn contract_set_verify_returns_only_failing_contracts() {
        // V anchor codon TAA → NoStopCodonInJunction violates.
        // V trim_3 = 0, V trim_5 = 0 → AnchorPreserved.V passes.
        // → exactly one violation.
        let cfg = make_vj_with_anchor_codons(b"TAA", b"GGG");
        let sim = make_assembled_sim_from_refdata(&cfg);

        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(NoStopCodonInJunction::new()))
            .with(Box::new(ProductiveJunctionFrame::new()));

        let violations = s.verify(&sim, Some(&cfg)).unwrap_err();
        assert_eq!(violations.len(), 1);
        assert_eq!(violations[0].contract_name, "no_stop_codon_in_junction");
    }

    #[test]
    fn contract_set_verify_collects_multiple_real_violations() {
        // Force two genuine violations: AnchorPreserved.V (trim_3=99)
        // AND a contract that always fails. Use a custom test contract.
        struct AlwaysFail;
        impl Contract for AlwaysFail {
            fn name(&self) -> &str {
                "always_fail"
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), ContractViolation> {
                Err(ContractViolation::new("always_fail", "always fails"))
            }
        }

        let cfg = make_v_anchor_at(10, Some(3));
        let bad_sim = Simulation::new().with_allele_assigned(
            Segment::V,
            AlleleInstance::new(AlleleId::new(0)).with_trim_5(5),
        );

        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(AlwaysFail));

        let violations = s.verify(&bad_sim, Some(&cfg)).unwrap_err();
        assert_eq!(violations.len(), 2);
        let names: Vec<&str> = violations
            .iter()
            .map(|v| v.contract_name.as_str())
            .collect();
        assert!(names.contains(&"anchor_preserved.v"));
        assert!(names.contains(&"always_fail"));
    }

    #[test]
    fn contract_set_admits_short_circuits_on_first_violator() {
        // Two contracts, both have admits returning Err. The set
        // should return the first one's violation.
        struct RejectAt(&'static str);
        impl Contract for RejectAt {
            fn name(&self) -> &str {
                self.0
            }
            fn verify(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
            ) -> Result<(), ContractViolation> {
                Ok(())
            }
            fn admits(
                &self,
                _sim: &Simulation,
                _refdata: Option<&RefDataConfig>,
                _address: &str,
                _candidate: &ChoiceValue,
            ) -> Result<(), ContractViolation> {
                Err(ContractViolation::new(self.name(), "rejected"))
            }
        }

        let s = ContractSet::new()
            .with(Box::new(RejectAt("first")))
            .with(Box::new(RejectAt("second")));

        let sim = Simulation::new();
        let err = s
            .admits(&sim, None, "any.address", &ChoiceValue::Int(0))
            .unwrap_err();
        assert_eq!(err.contract_name, "first");
    }

    #[test]
    fn contract_set_admits_succeeds_when_all_admit() {
        // Default `admits` is Ok, so a set of default contracts
        // admits everything.
        let s = ContractSet::new()
            .with(Box::new(AnchorPreserved::new(Segment::V)))
            .with(Box::new(ProductiveJunctionFrame::new()))
            .with(Box::new(NoStopCodonInJunction::new()));

        let sim = Simulation::new();
        assert!(s
            .admits(&sim, None, "trim.v_3", &ChoiceValue::Int(0))
            .is_ok());
    }

    // ── productive() bundle tests (D.5) ────────────────────────────

    #[test]
    fn productive_bundle_contains_four_contracts() {
        let s = productive();
        assert_eq!(s.len(), 4);
        let names: Vec<&str> = s.iter().map(|c| c.name()).collect();
        assert!(names.contains(&"productive_junction_frame"));
        assert!(names.contains(&"no_stop_codon_in_junction"));
        assert!(names.contains(&"anchor_preserved.v"));
        assert!(names.contains(&"anchor_preserved.j"));
    }

    #[test]
    fn productive_bundle_verifies_clean_sim() {
        // V anchor codon GGG (Gly), J anchor codon TTT (Phe). In-frame.
        let cfg = make_vj_with_anchor_codons(b"GGG", b"TTT");
        let sim = make_assembled_sim_from_refdata(&cfg);

        assert!(productive().verify(&sim, Some(&cfg)).is_ok());
    }

    #[test]
    fn productive_bundle_flags_stop_codon() {
        let cfg = make_vj_with_anchor_codons(b"TAA", b"GGG");
        let sim = make_assembled_sim_from_refdata(&cfg);

        let violations = productive().verify(&sim, Some(&cfg)).unwrap_err();
        let names: Vec<&str> = violations
            .iter()
            .map(|v| v.contract_name.as_str())
            .collect();
        assert!(names.contains(&"no_stop_codon_in_junction"));
    }

    #[test]
    fn productive_bundle_admits_returns_ok_for_default_admits() {
        // None of the productive bundle's contracts override admits
        // yet (D.6 will add it). All defaults → Ok.
        let s = productive();
        let sim = Simulation::new();
        assert!(s
            .admits(&sim, None, "np.np1.length", &ChoiceValue::Int(3))
            .is_ok());
    }
}
