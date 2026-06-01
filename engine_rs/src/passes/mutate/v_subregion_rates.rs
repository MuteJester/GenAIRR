//! `VSubregionRateWeights` — per-V-subregion SHM rate scalars.
//!
//! Slice B of the V-region substructure work
//! (`docs/v_subregion_shm_rate_design.md`). Layers per-V-subregion
//! rate scalars on top of the existing per-segment rates: a V site's
//! final weight is `base × segment_rate(V) × v_subregion_rate(label)`.
//! Non-V sites ignore the subregion rate entirely. V sites belonging
//! to a V allele without IMGT subregion annotation receive the
//! identity factor (1.0) — supports mixed cartridges.
//!
//! See the audit doc §1 (rate model shape) and §2 (pool position →
//! V-allele subregion mapping feasibility) for the architectural
//! rationale; the helper [`v_subregion_at_position`] is the audit's
//! sibling of [`super::segment_rates::segment_at_position`].

use crate::assignment::AlleleAssignments;
use crate::ir::{Segment, Sequence};
use crate::refdata::{RefDataConfig, VSubregionLabel};

/// Per-V-subregion SHM rate scalars in canonical IMGT order.
///
/// Fields are laid out FWR1 / CDR1 / FWR2 / CDR2 / FWR3 so positional
/// plumbing through PyO3 stays straightforward.
///
/// Construction invariants the caller must honour:
/// - Every field is finite and `>= 0.0`.
/// - At least one field is strictly positive (an all-zero vector is
///   unsatisfiable — every V site would have zero proposal weight).
///
/// The Python DSL boundary (`experiment._validate_v_subregion_rates`)
/// is the only caller that constructs non-default values; it
/// enforces these invariants before the tuple crosses the bridge.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct VSubregionRateWeights {
    pub fwr1: f64,
    pub cdr1: f64,
    pub fwr2: f64,
    pub cdr2: f64,
    pub fwr3: f64,
}

impl VSubregionRateWeights {
    /// Flat default — all five labels at `1.0`. Produces byte-
    /// identical sampling distributions to the pre-slice engine
    /// (both uniform and S5F passes detect this via
    /// [`Self::is_default`] and take the existing fast path).
    pub const fn default() -> Self {
        Self {
            fwr1: 1.0,
            cdr1: 1.0,
            fwr2: 1.0,
            cdr2: 1.0,
            fwr3: 1.0,
        }
    }

    /// Construct from the canonical DSL tuple ordering. Caller
    /// must have validated the inputs (no NaN/inf/negative); this
    /// constructor is a thin packaging step.
    pub const fn from_tuple(fwr1: f64, cdr1: f64, fwr2: f64, cdr2: f64, fwr3: f64) -> Self {
        Self {
            fwr1,
            cdr1,
            fwr2,
            cdr2,
            fwr3,
        }
    }

    /// Whether this rate vector matches the flat default. Used by
    /// both mutation passes to short-circuit the per-position
    /// subregion lookup and preserve the pre-slice fast path.
    #[inline]
    pub fn is_default(&self) -> bool {
        self.fwr1 == 1.0
            && self.cdr1 == 1.0
            && self.fwr2 == 1.0
            && self.cdr2 == 1.0
            && self.fwr3 == 1.0
    }

    /// Per-label rate lookup.
    #[inline]
    pub fn rate_for(&self, label: VSubregionLabel) -> f64 {
        match label {
            VSubregionLabel::Fwr1 => self.fwr1,
            VSubregionLabel::Cdr1 => self.cdr1,
            VSubregionLabel::Fwr2 => self.fwr2,
            VSubregionLabel::Cdr2 => self.cdr2,
            VSubregionLabel::Fwr3 => self.fwr3,
        }
    }
}

impl Default for VSubregionRateWeights {
    fn default() -> Self {
        Self::default()
    }
}

/// Locate the [`VSubregionLabel`] that a pool position `pos` falls
/// into on the **assigned V allele's** ungapped sequence. Returns:
///
/// - `None` if `pos` is outside the V region (segment is D/J/NP/C
///   or pre-assembly).
/// - `None` if the V region exists but no V allele is assigned (a
///   pre-sampling state — shouldn't happen during the mutate pass
///   but coded defensively).
/// - `None` if the assigned V allele has no subregion annotations
///   (legacy / user-authored V allele without `gapped_seq`). The
///   caller's policy is to apply factor `1.0` in that case.
/// - `None` if the allele-local position falls outside every
///   `[start, end)` interval — covers the V-side CDR3 contribution
///   (positions between `FWR3.end` and `len(allele.seq)`), which
///   the five-label annotation set deliberately leaves
///   unannotated.
/// - `Some(label)` for a position inside one of FWR1 / CDR1 /
///   FWR2 / CDR2 / FWR3.
///
/// Allele-local arithmetic (audit §2):
///
/// ```text
/// allele_local = (pos - v_region.start) + instance.trim_5
/// ```
///
/// Linear over `Allele.subregions.len()` (always ≤ 5) — same cost
/// shape as [`super::segment_rates::segment_at_position`]'s linear
/// scan over `sequence.regions` (always ≤ 5 regions in a VDJ pool).
#[inline]
pub fn v_subregion_at_position(
    sequence: &Sequence,
    assignments: &AlleleAssignments,
    refdata: &RefDataConfig,
    pos: u32,
) -> Option<VSubregionLabel> {
    // 1. Find the V region in the assembled pool.
    let v_region = sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V)?;
    // Defensive: ensure pos actually falls inside that region.
    if pos < v_region.start.index() || pos >= v_region.end.index() {
        return None;
    }
    // 2. Get the assigned V allele instance.
    let instance = assignments.get(Segment::V)?;
    let allele = refdata.v_pool.get(instance.allele_id)?;
    // 3. Compute allele-local offset.
    let pool_offset_in_v = pos - v_region.start.index();
    let allele_local: u32 = pool_offset_in_v + instance.trim_5 as u32;
    if allele_local >= allele.seq.len() as u32 {
        // Defence-in-depth: shouldn't happen for well-formed
        // assignments, but a bug upstream shouldn't crash sampling.
        return None;
    }
    // 4. Walk subregions for an interval containing allele_local.
    for sub in &allele.subregions {
        if (sub.start as u32) <= allele_local && allele_local < (sub.end as u32) {
            return Some(sub.label);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_all_ones() {
        let w = VSubregionRateWeights::default();
        assert!(w.is_default());
        assert_eq!(w.fwr1, 1.0);
        assert_eq!(w.cdr1, 1.0);
        assert_eq!(w.fwr2, 1.0);
        assert_eq!(w.cdr2, 1.0);
        assert_eq!(w.fwr3, 1.0);
    }

    #[test]
    fn rate_for_routes_each_label() {
        let w = VSubregionRateWeights::from_tuple(0.5, 2.0, 1.5, 3.0, 0.25);
        assert_eq!(w.rate_for(VSubregionLabel::Fwr1), 0.5);
        assert_eq!(w.rate_for(VSubregionLabel::Cdr1), 2.0);
        assert_eq!(w.rate_for(VSubregionLabel::Fwr2), 1.5);
        assert_eq!(w.rate_for(VSubregionLabel::Cdr2), 3.0);
        assert_eq!(w.rate_for(VSubregionLabel::Fwr3), 0.25);
        assert!(!w.is_default());
    }

    #[test]
    fn one_field_off_default_is_not_default() {
        let w = VSubregionRateWeights::from_tuple(1.0, 1.0, 1.0, 1.0, 0.5);
        assert!(!w.is_default());
    }
}
