//! Typed action-space supports returned by contracts under the v3.0
//! **constrain-before-propose** API.
//!
//! ## Architectural north star
//!
//! > **Support comes from contracts. Probabilities come from the
//! > natural pass distribution restricted to that support.**
//!
//! Each pass that has migrated to the constrain-before-propose hot
//! path queries the active [`super::ContractSet`] for the narrowed
//! action space at its current decision point, then samples from
//! the natural pass distribution restricted to that support. The
//! pass cannot — by construction — propose a contract-violating
//! action.
//!
//! ## Why split per-kind types instead of one umbrella enum
//!
//! An earlier draft used a single `AdmissibleSet` enum with mixed
//! variants (`BaseMask`, `Positions`, etc) dispatched on an
//! `ActionKind`. That conflated **candidate queries** ("is this
//! specific length admissible?") with **support queries** ("which
//! lengths are admissible?"). The per-kind split below makes the
//! support query the only shape the API can express; candidate
//! queries flow through the pre-existing `admits` predicate. Each
//! support carries the natural shape its consumer wants (bitmask
//! for bases, classified event for indels, length set for trim) —
//! no defensive shape-mixing.
//!
//! ## The hot path versus the pinned-value path
//!
//! [`BaseMask`] expresses the canonical 4-base A/C/G/T support
//! that S5F / uniform / PCR / quality / N-injection-via-A-C-G-T
//! all care about. It is intentionally a 4-bit value so the
//! hot-loop cost is minimal.
//!
//! Non-canonical writes — `N`-injection at lowercase positions,
//! IUPAC-ambiguity flags, etc — don't fit the 4-bit mask. Those
//! go through [`super::Contract::admits_fixed_base_at`] (a yes/no
//! candidate check) instead of the mask. This keeps the hot path
//! lean while keeping pinned-value writes contract-aware.

pub use crate::assignment::TrimEnd;
use crate::ir::Segment;

// ──────────────────────────────────────────────────────────────────
// Base substitution support
// ──────────────────────────────────────────────────────────────────

/// Bitmask over `A/C/G/T` for the set of admissible bases at a pool
/// site, under the constrain-before-propose API.
///
/// - bit 0 (`1`)  = `b'A'`
/// - bit 1 (`2`)  = `b'C'`
/// - bit 2 (`4`)  = `b'G'`
/// - bit 3 (`8`)  = `b'T'`
///
/// Mask of `0b1111` (`15`) is the default Unconstrained — all four
/// bases admissible. Mask of `0b0000` means no base admissible →
/// the pass should skip (permissive) or surface
/// `PassError::ConstraintSampling` (strict).
///
/// **Invariant**: a base `b` is admitted by [`super::Contract::admits`]
/// iff the contract's `admissible_bases_at` returns a mask whose
/// bit for `b` is set. Contract authors maintain this; the engine
/// doesn't enforce it at runtime (regression tests should).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BaseMask(pub u8);

impl BaseMask {
    /// All four bases admissible — the default for contracts with
    /// no opinion on a given site.
    pub const UNCONSTRAINED: Self = BaseMask(0b1111);

    /// No base admissible. Pass should skip / error.
    pub const EMPTY: Self = BaseMask(0b0000);

    /// Bit position (0..4) for the canonical base byte `b`. Returns
    /// `None` for non-A/C/G/T input (e.g. `N`, lowercase) — those
    /// flow through the fixed-base candidate API instead.
    #[inline]
    pub fn bit_for(base: u8) -> Option<u8> {
        match base {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }

    /// Returns the canonical base byte for bit index `i` (0..4).
    pub fn base_for_bit(i: u8) -> u8 {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        BASES[i as usize & 0b11]
    }

    /// Whether base `b` is admissible by this mask. Returns `false`
    /// for non-A/C/G/T values — those should go through
    /// `Contract::admits_fixed_base_at`.
    #[inline]
    pub fn admits(self, base: u8) -> bool {
        match Self::bit_for(base) {
            Some(bit) => (self.0 & (1 << bit)) != 0,
            None => false,
        }
    }

    /// Number of admissible bases (popcount of the mask).
    #[inline]
    pub fn count(self) -> u32 {
        self.0.count_ones()
    }

    /// `true` iff at least one base is admissible.
    #[inline]
    pub fn is_satisfiable(self) -> bool {
        self.0 != 0
    }

    /// Iterate over the canonical base bytes admissible under this
    /// mask, in A/C/G/T order. Used by hot loops that want to
    /// enumerate the admissible support.
    pub fn iter_bases(self) -> impl Iterator<Item = u8> {
        (0u8..4).filter_map(move |bit| {
            if self.0 & (1 << bit) != 0 {
                Some(Self::base_for_bit(bit))
            } else {
                None
            }
        })
    }
}

// ──────────────────────────────────────────────────────────────────
// Indel event class
// ──────────────────────────────────────────────────────────────────

/// The kind of indel candidate being classified — insertion or
/// deletion. Insertions have one extra dimension (the base) the
/// classification can ignore when the contract's verdict is
/// independent of the inserted base (e.g. frame contracts care
/// only about position, not value). Contracts that DO care about
/// the inserted base evaluate it via [`super::Contract::admits_post_event`]
/// in the post-event validator pass.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IndelKindHint {
    Insertion,
    Deletion,
}

/// Per-(site, kind) classification of an indel's impact on the
/// active contract bundle. Returned by
/// [`super::Contract::admissible_indel_class_at`] so the indel pass
/// can run a cross-slot coordination algorithm (mod-3 frame DP
/// for the productive contract) over a typed event vocabulary.
///
/// **The DP only solves the *frame* part of the productive
/// bundle.** A frame-preserving tuple from the DP can still
/// introduce a stop codon in the junction (e.g. an insertion at a
/// junction site with a base that completes a `TAA`). After the
/// pass samples a tuple from the DP-weighted distribution, it
/// must run a final exact validator against the assembled
/// post-event simulation (`ContractSet::admits_post_event`) and
/// retry / no-op / error according to the pass's empty-support
/// policy if the full bundle rejects.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IndelEventClass {
    /// The event has no effect on the active contract — either it
    /// touches a site the contract doesn't constrain (e.g.
    /// out-of-junction position under `ProductiveJunctionFrame`)
    /// or the event is symmetric in the relevant dimension. This
    /// is the identity element of [`Self::compose`].
    FrameNeutral,
    /// The event introduces a net frame delta of `delta` (in
    /// `{-1, +1}`) modulo 3. The pass-level coordinator must
    /// ensure the sum of `delta`s across the selected tuple is
    /// ≡ 0 (mod 3) to satisfy a frame-preservation contract.
    FrameDelta(i8),
    /// The event is forbidden outright by some contract — e.g. an
    /// indel that would delete the V-anchor codon, or an insertion
    /// at a position the active contract rules out unconditionally.
    /// This element absorbs in [`Self::compose`].
    Forbidden,
}

impl IndelEventClass {
    /// Compose two event classes under intersection semantics.
    ///
    /// - `Forbidden` absorbs (any-with-Forbidden = Forbidden).
    /// - `FrameNeutral` is the identity (any-with-Neutral = any).
    /// - Two compatible `FrameDelta` carry through (currently this
    ///   means: two contracts agreeing on the same delta. Internal
    ///   defensive convention: conflicting deltas collapse to
    ///   `Forbidden` so the pass treats the event as inadmissible
    ///   rather than silently picking one delta).
    pub fn compose(self, other: IndelEventClass) -> IndelEventClass {
        use IndelEventClass::*;
        match (self, other) {
            (Forbidden, _) | (_, Forbidden) => Forbidden,
            (FrameNeutral, x) | (x, FrameNeutral) => x,
            (FrameDelta(a), FrameDelta(b)) => {
                if a == b {
                    FrameDelta(a)
                } else {
                    // Two contracts that both classify the same
                    // event with different frame deltas are
                    // architecturally inconsistent — there's no
                    // single coherent action class for the pass.
                    // Treat as Forbidden so the pass can't propose
                    // it.
                    Forbidden
                }
            }
        }
    }

    /// The net frame delta this event contributes, in
    /// `{-1, 0, +1}`. `None` for `Forbidden`.
    pub fn frame_delta(self) -> Option<i8> {
        match self {
            IndelEventClass::FrameNeutral => Some(0),
            IndelEventClass::FrameDelta(d) => Some(d),
            IndelEventClass::Forbidden => None,
        }
    }
}

/// Which segment + end a trim is targeting. Anchor-preserving
/// contracts need both: V-3' trim guards the V anchor codon
/// position; J-5' trim guards the J anchor codon position;
/// D-5'/D-3' guard the D coding region but no anchor; etc.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TrimTarget {
    pub segment: Segment,
    pub end: TrimEnd,
}

// `TrimEnd` is re-exported from `crate::assignment` so the
// admit-set API and the simulation's trim metadata share the same
// type (no conversion churn at the pass / contract boundary).
// The re-export lives at module top via `use crate::assignment::TrimEnd`.

// ──────────────────────────────────────────────────────────────────
// Trim length support
// ──────────────────────────────────────────────────────────────────

/// Admissible trim lengths for an end-loss-style pass that wants
/// to drop `k` bases from one end of the pool. Returned by
/// [`super::Contract::admissible_trim_lengths`] given a max
/// `requested` length the pass would otherwise apply.
///
/// The compact `Full(max)` shape avoids materializing `0..=max`
/// vectors in the common "no contract narrowing" case.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum LengthSupport {
    /// Every length in `0..=max` is admissible (the unconstrained
    /// case). Encoded compactly so the contract bundle doesn't
    /// allocate a full vector when no narrowing happens.
    Full(u32),
    /// Explicit subset of admissible lengths. Always sorted
    /// ascending and unique.
    Subset(Vec<u32>),
    /// No length is admissible — pass should skip / error.
    Empty,
}

impl LengthSupport {
    /// True iff at least one length is admissible.
    pub fn is_satisfiable(&self) -> bool {
        match self {
            LengthSupport::Full(_) => true,
            LengthSupport::Subset(v) => !v.is_empty(),
            LengthSupport::Empty => false,
        }
    }

    /// Maximum admissible length, or `None` if empty.
    pub fn max(&self) -> Option<u32> {
        match self {
            LengthSupport::Full(m) => Some(*m),
            LengthSupport::Subset(v) => v.last().copied(),
            LengthSupport::Empty => None,
        }
    }

    /// Materialize the explicit set of admissible lengths.
    /// `Full(m)` expands to `0..=m`.
    pub fn to_vec(&self) -> Vec<u32> {
        match self {
            LengthSupport::Full(m) => (0..=*m).collect(),
            LengthSupport::Subset(v) => v.clone(),
            LengthSupport::Empty => Vec::new(),
        }
    }

    /// Compose two length supports under intersection semantics.
    /// `Empty` absorbs; `Full(a) & Full(b)` = `Full(min(a, b))`;
    /// otherwise expand to subsets and intersect.
    pub fn intersect(self, other: LengthSupport) -> LengthSupport {
        use LengthSupport::*;
        match (self, other) {
            (Empty, _) | (_, Empty) => Empty,
            (Full(a), Full(b)) => Full(a.min(b)),
            (Full(a), Subset(v)) | (Subset(v), Full(a)) => {
                let filtered: Vec<u32> = v.into_iter().filter(|x| *x <= a).collect();
                if filtered.is_empty() {
                    Empty
                } else {
                    Subset(filtered)
                }
            }
            (Subset(a), Subset(b)) => {
                use std::collections::BTreeSet;
                let set_b: BTreeSet<u32> = b.into_iter().collect();
                let combined: Vec<u32> = a.into_iter().filter(|x| set_b.contains(x)).collect();
                if combined.is_empty() {
                    Empty
                } else {
                    Subset(combined)
                }
            }
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // BaseMask ─────────────────────────────────────────────────

    #[test]
    fn base_mask_unconstrained_admits_all_canonical_bases() {
        let m = BaseMask::UNCONSTRAINED;
        for b in [b'A', b'C', b'G', b'T'] {
            assert!(m.admits(b), "unconstrained should admit {}", b as char);
        }
    }

    #[test]
    fn base_mask_empty_admits_nothing() {
        let m = BaseMask::EMPTY;
        for b in [b'A', b'C', b'G', b'T'] {
            assert!(!m.admits(b));
        }
        assert_eq!(m.count(), 0);
        assert!(!m.is_satisfiable());
    }

    #[test]
    fn base_mask_rejects_non_canonical_bytes() {
        let m = BaseMask::UNCONSTRAINED;
        assert!(!m.admits(b'N'));
        assert!(!m.admits(b'X'));
        assert!(!m.admits(0));
    }

    #[test]
    fn base_mask_case_insensitive() {
        let m = BaseMask(0b0001); // only A
        assert!(m.admits(b'A'));
        assert!(m.admits(b'a'));
        assert!(!m.admits(b'C'));
    }

    #[test]
    fn base_mask_iter_bases_in_acgt_order() {
        let m = BaseMask(0b1111);
        let bases: Vec<u8> = m.iter_bases().collect();
        assert_eq!(bases, vec![b'A', b'C', b'G', b'T']);

        let m = BaseMask(0b1010); // C, T
        let bases: Vec<u8> = m.iter_bases().collect();
        assert_eq!(bases, vec![b'C', b'T']);
    }

    // IndelEventClass ───────────────────────────────────────────

    #[test]
    fn indel_event_class_forbidden_absorbs() {
        use IndelEventClass::*;
        for other in [FrameNeutral, FrameDelta(1), FrameDelta(-1), Forbidden] {
            assert_eq!(Forbidden.compose(other), Forbidden);
            assert_eq!(other.compose(Forbidden), Forbidden);
        }
    }

    #[test]
    fn indel_event_class_neutral_is_identity() {
        use IndelEventClass::*;
        assert_eq!(FrameNeutral.compose(FrameNeutral), FrameNeutral);
        assert_eq!(FrameNeutral.compose(FrameDelta(1)), FrameDelta(1));
        assert_eq!(FrameDelta(-1).compose(FrameNeutral), FrameDelta(-1));
    }

    #[test]
    fn indel_event_class_matching_deltas_compose() {
        use IndelEventClass::*;
        assert_eq!(FrameDelta(1).compose(FrameDelta(1)), FrameDelta(1));
        assert_eq!(FrameDelta(-1).compose(FrameDelta(-1)), FrameDelta(-1));
    }

    #[test]
    fn indel_event_class_conflicting_deltas_become_forbidden() {
        use IndelEventClass::*;
        // Two contracts disagreeing on the frame delta of an event
        // is architecturally inconsistent — we collapse to
        // Forbidden defensively rather than silently picking one.
        assert_eq!(FrameDelta(1).compose(FrameDelta(-1)), Forbidden);
    }

    #[test]
    fn indel_event_class_frame_delta_accessor() {
        use IndelEventClass::*;
        assert_eq!(FrameNeutral.frame_delta(), Some(0));
        assert_eq!(FrameDelta(1).frame_delta(), Some(1));
        assert_eq!(FrameDelta(-1).frame_delta(), Some(-1));
        assert_eq!(Forbidden.frame_delta(), None);
    }

    // LengthSupport ─────────────────────────────────────────────

    #[test]
    fn length_support_full_satisfiable() {
        assert!(LengthSupport::Full(8).is_satisfiable());
        assert!(LengthSupport::Full(0).is_satisfiable()); // {0}
    }

    #[test]
    fn length_support_intersect_full_full_keeps_min() {
        let s = LengthSupport::Full(8).intersect(LengthSupport::Full(3));
        assert_eq!(s, LengthSupport::Full(3));
    }

    #[test]
    fn length_support_intersect_subset_full_filters() {
        let s = LengthSupport::Subset(vec![0, 3, 6, 9]).intersect(LengthSupport::Full(6));
        assert_eq!(s, LengthSupport::Subset(vec![0, 3, 6]));
    }

    #[test]
    fn length_support_intersect_disjoint_subsets_is_empty() {
        let s =
            LengthSupport::Subset(vec![0, 3, 6]).intersect(LengthSupport::Subset(vec![1, 4, 7]));
        assert_eq!(s, LengthSupport::Empty);
    }

    #[test]
    fn length_support_full_intersect_empty_is_empty() {
        assert_eq!(
            LengthSupport::Full(5).intersect(LengthSupport::Empty),
            LengthSupport::Empty
        );
    }

    #[test]
    fn length_support_to_vec_expands_full() {
        assert_eq!(LengthSupport::Full(3).to_vec(), vec![0, 1, 2, 3]);
        assert_eq!(LengthSupport::Subset(vec![0, 5]).to_vec(), vec![0, 5]);
        assert_eq!(LengthSupport::Empty.to_vec(), Vec::<u32>::new());
    }
}
