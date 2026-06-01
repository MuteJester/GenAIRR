//! Shared allele-call scoring kernel.
//!
//! Single source of truth for the engine's allele match semantics
//! and tie-set construction. Both the live-call walker and the
//! AIRR-record validator route their match decisions through here
//! so the two can never silently diverge.
//!
//! What this owns:
//!
//! - `classify_base` — A/C/G/T (case-insensitive) → `Canonical(upper)`,
//!   N/n → `Wildcard`, anything else → `Invalid`.
//! - `matches_observed` — does an observed pool byte match an allele's
//!   germline byte at the same position? `N` is wildcard; lowercase
//!   is upper-cased before comparison (NOT treated as wildcard);
//!   non-canonical observed bases never match.
//! - `score_alleles_in_region` — region-level scorer. Walks the pool
//!   region, counts matches per allele against each allele's germline
//!   byte at the matching `germline_pos`. Returns per-allele match
//!   counts.
//! - `tie_set_at_max_score` — select alleles tied at the maximum
//!   match count. Returns an `AlleleBitSet`. When max score is 0
//!   (no germline evidence) the caller chooses the semantics:
//!   walker returns the full allele universe (every allele equally
//!   consistent with absence of evidence); validator skips the
//!   check (no rescoring oracle).
//!
//! What this deliberately does NOT own:
//!
//! - Truth-allele favoritism. The scoring kernel is biology-only:
//!   it ranks by evidence with no knowledge of which allele was
//!   sampled. Truth-first CSV ordering is a projection convention
//!   (`airr_record::projection::live_call_name`), not scoring.
//! - The walker's inverted index. The walker uses `ReferenceMatchIndex`
//!   for O(1) per-byte lookups; that's a performance structure built
//!   on top of these primitives, not part of the scoring model.

use super::bitset::AlleleBitSet;
use crate::assignment::SegmentOrientation;
use crate::ir::{complement_base, NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::refdata::{Allele, AlleleId, RefDataConfig};

/// Classification of a single observed byte.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum BaseKind {
    /// One of A/C/G/T (upper-cased). The contained byte is the
    /// canonical (upper-case) form.
    Canonical(u8),
    /// `N` or `n`. Matches any canonical germline base.
    Wildcard,
    /// Anything else — never matches.
    Invalid,
}

/// Classify a byte under the engine's match semantics. Lowercase
/// `a/c/g/t` is upper-cased (NOT treated as wildcard — the
/// `sequencing_errors` lowercase convention is presentational, not
/// evidential). Only literal `N`/`n` is wildcard.
#[inline]
pub fn classify_base(b: u8) -> BaseKind {
    match b.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' => BaseKind::Canonical(b.to_ascii_uppercase()),
        b'N' => BaseKind::Wildcard,
        _ => BaseKind::Invalid,
    }
}

/// Transform an observed pool byte into the byte the orientation-
/// agnostic match rule should compare against the allele's germline.
///
/// For [`SegmentOrientation::Forward`]: identity — the observed byte
/// is already in the same orientation as the allele's reference.
///
/// For [`SegmentOrientation::ReverseComplement`]: pre-complement.
/// The assembly pass emits inverted-D bytes as
/// `complement(allele[allele_pos])` while keeping the
/// `germline_pos` pointing at the original allele coordinate.
/// Complementing the observed byte once at the comparison boundary
/// recovers the original allele byte, so the existing match rule
/// (and the inverted index keyed on canonical germline bytes) works
/// against inverted segments without any duplicate data structure.
///
/// `N` / `n` / invalid bytes pass through unchanged: complementing
/// them does not produce a more informative value, and the
/// downstream `classify_base` rule already handles wildcards.
#[inline]
pub fn observed_in_germline_orientation(observed: u8, orientation: SegmentOrientation) -> u8 {
    match orientation {
        SegmentOrientation::Forward => observed,
        SegmentOrientation::ReverseComplement => complement_base(observed),
    }
}

/// Orientation-aware variant of [`matches_observed`].
///
/// For [`SegmentOrientation::Forward`] the behaviour is identical
/// to [`matches_observed`].
///
/// For [`SegmentOrientation::ReverseComplement`] the comparison is
/// performed against `complement_base(observed)` so the rule
/// matches inverted-D evidence: the assembled pool byte at
/// `germline_pos = allele_pos` is the Watson-Crick complement of
/// `allele[allele_pos]`, so complementing the observed byte once
/// recovers the original allele coordinate.
///
/// Used by the walker observer and the AIRR validator's allele-call
/// oracle so both surfaces score inverted D evidence consistently.
#[inline]
pub fn matches_observed_with_orientation(
    observed: u8,
    germline: u8,
    orientation: SegmentOrientation,
) -> bool {
    matches_observed(observed_in_germline_orientation(observed, orientation), germline)
}

/// Does the observed pool byte match the allele's germline byte at
/// the same position?
///
/// Match rules:
/// - Both bytes are classified via [`classify_base`].
/// - `Wildcard` observed matches any `Canonical` germline.
/// - `Canonical` observed matches a `Canonical` germline iff the
///   canonical bytes are equal (case-insensitive).
/// - `Invalid` either side → no match.
/// - `Wildcard` germline (allele containing literal `N` in its
///   reference sequence) is an edge case the engine doesn't model
///   today; treated as wildcard for symmetry.
#[inline]
pub fn matches_observed(observed: u8, germline: u8) -> bool {
    match (classify_base(observed), classify_base(germline)) {
        (BaseKind::Wildcard, BaseKind::Canonical(_)) => true,
        (BaseKind::Canonical(_), BaseKind::Wildcard) => true,
        (BaseKind::Wildcard, BaseKind::Wildcard) => true,
        (BaseKind::Canonical(o), BaseKind::Canonical(g)) => o == g,
        _ => false,
    }
}

/// Per-allele match count over a region's pool bytes.
///
/// For each pool nucleotide at index `[region_start, region_end)`:
/// - if the nucleotide is an indel-insert (`germline_pos == None`),
///   it carries no allele-call evidence and is skipped;
/// - otherwise, the byte is compared (via [`matches_observed`])
///   against each allele's reference byte at the same
///   `germline_pos`, incrementing that allele's score when matching.
///
/// Returns one count per allele in `alleles`, in the same order.
///
/// This is the pure batch form of the per-byte scoring loop the
/// walker observer runs incrementally. The walker uses a
/// precomputed inverted index for O(1) per-byte lookups; this
/// function uses the straight `O(region_len * alleles_len)` scan.
/// Both produce the same scores by construction.
pub fn score_alleles_in_region(
    pool: &[Nucleotide],
    region_start: u32,
    region_end: u32,
    alleles: &[Allele],
) -> Vec<u32> {
    let mut scores = vec![0u32; alleles.len()];
    for nuc_idx in region_start..region_end {
        let Some(nuc) = pool.get(nuc_idx as usize) else {
            break;
        };
        let Some(gp) = nuc.germline_pos.get() else {
            continue;
        };
        for (allele_idx, allele) in alleles.iter().enumerate() {
            if let Some(&ref_byte) = allele.seq.get(gp as usize) {
                if matches_observed(nuc.base, ref_byte) {
                    scores[allele_idx] += 1;
                }
            }
        }
    }
    scores
}

/// Build the tie-set bitset of alleles at the maximum score.
///
/// If `max(scores) == 0` returns `None` — the caller decides what
/// "no evidence" means (the walker returns the full universe; the
/// validator skips the oracle check).
///
/// `universe_len` is the bitset capacity (typically `alleles.len()`).
/// Scores indices map to `AlleleId(i)`.
#[allow(dead_code)] // walker uses its own tie-set construction inline; this is the public form
pub fn tie_set_at_max_score(scores: &[u32], universe_len: usize) -> Option<AlleleBitSet> {
    let max_score = *scores.iter().max().unwrap_or(&0);
    if max_score == 0 {
        return None;
    }
    let mut set = AlleleBitSet::empty(universe_len);
    for (idx, &s) in scores.iter().enumerate() {
        if s == max_score {
            set.insert(AlleleId::new(idx as u32));
        }
    }
    Some(set)
}

/// Helper: list the AlleleIds at max score as a sorted Vec. Used by
/// the validator's oracle which works in vector form. Returns an
/// empty vec when no evidence exists.
pub fn tie_set_ids_at_max_score(scores: &[u32]) -> Vec<AlleleId> {
    let max_score = *scores.iter().max().unwrap_or(&0);
    if max_score == 0 {
        return Vec::new();
    }
    scores
        .iter()
        .enumerate()
        .filter(|(_, &s)| s == max_score)
        .map(|(i, _)| AlleleId::new(i as u32))
        .collect()
}

/// Allele pool lookup helper. Centralises the `match segment { V =>
/// refdata.v_pool, ... }` table so callers don't drift on which
/// segments are scoreable.
pub fn allele_pool_for_segment<'a>(
    refdata: &'a RefDataConfig,
    segment: Segment,
) -> Option<&'a [Allele]> {
    match segment {
        Segment::V => Some(refdata.v_pool.as_slice()),
        Segment::D => Some(refdata.d_pool.as_slice()),
        Segment::J => Some(refdata.j_pool.as_slice()),
        _ => None,
    }
}

// ──────────────────────────────────────────────────────────────────
// NP-extension-aware scoring (mirrors walker/extensions.rs)
//
// The live-call walker doesn't just score the structural V/D/J
// region — it walks BACKWARDS into the preceding NP region (capped
// by the assigned allele's `trim_5`) and FORWARDS into the
// following NP region (capped by `trim_3`), counting matches against
// each allele's reference at the extended ref positions. Each
// extension step is gated by `extension_narrows_tie_set`: extend
// only when the new byte strictly narrows the current max-score
// tie set (so NP bytes can discriminate but never widen).
//
// The validator's allele-call oracle (C4) must match the walker
// to avoid spurious tie-set mismatches under non-zero trim. Both
// paths now share the rules below.
// ──────────────────────────────────────────────────────────────────

/// Returns a boolean mask (one entry per allele in `alleles`)
/// indicating which alleles' reference byte at `ref_pos` matches
/// the `observed` pool byte under the shared
/// [`matches_observed_with_orientation`] semantics.
///
/// `orientation` is `Forward` for ordinary scoring; inverted-D
/// callers pass `ReverseComplement` so the observed byte is
/// pre-complemented before the comparison (see
/// [`matches_observed_with_orientation`] for the rationale).
pub fn alleles_compatible_at(
    alleles: &[Allele],
    ref_pos: u32,
    observed: u8,
    orientation: SegmentOrientation,
) -> Vec<bool> {
    alleles
        .iter()
        .map(|a| {
            a.seq
                .get(ref_pos as usize)
                .map(|&ref_byte| {
                    matches_observed_with_orientation(observed, ref_byte, orientation)
                })
                .unwrap_or(false)
        })
        .collect()
}

/// Mirror of `walker/extensions.rs::extension_narrows_tie_set`.
///
/// Returns `true` iff incrementing scores for the alleles in
/// `matched` would strictly narrow the current max-score tie set.
/// Three conditions must all hold:
///
/// 1. `pre_max > 0` — extensions don't operate from an empty
///    evidence floor (the structural walk seeds the scores).
/// 2. At least one allele tied at `pre_max` is in `matched` —
///    so a new `pre_max + 1` plateau exists after the increment.
/// 3. At least one allele tied at `pre_max` is NOT in `matched` —
///    so the new plateau excludes them (strict narrowing).
pub fn extension_narrows_tie_set(scores: &[u32], matched: &[bool]) -> bool {
    let pre_max = scores.iter().copied().max().unwrap_or(0);
    if pre_max == 0 {
        return false;
    }
    let mut any_matched = false;
    let mut any_missed = false;
    for (i, &score) in scores.iter().enumerate() {
        if score != pre_max {
            continue;
        }
        if matched.get(i).copied().unwrap_or(false) {
            any_matched = true;
        } else {
            any_missed = true;
        }
        if any_matched && any_missed {
            return true;
        }
    }
    false
}

/// Score every allele against a segment's structural region plus
/// the walker's NP-region extensions. Mirrors
/// `live_call::walker::call_from_region` (structural walk +
/// `walk_left_extension` + `walk_right_extension`) but uses
/// per-byte allele scans instead of the walker's inverted
/// `ReferenceMatchIndex`. The per-allele match counts match the
/// walker by construction (same match semantics via
/// [`matches_observed`], same extension gating via
/// [`extension_narrows_tie_set`], same trim caps).
///
/// Used by the validator's C4 allele-call oracle so the oracle
/// agrees with the walker under arbitrary trim_5/trim_3 caps.
///
/// Returns the per-allele score vector. The caller selects the
/// tie-set via [`tie_set_ids_at_max_score`].
pub fn score_alleles_with_extensions(
    sim: &Simulation,
    segment: Segment,
    alleles: &[Allele],
    structural_region: &Region,
    trim_5_cap: u32,
    trim_3_cap: u32,
    orientation: SegmentOrientation,
) -> Vec<u32> {
    let pool = sim.pool.as_slice();
    let mut scores = vec![0u32; alleles.len()];

    // ── Structural walk ─────────────────────────────────────────
    // Iterate the structural region's pool bytes, score each allele
    // against its reference at the byte's germline_pos. Indel-
    // inserted bytes (germline_pos == None) carry no evidence and
    // are skipped, matching walker.rs:72-74.
    let mut ref_start: Option<u32> = None;
    let mut next_ref_pos: Option<u32> = None;
    for seq_pos in structural_region.start.index()..structural_region.end.index() {
        let Some(nuc) = pool.get(seq_pos as usize) else {
            break;
        };
        if nuc.segment != segment {
            // walker.rs:64 returns unsupported_call here; for the
            // validator's oracle we abort scoring (no tie set will
            // be meaningful against a malformed region).
            return scores;
        }
        let Some(ref_pos) = nuc.germline_pos.get().map(|p| p as u32) else {
            continue;
        };
        // Orientation-aware motion. Under `Forward` we require
        // strictly-monotonic-increasing `germline_pos` (gap-up
        // permitted for deletions). Under `ReverseComplement`
        // the assembly emits bytes in reverse allele order so we
        // mirror the malformed check: backwards = increasing.
        // `ref_start` tracks the smallest seen position;
        // `next_ref_pos` tracks `max_seen + 1`. Same `[start,
        // end)` semantics as the walker; the iteration order is
        // the only orientation-dependent piece.
        let is_rc = matches!(
            orientation,
            SegmentOrientation::ReverseComplement
        );
        if let Some(expected) = next_ref_pos {
            let backwards = if is_rc {
                ref_pos > expected
            } else {
                ref_pos < expected
            };
            if backwards {
                return scores;
            }
        }
        ref_start = Some(match ref_start {
            Some(prev) => prev.min(ref_pos),
            None => ref_pos,
        });
        let new_end_candidate = ref_pos.saturating_add(1);
        next_ref_pos = Some(match next_ref_pos {
            Some(prev) => prev.max(new_end_candidate),
            None => new_end_candidate,
        });
        // Per-position score increment via the shared kernel rules.
        // Orientation pre-complements the observed byte when the
        // segment is in `ReverseComplement` (inverted-D path), so
        // the comparison sees the original allele coordinate.
        for (allele_idx, allele) in alleles.iter().enumerate() {
            if let Some(&ref_byte) = allele.seq.get(ref_pos as usize) {
                if matches_observed_with_orientation(nuc.base, ref_byte, orientation) {
                    scores[allele_idx] += 1;
                }
            }
        }
    }

    // No structural evidence (every byte was an indel-insert)
    // means no ref window to extend from — return the zero scores.
    let (Some(mut current_ref_start), Some(mut current_ref_end)) = (ref_start, next_ref_pos)
    else {
        return scores;
    };
    let mut current_seq_start = structural_region.start.index();
    let mut current_seq_end = structural_region.end.index();

    // Trim caps swap by allele direction under RC (audit §3).
    // The pool-direction iteration stays identical; what changes
    // is which cap bounds which walk and which `ref_*` boundary
    // moves on accept.
    let is_rc = matches!(orientation, SegmentOrientation::ReverseComplement);
    let (left_cap, right_cap) = if is_rc {
        (trim_3_cap, trim_5_cap)
    } else {
        (trim_5_cap, trim_3_cap)
    };

    // ── Left extension ──────────────────────────────────────────
    // Mirrors walker/extensions.rs::walk_left_extension. Look for an
    // NP region whose end == current_seq_start; walk backwards
    // through it, capped by `left_cap`. Under Forward, the
    // candidate allele position is `current_ref_start - 1` and an
    // accepted step decreases `current_ref_start`. Under RC, the
    // candidate is `current_ref_end` and an accepted step increases
    // `current_ref_end`.
    if let Some(np_region) = find_left_extension(sim, current_seq_start) {
        let lower_bound = sim
            .sequence
            .regions
            .iter()
            .find(|r| {
                matches!(r.segment, Segment::V | Segment::D | Segment::J)
                    && r.end == np_region.start
            })
            .map(|r| r.start.index())
            .unwrap_or(0);

        let mut steps_taken: u32 = 0;
        for seq_pos in (lower_bound..np_region.end.index()).rev() {
            if steps_taken >= left_cap {
                break;
            }
            let candidate_ref_pos = if is_rc {
                current_ref_end
            } else {
                if current_ref_start == 0 {
                    break;
                }
                current_ref_start - 1
            };
            let Some(nuc) = pool.get(seq_pos as usize) else {
                break;
            };
            let matched =
                alleles_compatible_at(alleles, candidate_ref_pos, nuc.base, orientation);
            if !matched.iter().any(|&m| m) {
                break;
            }
            if !extension_narrows_tie_set(&scores, &matched) {
                break;
            }
            for (idx, &m) in matched.iter().enumerate() {
                if m {
                    scores[idx] = scores[idx].saturating_add(1);
                }
            }
            current_seq_start = seq_pos;
            if is_rc {
                current_ref_end = candidate_ref_pos.saturating_add(1);
            } else {
                current_ref_start = candidate_ref_pos;
            }
            steps_taken = steps_taken.saturating_add(1);
        }
    }

    // ── Right extension ─────────────────────────────────────────
    // Mirrors walker/extensions.rs::walk_right_extension. Under
    // Forward the candidate is `current_ref_end`; under RC it's
    // `current_ref_start - 1`.
    if let Some(np_region) = find_right_extension(sim, current_seq_end) {
        let pool_len = pool.len() as u32;
        let upper_bound = sim
            .sequence
            .regions
            .iter()
            .find(|r| {
                matches!(r.segment, Segment::V | Segment::D | Segment::J)
                    && r.start == np_region.end
            })
            .map(|r| r.end.index())
            .unwrap_or(pool_len)
            .min(pool_len);

        let mut steps_taken: u32 = 0;
        for seq_pos in np_region.start.index()..upper_bound {
            if steps_taken >= right_cap {
                break;
            }
            let candidate_ref_pos = if is_rc {
                if current_ref_start == 0 {
                    break;
                }
                current_ref_start - 1
            } else {
                current_ref_end
            };
            let Some(nuc) = pool.get(seq_pos as usize) else {
                break;
            };
            let matched =
                alleles_compatible_at(alleles, candidate_ref_pos, nuc.base, orientation);
            if !matched.iter().any(|&m| m) {
                break;
            }
            if !extension_narrows_tie_set(&scores, &matched) {
                break;
            }
            for (idx, &m) in matched.iter().enumerate() {
                if m {
                    scores[idx] = scores[idx].saturating_add(1);
                }
            }
            current_seq_end = seq_pos.saturating_add(1);
            if is_rc {
                current_ref_start = candidate_ref_pos;
            } else {
                current_ref_end = candidate_ref_pos.saturating_add(1);
            }
            steps_taken = steps_taken.saturating_add(1);
        }
    }

    let _ = (current_seq_start, current_seq_end, NucHandle::new(0));
    scores
}

/// Find an NP region whose end touches `seq_start` (i.e. the NP is
/// the immediate left neighbor of a structural region starting at
/// `seq_start`).
fn find_left_extension(sim: &Simulation, seq_start: u32) -> Option<Region> {
    sim.sequence
        .regions
        .iter()
        .find(|r| {
            matches!(r.segment, Segment::Np1 | Segment::Np2) && r.end.index() == seq_start
        })
        .cloned()
}

/// Find an NP region whose start touches `seq_end` (i.e. the NP is
/// the immediate right neighbor of a structural region ending at
/// `seq_end`).
fn find_right_extension(sim: &Simulation, seq_end: u32) -> Option<Region> {
    sim.sequence
        .regions
        .iter()
        .find(|r| {
            matches!(r.segment, Segment::Np1 | Segment::Np2) && r.start.index() == seq_end
        })
        .cloned()
}

/// Convenience: rescore every allele in a segment's pool against
/// the segment's region in `sim`. Returns the per-allele scores
/// or `None` when the segment has no assembled region.
#[allow(dead_code)] // reserved for future validator/walker callers
pub fn score_segment_in_sim<'a>(
    sim: &Simulation,
    refdata: &'a RefDataConfig,
    segment: Segment,
) -> Option<(Vec<u32>, &'a [Allele])> {
    let region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == segment)?;
    let alleles = allele_pool_for_segment(refdata, segment)?;
    if alleles.is_empty() {
        return None;
    }
    let scores = score_alleles_in_region(
        sim.pool.as_slice(),
        region.start.index(),
        region.end.index(),
        alleles,
    );
    Some((scores, alleles))
}

#[cfg(test)]
mod tests {
    //! Tests live alongside the rest of the live_call tests at
    //! `engine_rs/src/live_call/tests/scoring.rs`.
}
