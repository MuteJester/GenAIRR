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
use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
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
/// [`matches_observed`] semantics.
pub fn alleles_compatible_at(alleles: &[Allele], ref_pos: u32, observed: u8) -> Vec<bool> {
    alleles
        .iter()
        .map(|a| {
            a.seq
                .get(ref_pos as usize)
                .map(|&ref_byte| matches_observed(observed, ref_byte))
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
        // Backwards ref_pos motion is unsupported (walker.rs:81).
        if let Some(expected) = next_ref_pos {
            if ref_pos < expected {
                return scores;
            }
        }
        next_ref_pos = Some(ref_pos.saturating_add(1));
        if ref_start.is_none() {
            ref_start = Some(ref_pos);
        }
        // Per-position score increment via the shared kernel rules.
        for (allele_idx, allele) in alleles.iter().enumerate() {
            if let Some(&ref_byte) = allele.seq.get(ref_pos as usize) {
                if matches_observed(nuc.base, ref_byte) {
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

    // ── Left extension ──────────────────────────────────────────
    // Mirrors walker/extensions.rs::walk_left_extension. Look for an
    // NP region whose end == current_seq_start; walk backwards
    // through it (and into the previous V/D/J region in the overlap
    // case), capped by trim_5_cap, extending only when the byte
    // strictly narrows the current max-score tie set.
    if let Some(np_region) = find_left_extension(sim, current_seq_start) {
        // Lower bound: the start of the V/D/J region whose end
        // touches np_region.start (overlap-into-neighbor case), or 0.
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
            if steps_taken >= trim_5_cap {
                break;
            }
            if current_ref_start == 0 {
                break;
            }
            let candidate_ref_pos = current_ref_start - 1;
            let Some(nuc) = pool.get(seq_pos as usize) else {
                break;
            };
            let matched = alleles_compatible_at(alleles, candidate_ref_pos, nuc.base);
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
            current_ref_start = candidate_ref_pos;
            steps_taken = steps_taken.saturating_add(1);
        }
    }

    // ── Right extension ─────────────────────────────────────────
    // Mirrors walker/extensions.rs::walk_right_extension.
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
            if steps_taken >= trim_3_cap {
                break;
            }
            let Some(nuc) = pool.get(seq_pos as usize) else {
                break;
            };
            let matched = alleles_compatible_at(alleles, current_ref_end, nuc.base);
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
            current_ref_end = current_ref_end.saturating_add(1);
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
