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
use crate::ir::{Nucleotide, Segment, Simulation};
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
