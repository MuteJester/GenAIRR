use crate::ir::{Region, Segment, Simulation};
use crate::live_call::SegmentLiveCall;
use crate::refdata::{Allele, RefDataConfig};

fn lookup_name(
    refdata: &RefDataConfig,
    segment: Segment,
    id: Option<crate::refdata::AlleleId>,
) -> String {
    id.and_then(|aid| refdata.get(segment, aid))
        .map(|a| a.name.clone())
        .unwrap_or_default()
}

/// pick the canonical allele id for AIRR projection on
/// a V/D/J segment. The first allele in the live-call's narrowed
/// `allele_call` set wins; fall back to the originally-sampled
/// allele (provenance) when the live layer is absent or the call
/// is `Unsupported` (empty candidate set).
///
/// This is the seam through which `germline_alignment` becomes
/// evidence-driven: when SHM mutations narrow the live call to a
/// different allele than the one originally sampled, the column
/// walker emits *that* allele's bytes (matching the `*_call` field
/// the record reports) instead of the historical provenance bytes.
pub(super) fn projected_allele_id(
    sim: &Simulation,
    segment: Segment,
) -> Option<crate::refdata::AlleleId> {
    // The sampled allele from recombination - used both as the
    // fallback when no live-call exists and as a preference among
    // alleles tied in the live-call's score-and-tie set. Several
    // alleles may share the same per-position score (the live-call
    // honestly reports the ambiguity in `v_call`) but for downstream
    // projection (germline_alignment, identity, CIGAR) we prefer the
    // truth allele when it sits inside the tie-set: it's the only
    // allele whose germline bytes really do match the pool bases. If
    // the truth isn't in the tie-set (mutations have shifted evidence
    // toward a different allele), fall back to the first tied id -
    // that's the genuine aligner-drift case the divergence narrative
    // captures.
    let truth_id = match segment {
        Segment::V => sim.assignments.get(Segment::V).copied().map(|i| i.allele_id),
        Segment::D => sim.assignments.get(Segment::D).copied().map(|i| i.allele_id),
        Segment::J => sim.assignments.get(Segment::J).copied().map(|i| i.allele_id),
        _ => None,
    };
    if let Some(call) = sim.segment_calls.get(segment) {
        if let Some(tid) = truth_id {
            if call.allele_call.contains(tid) {
                return Some(tid);
            }
        }
        if let Some(id) = call.allele_call.iter_ids().next() {
            return Some(id);
        }
    }
    truth_id
}

pub(super) fn projected_call_name(
    refdata: &RefDataConfig,
    sim: &Simulation,
    segment: Segment,
    origin_id: Option<crate::refdata::AlleleId>,
) -> String {
    // When the truth allele is inside the tie-set, list it FIRST in
    // the comma-separated `v_call` string. This keeps the first
    // allele in `v_call` aligned with `projected_allele_id`'s
    // truth-preference, so downstream consumers that read
    // `v_call.split(",")[0]` see the same allele that
    // `germline_alignment`, `v_identity`, and `v_cigar` are computed
    // against.
    live_call_name(
        refdata,
        sim.segment_calls.get(segment),
        segment,
        origin_id,
    )
    .unwrap_or_else(|| lookup_name(refdata, segment, origin_id))
}

/// Figure out which V/D/J segment (if any) has an extension
/// hypothesis that claims this NP pool position.
///
/// Returns `Some((segment, germline_byte))` when an adjacent V/D/J
/// live-call hypothesis extended past its structural region into
/// `pool_pos`. The germline byte comes from the segment's source
/// allele at the inferred reference position. Returns `None` when
/// no V/D/J claim exists - the position stays an NP column with `N`
/// in `germline_alignment`.
///
/// Tie-breaking when both V (right extension) and D (left extension)
/// claim the same NP1 position: V wins. Same for D vs J on NP2.
/// This matches the natural V->D->J order and is deterministic.
pub(super) fn np_claim_owner(
    sim: &Simulation,
    refdata: &RefDataConfig,
    pool_pos: usize,
) -> Option<(Segment, u8, u32)> {
    for &seg in Segment::assignable() {
        if let Some(claim) = check_segment_claim(sim, refdata, seg, pool_pos) {
            return Some(claim);
        }
    }
    None
}

/// build the AIRR `np1` / `np2` string from an NP
/// structural region by skipping any pool position that has been
/// claimed by a V/D/J live-call extension. The resulting string is
/// the "true" non-templated insertion as the evidence layer sees it.
/// pool-position bounds (`[start, end)`) of the contiguous unclaimed
/// interior of an NP region - the bytes that neither V's right
/// extension nor the adjacent segment's left extension claimed. Returns
/// `None` when every position in the region was claimed.
///
/// Used by `np1_aa` / `np2_aa` to slice `sequence_aa` at the codon
/// positions inside the unclaimed span, keeping the AA fields aligned
/// with `np_length` (which also counts only unclaimed bytes).
pub(super) fn unclaimed_np_bounds(
    sim: &Simulation,
    refdata: &RefDataConfig,
    region: &Region,
) -> Option<(usize, usize)> {
    let r_start = region.start.index() as usize;
    let r_end = region.end.index() as usize;
    let mut first: Option<usize> = None;
    let mut last: Option<usize> = None;
    for i in r_start..r_end {
        if np_claim_owner(sim, refdata, i).is_none() {
            if first.is_none() {
                first = Some(i);
            }
            last = Some(i);
        }
    }
    match (first, last) {
        (Some(s), Some(e)) => Some((s, e + 1)),
        _ => None,
    }
}

pub(super) fn unclaimed_np_string(
    sim: &Simulation,
    refdata: &RefDataConfig,
    sequence: &str,
    region: &Region,
) -> (String, i64) {
    let bytes = sequence.as_bytes();
    let r_start = region.start.index() as usize;
    let r_end = region.end.index() as usize;
    let mut out = String::with_capacity(r_end.saturating_sub(r_start));
    for i in r_start..r_end.min(bytes.len()) {
        if np_claim_owner(sim, refdata, i).is_none() {
            out.push(bytes[i] as char);
        }
    }
    let len = out.len() as i64;
    (out, len)
}

fn check_segment_claim(
    sim: &Simulation,
    refdata: &RefDataConfig,
    seg: Segment,
    pool_pos: usize,
) -> Option<(Segment, u8, u32)> {
    let call = sim.segment_calls.get(seg)?;
    let h = call.hypotheses.first()?;
    // NP-claim germline byte comes from the projected
    // allele (live-call's first allele, fallback to provenance) so
    // it matches `*_call`. See `projected_allele_id` for rationale.
    let allele_id = projected_allele_id(sim, seg)?;
    let allele = refdata.get(seg, allele_id)?;
    let region = sim.sequence.regions.iter().find(|r| r.segment == seg)?;
    let r_start = region.start.index() as usize;
    let r_end = region.end.index() as usize;
    let h_seq_start = h.seq_start as usize;
    let h_seq_end = h.seq_end as usize;
    // Position must be inside the live span...
    if pool_pos < h_seq_start || pool_pos >= h_seq_end {
        return None;
    }
    // ...and OUTSIDE the structural region (i.e. extension territory).
    if pool_pos >= r_start && pool_pos < r_end {
        return None;
    }
    // Linear interpolation: at `h.seq_start + k` the walker matched
    // ref pos `h.ref_start + k`. NP regions don't carry indels, so
    // the offset is monotonic.
    let offset = pool_pos - h_seq_start;
    let ref_pos = h.ref_start as usize + offset;
    if ref_pos >= allele.seq.len() {
        return None;
    }
    Some((seg, allele.seq[ref_pos], ref_pos as u32))
}

fn live_call_name(
    refdata: &RefDataConfig,
    call: Option<&SegmentLiveCall>,
    segment: Segment,
    truth_id: Option<crate::refdata::AlleleId>,
) -> Option<String> {
    let call = call?;
    // Order the tie-set so the truth allele appears first when it is
    // a member. The remaining alleles follow in iter_ids() (ascending
    // bitset) order. The first comma-separated entry of the resulting
    // string therefore matches `projected_allele_id`'s pick, which is
    // what `germline_alignment`, identity, and CIGAR are computed
    // against.
    let truth_in_tie = truth_id
        .map(|tid| call.allele_call.contains(tid))
        .unwrap_or(false);
    let mut ordered_ids: Vec<crate::refdata::AlleleId> = Vec::new();
    if let (Some(tid), true) = (truth_id, truth_in_tie) {
        ordered_ids.push(tid);
    }
    for id in call.allele_call.iter_ids() {
        if Some(id) == truth_id && truth_in_tie {
            continue;
        }
        ordered_ids.push(id);
    }
    let names: Vec<String> = ordered_ids
        .into_iter()
        .filter_map(|id| refdata.get(segment, id).map(|allele| allele.name.clone()))
        .collect();
    // Defensive fallback: under the score-and-tie caller, the bitset
    // should always have at least one allele (max=0 returns the full
    // pool). An empty bitset only arises today when an upstream
    // constructor builds an unsupported live-call directly. Treat
    // that as "no evidence-based call available" so projected_call_name's
    // unwrap_or_else fires and v_call falls back to the truth assignment,
    // staying consistent with projected_allele_id's existing fallback.
    if names.is_empty() {
        return None;
    }
    Some(names.join(","))
}

pub(super) fn lookup_allele(
    refdata: &RefDataConfig,
    segment: Segment,
    id: Option<crate::refdata::AlleleId>,
) -> Option<&Allele> {
    id.and_then(|aid| refdata.get(segment, aid))
}
