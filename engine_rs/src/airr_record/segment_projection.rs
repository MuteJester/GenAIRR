//! Segment projection kernel — "what does each V/D/J segment
//! project to?"
//!
//! Renders, per segment, the truth decisions that drive every
//! AIRR-record alignment field: which allele is projected, what
//! pool slice and ref slice the segment owns (structural region
//! plus any NP claims it won), the CIGAR run, identity ratio, and
//! alignment-column span. Separately surfaces the accepted NP
//! claims as a flat [`Vec<NpClaim>`] so callers don't have to
//! re-derive ownership later.
//!
//! ## Scope
//!
//! This module is the "segment truth" layer. It deliberately does
//! NOT render the assembled alignment strings (`sequence_alignment`,
//! `germline_alignment`, `dmask`): those remain in
//! [`super::walk`]. The two are intended to be consumed together —
//! the walker will eventually take per-segment decisions from
//! [`compute_segment_projections`] and render columns on top.
//!
//! Phase 1: kernel lives beside the existing walker, computes the
//! same numbers, and is verified by parity tests against
//! [`super::walk::walk_alignment_columns`] (see `tests::projection_parity`).
//! Phases 2+ swap the walker to consume this kernel; phase 4 deletes
//! the duplicated bookkeeping. No AIRR output change until then.
//!
//! ## Why this split exists
//!
//! `walk_alignment_columns` previously owned segment allele choice,
//! sequence/germline/ref ranges, CIGAR generation, identity counting,
//! deletion filling, insertion handling, NP claim ownership, NP claim
//! bounds, AND full alignment string rendering — all interleaved in
//! one ~400-line function. Splitting "segment truth" out lets the
//! column walker render against an immutable decision set rather than
//! deriving truth and rendering it in lockstep.

use crate::ir::{Nucleotide, PoolRange, RefRange, Segment, Simulation};
use crate::refdata::{AlleleId, RefDataConfig};

use super::projection::{np_claim_owner, projected_allele_id};
use super::sequence::pool_bases;
use super::walk::helpers::{
    eq_ascii_case_insensitive, extend_ref_range, push_cigar_op, ref_pos_already_covered,
    runlength_to_string,
};
use super::AirrRecord;

/// What one V/D/J segment projects to under the current simulation +
/// refdata. Fields are `None` when the segment didn't contribute to
/// the projection (no region present, or no allele assigned).
///
/// `segment` and `projected_allele_id` are self-describing fields
/// surfaced for downstream consumers (e.g. validators, future
/// projection rewrites); the alignment renderer reads only the
/// range / cigar / identity fields.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub(in crate::airr_record) struct SegmentProjection {
    pub segment: Segment,
    pub projected_allele_id: Option<AlleleId>,
    /// Pool-position span the segment owns (structural region plus
    /// any NP positions it won via extension claims). Half-open.
    pub sequence_range: Option<PoolRange>,
    /// Ref-position span the segment consumed (union of `M` and `D`
    /// ops). Half-open. `germline_span == M + D` by construction.
    pub ref_range: Option<RefRange>,
    /// CIGAR string in standard `<count><op>` form.
    pub cigar: String,
    /// Identity ratio (matches / total aligned positions), or
    /// `None` when the segment contributed no aligned columns.
    pub identity: Option<f64>,
    /// Alignment-column span (offsets into the assembled
    /// `sequence_alignment` string). Half-open.
    pub alignment_range: Option<(i64, i64)>,
}

/// One accepted NP-extension claim. Reflects the live-call
/// hypothesis on segment `owner` extending into an NP region and
/// winning the position (NP claims that were blocked by the
/// structural overlap / past-allele / already-covered guards do NOT
/// appear here — only emitted claims do).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(in crate::airr_record) struct NpClaim {
    pub pool_pos: u32,
    pub owner: Segment,
    pub ref_pos: u32,
    pub allele_byte: u8,
}

/// All three segment projections plus the flat NP-claim list.
#[derive(Clone, Debug)]
pub(in crate::airr_record) struct SegmentProjectionSet {
    pub v: SegmentProjection,
    pub d: SegmentProjection,
    pub j: SegmentProjection,
    pub np_claims: Vec<NpClaim>,
}

impl SegmentProjectionSet {
    #[allow(dead_code)]
    pub(in crate::airr_record) fn get(&self, seg: Segment) -> Option<&SegmentProjection> {
        match seg {
            Segment::V => Some(&self.v),
            Segment::D => Some(&self.d),
            Segment::J => Some(&self.j),
            _ => None,
        }
    }
}

fn seg_idx(seg: Segment) -> usize {
    match seg {
        Segment::V => 0,
        Segment::D => 1,
        Segment::J => 2,
        _ => unreachable!("seg_idx called on non-VDJ segment"),
    }
}

fn trims_for(rec: &AirrRecord, seg: Segment) -> (i64, i64) {
    match seg {
        Segment::V => (rec.v_trim_5, rec.v_trim_3),
        Segment::D => (rec.d_trim_5, rec.d_trim_3),
        Segment::J => (rec.j_trim_5, rec.j_trim_3),
        _ => (0, 0),
    }
}

/// Build the per-segment projection set in one pass over the
/// simulation regions.
///
/// Mirrors the bookkeeping in [`super::walk::walk_alignment_columns`]
/// exactly — same CIGAR rules, same NP claim guards, same
/// trim-derived structural bounds, same span-merge semantics across
/// NP and structural contributions to the segment. Differs only in
/// that it tracks a `col` counter in place of the walker's `sa`
/// vector so the assembled-alignment string is not constructed.
pub(in crate::airr_record) fn compute_segment_projections(
    sim: &Simulation,
    refdata: &RefDataConfig,
    rec: &AirrRecord,
) -> SegmentProjectionSet {
    let bases = pool_bases(sim);

    let mut cigar_runs: [Vec<(u32, u8)>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    let mut align_ranges: [Option<(i64, i64)>; 3] = [None, None, None];
    let mut seq_ranges: [Option<(i64, i64)>; 3] = [None, None, None];
    let mut ref_ranges: [Option<(i64, i64)>; 3] = [None, None, None];
    let mut id_counts: [[u64; 2]; 3] = [[0, 0], [0, 0], [0, 0]];
    let mut np_claims: Vec<NpClaim> = Vec::new();

    // Column counter — mirrors `sa.len()` in walk_alignment_columns.
    // Every position the assembled-alignment string would push, we
    // bump this by 1 instead.
    let mut col: i64 = 0;

    for region in &sim.sequence.regions {
        let seg = region.segment;
        let r_start = region.start.index() as usize;
        let r_end = region.end.index() as usize;

        match seg {
            Segment::V | Segment::D | Segment::J => {
                let allele_id = projected_allele_id(sim, seg);
                let allele_seq: Option<&[u8]> = allele_id
                    .and_then(|aid| refdata.get(seg, aid))
                    .map(|a| a.seq.as_slice());
                let (trim_5, trim_3) = trims_for(rec, seg);
                let i = seg_idx(seg);

                let span_start = col;

                if let Some(allele_seq) = allele_seq {
                    let np_covered_end =
                        ref_ranges[i].map(|(_, e)| e as usize).unwrap_or(0);
                    let mut expected_pos: usize =
                        (trim_5.max(0) as usize).max(np_covered_end);
                    let end_germ: usize = (allele_seq.len() as i64 - trim_3).max(0) as usize;
                    for p in r_start..r_end {
                        let nuc: &Nucleotide = &sim.pool.as_slice()[p];
                        let base_char = bases[p];
                        let Some(germ_pos) = nuc.germline_pos.get() else {
                            // Indel insertion: gap in germline. One column.
                            col += 1;
                            push_cigar_op(&mut cigar_runs[i], b'I');
                            id_counts[i][1] += 1;
                            continue;
                        };
                        let germ_pos = germ_pos as usize;
                        // Fill preceding deletion gaps.
                        while expected_pos < germ_pos && expected_pos < allele_seq.len() {
                            col += 1;
                            push_cigar_op(&mut cigar_runs[i], b'D');
                            id_counts[i][1] += 1;
                            extend_ref_range(&mut ref_ranges, i, expected_pos as i64);
                            expected_pos += 1;
                        }
                        // Match/mismatch column.
                        col += 1;
                        let g = if germ_pos < allele_seq.len() {
                            allele_seq[germ_pos]
                        } else {
                            nuc.germline
                        };
                        push_cigar_op(&mut cigar_runs[i], b'M');
                        id_counts[i][1] += 1;
                        if eq_ascii_case_insensitive(base_char, g) {
                            id_counts[i][0] += 1;
                        }
                        extend_ref_range(&mut ref_ranges, i, germ_pos as i64);
                        expected_pos = germ_pos + 1;
                    }
                    // Trailing deletion gaps.
                    while expected_pos < end_germ && expected_pos < allele_seq.len() {
                        col += 1;
                        push_cigar_op(&mut cigar_runs[i], b'D');
                        id_counts[i][1] += 1;
                        extend_ref_range(&mut ref_ranges, i, expected_pos as i64);
                        expected_pos += 1;
                    }
                } else {
                    // No source allele: bases pass through unaligned.
                    col += (r_end - r_start) as i64;
                }

                let span_end = col;
                if span_end > span_start {
                    let new_start = match align_ranges[i] {
                        Some((s, _)) => s.min(span_start),
                        None => span_start,
                    };
                    align_ranges[i] = Some((new_start, span_end));
                    let pool_start = r_start as i64;
                    let pool_end = r_end as i64;
                    let new_pool_start = match seq_ranges[i] {
                        Some((s, _)) => s.min(pool_start),
                        None => pool_start,
                    };
                    seq_ranges[i] = Some((new_pool_start, pool_end));
                }
            }
            Segment::Np1 | Segment::Np2 => {
                let mut np_blocked: [bool; 3] = [false; 3];
                for p in r_start..r_end {
                    let base_char = bases[p];
                    col += 1;
                    let Some((claim_seg, allele_byte_proj, ref_pos_proj)) =
                        np_claim_owner(sim, refdata, p)
                    else {
                        continue;
                    };
                    let claim_idx = seg_idx(claim_seg);
                    let (ref_pos, allele_byte) = match ref_ranges[claim_idx] {
                        Some((_, e)) => {
                            let alid = projected_allele_id(sim, claim_seg);
                            let allele = alid.and_then(|aid| refdata.get(claim_seg, aid));
                            let byte = allele
                                .and_then(|a| a.seq.get(e as usize).copied())
                                .unwrap_or(b'N');
                            (e as u32, byte)
                        }
                        None => (ref_pos_proj, allele_byte_proj),
                    };
                    // Extension-territory bounds — see walk/mod.rs
                    // for the long-form story. Boundary semantics
                    // live in the [`RefRange::contains`] check below.
                    let (claim_in_structural, claim_past_allele) = {
                        let claim_alid = projected_allele_id(sim, claim_seg);
                        let claim_allele =
                            claim_alid.and_then(|aid| refdata.get(claim_seg, aid));
                        let claim_allele_len = claim_allele
                            .map(|a| a.seq.len() as i64)
                            .unwrap_or(0);
                        let (claim_t5, claim_t3) = trims_for(rec, claim_seg);
                        let structural = RefRange::new(
                            claim_t5,
                            (claim_allele_len - claim_t3).max(claim_t5),
                        );
                        let rp = ref_pos as i64;
                        (structural.contains(rp), rp >= claim_allele_len)
                    };
                    if np_blocked[claim_idx]
                        || ref_pos_already_covered(&ref_ranges, claim_idx, ref_pos as i64)
                        || claim_in_structural
                        || claim_past_allele
                    {
                        np_blocked[claim_idx] = true;
                        continue;
                    }
                    push_cigar_op(&mut cigar_runs[claim_idx], b'M');
                    id_counts[claim_idx][1] += 1;
                    if eq_ascii_case_insensitive(base_char, allele_byte) {
                        id_counts[claim_idx][0] += 1;
                    }
                    extend_ref_range(&mut ref_ranges, claim_idx, ref_pos as i64);
                    let col_pos = col - 1;
                    align_ranges[claim_idx] = Some(match align_ranges[claim_idx] {
                        Some((s, _)) => (s, col_pos + 1),
                        None => (col_pos, col_pos + 1),
                    });
                    let pool_pos = p as i64;
                    seq_ranges[claim_idx] = Some(match seq_ranges[claim_idx] {
                        Some((s, e)) => (s.min(pool_pos), e.max(pool_pos + 1)),
                        None => (pool_pos, pool_pos + 1),
                    });
                    np_claims.push(NpClaim {
                        pool_pos: p as u32,
                        owner: claim_seg,
                        ref_pos,
                        allele_byte,
                    });
                }
            }
        }
    }

    // Tail-end indel insertions (positions outside every region)
    // don't affect any per-segment range — they only bump `col` in
    // the walker, which we ignore here since the kernel doesn't
    // surface a global column count.

    let identity = |idx: usize| -> Option<f64> {
        if id_counts[idx][1] > 0 {
            Some(id_counts[idx][0] as f64 / id_counts[idx][1] as f64)
        } else {
            None
        }
    };

    let build = |seg: Segment| -> SegmentProjection {
        let i = seg_idx(seg);
        SegmentProjection {
            segment: seg,
            projected_allele_id: projected_allele_id(sim, seg),
            sequence_range: seq_ranges[i]
                .map(|(s, e)| PoolRange::new(s as u32, e as u32)),
            ref_range: ref_ranges[i].map(|(s, e)| RefRange::new(s, e)),
            cigar: runlength_to_string(&cigar_runs[i]),
            identity: identity(i),
            alignment_range: align_ranges[i],
        }
    };

    SegmentProjectionSet {
        v: build(Segment::V),
        d: build(Segment::D),
        j: build(Segment::J),
        np_claims,
    }
}

