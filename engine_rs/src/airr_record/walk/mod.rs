//! Alignment-column renderer.
//!
//! Builds the assembled `sequence_alignment`, `germline_alignment`,
//! and `germline_alignment_d_mask` strings for an AIRR record by
//! walking every region in the simulation once and emitting one
//! column per assembled position.
//!
//! ## Ownership
//!
//! Segment-level truth (allele choice, sequence/ref ranges, CIGAR,
//! identity, accepted NP claims) lives in
//! [`super::segment_projection::compute_segment_projections`]. The
//! renderer **consumes** that truth and uses it to decide which
//! germline byte (or `-` gap, or `N`) each column emits. It owns no
//! per-segment bookkeeping of its own — the only state it tracks is
//! the three output strings and a pool-coverage bitmap used to
//! locate tail-end indel insertions.
//!
//! This split was driven by historical drift: the previous walker
//! interleaved segment truth and rendering in one 400-line function,
//! so NP-claim bounds, segment range arithmetic, and CIGAR rules
//! lived in two places (here and the validator) and could disagree.
//! The kernel collapses them.

pub(in crate::airr_record) mod helpers;
mod types;

pub(in crate::airr_record) use types::AlignmentWalk;

use helpers::{bytes_uppercase_in_place, push_dmask_for_seg};

use super::projection::projected_allele_id;
use super::segment_projection::{compute_segment_projections, NpClaim};
use super::sequence::{bytes_to_string, pool_bases};
use super::AirrRecord;
use crate::ir::{Nucleotide, PoolRange, RefRange, Segment, Simulation};
use crate::refdata::RefDataConfig;

pub(super) fn walk_alignment_columns(
    sim: &Simulation,
    refdata: &RefDataConfig,
    rec: &AirrRecord,
) -> AlignmentWalk {
    // Phase 3: pre-compute segment truth. The renderer reads the
    // kernel's decisions and never recomputes them.
    let projections = compute_segment_projections(sim, refdata, rec);

    let bases = pool_bases(sim);
    let pool_len = bases.len();

    // NP-claim lookup keyed by pool position. Built from the
    // kernel's accepted claims (positions blocked by the kernel's
    // structural / past-allele / already-covered guards never make
    // it into this list — those columns render as plain NP).
    let np_claim_by_pool = build_np_claim_lookup(&projections.np_claims, pool_len);

    let mut sa: Vec<u8> = Vec::with_capacity(pool_len + 8);
    let mut galn: Vec<u8> = Vec::with_capacity(pool_len + 8);
    let mut dmask: Vec<u8> = Vec::with_capacity(pool_len + 8);
    let mut covered: Vec<bool> = vec![false; pool_len];

    for region in &sim.sequence.regions {
        let seg = region.segment;
        let r_start = region.start.index() as usize;
        let r_end = region.end.index() as usize;
        for idx in r_start..r_end {
            covered[idx] = true;
        }

        match seg {
            Segment::V | Segment::D | Segment::J => {
                // `np_covered_end` for this segment must reflect ONLY
                // the NP claims processed BEFORE this structural region
                // — i.e. claims whose pool position is strictly less
                // than `r_start`. The kernel's running `ref_ranges[seg]`
                // at the moment this region begins iteration carries
                // exactly that value; here we recover it by filtering
                // the flat claim list on pool position.
                //
                // Without this filter, V-right-extension claims (which
                // sit in NP1, *after* V's structural region) would
                // bleed back into V's `expected_pos`, shifting V's
                // D-fill column count by one and corrupting the
                // assembled `germline_alignment`.
                let np_covered_end = projections
                    .np_claims
                    .iter()
                    .filter(|c| c.owner == seg && (c.pool_pos as usize) < r_start)
                    .map(|c| c.ref_pos as usize + 1)
                    .max()
                    .unwrap_or(0);
                render_structural_region(
                    sim, refdata, rec, &bases, seg, r_start, r_end,
                    np_covered_end, &mut sa, &mut galn, &mut dmask,
                );
            }
            Segment::Np1 | Segment::Np2 => {
                for p in r_start..r_end {
                    sa.push(bases[p]);
                    match np_claim_by_pool[p] {
                        Some((claim_seg, allele_byte)) => {
                            galn.push(allele_byte);
                            push_dmask_for_seg(&mut dmask, claim_seg, allele_byte);
                        }
                        None => {
                            galn.push(b'N');
                            dmask.push(b'N');
                        }
                    }
                }
            }
        }
    }

    // Tail-end indel insertions outside every region.
    for i in 0..pool_len {
        if !covered[i] {
            sa.push(bases[i]);
            galn.push(b'-');
            dmask.push(b'-');
        }
    }

    bytes_uppercase_in_place(&mut sa);
    bytes_uppercase_in_place(&mut galn);
    bytes_uppercase_in_place(&mut dmask);

    // Per-segment outputs all come from the kernel.
    let seq_range_tuple = |r: Option<PoolRange>| r.map(|p| (p.start as i64, p.end as i64));
    let ref_range_tuple = |r: Option<RefRange>| r.map(|p| (p.start, p.end));

    AlignmentWalk {
        sa: bytes_to_string(&sa),
        galn: bytes_to_string(&galn),
        dmask: bytes_to_string(&dmask),
        cigars: [
            projections.v.cigar.clone(),
            projections.d.cigar.clone(),
            projections.j.cigar.clone(),
        ],
        align_ranges: [
            projections.v.alignment_range,
            projections.d.alignment_range,
            projections.j.alignment_range,
        ],
        seq_ranges: [
            seq_range_tuple(projections.v.sequence_range),
            seq_range_tuple(projections.d.sequence_range),
            seq_range_tuple(projections.j.sequence_range),
        ],
        ref_ranges: [
            ref_range_tuple(projections.v.ref_range),
            ref_range_tuple(projections.d.ref_range),
            ref_range_tuple(projections.j.ref_range),
        ],
        identities: [
            projections.v.identity,
            projections.d.identity,
            projections.j.identity,
        ],
    }
}

fn build_np_claim_lookup(claims: &[NpClaim], pool_len: usize) -> Vec<Option<(Segment, u8)>> {
    let mut lookup: Vec<Option<(Segment, u8)>> = vec![None; pool_len];
    for c in claims {
        let p = c.pool_pos as usize;
        if p < pool_len {
            lookup[p] = Some((c.owner, c.allele_byte));
        }
    }
    lookup
}

#[allow(clippy::too_many_arguments)]
fn render_structural_region(
    sim: &Simulation,
    refdata: &RefDataConfig,
    rec: &AirrRecord,
    bases: &[u8],
    seg: Segment,
    r_start: usize,
    r_end: usize,
    np_covered_end: usize,
    sa: &mut Vec<u8>,
    galn: &mut Vec<u8>,
    dmask: &mut Vec<u8>,
) {
    let allele_id = projected_allele_id(sim, seg);
    let allele_seq: Option<&[u8]> = allele_id
        .and_then(|aid| refdata.get(seg, aid))
        .map(|a| a.seq.as_slice());
    let (trim_5, trim_3) = match seg {
        Segment::V => (rec.v_trim_5, rec.v_trim_3),
        Segment::D => (rec.d_trim_5, rec.d_trim_3),
        Segment::J => (rec.j_trim_5, rec.j_trim_3),
        _ => unreachable!(),
    };

    let Some(allele_seq) = allele_seq else {
        // No source allele assigned; emit raw bases with 'N'
        // germline so column lengths stay consistent.
        for i in r_start..r_end {
            sa.push(bases[i]);
            galn.push(b'N');
            push_dmask_for_seg(dmask, seg, b'N');
        }
        return;
    };

    // NP-region left-extension claims may have already covered some
    // ref positions for this segment. Start `expected_pos` past
    // them so a structural-indel deletion at the leftmost germline_pos
    // doesn't trigger a D-fill at a ref position the NP claim
    // already counted. `np_covered_end` is computed by the caller
    // from the subset of claims processed before this region begins.
    let mut expected_pos: usize = (trim_5.max(0) as usize).max(np_covered_end);
    let end_germ: usize = (allele_seq.len() as i64 - trim_3).max(0) as usize;

    for i in r_start..r_end {
        let nuc: &Nucleotide = &sim.pool.as_slice()[i];
        let base_char = bases[i];
        let Some(germ_pos) = nuc.germline_pos.get() else {
            // Indel insertion: gap in germline. The walker emits one
            // I column. The kernel already counted the I op in the
            // segment's CIGAR.
            sa.push(base_char);
            galn.push(b'-');
            push_dmask_for_seg(dmask, seg, b'-');
            continue;
        };
        let germ_pos = germ_pos as usize;
        // Fill any preceding deletion gap.
        while expected_pos < germ_pos && expected_pos < allele_seq.len() {
            sa.push(b'-');
            let g = allele_seq[expected_pos];
            galn.push(g);
            push_dmask_for_seg(dmask, seg, g);
            expected_pos += 1;
        }
        // Match/mismatch column. The germline byte comes from the
        // projected allele's sequence so it matches `*_call` (live-
        // call narrowed to this allele) rather than the originally-
        // sampled `nuc.germline`.
        sa.push(base_char);
        let g = if germ_pos < allele_seq.len() {
            allele_seq[germ_pos]
        } else {
            nuc.germline
        };
        galn.push(g);
        push_dmask_for_seg(dmask, seg, g);
        expected_pos = germ_pos + 1;
    }
    // Trailing deletion gaps: positions in `[expected_pos, end_germ)`
    // were trimmed off by an indel pass, not by recombination.
    while expected_pos < end_germ && expected_pos < allele_seq.len() {
        sa.push(b'-');
        let g = allele_seq[expected_pos];
        galn.push(g);
        push_dmask_for_seg(dmask, seg, g);
        expected_pos += 1;
    }
}
