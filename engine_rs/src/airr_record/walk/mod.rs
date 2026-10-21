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
//! identity) lives in
//! [`super::segment_projection::compute_segment_projections`]. The
//! renderer **consumes** that truth and uses it to decide which
//! germline byte (or `-` gap, or `N`) each column emits. It owns no
//! per-segment bookkeeping of its own — it tracks the three output strings
//! and a pool cursor so junction P bases and unowned insertions are emitted
//! at their actual positions.
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

use super::projection::{projected_reference, visit_structural_alignment};
use super::segment_projection::compute_segment_projections;
use super::sequence::{bytes_to_string, pool_bases};
use super::AirrRecord;
use crate::ir::{flag, PoolRange, RefRange, Segment, Simulation};
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

    let mut sa: Vec<u8> = Vec::with_capacity(pool_len + 8);
    let mut galn: Vec<u8> = Vec::with_capacity(pool_len + 8);
    let mut dmask: Vec<u8> = Vec::with_capacity(pool_len + 8);
    let mut covered_until = 0usize;

    for region in &sim.sequence.regions {
        let seg = region.segment;
        let r_start = region.start.index() as usize;
        let r_end = region.end.index() as usize;
        // P additions are intentionally outside structural regions. Render
        // these intervening pool positions here, not after the final J base.
        render_uncovered(
            sim,
            &bases,
            covered_until..r_start,
            &mut sa,
            &mut galn,
            &mut dmask,
        );
        covered_until = r_end;

        match seg {
            Segment::V | Segment::D | Segment::J => {
                render_structural_region(
                    sim, refdata, rec, &bases, seg, r_start, r_end, &mut sa, &mut galn, &mut dmask,
                );
            }
            Segment::Np1 | Segment::Np2 => {
                for p in r_start..r_end {
                    sa.push(bases[p]);
                    galn.push(b'N');
                    dmask.push(b'N');
                }
            }
        }
    }

    render_uncovered(
        sim,
        &bases,
        covered_until..pool_len,
        &mut sa,
        &mut galn,
        &mut dmask,
    );

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

fn render_uncovered(
    sim: &Simulation,
    bases: &[u8],
    range: std::ops::Range<usize>,
    sa: &mut Vec<u8>,
    galn: &mut Vec<u8>,
    dmask: &mut Vec<u8>,
) {
    for i in range {
        sa.push(bases[i]);
        // P bases belong to the combined N/P junction, not to a germline
        // segment or an indel. Other unowned bases are tail indel insertions.
        let germline = if sim.pool.as_slice()[i].flags.contains(flag::P_NUC) {
            b'N'
        } else {
            b'-'
        };
        galn.push(germline);
        dmask.push(germline);
    }
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
    sa: &mut Vec<u8>,
    galn: &mut Vec<u8>,
    dmask: &mut Vec<u8>,
) {
    let reference = projected_reference(sim, refdata, seg);
    let (trim_5, trim_3) = match seg {
        Segment::V => (rec.v_trim_5, rec.v_trim_3),
        Segment::D => (rec.d_trim_5, rec.d_trim_3),
        Segment::J => (rec.j_trim_5, rec.j_trim_3),
        _ => unreachable!(),
    };

    let Some(reference) = reference else {
        // No source allele assigned; emit raw bases with 'N'
        // germline so column lengths stay consistent.
        for i in r_start..r_end {
            sa.push(bases[i]);
            galn.push(b'N');
            push_dmask_for_seg(dmask, seg, b'N');
        }
        return;
    };

    visit_structural_alignment(
        sim,
        &reference,
        r_start,
        r_end,
        trim_5,
        trim_3,
        |_, _, sequence_base, germline_base| {
            sa.push(sequence_base);
            galn.push(germline_base);
            push_dmask_for_seg(dmask, seg, germline_base);
        },
    );
}
