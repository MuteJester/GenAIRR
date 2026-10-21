//! Per-segment generating alignment coordinates, CIGARs and identities.
//!
//! The renderer consumes these projections. Generating assignments and
//! nucleotide origin determine segment membership; live re-call hypotheses
//! never reassign junction bases or choose the alignment reference.

use crate::ir::{PoolRange, RefRange, Segment, Simulation};
use crate::refdata::{AlleleId, RefDataConfig};

use super::projection::{projected_allele_id, projected_reference, visit_structural_alignment};
use super::walk::helpers::{
    eq_ascii_case_insensitive, extend_ref_range, push_cigar_op, runlength_to_string,
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
    /// Pool-position span of the generating structural region. Half-open.
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

/// All three generating segment projections.
#[derive(Clone, Debug)]
pub(in crate::airr_record) struct SegmentProjectionSet {
    pub v: SegmentProjection,
    pub d: SegmentProjection,
    pub j: SegmentProjection,
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

/// Calculate alignment-column spans and per-segment summaries using the
/// same structural traversal that renders the actual alignment strings.
pub(in crate::airr_record) fn compute_segment_projections(
    sim: &Simulation,
    refdata: &RefDataConfig,
    rec: &AirrRecord,
) -> SegmentProjectionSet {
    let mut cigar_runs: [Vec<(u32, u8)>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    let mut align_ranges: [Option<(i64, i64)>; 3] = [None, None, None];
    let mut seq_ranges: [Option<(i64, i64)>; 3] = [None, None, None];
    let mut ref_ranges: [Option<(i64, i64)>; 3] = [None, None, None];
    let mut id_counts: [[u64; 2]; 3] = [[0, 0], [0, 0], [0, 0]];

    // Column counter — mirrors `sa.len()` in walk_alignment_columns.
    // Every position the assembled-alignment string would push, we
    // bump this by 1 instead.
    let mut col: i64 = 0;
    let mut covered_until = 0usize;

    for region in &sim.sequence.regions {
        let seg = region.segment;
        let r_start = region.start.index() as usize;
        let r_end = region.end.index() as usize;
        // P additions occupy pool positions between structural regions.
        // They contribute alignment columns even though no V/D/J owns them.
        col += r_start.saturating_sub(covered_until) as i64;
        covered_until = r_end;

        match seg {
            Segment::V | Segment::D | Segment::J => {
                let reference = projected_reference(sim, refdata, seg);
                let (trim_5, trim_3) = trims_for(rec, seg);
                let i = seg_idx(seg);

                let span_start = col;

                if let Some(reference) = reference {
                    visit_structural_alignment(
                        sim,
                        &reference,
                        r_start,
                        r_end,
                        trim_5,
                        trim_3,
                        |pool_pos, ref_pos, base, germline| {
                            col += 1;
                            let op = match (pool_pos, ref_pos) {
                                (None, _) => b'D',
                                (_, None) => b'I',
                                _ => b'M',
                            };
                            push_cigar_op(&mut cigar_runs[i], op);
                            id_counts[i][1] += 1;
                            if op == b'M' && eq_ascii_case_insensitive(base, germline) {
                                id_counts[i][0] += 1;
                            }
                            if let Some(ref_pos) = ref_pos {
                                extend_ref_range(&mut ref_ranges, i, ref_pos as i64);
                            }
                        },
                    );
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
                // Sequence similarity does not change a generated N/P base's
                // origin. Re-call extension hypotheses remain in the live IR.
                col += (r_end - r_start) as i64;
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
            sequence_range: seq_ranges[i].map(|(s, e)| PoolRange::new(s as u32, e as u32)),
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
    }
}
