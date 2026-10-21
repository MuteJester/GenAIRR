use crate::ir::{complement_base, Segment, Simulation};
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

/// Standard AIRR alignment fields describe the generating assignment.
/// Mutational evidence can change a re-call, but cannot change the origin
/// reference frame carried by the nucleotide pool's germline coordinates.
pub(super) fn projected_allele_id(
    sim: &Simulation,
    segment: Segment,
) -> Option<crate::refdata::AlleleId> {
    sim.assignments
        .get(segment)
        .map(|instance| instance.allele_id)
}

/// The generating reference in the segment's assembled orientation.
///
/// Germline coordinates always refer to the forward reference; only traversal
/// order and rendered bases change for inverted D segments.
pub(super) struct ProjectedReference<'a> {
    pub seq: &'a [u8],
    pub reverse: bool,
    observed_start: Option<usize>,
    observed_end: Option<usize>,
}

impl ProjectedReference<'_> {
    pub fn oriented_position(&self, reference_pos: usize) -> Option<usize> {
        if reference_pos >= self.seq.len() {
            return None;
        }
        Some(if self.reverse {
            self.seq.len() - 1 - reference_pos
        } else {
            reference_pos
        })
    }

    pub fn reference_position(&self, oriented_pos: usize) -> usize {
        if self.reverse {
            self.seq.len() - 1 - oriented_pos
        } else {
            oriented_pos
        }
    }

    pub fn base(&self, reference_pos: usize) -> Option<u8> {
        self.seq.get(reference_pos).map(|&base| {
            if self.reverse {
                complement_base(base)
            } else {
                base
            }
        })
    }

    /// Reference interval retained by the generating recombination.
    pub fn structural_bounds(&self, trim_5: i64, trim_3: i64) -> (usize, usize) {
        let end = self
            .seq
            .len()
            .saturating_sub(trim_3.max(0) as usize)
            .min(self.seq.len());
        ((trim_5.max(0) as usize).min(end), end)
    }

    /// Reference interval covered by this observation. Explicit clipping
    /// removes reference coverage; structural deletions within the remaining
    /// interval still contribute deletion columns.
    pub fn observed_bounds(&self, trim_5: i64, trim_3: i64) -> (usize, usize) {
        let (start, end) = self.structural_bounds(trim_5, trim_3);
        let end = end.min(self.observed_end.unwrap_or(end));
        (
            start.max(self.observed_start.unwrap_or(start)).min(end),
            end,
        )
    }
}

pub(super) fn projected_reference<'a>(
    sim: &Simulation,
    refdata: &'a RefDataConfig,
    segment: Segment,
) -> Option<ProjectedReference<'a>> {
    let allele = refdata.get(segment, projected_allele_id(sim, segment)?)?;
    let reverse = segment == Segment::D
        && sim
            .assignments
            .get(segment)
            .map(|instance| instance.orientation.is_reverse())
            .unwrap_or(false);
    Some(ProjectedReference {
        seq: &allele.seq,
        reverse,
        observed_start: sim
            .assignments
            .get(segment)
            .and_then(|instance| instance.observed_ref_start)
            .map(|n| n as usize),
        observed_end: sim
            .assignments
            .get(segment)
            .and_then(|instance| instance.observed_ref_end)
            .map(|n| n as usize),
    })
}

/// Emit a structural segment's alignment columns in query order.
///
/// The shared traversal keeps CIGAR/identity/range calculation and string
/// rendering on the same reference. Coordinates in each emitted column are
/// `(pool_position, forward_reference_position)`; a missing position denotes an
/// insertion/deletion gap. Reference bounds and bases come from the same
/// generating allele as the pool's germline coordinates.
#[allow(clippy::too_many_arguments)]
pub(super) fn visit_structural_alignment(
    sim: &Simulation,
    reference: &ProjectedReference<'_>,
    pool_start: usize,
    pool_end: usize,
    trim_5: i64,
    trim_3: i64,
    mut emit: impl FnMut(Option<usize>, Option<usize>, u8, u8),
) {
    let len = reference.seq.len();
    let (start, end) = reference.observed_bounds(trim_5, trim_3);
    let (start, end) = if reference.reverse {
        (len - end, len - start)
    } else {
        (start, end)
    };
    let mut expected = start.min(len);

    for pool_pos in pool_start..pool_end {
        let nucleotide = &sim.pool.as_slice()[pool_pos];
        let position = nucleotide
            .germline_pos
            .get()
            .and_then(|position| reference.oriented_position(position as usize));
        let Some(position) = position else {
            emit(Some(pool_pos), None, nucleotide.base, b'-');
            continue;
        };
        while expected < position {
            let ref_pos = reference.reference_position(expected);
            emit(None, Some(ref_pos), b'-', reference.base(ref_pos).unwrap());
            expected += 1;
        }
        let ref_pos = reference.reference_position(position);
        emit(
            Some(pool_pos),
            Some(ref_pos),
            nucleotide.base,
            reference.base(ref_pos).unwrap(),
        );
        expected = position + 1;
    }
    while expected < end {
        let ref_pos = reference.reference_position(expected);
        emit(None, Some(ref_pos), b'-', reference.base(ref_pos).unwrap());
        expected += 1;
    }
}

pub(super) fn projected_call_name(
    refdata: &RefDataConfig,
    _sim: &Simulation,
    segment: Segment,
    origin_id: Option<crate::refdata::AlleleId>,
) -> String {
    lookup_name(refdata, segment, origin_id)
}

/// Evidence-based re-calls remain available separately from the standard
/// generating calls and alignment reference.
pub(crate) fn recalled_call_name(
    refdata: &RefDataConfig,
    sim: &Simulation,
    segment: Segment,
    origin_id: Option<crate::refdata::AlleleId>,
) -> String {
    live_call_name(refdata, sim.segment_calls.get(segment), segment, origin_id).unwrap_or_default()
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
    // bitset) order. These are evidence-based calls, not alignment
    // reference choices.
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
    // An observed all-zero-score call still contains the supported-pool tie.
    // An absent/unsupported call has no inference to export; generating
    // identity is already available separately in the standard call.
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
