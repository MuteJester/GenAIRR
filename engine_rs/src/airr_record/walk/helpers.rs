use crate::ir::Segment;

pub(in crate::airr_record) fn bytes_uppercase_in_place(bytes: &mut [u8]) {
    for b in bytes.iter_mut() {
        if (*b).is_ascii_lowercase() {
            *b = (*b).to_ascii_uppercase();
        }
    }
}

pub(in crate::airr_record) fn eq_ascii_case_insensitive(a: u8, b: u8) -> bool {
    a.to_ascii_uppercase() == b.to_ascii_uppercase()
}

pub(in crate::airr_record) fn push_cigar_op(runs: &mut Vec<(u32, u8)>, op: u8) {
    if let Some(last) = runs.last_mut() {
        if last.1 == op {
            last.0 += 1;
            return;
        }
    }
    runs.push((1, op));
}

pub(in crate::airr_record) fn runlength_to_string(runs: &[(u32, u8)]) -> String {
    let mut s = String::with_capacity(runs.len() * 4);
    for (count, op) in runs {
        s.push_str(&count.to_string());
        s.push(*op as char);
    }
    s
}

pub(in crate::airr_record) fn push_dmask_for_seg(dmask: &mut Vec<u8>, seg: Segment, ga_char: u8) {
    if seg == Segment::D && ga_char != b'-' {
        dmask.push(b'N');
    } else {
        dmask.push(ga_char);
    }
}

/// Extend `ranges[idx]` so it covers the (single) ref position
/// `ref_pos`. Used for `ref_ranges` in the column walker — every
/// `M` and `D` op consumes one ref position; this helper folds that
/// into the per-segment span.
pub(in crate::airr_record) fn extend_ref_range(ranges: &mut [Option<(i64, i64)>; 3], idx: usize, ref_pos: i64) {
    ranges[idx] = Some(match ranges[idx] {
        Some((s, e)) => (s.min(ref_pos), e.max(ref_pos + 1)),
        None => (ref_pos, ref_pos + 1),
    });
}

/// returns true when `ref_pos` is already inside
/// `ranges[idx]`. Used by the NP-claim path to detect when a
/// hypothesis's `pool_pos → ref_pos` projection collides with a
/// ref position the structural walker already accounted for. This
/// happens when the segment contains a structural-indel deletion
/// (the live-call hypothesis tracks pool↔ref linearly, but a
/// deletion breaks that mapping). When detected, the NP claim is
/// skipped — structural CIGAR ops take precedence over extension
/// ops on overlap.
pub(in crate::airr_record) fn ref_pos_already_covered(
    ranges: &[Option<(i64, i64)>; 3],
    idx: usize,
    ref_pos: i64,
) -> bool {
    match ranges[idx] {
        Some((s, e)) => ref_pos >= s && ref_pos < e,
        None => false,
    }
}
