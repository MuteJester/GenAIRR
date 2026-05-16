//! Single-pass alignment column walker for AIRR record building.
//!
//! `walk_alignment_columns` reads every region in the simulated
//! sequence once and produces every per-segment alignment artefact
//! the AIRR record needs: the sequence/germline/dmask alignment
//! strings, per-segment CIGAR strings, alignment column spans, pool
//! position spans, ref position spans, and per-segment identity
//! ratios. The Python builder did four passes over the IR for these
//! fields; this is one of the main reasons we expect a >5× speedup
//! at scale.

use super::{bytes_to_string, np_claim_owner, pool_bases, projected_allele_id, AirrRecord};
use crate::ir::{Nucleotide, Segment, Simulation};
use crate::refdata::RefDataConfig;

pub(super) struct AlignmentWalk {
    pub(super) sa: String,
    pub(super) galn: String,
    pub(super) dmask: String,
    /// One CIGAR string per V/D/J segment (in that order).
    pub(super) cigars: [String; 3],
    /// Per-segment alignment-string spans `(start, end)` (0-based
    /// half-open). `None` when the segment contributed no columns.
    pub(super) align_ranges: [Option<(i64, i64)>; 3],
    /// Per-segment pool-position spans `(start, end)` (0-based
    /// half-open). Tracks the structural region plus any NP columns
    /// the column-walker claimed for the segment — i.e. the slice of
    /// the sequence that this segment "owns" after live-call overlap
    /// resolution. Used for `*_sequence_start/end`.
    pub(super) seq_ranges: [Option<(i64, i64)>; 3],
    /// Per-segment allele-position spans `(start, end)` (0-based
    /// half-open). Mirrors `seq_ranges` but in reference space —
    /// the union of ref positions consumed by `M` and `D` ops on
    /// the segment. Used for `*_germline_start/end` so the AIRR
    /// record's germline span matches the CIGAR exactly:
    /// `germline_span == M + D` by construction.
    pub(super) ref_ranges: [Option<(i64, i64)>; 3],
    /// Per-segment identity (matches / total), or `None` when no
    /// columns. Indexed V=0, D=1, J=2.
    pub(super) identities: [Option<f64>; 3],
}

pub(super) fn walk_alignment_columns(
    sim: &Simulation,
    refdata: &RefDataConfig,
    rec: &AirrRecord,
) -> AlignmentWalk {
    // We read germline bytes directly off each `Nucleotide` inside the
    // walk loop rather than gathering the whole germline column up front
    // — saves an extra Vec allocation per record.
    let bases = pool_bases(sim);
    let pool_len = bases.len();

    // Grow incrementally to typical sequence size.
    let mut sa: Vec<u8> = Vec::with_capacity(pool_len + 8);
    let mut galn: Vec<u8> = Vec::with_capacity(pool_len + 8);
    let mut dmask: Vec<u8> = Vec::with_capacity(pool_len + 8);

    // Per-segment CIGAR run-length state.
    // [V, D, J] — accumulated `(count, op)` pairs for the segment.
    let mut cigar_runs: [Vec<(u32, u8)>; 3] = [Vec::new(), Vec::new(), Vec::new()];

    // Per-segment alignment-string spans (column offsets).
    let mut align_ranges: [Option<(i64, i64)>; 3] = [None, None, None];

    // Per-segment pool-position spans. Mirrors `align_ranges`, but
    // in pool coordinates so it's usable for AIRR `*_sequence_*`.
    let mut seq_ranges: [Option<(i64, i64)>; 3] = [None, None, None];

    // Per-segment ref-position spans. Mirrors `seq_ranges` but in
    // allele coordinates — the union of ref positions consumed by
    // `M` and `D` ops on the segment. Used for `*_germline_*`.
    let mut ref_ranges: [Option<(i64, i64)>; 3] = [None, None, None];

    // Per-segment identity counters [matches, total].
    let mut id_counts: [[u64; 2]; 3] = [[0, 0], [0, 0], [0, 0]];

    let mut covered: Vec<bool> = vec![false; pool_len];

    // Iterate regions in the order they appear in the sequence.
    for region in &sim.sequence.regions {
        let seg = region.segment;
        let r_start = region.start.index() as usize;
        let r_end = region.end.index() as usize;

        // Mark every pool position in the region as covered.
        for idx in r_start..r_end {
            covered[idx] = true;
        }

        match seg {
            Segment::V | Segment::D | Segment::J => {
                // pick the canonical allele from the live
                // call (with a provenance fallback). All `germline_*`
                // outputs in this branch read from this allele's
                // bytes, so they match the `*_call` field the AIRR
                // record reports — even when SHM narrowed the live
                // call to a different allele than originally sampled.
                let allele_id = projected_allele_id(sim, seg);
                let allele_seq: Option<&[u8]> = allele_id
                    .and_then(|aid| refdata.get(seg, aid))
                    .map(|a| a.seq.as_slice());
                let trim_5: i64 = match seg {
                    Segment::V => rec.v_trim_5,
                    Segment::D => rec.d_trim_5,
                    Segment::J => rec.j_trim_5,
                    _ => unreachable!(),
                };
                let trim_3: i64 = match seg {
                    Segment::V => rec.v_trim_3,
                    Segment::D => rec.d_trim_3,
                    Segment::J => rec.j_trim_3,
                    _ => unreachable!(),
                };
                let seg_idx = match seg {
                    Segment::V => 0,
                    Segment::D => 1,
                    Segment::J => 2,
                    _ => unreachable!(),
                };

                // The column-span start for this segment.
                let span_start = sa.len() as i64;

                if let Some(allele_seq) = allele_seq {
                    // NP1 left-extension might have already
                    // claimed ref positions up to (or past) `trim_5`.
                    // Start `expected_pos` past whatever NP claims
                    // covered, so a structural-indel deletion at the
                    // leftmost germline_pos doesn't trigger a D-fill
                    // at a ref position the NP claim already counted.
                    let np_covered_end = ref_ranges[seg_idx]
                        .map(|(_, e)| e as usize)
                        .unwrap_or(0);
                    let mut expected_pos: usize =
                        (trim_5.max(0) as usize).max(np_covered_end);
                    let end_germ: usize = (allele_seq.len() as i64 - trim_3).max(0) as usize;
                    for i in r_start..r_end {
                        let nuc: &Nucleotide = &sim.pool.as_slice()[i];
                        let base_char = bases[i];
                        let is_synthetic = nuc.germline_pos.is_none();
                        if is_synthetic {
                            // Indel insertion: gap in germline.
                            sa.push(base_char);
                            galn.push(b'-');
                            push_dmask_for_seg(&mut dmask, seg, b'-');
                            push_cigar_op(&mut cigar_runs[seg_idx], b'I');
                            id_counts[seg_idx][1] += 1;
                            // No match credit (germline gap).
                        } else {
                            // SAFETY: is_synthetic == false ⇒ germline_pos.is_some().
                            let germ_pos = nuc.germline_pos.get().expect("checked above") as usize;
                            // Fill any preceding deletion gap.
                            while expected_pos < germ_pos && expected_pos < allele_seq.len() {
                                sa.push(b'-');
                                let g = allele_seq[expected_pos];
                                galn.push(g);
                                push_dmask_for_seg(&mut dmask, seg, g);
                                push_cigar_op(&mut cigar_runs[seg_idx], b'D');
                                id_counts[seg_idx][1] += 1;
                                extend_ref_range(&mut ref_ranges, seg_idx, expected_pos as i64);
                                expected_pos += 1;
                            }
                            // Emit the matched/mismatched column. Phase
                            // 12.C: read the germline byte from the
                            // projected allele's sequence (matching
                            // `*_call`), not from `nuc.germline` which
                            // captures the originally-sampled byte and
                            // can diverge when the live-call narrowed
                            // to a different allele. Bounds-fall-back
                            // to `nuc.germline` if `germ_pos` exceeds
                            // the projected allele length (defensive,
                            // shouldn't happen in well-formed data).
                            sa.push(base_char);
                            let g = if germ_pos < allele_seq.len() {
                                allele_seq[germ_pos]
                            } else {
                                nuc.germline
                            };
                            galn.push(g);
                            push_dmask_for_seg(&mut dmask, seg, g);
                            push_cigar_op(&mut cigar_runs[seg_idx], b'M');
                            id_counts[seg_idx][1] += 1;
                            if eq_ascii_case_insensitive(base_char, g) {
                                id_counts[seg_idx][0] += 1;
                            }
                            extend_ref_range(&mut ref_ranges, seg_idx, germ_pos as i64);
                            expected_pos = germ_pos + 1;
                        }
                    }
                    // Trailing deletion gaps: positions in
                    // `[expected_pos, end_germ)` were trimmed off by
                    // an indel pass, not by recombination.
                    while expected_pos < end_germ && expected_pos < allele_seq.len() {
                        sa.push(b'-');
                        let g = allele_seq[expected_pos];
                        galn.push(g);
                        push_dmask_for_seg(&mut dmask, seg, g);
                        push_cigar_op(&mut cigar_runs[seg_idx], b'D');
                        id_counts[seg_idx][1] += 1;
                        extend_ref_range(&mut ref_ranges, seg_idx, expected_pos as i64);
                        expected_pos += 1;
                    }
                } else {
                    // No source allele assigned; emit raw bases with
                    // 'N' germline so lengths stay consistent.
                    for i in r_start..r_end {
                        sa.push(bases[i]);
                        galn.push(b'N');
                        push_dmask_for_seg(&mut dmask, seg, b'N');
                        // No CIGAR / identity entries — unaligned.
                    }
                }

                let span_end = sa.len() as i64;
                if span_end > span_start {
                    // a prior NP-side claim (e.g. J's
                    // left-extension into NP2) may have already set
                    // `align_ranges[seg_idx]` to a span that begins
                    // *before* this structural region. Preserve that
                    // earlier start instead of overwriting.
                    let new_start = match align_ranges[seg_idx] {
                        Some((s, _)) => s.min(span_start),
                        None => span_start,
                    };
                    align_ranges[seg_idx] = Some((new_start, span_end));
                    // Mirror in pool space.
                    let pool_start = r_start as i64;
                    let pool_end = r_end as i64;
                    let new_pool_start = match seq_ranges[seg_idx] {
                        Some((s, _)) => s.min(pool_start),
                        None => pool_start,
                    };
                    seq_ranges[seg_idx] = Some((new_pool_start, pool_end));
                }
            }
            Segment::Np1 | Segment::Np2 => {
                // NP positions claimed by an adjacent V/D/J live-call
                // extension are relabelled. The claimed column emits
                // the source allele's germline byte (instead of `N`),
                // pushes an `M` op onto that segment's CIGAR,
                // contributes to its identity counter, and extends its
                // `align_ranges`. Unclaimed positions remain plain NP
                // columns.
                //
                // a structural-indel deletion inside V/D/J
                // breaks the hypothesis's `pool_pos → ref_pos` linear
                // projection. When that projection points at a ref
                // position the structural walker already counted, we
                // drop the claim — and lock the segment out for the
                // rest of this NP region, since later positions in
                // the same projection would be similarly inconsistent.
                let mut np_blocked: [bool; 3] = [false; 3];
                for i in r_start..r_end {
                    let base_char = bases[i];
                    sa.push(base_char);
                    if let Some((claim_seg, allele_byte_proj, ref_pos_proj)) =
                        np_claim_owner(sim, refdata, i)
                    {
                        let claim_idx = match claim_seg {
                            Segment::V => 0,
                            Segment::D => 1,
                            Segment::J => 2,
                            _ => unreachable!(),
                        };
                        // pick the canonical ref_pos for
                        // this claim. For an NP2-side extension (the
                        // segment's structural region has already been
                        // walked, so `ref_ranges[claim_idx]` is set),
                        // use `ref_ranges.end` as the next ref slot —
                        // this sidesteps the hypothesis's linear
                        // pool→ref projection, which mis-shifts under
                        // synthetic insertions in the structural
                        // region. For NP1-side extensions
                        // (`ref_ranges[claim_idx]` is still empty),
                        // keep the hypothesis projection — NP regions
                        // themselves don't carry indels, so the linear
                        // formula is exact there.
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
                        if np_blocked[claim_idx]
                            || ref_pos_already_covered(&ref_ranges, claim_idx, ref_pos as i64)
                        {
                            np_blocked[claim_idx] = true;
                            galn.push(b'N');
                            dmask.push(b'N');
                            continue;
                        }
                        galn.push(allele_byte);
                        push_dmask_for_seg(&mut dmask, claim_seg, allele_byte);
                        push_cigar_op(&mut cigar_runs[claim_idx], b'M');
                        id_counts[claim_idx][1] += 1;
                        if eq_ascii_case_insensitive(base_char, allele_byte) {
                            id_counts[claim_idx][0] += 1;
                        }
                        extend_ref_range(&mut ref_ranges, claim_idx, ref_pos as i64);
                        let col = (sa.len() as i64) - 1;
                        align_ranges[claim_idx] = Some(match align_ranges[claim_idx] {
                            Some((s, _)) => (s, col + 1),
                            None => (col, col + 1),
                        });
                        let pool_pos = i as i64;
                        seq_ranges[claim_idx] = Some(match seq_ranges[claim_idx] {
                            Some((s, e)) => (s.min(pool_pos), e.max(pool_pos + 1)),
                            None => (pool_pos, pool_pos + 1),
                        });
                    } else {
                        galn.push(b'N');
                        dmask.push(b'N');
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

    // Uppercase strings for AIRR-tooling compatibility.
    bytes_uppercase_in_place(&mut sa);
    bytes_uppercase_in_place(&mut galn);
    bytes_uppercase_in_place(&mut dmask);

    let cigars = [
        runlength_to_string(&cigar_runs[0]),
        runlength_to_string(&cigar_runs[1]),
        runlength_to_string(&cigar_runs[2]),
    ];

    let identities = [
        if id_counts[0][1] > 0 {
            Some(id_counts[0][0] as f64 / id_counts[0][1] as f64)
        } else {
            None
        },
        if id_counts[1][1] > 0 {
            Some(id_counts[1][0] as f64 / id_counts[1][1] as f64)
        } else {
            None
        },
        if id_counts[2][1] > 0 {
            Some(id_counts[2][0] as f64 / id_counts[2][1] as f64)
        } else {
            None
        },
    ];

    AlignmentWalk {
        sa: bytes_to_string(&sa),
        galn: bytes_to_string(&galn),
        dmask: bytes_to_string(&dmask),
        cigars,
        align_ranges,
        seq_ranges,
        ref_ranges,
        identities,
    }
}

fn bytes_uppercase_in_place(bytes: &mut [u8]) {
    for b in bytes.iter_mut() {
        if (*b).is_ascii_lowercase() {
            *b = (*b).to_ascii_uppercase();
        }
    }
}

fn eq_ascii_case_insensitive(a: u8, b: u8) -> bool {
    a.to_ascii_uppercase() == b.to_ascii_uppercase()
}

fn push_cigar_op(runs: &mut Vec<(u32, u8)>, op: u8) {
    if let Some(last) = runs.last_mut() {
        if last.1 == op {
            last.0 += 1;
            return;
        }
    }
    runs.push((1, op));
}

pub(super) fn runlength_to_string(runs: &[(u32, u8)]) -> String {
    let mut s = String::with_capacity(runs.len() * 4);
    for (count, op) in runs {
        s.push_str(&count.to_string());
        s.push(*op as char);
    }
    s
}

fn push_dmask_for_seg(dmask: &mut Vec<u8>, seg: Segment, ga_char: u8) {
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
fn extend_ref_range(ranges: &mut [Option<(i64, i64)>; 3], idx: usize, ref_pos: i64) {
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
fn ref_pos_already_covered(
    ranges: &[Option<(i64, i64)>; 3],
    idx: usize,
    ref_pos: i64,
) -> bool {
    match ranges[idx] {
        Some((s, e)) => ref_pos >= s && ref_pos < e,
        None => false,
    }
}
