use crate::address;
use crate::codon::translate_seq;
use crate::ir::Segment;
use crate::pass::Outcome;
use crate::refdata::RefDataConfig;

use super::junction::{
    aa_slice_for_region, anchor_amino_acid_preserved, anchor_pool_position, junction_has_stop,
};
use super::locus::{derive_locus, locus_from_refdata};
use super::projection::{
    lookup_allele, projected_call_name, unclaimed_np_bounds, unclaimed_np_string,
};
use super::sequence::{apply_rev_comp_projection, bytes_to_string, pool_bases};
use super::trace_fields::{trace_bool, trace_int};
use super::walk::walk_alignment_columns;
use super::AirrRecord;

/// Build an AIRR record from an outcome + refdata.
///
/// `sequence_id` is set on the record verbatim; pass an empty
/// string for unset.
pub fn build_airr_record(
    outcome: &Outcome,
    refdata: &RefDataConfig,
    sequence_id: &str,
) -> AirrRecord {
    let sim = outcome.final_simulation();
    let trace = &outcome.trace;

    let mut rec = AirrRecord::default();
    rec.sequence_id = sequence_id.to_string();
    rec.rev_comp = false;
    rec.c_call = String::new();

    // Bases as a string (preserving case from IR).
    rec.sequence = bytes_to_string(&pool_bases(sim));
    rec.sequence_length = rec.sequence.len() as i64;

    if rec.sequence.is_empty() {
        // Pre-recombine sim: nothing else to fill in.
        return rec;
    }

    // Trim values from trace. Our DSL records v_3, d_5, d_3, j_5;
    // v_5 and j_3 stay 0.
    rec.v_trim_5 = 0;
    rec.v_trim_3 = trace_int(trace, address::TRIM_V_3);
    rec.d_trim_5 = trace_int(trace, address::TRIM_D_5);
    rec.d_trim_3 = trace_int(trace, address::TRIM_D_3);
    rec.j_trim_5 = trace_int(trace, address::TRIM_J_5);
    rec.j_trim_3 = 0;

    // Calls + locus.
    let v_id = sim.assignments.v.map(|i| i.allele_id);
    let d_id = sim.assignments.d.map(|i| i.allele_id);
    let j_id = sim.assignments.j.map(|i| i.allele_id);

    rec.v_call = projected_call_name(refdata, sim, Segment::V, v_id);
    rec.d_call = projected_call_name(refdata, sim, Segment::D, d_id);
    rec.j_call = projected_call_name(refdata, sim, Segment::J, j_id);
    rec.locus = derive_locus(&rec.v_call, &rec.j_call, &rec.d_call);
    if rec.locus.is_empty() {
        // Under heavy corruption the live-call layer can wipe every
        // V/D/J call (no allele supports the mutated sequence),
        // leaving `derive_locus` with nothing to parse. Fall back
        // to the refdata's pool - any allele name gives us a locus
        // prefix, which is still meaningful AIRR output (the chain
        // didn't change, only the call evidence is absent).
        rec.locus = locus_from_refdata(refdata);
    }

    // Sequence-coord regions (raw pool start/end).
    let v_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V);
    let d_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::D);
    let j_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::J);
    let np1_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::Np1);
    let np2_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::Np2);

    // sequence coords come from the column-walker's
    // pool-range output, which honours the same overlap-resolution
    // rule the column relabel uses. The fallback is the structural
    // region (raw RefDataConfig runs, no live calls). Note: the
    // walker output itself is computed below; placeholder values
    // from the structural region are written here and then refined
    // post-walk to keep the data flow obvious.
    if let Some(r) = v_region {
        rec.v_sequence_start = Some(r.start.index() as i64);
        rec.v_sequence_end = Some(r.end.index() as i64);
    }
    if let Some(r) = d_region {
        rec.d_sequence_start = Some(r.start.index() as i64);
        rec.d_sequence_end = Some(r.end.index() as i64);
    }
    if let Some(r) = j_region {
        rec.j_sequence_start = Some(r.start.index() as i64);
        rec.j_sequence_end = Some(r.end.index() as i64);
    }

    // NP region nucleotide strings.
    //
    // when a V/D/J live-call extension claims one or more
    // NP positions (e.g. NP1[0..3] reabsorbed back into V), those
    // positions are no longer "non-templated" - they belong to a
    // segment. The AIRR `np1` / `np2` strings drop them.
    if let Some(r) = np1_region {
        let (s, n) = unclaimed_np_string(sim, refdata, &rec.sequence, r);
        rec.np1 = s;
        rec.np1_length = n;
    }
    if let Some(r) = np2_region {
        let (s, n) = unclaimed_np_string(sim, refdata, &rec.sequence, r);
        rec.np2 = s;
        rec.np2_length = n;
    }

    // Mutation + corruption counters. `n_mutations` reads from
    // `LiveCallState`, stashed by the S5F / Uniform passes at seal
    // time.
    rec.n_mutations = sim
        .live_calls
        .as_ref()
        .map(|s| s.mutation_count as i64)
        .unwrap_or(0);
    rec.mutation_rate = if rec.sequence_length > 0 {
        rec.n_mutations as f64 / rec.sequence_length as f64
    } else {
        0.0
    };
    rec.n_pcr_errors = trace_int(trace, address::CORRUPT_PCR_COUNT);
    rec.n_quality_errors = trace_int(trace, address::CORRUPT_QUALITY_COUNT);
    rec.n_indels = trace_int(trace, address::CORRUPT_INDEL_COUNT);
    rec.is_contaminant = trace_bool(trace, address::CORRUPT_CONTAMINANT_APPLIED);

    // Junction window: permissive semantics - compute as long as
    // both V/J anchors and regions exist, even when ``j_trim_5 >
    // j_anchor`` (which the strict `compute_junction` rejects).
    // Accept possibly-unphysical coordinates here; the
    // AnchorPreserved contract is the place to enforce strictness.
    //
    // anchor pool positions come from a `germline_pos`
    // scan over each segment's structural region rather than from
    // a `region.start + anchor` arithmetic offset. The scan finds
    // the actual pool position of the anchor codon's first base,
    // which is correct under indel passes (insertions / deletions
    // between region.start and the anchor shift the offset, but
    // each surviving germline node still carries its original
    // `germline_pos`).
    let v_anchor = lookup_allele(refdata, Segment::V, v_id).and_then(|a| a.anchor);
    let j_anchor = lookup_allele(refdata, Segment::J, j_id).and_then(|a| a.anchor);
    if let (Some(vr), Some(jr), Some(va), Some(ja)) = (v_region, j_region, v_anchor, j_anchor) {
        let v_anchor_pool = anchor_pool_position(sim, vr, va as u32);
        let j_anchor_pool = anchor_pool_position(sim, jr, ja as u32);
        if let (Some(vap), Some(jap)) = (v_anchor_pool, j_anchor_pool) {
            let v_anchor_in_pool: i64 = vap as i64;
            let j_anchor_in_pool: i64 = jap as i64;
            if j_anchor_in_pool + 3 > v_anchor_in_pool {
                let jstart = v_anchor_in_pool;
                let jend = j_anchor_in_pool + 3;
                // Slice into the assembled sequence using saturating
                // bounds - Python's slicing collapses out-of-range
                // negatives gracefully; we mirror that behaviour by
                // clamping to [0, sequence_length].
                let seq_len = rec.sequence_length;
                let safe_start = jstart.clamp(0, seq_len) as usize;
                let safe_end = jend.clamp(0, seq_len) as usize;
                let junction_nt: &str = if safe_end > safe_start {
                    &rec.sequence[safe_start..safe_end]
                } else {
                    ""
                };
                rec.junction = junction_nt.to_string();
                rec.junction_start = Some(jstart);
                rec.junction_end = Some(jend);
                // Length is the explicit nt length the Python builder
                // emitted, equal to len(junction_nt) (after slicing).
                rec.junction_length = Some(junction_nt.len() as i64);
                let in_frame = junction_nt.len() % 3 == 0;
                rec.vj_in_frame = Some(in_frame);
                if in_frame {
                    let has_stop = junction_has_stop(junction_nt);
                    let anchors_preserved = anchor_amino_acid_preserved(
                        sim,
                        refdata,
                        Segment::V,
                        vr,
                        v_id,
                        rec.v_trim_5,
                    ) && anchor_amino_acid_preserved(
                        sim,
                        refdata,
                        Segment::J,
                        jr,
                        j_id,
                        rec.j_trim_5,
                    );
                    rec.stop_codon = Some(has_stop);
                    rec.junction_aa = translate_seq(junction_nt);
                    rec.productive = Some(!has_stop && anchors_preserved);
                } else {
                    rec.stop_codon = Some(false);
                    rec.junction_aa = String::new();
                    rec.productive = Some(false);
                }
            }
        }
    }

    // sequence_aa / np_aa: V-anchor reading frame.
    //
    // `np1_aa` / `np2_aa` cover only the *unclaimed* span of the NP
    // region - the contiguous interior left after V's right-extension
    // and (for NP1) D's left-extension claim adjacent positions as
    // their own segment. Translating that span via `aa_slice_for_region`
    // gives an in-frame slice of `sequence_aa`, keeping `np2_aa`
    // consistent with both `np2_length` (unclaimed byte count) and
    // the codon-aligned position of NP inside the full sequence.
    if let Some(jstart) = rec.junction_start {
        let frame_offset = (jstart % 3) as usize;
        rec.sequence_aa = translate_seq(&rec.sequence[frame_offset..]);
        if let Some(r) = np1_region {
            if let Some((s, e)) = unclaimed_np_bounds(sim, refdata, r) {
                rec.np1_aa = aa_slice_for_region(s, e, frame_offset, &rec.sequence_aa);
            }
        }
        if let Some(r) = np2_region {
            if let Some((s, e)) = unclaimed_np_bounds(sim, refdata, r) {
                rec.np2_aa = aa_slice_for_region(s, e, frame_offset, &rec.sequence_aa);
            }
        }
    }

    // Single column walk: fills sequence_alignment, germline_alignment,
    // germline_alignment_d_mask, per-segment CIGARs, identity counts,
    // and alignment-coord pairs.
    let walk = walk_alignment_columns(sim, refdata, &rec);

    rec.sequence_alignment = walk.sa;
    rec.germline_alignment = walk.galn;
    rec.germline_alignment_d_mask = walk.dmask;

    rec.v_cigar = walk.cigars[0].clone();
    rec.d_cigar = walk.cigars[1].clone();
    rec.j_cigar = walk.cigars[2].clone();

    rec.v_alignment_start = walk.align_ranges[0].map(|(s, _)| s);
    rec.v_alignment_end = walk.align_ranges[0].map(|(_, e)| e);
    rec.d_alignment_start = walk.align_ranges[1].map(|(s, _)| s);
    rec.d_alignment_end = walk.align_ranges[1].map(|(_, e)| e);
    rec.j_alignment_start = walk.align_ranges[2].map(|(s, _)| s);
    rec.j_alignment_end = walk.align_ranges[2].map(|(_, e)| e);

    // refine sequence coords from the walker output.
    // The walker resolved any overlap between V/D/J live-call
    // hypotheses, so its `seq_ranges` reflects each segment's
    // *effective* span - guaranteeing CIGAR M+I+D == seq_len + D_ops.
    if let Some((s, e)) = walk.seq_ranges[0] {
        rec.v_sequence_start = Some(s);
        rec.v_sequence_end = Some(e);
    }
    if let Some((s, e)) = walk.seq_ranges[1] {
        rec.d_sequence_start = Some(s);
        rec.d_sequence_end = Some(e);
    }
    if let Some((s, e)) = walk.seq_ranges[2] {
        rec.j_sequence_start = Some(s);
        rec.j_sequence_end = Some(e);
    }

    // Germline coord pairs.
    //
    // read from the column walker's `ref_ranges` output
    // - the union of ref positions consumed by `M` and `D` ops on
    // each segment. This guarantees `germline_span == M + D` by
    // construction, closing the H.5 invariant gap that a raw
    // live-call-hypothesis read had (the hypothesis can extend past
    // what the column walker actually claims when an adjacent segment
    // wins an NP-position tiebreak).
    //
    // Fall back to the structural trim-derived range when there is
    // no allele assigned (raw RefDataConfig runs).
    if let Some((g_start, g_end)) = walk.ref_ranges[0] {
        rec.v_germline_start = Some(g_start);
        rec.v_germline_end = Some(g_end);
    } else if let Some(allele) = lookup_allele(refdata, Segment::V, v_id) {
        rec.v_germline_start = Some(rec.v_trim_5);
        rec.v_germline_end = Some(allele.seq.len() as i64 - rec.v_trim_3);
    }
    if let Some((g_start, g_end)) = walk.ref_ranges[1] {
        rec.d_germline_start = Some(g_start);
        rec.d_germline_end = Some(g_end);
    } else if let Some(allele) = lookup_allele(refdata, Segment::D, d_id) {
        rec.d_germline_start = Some(rec.d_trim_5);
        rec.d_germline_end = Some(allele.seq.len() as i64 - rec.d_trim_3);
    }
    if let Some((g_start, g_end)) = walk.ref_ranges[2] {
        rec.j_germline_start = Some(g_start);
        rec.j_germline_end = Some(g_end);
    } else if let Some(allele) = lookup_allele(refdata, Segment::J, j_id) {
        rec.j_germline_start = Some(rec.j_trim_5);
        rec.j_germline_end = Some(allele.seq.len() as i64 - rec.j_trim_3);
    }

    rec.v_identity = walk.identities[0];
    rec.d_identity = walk.identities[1];
    rec.j_identity = walk.identities[2];

    // reverse-complement projection. The RevCompPass
    // records `corrupt.rev_comp.applied` in the trace without
    // touching the IR (subsequent passes still see the forward
    // sequence). Here, *after* the forward record is fully built,
    // we flip the antisense-orientation fields per AIRR spec:
    //   - `sequence`, `np1`, `np2`, `junction` -> reverse complement
    //   - `*_sequence_start/end`, `junction_start/end` -> flipped
    //     so they index into the new (antisense) sequence
    //   - `*_aa` (sequence_aa, junction_aa, np1_aa, np2_aa) ->
    //     re-translated from the flipped strings
    //   - `sequence_alignment`, `germline_alignment`, CIGAR,
    //     `*_alignment_start/end`, `*_germline_start/end`,
    //     identity -> unchanged (forward orientation per spec)
    if trace_bool(trace, address::CORRUPT_REV_COMP_APPLIED) {
        apply_rev_comp_projection(&mut rec);
    }

    rec
}
