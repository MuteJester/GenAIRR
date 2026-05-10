//! AIRR Rearrangement record builder.
//!
//! Builds a fully-populated AIRR-format record from an `Outcome` and
//! its `RefDataConfig` in one walk. This is the Rust replacement for
//! the Python `_airr_record.outcome_to_airr_record` builder; the
//! Python side is now a thin wrapper around the PyO3 method this
//! module backs.
//!
//! Field semantics match exactly what the Python code produced (the
//! Phase H tests are the regression suite). Convention is **0-based
//! half-open** for every coordinate field; the `airr_strict=True`
//! export flag in `result.py` does the 1-based-inclusive conversion
//! at TSV/CSV/DataFrame time.
//!
//! Walks the IR exactly once, accumulating all per-column data
//! (alignment strings, CIGAR ops, identity counts, segment spans)
//! in a single pass. The Python builder did four passes; this is
//! one of the main reasons we expect a >5× speedup at scale.

use crate::ir::{NucHandle, Nucleotide, Region, Segment, Simulation};
use crate::live_call::SegmentLiveCall;
use crate::pass::Outcome;
use crate::refdata::{Allele, RefDataConfig};
use crate::trace::{ChoiceValue, Trace};

// ──────────────────────────────────────────────────────────────────
// Public record struct
// ──────────────────────────────────────────────────────────────────

/// One AIRR Rearrangement record. All ~50 fields populated as
/// ground truth from the IR, refdata, and trace — no aligner.
///
/// `Option<i64>` for coordinate fields lets us emit `null` /
/// missing values for VJ chains' D coords or pre-recombine sims.
/// Strings are owned `String`s; per-record allocation cost is
/// dwarfed by the actual sequence string lengths anyway.
#[derive(Debug, Clone, Default)]
pub struct AirrRecord {
    // AIRR metadata
    pub sequence_id: String,
    pub sequence: String,
    pub sequence_aa: String,
    pub sequence_alignment: String,
    pub germline_alignment: String,
    pub germline_alignment_d_mask: String,
    pub sequence_length: i64,
    pub rev_comp: bool,
    pub locus: String,

    // V
    pub v_call: String,
    pub v_cigar: String,
    pub v_score: Option<f64>,
    pub v_identity: Option<f64>,
    pub v_support: Option<f64>,
    pub v_sequence_start: Option<i64>,
    pub v_sequence_end: Option<i64>,
    pub v_alignment_start: Option<i64>,
    pub v_alignment_end: Option<i64>,
    pub v_germline_start: Option<i64>,
    pub v_germline_end: Option<i64>,
    pub v_trim_5: i64,
    pub v_trim_3: i64,

    // D
    pub d_call: String,
    pub d_cigar: String,
    pub d_score: Option<f64>,
    pub d_identity: Option<f64>,
    pub d_support: Option<f64>,
    pub d_sequence_start: Option<i64>,
    pub d_sequence_end: Option<i64>,
    pub d_alignment_start: Option<i64>,
    pub d_alignment_end: Option<i64>,
    pub d_germline_start: Option<i64>,
    pub d_germline_end: Option<i64>,
    pub d_trim_5: i64,
    pub d_trim_3: i64,

    // J
    pub j_call: String,
    pub j_cigar: String,
    pub j_score: Option<f64>,
    pub j_identity: Option<f64>,
    pub j_support: Option<f64>,
    pub j_sequence_start: Option<i64>,
    pub j_sequence_end: Option<i64>,
    pub j_alignment_start: Option<i64>,
    pub j_alignment_end: Option<i64>,
    pub j_germline_start: Option<i64>,
    pub j_germline_end: Option<i64>,
    pub j_trim_5: i64,
    pub j_trim_3: i64,

    // C
    pub c_call: String,

    // Junction
    pub junction: String,
    pub junction_aa: String,
    pub junction_start: Option<i64>,
    pub junction_end: Option<i64>,
    pub junction_length: Option<i64>,

    // NP regions
    pub np1: String,
    pub np1_aa: String,
    pub np1_length: i64,
    pub np2: String,
    pub np2_aa: String,
    pub np2_length: i64,

    // Functionality
    pub productive: Option<bool>,
    pub vj_in_frame: Option<bool>,
    pub stop_codon: Option<bool>,

    // SHM + corruption (non-AIRR; GenAIRR additions)
    pub n_mutations: i64,
    pub mutation_rate: f64,
    pub n_pcr_errors: i64,
    pub n_quality_errors: i64,
    pub n_indels: i64,
    pub is_contaminant: bool,
}

// ──────────────────────────────────────────────────────────────────
// Builder entry point
// ──────────────────────────────────────────────────────────────────

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
    rec.v_trim_3 = trace_int(trace, "trim.v_3");
    rec.d_trim_5 = trace_int(trace, "trim.d_5");
    rec.d_trim_3 = trace_int(trace, "trim.d_3");
    rec.j_trim_5 = trace_int(trace, "trim.j_5");
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
        // Phase 12.C follow-up: under heavy corruption the live-call
        // layer can wipe every V/D/J call (no allele supports the
        // mutated sequence), leaving `derive_locus` with nothing to
        // parse. Fall back to the refdata's pool — any allele name
        // gives us a locus prefix, which is still meaningful AIRR
        // output (the chain didn't change, only the call evidence
        // is absent).
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

    // Phase 11.4: sequence coords come from the column-walker's
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
    // Phase 11.4: when a V/D/J live-call extension claims one or more
    // NP positions (e.g. NP1[0..3] reabsorbed back into V), those
    // positions are no longer "non-templated" — they belong to a
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

    // Mutation + corruption counters.
    rec.n_mutations = mutation_count(trace);
    rec.mutation_rate = if rec.sequence_length > 0 {
        rec.n_mutations as f64 / rec.sequence_length as f64
    } else {
        0.0
    };
    rec.n_pcr_errors = trace_int(trace, "corrupt.pcr.count");
    rec.n_quality_errors = trace_int(trace, "corrupt.quality.count");
    rec.n_indels = trace_int(trace, "corrupt.indel.count");
    rec.is_contaminant = trace_bool(trace, "corrupt.contaminant.applied");

    // Junction window: mirror the Python builder's permissive
    // semantics — compute as long as both V/J anchors and regions
    // exist, even when ``j_trim_5 > j_anchor`` (which the strict
    // `compute_junction` rejects). Keeping bit-parity with the
    // Phase H Python builder means accepting possibly-unphysical
    // coordinates here; the AnchorPreserved contract is the place
    // to enforce strictness.
    //
    // Phase 11.6: anchor pool positions come from a `germline_pos`
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
            // bounds — Python's slicing collapses out-of-range
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
                let anchors_preserved =
                    anchor_amino_acid_preserved(sim, refdata, Segment::V, vr, v_id, rec.v_trim_5)
                        && anchor_amino_acid_preserved(
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
    if let Some(jstart) = rec.junction_start {
        let frame_offset = (jstart % 3) as usize;
        rec.sequence_aa = translate_seq(&rec.sequence[frame_offset..]);
        if let Some(r) = np1_region {
            rec.np1_aa = aa_slice_for_region(
                r.start.index() as usize,
                r.end.index() as usize,
                frame_offset,
                &rec.sequence_aa,
            );
        }
        if let Some(r) = np2_region {
            rec.np2_aa = aa_slice_for_region(
                r.start.index() as usize,
                r.end.index() as usize,
                frame_offset,
                &rec.sequence_aa,
            );
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

    // Phase 11.4: refine sequence coords from the walker output.
    // The walker resolved any overlap between V/D/J live-call
    // hypotheses, so its `seq_ranges` reflects each segment's
    // *effective* span — guaranteeing CIGAR M+I+D == seq_len + D_ops.
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
    // Phase 11.7: read from the column walker's `ref_ranges` output
    // — the union of ref positions consumed by `M` and `D` ops on
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

    // Phase 12.D: reverse-complement projection. The RevCompPass
    // records `corrupt.rev_comp.applied` in the trace without
    // touching the IR (subsequent passes still see the forward
    // sequence). Here, *after* the forward record is fully built,
    // we flip the antisense-orientation fields per AIRR spec:
    //   - `sequence`, `np1`, `np2`, `junction` → reverse complement
    //   - `*_sequence_start/end`, `junction_start/end` → flipped
    //     so they index into the new (antisense) sequence
    //   - `*_aa` (sequence_aa, junction_aa, np1_aa, np2_aa) →
    //     re-translated from the flipped strings
    //   - `sequence_alignment`, `germline_alignment`, CIGAR,
    //     `*_alignment_start/end`, `*_germline_start/end`,
    //     identity → unchanged (forward orientation per spec)
    if trace_bool(trace, "corrupt.rev_comp.applied") {
        apply_rev_comp_projection(&mut rec);
    }

    rec
}

fn apply_rev_comp_projection(rec: &mut AirrRecord) {
    rec.rev_comp = true;
    let seq_len = rec.sequence_length;
    rec.sequence = reverse_complement(&rec.sequence);
    rec.np1 = reverse_complement(&rec.np1);
    rec.np2 = reverse_complement(&rec.np2);
    rec.junction = reverse_complement(&rec.junction);

    // Flip pool-position coords: new_start = seq_len - old_end,
    // new_end = seq_len - old_start.
    flip_coord_pair(&mut rec.v_sequence_start, &mut rec.v_sequence_end, seq_len);
    flip_coord_pair(&mut rec.d_sequence_start, &mut rec.d_sequence_end, seq_len);
    flip_coord_pair(&mut rec.j_sequence_start, &mut rec.j_sequence_end, seq_len);
    flip_coord_pair(&mut rec.junction_start, &mut rec.junction_end, seq_len);

    // Re-translate AA fields from the flipped strings.
    rec.junction_aa = translate_seq(&rec.junction);
    if let Some(jstart) = rec.junction_start {
        let frame_offset = (jstart % 3) as usize;
        rec.sequence_aa = translate_seq(&rec.sequence[frame_offset..]);
        // np1_aa / np2_aa come from `aa_slice_for_region`, but the
        // np region pool offsets aren't carried on the AirrRecord
        // (only the strings + lengths are). Re-derive by walking
        // the new sequence: find the np1 / np2 substring positions
        // and slice. Simpler: just translate each string in its
        // own frame, which is what the forward path effectively did
        // for these short np regions when the frame happened to
        // align. AIRR consumers that care use `junction_aa`; np_aa
        // is a convenience field. Translate plain:
        rec.np1_aa = translate_seq(&rec.np1);
        rec.np2_aa = translate_seq(&rec.np2);
    } else {
        rec.sequence_aa = String::new();
        rec.np1_aa = String::new();
        rec.np2_aa = String::new();
    }
}

fn flip_coord_pair(start: &mut Option<i64>, end: &mut Option<i64>, total_len: i64) {
    if let (Some(s), Some(e)) = (*start, *end) {
        *start = Some(total_len - e);
        *end = Some(total_len - s);
    }
}

fn reverse_complement(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars().rev() {
        out.push(match c {
            'A' => 'T',
            'T' => 'A',
            'U' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'u' => 'a',
            'c' => 'g',
            'g' => 'c',
            'N' => 'N',
            'n' => 'n',
            other => other,
        });
    }
    out
}

// ──────────────────────────────────────────────────────────────────
// Single-pass column walker
// ──────────────────────────────────────────────────────────────────

/// Result of one walk over the alignment columns.
struct AlignmentWalk {
    sa: String,
    galn: String,
    dmask: String,
    /// One CIGAR string per V/D/J segment (in that order).
    cigars: [String; 3],
    /// Per-segment alignment-string spans `(start, end)` (0-based
    /// half-open). `None` when the segment contributed no columns.
    align_ranges: [Option<(i64, i64)>; 3],
    /// Per-segment pool-position spans `(start, end)` (0-based
    /// half-open). Tracks the structural region plus any NP columns
    /// the column-walker claimed for the segment — i.e. the slice of
    /// the sequence that this segment "owns" after live-call overlap
    /// resolution. Used for `*_sequence_start/end`.
    seq_ranges: [Option<(i64, i64)>; 3],
    /// Per-segment allele-position spans `(start, end)` (0-based
    /// half-open). Mirrors `seq_ranges` but in reference space —
    /// the union of ref positions consumed by `M` and `D` ops on
    /// the segment. Used for `*_germline_start/end` so the AIRR
    /// record's germline span matches the CIGAR exactly:
    /// `germline_span == M + D` by construction.
    ref_ranges: [Option<(i64, i64)>; 3],
    /// Per-segment identity (matches / total), or `None` when no
    /// columns. Indexed V=0, D=1, J=2.
    identities: [Option<f64>; 3],
}

fn walk_alignment_columns(
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
                // Phase 12.C: pick the canonical allele from the live
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
                    // Phase 12.C: NP1 left-extension might have already
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
                        let is_synthetic = nuc.germline_pos == Nucleotide::NO_GERMLINE_POS;
                        if is_synthetic {
                            // Indel insertion: gap in germline.
                            sa.push(base_char);
                            galn.push(b'-');
                            push_dmask_for_seg(&mut dmask, seg, b'-');
                            push_cigar_op(&mut cigar_runs[seg_idx], b'I');
                            id_counts[seg_idx][1] += 1;
                            // No match credit (germline gap).
                        } else {
                            let germ_pos = nuc.germline_pos as usize;
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
                    // Phase 11.4: a prior NP-side claim (e.g. J's
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
                // Phase 11.3 + 11.4 + 11.5: NP positions claimed by an
                // adjacent V/D/J live-call extension are relabelled.
                // The claimed column emits the source allele's germline
                // byte (instead of `N`), pushes an `M` op onto that
                // segment's CIGAR, contributes to its identity counter,
                // and extends its `align_ranges`. Unclaimed positions
                // remain plain NP columns.
                //
                // Phase 11.7: a structural-indel deletion inside V/D/J
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
                        // Phase 12.C: pick the canonical ref_pos for
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

// ──────────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────────

fn pool_bases(sim: &Simulation) -> Vec<u8> {
    sim.pool.as_slice().iter().map(|n| n.base).collect()
}

fn bytes_to_string(bytes: &[u8]) -> String {
    // Bases are ASCII; lossy decode is safe.
    String::from_utf8_lossy(bytes).into_owned()
}

/// Phase 11.6: locate the pool position where an allele's anchor
/// codon resides in the assembled sequence. Scans the segment's
/// structural region for the first node whose `germline_pos`
/// matches the anchor's allele coordinate. Returns `None` when
/// the anchor was trimmed off, deleted by a structural-indel
/// pass, or otherwise not present in the live data.
///
/// Biologically correct under indels: each surviving germline
/// node carries its original `germline_pos`, so the anchor codon
/// can be found wherever the actual base ended up — even after
/// insertions or deletions shifted things around inside the
/// segment.
fn anchor_pool_position(
    sim: &Simulation,
    region: &Region,
    anchor_ref: u32,
) -> Option<u32> {
    let r_start = region.start.index();
    let r_end = region.end.index();
    let pool = sim.pool.as_slice();
    for i in r_start..r_end {
        let nuc = &pool[i as usize];
        if nuc.germline_pos != Nucleotide::NO_GERMLINE_POS
            && nuc.germline_pos as u32 == anchor_ref
        {
            return Some(i);
        }
    }
    None
}

fn anchor_amino_acid_preserved(
    sim: &Simulation,
    refdata: &RefDataConfig,
    segment: Segment,
    region: &Region,
    allele_id: Option<crate::refdata::AlleleId>,
    _trim_5: i64,
) -> bool {
    let Some(allele) = lookup_allele(refdata, segment, allele_id) else {
        return true;
    };
    let Some(anchor) = allele.anchor else {
        return true;
    };

    let anchor = anchor as i64;
    if anchor < 0 || anchor + 3 > allele.seq.len() as i64 {
        return false;
    }

    // Phase 11.6: locate the anchor by `germline_pos`, not by a
    // structural offset. Same rationale as the junction-window
    // scanner — a structural offset miscomputes the anchor pool
    // position when indels disturb V/J's coding region.
    let Some(pool_start) = anchor_pool_position(sim, region, anchor as u32) else {
        return false;
    };
    if pool_start + 3 > region.end.index() {
        return false;
    }

    let mut live = [b'N'; 3];
    for offset in 0..3 {
        let Some(nuc) = sim.pool.get(NucHandle::new(pool_start + offset)) else {
            return false;
        };
        live[offset as usize] = nuc.base;
    }

    let ref_start = anchor as usize;
    let reference = [
        allele.seq[ref_start],
        allele.seq[ref_start + 1],
        allele.seq[ref_start + 2],
    ];
    translate_codon(&live) == translate_codon(&reference)
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

fn runlength_to_string(runs: &[(u32, u8)]) -> String {
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

fn lookup_name(
    refdata: &RefDataConfig,
    segment: Segment,
    id: Option<crate::refdata::AlleleId>,
) -> String {
    id.and_then(|aid| refdata.get(segment, aid))
        .map(|a| a.name.clone())
        .unwrap_or_default()
}

/// Phase 12.C: pick the canonical allele id for AIRR projection on
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
fn projected_allele_id(
    sim: &Simulation,
    segment: Segment,
) -> Option<crate::refdata::AlleleId> {
    if let Some(state) = &sim.live_calls {
        if let Some(call) = state.get(segment) {
            if let Some(id) = call.allele_call.iter_ids().next() {
                return Some(id);
            }
        }
    }
    match segment {
        Segment::V => sim.assignments.v.map(|i| i.allele_id),
        Segment::D => sim.assignments.d.map(|i| i.allele_id),
        Segment::J => sim.assignments.j.map(|i| i.allele_id),
        _ => None,
    }
}

fn projected_call_name(
    refdata: &RefDataConfig,
    sim: &Simulation,
    segment: Segment,
    origin_id: Option<crate::refdata::AlleleId>,
) -> String {
    live_call_name(
        refdata,
        sim.live_calls.as_ref().and_then(|state| state.get(segment)),
        segment,
    )
    .unwrap_or_else(|| lookup_name(refdata, segment, origin_id))
}

/// Phase 11.3 / 11.4: figure out which V/D/J segment (if any) has
/// an extension hypothesis that claims this NP pool position.
///
/// Returns `Some((segment, germline_byte))` when an adjacent V/D/J
/// live-call hypothesis extended past its structural region into
/// `pool_pos`. The germline byte comes from the segment's source
/// allele at the inferred reference position. Returns `None` when
/// no V/D/J claim exists — the position stays an NP column with `N`
/// in `germline_alignment`.
///
/// Tie-breaking when both V (right extension) and D (left extension)
/// claim the same NP1 position: V wins. Same for D vs J on NP2.
/// This matches the natural V→D→J order and is deterministic.
fn np_claim_owner(
    sim: &Simulation,
    refdata: &RefDataConfig,
    pool_pos: usize,
) -> Option<(Segment, u8, u32)> {
    for seg in [Segment::V, Segment::D, Segment::J] {
        if let Some(claim) = check_segment_claim(sim, refdata, seg, pool_pos) {
            return Some(claim);
        }
    }
    None
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

/// Phase 11.7: returns true when `ref_pos` is already inside
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



/// Phase 11.4: build the AIRR `np1` / `np2` string from an NP
/// structural region by skipping any pool position that has been
/// claimed by a V/D/J live-call extension. The resulting string is
/// the "true" non-templated insertion as the evidence layer sees it.
fn unclaimed_np_string(
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
    let call = sim.live_calls.as_ref()?.get(seg)?;
    let h = call.hypotheses.first()?;
    // Phase 12.C: NP-claim germline byte comes from the projected
    // allele (live-call's first allele, fallback to provenance) so
    // it matches `*_call`. See `projected_allele_id` for rationale.
    let allele_id = projected_allele_id(sim, seg)?;
    let allele = refdata.get(seg, allele_id)?;
    let region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == seg)?;
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
) -> Option<String> {
    let call = call?;
    let names: Vec<String> = call
        .allele_call
        .iter_ids()
        .filter_map(|id| refdata.get(segment, id).map(|allele| allele.name.clone()))
        .collect();
    Some(names.join(","))
}

fn lookup_allele<'a>(
    refdata: &'a RefDataConfig,
    segment: Segment,
    id: Option<crate::refdata::AlleleId>,
) -> Option<&'a Allele> {
    id.and_then(|aid| refdata.get(segment, aid))
}

fn trace_int(trace: &Trace, address: &str) -> i64 {
    match trace.find(address) {
        Some(rec) => match rec.value {
            ChoiceValue::Int(v) => v,
            _ => 0,
        },
        None => 0,
    }
}

fn trace_bool(trace: &Trace, address: &str) -> bool {
    match trace.find(address) {
        Some(rec) => match rec.value {
            ChoiceValue::Bool(v) => v,
            _ => false,
        },
        None => false,
    }
}

fn mutation_count(trace: &Trace) -> i64 {
    if let Some(r) = trace.find("mutate.s5f.count") {
        if let ChoiceValue::Int(v) = r.value {
            return v;
        }
    }
    if let Some(r) = trace.find("mutate.uniform.count") {
        if let ChoiceValue::Int(v) = r.value {
            return v;
        }
    }
    0
}

const AIRR_LOCI: [&str; 7] = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"];

/// Phase 12.C: refdata-driven locus fallback. Walks each pool's
/// first allele in turn and returns the locus prefix when the name
/// starts with one of the AIRR loci. Used when live-call evidence
/// has wiped every `*_call` (heavy SHM under corruption stack) but
/// the chain is still well-defined by the source data.
fn locus_from_refdata(refdata: &RefDataConfig) -> String {
    for entry in [
        refdata.v_pool.iter().next(),
        refdata.j_pool.iter().next(),
        refdata.d_pool.iter().next(),
    ]
    .into_iter()
    .flatten()
    {
        let candidate = locus_prefix(&entry.1.name);
        if !candidate.is_empty() {
            return candidate;
        }
    }
    String::new()
}

fn locus_prefix(name: &str) -> String {
    if name.len() < 3 {
        return String::new();
    }
    let mut prefix = String::with_capacity(3);
    for c in name.chars().take(3) {
        prefix.push(c.to_ascii_uppercase());
    }
    if AIRR_LOCI.contains(&prefix.as_str()) {
        prefix
    } else {
        String::new()
    }
}

fn derive_locus(v_call: &str, j_call: &str, d_call: &str) -> String {
    for name in [v_call, j_call, d_call] {
        if name.is_empty() {
            continue;
        }
        let candidate = locus_prefix(name);
        if !candidate.is_empty() {
            return candidate;
        }
    }
    String::new()
}

// ── Translation ────────────────────────────────────────────────────

const GENETIC_CODE: &[(&[u8; 3], char)] = &[
    (b"TTT", 'F'),
    (b"TTC", 'F'),
    (b"TTA", 'L'),
    (b"TTG", 'L'),
    (b"CTT", 'L'),
    (b"CTC", 'L'),
    (b"CTA", 'L'),
    (b"CTG", 'L'),
    (b"ATT", 'I'),
    (b"ATC", 'I'),
    (b"ATA", 'I'),
    (b"ATG", 'M'),
    (b"GTT", 'V'),
    (b"GTC", 'V'),
    (b"GTA", 'V'),
    (b"GTG", 'V'),
    (b"TCT", 'S'),
    (b"TCC", 'S'),
    (b"TCA", 'S'),
    (b"TCG", 'S'),
    (b"AGT", 'S'),
    (b"AGC", 'S'),
    (b"CCT", 'P'),
    (b"CCC", 'P'),
    (b"CCA", 'P'),
    (b"CCG", 'P'),
    (b"ACT", 'T'),
    (b"ACC", 'T'),
    (b"ACA", 'T'),
    (b"ACG", 'T'),
    (b"GCT", 'A'),
    (b"GCC", 'A'),
    (b"GCA", 'A'),
    (b"GCG", 'A'),
    (b"TAT", 'Y'),
    (b"TAC", 'Y'),
    (b"TAA", '*'),
    (b"TAG", '*'),
    (b"TGA", '*'),
    (b"CAT", 'H'),
    (b"CAC", 'H'),
    (b"CAA", 'Q'),
    (b"CAG", 'Q'),
    (b"AAT", 'N'),
    (b"AAC", 'N'),
    (b"AAA", 'K'),
    (b"AAG", 'K'),
    (b"GAT", 'D'),
    (b"GAC", 'D'),
    (b"GAA", 'E'),
    (b"GAG", 'E'),
    (b"TGT", 'C'),
    (b"TGC", 'C'),
    (b"TGG", 'W'),
    (b"CGT", 'R'),
    (b"CGC", 'R'),
    (b"CGA", 'R'),
    (b"CGG", 'R'),
    (b"AGA", 'R'),
    (b"AGG", 'R'),
    (b"GGT", 'G'),
    (b"GGC", 'G'),
    (b"GGA", 'G'),
    (b"GGG", 'G'),
];

fn translate_codon(codon: &[u8]) -> char {
    if codon.len() != 3 {
        return 'X';
    }
    let upper = [
        match codon[0] {
            b'U' | b'u' => b'T',
            c => c.to_ascii_uppercase(),
        },
        match codon[1] {
            b'U' | b'u' => b'T',
            c => c.to_ascii_uppercase(),
        },
        match codon[2] {
            b'U' | b'u' => b'T',
            c => c.to_ascii_uppercase(),
        },
    ];
    for (k, v) in GENETIC_CODE {
        if k[0] == upper[0] && k[1] == upper[1] && k[2] == upper[2] {
            return *v;
        }
    }
    'X'
}

fn translate_seq(seq: &str) -> String {
    let bytes = seq.as_bytes();
    let n_codons = bytes.len() / 3;
    let mut out = String::with_capacity(n_codons);
    for i in 0..n_codons {
        let start = i * 3;
        out.push(translate_codon(&bytes[start..start + 3]));
    }
    out
}

fn aa_slice_for_region(
    region_start: usize,
    region_end: usize,
    frame_offset: usize,
    sequence_aa: &str,
) -> String {
    if region_end <= region_start {
        return String::new();
    }
    // First codon-aligned position at or after `region_start`. Mirror
    // the Python `(frame_offset - region_start) % 3` semantics with
    // `rem_euclid` so a negative arithmetic difference still produces
    // a non-negative remainder in `[0, 3)`.
    let delta = (frame_offset as i64 - region_start as i64).rem_euclid(3) as usize;
    let codon_start = region_start + delta;
    if codon_start >= region_end {
        return String::new();
    }
    let n_codons = (region_end - codon_start) / 3;
    if n_codons == 0 {
        return String::new();
    }
    let aa_index = (codon_start - frame_offset) / 3;
    let aa_bytes = sequence_aa.as_bytes();
    if aa_index >= aa_bytes.len() {
        return String::new();
    }
    let end = (aa_index + n_codons).min(aa_bytes.len());
    String::from_utf8_lossy(&aa_bytes[aa_index..end]).into_owned()
}

fn junction_has_stop(seq: &str) -> bool {
    let bytes = seq.as_bytes();
    let n = bytes.len() / 3;
    for i in 0..n {
        let start = i * 3;
        let mut codon = [0u8; 3];
        for j in 0..3 {
            codon[j] = match bytes[start + j] {
                b'U' | b'u' => b'T',
                c => c.to_ascii_uppercase(),
            };
        }
        if (codon == *b"TAA") || (codon == *b"TAG") || (codon == *b"TGA") {
            return true;
        }
    }
    false
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assignment::AlleleInstance;
    use crate::live_call::{
        AlleleBitSet, EvidenceScore, HypothesisFlags, LiveCallState, PlacementHypothesis,
        SegmentLiveCall,
    };
    use crate::refdata::{Allele, AlleleId, ChainType};
    use crate::trace::Trace;

    #[test]
    fn translate_codon_basic() {
        assert_eq!(translate_codon(b"ATG"), 'M');
        assert_eq!(translate_codon(b"taa"), '*');
        assert_eq!(translate_codon(b"NNN"), 'X');
    }

    #[test]
    fn translate_seq_truncates_partial_codon() {
        assert_eq!(translate_seq("ATGCCC"), "MP");
        assert_eq!(translate_seq("ATGCCCA"), "MP"); // partial codon dropped
    }

    #[test]
    fn derive_locus_picks_first_known_prefix() {
        assert_eq!(derive_locus("IGHV1-2*01", "", ""), "IGH");
        assert_eq!(derive_locus("", "IGKJ4*01", ""), "IGK");
        assert_eq!(derive_locus("ighv1-2*01", "", ""), "IGH"); // case-insensitive
        assert_eq!(derive_locus("", "", ""), "");
        assert_eq!(derive_locus("XYZ1*01", "", ""), "");
    }

    #[test]
    fn locus_from_refdata_falls_back_to_pool_allele_names() {
        // Phase 12.C: when every live call has been wiped, the
        // refdata's pool allele names still tell us the locus.
        let mut cfg = RefDataConfig::empty(ChainType::Vdj);
        let _ = cfg.v_pool.push(Allele {
            name: "IGHV1-2*01".into(),
            gene: "IGHV1-2".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        assert_eq!(locus_from_refdata(&cfg), "IGH");

        // Non-AIRR-prefixed names should yield empty (the helper
        // requires a recognisable AIRR locus).
        let mut alien = RefDataConfig::empty(ChainType::Vdj);
        let _ = alien.v_pool.push(Allele {
            name: "XYZV1*01".into(),
            gene: "XYZV1".into(),
            seq: b"AAA".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        assert_eq!(locus_from_refdata(&alien), "");

        // Empty pool → empty locus.
        let empty = RefDataConfig::empty(ChainType::Vj);
        assert_eq!(locus_from_refdata(&empty), "");
    }

    #[test]
    fn runlength_collapses_repeated_ops() {
        let runs = vec![(5, b'M'), (2, b'I'), (3, b'M'), (1, b'D')];
        assert_eq!(runlength_to_string(&runs), "5M2I3M1D");
    }

    #[test]
    fn runlength_empty_is_empty_string() {
        assert_eq!(runlength_to_string(&[]), "");
    }

    #[test]
    fn aa_slice_for_region_basic() {
        // sequence_aa "ABCDE", frame_offset 0, region [3, 12) → "BCD".
        assert_eq!(aa_slice_for_region(3, 12, 0, "ABCDE"), "BCD");
        // Region [4, 12) starts mid-codon → first complete codon at 6.
        assert_eq!(aa_slice_for_region(4, 12, 0, "ABCDE"), "CD");
        // Empty region.
        assert_eq!(aa_slice_for_region(5, 5, 0, "ABCDE"), "");
        // Region too short for a complete codon.
        assert_eq!(aa_slice_for_region(4, 5, 0, "ABCDE"), "");
    }

    #[test]
    fn junction_has_stop_detects_stops() {
        assert!(junction_has_stop("TAA"));
        assert!(junction_has_stop("ATGTAA"));
        assert!(junction_has_stop("ATGTAG"));
        assert!(junction_has_stop("ATGTGA"));
        assert!(!junction_has_stop("ATG"));
        assert!(!junction_has_stop("ATGCCC"));
    }

    fn anchor_record_fixture() -> (RefDataConfig, Simulation) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_airr*01".into(),
            gene: "v_airr".into(),
            seq: b"TGT".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });

        let mut sim = Simulation::new();
        for (i, &b) in b"TGT".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(b, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(3))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim.with_region_added(region).with_allele_assigned(
            Segment::V,
            crate::assignment::AlleleInstance::new(AlleleId::new(0)),
        );

        (cfg, sim)
    }

    fn outcome_from_sim(sim: Simulation) -> Outcome {
        Outcome {
            revisions: vec![sim],
            pass_names: Vec::new(),
            trace: Trace::new(),
            events: Vec::new(),
        }
    }

    fn call_projection_fixture() -> (RefDataConfig, Simulation, AlleleId, AlleleId) {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let v0 = cfg.v_pool.push(Allele {
            name: "IGHV1-1*01".into(),
            gene: "IGHV1-1".into(),
            seq: b"AAACCC".to_vec(),
            segment: Segment::V,
            anchor: None,
        });
        let v1 = cfg.v_pool.push(Allele {
            name: "IGHV1-1*02".into(),
            gene: "IGHV1-1".into(),
            seq: b"AAACCC".to_vec(),
            segment: Segment::V,
            anchor: None,
        });

        let mut sim = Simulation::new();
        for (i, &base) in b"AAACCC".iter().enumerate() {
            let (next, _) =
                sim.with_nucleotide_pushed(Nucleotide::germline(base, i as u16, Segment::V));
            sim = next;
        }
        let region = Region::new(Segment::V, NucHandle::new(0), NucHandle::new(6))
            .with_codon_rail_recomputed(&sim.pool);
        sim = sim
            .with_region_added(region)
            .with_allele_assigned(Segment::V, AlleleInstance::new(v0));

        (cfg, sim, v0, v1)
    }

    #[test]
    fn airr_call_projection_falls_back_to_origin_assignment_without_live_call() {
        let (cfg, sim, _v0, _v1) = call_projection_fixture();
        let outcome = outcome_from_sim(sim);
        let rec = build_airr_record(&outcome, &cfg, "fallback");

        assert_eq!(rec.v_call, "IGHV1-1*01");
        assert_eq!(rec.locus, "IGH");
    }

    #[test]
    fn airr_call_projection_prefers_live_call_allele_set() {
        let (cfg, sim, v0, v1) = call_projection_fixture();
        let hypothesis = PlacementHypothesis::new(
            Segment::V,
            0,
            6,
            0,
            6,
            AlleleBitSet::from_ids(cfg.v_pool.len(), [v0, v1]),
            EvidenceScore::exact(6, 0),
            HypothesisFlags::EMPTY,
        );
        let call =
            SegmentLiveCall::from_hypotheses(Segment::V, cfg.v_pool.len(), vec![hypothesis], 1);
        let live = LiveCallState::empty().with_segment_call(call);
        let outcome = outcome_from_sim(sim.with_live_calls(live));

        let rec = build_airr_record(&outcome, &cfg, "live");

        assert_eq!(rec.v_call, "IGHV1-1*01,IGHV1-1*02");
        assert_eq!(rec.locus, "IGH");
    }

    #[test]
    fn airr_call_projection_uses_empty_call_when_live_call_is_unsupported() {
        let (cfg, sim, _v0, _v1) = call_projection_fixture();
        let call = SegmentLiveCall::from_hypotheses(Segment::V, cfg.v_pool.len(), Vec::new(), 1);
        let live = LiveCallState::empty().with_segment_call(call);
        let outcome = outcome_from_sim(sim.with_live_calls(live));

        let rec = build_airr_record(&outcome, &cfg, "unsupported");

        assert_eq!(rec.v_call, "");
    }

    #[test]
    fn anchor_amino_acid_preserved_allows_synonymous_change() {
        let (cfg, sim) = anchor_record_fixture();
        let changed = sim.with_base_changed(NucHandle::new(2), b'C');
        let region = changed
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)
            .unwrap();

        assert!(anchor_amino_acid_preserved(
            &changed,
            &cfg,
            Segment::V,
            region,
            Some(AlleleId::new(0)),
            0,
        ));
    }

    #[test]
    fn anchor_amino_acid_preserved_rejects_nonsynonymous_change() {
        let (cfg, sim) = anchor_record_fixture();
        let changed = sim.with_base_changed(NucHandle::new(0), b'A');
        let region = changed
            .sequence
            .regions
            .iter()
            .find(|r| r.segment == Segment::V)
            .unwrap();

        assert!(!anchor_amino_acid_preserved(
            &changed,
            &cfg,
            Segment::V,
            region,
            Some(AlleleId::new(0)),
            0,
        ));
    }
}
