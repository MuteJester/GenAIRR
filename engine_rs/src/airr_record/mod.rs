//! AIRR Rearrangement record builder.
//!
//! Builds a fully-populated AIRR-format record from an `Outcome` and
//! its `RefDataConfig` in one walk. This is the Rust replacement for
//! the Python `_airr_record.outcome_to_airr_record` builder; the
//! Python side is now a thin wrapper around the PyO3 method this
//! module backs.
//!
//! Field semantics match exactly what the Python code produced.
//! Convention is **0-based half-open** for every coordinate field;
//! the `airr_strict=True` export flag in `result.py` does the
//! 1-based-inclusive conversion at TSV/CSV/DataFrame time.
//!
//! Walks the IR exactly once, accumulating all per-column data
//! (alignment strings, CIGAR ops, identity counts, segment spans)
//! in a single pass. The Python builder did four passes; this is
//! one of the main reasons we expect a >5× speedup at scale.

use crate::codon::{translate_codon_slice, translate_seq};
use crate::ir::{NucHandle, Region, Segment, Simulation};
#[cfg(test)]
use crate::ir::Nucleotide;
use crate::live_call::SegmentLiveCall;
use crate::pass::Outcome;
use crate::refdata::{Allele, RefDataConfig};
use crate::trace::{ChoiceValue, Trace};

mod walk;
use walk::walk_alignment_columns;

// Test-only re-import — `airr_record_tests.rs` pulls `runlength_to_string`
// through `use super::*`.
#[cfg(test)]
use walk::runlength_to_string;

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
        // Under heavy corruption the live-call layer can wipe every
        // V/D/J call (no allele supports the mutated sequence),
        // leaving `derive_locus` with nothing to parse. Fall back
        // to the refdata's pool — any allele name gives us a locus
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

    // Junction window: permissive semantics — compute as long as
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
    //
    // `np1_aa` / `np2_aa` cover only the *unclaimed* span of the NP
    // region — the contiguous interior left after V's right-extension
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
    // read from the column walker's `ref_ranges` output
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

    // reverse-complement projection. The RevCompPass
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
// Helpers
// ──────────────────────────────────────────────────────────────────

fn pool_bases(sim: &Simulation) -> Vec<u8> {
    sim.pool.as_slice().iter().map(|n| n.base).collect()
}

fn bytes_to_string(bytes: &[u8]) -> String {
    // Bases are ASCII; lossy decode is safe.
    String::from_utf8_lossy(bytes).into_owned()
}

/// locate the pool position where an allele's anchor
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
        if nuc.germline_pos.get().map(|p| p as u32) == Some(anchor_ref) {
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

    // locate the anchor by `germline_pos`, not by a
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
    translate_codon_slice(&live) == translate_codon_slice(&reference)
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
fn projected_allele_id(
    sim: &Simulation,
    segment: Segment,
) -> Option<crate::refdata::AlleleId> {
    // The sampled allele from recombination — used both as the
    // fallback when no live-call exists and as a preference among
    // alleles tied in the live-call's score-and-tie set. Several
    // alleles may share the same per-position score (the live-call
    // honestly reports the ambiguity in `v_call`) but for downstream
    // projection (germline_alignment, identity, CIGAR) we prefer the
    // truth allele when it sits inside the tie-set: it's the only
    // allele whose germline bytes really do match the pool bases. If
    // the truth isn't in the tie-set (mutations have shifted evidence
    // toward a different allele), fall back to the first tied id —
    // that's the genuine aligner-drift case the divergence narrative
    // captures.
    let truth_id = match segment {
        Segment::V => sim.assignments.v.map(|i| i.allele_id),
        Segment::D => sim.assignments.d.map(|i| i.allele_id),
        Segment::J => sim.assignments.j.map(|i| i.allele_id),
        _ => None,
    };
    if let Some(state) = &sim.live_calls {
        if let Some(call) = state.get(segment) {
            if let Some(tid) = truth_id {
                if call.allele_call.contains(tid) {
                    return Some(tid);
                }
            }
            if let Some(id) = call.allele_call.iter_ids().next() {
                return Some(id);
            }
        }
    }
    truth_id
}

fn projected_call_name(
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
        sim.live_calls.as_ref().and_then(|state| state.get(segment)),
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



/// build the AIRR `np1` / `np2` string from an NP
/// structural region by skipping any pool position that has been
/// claimed by a V/D/J live-call extension. The resulting string is
/// the "true" non-templated insertion as the evidence layer sees it.
/// pool-position bounds (`[start, end)`) of the contiguous unclaimed
/// interior of an NP region — the bytes that neither V's right
/// extension nor the adjacent segment's left extension claimed. Returns
/// `None` when every position in the region was claimed.
///
/// Used by `np1_aa` / `np2_aa` to slice `sequence_aa` at the codon
/// positions inside the unclaimed span, keeping the AA fields aligned
/// with `np_length` (which also counts only unclaimed bytes).
fn unclaimed_np_bounds(
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
    // NP-claim germline byte comes from the projected
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

/// refdata-driven locus fallback. Walks each pool's
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
#[path = "../airr_record_tests.rs"]
mod tests;
