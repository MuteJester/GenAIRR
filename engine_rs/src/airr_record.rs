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

use crate::ir::{Nucleotide, Segment, Simulation};
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

    rec.v_call = lookup_name(refdata, Segment::V, v_id);
    rec.d_call = lookup_name(refdata, Segment::D, d_id);
    rec.j_call = lookup_name(refdata, Segment::J, j_id);
    rec.locus = derive_locus(&rec.v_call, &rec.j_call, &rec.d_call);

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

    rec.v_sequence_start = v_region.map(|r| r.start.index() as i64);
    rec.v_sequence_end = v_region.map(|r| r.end.index() as i64);
    rec.d_sequence_start = d_region.map(|r| r.start.index() as i64);
    rec.d_sequence_end = d_region.map(|r| r.end.index() as i64);
    rec.j_sequence_start = j_region.map(|r| r.start.index() as i64);
    rec.j_sequence_end = j_region.map(|r| r.end.index() as i64);

    // NP region nucleotide strings.
    if let Some(r) = np1_region {
        rec.np1 = rec.sequence[r.start.index() as usize..r.end.index() as usize].to_string();
        rec.np1_length = (r.end.index() - r.start.index()) as i64;
    }
    if let Some(r) = np2_region {
        rec.np2 = rec.sequence[r.start.index() as usize..r.end.index() as usize].to_string();
        rec.np2_length = (r.end.index() - r.start.index()) as i64;
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
    let v_anchor = lookup_allele(refdata, Segment::V, v_id).and_then(|a| a.anchor);
    let j_anchor = lookup_allele(refdata, Segment::J, j_id).and_then(|a| a.anchor);
    if let (Some(vr), Some(jr), Some(va), Some(ja)) = (v_region, j_region, v_anchor, j_anchor) {
        let v_anchor_in_pool: i64 = vr.start.index() as i64 + (va as i64 - 0); // v_trim_5 always 0 in our DSL
        let j_anchor_in_pool: i64 = jr.start.index() as i64 + (ja as i64 - rec.j_trim_5);
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
                rec.stop_codon = Some(has_stop);
                rec.junction_aa = translate_seq(junction_nt);
                rec.productive = Some(!has_stop);
            } else {
                rec.stop_codon = Some(false);
                rec.junction_aa = String::new();
                rec.productive = Some(false);
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

    // Germline coord pairs: derived from trim values + allele length.
    if let Some(allele) = lookup_allele(refdata, Segment::V, v_id) {
        if rec.v_alignment_start.is_some() {
            rec.v_germline_start = Some(rec.v_trim_5);
            rec.v_germline_end = Some(allele.seq.len() as i64 - rec.v_trim_3);
        }
    }
    if let Some(allele) = lookup_allele(refdata, Segment::D, d_id) {
        if rec.d_alignment_start.is_some() {
            rec.d_germline_start = Some(rec.d_trim_5);
            rec.d_germline_end = Some(allele.seq.len() as i64 - rec.d_trim_3);
        }
    }
    if let Some(allele) = lookup_allele(refdata, Segment::J, j_id) {
        if rec.j_alignment_start.is_some() {
            rec.j_germline_start = Some(rec.j_trim_5);
            rec.j_germline_end = Some(allele.seq.len() as i64 - rec.j_trim_3);
        }
    }

    rec.v_identity = walk.identities[0];
    rec.d_identity = walk.identities[1];
    rec.j_identity = walk.identities[2];

    rec
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
                let allele_id = match seg {
                    Segment::V => sim.assignments.v.map(|i| i.allele_id),
                    Segment::D => sim.assignments.d.map(|i| i.allele_id),
                    Segment::J => sim.assignments.j.map(|i| i.allele_id),
                    _ => unreachable!(),
                };
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
                    let mut expected_pos: usize = trim_5.max(0) as usize;
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
                                expected_pos += 1;
                            }
                            // Emit the matched/mismatched column.
                            sa.push(base_char);
                            let g = nuc.germline;
                            galn.push(g);
                            push_dmask_for_seg(&mut dmask, seg, g);
                            push_cigar_op(&mut cigar_runs[seg_idx], b'M');
                            id_counts[seg_idx][1] += 1;
                            if eq_ascii_case_insensitive(base_char, g) {
                                id_counts[seg_idx][0] += 1;
                            }
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
                    align_ranges[seg_idx] = Some((span_start, span_end));
                }
            }
            Segment::Np1 | Segment::Np2 => {
                for i in r_start..r_end {
                    sa.push(bases[i]);
                    galn.push(b'N');
                    dmask.push(b'N');
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

fn derive_locus(v_call: &str, j_call: &str, d_call: &str) -> String {
    for name in [v_call, j_call, d_call] {
        if name.is_empty() {
            continue;
        }
        if name.len() < 3 {
            continue;
        }
        let mut prefix = String::with_capacity(3);
        for c in name.chars().take(3) {
            prefix.push(c.to_ascii_uppercase());
        }
        if AIRR_LOCI.contains(&prefix.as_str()) {
            return prefix;
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
}
