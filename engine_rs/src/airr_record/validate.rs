//! AIRR-record postcondition validator.
//!
//! Given a final `(Simulation, Outcome, RefDataConfig)` and the
//! `AirrRecord` produced from it, this module re-derives every
//! field from independent sources and reports any divergence as a
//! `RecordValidationIssue`. The validator is **read-only**: it
//! never mutates the record, the outcome, or refdata.
//!
//! See [`docs/airr_record_validator.md`](../../../docs/airr_record_validator.md)
//! for the full check catalogue and the design discussion.

use crate::address::{ChoiceAddress, PrimeEnd, VdjSegment};
use crate::ir::{Segment, SimulationEvent};
use crate::live_call::scoring::{
    allele_pool_for_segment, score_alleles_with_extensions, tie_set_ids_at_max_score,
};
use crate::pass::Outcome;
use crate::refdata::RefDataConfig;
use crate::trace::ChoiceValue;

use super::junction::{anchor_amino_acid_preserved, anchor_pool_position, junction_has_stop};
use super::projection::lookup_allele;
use super::record::AirrRecord;
use super::sequence::pool_bases;

/// One reported issue. Variants are stable identifiers; downstream
/// CI / dashboards can pattern-match on them.
#[derive(Clone, Debug, PartialEq)]
pub enum RecordValidationIssue {
    // ── C1: Structural record invariants ────────────────────────
    SequenceLengthMismatch {
        reported: i64,
        actual_bytes: usize,
    },
    SequenceContentMismatch {
        reported_prefix: String,
        actual_prefix: String,
    },
    SegmentCoordinatesOutOfOrder {
        segment: Segment,
        start: i64,
        end: i64,
    },
    GermlineCoordinatesOutOfOrder {
        segment: Segment,
        start: i64,
        end: i64,
    },
    CigarReadsInvalid {
        segment: Segment,
        reason: String,
    },
    CigarSpanMismatch {
        segment: Segment,
        cigar_query_span: usize,
        sequence_span: usize,
    },

    // ── C2: Counter provenance ──────────────────────────────────
    NMutationsMismatch {
        reported: i64,
        sim_count: i64,
    },
    NPcrErrorsMismatch {
        reported: i64,
        trace_count: i64,
    },
    NQualityErrorsMismatch {
        reported: i64,
        trace_count: i64,
    },
    NIndelsMismatch {
        reported: i64,
        event_count: i64,
    },
    NSegmentIndelsMismatch {
        segment: Segment,
        reported: i64,
        event_count: i64,
    },
    EndLossLengthMismatch {
        side: PrimeEnd,
        reported: i64,
        trace_count: i64,
    },

    // ── C3: Junction truth ─────────────────────────────────────
    JunctionLengthMismatch {
        reported: Option<i64>,
        recomputed: u32,
    },
    JunctionContentMismatch {
        reported: String,
        recomputed: String,
    },
    JunctionAaMismatch {
        reported: String,
        recomputed: String,
    },
    VjInFrameMismatch {
        reported: Option<bool>,
        recomputed: bool,
    },
    StopCodonMismatch {
        reported: Option<bool>,
        recomputed: bool,
    },
    ProductiveMismatch {
        reported: Option<bool>,
        recomputed: bool,
        reason: ProductiveDecidedBy,
    },

    // ── C4: Allele-call oracle ─────────────────────────────────
    AlleleCallTieSetMismatch {
        segment: Segment,
        reported: Vec<String>,
        recomputed: Vec<String>,
    },
    AlleleCallOrderMismatch {
        segment: Segment,
        reported_first: String,
        expected_first: String,
        reason: AlleleOrderReason,
    },

    // ── C5: Region / live-call structural invariants ───────────
    MultipleRegionsForSegment {
        segment: Segment,
        count: usize,
    },
    MultipleHypothesesInLiveCall {
        segment: Segment,
        count: usize,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub enum ProductiveDecidedBy {
    OutOfFrame,
    JunctionStopCodon,
    VAnchorAaChanged,
    JAnchorAaChanged,
    InFrameAndAnchorsPreserved,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AlleleOrderReason {
    TruthFirstIfInTieSet,
    AscendingAlleleIdOtherwise,
}

/// Validate the AIRR record against the simulation outcome that
/// produced it. Returns an empty vector when the record passes
/// every check.
pub fn validate_airr_record(
    record: &AirrRecord,
    outcome: &Outcome,
    refdata: &RefDataConfig,
) -> Vec<RecordValidationIssue> {
    let mut issues = Vec::new();
    let sim = outcome.final_simulation();

    check_structural(record, sim, &mut issues);
    check_counters(record, outcome, &mut issues);
    check_junction(record, sim, refdata, &mut issues);
    check_allele_oracle(record, outcome, refdata, &mut issues);
    check_region_and_hypothesis_invariants(sim, &mut issues);

    issues
}

// ──────────────────────────────────────────────────────────────────
// C1: Structural record invariants
// ──────────────────────────────────────────────────────────────────

fn check_structural(
    record: &AirrRecord,
    sim: &crate::ir::Simulation,
    issues: &mut Vec<RecordValidationIssue>,
) {
    // sequence_length matches actual sequence string bytes.
    let actual_bytes = record.sequence.len();
    if record.sequence_length as usize != actual_bytes {
        issues.push(RecordValidationIssue::SequenceLengthMismatch {
            reported: record.sequence_length,
            actual_bytes,
        });
    }

    // sequence content matches the pool's bases (case-folded).
    // Compare case-insensitively because sequencing errors lowercase
    // their substitutions, and the AIRR record may carry rev-comp.
    let pool_seq = String::from_utf8(pool_bases(sim)).unwrap_or_default();
    if !record.rev_comp && !record.sequence.eq_ignore_ascii_case(&pool_seq) {
        let prefix_len = record.sequence.len().min(pool_seq.len()).min(40);
        issues.push(RecordValidationIssue::SequenceContentMismatch {
            reported_prefix: record.sequence.chars().take(prefix_len).collect(),
            actual_prefix: pool_seq.chars().take(prefix_len).collect(),
        });
    }

    // Coordinate ordering for V/D/J segments and germline ranges.
    check_coord_pair(
        Segment::V,
        record.v_sequence_start,
        record.v_sequence_end,
        false,
        issues,
    );
    check_coord_pair(
        Segment::V,
        record.v_germline_start,
        record.v_germline_end,
        true,
        issues,
    );
    check_coord_pair(
        Segment::D,
        record.d_sequence_start,
        record.d_sequence_end,
        false,
        issues,
    );
    check_coord_pair(
        Segment::D,
        record.d_germline_start,
        record.d_germline_end,
        true,
        issues,
    );
    check_coord_pair(
        Segment::J,
        record.j_sequence_start,
        record.j_sequence_end,
        false,
        issues,
    );
    check_coord_pair(
        Segment::J,
        record.j_germline_start,
        record.j_germline_end,
        true,
        issues,
    );

    // CIGAR span sanity: M+I ops must equal the sequence-side span.
    for (seg, start, end, cigar) in [
        (
            Segment::V,
            record.v_sequence_start,
            record.v_sequence_end,
            &record.v_cigar,
        ),
        (
            Segment::D,
            record.d_sequence_start,
            record.d_sequence_end,
            &record.d_cigar,
        ),
        (
            Segment::J,
            record.j_sequence_start,
            record.j_sequence_end,
            &record.j_cigar,
        ),
    ] {
        if cigar.is_empty() {
            continue;
        }
        match parse_cigar(cigar) {
            Ok(ops) => {
                let query_span: usize = ops
                    .iter()
                    .filter(|(_, op)| *op == b'M' || *op == b'I')
                    .map(|(n, _)| *n)
                    .sum();
                if let (Some(s), Some(e)) = (start, end) {
                    let sequence_span = (e - s).max(0) as usize;
                    if query_span != sequence_span {
                        issues.push(RecordValidationIssue::CigarSpanMismatch {
                            segment: seg,
                            cigar_query_span: query_span,
                            sequence_span,
                        });
                    }
                }
            }
            Err(reason) => issues.push(RecordValidationIssue::CigarReadsInvalid {
                segment: seg,
                reason,
            }),
        }
    }
}

fn check_coord_pair(
    segment: Segment,
    start: Option<i64>,
    end: Option<i64>,
    is_germline: bool,
    issues: &mut Vec<RecordValidationIssue>,
) {
    if let (Some(s), Some(e)) = (start, end) {
        // Half-open: end > start is required (or both 0 for empty).
        if s < 0 || e < s {
            let issue = if is_germline {
                RecordValidationIssue::GermlineCoordinatesOutOfOrder {
                    segment,
                    start: s,
                    end: e,
                }
            } else {
                RecordValidationIssue::SegmentCoordinatesOutOfOrder {
                    segment,
                    start: s,
                    end: e,
                }
            };
            issues.push(issue);
        }
    }
}

fn parse_cigar(cigar: &str) -> Result<Vec<(usize, u8)>, String> {
    let mut ops = Vec::new();
    let mut digits = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            digits.push(ch);
        } else {
            let n: usize = digits
                .parse()
                .map_err(|_| format!("non-numeric op length in CIGAR {cigar:?}"))?;
            if !matches!(ch as u8, b'M' | b'I' | b'D' | b'S' | b'N' | b'P' | b'X' | b'=') {
                return Err(format!("unrecognized CIGAR op {ch:?} in {cigar:?}"));
            }
            ops.push((n, ch as u8));
            digits.clear();
        }
    }
    if !digits.is_empty() {
        return Err(format!("trailing digits in CIGAR {cigar:?}"));
    }
    Ok(ops)
}

// ──────────────────────────────────────────────────────────────────
// C2: Counter provenance
// ──────────────────────────────────────────────────────────────────

fn check_counters(
    record: &AirrRecord,
    outcome: &Outcome,
    issues: &mut Vec<RecordValidationIssue>,
) {
    let sim = outcome.final_simulation();

    // n_mutations comes from sim.mutation_count (set by S5F / Uniform
    // at seal). Direct equality.
    if record.n_mutations != sim.mutation_count as i64 {
        issues.push(RecordValidationIssue::NMutationsMismatch {
            reported: record.n_mutations,
            sim_count: sim.mutation_count as i64,
        });
    }

    // n_pcr_errors / n_quality_errors from trace.
    let trace_int = |addr: ChoiceAddress| -> i64 {
        outcome
            .trace
            .find_choice(addr)
            .and_then(|rec| match rec.value {
                ChoiceValue::Int(n) => Some(n),
                _ => None,
            })
            .unwrap_or(0)
    };

    let pcr_attempts = trace_int(ChoiceAddress::CorruptPcrCount);
    if record.n_pcr_errors != pcr_attempts {
        issues.push(RecordValidationIssue::NPcrErrorsMismatch {
            reported: record.n_pcr_errors,
            trace_count: pcr_attempts,
        });
    }

    let quality_attempts = trace_int(ChoiceAddress::CorruptQualityCount);
    if record.n_quality_errors != quality_attempts {
        issues.push(RecordValidationIssue::NQualityErrorsMismatch {
            reported: record.n_quality_errors,
            trace_count: quality_attempts,
        });
    }

    // n_indels and per-segment indel counters from event ledger of the
    // corrupt.indel pass only (per audit §6.1 / §6.2).
    let mut total = 0i64;
    let mut by_segment = [0i64; Segment::COUNT];
    for er in outcome.events() {
        if er.pass_name != crate::address::CORRUPT_INDEL {
            continue;
        }
        for ev in &er.simulation_events {
            let segment = match ev {
                SimulationEvent::IndelInserted { segment, .. } => *segment,
                SimulationEvent::IndelDeleted { segment, .. } => *segment,
                _ => continue,
            };
            total += 1;
            by_segment[segment as usize] += 1;
        }
    }
    if record.n_indels != total {
        issues.push(RecordValidationIssue::NIndelsMismatch {
            reported: record.n_indels,
            event_count: total,
        });
    }
    for (seg, reported) in [
        (Segment::V, record.n_v_indels),
        (Segment::D, record.n_d_indels),
        (Segment::J, record.n_j_indels),
    ] {
        let expected = by_segment[seg as usize];
        if reported != expected {
            issues.push(RecordValidationIssue::NSegmentIndelsMismatch {
                segment: seg,
                reported,
                event_count: expected,
            });
        }
    }

    // End-loss lengths from trace.
    let el5 = trace_int(ChoiceAddress::CorruptEndLoss(PrimeEnd::Five));
    if record.end_loss_5_length != el5 {
        issues.push(RecordValidationIssue::EndLossLengthMismatch {
            side: PrimeEnd::Five,
            reported: record.end_loss_5_length,
            trace_count: el5,
        });
    }
    let el3 = trace_int(ChoiceAddress::CorruptEndLoss(PrimeEnd::Three));
    if record.end_loss_3_length != el3 {
        issues.push(RecordValidationIssue::EndLossLengthMismatch {
            side: PrimeEnd::Three,
            reported: record.end_loss_3_length,
            trace_count: el3,
        });
    }
}

// ──────────────────────────────────────────────────────────────────
// C3: Junction truth
// ──────────────────────────────────────────────────────────────────

fn check_junction(
    record: &AirrRecord,
    sim: &crate::ir::Simulation,
    refdata: &RefDataConfig,
    issues: &mut Vec<RecordValidationIssue>,
) {
    // Junction recomputation only meaningful for non-rev-comp records;
    // rev-comp flips coordinates and re-translates AA after the
    // junction is sliced. Skip when record is reverse-complemented —
    // the §5 audit covers that path with dedicated tests.
    if record.rev_comp {
        return;
    }

    // Re-derive the junction window the same way builder.rs does:
    // germline_pos scan over V/J regions, not offset arithmetic.
    // (The contract-side `crate::junction::compute_junction` uses
    // offsets and diverges under indels/end-loss; the AIRR builder
    // uses the scan, so the validator must match the builder.)
    let v_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V);
    let j_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::J);
    let v_id = sim.assignments.get(Segment::V).map(|i| i.allele_id);
    let j_id = sim.assignments.get(Segment::J).map(|i| i.allele_id);
    let v_anchor = lookup_allele(refdata, Segment::V, v_id).and_then(|a| a.anchor);
    let j_anchor = lookup_allele(refdata, Segment::J, j_id).and_then(|a| a.anchor);

    let (Some(vr), Some(jr), Some(va), Some(ja)) = (v_region, j_region, v_anchor, j_anchor) else {
        // Junction undefined; builder leaves fields as defaults.
        return;
    };
    let (Some(vap), Some(jap)) = (
        anchor_pool_position(sim, vr, va as u32),
        anchor_pool_position(sim, jr, ja as u32),
    ) else {
        return;
    };
    let v_anchor_in_pool = vap as i64;
    let j_anchor_in_pool = jap as i64;
    if j_anchor_in_pool + 3 <= v_anchor_in_pool {
        return;
    }
    let jstart = v_anchor_in_pool;
    let jend = j_anchor_in_pool + 3;
    let seq_len = record.sequence_length;
    let safe_start = jstart.clamp(0, seq_len) as usize;
    let safe_end = jend.clamp(0, seq_len) as usize;
    let recomputed_content: String = if safe_end > safe_start {
        record.sequence[safe_start..safe_end].to_string()
    } else {
        String::new()
    };
    let recomputed_len = recomputed_content.len() as u32;

    let reported_len = record.junction_length;
    if reported_len != Some(recomputed_len as i64) {
        issues.push(RecordValidationIssue::JunctionLengthMismatch {
            reported: reported_len,
            recomputed: recomputed_len,
        });
    }

    if !record.junction.eq_ignore_ascii_case(&recomputed_content) {
        issues.push(RecordValidationIssue::JunctionContentMismatch {
            reported: record.junction.clone(),
            recomputed: recomputed_content.clone(),
        });
    }

    // Frame.
    let recomputed_in_frame = recomputed_len % 3 == 0;
    if record.vj_in_frame != Some(recomputed_in_frame) {
        issues.push(RecordValidationIssue::VjInFrameMismatch {
            reported: record.vj_in_frame,
            recomputed: recomputed_in_frame,
        });
    }

    // Stop codon (only meaningful when in-frame).
    let recomputed_stop = recomputed_in_frame && junction_has_stop(&record.junction);
    if record.stop_codon != Some(recomputed_stop) {
        issues.push(RecordValidationIssue::StopCodonMismatch {
            reported: record.stop_codon,
            recomputed: recomputed_stop,
        });
    }

    // Productive triad: in-frame ∧ no stop ∧ V/J anchor amino acids preserved.
    let v_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::V);
    let j_region = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == Segment::J);
    let v_anchor_ok = v_region
        .map(|r| {
            anchor_amino_acid_preserved(
                sim,
                refdata,
                Segment::V,
                r,
                sim.assignments.get(Segment::V).map(|i| i.allele_id),
                record.v_trim_5,
            )
        })
        .unwrap_or(true);
    let j_anchor_ok = j_region
        .map(|r| {
            anchor_amino_acid_preserved(
                sim,
                refdata,
                Segment::J,
                r,
                sim.assignments.get(Segment::J).map(|i| i.allele_id),
                record.j_trim_5,
            )
        })
        .unwrap_or(true);

    let (recomputed_productive, reason) = if !recomputed_in_frame {
        (false, ProductiveDecidedBy::OutOfFrame)
    } else if recomputed_stop {
        (false, ProductiveDecidedBy::JunctionStopCodon)
    } else if !v_anchor_ok {
        (false, ProductiveDecidedBy::VAnchorAaChanged)
    } else if !j_anchor_ok {
        (false, ProductiveDecidedBy::JAnchorAaChanged)
    } else {
        (true, ProductiveDecidedBy::InFrameAndAnchorsPreserved)
    };
    if record.productive != Some(recomputed_productive) {
        issues.push(RecordValidationIssue::ProductiveMismatch {
            reported: record.productive,
            recomputed: recomputed_productive,
            reason,
        });
    }
}

// ──────────────────────────────────────────────────────────────────
// C4: Allele-call oracle
//
// Independent reimplementation of the walker's max-match-count
// tie-set selection. We walk the segment's region bytes, count
// matches per allele, and pick the alleles at the max score.
// Then we re-derive the projected CSV order (truth first if in
// tie-set, otherwise ascending allele id).
// ──────────────────────────────────────────────────────────────────

fn check_allele_oracle(
    record: &AirrRecord,
    outcome: &Outcome,
    refdata: &RefDataConfig,
    issues: &mut Vec<RecordValidationIssue>,
) {
    let sim = outcome.final_simulation();
    if record.rev_comp {
        // Rev-comp flips the sequence; the allele call was computed
        // pre-flip. Skip the oracle here; the rev-comp audit covers
        // it with dedicated tests.
        return;
    }
    for (segment, vdj, reported_call) in [
        (Segment::V, VdjSegment::V, &record.v_call),
        (Segment::D, VdjSegment::D, &record.d_call),
        (Segment::J, VdjSegment::J, &record.j_call),
    ] {
        let _ = vdj;
        oracle_check_segment(sim, refdata, segment, reported_call, issues);
    }
}

fn oracle_check_segment(
    sim: &crate::ir::Simulation,
    refdata: &RefDataConfig,
    segment: Segment,
    reported_call: &str,
    issues: &mut Vec<RecordValidationIssue>,
) {
    let Some(region) = sim
        .sequence
        .regions
        .iter()
        .find(|r| r.segment == segment)
    else {
        return;
    };
    let Some(allele_pool) = allele_pool_for_segment(refdata, segment) else {
        return;
    };
    if allele_pool.is_empty() {
        return;
    }
    let assignment = sim.assignments.get(segment);
    let truth_id = assignment.map(|a| a.allele_id);
    let trim_5_cap = assignment.map(|a| a.trim_5 as u32).unwrap_or(0);
    let trim_3_cap = assignment.map(|a| a.trim_3 as u32).unwrap_or(0);

    // Independent rescore via the shared scoring kernel: structural
    // region + NP-region extensions under the assigned allele's trim
    // caps. Mirrors `live_call::walker::call_from_region` so the
    // oracle and the walker agree on the tie-set under arbitrary
    // trim.
    let scores =
        score_alleles_with_extensions(sim, segment, allele_pool, region, trim_5_cap, trim_3_cap);
    let tied_ids = tie_set_ids_at_max_score(&scores);
    if tied_ids.is_empty() {
        return; // No germline evidence; oracle abstains.
    }
    let tied_indices: Vec<usize> = tied_ids.iter().map(|id| id.as_usize()).collect();

    // Expected CSV order: truth allele first when in tie-set,
    // otherwise ascending by allele id (already sorted by the
    // kernel's iteration order).
    let truth_idx = truth_id
        .map(|id| id.as_usize())
        .filter(|&i| i < allele_pool.len());

    let mut expected_order = tied_indices.clone();
    if let Some(t) = truth_idx {
        if let Some(pos) = expected_order.iter().position(|&i| i == t) {
            let truth_first = expected_order.remove(pos);
            expected_order.insert(0, truth_first);
        }
    }

    let expected_names: Vec<String> = expected_order
        .iter()
        .map(|&i| allele_pool[i].name.clone())
        .collect();
    let reported_names: Vec<String> = reported_call
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect();

    // Tie-set equality (order-insensitive).
    let mut reported_sorted = reported_names.clone();
    reported_sorted.sort();
    let mut expected_sorted = expected_names.clone();
    expected_sorted.sort();
    if reported_sorted != expected_sorted {
        issues.push(RecordValidationIssue::AlleleCallTieSetMismatch {
            segment,
            reported: reported_names.clone(),
            recomputed: expected_names.clone(),
        });
        return; // Order check is meaningless if the sets differ.
    }

    // Order check: first element must match expected_order's first.
    if let (Some(reported_first), Some(expected_first)) =
        (reported_names.first(), expected_names.first())
    {
        if reported_first != expected_first {
            let reason = if truth_idx.is_some()
                && expected_order
                    .first()
                    .map(|&i| i == truth_idx.unwrap())
                    .unwrap_or(false)
            {
                AlleleOrderReason::TruthFirstIfInTieSet
            } else {
                AlleleOrderReason::AscendingAlleleIdOtherwise
            };
            issues.push(RecordValidationIssue::AlleleCallOrderMismatch {
                segment,
                reported_first: reported_first.clone(),
                expected_first: expected_first.clone(),
                reason,
            });
        }
    }
}

// Match semantics live in crate::live_call::scoring; this module
// uses them via score_alleles_in_region / tie_set_ids_at_max_score.

// ──────────────────────────────────────────────────────────────────
// C5: Region / live-call structural invariants
//
// Per the §5 architectural audit:
//   - Each Segment (V/D/J) appears at most once in
//     sim.sequence.regions. Live-call and AIRR projection pick the
//     "latest" and "first" respectively; identical-by-construction
//     today, drift-able if invariant breaks.
//   - SegmentLiveCall.hypotheses has length <= 1 in production runs.
//     Projection silently uses hypotheses[0]; multi-hypothesis would
//     be lossy.
// ──────────────────────────────────────────────────────────────────

fn check_region_and_hypothesis_invariants(
    sim: &crate::ir::Simulation,
    issues: &mut Vec<RecordValidationIssue>,
) {
    for seg in [Segment::V, Segment::D, Segment::J] {
        let count = sim
            .sequence
            .regions
            .iter()
            .filter(|r| r.segment == seg)
            .count();
        if count > 1 {
            issues.push(RecordValidationIssue::MultipleRegionsForSegment {
                segment: seg,
                count,
            });
        }
    }

    for seg in [Segment::V, Segment::D, Segment::J] {
        if let Some(call) = sim.segment_calls.get(seg) {
            if call.hypotheses.len() > 1 {
                issues.push(RecordValidationIssue::MultipleHypothesesInLiveCall {
                    segment: seg,
                    count: call.hypotheses.len(),
                });
            }
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    //! Unit tests live alongside the integration tests in
    //! `engine_rs/src/airr_record/tests/validate.rs`. Keeping them
    //! out of this file makes the module easier to skim.
}
