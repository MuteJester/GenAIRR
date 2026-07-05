//! C2: Counter provenance.

use super::*;

use crate::address::{ChoiceAddress, PrimeEnd};
use crate::ir::{Segment, SimulationEvent};
use crate::trace::ChoiceValue;

pub(super) fn check_counters(
    record: &AirrRecord,
    outcome: &Outcome,
    refdata: &RefDataConfig,
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

    // Per-segment SHM counters from the event ledger of the
    // mutate.{uniform,s5f} passes only. Mirrors the indel walk
    // above but counts `BaseChanged` events instead of indel
    // events, and rolls NP1+NP2 into the single NP bucket. The
    // four per-bucket checks fire `N*MutationsMismatch` on
    // disagreement; the sum-invariant cross-check fires
    // `MutationCountSumMismatch`.
    //
    // Same walk also produces the V-subregion partition (Slice —
    // V-Subregion Mutation Counters): each V `BaseChanged.germline_pos`
    // is matched against the assigned V allele's `subregions` table
    // from scratch, independent of the projection's bucketing. Six
    // per-bucket mismatches plus a cross-field sum invariant
    // (`VSubregionMutationCountSumMismatch`) fire on disagreement.
    let mut shm_by_segment = [0i64; Segment::COUNT];
    let mut shm_by_subregion = [0i64; 5];
    let mut shm_v_unannotated = 0i64;
    let v_subregions_for_validate: Option<&[crate::refdata::VSubregion]> = outcome
        .final_simulation()
        .assignments
        .get(Segment::V)
        .and_then(|inst| refdata.v_pool.get(inst.allele_id))
        .map(|allele| allele.subregions.as_slice());
    for er in outcome.events() {
        if er.pass_name != crate::address::MUTATE_UNIFORM
            && er.pass_name != crate::address::MUTATE_S5F
        {
            continue;
        }
        for ev in &er.simulation_events {
            let SimulationEvent::BaseChanged {
                segment,
                germline_pos,
                ..
            } = ev
            else {
                continue;
            };
            shm_by_segment[*segment as usize] += 1;
            if *segment == Segment::V {
                let label = v_subregions_for_validate.and_then(|subs| {
                    germline_pos.and_then(|pos| {
                        subs.iter()
                            .find(|s| s.start <= pos && pos < s.end)
                            .map(|s| s.label)
                    })
                });
                match label {
                    Some(crate::refdata::VSubregionLabel::Fwr1) => {
                        shm_by_subregion[0] += 1
                    }
                    Some(crate::refdata::VSubregionLabel::Cdr1) => {
                        shm_by_subregion[1] += 1
                    }
                    Some(crate::refdata::VSubregionLabel::Fwr2) => {
                        shm_by_subregion[2] += 1
                    }
                    Some(crate::refdata::VSubregionLabel::Cdr2) => {
                        shm_by_subregion[3] += 1
                    }
                    Some(crate::refdata::VSubregionLabel::Fwr3) => {
                        shm_by_subregion[4] += 1
                    }
                    None => shm_v_unannotated += 1,
                }
            }
        }
    }
    let expected_v = shm_by_segment[Segment::V as usize];
    let expected_d = shm_by_segment[Segment::D as usize];
    let expected_j = shm_by_segment[Segment::J as usize];
    let expected_np = shm_by_segment[Segment::Np1 as usize]
        + shm_by_segment[Segment::Np2 as usize];
    if record.n_v_mutations != expected_v {
        issues.push(RecordValidationIssue::NVMutationsMismatch {
            reported: record.n_v_mutations,
            event_count: expected_v,
        });
    }
    if record.n_d_mutations != expected_d {
        issues.push(RecordValidationIssue::NDMutationsMismatch {
            reported: record.n_d_mutations,
            event_count: expected_d,
        });
    }
    if record.n_j_mutations != expected_j {
        issues.push(RecordValidationIssue::NJMutationsMismatch {
            reported: record.n_j_mutations,
            event_count: expected_j,
        });
    }
    if record.n_np_mutations != expected_np {
        issues.push(RecordValidationIssue::NNpMutationsMismatch {
            reported: record.n_np_mutations,
            event_count: expected_np,
        });
    }
    // Sum-invariant: the four per-bucket fields must add up to
    // ``n_mutations``. Validates the consistency of any consumer-
    // supplied record dict; the engine-projected record satisfies
    // it by construction.
    let sum_of_buckets = record
        .n_v_mutations
        .saturating_add(record.n_d_mutations)
        .saturating_add(record.n_j_mutations)
        .saturating_add(record.n_np_mutations);
    if record.n_mutations != sum_of_buckets {
        issues.push(RecordValidationIssue::MutationCountSumMismatch {
            reported_total: record.n_mutations,
            sum_of_buckets,
        });
    }
    // V-subregion partition mismatch checks. Each per-bucket
    // mismatch fires `N<Region>MutationsMismatch`; the sum
    // invariant fires `VSubregionMutationCountSumMismatch`.
    if record.n_fwr1_mutations != shm_by_subregion[0] {
        issues.push(RecordValidationIssue::NFwr1MutationsMismatch {
            reported: record.n_fwr1_mutations,
            event_count: shm_by_subregion[0],
        });
    }
    if record.n_cdr1_mutations != shm_by_subregion[1] {
        issues.push(RecordValidationIssue::NCdr1MutationsMismatch {
            reported: record.n_cdr1_mutations,
            event_count: shm_by_subregion[1],
        });
    }
    if record.n_fwr2_mutations != shm_by_subregion[2] {
        issues.push(RecordValidationIssue::NFwr2MutationsMismatch {
            reported: record.n_fwr2_mutations,
            event_count: shm_by_subregion[2],
        });
    }
    if record.n_cdr2_mutations != shm_by_subregion[3] {
        issues.push(RecordValidationIssue::NCdr2MutationsMismatch {
            reported: record.n_cdr2_mutations,
            event_count: shm_by_subregion[3],
        });
    }
    if record.n_fwr3_mutations != shm_by_subregion[4] {
        issues.push(RecordValidationIssue::NFwr3MutationsMismatch {
            reported: record.n_fwr3_mutations,
            event_count: shm_by_subregion[4],
        });
    }
    if record.n_v_unannotated_mutations != shm_v_unannotated {
        issues.push(RecordValidationIssue::NVUnannotatedMutationsMismatch {
            reported: record.n_v_unannotated_mutations,
            event_count: shm_v_unannotated,
        });
    }
    let sum_of_subregion_buckets = record
        .n_fwr1_mutations
        .saturating_add(record.n_cdr1_mutations)
        .saturating_add(record.n_fwr2_mutations)
        .saturating_add(record.n_cdr2_mutations)
        .saturating_add(record.n_fwr3_mutations)
        .saturating_add(record.n_v_unannotated_mutations);
    if record.n_v_mutations != sum_of_subregion_buckets {
        issues.push(
            RecordValidationIssue::VSubregionMutationCountSumMismatch {
                reported_v_total: record.n_v_mutations,
                sum_of_subregion_buckets,
            },
        );
    }

    // Per-end P-nucleotide length counters (Slice —
    // P-nucleotide v1). Independent event-ledger recompute:
    // walk `PRegionAdded { end, region }` events from the
    // matching `p_addition.*` passes and sum `region.len()`
    // per end. Catches downstream consumers that tamper with
    // the four `p_*_length` fields (record edits, fork-
    // patched builders, deserialised dicts).
    let mut p_v_3_recompute = 0i64;
    let mut p_d_5_recompute = 0i64;
    let mut p_d_3_recompute = 0i64;
    let mut p_j_5_recompute = 0i64;
    for ev_record in outcome.events() {
        let is_p_addition = ev_record.pass_name == crate::address::P_ADDITION_V_3
            || ev_record.pass_name == crate::address::P_ADDITION_D_5
            || ev_record.pass_name == crate::address::P_ADDITION_D_3
            || ev_record.pass_name == crate::address::P_ADDITION_J_5;
        if !is_p_addition {
            continue;
        }
        for ev in &ev_record.simulation_events {
            if let crate::ir::SimulationEvent::PRegionAdded { end, region } = ev {
                let len = region.len() as i64;
                match end {
                    crate::address::PEnd::V3 => p_v_3_recompute += len,
                    crate::address::PEnd::D5 => p_d_5_recompute += len,
                    crate::address::PEnd::D3 => p_d_3_recompute += len,
                    crate::address::PEnd::J5 => p_j_5_recompute += len,
                }
            }
        }
    }
    if record.p_v_3_length != p_v_3_recompute {
        issues.push(RecordValidationIssue::PLengthMismatch {
            end: crate::address::PEnd::V3,
            reported: record.p_v_3_length,
            event_count: p_v_3_recompute,
        });
    }
    if record.p_d_5_length != p_d_5_recompute {
        issues.push(RecordValidationIssue::PLengthMismatch {
            end: crate::address::PEnd::D5,
            reported: record.p_d_5_length,
            event_count: p_d_5_recompute,
        });
    }
    if record.p_d_3_length != p_d_3_recompute {
        issues.push(RecordValidationIssue::PLengthMismatch {
            end: crate::address::PEnd::D3,
            reported: record.p_d_3_length,
            event_count: p_d_3_recompute,
        });
    }
    if record.p_j_5_length != p_j_5_recompute {
        issues.push(RecordValidationIssue::PLengthMismatch {
            end: crate::address::PEnd::J5,
            reported: record.p_j_5_length,
            event_count: p_j_5_recompute,
        });
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

    // D inversion provenance (Slice E). Expected value reads from
    // the simulation's final D assignment; defaults to `false` when
    // D is absent (VJ chains) — matching the builder's `unwrap_or`.
    let expected_inverted = sim
        .assignments
        .get(Segment::D)
        .map(|inst| inst.orientation.is_reverse())
        .unwrap_or(false);
    if record.d_inverted != expected_inverted {
        issues.push(RecordValidationIssue::DInvertedMismatch {
            reported: record.d_inverted,
            expected: expected_inverted,
        });
    }

    // Receptor revision provenance — IR-sourced (Bug D fix).
    // Originally trace-sourced; the descendant trace omits pre-fork
    // choices, which made every clonal descendant's projection
    // disagree with the parent's actual revision state. The
    // assignments slot persists across the parent→descendant
    // boundary, so reading from it produces identical behaviour
    // for non-clonal and clonal pipelines.
    let v_inst = sim.assignments.get(Segment::V);
    let expected_applied = v_inst
        .map(|inst| inst.receptor_revision_original_id.is_some())
        .unwrap_or(false);
    if record.receptor_revision_applied != expected_applied {
        issues.push(RecordValidationIssue::ReceptorRevisionAppliedMismatch {
            reported: record.receptor_revision_applied,
            expected: expected_applied,
        });
    }

    let expected_original_v_call = if expected_applied {
        original_v_name_from_assignment(
            v_inst.expect("expected_applied implies v_inst is Some"),
            refdata,
        )
    } else {
        String::new()
    };
    if record.original_v_call != expected_original_v_call {
        issues.push(RecordValidationIssue::OriginalVCallMismatch {
            reported: record.original_v_call.clone(),
            expected: expected_original_v_call,
        });
    }
}

/// IR-sourced counterpart of the (now-removed) trace-based helper.
/// Resolves the pre-revision V allele's refdata name from the
/// persistent provenance slot the receptor-revision pass installs.
/// Mirrors `original_v_call_from_assignment` in
/// `airr_record::builder` so the validator stays self-contained.
fn original_v_name_from_assignment(
    v_inst: &crate::assignment::AlleleInstance,
    refdata: &RefDataConfig,
) -> String {
    let Some(original_id) = v_inst.receptor_revision_original_id else {
        return String::new();
    };
    refdata
        .get(Segment::V, original_id)
        .map(|a| a.name.clone())
        .unwrap_or_default()
}
