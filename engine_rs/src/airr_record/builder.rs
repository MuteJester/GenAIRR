use crate::address::{self, ChoiceAddress, PrimeEnd};
use crate::codon::translate_seq;
use crate::ir::{Segment, SimulationEvent};
use crate::pass::Outcome;
use crate::refdata::RefDataConfig;

use super::junction::aa_slice_for_region;
use super::locus::{derive_locus, locus_from_refdata};
use super::projection::{projected_call_name, recalled_call_name};
use super::sequence::{apply_rev_comp_projection, bytes_to_string, pool_bases};
use super::trace_fields::{trace_bool_choice, trace_int_choice};
use super::walk::walk_alignment_columns;
use super::AirrRecord;

const DUPLEX_UNSUPPORTED: &str =
    "legacy AIRR records cannot represent unknown duplex P/NP metadata; use the duplex projection API";

/// Checked legacy-record boundary. Duplex descendants require the separate
/// nullable duplex projector, even when their sequence layout is representable.
pub fn try_build_airr_record(
    outcome: &Outcome,
    refdata: &RefDataConfig,
    sequence_id: &str,
) -> Result<AirrRecord, String> {
    if outcome
        .revisions
        .iter()
        .any(|simulation| simulation.duplex_provenance.is_some())
    {
        return Err(DUPLEX_UNSUPPORTED.into());
    }
    Ok(build_airr_record(outcome, refdata, sequence_id))
}

/// Build an AIRR record from an outcome + refdata.
///
/// `sequence_id` is set on the record verbatim; pass an empty
/// string for unset.
///
/// # Panics
///
/// The legacy schema cannot represent duplex uncertainty. Call
/// [`try_build_airr_record`] at fallible boundaries, or use the duplex
/// projector for duplex outcomes. This legacy entrypoint fails closed rather
/// than substituting zero for an unknown count.
pub fn build_airr_record(
    outcome: &Outcome,
    refdata: &RefDataConfig,
    sequence_id: &str,
) -> AirrRecord {
    assert!(
        outcome
            .revisions
            .iter()
            .all(|simulation| simulation.duplex_provenance.is_none()),
        "{DUPLEX_UNSUPPORTED}"
    );
    let sim = outcome.final_simulation();
    let trace = &outcome.trace;

    let mut rec = AirrRecord::default();
    let assessment = crate::productivity::assess(sim, refdata);
    rec.productivity_policy = crate::productivity::POLICY.into();
    rec.productivity_scope = crate::productivity::SCOPE.into();
    rec.productivity_reasons = assessment.reasons.iter().map(|r| (*r).into()).collect();
    rec.productive = assessment.productive;
    rec.vj_in_frame = assessment.vj_in_frame;
    rec.stop_codon = assessment.stop_codon;
    rec.coding_end_processing = super::coding_end::project(sim, refdata);
    rec.sequence_id = sequence_id.to_string();
    rec.rev_comp = false;
    rec.c_call = String::new();

    // Bases as a string (preserving case from IR).
    rec.sequence = bytes_to_string(&pool_bases(sim));
    rec.sequence_length = rec.sequence.len() as i64;

    // D inversion (Slice E). Sourced from the final IR's
    // `AlleleAssignments.D.orientation` — independent of the pool's
    // post-assembly bytes. Set BEFORE the pre-recombine early return
    // so a sim that has D assigned but hasn't run assembly yet
    // still surfaces a correct `d_inverted` value (instrumentation
    // and unit-test fixtures both rely on this).
    rec.d_inverted = sim
        .assignments
        .get(Segment::D)
        .map(|inst| inst.orientation.is_reverse())
        .unwrap_or(false);

    // Receptor revision provenance — IR-sourced (Bug D fix).
    // Originally this read from the descendant's trace, but the
    // descendant trace omits pre-fork choices, which silently
    // dropped revision provenance on clonal descendants whose
    // receptor revision had run in the parent. The assignments
    // slot is part of the persistent IR and survives the
    // parent→descendant boundary, so reading from it gives
    // identical results for non-clonal and clonal pipelines.
    // Same pattern as `d_inverted` above.
    let v_inst = sim.assignments.get(Segment::V);
    rec.receptor_revision_applied = v_inst
        .map(|inst| inst.receptor_revision_original_id.is_some())
        .unwrap_or(false);
    rec.original_v_call = if rec.receptor_revision_applied {
        original_v_call_from_assignment(
            v_inst.expect("receptor_revision_applied implies v_inst is Some"),
            refdata,
        )
    } else {
        String::new()
    };

    // Assignments are the current biological state: they survive a clonal
    // fork and receptor revision updates their trims. A descendant trace
    // can omit the founding recombination, and an original trim choice
    // does not describe a replacement allele.
    if let Some(inst) = sim.assignments.get(Segment::V) {
        rec.v_trim_5 = inst.trim_5 as i64;
        rec.v_trim_3 = inst.trim_3 as i64;
    }
    if let Some(inst) = sim.assignments.get(Segment::D) {
        rec.d_trim_5 = inst.trim_5 as i64;
        rec.d_trim_3 = inst.trim_3 as i64;
    }
    if let Some(inst) = sim.assignments.get(Segment::J) {
        rec.j_trim_5 = inst.trim_5 as i64;
        rec.j_trim_3 = inst.trim_3 as i64;
    }

    // Calls + locus.
    let v_id = sim
        .assignments
        .get(Segment::V)
        .copied()
        .map(|i| i.allele_id);
    let d_id = sim
        .assignments
        .get(Segment::D)
        .copied()
        .map(|i| i.allele_id);
    let j_id = sim
        .assignments
        .get(Segment::J)
        .copied()
        .map(|i| i.allele_id);

    rec.v_call = projected_call_name(refdata, sim, Segment::V, v_id);
    rec.d_call = projected_call_name(refdata, sim, Segment::D, d_id);
    rec.j_call = projected_call_name(refdata, sim, Segment::J, j_id);
    rec.recalled_v_call = recalled_call_name(refdata, sim, Segment::V, v_id);
    rec.recalled_d_call = recalled_call_name(refdata, sim, Segment::D, d_id);
    rec.recalled_j_call = recalled_call_name(refdata, sim, Segment::J, j_id);
    rec.locus = derive_locus(&rec.v_call, &rec.j_call, &rec.d_call);
    if rec.locus.is_empty() {
        // Synthetic references need not use locus-prefixed allele names.
        // Infer a locus from the cartridge when assignments cannot name one.
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
    // Initialize query spans from structural regions, then use the alignment
    // walker's actual observed spans below. Re-calls do not move these bounds.
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

    // Mutation + corruption counters. `n_mutations` reads from
    // `Simulation.mutation_count`, stashed by the S5F / Uniform
    // passes at seal time.
    rec.n_mutations = sim.mutation_count as i64;
    // Each mutation stage records its own input molecule length. Later read
    // loss must not inflate this cumulative event burden.
    rec.mutation_rate = sim.mutation_burden;
    rec.n_pcr_errors = trace_int_choice(trace, ChoiceAddress::CorruptPcrCount);
    rec.n_quality_errors = trace_int_choice(trace, ChoiceAddress::CorruptQualityCount);
    // Count realized structural indels from the event ledger, not the
    // trace's `corrupt.indel.count` (which is the *attempted* count
    // and includes pool-empty / contract no-op sentinels). Scoped to
    // the IndelPass: `EndLossPass` reuses `IndelDeleted` as its IR
    // primitive but isn't a structural indel. See
    // `docs/indel_provenance_audit.md` §6.1 / §6.2.
    let mut n_indels = 0i64;
    let mut n_v_indels = 0i64;
    let mut n_d_indels = 0i64;
    let mut n_j_indels = 0i64;
    for record in outcome.events() {
        if record.pass_name != address::CORRUPT_INDEL {
            continue;
        }
        for ev in &record.simulation_events {
            let segment = match ev {
                SimulationEvent::IndelInserted { segment, .. } => *segment,
                SimulationEvent::IndelDeleted { segment, .. } => *segment,
                _ => continue,
            };
            n_indels += 1;
            match segment {
                Segment::V => n_v_indels += 1,
                Segment::D => n_d_indels += 1,
                Segment::J => n_j_indels += 1,
                // NP1 / NP2 events count toward `n_indels` only.
                Segment::Np1 | Segment::Np2 => {}
            }
        }
    }
    rec.n_indels = n_indels;
    rec.n_v_indels = n_v_indels;
    rec.n_d_indels = n_d_indels;
    rec.n_j_indels = n_j_indels;

    // Per-segment biological-SHM mutation counters. Same shape as
    // the indel walk above: scan the per-pass event ledger, filter
    // to the two SHM passes by their canonical pass-name constants,
    // bucket every `BaseChanged` event by its carried segment. NP1
    // + NP2 roll into a single `n_np_mutations` bucket — same
    // collapse the `segment_rates` DSL uses (per the audit's
    // `n_v/d/j_mutations` + `n_np_mutations` recommendation).
    // Corruption-stage `BaseChanged` events (PCR, quality, contaminant)
    // are excluded by the pass-name filter. By construction the
    // four fields partition `n_mutations`; the validator surfaces
    // any divergence via `MutationCountSumMismatch`.
    let mut n_v_mutations = 0i64;
    let mut n_d_mutations = 0i64;
    let mut n_j_mutations = 0i64;
    let mut n_np_mutations = 0i64;
    // V-subregion partition of `n_v_mutations` (Slice — V-Subregion
    // Mutation Counters). Aggregated in the same pass as the
    // per-segment counters: each V `BaseChanged` event's
    // `germline_pos: Option<u16>` is matched against the assigned V
    // allele's `subregions` table; events landing in one of the five
    // canonical IMGT intervals route to the matching bucket,
    // everything else (missing assignment, empty annotations,
    // `germline_pos = None`, V-side CDR3 stretch outside FWR1..FWR3)
    // routes to `n_v_unannotated_mutations`. See
    // `docs/v_subregion_mutation_counters_audit.md`.
    let mut n_fwr1_mutations = 0i64;
    let mut n_cdr1_mutations = 0i64;
    let mut n_fwr2_mutations = 0i64;
    let mut n_cdr2_mutations = 0i64;
    let mut n_fwr3_mutations = 0i64;
    let mut n_v_unannotated_mutations = 0i64;
    // Hoist the assigned V allele's subregion table out of the
    // event loop — it's invariant within a record. `None` is the
    // pre-recombine state (no V assignment); under that condition
    // every V `BaseChanged` (rare but possible in pathological
    // plans) routes to the unannotated bucket.
    let v_subregions: Option<&[crate::refdata::VSubregion]> = sim
        .assignments
        .get(Segment::V)
        .and_then(|inst| refdata.v_pool.get(inst.allele_id))
        .map(|allele| allele.subregions.as_slice());
    for record in outcome.events() {
        if record.pass_name != address::MUTATE_UNIFORM && record.pass_name != address::MUTATE_S5F {
            continue;
        }
        for ev in &record.simulation_events {
            let SimulationEvent::BaseChanged {
                segment,
                germline_pos,
                old_base,
                new_base,
                ..
            } = ev
            else {
                continue;
            };
            // A write whose drawn base equalled the base already there
            // changed nothing, so it is not a mutation and must not reach
            // these counters either. Skipping it here is what keeps the
            // four fields a partition of `n_mutations`: the transaction
            // subtracts the same writes from the total, and omitting this
            // half broke the sum invariant the comment above asserts --
            // caught by test_per_segment_counts_sum_to_the_total.
            //
            // Case-insensitive for the same reason as the transaction's
            // check: a templated base is stored lowercase and an SHM write
            // emits uppercase, so a byte comparison would see a change
            // where a consumer comparing to the germline sees none.
            if old_base.eq_ignore_ascii_case(new_base) {
                continue;
            }
            match segment {
                Segment::V => {
                    n_v_mutations += 1;
                    // Dispatch into the V-subregion partition.
                    // Three "unannotated" cases collapse into the
                    // same arm:
                    //   1. No V allele assigned (v_subregions is None).
                    //   2. germline_pos is None (indel-inserted base).
                    //   3. The position falls outside every subregion
                    //      interval (V-side CDR3 stretch, or an
                    //      allele with empty `subregions`).
                    let label = v_subregions.and_then(|subs| {
                        germline_pos.and_then(|pos| {
                            subs.iter()
                                .find(|s| s.start <= pos && pos < s.end)
                                .map(|s| s.label)
                        })
                    });
                    match label {
                        Some(crate::refdata::VSubregionLabel::Fwr1) => n_fwr1_mutations += 1,
                        Some(crate::refdata::VSubregionLabel::Cdr1) => n_cdr1_mutations += 1,
                        Some(crate::refdata::VSubregionLabel::Fwr2) => n_fwr2_mutations += 1,
                        Some(crate::refdata::VSubregionLabel::Cdr2) => n_cdr2_mutations += 1,
                        Some(crate::refdata::VSubregionLabel::Fwr3) => n_fwr3_mutations += 1,
                        None => n_v_unannotated_mutations += 1,
                    }
                }
                Segment::D => n_d_mutations += 1,
                Segment::J => n_j_mutations += 1,
                Segment::Np1 | Segment::Np2 => n_np_mutations += 1,
            }
        }
    }
    rec.n_v_mutations = n_v_mutations;
    rec.n_d_mutations = n_d_mutations;
    rec.n_j_mutations = n_j_mutations;
    rec.n_np_mutations = n_np_mutations;
    rec.n_fwr1_mutations = n_fwr1_mutations;
    rec.n_cdr1_mutations = n_cdr1_mutations;
    rec.n_fwr2_mutations = n_fwr2_mutations;
    rec.n_cdr2_mutations = n_cdr2_mutations;
    rec.n_fwr3_mutations = n_fwr3_mutations;
    rec.n_v_unannotated_mutations = n_v_unannotated_mutations;

    // Persistent P-addition provenance survives clonal forks even when the
    // descendant ledger contains only observation passes. Retain the event
    // fallback for synthetic outcomes with an event but no stamped assignment.
    let mut p_v_3_length = 0i64;
    let mut p_d_5_length = 0i64;
    let mut p_d_3_length = 0i64;
    let mut p_j_5_length = 0i64;
    for record in outcome.events() {
        let is_p_addition = record.pass_name == crate::address::P_ADDITION_V_3
            || record.pass_name == crate::address::P_ADDITION_D_5
            || record.pass_name == crate::address::P_ADDITION_D_3
            || record.pass_name == crate::address::P_ADDITION_J_5;
        if !is_p_addition {
            continue;
        }
        for ev in &record.simulation_events {
            if let SimulationEvent::PRegionAdded { end, region } = ev {
                let len = region.len() as i64;
                match end {
                    crate::address::PEnd::V3 => p_v_3_length += len,
                    crate::address::PEnd::D5 => p_d_5_length += len,
                    crate::address::PEnd::D3 => p_d_3_length += len,
                    crate::address::PEnd::J5 => p_j_5_length += len,
                }
            }
        }
    }
    if let Some(inst) = sim.assignments.get(Segment::V) {
        p_v_3_length = p_v_3_length.max(inst.p_length_3 as i64);
    }
    if let Some(inst) = sim.assignments.get(Segment::D) {
        p_d_5_length = p_d_5_length.max(inst.p_length_5 as i64);
        p_d_3_length = p_d_3_length.max(inst.p_length_3 as i64);
    }
    if let Some(inst) = sim.assignments.get(Segment::J) {
        p_j_5_length = p_j_5_length.max(inst.p_length_5 as i64);
    }
    rec.p_v_3_length = p_v_3_length;
    rec.p_d_5_length = p_d_5_length;
    rec.p_d_3_length = p_d_3_length;
    rec.p_j_5_length = p_j_5_length;

    // Observation-stage end-loss / primer-trim amounts. Distinct
    // provenance from recombination-stage `v_trim_5/3` /
    // `j_trim_5/3`; carries the *realized* (post-clamp) byte count
    // the `EndLossPass` recorded. Defaults to 0 when no end-loss
    // pass ran. See `docs/primer_trim_end_loss_audit.md`.
    rec.end_loss_5_length = trace_int_choice(trace, ChoiceAddress::CorruptEndLoss(PrimeEnd::Five));
    rec.end_loss_3_length = trace_int_choice(trace, ChoiceAddress::CorruptEndLoss(PrimeEnd::Three));
    rec.is_contaminant = trace_bool_choice(trace, ChoiceAddress::CorruptContaminantApplied);
    // Empty observations can result from complete end loss after
    // recombination. Keep their realized counters and assignments.
    if rec.sequence.is_empty() {
        rec.rev_comp = trace_bool_choice(trace, ChoiceAddress::CorruptRevCompApplied);
        return rec;
    }
    // `rec.d_inverted` / `rec.receptor_revision_applied` /
    // `rec.original_v_call` are populated earlier (before the
    // pre-recombine early return) — `d_inverted` is sim-level
    // metadata and the receptor-revision fields are trace-sourced;
    // neither needs the pool to be assembled to read its value.

    // Translate using trusted coding-phase evidence; landmark uncertainty
    // must not invent a frame or a canonical junction window.
    if let Some(frame) = assessment.frame {
        rec.sequence_aa = translate_seq(&rec.sequence[frame.phase.min(rec.sequence.len())..]);
    }
    if let Some((start, end)) = assessment.junction {
        rec.junction = rec.sequence[start..end].to_string();
        rec.junction_aa = translate_seq(&rec.junction);
        rec.junction_start = Some(start as i64);
        rec.junction_end = Some(end as i64);
        rec.junction_length = Some((end - start) as i64);
    }

    // Single column walk: renders sequence_alignment, germline_alignment,
    // and germline_alignment_d_mask. Per-segment CIGARs, identities,
    // and ranges come from the segment-projection kernel that the
    // walker internally pre-calls.
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

    // Refine query coordinates from the generating-reference walk.
    // Each segment's M+I span equals its observed query span.
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

    // AIRR NP fields describe the combined N/P interval between the final
    // segment alignments. P bases have no structural NP region, so deriving
    // these strings from NP regions alone silently drops their contribution.
    // Generating boundaries preserve sampled NP bases even when they happen
    // to match adjacent trimmed germline sequence.
    let np1_bounds = np_interval(
        rec.v_sequence_end,
        if refdata.chain_type.has_d() {
            rec.d_sequence_start
        } else {
            rec.j_sequence_start
        },
        rec.sequence.len(),
    );
    let np2_bounds = if refdata.chain_type.has_d() {
        np_interval(rec.d_sequence_end, rec.j_sequence_start, rec.sequence.len())
    } else {
        None
    };
    if let Some((start, end)) = np1_bounds {
        rec.np1 = rec.sequence[start..end].to_string();
        rec.np1_length = (end - start) as i64;
    }
    if let Some((start, end)) = np2_bounds {
        rec.np2 = rec.sequence[start..end].to_string();
        rec.np2_length = (end - start) as i64;
    }
    if let Some(frame) = crate::coding::coding_frame(sim, refdata) {
        if let Some((start, end)) = np1_bounds {
            rec.np1_aa = aa_slice_for_region(start, end, frame.phase, &rec.sequence_aa);
        }
        if let Some((start, end)) = np2_bounds {
            rec.np2_aa = aa_slice_for_region(start, end, frame.phase, &rec.sequence_aa);
        }
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
    // Partial-state outcomes can have an assignment but no alignment
    // columns. Their fallback must obey the same observed coverage bounds;
    // a completely clipped segment has no germline alignment interval.
    if let Some((g_start, g_end)) = walk.ref_ranges[0] {
        rec.v_germline_start = Some(g_start);
        rec.v_germline_end = Some(g_end);
    } else if let Some(reference) = super::projection::projected_reference(sim, refdata, Segment::V)
    {
        let (start, end) = reference.observed_bounds(rec.v_trim_5, rec.v_trim_3);
        if start < end {
            rec.v_germline_start = Some(start as i64);
            rec.v_germline_end = Some(end as i64);
        }
    }
    if let Some((g_start, g_end)) = walk.ref_ranges[1] {
        rec.d_germline_start = Some(g_start);
        rec.d_germline_end = Some(g_end);
    } else if let Some(reference) = super::projection::projected_reference(sim, refdata, Segment::D)
    {
        let (start, end) = reference.observed_bounds(rec.d_trim_5, rec.d_trim_3);
        if start < end {
            rec.d_germline_start = Some(start as i64);
            rec.d_germline_end = Some(end as i64);
        }
    }
    if let Some((g_start, g_end)) = walk.ref_ranges[2] {
        rec.j_germline_start = Some(g_start);
        rec.j_germline_end = Some(g_end);
    } else if let Some(reference) = super::projection::projected_reference(sim, refdata, Segment::J)
    {
        let (start, end) = reference.observed_bounds(rec.j_trim_5, rec.j_trim_3);
        if start < end {
            rec.j_germline_start = Some(start as i64);
            rec.j_germline_end = Some(end as i64);
        }
    }

    rec.v_identity = walk.identities[0];
    rec.d_identity = walk.identities[1];
    rec.j_identity = walk.identities[2];

    // RevCompPass records the observation strand without touching the IR.
    // Only the emitted sequence flips; biological annotations continue to
    // index its coding-oriented reverse complement, as AIRR requires.
    if trace_bool_choice(trace, ChoiceAddress::CorruptRevCompApplied) {
        apply_rev_comp_projection(&mut rec);
    }

    // Paired-end / read-layout projection (Slice C). The
    // `PairedEndSamplingPass` writes three Int records on the
    // trace; if all three are present, the kernel populates the
    // eight AIRR fields. Absent records leave the fields at their
    // Slice A defaults (`""` / `None` / `0`) — the canonical
    // no-paired-end baseline path.
    //
    // Runs AFTER `apply_rev_comp_projection` per the audit §7
    // composition order: paired-end windows are drawn against the
    // finalised, possibly rev-comped, molecule.
    apply_paired_end_projection_from_trace(&mut rec, trace);

    rec
}

/// A non-overlapping inter-segment interval in the emitted sequence.
/// Missing alignments do not establish an NP interval, and overlapping
/// segment hypotheses have no intervening nucleotides.
fn np_interval(
    left_end: Option<i64>,
    right_start: Option<i64>,
    sequence_len: usize,
) -> Option<(usize, usize)> {
    let start = usize::try_from(left_end?).ok()?;
    let end = usize::try_from(right_start?).ok()?;
    (start <= end && end <= sequence_len).then_some((start, end))
}

/// Read the three paired-end trace records and apply the
/// projection kernel when all three are present.
///
/// - All three absent → no-paired-end baseline; `rec` stays at
///   Slice A defaults.
/// - All three present + kernel succeeds → populate the eight
///   fields.
/// - All three present + kernel rejects (e.g. `insert_size >
///   sequence_length` from a refdata swap) → set
///   `read_layout = "paired_end"` and leave coords / strings
///   empty so the Slice B validator surfaces the divergence as
///   `ReadWindowOutOfBounds` rather than silently producing
///   defaults that look like the no-layout baseline.
/// - Partial (one or two of three) → treated as not-opted-in;
///   the `PairedEndSamplingPass` validates atomicity before
///   recording, so this case can only arise from a hand-edited
///   trace. Leaving the fields at their Slice A defaults makes
///   the no-layout validator catch the partial trace through the
///   `PairedEndFieldWithoutLayout` channel.
fn apply_paired_end_projection_from_trace(
    rec: &mut super::AirrRecord,
    trace: &crate::trace::Trace,
) {
    use crate::airr_record::sequence::project_paired_end_layout;
    use crate::airr_record::trace_fields::trace_int_choice_opt;

    let r1_length = trace_int_choice_opt(trace, ChoiceAddress::PairedEndR1Length);
    let r2_length = trace_int_choice_opt(trace, ChoiceAddress::PairedEndR2Length);
    let insert_size = trace_int_choice_opt(trace, ChoiceAddress::PairedEndInsertSize);
    let (Some(r1), Some(r2), Some(insert)) = (r1_length, r2_length, insert_size) else {
        return;
    };

    match project_paired_end_layout(rec, r1, r2, insert) {
        Ok(p) => {
            rec.read_layout = p.read_layout;
            rec.r1_sequence = p.r1_sequence;
            rec.r2_sequence = p.r2_sequence;
            rec.r1_start = p.r1_start;
            rec.r1_end = p.r1_end;
            rec.r2_start = p.r2_start;
            rec.r2_end = p.r2_end;
            rec.insert_size = p.insert_size;
        }
        Err(_) => {
            // The kernel rejected the geometry against the
            // finalised `rec.sequence` (e.g. `insert_size >
            // sequence_length`). Surface the failure to the
            // validator by tagging the layout without populating
            // coords/strings — `ReadWindowOutOfBounds` fires for
            // the missing coord pair.
            rec.read_layout = "paired_end".to_string();
        }
    }
}

/// Resolve the pre-revision V allele's refdata name from the IR
/// assignments. The caller has already established that receptor
/// revision applied (i.e. `inst.receptor_revision_original_id` is
/// `Some`); a missing or unresolvable id at that point is a
/// structural bug.
///
/// Returns an empty string when the recorded id doesn't resolve in
/// refdata — the same documented fallback the old
/// trace-sourced helper used for refdata swaps between record and
/// replay. The validator's `OriginalVCallMismatch` continues to
/// catch divergence between the record value and the expected.
fn original_v_call_from_assignment(
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
