use crate::address::{self, ChoiceAddress, PrimeEnd, VdjSegment};
use crate::assignment::TrimEnd;
use crate::codon::translate_seq;
use crate::ir::{Segment, SimulationEvent};
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
use super::trace_fields::{trace_bool_choice, trace_int_choice};
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

    if rec.sequence.is_empty() {
        // Pre-recombine sim: nothing else to fill in.
        return rec;
    }

    // Trim values from trace. Our DSL records v_3, d_5, d_3, j_5;
    // v_5 and j_3 stay 0.
    rec.v_trim_5 = 0;
    rec.v_trim_3 = trace_int_choice(
        trace,
        ChoiceAddress::Trim {
            segment: VdjSegment::V,
            end: TrimEnd::Three,
        },
    );
    rec.d_trim_5 = trace_int_choice(
        trace,
        ChoiceAddress::Trim {
            segment: VdjSegment::D,
            end: TrimEnd::Five,
        },
    );
    rec.d_trim_3 = trace_int_choice(
        trace,
        ChoiceAddress::Trim {
            segment: VdjSegment::D,
            end: TrimEnd::Three,
        },
    );
    rec.j_trim_5 = trace_int_choice(
        trace,
        ChoiceAddress::Trim {
            segment: VdjSegment::J,
            end: TrimEnd::Five,
        },
    );
    rec.j_trim_3 = 0;

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
    // `Simulation.mutation_count`, stashed by the S5F / Uniform
    // passes at seal time.
    rec.n_mutations = sim.mutation_count as i64;
    rec.mutation_rate = if rec.sequence_length > 0 {
        rec.n_mutations as f64 / rec.sequence_length as f64
    } else {
        0.0
    };
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
        if record.pass_name != address::MUTATE_UNIFORM
            && record.pass_name != address::MUTATE_S5F
        {
            continue;
        }
        for ev in &record.simulation_events {
            let SimulationEvent::BaseChanged {
                segment,
                germline_pos,
                ..
            } = ev
            else {
                continue;
            };
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
                                .find(|s| {
                                    s.start <= pos && pos < s.end
                                })
                                .map(|s| s.label)
                        })
                    });
                    match label {
                        Some(crate::refdata::VSubregionLabel::Fwr1) => {
                            n_fwr1_mutations += 1
                        }
                        Some(crate::refdata::VSubregionLabel::Cdr1) => {
                            n_cdr1_mutations += 1
                        }
                        Some(crate::refdata::VSubregionLabel::Fwr2) => {
                            n_fwr2_mutations += 1
                        }
                        Some(crate::refdata::VSubregionLabel::Cdr2) => {
                            n_cdr2_mutations += 1
                        }
                        Some(crate::refdata::VSubregionLabel::Fwr3) => {
                            n_fwr3_mutations += 1
                        }
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

    // Per-end P-nucleotide length counters (Slice — P-nucleotide
    // v1). Walk the event ledger for `PRegionAdded` events and
    // sum `region.len()` per `end`. Symmetric with the
    // per-segment indel / mutation walks above. The pass-name
    // filter selects the four `p_addition.*` passes by the
    // canonical name constants. Records from cartridges with
    // no P-plane authored carry zeros (no events emitted).
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
    rec.p_v_3_length = p_v_3_length;
    rec.p_d_5_length = p_d_5_length;
    rec.p_d_3_length = p_d_3_length;
    rec.p_j_5_length = p_j_5_length;

    // Observation-stage end-loss / primer-trim amounts. Distinct
    // provenance from recombination-stage `v_trim_5/3` /
    // `j_trim_5/3`; carries the *realized* (post-clamp) byte count
    // the `EndLossPass` recorded. Defaults to 0 when no end-loss
    // pass ran. See `docs/primer_trim_end_loss_audit.md`.
    rec.end_loss_5_length =
        trace_int_choice(trace, ChoiceAddress::CorruptEndLoss(PrimeEnd::Five));
    rec.end_loss_3_length =
        trace_int_choice(trace, ChoiceAddress::CorruptEndLoss(PrimeEnd::Three));
    rec.is_contaminant = trace_bool_choice(trace, ChoiceAddress::CorruptContaminantApplied);
    // `rec.d_inverted` / `rec.receptor_revision_applied` /
    // `rec.original_v_call` are populated earlier (before the
    // pre-recombine early return) — `d_inverted` is sim-level
    // metadata and the receptor-revision fields are trace-sourced;
    // neither needs the pool to be assembled to read its value.

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
        // Trim-only fallback for partial-state outcomes where the
        // walker couldn't compute `ref_ranges[V]` (e.g. no V
        // region present, allele assigned for projection only).
        //
        // **End-loss interaction.** The walker encodes a 5'
        // end-loss as `D` (deletion) ops at the start of V's
        // CIGAR — those positions are "covered by the
        // alignment" in the AIRR sense even though they were
        // removed from the read. The walker therefore reports
        // `v_germline_start = v_trim_5` regardless of end-loss.
        // This fallback matches that semantic intentionally: a
        // shift-by-end-loss formulation would diverge from the
        // walker and break AIRR-tool parity.
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
        // Mirror of the V fallback: walker encodes 3' end-loss as
        // trailing `D` ops, keeping `j_germline_end` at
        // `allele_len - j_trim_3`. Fallback matches.
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
