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
    /// Per-segment SHM counter divergence. Re-derived by walking
    /// `outcome.events()`, filtering to the SHM passes
    /// (`mutate.uniform` / `mutate.s5f`), and counting
    /// `SimulationEvent::BaseChanged` by carried segment. NP1+NP2
    /// roll into the NP bucket. See
    /// `docs/mutation_provenance_audit.md`.
    NVMutationsMismatch { reported: i64, event_count: i64 },
    NDMutationsMismatch { reported: i64, event_count: i64 },
    NJMutationsMismatch { reported: i64, event_count: i64 },
    NNpMutationsMismatch { reported: i64, event_count: i64 },
    /// Sum-invariant cross-check: the four per-segment SHM
    /// counters must add up to `n_mutations`. Fires when a
    /// downstream record's reported fields drift from this
    /// arithmetic identity (e.g. someone bumped one bucket
    /// without bumping the global).
    MutationCountSumMismatch {
        reported_total: i64,
        sum_of_buckets: i64,
    },
    /// Per-V-subregion biological-SHM mutation counter mismatches.
    /// Each variant compares the AIRR record's reported field
    /// against an independent recompute over `outcome.events()`
    /// — the validator re-walks the same event ledger the
    /// builder did, applying the same pass-name filter
    /// (`mutate.uniform` + `mutate.s5f`), and matches each V
    /// `BaseChanged.germline_pos` against the assigned V allele's
    /// `subregions` table from scratch. Surfaces when a downstream
    /// record's reported per-subregion fields drift from the
    /// event-derived ground truth. Mirrors the per-segment
    /// `N{V,D,J,Np}MutationsMismatch` shape.
    NFwr1MutationsMismatch { reported: i64, event_count: i64 },
    NCdr1MutationsMismatch { reported: i64, event_count: i64 },
    NFwr2MutationsMismatch { reported: i64, event_count: i64 },
    NCdr2MutationsMismatch { reported: i64, event_count: i64 },
    NFwr3MutationsMismatch { reported: i64, event_count: i64 },
    NVUnannotatedMutationsMismatch { reported: i64, event_count: i64 },
    /// Sum-invariant cross-check for the V-subregion partition:
    /// the five per-subregion buckets plus the unannotated
    /// bucket must add up to `n_v_mutations`. Fires when a
    /// downstream record's reported fields drift from the
    /// arithmetic identity. See
    /// `docs/v_subregion_mutation_counters_audit.md`.
    VSubregionMutationCountSumMismatch {
        reported_v_total: i64,
        sum_of_subregion_buckets: i64,
    },
    /// Per-end P-nucleotide length (`p_v_3_length` /
    /// `p_d_5_length` / `p_d_3_length` / `p_j_5_length`)
    /// disagrees with the recompute from the event ledger.
    /// `end` discriminates which side fired; the validator
    /// recomputes by summing `region.len()` over
    /// `PRegionAdded { end, region }` events emitted by the
    /// matching `p_addition.*` pass. Slice — P-nucleotide v1.
    PLengthMismatch {
        end: crate::address::PEnd,
        reported: i64,
        event_count: i64,
    },
    EndLossLengthMismatch {
        side: PrimeEnd,
        reported: i64,
        trace_count: i64,
    },
    /// `record.d_inverted` disagrees with the simulation's final D
    /// orientation. `expected` is read from
    /// `Simulation.assignments.get(Segment::D).orientation.is_reverse()`
    /// (defaulting to `false` when D is absent). Surfaces when a
    /// downstream consumer manually edits the AIRR record or when
    /// a fork of the AIRR builder forgets to populate the field.
    DInvertedMismatch {
        reported: bool,
        expected: bool,
    },
    /// `record.receptor_revision_applied` disagrees with the trace.
    /// `expected` is the Bool at `receptor_revision.applied`,
    /// defaulting to `false` when the address is absent (no revision
    /// step ran).
    ReceptorRevisionAppliedMismatch {
        reported: bool,
        expected: bool,
    },
    /// `record.original_v_call` disagrees with the trace+refdata.
    /// `expected` is the V allele name from the trace's first
    /// `sample_allele.v` record when
    /// `receptor_revision_applied=true`, and the empty string when
    /// the revision didn't fire. Surfaces when a fork of the AIRR
    /// builder forgets to populate the field, or when refdata
    /// drift between record-time and replay leaves the recorded
    /// allele id unresolvable.
    OriginalVCallMismatch {
        reported: String,
        expected: String,
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

    // ── C6: Paired-end / read-layout invariants ────────────────
    //
    // Slice A of the paired-end roadmap. Five variants land now;
    // only `PairedEndFieldWithoutLayout` fires in Slice A (the
    // "fields default when read_layout is empty" invariant).
    // The other four are reserved scaffolding for Slice B/C,
    // where projection logic provides values the validator can
    // re-derive and compare against.

    /// A record carries `read_layout == ""` (no paired-end
    /// projection requested) but one of the eight paired-end
    /// fields is set to a non-default value. The `reported`
    /// value is the offending field name (mirrors
    /// `OriginalVCallMismatch`'s string-payload shape); the
    /// `expected` value is the literal `"<default>"` token so a
    /// future v2 variant that wants to attach the actual
    /// default representation can extend the payload without
    /// breaking string consumers.
    PairedEndFieldWithoutLayout { field: PairedEndField },
    /// Reserved for Slice B: an R1/R2 window's `[start, end)`
    /// range falls outside `[0, sequence_length]` or is
    /// inverted. Not fired in Slice A.
    ReadWindowOutOfBounds {
        side: PairedEndRead,
        start: i64,
        end: i64,
        sequence_length: i64,
    },
    /// Reserved for Slice B: `r1_sequence` /
    /// `r2_sequence` disagrees with the re-derived window
    /// substring. Not fired in Slice A.
    ReadSequenceMismatch {
        side: PairedEndRead,
        reported: String,
        expected: String,
    },
    /// Reserved for Slice B/C: `insert_size` disagrees with the
    /// window geometry (audit §8 pins
    /// `insert_size == r2_end`). Not fired in Slice A.
    ReadInsertSizeMismatch { reported: i64, expected: i64 },
    /// Reserved for Slice B/C: `read_layout` carries an
    /// unsupported value (not `""` / `"paired_end"` /
    /// `"single_end"`). Not fired in Slice A.
    ReadLayoutMismatch { reported: String, expected: String },
}

/// Which paired-end field tripped a structured issue. Used by
/// `PairedEndFieldWithoutLayout` (Slice A) and reserved for the
/// per-field geometry checks Slice B/C will add.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PairedEndField {
    ReadLayout,
    R1Sequence,
    R2Sequence,
    R1Start,
    R1End,
    R2Start,
    R2End,
    InsertSize,
}

impl PairedEndField {
    /// Snake-case field name as it appears in the AIRR record
    /// dict; pinned by the Slice A contract test so a future
    /// rename has to come with an explicit `pin_scaffold_*` flip.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ReadLayout => "read_layout",
            Self::R1Sequence => "r1_sequence",
            Self::R2Sequence => "r2_sequence",
            Self::R1Start => "r1_start",
            Self::R1End => "r1_end",
            Self::R2Start => "r2_start",
            Self::R2End => "r2_end",
            Self::InsertSize => "insert_size",
        }
    }
}

/// Which side of a paired-end read tripped a structured issue.
/// Reserved for Slice B/C; declared in Slice A so the variant
/// shapes are stable from day one.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PairedEndRead {
    R1,
    R2,
}

impl PairedEndRead {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::R1 => "r1",
            Self::R2 => "r2",
        }
    }
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
    check_counters(record, outcome, refdata, &mut issues);
    check_junction(record, sim, refdata, &mut issues);
    check_allele_oracle(record, outcome, refdata, &mut issues);
    check_region_and_hypothesis_invariants(sim, &mut issues);
    check_paired_end_defaults(record, &mut issues);

    issues
}

// ──────────────────────────────────────────────────────────────────
// C6: Paired-end / read-layout invariants
//
// Slice A: enforces the no-layout default invariant (`read_layout
// == ""` ⇒ every paired-end field is at its default).
//
// Slice B: extends the dispatch with per-layout geometry checks.
// When `read_layout == "paired_end"`, the validator re-derives
// R1/R2 windows from the rules in `docs/paired_end_design.md` §8
// and surfaces any divergence as one of the four reserved variants
// (`ReadWindowOutOfBounds`, `ReadSequenceMismatch`,
// `ReadInsertSizeMismatch`). Unknown non-empty layouts surface as
// `ReadLayoutMismatch`. `"single_end"` is documented as reserved
// (§2.2) and treated as a no-op for now.
// ──────────────────────────────────────────────────────────────────

fn check_paired_end_defaults(
    record: &AirrRecord,
    issues: &mut Vec<RecordValidationIssue>,
) {
    match record.read_layout.as_str() {
        // Slice A: no-layout default invariant.
        "" => check_paired_end_default_values(record, issues),
        // Slice B: geometry against the projection rules.
        "paired_end" => check_paired_end_geometry(record, issues),
        // Reserved (per §2.2). A future slice may attach a
        // single-read geometry check; today we don't validate
        // anything beyond accepting the layout string.
        "single_end" => {}
        // Anything else is an unsupported layout value. Surface
        // the structured mismatch so a typo (`"pair_end"`) or a
        // refactor that introduced a new layout without wiring
        // it through the validator fails closed.
        _ => issues.push(RecordValidationIssue::ReadLayoutMismatch {
            reported: record.read_layout.clone(),
            expected: r#"one of: "", "paired_end", "single_end""#.to_string(),
        }),
    }
}

fn check_paired_end_default_values(
    record: &AirrRecord,
    issues: &mut Vec<RecordValidationIssue>,
) {
    if !record.r1_sequence.is_empty() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R1Sequence,
        });
    }
    if !record.r2_sequence.is_empty() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R2Sequence,
        });
    }
    if record.r1_start.is_some() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R1Start,
        });
    }
    if record.r1_end.is_some() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R1End,
        });
    }
    if record.r2_start.is_some() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R2Start,
        });
    }
    if record.r2_end.is_some() {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::R2End,
        });
    }
    if record.insert_size != 0 {
        issues.push(RecordValidationIssue::PairedEndFieldWithoutLayout {
            field: PairedEndField::InsertSize,
        });
    }
}

fn check_paired_end_geometry(
    record: &AirrRecord,
    issues: &mut Vec<RecordValidationIssue>,
) {
    let seq_len = record.sequence_length;
    let seq = &record.sequence;

    // R1 window: bounds + byte equality.
    match resolve_window(record.r1_start, record.r1_end, seq_len) {
        WindowResolution::Valid { start, end } => {
            let expected = seq[start as usize..end as usize].to_string();
            if record.r1_sequence != expected {
                issues.push(RecordValidationIssue::ReadSequenceMismatch {
                    side: PairedEndRead::R1,
                    reported: record.r1_sequence.clone(),
                    expected,
                });
            }
        }
        WindowResolution::OutOfBounds { start, end } => {
            issues.push(RecordValidationIssue::ReadWindowOutOfBounds {
                side: PairedEndRead::R1,
                start,
                end,
                sequence_length: seq_len,
            });
        }
    }

    // R2 window: bounds + reverse-complement byte equality + insert
    // size consistency. The audit pins
    // `insert_size == r2_end` (§8) — the only insert-size-mismatch
    // surface today.
    match resolve_window(record.r2_start, record.r2_end, seq_len) {
        WindowResolution::Valid { start, end } => {
            let r2_inner = &seq[start as usize..end as usize];
            let expected = super::sequence::reverse_complement(r2_inner);
            if record.r2_sequence != expected {
                issues.push(RecordValidationIssue::ReadSequenceMismatch {
                    side: PairedEndRead::R2,
                    reported: record.r2_sequence.clone(),
                    expected,
                });
            }
            if record.insert_size != end {
                issues.push(RecordValidationIssue::ReadInsertSizeMismatch {
                    reported: record.insert_size,
                    expected: end,
                });
            }
        }
        WindowResolution::OutOfBounds { start, end } => {
            issues.push(RecordValidationIssue::ReadWindowOutOfBounds {
                side: PairedEndRead::R2,
                start,
                end,
                sequence_length: seq_len,
            });
            // With R2 bounds unresolved we can't recompute the
            // expected insert size; skip the insert-size check
            // rather than fabricate a sentinel. The window
            // out-of-bounds issue is the actionable signal.
        }
    }
}

/// Result of resolving a paired-end window `(start, end)` against
/// the projected `sequence_length`. `OutOfBounds` collapses missing
/// coords (`None`) and out-of-range coords (negative, swapped,
/// past the end) into one variant whose `start`/`end` fields use
/// the sentinel `-1` for missing values — surfaces cleanly in the
/// structured-issue payload without inventing a new variant.
enum WindowResolution {
    Valid { start: i64, end: i64 },
    OutOfBounds { start: i64, end: i64 },
}

fn resolve_window(
    start: Option<i64>,
    end: Option<i64>,
    sequence_length: i64,
) -> WindowResolution {
    let start_val = start.unwrap_or(-1);
    let end_val = end.unwrap_or(-1);
    match (start, end) {
        (Some(s), Some(e)) if s >= 0 && e >= s && e <= sequence_length => {
            WindowResolution::Valid { start: s, end: e }
        }
        _ => WindowResolution::OutOfBounds {
            start: start_val,
            end: end_val,
        },
    }
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
    // Orientation drives the per-byte comparison rule via the
    // shared `matches_observed_with_orientation` primitive: under
    // `ReverseComplement` the observed byte is pre-complemented
    // before matching the allele's germline byte at the same
    // `germline_pos`. Defaults to `Forward` when the segment is
    // unassigned. See `scoring::observed_in_germline_orientation`
    // for the rationale.
    let orientation = assignment
        .map(|a| a.orientation)
        .unwrap_or(crate::assignment::SegmentOrientation::Forward);

    // Independent rescore via the shared scoring kernel: structural
    // region + NP-region extensions under the assigned allele's trim
    // caps. Mirrors `live_call::walker::call_from_region` so the
    // oracle and the walker agree on the tie-set under arbitrary
    // trim.
    let scores = score_alleles_with_extensions(
        sim,
        segment,
        allele_pool,
        region,
        trim_5_cap,
        trim_3_cap,
        orientation,
    );
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
