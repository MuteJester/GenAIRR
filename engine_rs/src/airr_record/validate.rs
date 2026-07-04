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

use crate::address::PrimeEnd;
use crate::ir::Segment;
use crate::pass::Outcome;
use crate::refdata::RefDataConfig;

use super::record::AirrRecord;

mod allele_oracle;
mod counters;
mod junction;
mod paired_end;
mod region_invariants;
mod structural;

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

    structural::check_structural(record, sim, &mut issues);
    counters::check_counters(record, outcome, refdata, &mut issues);
    junction::check_junction(record, sim, refdata, &mut issues);
    allele_oracle::check_allele_oracle(record, outcome, refdata, &mut issues);
    region_invariants::check_region_and_hypothesis_invariants(sim, &mut issues);
    paired_end::check_paired_end_defaults(record, &mut issues);

    issues
}
