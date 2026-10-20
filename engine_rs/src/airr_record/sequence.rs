use super::AirrRecord;
use crate::ir::Simulation;

pub(super) fn pool_bases(sim: &Simulation) -> Vec<u8> {
    sim.pool.as_slice().iter().map(|n| n.base).collect()
}

pub(super) fn bytes_to_string(bytes: &[u8]) -> String {
    // Bases are ASCII; lossy decode is safe.
    String::from_utf8_lossy(bytes).into_owned()
}

pub(super) fn apply_rev_comp_projection(rec: &mut AirrRecord) {
    // AIRR annotations are always expressed in the coding orientation.
    // rev_comp describes the relationship between that orientation and
    // the observed sequence, not a different biological rearrangement.
    rec.rev_comp = true;
    rec.sequence = reverse_complement(&rec.sequence);
}

pub(super) fn reverse_complement(s: &str) -> String {
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
// Paired-end / read-layout projection kernel (Slice B)
//
// Pure function: takes an already-built `AirrRecord` plus explicit
// layout parameters and returns the projected field values as a
// `PairedEndProjection`. No mutation, no RNG, no trace access. The
// caller is responsible for applying the projection to a record
// (Slice C/D wires the per-record parameter source through the
// trace).
//
// Coordinate space: `record.sequence` is the only substrate. End-
// loss + `random_strand_orientation` have already finalised the
// molecule before paired-end window selection — paired-end is the
// last projection layer per `docs/paired_end_design.md` §6 / §7.
// ──────────────────────────────────────────────────────────────────

use super::validate::PairedEndRead;

/// Projected field values for the `read_layout = "paired_end"`
/// case. Returned by [`project_paired_end_layout`]; the caller
/// applies these to an `AirrRecord` to convert it from the
/// no-layout default state to a populated paired-end record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PairedEndProjection {
    pub read_layout: String,
    pub r1_sequence: String,
    pub r2_sequence: String,
    pub r1_start: Option<i64>,
    pub r1_end: Option<i64>,
    pub r2_start: Option<i64>,
    pub r2_end: Option<i64>,
    pub insert_size: i64,
}

/// Structured projection error. Each variant names a specific
/// invariant violation so a future caller (Slice C+) can surface
/// the failure as a `PassError` without re-deriving the cause.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum PairedEndProjectionError {
    /// `r1_length` or `r2_length` is non-positive. Read lengths
    /// must be `> 0`; the no-layout default (`""`) is the
    /// canonical "no reads" representation.
    NonPositiveLength { side: PairedEndRead, length: i64 },
    /// `insert_size` falls outside `[0, sequence_length]`. The
    /// audit's strict-truncation rule (§6.1): a sequencer can't
    /// read past the fragment 3' end.
    InsertSizeOutOfBounds {
        insert_size: i64,
        sequence_length: i64,
    },
    /// `r1_length > insert_size` or `r2_length > insert_size`.
    /// Either alone would run past the fragment 3' end; the
    /// degenerate case the audit §6.2 explicitly rejects.
    /// (`r1_length + r2_length > insert_size` — read overlap —
    /// is *allowed* per §6.2 and not surfaced here.)
    ReadExceedsInsert {
        side: PairedEndRead,
        length: i64,
        insert_size: i64,
    },
}

/// Project a paired-end read layout onto the supplied record's
/// `sequence`.
///
/// Semantics (per `docs/paired_end_design.md` §6 + §8):
///
/// - `read_layout = "paired_end"`
/// - `r1_start = 0`, `r1_end = r1_length`
/// - `r2_start = insert_size - r2_length`, `r2_end = insert_size`
/// - `r1_sequence = sequence[r1_start..r1_end]`
/// - `r2_sequence = reverse_complement(sequence[r2_start..r2_end])`
/// - `insert_size` is the supplied parameter, unchanged
///
/// The helper validates input bounds before slicing; out-of-range
/// inputs surface as a structured `PairedEndProjectionError` so
/// the caller (Slice C+) can map them onto a `PassError`.
///
/// **No `Record` mutation.** Slice B is a pure projection kernel
/// — the caller applies the returned struct to a record (or
/// compares against an already-applied record in the validator).
pub(crate) fn project_paired_end_layout(
    record: &super::AirrRecord,
    r1_length: i64,
    r2_length: i64,
    insert_size: i64,
) -> Result<PairedEndProjection, PairedEndProjectionError> {
    let seq_len = record.sequence_length;

    if r1_length <= 0 {
        return Err(PairedEndProjectionError::NonPositiveLength {
            side: PairedEndRead::R1,
            length: r1_length,
        });
    }
    if r2_length <= 0 {
        return Err(PairedEndProjectionError::NonPositiveLength {
            side: PairedEndRead::R2,
            length: r2_length,
        });
    }
    if insert_size < 0 || insert_size > seq_len {
        return Err(PairedEndProjectionError::InsertSizeOutOfBounds {
            insert_size,
            sequence_length: seq_len,
        });
    }
    if r1_length > insert_size {
        return Err(PairedEndProjectionError::ReadExceedsInsert {
            side: PairedEndRead::R1,
            length: r1_length,
            insert_size,
        });
    }
    if r2_length > insert_size {
        return Err(PairedEndProjectionError::ReadExceedsInsert {
            side: PairedEndRead::R2,
            length: r2_length,
            insert_size,
        });
    }

    let seq = &record.sequence;
    let r1_start: i64 = 0;
    let r1_end: i64 = r1_length;
    let r2_start: i64 = insert_size - r2_length;
    let r2_end: i64 = insert_size;

    let r1_sequence = seq[r1_start as usize..r1_end as usize].to_string();
    let r2_inner = &seq[r2_start as usize..r2_end as usize];
    let r2_sequence = reverse_complement(r2_inner);

    Ok(PairedEndProjection {
        read_layout: "paired_end".to_string(),
        r1_sequence,
        r2_sequence,
        r1_start: Some(r1_start),
        r1_end: Some(r1_end),
        r2_start: Some(r2_start),
        r2_end: Some(r2_end),
        insert_size,
    })
}
