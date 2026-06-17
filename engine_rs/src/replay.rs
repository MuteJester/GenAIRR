//! Trace consumer core for trace-injected replay (Option B, first
//! slice).
//!
//! ## What this provides
//!
//! - [`TraceCursor`]: an ordered cursor over `&[ChoiceRecord]` that
//!   passes consume from when running in **replay** mode.
//! - [`ReplayError`]: a closed enum of every way replay can fail
//!   (address mismatch, value-kind mismatch, exhausted trace,
//!   unused trailing records). Each variant carries enough context
//!   to produce a readable diagnostic without re-reading the trace.
//! - Typed consumption helpers: [`TraceCursor::expect_int`],
//!   [`TraceCursor::expect_base`], [`TraceCursor::expect_bases`],
//!   [`TraceCursor::expect_allele_id`], [`TraceCursor::expect_bool`].
//!   Each verifies the expected address against the cursor's next
//!   record, advances the cursor on success, and returns the typed
//!   payload.
//!
//! ## Architectural contract
//!
//! Replay mode is **strictly positional**: the i-th `record_choice`
//! a pass would have made in fresh-RNG mode must correspond to the
//! i-th record in the input trace. The cursor enforces this:
//!
//! - The expected address is the typed [`ChoiceAddress`] the pass
//!   would have recorded. The cursor compares it (via `Display`)
//!   to the recorded address string and refuses on mismatch.
//! - The expected value kind is encoded in the helper name
//!   (`expect_int` etc). The cursor refuses on kind mismatch.
//! - Sites that the pass would have skipped in fresh-RNG mode
//!   (e.g. permissive empty-support no-ops) must NOT appear in the
//!   trace, because the original run did not record them. Replay
//!   does not synthesise these slots.
//!
//! After execution completes the runtime calls
//! [`TraceCursor::assert_drained`] to detect a trace that's longer
//! than the plan consumes — typically caused by replaying against a
//! different plan or rerunning after a pass changed how many
//! records it emits.
//!
//! ## What this is *not*
//!
//! Not a sampler. The cursor reads from a pre-recorded buffer; it
//! has no view of the contract set, the RNG, or the simulation
//! state. Passes that want "consume from trace if available, else
//! sample fresh" branch on `ctx.replay_cursor.is_some()` at their
//! top-level execute path.

use crate::address::ChoiceAddress;
use crate::trace::{ChoiceRecord, ChoiceValue};

/// Closed enum of replay failure modes.
///
/// Variants are designed to render cleanly via `Display` without
/// needing the consumer to introspect the variant.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReplayError {
    /// Cursor was empty when a pass asked for the next record.
    Exhausted {
        position: usize,
        expected_address: String,
    },
    /// Cursor's next record had a different address than the pass
    /// expected. `position` is the 0-based index of the offending
    /// record in the input trace.
    AddressMismatch {
        position: usize,
        expected: String,
        got: String,
    },
    /// Cursor's next record had the expected address but a value
    /// of the wrong [`ChoiceValue`] variant.
    ValueKindMismatch {
        position: usize,
        address: String,
        expected_kind: &'static str,
        got_kind: &'static str,
    },
    /// Cursor still had records after the plan finished. The plan
    /// either changed since the trace was recorded, or replay
    /// targeted the wrong simulator. Returned by
    /// [`TraceCursor::assert_drained`]; never raised mid-pass.
    UnusedTrailingRecords {
        consumed: usize,
        remaining: usize,
    },
}

impl std::fmt::Display for ReplayError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Exhausted {
                position,
                expected_address,
            } => write!(
                f,
                "replay: trace exhausted at position {position}; \
                 pass asked for address {expected_address:?}",
            ),
            Self::AddressMismatch {
                position,
                expected,
                got,
            } => write!(
                f,
                "replay: address mismatch at position {position};\n  \
                 expected: {expected}\n  got:      {got}",
            ),
            Self::ValueKindMismatch {
                position,
                address,
                expected_kind,
                got_kind,
            } => write!(
                f,
                "replay: value-kind mismatch at position {position} \
                 (address {address:?}); expected {expected_kind}, got {got_kind}",
            ),
            Self::UnusedTrailingRecords { consumed, remaining } => write!(
                f,
                "replay: {remaining} unused record{s} remained after the \
                 plan consumed {consumed}",
                s = if *remaining == 1 { "" } else { "s" },
            ),
        }
    }
}

impl std::error::Error for ReplayError {}

/// Stable, human-readable kind name for a [`ChoiceValue`] variant.
/// Used in [`ReplayError::ValueKindMismatch`].
pub fn choice_value_kind(value: &ChoiceValue) -> &'static str {
    match value {
        ChoiceValue::Int(_) => "Int",
        ChoiceValue::Base(_) => "Base",
        ChoiceValue::Bases(_) => "Bases",
        ChoiceValue::AlleleId(_) => "AlleleId",
        ChoiceValue::Bool(_) => "Bool",
        ChoiceValue::Haplotype(_) => "Haplotype",
        ChoiceValue::GeneId(_) => "GeneId",
    }
}

/// Positional cursor over a pre-recorded sequence of
/// [`ChoiceRecord`]s.
///
/// Owns the records by clone (rather than borrowing) so the cursor
/// lifetime is unentangled from the records' source. This costs one
/// `Vec` clone at the start of a replay run — negligible compared to
/// the simulation itself — and avoids two lifetime parameters on
/// every `PassContext` borrow downstream.
///
/// Passes consume from this cursor through the typed `expect_*`
/// helpers; each consumption advances the position.
#[derive(Debug)]
pub struct TraceCursor {
    records: Vec<ChoiceRecord>,
    position: usize,
}

impl TraceCursor {
    /// Build a fresh cursor by cloning `records`. Position starts at
    /// zero.
    pub fn new(records: &[ChoiceRecord]) -> Self {
        Self {
            records: records.to_vec(),
            position: 0,
        }
    }

    /// Build a fresh cursor from an owned record vector. Used by the
    /// runtime when it already has ownership of the trace records.
    pub fn from_owned(records: Vec<ChoiceRecord>) -> Self {
        Self {
            records,
            position: 0,
        }
    }

    /// Number of records already consumed.
    pub fn position(&self) -> usize {
        self.position
    }

    /// Number of records remaining.
    pub fn remaining(&self) -> usize {
        self.records.len().saturating_sub(self.position)
    }

    /// `true` iff every record has been consumed.
    pub fn is_drained(&self) -> bool {
        self.position == self.records.len()
    }

    /// Read the next record's address without advancing the cursor.
    /// Returns `None` when the cursor is drained.
    ///
    /// Used by passes whose iteration shape doesn't 1:1 align with
    /// trace positions — e.g. `ContaminantPass` walks every pool
    /// site but only records (and consumes) at sites the original
    /// contract bundle admitted. The pass peeks the cursor at each
    /// iteration and decides whether to consume.
    pub fn peek_address(&self) -> Option<&str> {
        self.records.get(self.position).map(|r| r.address.as_str())
    }

    /// Verify the cursor was fully drained by the run. Returns
    /// `Err(ReplayError::UnusedTrailingRecords)` if records remain.
    /// Called by the runtime after the plan completes.
    pub fn assert_drained(&self) -> Result<(), ReplayError> {
        if self.is_drained() {
            Ok(())
        } else {
            Err(ReplayError::UnusedTrailingRecords {
                consumed: self.position,
                remaining: self.remaining(),
            })
        }
    }

    /// Internal: consume the next record, verifying it matches the
    /// expected address. On success the cursor advances and returns
    /// `(consumed_position, cloned_address, cloned_value)`. On
    /// failure the cursor stays put so the error diagnostic can
    /// identify the bad record.
    ///
    /// Returns cloned data (rather than a borrow of the cursor's
    /// stored record) so callers can read `self.position` /
    /// `self.records` without overlapping borrows. The clones are
    /// small (`String` + `ChoiceValue`, where the largest variant —
    /// `Bases` — is a few bases of `Vec<u8>` in practice).
    fn advance_with_address(
        &mut self,
        expected: ChoiceAddress,
    ) -> Result<(usize, String, ChoiceValue), ReplayError> {
        let idx = self.position;
        let (rec_address, rec_value) = match self.records.get(idx) {
            Some(rec) => (rec.address.clone(), rec.value.clone()),
            None => {
                return Err(ReplayError::Exhausted {
                    position: idx,
                    expected_address: expected.to_string(),
                });
            }
        };
        let expected_str = expected.to_string();
        if rec_address != expected_str {
            return Err(ReplayError::AddressMismatch {
                position: idx,
                expected: expected_str,
                got: rec_address,
            });
        }
        self.position = idx + 1;
        Ok((idx, rec_address, rec_value))
    }

    /// Consume the next record as a `ChoiceValue::Int`. Returns the
    /// raw `i64`; the caller is responsible for any range / sign
    /// validation domain-specific to the consuming pass.
    pub fn expect_int(&mut self, address: ChoiceAddress) -> Result<i64, ReplayError> {
        let (position, address_str, value) = self.advance_with_address(address)?;
        match value {
            ChoiceValue::Int(n) => Ok(n),
            other => Err(kind_mismatch(position, &address_str, "Int", &other)),
        }
    }

    /// Consume the next record as a `ChoiceValue::Base`.
    pub fn expect_base(&mut self, address: ChoiceAddress) -> Result<u8, ReplayError> {
        let (position, address_str, value) = self.advance_with_address(address)?;
        match value {
            ChoiceValue::Base(b) => Ok(b),
            other => Err(kind_mismatch(position, &address_str, "Base", &other)),
        }
    }

    /// Consume the next record as a `ChoiceValue::Bases`. Returns an
    /// owned `Vec<u8>`.
    pub fn expect_bases(&mut self, address: ChoiceAddress) -> Result<Vec<u8>, ReplayError> {
        let (position, address_str, value) = self.advance_with_address(address)?;
        match value {
            ChoiceValue::Bases(bs) => Ok(bs),
            other => Err(kind_mismatch(position, &address_str, "Bases", &other)),
        }
    }

    /// Consume the next record as a `ChoiceValue::AlleleId`. Returns
    /// the raw `u32` pool index.
    pub fn expect_allele_id(&mut self, address: ChoiceAddress) -> Result<u32, ReplayError> {
        let (position, address_str, value) = self.advance_with_address(address)?;
        match value {
            ChoiceValue::AlleleId(id) => Ok(id),
            other => Err(kind_mismatch(position, &address_str, "AlleleId", &other)),
        }
    }

    /// Consume the next record as a `ChoiceValue::Bool`.
    pub fn expect_bool(&mut self, address: ChoiceAddress) -> Result<bool, ReplayError> {
        let (position, address_str, value) = self.advance_with_address(address)?;
        match value {
            ChoiceValue::Bool(b) => Ok(b),
            other => Err(kind_mismatch(position, &address_str, "Bool", &other)),
        }
    }

    /// Consume the next record as a `ChoiceValue::Haplotype`.
    pub fn expect_haplotype(&mut self, address: ChoiceAddress) -> Result<u8, ReplayError> {
        let (position, address_str, value) = self.advance_with_address(address)?;
        match value {
            ChoiceValue::Haplotype(h) => Ok(h),
            other => Err(kind_mismatch(position, &address_str, "Haplotype", &other)),
        }
    }

    /// Consume the next record as a `ChoiceValue::GeneId`.
    pub fn expect_gene_id(&mut self, address: ChoiceAddress) -> Result<u32, ReplayError> {
        let (position, address_str, value) = self.advance_with_address(address)?;
        match value {
            ChoiceValue::GeneId(g) => Ok(g),
            other => Err(kind_mismatch(position, &address_str, "GeneId", &other)),
        }
    }
}

/// Build a `ValueKindMismatch` without re-borrowing the cursor. Free
/// function so the expect_* helpers can call it after the
/// `advance_with_address` borrow has already produced the record.
fn kind_mismatch(
    position: usize,
    address: &str,
    expected_kind: &'static str,
    got: &ChoiceValue,
) -> ReplayError {
    ReplayError::ValueKindMismatch {
        position,
        address: address.to_string(),
        expected_kind,
        got_kind: choice_value_kind(got),
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::address::{NpSegment, VdjSegment};
    use crate::assignment::TrimEnd;
    use crate::trace::Trace;

    fn make_records(entries: &[(ChoiceAddress, ChoiceValue)]) -> Trace {
        let mut t = Trace::new();
        for (addr, value) in entries {
            t.record_choice(*addr, value.clone());
        }
        t
    }

    // ── Position / remaining / drained ──────────────────────────────

    #[test]
    fn cursor_starts_at_position_zero() {
        let records: Vec<ChoiceRecord> = Vec::new();
        let cursor = TraceCursor::new(&records);
        assert_eq!(cursor.position(), 0);
        assert_eq!(cursor.remaining(), 0);
        assert!(cursor.is_drained());
        assert!(cursor.assert_drained().is_ok());
    }

    #[test]
    fn cursor_advances_after_successful_consume() {
        let trace = make_records(&[
            (
                ChoiceAddress::SampleAllele(VdjSegment::V),
                ChoiceValue::AlleleId(7),
            ),
            (
                ChoiceAddress::SampleAllele(VdjSegment::J),
                ChoiceValue::AlleleId(3),
            ),
        ]);
        let mut cursor = TraceCursor::new(trace.choices());
        assert_eq!(cursor.position(), 0);

        let id = cursor
            .expect_allele_id(ChoiceAddress::SampleAllele(VdjSegment::V))
            .unwrap();
        assert_eq!(id, 7);
        assert_eq!(cursor.position(), 1);
        assert_eq!(cursor.remaining(), 1);

        let id = cursor
            .expect_allele_id(ChoiceAddress::SampleAllele(VdjSegment::J))
            .unwrap();
        assert_eq!(id, 3);
        assert!(cursor.is_drained());
    }

    // ── Each typed consumer ─────────────────────────────────────────

    #[test]
    fn expect_int_returns_value() {
        let trace = make_records(&[(
            ChoiceAddress::NpLength(NpSegment::Np1),
            ChoiceValue::Int(5),
        )]);
        let mut cursor = TraceCursor::new(trace.choices());
        assert_eq!(
            cursor.expect_int(ChoiceAddress::NpLength(NpSegment::Np1)).unwrap(),
            5,
        );
    }

    #[test]
    fn expect_base_returns_value() {
        let trace = make_records(&[(
            ChoiceAddress::NpBase {
                segment: NpSegment::Np1,
                index: 0,
            },
            ChoiceValue::Base(b'A'),
        )]);
        let mut cursor = TraceCursor::new(trace.choices());
        assert_eq!(
            cursor
                .expect_base(ChoiceAddress::NpBase {
                    segment: NpSegment::Np1,
                    index: 0,
                })
                .unwrap(),
            b'A',
        );
    }

    #[test]
    fn expect_allele_id_returns_value() {
        let trace = make_records(&[(
            ChoiceAddress::SampleAllele(VdjSegment::V),
            ChoiceValue::AlleleId(42),
        )]);
        let mut cursor = TraceCursor::new(trace.choices());
        assert_eq!(
            cursor
                .expect_allele_id(ChoiceAddress::SampleAllele(VdjSegment::V))
                .unwrap(),
            42,
        );
    }

    #[test]
    fn expect_bases_returns_clone() {
        let trace = make_records(&[(
            ChoiceAddress::NpBase {
                segment: NpSegment::Np2,
                index: 0,
            },
            ChoiceValue::Bases(b"ACGT".to_vec()),
        )]);
        let mut cursor = TraceCursor::new(trace.choices());
        let bs = cursor
            .expect_bases(ChoiceAddress::NpBase {
                segment: NpSegment::Np2,
                index: 0,
            })
            .unwrap();
        assert_eq!(bs, b"ACGT");
    }

    #[test]
    fn expect_bool_returns_value() {
        let trace = make_records(&[(
            ChoiceAddress::CorruptRevCompApplied,
            ChoiceValue::Bool(true),
        )]);
        let mut cursor = TraceCursor::new(trace.choices());
        assert!(cursor
            .expect_bool(ChoiceAddress::CorruptRevCompApplied)
            .unwrap());
    }

    // ── Error variants ─────────────────────────────────────────────

    #[test]
    fn expect_int_on_empty_cursor_returns_exhausted() {
        let records: Vec<ChoiceRecord> = Vec::new();
        let mut cursor = TraceCursor::new(&records);
        let err = cursor
            .expect_int(ChoiceAddress::NpLength(NpSegment::Np1))
            .unwrap_err();
        match err {
            ReplayError::Exhausted {
                position,
                expected_address,
            } => {
                assert_eq!(position, 0);
                assert_eq!(expected_address, "np.np1.length");
            }
            other => panic!("expected Exhausted, got {other:?}"),
        }
    }

    #[test]
    fn expect_with_wrong_address_returns_address_mismatch() {
        let trace = make_records(&[(
            ChoiceAddress::SampleAllele(VdjSegment::V),
            ChoiceValue::AlleleId(0),
        )]);
        let mut cursor = TraceCursor::new(trace.choices());
        let err = cursor
            .expect_allele_id(ChoiceAddress::SampleAllele(VdjSegment::J))
            .unwrap_err();
        match err {
            ReplayError::AddressMismatch {
                position,
                expected,
                got,
            } => {
                assert_eq!(position, 0);
                assert_eq!(expected, "sample_allele.j");
                assert_eq!(got, "sample_allele.v");
            }
            other => panic!("expected AddressMismatch, got {other:?}"),
        }
        // Cursor must not have advanced on the failed call.
        assert_eq!(cursor.position(), 0);
    }

    #[test]
    fn expect_with_wrong_kind_returns_value_kind_mismatch() {
        let trace = make_records(&[(
            ChoiceAddress::SampleAllele(VdjSegment::V),
            ChoiceValue::AlleleId(7),
        )]);
        let mut cursor = TraceCursor::new(trace.choices());
        let err = cursor
            .expect_int(ChoiceAddress::SampleAllele(VdjSegment::V))
            .unwrap_err();
        match err {
            ReplayError::ValueKindMismatch {
                position,
                address,
                expected_kind,
                got_kind,
            } => {
                assert_eq!(position, 0);
                assert_eq!(address, "sample_allele.v");
                assert_eq!(expected_kind, "Int");
                assert_eq!(got_kind, "AlleleId");
            }
            other => panic!("expected ValueKindMismatch, got {other:?}"),
        }
    }

    #[test]
    fn assert_drained_with_remaining_returns_unused_trailing_records() {
        let trace = make_records(&[
            (
                ChoiceAddress::SampleAllele(VdjSegment::V),
                ChoiceValue::AlleleId(0),
            ),
            (
                ChoiceAddress::SampleAllele(VdjSegment::J),
                ChoiceValue::AlleleId(1),
            ),
        ]);
        let mut cursor = TraceCursor::new(trace.choices());
        let _ = cursor
            .expect_allele_id(ChoiceAddress::SampleAllele(VdjSegment::V))
            .unwrap();
        let err = cursor.assert_drained().unwrap_err();
        match err {
            ReplayError::UnusedTrailingRecords { consumed, remaining } => {
                assert_eq!(consumed, 1);
                assert_eq!(remaining, 1);
            }
            other => panic!("expected UnusedTrailingRecords, got {other:?}"),
        }
    }

    // ── Display ────────────────────────────────────────────────────

    #[test]
    fn display_renders_each_variant() {
        let e1 = ReplayError::Exhausted {
            position: 3,
            expected_address: "trim.v_3".into(),
        };
        assert!(format!("{e1}").contains("exhausted"));
        assert!(format!("{e1}").contains("trim.v_3"));

        let e2 = ReplayError::AddressMismatch {
            position: 1,
            expected: "sample_allele.v".into(),
            got: "sample_allele.j".into(),
        };
        let s = format!("{e2}");
        assert!(s.contains("expected"));
        assert!(s.contains("got"));

        let e3 = ReplayError::ValueKindMismatch {
            position: 0,
            address: "trim.v_3".into(),
            expected_kind: "Int",
            got_kind: "Base",
        };
        assert!(format!("{e3}").contains("value-kind mismatch"));

        let e4 = ReplayError::UnusedTrailingRecords {
            consumed: 5,
            remaining: 3,
        };
        let s = format!("{e4}");
        assert!(s.contains("3 unused"));
    }

    // Compile-touch the trim end address spelling so the test
    // module's `TrimEnd` import stays load-bearing if the surface
    // grows.
    #[test]
    fn expect_address_for_trim_uses_canonical_spelling() {
        let trace = make_records(&[(
            ChoiceAddress::Trim {
                segment: VdjSegment::V,
                end: TrimEnd::Three,
            },
            ChoiceValue::Int(2),
        )]);
        let mut cursor = TraceCursor::new(trace.choices());
        let n = cursor
            .expect_int(ChoiceAddress::Trim {
                segment: VdjSegment::V,
                end: TrimEnd::Three,
            })
            .unwrap();
        assert_eq!(n, 2);
    }
}
