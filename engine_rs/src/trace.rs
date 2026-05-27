//! Addressed stochastic trace — the record of every random choice
//! the engine made during a simulation.
//!
//! ## What this is
//!
//! Per D3, every random draw the engine makes carries an **address**:
//! a hierarchical dotted string (e.g. `"trim.v_3"`,
//! `"np.np1.length"`, `"mutate.s5f.site[14]"`) that identifies the
//! choice independently of physical position. Per D9, the history
//! exposed to introspection includes this trace — replay, diff, and
//! debugging all operate on it.
//!
//! ## What it is not
//!
//! Not part of `Simulation`. Trace lives as a side-channel maintained
//! by the runtime alongside the IR revision chain. A persistent
//! Trace-on-Simulation would have to clone a growing Vec at every
//! pass — quadratic cost for no architectural benefit. Keeping the
//! trace separate also matches Gen.jl's separation of "the generative
//! function's IR" from "the trace of choices made during one run."
//!
//! ## Scope
//!
//! Just the data structures. No pass integration yet (that's B.2).
//! No serialization — trace is in-memory only at this stage.

use crate::address::ChoiceAddress;
use serde::{Deserialize, Serialize};

/// The typed value of a single random choice.
///
/// Intentionally a closed enum so contracts and replayers can match
/// exhaustively on the cases. New variants get added when new choice
/// kinds appear in the engine; the compiler refuses to build until
/// every consumer handles them.
///
/// # Serialized form
///
/// Each variant is tagged with `"kind"` and carries its data under
/// `"value"`. Byte payloads serialise as ASCII strings rather than
/// raw integers so JSON dumps stay readable — `Base(b'A')` lands as
/// `{"kind":"Base","value":"A"}`, not `{"kind":"Base","value":65}`.
/// The custom impls handle this; the `derive` would emit raw bytes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ChoiceValue {
    /// A signed integer. Used for length samples (NP length, trim
    /// amount, indel count, etc.) and for any choice naturally
    /// expressed as a number.
    Int(i64),

    /// A single nucleotide base byte (`b'A'`, `b'C'`, `b'G'`,
    /// `b'T'`, plus lowercase / `b'N'` variants).
    Base(u8),

    /// A sequence of nucleotide base bytes — used when a single
    /// pass samples a multi-base block (e.g., NP region bases).
    Bases(Vec<u8>),

    /// An allele identifier, indexed into the `RefDataConfig`'s
    /// allele table.
    AlleleId(u32),

    /// A boolean choice (e.g., D inversion: yes/no, receptor
    /// revision: yes/no, contaminant injection: yes/no).
    Bool(bool),
}

/// On-disk discriminant tag for `ChoiceValue`. Lives as a separate
/// enum so the wire format is a stable string set independent of the
/// Rust variant names — renaming a Rust variant cannot accidentally
/// change the JSON shape; only changing this enum can.
#[derive(Serialize, Deserialize)]
enum ChoiceValueWire {
    Int(i64),
    /// Single ASCII character.
    Base(char),
    /// ASCII string. Validated as ASCII at parse time.
    Bases(String),
    AlleleId(u32),
    Bool(bool),
}

impl Serialize for ChoiceValue {
    fn serialize<S: serde::Serializer>(&self, ser: S) -> Result<S::Ok, S::Error> {
        let wire = match *self {
            ChoiceValue::Int(n) => ChoiceValueWire::Int(n),
            ChoiceValue::Base(b) => ChoiceValueWire::Base(b as char),
            ChoiceValue::Bases(ref bs) => ChoiceValueWire::Bases(
                std::str::from_utf8(bs)
                    .map_err(|e| serde::ser::Error::custom(format!("non-ASCII bases: {e}")))?
                    .to_string(),
            ),
            ChoiceValue::AlleleId(id) => ChoiceValueWire::AlleleId(id),
            ChoiceValue::Bool(b) => ChoiceValueWire::Bool(b),
        };
        wire.serialize(ser)
    }
}

impl<'de> Deserialize<'de> for ChoiceValue {
    fn deserialize<D: serde::Deserializer<'de>>(de: D) -> Result<Self, D::Error> {
        let wire = ChoiceValueWire::deserialize(de)?;
        Ok(match wire {
            ChoiceValueWire::Int(n) => ChoiceValue::Int(n),
            ChoiceValueWire::Base(c) => {
                if !c.is_ascii() {
                    return Err(serde::de::Error::custom(format!(
                        "ChoiceValue::Base must be ASCII, got {c:?}"
                    )));
                }
                ChoiceValue::Base(c as u8)
            }
            ChoiceValueWire::Bases(s) => {
                if !s.is_ascii() {
                    return Err(serde::de::Error::custom(
                        "ChoiceValue::Bases must be ASCII",
                    ));
                }
                ChoiceValue::Bases(s.into_bytes())
            }
            ChoiceValueWire::AlleleId(id) => ChoiceValue::AlleleId(id),
            ChoiceValueWire::Bool(b) => ChoiceValue::Bool(b),
        })
    }
}

/// One entry in the trace: the address at which a choice was made
/// and the value that was sampled there.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ChoiceRecord {
    pub address: String,
    pub value: ChoiceValue,
}

impl ChoiceRecord {
    pub fn new(address: impl Into<String>, value: ChoiceValue) -> Self {
        Self {
            address: address.into(),
            value,
        }
    }

    /// Parse this record's persisted string address into the typed
    /// built-in address form, if it is one of the engine's canonical
    /// stochastic choices.
    pub fn choice_address(&self) -> Option<ChoiceAddress> {
        ChoiceAddress::parse(&self.address)
    }
}

/// The append-only record of every choice made during one simulation.
///
/// Maintained by the runtime as a side-channel; not part of any
/// `Simulation` IR revision. Once a simulation completes the trace
/// is delivered to the caller (via the `History` API in B.2 or via
/// `with_history=True` in the public `.run()` from D9) and dropped
/// otherwise.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Trace {
    choices: Vec<ChoiceRecord>,
}

impl Trace {
    /// Empty trace.
    pub fn new() -> Self {
        Self::default()
    }

    /// Empty trace pre-allocated for an expected number of choices.
    /// Used to avoid reallocations when the rough count is known
    /// from the active pass plan.
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            choices: Vec::with_capacity(cap),
        }
    }

    /// Number of choices recorded so far.
    pub fn len(&self) -> usize {
        self.choices.len()
    }

    /// Whether no choices have been recorded yet.
    pub fn is_empty(&self) -> bool {
        self.choices.is_empty()
    }

    /// Append a choice to the trace. The runtime calls this once
    /// per random draw inside a sampling pass. Address ownership
    /// is taken by the trace (the runtime usually has a `&'static
    /// str` constant for static addresses or a freshly-formatted
    /// `String` for indexed ones).
    pub fn record(&mut self, address: impl Into<String>, value: ChoiceValue) {
        self.choices.push(ChoiceRecord::new(address, value));
    }

    /// Append a choice at a typed built-in address.
    ///
    /// This preserves the persisted trace representation as the
    /// canonical string spelling while letting new internal call
    /// sites avoid hand-built address strings.
    pub fn record_choice(&mut self, address: ChoiceAddress, value: ChoiceValue) {
        self.record(address, value);
    }

    /// Atomically append a completed per-pass trace delta.
    ///
    /// The compiled simulator builds trace entries in a local
    /// per-pass buffer. Only after the pass result has passed all
    /// execution fences does it move that delta into the run trace.
    /// This keeps trace/history/state commits aligned.
    pub fn append_delta(&mut self, mut delta: Trace) {
        self.choices.append(&mut delta.choices);
    }

    /// All recorded choices in chronological order.
    pub fn choices(&self) -> &[ChoiceRecord] {
        &self.choices
    }

    /// Find the first choice recorded at the given address, if any.
    /// Address comparison is exact-string equality.
    pub fn find(&self, address: &str) -> Option<&ChoiceRecord> {
        self.choices.iter().find(|c| c.address == address)
    }

    /// Find the first choice recorded at the given typed built-in
    /// address, if any.
    pub fn find_choice(&self, address: ChoiceAddress) -> Option<&ChoiceRecord> {
        let address = address.to_string();
        self.find(&address)
    }

    /// Find every choice recorded at addresses with the given prefix.
    /// Useful for queries like "all S5F choices" via
    /// `prefix_query("mutate.s5f.")`. Returns choices in chronological
    /// order, not de-duplicated.
    ///
    /// Lifetime note: the returned iterator borrows from both `self`
    /// and `prefix`, so both must outlive the iterator. The named
    /// lifetime `'a` makes that explicit.
    pub fn prefix_query<'a>(
        &'a self,
        prefix: &'a str,
    ) -> impl Iterator<Item = &'a ChoiceRecord> + 'a {
        self.choices
            .iter()
            .filter(move |c| c.address.starts_with(prefix))
    }

    /// Number of choices made at addresses with the given prefix.
    pub fn prefix_count(&self, prefix: &str) -> usize {
        self.prefix_query(prefix).count()
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::address::{ChoiceAddress, NpSegment, VdjSegment};
    use crate::assignment::TrimEnd;

    #[test]
    fn trace_starts_empty() {
        let t = Trace::new();
        assert!(t.is_empty());
        assert_eq!(t.len(), 0);
        assert!(t.choices().is_empty());
        assert!(t.find("any.address").is_none());
        assert_eq!(t.prefix_count("anything."), 0);
    }

    #[test]
    fn trace_record_appends_in_order() {
        let mut t = Trace::new();
        t.record("trim.v_3", ChoiceValue::Int(2));
        t.record("np.np1.length", ChoiceValue::Int(7));
        t.record("np.np1.bases", ChoiceValue::Bases(b"acgtcga".to_vec()));

        assert_eq!(t.len(), 3);
        assert!(!t.is_empty());

        assert_eq!(t.choices()[0].address, "trim.v_3");
        assert_eq!(t.choices()[0].value, ChoiceValue::Int(2));

        assert_eq!(t.choices()[1].address, "np.np1.length");
        assert_eq!(t.choices()[1].value, ChoiceValue::Int(7));

        assert_eq!(t.choices()[2].address, "np.np1.bases");
        assert_eq!(
            t.choices()[2].value,
            ChoiceValue::Bases(b"acgtcga".to_vec())
        );
    }

    #[test]
    fn trace_append_delta_moves_records_in_order() {
        let mut run_trace = Trace::new();
        run_trace.record("first.choice", ChoiceValue::Int(1));

        let mut delta = Trace::new();
        delta.record("second.choice", ChoiceValue::Base(b'A'));
        delta.record("third.choice", ChoiceValue::Bool(true));

        run_trace.append_delta(delta);

        assert_eq!(run_trace.len(), 3);
        assert_eq!(run_trace.choices()[0].address, "first.choice");
        assert_eq!(run_trace.choices()[1].address, "second.choice");
        assert_eq!(run_trace.choices()[2].address, "third.choice");
    }

    #[test]
    fn trace_find_by_exact_address() {
        let mut t = Trace::new();
        t.record("trim.v_3", ChoiceValue::Int(4));
        t.record("trim.j_5", ChoiceValue::Int(1));

        let v3 = t.find("trim.v_3").expect("trim.v_3 should exist");
        assert_eq!(v3.value, ChoiceValue::Int(4));

        let j5 = t.find("trim.j_5").expect("trim.j_5 should exist");
        assert_eq!(j5.value, ChoiceValue::Int(1));

        assert!(t.find("trim.v_5").is_none());
    }

    #[test]
    fn trace_find_returns_first_match_only() {
        // Some addresses (e.g. retried draws in future phases) may end
        // up recorded twice. `find` returns the first chronologically.
        let mut t = Trace::new();
        t.record("retry.same_addr", ChoiceValue::Int(1));
        t.record("retry.same_addr", ChoiceValue::Int(2));

        let first = t.find("retry.same_addr").unwrap();
        assert_eq!(first.value, ChoiceValue::Int(1));
    }

    #[test]
    fn trace_typed_address_helpers_preserve_string_storage() {
        let mut t = Trace::new();
        t.record_choice(
            ChoiceAddress::Trim {
                segment: VdjSegment::V,
                end: TrimEnd::Three,
            },
            ChoiceValue::Int(4),
        );
        t.record_choice(
            ChoiceAddress::NpBase {
                segment: NpSegment::Np1,
                index: 2,
            },
            ChoiceValue::Base(b'G'),
        );

        assert_eq!(t.choices()[0].address, "trim.v_3");
        assert_eq!(
            t.choices()[0].choice_address(),
            Some(ChoiceAddress::Trim {
                segment: VdjSegment::V,
                end: TrimEnd::Three,
            })
        );

        let rec = t
            .find_choice(ChoiceAddress::NpBase {
                segment: NpSegment::Np1,
                index: 2,
            })
            .expect("typed NP base address should be present");
        assert_eq!(rec.address, "np.np1.bases[2]");
        assert_eq!(rec.value, ChoiceValue::Base(b'G'));
    }

    #[test]
    fn choice_record_typed_address_returns_none_for_custom_addresses() {
        let rec = ChoiceRecord::new("custom.choice", ChoiceValue::Bool(true));
        assert_eq!(rec.choice_address(), None);
    }

    #[test]
    fn trace_prefix_query_iterates_matching_choices() {
        let mut t = Trace::new();
        t.record("mutate.s5f.site[0]", ChoiceValue::Int(12));
        t.record("mutate.s5f.base[0]", ChoiceValue::Base(b'G'));
        t.record("mutate.s5f.site[1]", ChoiceValue::Int(45));
        t.record("mutate.s5f.base[1]", ChoiceValue::Base(b'A'));
        t.record("trim.v_3", ChoiceValue::Int(3));

        let s5f: Vec<&ChoiceRecord> = t.prefix_query("mutate.s5f.").collect();
        assert_eq!(s5f.len(), 4);

        let s5f_sites: Vec<&ChoiceRecord> = t.prefix_query("mutate.s5f.site").collect();
        assert_eq!(s5f_sites.len(), 2);

        let trims: Vec<&ChoiceRecord> = t.prefix_query("trim.").collect();
        assert_eq!(trims.len(), 1);

        let nothing: Vec<&ChoiceRecord> = t.prefix_query("corrupt.").collect();
        assert_eq!(nothing.len(), 0);
    }

    #[test]
    fn trace_prefix_count_matches_iter_len() {
        let mut t = Trace::new();
        t.record("mutate.s5f.site[0]", ChoiceValue::Int(0));
        t.record("mutate.s5f.site[1]", ChoiceValue::Int(1));
        t.record("mutate.s5f.site[2]", ChoiceValue::Int(2));

        assert_eq!(t.prefix_count("mutate.s5f."), 3);
        assert_eq!(t.prefix_count("mutate."), 3);
        assert_eq!(t.prefix_count(""), 3); // empty prefix matches all
        assert_eq!(t.prefix_count("none."), 0);
    }

    #[test]
    fn trace_with_capacity_does_not_change_observable_state() {
        let t = Trace::with_capacity(64);
        assert!(t.is_empty());
        assert_eq!(t.len(), 0);
    }

    #[test]
    fn choice_value_variants_compare_structurally() {
        // Pin the exhaustiveness contract for the typed enum: every
        // variant must roundtrip through PartialEq cleanly.
        assert_eq!(ChoiceValue::Int(7), ChoiceValue::Int(7));
        assert_ne!(ChoiceValue::Int(7), ChoiceValue::Int(8));

        assert_eq!(ChoiceValue::Base(b'C'), ChoiceValue::Base(b'C'));
        assert_ne!(ChoiceValue::Base(b'C'), ChoiceValue::Base(b'G'));

        assert_eq!(
            ChoiceValue::Bases(b"acgt".to_vec()),
            ChoiceValue::Bases(b"acgt".to_vec())
        );
        assert_ne!(
            ChoiceValue::Bases(b"acgt".to_vec()),
            ChoiceValue::Bases(b"acg".to_vec())
        );

        assert_eq!(ChoiceValue::AlleleId(42), ChoiceValue::AlleleId(42));
        assert_ne!(ChoiceValue::AlleleId(42), ChoiceValue::AlleleId(43));

        assert_eq!(ChoiceValue::Bool(true), ChoiceValue::Bool(true));
        assert_ne!(ChoiceValue::Bool(true), ChoiceValue::Bool(false));

        // Cross-variant inequality.
        assert_ne!(ChoiceValue::Int(0), ChoiceValue::Bool(false));
        assert_ne!(ChoiceValue::Base(b'A'), ChoiceValue::AlleleId(b'A' as u32));
    }

    #[test]
    fn trace_supports_dynamically_built_indexed_addresses() {
        // D3 calls out that indexed addresses use `format!` at
        // construction time. Confirm that pattern works through the
        // `impl Into<String>` API.
        let mut t = Trace::new();
        for i in 0..5 {
            t.record(
                format!("mutate.s5f.site[{}]", i),
                ChoiceValue::Int(i as i64),
            );
        }
        assert_eq!(t.prefix_count("mutate.s5f.site["), 5);
        for i in 0..5 {
            let addr = format!("mutate.s5f.site[{}]", i);
            let rec = t.find(&addr).expect("indexed address should be present");
            assert_eq!(rec.value, ChoiceValue::Int(i as i64));
        }
    }

    // ── serde round-trip ──────────────────────────────────────────

    #[test]
    fn choice_value_int_round_trips_through_json() {
        let v = ChoiceValue::Int(42);
        let s = serde_json::to_string(&v).unwrap();
        assert_eq!(s, r#"{"Int":42}"#);
        let back: ChoiceValue = serde_json::from_str(&s).unwrap();
        assert_eq!(back, v);
    }

    #[test]
    fn choice_value_base_serializes_as_character_string() {
        let v = ChoiceValue::Base(b'A');
        let s = serde_json::to_string(&v).unwrap();
        assert_eq!(s, r#"{"Base":"A"}"#);
        let back: ChoiceValue = serde_json::from_str(&s).unwrap();
        assert_eq!(back, v);
    }

    #[test]
    fn choice_value_lowercase_base_round_trips() {
        // Quality errors write lowercase bases — must round-trip.
        let v = ChoiceValue::Base(b'a');
        let s = serde_json::to_string(&v).unwrap();
        assert_eq!(s, r#"{"Base":"a"}"#);
        let back: ChoiceValue = serde_json::from_str(&s).unwrap();
        assert_eq!(back, v);
    }

    #[test]
    fn choice_value_bases_serializes_as_string() {
        let v = ChoiceValue::Bases(b"ACGTacgt".to_vec());
        let s = serde_json::to_string(&v).unwrap();
        assert_eq!(s, r#"{"Bases":"ACGTacgt"}"#);
        let back: ChoiceValue = serde_json::from_str(&s).unwrap();
        assert_eq!(back, v);
    }

    #[test]
    fn choice_value_allele_id_and_bool_round_trip() {
        for v in [
            ChoiceValue::AlleleId(7),
            ChoiceValue::Bool(true),
            ChoiceValue::Bool(false),
        ] {
            let s = serde_json::to_string(&v).unwrap();
            let back: ChoiceValue = serde_json::from_str(&s).unwrap();
            assert_eq!(back, v);
        }
    }

    #[test]
    fn choice_record_round_trips_through_json() {
        let rec = ChoiceRecord::new("mutate.uniform.base[0]", ChoiceValue::Base(b'G'));
        let s = serde_json::to_string(&rec).unwrap();
        let back: ChoiceRecord = serde_json::from_str(&s).unwrap();
        assert_eq!(back, rec);
    }

    #[test]
    fn trace_round_trips_through_json_preserving_order() {
        let mut t = Trace::new();
        t.record("trim.v_3", ChoiceValue::Int(2));
        t.record("np.np1.length", ChoiceValue::Int(5));
        t.record("np.np1.bases[0]", ChoiceValue::Base(b'A'));
        t.record("np.np1.bases[1]", ChoiceValue::Base(b'C'));
        t.record("sample_allele.v", ChoiceValue::AlleleId(12));

        let s = serde_json::to_string(&t).unwrap();
        let back: Trace = serde_json::from_str(&s).unwrap();

        assert_eq!(back.len(), t.len());
        for (a, b) in t.choices().iter().zip(back.choices()) {
            assert_eq!(a.address, b.address);
            assert_eq!(a.value, b.value);
        }
    }
}
