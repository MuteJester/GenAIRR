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
//! ## Phase B.1 scope
//!
//! Just the data structures. No pass integration yet (that's B.2).
//! No serialization — trace is in-memory only at this stage.

/// The typed value of a single random choice.
///
/// Intentionally a closed enum so contracts and replayers can match
/// exhaustively on the cases. New variants get added when new choice
/// kinds appear in the engine; the compiler refuses to build until
/// every consumer handles them.
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
    /// allele table. Concrete shape arrives in Phase C.
    AlleleId(u32),

    /// A boolean choice (e.g., D inversion: yes/no, receptor
    /// revision: yes/no, contaminant injection: yes/no).
    Bool(bool),
}

/// One entry in the trace: the address at which a choice was made
/// and the value that was sampled there.
#[derive(Clone, Debug, PartialEq, Eq)]
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
}

/// The append-only record of every choice made during one simulation.
///
/// Maintained by the runtime as a side-channel; not part of any
/// `Simulation` IR revision. Once a simulation completes the trace
/// is delivered to the caller (via the `History` API in B.2 or via
/// `with_history=True` in the public `.run()` from D9) and dropped
/// otherwise.
#[derive(Clone, Debug, Default)]
pub struct Trace {
    choices: Vec<ChoiceRecord>,
}

impl Trace {
    /// Empty trace.
    pub fn new() -> Self {
        Self::default()
    }

    /// Empty trace pre-allocated for an expected number of choices.
    /// Phase E will know the rough count from the active pass plan
    /// and call this to avoid reallocations.
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

    /// All recorded choices in chronological order.
    pub fn choices(&self) -> &[ChoiceRecord] {
        &self.choices
    }

    /// Find the first choice recorded at the given address, if any.
    /// Address comparison is exact-string equality.
    pub fn find(&self, address: &str) -> Option<&ChoiceRecord> {
        self.choices.iter().find(|c| c.address == address)
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
}
