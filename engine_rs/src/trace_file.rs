//! Durable trace artifact — the on-disk form of a simulation's
//! [`Trace`] paired with the metadata needed to interpret and
//! rerun it.
//!
//! # Schema versioning policy
//!
//! [`TraceFile`] is now a **public on-disk contract**, not a debug
//! dump. Once the engine has emitted a file at a given
//! `schema_version`, every future engine must continue to load and
//! correctly replay that file — or explicitly reject it with a
//! versioned error. Drift here translates directly into lost
//! reproducibility for downstream science artifacts.
//!
//! ## Versions in flight
//!
//! - **v1** (legacy): `schema_version`, `engine_version`, `seed`,
//!   `pass_plan_signature`, `refdata_signature` (structural only —
//!   chain type + per-pool counts), `trace`. No content hash. No
//!   address-vocabulary version field.
//! - **v2** (current): adds `address_schema_version` (the revision
//!   of the persisted [`crate::address::ChoiceAddress`] vocabulary)
//!   and `refdata_content_hash` (`"sha256:{hex}"` over the
//!   canonical chain+pools+alleles byte stream — see
//!   [`refdata_content_hash`]). Both fields are
//!   `#[serde(default)]` so v1 files still deserialise; the loader
//!   skips the corresponding validation gates when those fields
//!   are absent.
//!
//! [`KNOWN_TRACE_FILE_SCHEMA_VERSIONS`] enumerates every version
//! this engine accepts. [`TRACE_FILE_SCHEMA_VERSION`] is the
//! version every fresh emit carries.
//!
//! ## What triggers a `schema_version` bump
//!
//! Any change to the wire format that an old loader cannot
//! correctly handle. Concretely:
//!
//! - **New required field added.** v1 loaders cannot ignore it
//!   safely; they must reject.
//! - **Existing field semantics changed.** e.g. `refdata_signature`
//!   switched from structural to content-hashed (would be a
//!   semantic break, not just additive).
//! - **Trace record ordering rule changed.** The strict-positional
//!   cursor is part of the contract; reordering breaks replay.
//! - **A `ChoiceValue` variant gained / lost / renamed.** v1
//!   loaders can't deserialise the new tag.
//!
//! ## What does *not* trigger a bump
//!
//! - **New optional field added** (carries a serde `#[default]`,
//!   loadable as absent by older readers — but older readers won't
//!   gain whatever guarantee the new field provides).
//! - **New `ChoiceAddress` variant added** (covered by
//!   [`crate::address::ADDRESS_SCHEMA_VERSION`] instead — bump
//!   that, not `schema_version`).
//! - **A new built-in pass added** with new declared addresses
//!   (same — address-vocabulary version covers it).
//!
//! ## Why two version fields
//!
//! `schema_version` describes the **JSON wire shape**.
//! `address_schema_version` describes the **string spellings**
//! recorded inside `trace.choices[*].address`. They evolve
//! independently: adding a new sampling pass with a fresh address
//! bumps `ADDRESS_SCHEMA_VERSION` but leaves `schema_version`
//! alone; restructuring the top-level JSON bumps `schema_version`
//! without touching addresses. The two-axis design means the
//! engine can grow without forcing every fixture to be
//! regenerated.
//!
//! ## Signature design
//!
//! Three identity surfaces sit inside a v2 file, in increasing
//! strength:
//!
//! 1. **`pass_plan_signature`** — `"name1|name2|…"` of the pass
//!    names in plan order. Catches plan-changed regressions with
//!    a readable diff.
//! 2. **`refdata_signature`** — structural readable string. Cheap
//!    first-line identity check; insufficient on its own across
//!    machines.
//! 3. **`refdata_content_hash`** — SHA-256 over the canonical
//!    refdata byte stream. Two refdatas with matching structural
//!    signatures but different allele bytes produce different
//!    content hashes; this is the strong identity gate.
//!
//! [`Self::validate_against`] runs all three checks in order, so
//! the most informative error is reported.

use std::path::Path;

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::address::ChoiceAddress;
use crate::pass::PassPlan;
use crate::refdata::RefDataConfig;
use crate::trace::Trace;

/// Current on-disk schema version the engine emits. Bumped when a
/// new required field or wire-format change ships. See module docs
/// for the bump policy.
pub const TRACE_FILE_SCHEMA_VERSION: u32 = 2;

/// Schema versions this engine knows how to load. Loading a file
/// whose `schema_version` is outside this set surfaces as
/// [`TraceFileError::UnsupportedSchemaVersion`].
///
/// v1 → structural-only refdata signature, no address schema
/// version, no content hash.
/// v2 → adds `address_schema_version` + `refdata_content_hash`.
pub const KNOWN_TRACE_FILE_SCHEMA_VERSIONS: &[u32] = &[1, 2];

/// Durable artifact bundling a simulation's [`Trace`] with the
/// metadata needed to verify and rerun it. See the module docs.
///
/// **v1 → v2 evolution**: v2 adds two fields, both marked
/// `#[serde(default)]` so v1 files (which lack them) still
/// deserialize. The loader interprets a missing
/// `address_schema_version` as `0` and a missing
/// `refdata_content_hash` as `None`. Validation logic falls back
/// to the structural signature when no content hash is present;
/// when a hash IS present it must match the live refdata exactly.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TraceFile {
    pub schema_version: u32,
    pub engine_version: String,
    pub seed: u64,
    pub pass_plan_signature: String,
    pub refdata_signature: String,
    /// Address vocabulary revision the trace was recorded against.
    /// v1 traces predate this field and default to `0`; v2 traces
    /// always emit the current
    /// [`crate::address::ADDRESS_SCHEMA_VERSION`].
    #[serde(default)]
    pub address_schema_version: u32,
    /// Cryptographic identity of the refdata used to produce the
    /// trace. v1 traces predate this field; v2 traces always emit
    /// a `"sha256:{hex}"` string. See
    /// [`refdata_content_hash`].
    #[serde(default)]
    pub refdata_content_hash: Option<String>,
    pub trace: Trace,
}

/// Errors that can occur while loading, validating, or saving a
/// [`TraceFile`]. Distinct variants so callers can show a useful
/// diagnostic for each failure mode without string-parsing a
/// generic error.
#[derive(Debug)]
pub enum TraceFileError {
    /// On-disk schema version is not one this engine can read.
    UnsupportedSchemaVersion { got: u32, expected: u32 },
    /// The loaded plan does not match the trace's recorded plan.
    PassPlanSignatureMismatch { expected: String, got: String },
    /// The loaded refdata does not match the trace's recorded refdata.
    RefdataSignatureMismatch { expected: String, got: String },
    /// The loaded refdata's content hash doesn't match the trace's
    /// recorded hash. Stronger than `RefdataSignatureMismatch`:
    /// two refdatas with identical structural shape but different
    /// allele sequences would pass the structural check and fail
    /// here.
    RefdataContentHashMismatch { expected: String, got: String },
    /// The address vocabulary revision the trace was recorded
    /// against doesn't match what this engine speaks. A trace
    /// emitted under `address_schema_version = N` cannot be
    /// loaded by an engine that ships `M ≠ N` — the spelling /
    /// parse pair has shifted underneath it.
    AddressSchemaVersionMismatch { got: u32, expected: u32 },
    /// A recorded address doesn't parse to a built-in
    /// [`ChoiceAddress`] and the loader was asked to be strict.
    UnknownAddress { address: String },
    /// Filesystem I/O failure.
    Io(std::io::Error),
    /// JSON parse / serialise failure.
    Json(serde_json::Error),
}

impl std::fmt::Display for TraceFileError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedSchemaVersion { got, expected } => write!(
                f,
                "unsupported trace-file schema version: got {got}, this engine supports {expected}",
            ),
            Self::PassPlanSignatureMismatch { expected, got } => write!(
                f,
                "trace-file pass plan signature mismatch:\n  expected: {expected}\n  got:      {got}",
            ),
            Self::RefdataSignatureMismatch { expected, got } => write!(
                f,
                "trace-file refdata signature mismatch:\n  expected: {expected}\n  got:      {got}",
            ),
            Self::RefdataContentHashMismatch { expected, got } => write!(
                f,
                "trace-file refdata content-hash mismatch:\n  expected: {expected}\n  got:      {got}",
            ),
            Self::AddressSchemaVersionMismatch { got, expected } => write!(
                f,
                "trace-file address schema version mismatch: got {got}, this engine speaks {expected}",
            ),
            Self::UnknownAddress { address } => {
                write!(f, "trace-file contains non-built-in address {address:?}; pass allow_custom_addresses=true to accept it")
            }
            Self::Io(e) => write!(f, "trace-file I/O error: {e}"),
            Self::Json(e) => write!(f, "trace-file JSON error: {e}"),
        }
    }
}

impl std::error::Error for TraceFileError {}

impl From<std::io::Error> for TraceFileError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<serde_json::Error> for TraceFileError {
    fn from(e: serde_json::Error) -> Self {
        Self::Json(e)
    }
}

/// Canonical signature of a [`PassPlan`]. Plan order matters; the
/// signature is the pass names in insertion order joined by `|`.
pub fn pass_plan_signature(plan: &PassPlan) -> String {
    plan.passes()
        .iter()
        .map(|p| p.name())
        .collect::<Vec<_>>()
        .join("|")
}

/// Structural signature of a [`RefDataConfig`] — readable
/// chain-type + per-pool allele counts. Useful for debug output
/// and as a fast first-line identity check. Not sufficient for
/// cross-machine reproducibility (two refdatas with identical
/// counts but different allele sequences hash to the same
/// structural signature); pair it with [`refdata_content_hash`]
/// for content-level identity.
pub fn refdata_signature(refdata: &RefDataConfig) -> String {
    format!(
        "chain={:?};v={};d={};j={};c={}",
        refdata.chain_type,
        refdata.v_pool.len(),
        refdata.d_pool.len(),
        refdata.j_pool.len(),
        refdata.c_pool.len(),
    )
}

/// Deterministic content hash of a [`RefDataConfig`]. Two
/// `RefDataConfig`s produce the same hash iff they carry the same
/// chain type, the same pools in the same order, and the same
/// per-allele `(name, gene, segment, seq, anchor)` tuples.
///
/// The hash is a hex-encoded SHA-256 over a canonical byte stream:
///
/// ```text
/// chain={chain_type}\n
/// v_pool:{count}\n
///   v[{i}]:name={name}|gene={gene}|seg={segment}|seq={ascii}|anchor={anchor_or_NA}\n
///   …
/// d_pool:{count}\n
///   …
/// j_pool:{count}\n
///   …
/// c_pool:{count}\n
///   …
/// ```
///
/// **Ordering**: alleles are emitted in `AllelePool::iter` order,
/// which is the order they were `push`'d. The Python loader
/// constructs refdata in a deterministic order from the bundled
/// `DataConfig`, so a given config name always produces the same
/// hash. Two engine builds with different bundled datasets produce
/// different hashes — exactly what cross-machine identity needs.
///
/// **Hash strength**: SHA-256 is cryptographic, well beyond what
/// is needed to distinguish refdatas. The choice is about
/// determinism and freedom from collision concerns under any
/// realistic load. The on-disk representation is `"sha256:{hex}"`
/// so future hash algorithms can be substituted with an explicit
/// prefix change (and a corresponding schema bump).
pub fn refdata_content_hash(refdata: &RefDataConfig) -> String {
    let mut hasher = Sha256::new();
    let write_field = |buf: &mut Sha256, label: &str, value: &[u8]| {
        buf.update(label.as_bytes());
        buf.update(b"=");
        buf.update(value);
    };

    hasher.update(format!("chain={:?}\n", refdata.chain_type).as_bytes());

    for (label, pool) in [
        ("v_pool", &refdata.v_pool),
        ("d_pool", &refdata.d_pool),
        ("j_pool", &refdata.j_pool),
        ("c_pool", &refdata.c_pool),
    ] {
        hasher.update(format!("{label}:{}\n", pool.len()).as_bytes());
        for (id, allele) in pool.iter() {
            hasher.update(format!("  {label}[{}]:", id.index()).as_bytes());
            write_field(&mut hasher, "name", allele.name.as_bytes());
            hasher.update(b"|");
            write_field(&mut hasher, "gene", allele.gene.as_bytes());
            hasher.update(b"|");
            write_field(
                &mut hasher,
                "seg",
                format!("{:?}", allele.segment).as_bytes(),
            );
            hasher.update(b"|");
            write_field(&mut hasher, "seq", &allele.seq);
            hasher.update(b"|");
            let anchor_str = allele
                .anchor
                .map(|a| a.to_string())
                .unwrap_or_else(|| "NA".to_string());
            write_field(&mut hasher, "anchor", anchor_str.as_bytes());
            hasher.update(b"\n");
        }
    }

    let digest = hasher.finalize();
    format!("sha256:{:x}", digest)
}

impl TraceFile {
    /// Build a [`TraceFile`] from a fresh `(plan, refdata, seed, trace)`
    /// tuple. Stamps the current schema version, engine version,
    /// address schema version, and refdata content hash onto the
    /// file — i.e. emits v2.
    pub fn build(plan: &PassPlan, refdata: &RefDataConfig, seed: u64, trace: Trace) -> Self {
        Self {
            schema_version: TRACE_FILE_SCHEMA_VERSION,
            engine_version: env!("CARGO_PKG_VERSION").to_string(),
            seed,
            pass_plan_signature: pass_plan_signature(plan),
            refdata_signature: refdata_signature(refdata),
            address_schema_version: crate::address::ADDRESS_SCHEMA_VERSION,
            refdata_content_hash: Some(refdata_content_hash(refdata)),
            trace,
        }
    }

    /// Deserialise from a JSON string. Accepts any version listed
    /// in [`KNOWN_TRACE_FILE_SCHEMA_VERSIONS`]; rejects others.
    /// Does not check pass/refdata signatures — that's
    /// [`Self::validate_against`].
    pub fn from_json(s: &str) -> Result<Self, TraceFileError> {
        let f: TraceFile = serde_json::from_str(s)?;
        if !KNOWN_TRACE_FILE_SCHEMA_VERSIONS.contains(&f.schema_version) {
            return Err(TraceFileError::UnsupportedSchemaVersion {
                got: f.schema_version,
                expected: TRACE_FILE_SCHEMA_VERSION,
            });
        }
        Ok(f)
    }

    /// Serialise to a pretty-printed JSON string.
    pub fn to_json_pretty(&self) -> Result<String, TraceFileError> {
        Ok(serde_json::to_string_pretty(self)?)
    }

    /// Load a trace file from disk. Equivalent to reading the file
    /// then calling [`Self::from_json`].
    pub fn read_from(path: impl AsRef<Path>) -> Result<Self, TraceFileError> {
        let s = std::fs::read_to_string(path.as_ref())?;
        Self::from_json(&s)
    }

    /// Write a trace file to disk as pretty-printed JSON.
    pub fn write_to(&self, path: impl AsRef<Path>) -> Result<(), TraceFileError> {
        let s = self.to_json_pretty()?;
        std::fs::write(path.as_ref(), s)?;
        Ok(())
    }

    /// Verify that this trace file's plan and refdata signatures
    /// match the currently-loaded plan and refdata. Returns
    /// `Ok(())` on match, otherwise the first mismatch encountered.
    ///
    /// Checks, in order:
    ///   1. Pass-plan signature (always).
    ///   2. Refdata structural signature (always).
    ///   3. Refdata content hash — when the trace file carries one
    ///      (v2 and later). v1 traces predate the hash and skip
    ///      this check, falling back to the structural signature
    ///      as the only refdata-identity gate.
    ///   4. Address-vocabulary schema version — when the trace
    ///      file carries one (v2). A trace recorded against a
    ///      different address-spelling generation can't safely
    ///      replay through this engine's cursor.
    pub fn validate_against(
        &self,
        plan: &PassPlan,
        refdata: &RefDataConfig,
    ) -> Result<(), TraceFileError> {
        let live_plan_sig = pass_plan_signature(plan);
        if live_plan_sig != self.pass_plan_signature {
            return Err(TraceFileError::PassPlanSignatureMismatch {
                expected: self.pass_plan_signature.clone(),
                got: live_plan_sig,
            });
        }
        let live_refdata_sig = refdata_signature(refdata);
        if live_refdata_sig != self.refdata_signature {
            return Err(TraceFileError::RefdataSignatureMismatch {
                expected: self.refdata_signature.clone(),
                got: live_refdata_sig,
            });
        }
        if let Some(ref recorded_hash) = self.refdata_content_hash {
            let live_hash = refdata_content_hash(refdata);
            if *recorded_hash != live_hash {
                return Err(TraceFileError::RefdataContentHashMismatch {
                    expected: recorded_hash.clone(),
                    got: live_hash,
                });
            }
        }
        // v1 traces report `address_schema_version = 0` (the
        // serde `#[default]`); accept those without a version
        // gate. v2+ traces emit the live constant and must match.
        if self.address_schema_version != 0
            && self.address_schema_version != crate::address::ADDRESS_SCHEMA_VERSION
        {
            return Err(TraceFileError::AddressSchemaVersionMismatch {
                got: self.address_schema_version,
                expected: crate::address::ADDRESS_SCHEMA_VERSION,
            });
        }
        Ok(())
    }

    /// Walk every recorded address and confirm it parses to a
    /// built-in [`ChoiceAddress`]. Returns the first non-parsing
    /// address found, if any. `allow_custom_addresses=true` skips
    /// the check (the trace may carry addresses from user-defined
    /// passes).
    pub fn validate_addresses(&self, allow_custom_addresses: bool) -> Result<(), TraceFileError> {
        if allow_custom_addresses {
            return Ok(());
        }
        for rec in self.trace.choices() {
            if ChoiceAddress::parse(&rec.address).is_none() {
                return Err(TraceFileError::UnknownAddress {
                    address: rec.address.clone(),
                });
            }
        }
        Ok(())
    }
}

// ──────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pass::{PassPlan, Pass, PassContext, PassEffect};
    use crate::refdata::{Allele, ChainType, RefDataConfig};
    use crate::ir::{Segment, Simulation};
    use crate::address;
    use crate::trace::{ChoiceRecord, ChoiceValue};

    /// Minimal dummy pass for plan-signature tests. Does nothing on
    /// execute; its `name()` is the only thing the signature reads.
    struct NamedPass(&'static str);
    impl Pass for NamedPass {
        fn name(&self) -> &str {
            self.0
        }
        fn execute(&self, sim: &Simulation, _ctx: &mut PassContext) -> Simulation {
            sim.clone()
        }
        fn declared_choice_patterns(&self) -> Vec<address::ChoiceAddressPattern> {
            vec![]
        }
        fn effects(&self) -> Vec<PassEffect> {
            vec![]
        }
    }

    fn make_plan(names: &[&'static str]) -> PassPlan {
        let mut plan = PassPlan::new();
        for n in names {
            plan.push(Box::new(NamedPass(n)));
        }
        plan
    }

    fn make_refdata() -> RefDataConfig {
        let mut cfg = RefDataConfig::empty(ChainType::Vj);
        let _ = cfg.v_pool.push(Allele {
            name: "v_a*01".into(),
            gene: "v_a".into(),
            seq: b"ACGT".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = cfg.j_pool.push(Allele {
            name: "j_a*01".into(),
            gene: "j_a".into(),
            seq: b"TGGA".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });
        cfg
    }

    fn make_trace() -> Trace {
        let mut t = Trace::new();
        t.record("sample_allele.v", ChoiceValue::AlleleId(0));
        t.record("sample_allele.j", ChoiceValue::AlleleId(0));
        t.record("np.np1.length", ChoiceValue::Int(3));
        t.record("np.np1.bases[0]", ChoiceValue::Base(b'A'));
        t.record("np.np1.bases[1]", ChoiceValue::Base(b'C'));
        t.record("np.np1.bases[2]", ChoiceValue::Base(b'G'));
        t
    }

    #[test]
    fn pass_plan_signature_is_pipe_joined_names_in_plan_order() {
        let plan = make_plan(&["sample_allele.v", "sample_allele.j", "generate_np.np1"]);
        assert_eq!(
            pass_plan_signature(&plan),
            "sample_allele.v|sample_allele.j|generate_np.np1",
        );
    }

    #[test]
    fn refdata_signature_is_structural() {
        let cfg = make_refdata();
        let sig = refdata_signature(&cfg);
        assert_eq!(sig, "chain=Vj;v=1;d=0;j=1;c=0");
    }

    #[test]
    fn trace_file_round_trips_through_json() {
        let plan = make_plan(&["sample_allele.v", "sample_allele.j"]);
        let refdata = make_refdata();
        let trace = make_trace();
        let original = TraceFile::build(&plan, &refdata, 0xdead_beef, trace);

        let s = original.to_json_pretty().unwrap();
        let back = TraceFile::from_json(&s).unwrap();

        assert_eq!(back.schema_version, original.schema_version);
        assert_eq!(back.engine_version, original.engine_version);
        assert_eq!(back.seed, original.seed);
        assert_eq!(back.pass_plan_signature, original.pass_plan_signature);
        assert_eq!(back.refdata_signature, original.refdata_signature);
        assert_eq!(back.trace.len(), original.trace.len());
        for (a, b) in original.trace.choices().iter().zip(back.trace.choices()) {
            assert_eq!(a.address, b.address);
            assert_eq!(a.value, b.value);
        }
    }

    #[test]
    fn validate_against_same_plan_and_refdata_succeeds() {
        let plan = make_plan(&["a", "b"]);
        let refdata = make_refdata();
        let file = TraceFile::build(&plan, &refdata, 0, Trace::new());
        assert!(file.validate_against(&plan, &refdata).is_ok());
    }

    #[test]
    fn validate_against_mismatched_plan_returns_mismatch_error() {
        let plan_a = make_plan(&["a", "b"]);
        let refdata = make_refdata();
        let file = TraceFile::build(&plan_a, &refdata, 0, Trace::new());

        let plan_b = make_plan(&["a", "c"]);
        let err = file.validate_against(&plan_b, &refdata).unwrap_err();
        match err {
            TraceFileError::PassPlanSignatureMismatch { expected, got } => {
                assert_eq!(expected, "a|b");
                assert_eq!(got, "a|c");
            }
            other => panic!("expected PassPlanSignatureMismatch, got {other:?}"),
        }
    }

    #[test]
    fn validate_against_mismatched_refdata_returns_mismatch_error() {
        let plan = make_plan(&["a"]);
        let refdata_a = make_refdata();
        let file = TraceFile::build(&plan, &refdata_a, 0, Trace::new());

        let mut refdata_b = make_refdata();
        let _ = refdata_b.v_pool.push(Allele {
            name: "v_b*01".into(),
            gene: "v_b".into(),
            seq: b"AAAA".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });

        let err = file.validate_against(&plan, &refdata_b).unwrap_err();
        match err {
            TraceFileError::RefdataSignatureMismatch { expected, got } => {
                assert_eq!(expected, "chain=Vj;v=1;d=0;j=1;c=0");
                assert_eq!(got, "chain=Vj;v=2;d=0;j=1;c=0");
            }
            other => panic!("expected RefdataSignatureMismatch, got {other:?}"),
        }
    }

    #[test]
    fn from_json_rejects_future_schema_version() {
        // Hand-craft a JSON document with a forward schema version.
        let raw = serde_json::json!({
            "schema_version": 999,
            "engine_version": "0.0.0",
            "seed": 0,
            "pass_plan_signature": "",
            "refdata_signature": "",
            "trace": { "choices": [] },
        });
        let err = TraceFile::from_json(&raw.to_string()).unwrap_err();
        match err {
            TraceFileError::UnsupportedSchemaVersion { got, expected } => {
                assert_eq!(got, 999);
                assert_eq!(expected, TRACE_FILE_SCHEMA_VERSION);
            }
            other => panic!("expected UnsupportedSchemaVersion, got {other:?}"),
        }
    }

    #[test]
    fn validate_addresses_strict_rejects_unknown_address() {
        let plan = make_plan(&["a"]);
        let refdata = make_refdata();
        let mut trace = Trace::new();
        trace.record("user.custom.pass", ChoiceValue::Int(1));
        let file = TraceFile::build(&plan, &refdata, 0, trace);

        let err = file.validate_addresses(false).unwrap_err();
        match err {
            TraceFileError::UnknownAddress { address } => {
                assert_eq!(address, "user.custom.pass");
            }
            other => panic!("expected UnknownAddress, got {other:?}"),
        }
    }

    #[test]
    fn validate_addresses_permissive_accepts_unknown_address() {
        let plan = make_plan(&["a"]);
        let refdata = make_refdata();
        let mut trace = Trace::new();
        trace.record("user.custom.pass", ChoiceValue::Int(1));
        let file = TraceFile::build(&plan, &refdata, 0, trace);

        assert!(file.validate_addresses(true).is_ok());
    }

    #[test]
    fn validate_addresses_accepts_every_built_in_address() {
        let plan = make_plan(&["a"]);
        let refdata = make_refdata();
        let trace = make_trace();
        let file = TraceFile::build(&plan, &refdata, 0, trace);
        assert!(file.validate_addresses(false).is_ok());
    }

    #[test]
    fn write_and_read_round_trip_through_disk() {
        let plan = make_plan(&["sample_allele.v", "sample_allele.j"]);
        let refdata = make_refdata();
        let trace = make_trace();
        let original = TraceFile::build(&plan, &refdata, 42, trace);

        let tmp = std::env::temp_dir().join("genairr_trace_file_disk_round_trip.json");
        original.write_to(&tmp).unwrap();
        let back = TraceFile::read_from(&tmp).unwrap();
        std::fs::remove_file(&tmp).ok();

        assert_eq!(back.seed, original.seed);
        assert_eq!(back.trace.len(), original.trace.len());
        for (a, b) in original.trace.choices().iter().zip(back.trace.choices()) {
            assert_eq!(a.address, b.address);
            assert_eq!(a.value, b.value);
        }
    }

    #[test]
    fn unused_choice_record_field_compiles_through_re_export() {
        // Pin that ChoiceRecord remains the public companion of Trace.
        // (Compile-fence — the type still being a stable export is what
        // matters for downstream callers building TraceFiles by hand.)
        let _ = ChoiceRecord::new("trim.v_3", ChoiceValue::Int(0));
    }

    // ── Schema policy compatibility tests ─────────────────────────

    #[test]
    fn schema_version_constant_is_2() {
        // Pinned so a downstream consumer adding a new optional
        // field can't accidentally bump the version. Bumping must
        // come with explicit documentation + KNOWN list update.
        assert_eq!(TRACE_FILE_SCHEMA_VERSION, 2);
        assert!(KNOWN_TRACE_FILE_SCHEMA_VERSIONS.contains(&1));
        assert!(KNOWN_TRACE_FILE_SCHEMA_VERSIONS.contains(&2));
    }

    #[test]
    fn fresh_emit_carries_v2_metadata() {
        let plan = make_plan(&["sample_allele.v"]);
        let refdata = make_refdata();
        let file = TraceFile::build(&plan, &refdata, 0, Trace::new());

        assert_eq!(file.schema_version, 2);
        assert_eq!(
            file.address_schema_version,
            crate::address::ADDRESS_SCHEMA_VERSION
        );
        let hash = file.refdata_content_hash.as_ref().expect("v2 emits hash");
        assert!(hash.starts_with("sha256:"));
        // Hash is 64 hex chars after the prefix.
        assert_eq!(hash.len(), "sha256:".len() + 64);
    }

    #[test]
    fn loader_accepts_v1_trace_with_missing_new_fields() {
        // Hand-craft a v1-shape JSON: no address_schema_version,
        // no refdata_content_hash. Loading via from_json must
        // succeed; the missing fields take their serde defaults.
        let raw = serde_json::json!({
            "schema_version": 1,
            "engine_version": "0.0.1",
            "seed": 7,
            "pass_plan_signature": "a|b",
            "refdata_signature": "chain=Vj;v=1;d=0;j=1;c=0",
            "trace": { "choices": [] },
        });
        let f = TraceFile::from_json(&raw.to_string()).unwrap();
        assert_eq!(f.schema_version, 1);
        assert_eq!(f.address_schema_version, 0);
        assert!(f.refdata_content_hash.is_none());
    }

    #[test]
    fn v1_trace_validates_via_structural_signature_only() {
        // v1 file → validation falls back to structural signature.
        // The content-hash check is skipped because the file
        // doesn't have one.
        let plan = make_plan(&["a"]);
        let refdata = make_refdata();
        // Manually construct a v1-shape file pointing at the same
        // refdata's structural signature.
        let v1 = TraceFile {
            schema_version: 1,
            engine_version: "0.0.1".into(),
            seed: 0,
            pass_plan_signature: pass_plan_signature(&plan),
            refdata_signature: refdata_signature(&refdata),
            address_schema_version: 0,
            refdata_content_hash: None,
            trace: Trace::new(),
        };
        assert!(v1.validate_against(&plan, &refdata).is_ok());
    }

    #[test]
    fn v2_content_hash_match_succeeds() {
        let plan = make_plan(&["a"]);
        let refdata = make_refdata();
        let file = TraceFile::build(&plan, &refdata, 0, Trace::new());
        assert!(file.validate_against(&plan, &refdata).is_ok());
    }

    #[test]
    fn v2_content_hash_mismatch_returns_structured_error() {
        // Construct a v2 file with a manually-bogus hash so the
        // structural signature matches but the content hash does
        // not — proves the content check fires AFTER the
        // structural check and surfaces its own error variant.
        let plan = make_plan(&["a"]);
        let refdata = make_refdata();
        let mut file = TraceFile::build(&plan, &refdata, 0, Trace::new());
        file.refdata_content_hash = Some("sha256:deadbeef".into());

        let err = file.validate_against(&plan, &refdata).unwrap_err();
        match err {
            TraceFileError::RefdataContentHashMismatch { expected, .. } => {
                assert_eq!(expected, "sha256:deadbeef");
            }
            other => panic!("expected RefdataContentHashMismatch, got {other:?}"),
        }
    }

    #[test]
    fn refdata_content_hash_is_stable_under_repeated_calls() {
        let refdata = make_refdata();
        let a = refdata_content_hash(&refdata);
        let b = refdata_content_hash(&refdata);
        assert_eq!(a, b);
    }

    #[test]
    fn refdata_content_hash_distinguishes_seq_only_differences() {
        // Two refdatas with identical structural signature
        // (same chain, same counts) but different allele bytes
        // must hash differently. This is the core property that
        // makes the content hash worth carrying over the
        // structural signature.
        let mut a = RefDataConfig::empty(ChainType::Vj);
        let _ = a.v_pool.push(Allele {
            name: "v*01".into(),
            gene: "v".into(),
            seq: b"AAAA".to_vec(),
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = a.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TTTT".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        let mut b = RefDataConfig::empty(ChainType::Vj);
        let _ = b.v_pool.push(Allele {
            name: "v*01".into(),
            gene: "v".into(),
            seq: b"CCCC".to_vec(), // ← only difference
            segment: Segment::V,
            anchor: Some(0),
        });
        let _ = b.j_pool.push(Allele {
            name: "j*01".into(),
            gene: "j".into(),
            seq: b"TTTT".to_vec(),
            segment: Segment::J,
            anchor: Some(0),
        });

        // Structural signatures match → would pass v1 check.
        assert_eq!(refdata_signature(&a), refdata_signature(&b));
        // Content hashes differ → v2 catches the divergence.
        assert_ne!(refdata_content_hash(&a), refdata_content_hash(&b));
    }

    #[test]
    fn address_schema_version_mismatch_rejects_at_validate() {
        // Hand-craft a file with a forward address_schema_version.
        let plan = make_plan(&["a"]);
        let refdata = make_refdata();
        let mut file = TraceFile::build(&plan, &refdata, 0, Trace::new());
        file.address_schema_version = 99;

        let err = file.validate_against(&plan, &refdata).unwrap_err();
        match err {
            TraceFileError::AddressSchemaVersionMismatch { got, expected } => {
                assert_eq!(got, 99);
                assert_eq!(expected, crate::address::ADDRESS_SCHEMA_VERSION);
            }
            other => panic!("expected AddressSchemaVersionMismatch, got {other:?}"),
        }
    }
}
