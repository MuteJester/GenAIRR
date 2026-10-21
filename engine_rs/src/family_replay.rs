//! Versioned native foundations for self-contained family choice replay.
//! Only audited constructors produce executable state; this module never imports IR.
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::cell::RefCell;
use std::collections::HashSet;
use std::sync::{Arc, Mutex};

pub mod plan;
pub mod reference;
mod state;
pub mod strict_json;
pub use state::{outcome_events_digest, semantic_state_digest};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct F64Bits(pub u64);
impl F64Bits {
    pub fn new(value: f64) -> Result<Self, String> {
        if !value.is_finite() {
            return Err("nonfinite binary64 value".into());
        }
        Ok(Self(value.to_bits()))
    }
    pub fn from_f64(value: f64) -> Result<Self, String> {
        Self::new(value)
    }
    pub fn get(self) -> f64 {
        f64::from_bits(self.0)
    }
}
impl Serialize for F64Bits {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        if !self.get().is_finite() {
            return Err(serde::ser::Error::custom("nonfinite binary64"));
        }
        serializer.serialize_str(&format!("{:016x}", self.0))
    }
}
impl<'de> Deserialize<'de> for F64Bits {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let value = String::deserialize(deserializer)?;
        if value.len() != 16
            || !value
                .bytes()
                .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
        {
            return Err(serde::de::Error::custom(
                "binary64 requires 16 lowercase hex digits",
            ));
        }
        let bits = u64::from_str_radix(&value, 16).map_err(serde::de::Error::custom)?;
        if !f64::from_bits(bits).is_finite() {
            return Err(serde::de::Error::custom("nonfinite binary64"));
        }
        Ok(Self(bits))
    }
}

/// Length framed and type tagged; strings never alias labels or adjacent fields.
pub fn canonical_digest(domain: &str, value: &Value) -> Result<String, String> {
    fn bytes(hash: &mut Sha256, data: &[u8]) {
        hash.update((data.len() as u64).to_be_bytes());
        hash.update(data);
    }
    fn visit(hash: &mut Sha256, value: &Value) -> Result<(), String> {
        match value {
            Value::Null => hash.update(b"n"),
            Value::Bool(v) => hash.update(if *v { b"t" } else { b"f" }),
            Value::String(v) => {
                hash.update(b"s");
                bytes(hash, v.as_bytes());
            }
            Value::Number(v) => {
                if let Some(n) = v.as_i64() {
                    hash.update(b"i");
                    bytes(hash, n.to_string().as_bytes());
                } else if let Some(n) = v.as_u64() {
                    hash.update(b"i");
                    bytes(hash, n.to_string().as_bytes());
                } else {
                    let n = v.as_f64().ok_or("invalid JSON number")?;
                    if !n.is_finite() {
                        return Err("nonfinite canonical number".into());
                    }
                    hash.update(b"d");
                    hash.update(n.to_bits().to_be_bytes());
                }
            }
            Value::Array(values) => {
                hash.update(b"a");
                hash.update((values.len() as u64).to_be_bytes());
                for value in values {
                    visit(hash, value)?;
                }
            }
            Value::Object(values) => {
                hash.update(b"o");
                hash.update((values.len() as u64).to_be_bytes());
                let mut entries: Vec<_> = values.iter().collect();
                entries.sort_unstable_by(|a, b| a.0.as_bytes().cmp(b.0.as_bytes()));
                for (key, value) in entries {
                    bytes(hash, key.as_bytes());
                    visit(hash, value)?;
                }
            }
        }
        Ok(())
    }
    let mut hash = Sha256::new();
    bytes(&mut hash, b"genairr_typed_canonical_v1");
    bytes(&mut hash, domain.as_bytes());
    visit(&mut hash, value)?;
    Ok(format!("sha256:{:x}", hash.finalize()))
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct FamilyReplayLimits {
    pub max_file_bytes: u64,
    pub max_depth: u64,
    pub max_nodes: u64,
    pub max_choices: u64,
    pub max_bases: u64,
    pub max_work: u64,
    pub max_attempts: u64,
    pub max_metadata_bytes: u64,
    pub max_string_bytes: u64,
    pub max_components: u64,
    pub max_values: u64,
    pub max_history_entries: u64,
    pub max_retained_bytes: u64,
}
impl Default for FamilyReplayLimits {
    fn default() -> Self {
        Self {
            max_file_bytes: 256 * 1024 * 1024,
            max_depth: 64,
            max_nodes: 1_000_000,
            max_choices: 5_000_000,
            max_bases: 50_000_000,
            max_work: 10_000_000,
            max_attempts: 10_000,
            max_metadata_bytes: 16 * 1024 * 1024,
            max_string_bytes: 16 * 1024 * 1024,
            max_components: 100_000,
            max_values: 5_000_000,
            max_history_entries: 5_000_000,
            max_retained_bytes: 256 * 1024 * 1024,
        }
    }
}
#[derive(Debug, Default, Serialize)]
struct Counters {
    work: u64,
    bases: u64,
    nodes: u64,
    choices: u64,
    attempts: u64,
    history_entries: u64,
    retained_bytes: u64,
    #[serde(skip)]
    founders: HashSet<usize>,
}
#[derive(Clone, Debug)]
pub struct FamilyReplayBudget {
    pub limits: FamilyReplayLimits,
    counters: Arc<Mutex<Counters>>,
}
impl Default for FamilyReplayBudget {
    fn default() -> Self {
        Self::new(FamilyReplayLimits::default())
    }
}
impl FamilyReplayBudget {
    pub fn new(limits: FamilyReplayLimits) -> Self {
        Self {
            limits,
            counters: Arc::new(Mutex::new(Counters::default())),
        }
    }
    fn increment(value: &mut u64, amount: u64, limit: u64, label: &str) -> Result<(), String> {
        let next = value
            .checked_add(amount)
            .ok_or_else(|| format!("resource limit: {label} overflow"))?;
        if next > limit {
            return Err(format!("resource limit: {label}: {next} > {limit}"));
        }
        *value = next;
        Ok(())
    }
    pub fn charge_work(&mut self, n: u64) -> Result<(), String> {
        Self::increment(
            &mut self
                .counters
                .lock()
                .map_err(|_| "budget lock poisoned")?
                .work,
            n,
            self.limits.max_work,
            "work",
        )
    }
    pub fn charge_bases(&mut self, n: u64) -> Result<(), String> {
        Self::increment(
            &mut self
                .counters
                .lock()
                .map_err(|_| "budget lock poisoned")?
                .bases,
            n,
            self.limits.max_bases,
            "bases",
        )
    }
    pub fn charge_nodes(&mut self, n: u64) -> Result<(), String> {
        Self::increment(
            &mut self
                .counters
                .lock()
                .map_err(|_| "budget lock poisoned")?
                .nodes,
            n,
            self.limits.max_nodes,
            "nodes",
        )
    }
    pub fn charge_choices(&mut self, n: u64) -> Result<(), String> {
        Self::increment(
            &mut self
                .counters
                .lock()
                .map_err(|_| "budget lock poisoned")?
                .choices,
            n,
            self.limits.max_choices,
            "choices",
        )
    }
    /// Charge non-provenance trace/event containers before an outcome copy.
    /// Provenance vectors are charged by their actual Clone implementation.
    pub fn charge_history_copy(&mut self, entries: u64, bytes: u64) -> Result<(), String> {
        let mut counters = self.counters.lock().map_err(|_| "budget lock poisoned")?;
        Self::increment(
            &mut counters.history_entries,
            entries,
            self.limits.max_history_entries,
            "history_entries",
        )?;
        Self::increment(
            &mut counters.retained_bytes,
            bytes,
            self.limits.max_retained_bytes,
            "retained_bytes",
        )
    }
    pub fn charge_attempts(&mut self, n: u64) -> Result<(), String> {
        Self::increment(
            &mut self
                .counters
                .lock()
                .map_err(|_| "budget lock poisoned")?
                .attempts,
            n,
            self.limits.max_attempts,
            "attempts",
        )
    }
    pub fn charge_history(&mut self, sim: &crate::ir::Simulation) -> Result<(), String> {
        if let Some(p) = &sim.duplex_provenance {
            self.charge_provenance(p)?;
        }
        Ok(())
    }
    pub fn charge_provenance(&mut self, p: &crate::ir::DuplexProvenance) -> Result<(), String> {
        let mut count = self.counters.lock().map_err(|_| "budget lock poisoned")?;
        let entries = p
            .links
            .capacity()
            .checked_add(p.removed.capacity())
            .and_then(|n| n.checked_add(p.edits.capacity()))
            .and_then(|n| n.checked_add(p.attempt_stages.capacity()))
            .ok_or("history capacity overflow")? as u64;
        Self::increment(
            &mut count.history_entries,
            entries,
            self.limits.max_history_entries,
            "history_entries",
        )?;
        let bytes = (p.links.capacity() as u64)
            .checked_mul(std::mem::size_of::<crate::ir::DuplexQueryLink>() as u64)
            .and_then(|n| {
                n.checked_add(
                    p.removed.capacity() as u64
                        * std::mem::size_of::<crate::ir::DuplexRemovedLink>() as u64,
                )
            })
            .and_then(|n| {
                n.checked_add(
                    p.edits.capacity() as u64 * std::mem::size_of::<crate::ir::DuplexEdit>() as u64,
                )
            })
            .and_then(|n| {
                n.checked_add(
                    p.attempt_stages.capacity() as u64
                        * std::mem::size_of::<crate::ir::DuplexAttemptSummary>() as u64,
                )
            })
            .ok_or("history byte accounting overflow")?;
        let strings = p
            .edits
            .iter()
            .try_fold(0u64, |sum, edit| {
                let segment = match &edit.action {
                    crate::ir::DuplexEditAction::Substitution { segment, .. }
                    | crate::ir::DuplexEditAction::Insertion { segment, .. }
                    | crate::ir::DuplexEditAction::Deletion { segment, .. } => segment.capacity(),
                    _ => 0,
                };
                sum.checked_add(edit.source_pass.as_ref().map_or(0, |s| s.capacity()) as u64)
                    .and_then(|value| value.checked_add(segment as u64))
            })
            .ok_or("history strings overflow")?;
        Self::increment(
            &mut count.retained_bytes,
            bytes
                .checked_add(strings)
                .ok_or("history byte accounting overflow")?,
            self.limits.max_retained_bytes,
            "retained_bytes",
        )?;
        let key = Arc::as_ptr(&p.founder) as usize;
        if !count.founders.contains(&key) {
            let data=serde_json::to_vec(&json!({"annotated":p.founder.annotated.as_ref(),"bindings":p.founder.bindings.as_ref(),"body_claims":p.founder.body_claims,"p_metrics":p.founder.p_metrics,"selection":p.founder.selection})).map_err(|e|e.to_string())?;
            Self::increment(
                &mut count.retained_bytes,
                data.len() as u64,
                self.limits.max_retained_bytes,
                "retained_bytes",
            )?;
            count.founders.insert(key);
        }
        Ok(())
    }
    pub fn snapshot(&self) -> Result<String, String> {
        serde_json::to_string(&*self.counters.lock().map_err(|_| "budget lock poisoned")?)
            .map_err(|e| e.to_string())
    }
}

// Primitive mutation/RNG APIs are infallible. A private unwind aborts ONLY an
// opted-in bounded operation; its boundary catches this payload and returns an
// ordinary resource error. Other panics retain their normal behavior.
#[derive(Debug)]
struct LimitExceeded(String);
thread_local! { static ACTIVE:RefCell<Vec<FamilyReplayBudget>>=const{RefCell::new(Vec::new())}; }
thread_local! { static REPLAY:std::cell::Cell<u32>=const{std::cell::Cell::new(0)}; }
pub(crate) fn is_budgeted() -> bool {
    ACTIVE.with(|stack| !stack.borrow().is_empty())
}
pub fn with_budget_scope<T>(
    budget: &mut FamilyReplayBudget,
    f: impl FnOnce(&mut FamilyReplayBudget) -> Result<T, String>,
) -> Result<T, String> {
    ACTIVE.with(|s| s.borrow_mut().push(budget.clone()));
    struct Pop;
    impl Drop for Pop {
        fn drop(&mut self) {
            ACTIVE.with(|s| {
                s.borrow_mut().pop();
            });
        }
    }
    let _pop = Pop;
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| f(budget))) {
        Ok(v) => v,
        Err(payload) => match payload.downcast::<LimitExceeded>() {
            Ok(v) => Err(v.0),
            Err(other) => std::panic::resume_unwind(other),
        },
    }
}
pub(crate) fn primitive_work(n: u64) {
    ACTIVE.with(|s| {
        if let Some(b) = s.borrow_mut().last_mut() {
            if let Err(e) = b.charge_work(n) {
                std::panic::resume_unwind(Box::new(LimitExceeded(e)));
            }
        }
    });
}
pub fn with_replay_scope<T>(f: impl FnOnce() -> Result<T, String>) -> Result<T, String> {
    REPLAY.with(|depth| depth.set(depth.get() + 1));
    struct Pop;
    impl Drop for Pop {
        fn drop(&mut self) {
            REPLAY.with(|depth| depth.set(depth.get() - 1));
        }
    }
    let _pop = Pop;
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)) {
        Ok(v) => v,
        Err(payload) => match payload.downcast::<LimitExceeded>() {
            Ok(v) => Err(v.0),
            Err(other) => std::panic::resume_unwind(other),
        },
    }
}
pub fn is_replaying() -> bool {
    REPLAY.with(|depth| depth.get() > 0)
}
#[inline]
pub(crate) fn primitive_random_word() {
    REPLAY.with(|depth| {
        if depth.get() > 0 {
            std::panic::resume_unwind(Box::new(LimitExceeded(
                "unexpected fresh RNG during family semantic replay".into(),
            )));
        }
    });
    primitive_work(1);
}
pub(crate) fn primitive_provenance(p: &crate::ir::DuplexProvenance) {
    ACTIVE.with(|s| {
        if let Some(b) = s.borrow_mut().last_mut() {
            if let Err(e) = b.charge_provenance(p) {
                std::panic::resume_unwind(Box::new(LimitExceeded(e)));
            }
        }
    });
}
#[inline]
pub(crate) fn primitive_bases(n: u64) {
    ACTIVE.with(|s| {
        if let Some(b) = s.borrow_mut().last_mut() {
            if let Err(e) = b.charge_bases(n) {
                std::panic::resume_unwind(Box::new(LimitExceeded(e)));
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bits_are_lossless_and_strict() {
        for v in [0.0, -0.0, f64::from_bits(1), 1.0, 1e-250] {
            let b = F64Bits::new(v).unwrap();
            let s = serde_json::to_string(&b).unwrap();
            assert_eq!(serde_json::from_str::<F64Bits>(&s).unwrap().0, v.to_bits());
        }
        assert!(serde_json::from_str::<F64Bits>("\"7ff0000000000000\"").is_err());
        assert!(serde_json::from_str::<F64Bits>("\"000000000000000A\"").is_err());
    }
    #[test]
    fn canonical_types_and_delimiters_are_distinct() {
        assert_ne!(
            canonical_digest("x", &json!(false)).unwrap(),
            canonical_digest("x", &json!(0)).unwrap()
        );
        assert_ne!(
            canonical_digest("x", &json!(0)).unwrap(),
            canonical_digest("x", &json!(0.0)).unwrap()
        );
        assert_ne!(
            canonical_digest("x", &json!(["a|b", "c"])).unwrap(),
            canonical_digest("x", &json!(["a", "b|c"])).unwrap()
        );
    }
    #[test]
    fn scope_fences_infallible_work_and_restores_context() {
        let mut b = FamilyReplayBudget::new(FamilyReplayLimits {
            max_work: 1,
            ..Default::default()
        });
        assert!(with_budget_scope(&mut b, |_| {
            primitive_work(2);
            Ok(())
        })
        .unwrap_err()
        .contains("resource limit"));
        primitive_work(100); // no scoped budget leaked
    }
}
