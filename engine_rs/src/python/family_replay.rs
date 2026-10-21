//! Native lossless value codec, bounded parser and shared resource handles.
use crate::family_replay::{
    canonical_digest, strict_json, F64Bits, FamilyReplayBudget, FamilyReplayLimits,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBool, PyDict, PyFloat, PyInt, PyList, PyString};
use serde_json::{json, Value};
use std::collections::HashSet;
use std::sync::Mutex;

fn err(e: impl ToString) -> PyErr {
    PyValueError::new_err(e.to_string())
}
pub(crate) fn limits(text: Option<&str>) -> PyResult<FamilyReplayLimits> {
    match text {
        Some(s) => strict_json::from_str(s, &FamilyReplayLimits::default()).map_err(err),
        None => Ok(FamilyReplayLimits::default()),
    }
}
#[pyclass(name = "FamilyReplayBudget", module = "GenAIRR._engine")]
pub struct PyFamilyReplayBudget {
    pub inner: Mutex<FamilyReplayBudget>,
}
#[pymethods]
impl PyFamilyReplayBudget {
    #[new]
    #[pyo3(signature=(limits_json=None))]
    fn new(limits_json: Option<&str>) -> PyResult<Self> {
        Ok(Self {
            inner: Mutex::new(FamilyReplayBudget::new(limits(limits_json)?)),
        })
    }
    fn charge_work(&self, n: u64) -> PyResult<()> {
        self.inner.lock().map_err(err)?.charge_work(n).map_err(err)
    }
    fn charge_bases(&self, n: u64) -> PyResult<()> {
        self.inner.lock().map_err(err)?.charge_bases(n).map_err(err)
    }
    fn charge_nodes(&self, n: u64) -> PyResult<()> {
        self.inner.lock().map_err(err)?.charge_nodes(n).map_err(err)
    }
    fn charge_choices(&self, n: u64) -> PyResult<()> {
        self.inner
            .lock()
            .map_err(err)?
            .charge_choices(n)
            .map_err(err)
    }
    fn charge_attempts(&self, n: u64) -> PyResult<()> {
        self.inner
            .lock()
            .map_err(err)?
            .charge_attempts(n)
            .map_err(err)
    }
    fn snapshot(&self) -> PyResult<String> {
        self.inner.lock().map_err(err)?.snapshot().map_err(err)
    }
    fn limits_json(&self) -> PyResult<String> {
        serde_json::to_string(&self.inner.lock().map_err(err)?.limits).map_err(err)
    }
}

struct Pack<'a> {
    limits: &'a FamilyReplayLimits,
    values: u64,
    ancestors: HashSet<usize>,
}
impl Pack<'_> {
    fn value(&mut self, v: &Bound<'_, PyAny>, depth: u64, path: &str) -> PyResult<Value> {
        if depth > self.limits.max_depth {
            return Err(err(format!("{path}: resource limit: value depth")));
        }
        self.values = self
            .values
            .checked_add(1)
            .ok_or_else(|| err("value count overflow"))?;
        if self.values > self.limits.max_values {
            return Err(err(format!("{path}: resource limit: values")));
        }
        if v.is_none() {
            return Ok(json!({"t":"null"}));
        }
        if v.is_instance_of::<PyBool>() {
            return Ok(json!({"t":"bool","v":v.extract::<bool>()?}));
        }
        if v.is_instance_of::<PyInt>() {
            let decimal = if let Ok(n) = v.extract::<i64>() {
                n.to_string()
            } else {
                v.extract::<u64>()
                    .map_err(|_| err(format!("{path}: integer outside [-2**63,2**64-1]")))?
                    .to_string()
            };
            return Ok(json!({"t":"int","v":decimal}));
        }
        if v.is_instance_of::<PyFloat>() {
            return Ok(
                json!({"t":"f64","v":F64Bits::new(v.extract::<f64>()?).map_err(|e|err(format!("{path}: {e}")))?}),
            );
        }
        if let Ok(s) = v.downcast::<PyString>() {
            let s = s.extract::<String>()?;
            if s.len() as u64 > self.limits.max_string_bytes {
                return Err(err(format!("{path}: resource limit: string")));
            }
            return Ok(json!({"t":"str","v":s}));
        }
        let pointer = v.as_ptr() as usize;
        if !self.ancestors.insert(pointer) {
            return Err(err(format!("{path}: cyclic value")));
        }
        let result = (|| {
            if let Ok(d) = v.downcast::<PyDict>() {
                let mut entries = Vec::new();
                for (key, value) in d.iter() {
                    if !key.is_instance_of::<PyString>() {
                        return Err(err(format!("{path}: object keys must be strings")));
                    }
                    let key = key.extract::<String>()?;
                    if key.len() as u64 > self.limits.max_string_bytes {
                        return Err(err(format!("{path}: resource limit: key")));
                    }
                    let child = self.value(&value, depth + 1, &format!("{path}.{key}"))?;
                    entries.push((key, child));
                }
                entries.sort_unstable_by(|a, b| a.0.as_bytes().cmp(b.0.as_bytes()));
                return Ok(json!({"t":"object","v":entries}));
            }
            // JSON arrays are lists. Tuples are deliberately not silently
            // coerced: caller metadata must already use the declared subset.
            if let Ok(a) = v.downcast::<PyList>() {
                let mut values = Vec::new();
                for (i, value) in a.iter().enumerate() {
                    values.push(self.value(&value, depth + 1, &format!("{path}[{i}]"))?);
                }
                return Ok(json!({"t":"array","v":values}));
            }
            Err(err(format!("{path}: unsupported JSON value type")))
        })();
        self.ancestors.remove(&pointer);
        result
    }
}
#[pyfunction]
#[pyo3(signature=(value,limits_json=None))]
pub fn _family_replay_pack(
    value: &Bound<'_, PyAny>,
    limits_json: Option<&str>,
) -> PyResult<String> {
    let limits = limits(limits_json)?;
    let packed = Pack {
        limits: &limits,
        values: 0,
        ancestors: HashSet::new(),
    }
    .value(value, 0, "$")?;
    let text = serde_json::to_string(&packed).map_err(err)?;
    if text.len() as u64 > limits.max_file_bytes {
        return Err(err("resource limit: encoded bytes"));
    }
    strict_json::parse(&text, &limits).map_err(err)?;
    Ok(text)
}
fn unpack(
    py: Python<'_>,
    v: &Value,
    depth: u64,
    limits: &FamilyReplayLimits,
) -> PyResult<PyObject> {
    if depth > limits.max_depth {
        return Err(err("resource limit: dynamic value depth"));
    }
    let object = v
        .as_object()
        .ok_or_else(|| err("tagged value must be an object"))?;
    let tag = object
        .get("t")
        .and_then(Value::as_str)
        .ok_or_else(|| err("missing dynamic type tag"))?;
    let expected = if tag == "null" { 1 } else { 2 };
    if object.len() != expected || (expected == 2 && !object.contains_key("v")) {
        return Err(err("unknown/missing tagged value fields"));
    }
    let value = object.get("v").unwrap_or(&Value::Null);
    match tag {
        "null" => Ok(py.None()),
        "bool" => Ok(value
            .as_bool()
            .ok_or_else(|| err("invalid tagged Boolean"))?
            .into_py(py)),
        "str" => Ok(value
            .as_str()
            .ok_or_else(|| err("invalid tagged string"))?
            .into_py(py)),
        "int" => {
            let s = value
                .as_str()
                .ok_or_else(|| err("integer needs decimal string"))?;
            if let Ok(n) = s.parse::<i64>() {
                if n.to_string() != s {
                    return Err(err("noncanonical integer"));
                }
                Ok(n.into_py(py))
            } else {
                let n = s
                    .parse::<u64>()
                    .map_err(|_| err("integer outside 64-bit range"))?;
                if n.to_string() != s {
                    return Err(err("noncanonical integer"));
                }
                Ok(n.into_py(py))
            }
        }
        "f64" => {
            let bits: F64Bits = serde_json::from_value(value.clone()).map_err(err)?;
            Ok(bits.get().into_py(py))
        }
        "array" => {
            let a = value
                .as_array()
                .ok_or_else(|| err("invalid tagged array"))?;
            let mut values = Vec::new();
            for item in a {
                values.push(unpack(py, item, depth + 1, limits)?);
            }
            Ok(PyList::new_bound(py, values).into_py(py))
        }
        "object" => {
            let a = value
                .as_array()
                .ok_or_else(|| err("invalid tagged object"))?;
            let d = PyDict::new_bound(py);
            let mut keys = HashSet::new();
            for pair in a {
                let pair = pair
                    .as_array()
                    .filter(|p| p.len() == 2)
                    .ok_or_else(|| err("invalid tagged object entry"))?;
                let key = pair[0]
                    .as_str()
                    .ok_or_else(|| err("object key must be a string"))?;
                if !keys.insert(key) {
                    return Err(err("duplicate tagged object key"));
                }
                d.set_item(key, unpack(py, &pair[1], depth + 1, limits)?)?;
            }
            Ok(d.into_py(py))
        }
        _ => Err(err(format!("unsupported dynamic type tag: {tag}"))),
    }
}
#[pyfunction]
#[pyo3(signature=(text,limits_json=None))]
pub fn _family_replay_unpack(
    py: Python<'_>,
    text: &str,
    limits_json: Option<&str>,
) -> PyResult<PyObject> {
    let limits = limits(limits_json)?;
    let value = strict_json::parse(text, &limits).map_err(err)?;
    unpack(py, &value, 0, &limits)
}
#[pyfunction]
#[pyo3(signature=(text,limits_json=None))]
pub fn _family_replay_parse_json(text: &str, limits_json: Option<&str>) -> PyResult<String> {
    let value = strict_json::parse(text, &limits(limits_json)?).map_err(err)?;
    serde_json::to_string(&value).map_err(err)
}
#[pyfunction]
pub fn _family_replay_digest(text: &str, domain: &str) -> PyResult<String> {
    let value = strict_json::parse(text, &FamilyReplayLimits::default()).map_err(err)?;
    canonical_digest(domain, &value).map_err(err)
}
#[pyfunction]
pub fn _family_replay_numeric_profile() -> PyResult<String> {
    let probes: Vec<_> = [-700.0_f64, -64.0, -1.0, -0.0, 0.0, 0.01, 1.0, 5.0, 64.0]
        .into_iter()
        .map(|n| {
            json!([
                format!("{:016x}", n.to_bits()),
                format!("{:016x}", n.exp().to_bits())
            ])
        })
        .collect();
    serde_json::to_string(&json!({"id":"ieee754_binary64_exp_v1","architecture":std::env::consts::ARCH,"os":std::env::consts::OS,"probes":probes})).map_err(err)
}
// Output materialization is part of the same caller-owned budget as replay.
// Simulations share Arc-backed pools/provenance. Charge the containers, strings
// and event ledgers that Outcome::clone actually allocates, not shared payloads.
#[derive(Default)]
struct OutcomeCopyCharge {
    entries: u64,
    bytes: u64,
    work: u64,
}
impl OutcomeCopyCharge {
    fn array<T>(&mut self, len: usize) -> Result<(), String> {
        let n = u64::try_from(len).map_err(|_| "history size overflow")?;
        self.entries = self
            .entries
            .checked_add(n)
            .ok_or("history entry overflow")?;
        self.bytes = self
            .bytes
            .checked_add(
                n.checked_mul(std::mem::size_of::<T>() as u64)
                    .ok_or("history byte overflow")?,
            )
            .ok_or("history byte overflow")?;
        self.work = self.work.checked_add(n).ok_or("history work overflow")?;
        Ok(())
    }
    fn string(&mut self, len: usize) -> Result<(), String> {
        self.bytes = self
            .bytes
            .checked_add(len as u64)
            .ok_or("history string overflow")?;
        self.work = self
            .work
            .checked_add(len as u64 + 1)
            .ok_or("history work overflow")?;
        Ok(())
    }
    fn outcome(&mut self, outcome: &crate::pass::Outcome) -> Result<(), String> {
        self.array::<crate::ir::Simulation>(outcome.revisions.len())?;
        self.array::<String>(outcome.pass_names.len())?;
        for name in &outcome.pass_names {
            self.string(name.len())?;
        }
        self.array::<crate::trace::ChoiceRecord>(outcome.trace.len())?;
        for choice in outcome.trace.choices() {
            self.string(choice.address.len())?;
            if let crate::trace::ChoiceValue::Bases(value) = &choice.value {
                self.string(value.len())?;
            }
        }
        self.array::<crate::trace::SamplingDiagnostic>(outcome.trace.diagnostics().len())?;
        for diagnostic in outcome.trace.diagnostics() {
            for value in [
                &diagnostic.pass_name,
                &diagnostic.model,
                &diagnostic.strategy,
                &diagnostic.status,
            ] {
                self.string(value.len())?;
            }
        }
        self.array::<crate::event::EventRecord>(outcome.events.len())?;
        for event in &outcome.events {
            self.string(event.pass_name.len())?;
            self.array::<crate::pass::PassCompileEffect>(event.compile_effects.len())?;
            self.array::<crate::ir::SimulationEvent>(event.simulation_events.len())?;
            for change in &event.simulation_events {
                if let crate::ir::SimulationEvent::DuplexEditContextRecorded {
                    source_pass: Some(value),
                } = change
                {
                    self.string(value.len())?;
                }
            }
        }
        Ok(())
    }
    fn charge(self, budget: &mut FamilyReplayBudget, multiplier: u64) -> Result<(), String> {
        budget.charge_work(
            self.work
                .checked_mul(multiplier)
                .ok_or("history work overflow")?,
        )?;
        budget.charge_history_copy(
            self.entries
                .checked_mul(multiplier)
                .ok_or("history entry overflow")?,
            self.bytes
                .checked_mul(multiplier)
                .ok_or("history byte overflow")?,
        )
    }
}
fn scoped<T>(
    budget: &PyFamilyReplayBudget,
    f: impl FnOnce(&mut FamilyReplayBudget) -> Result<T, String>,
) -> PyResult<T> {
    let mut guard = budget
        .inner
        .lock()
        .map_err(|_| err("budget lock poisoned"))?;
    crate::family_replay::with_budget_scope(&mut guard, f).map_err(err)
}
#[pyfunction]
pub fn _family_replay_has_observed(
    family: &super::lineage::PyFamilyOutcome,
    budget: &PyFamilyReplayBudget,
) -> PyResult<bool> {
    scoped(budget, |b| {
        b.charge_work(family.node_outcomes.len() as u64)?;
        Ok(family.node_outcomes.iter().any(Option::is_some))
    })
}
#[pyfunction]
pub fn _family_replay_observed_outcomes(
    family: &super::lineage::PyFamilyOutcome,
    budget: &PyFamilyReplayBudget,
) -> PyResult<Vec<super::outcome::PyOutcome>> {
    scoped(budget, |b| {
        let mut charge = OutcomeCopyCharge::default();
        charge.array::<super::outcome::PyOutcome>(family.node_outcomes.iter().flatten().count())?;
        for outcome in family.node_outcomes.iter().flatten() {
            charge.outcome(outcome)?;
        }
        charge.charge(b, 1)?;
        Ok(family
            .node_outcomes
            .iter()
            .flatten()
            .cloned()
            .map(super::outcome::PyOutcome::new)
            .collect())
    })
}
#[pyfunction]
pub fn _family_replay_node_outcomes(
    family: &super::lineage::PyFamilyOutcome,
    budget: &PyFamilyReplayBudget,
) -> PyResult<Vec<Option<super::outcome::PyOutcome>>> {
    scoped(budget, |b| {
        let mut charge = OutcomeCopyCharge::default();
        charge.array::<Option<super::outcome::PyOutcome>>(family.node_outcomes.len())?;
        for outcome in family.node_outcomes.iter().flatten() {
            charge.outcome(outcome)?;
        }
        charge.charge(b, 1)?;
        Ok(family
            .node_outcomes
            .iter()
            .map(|o| o.clone().map(super::outcome::PyOutcome::new))
            .collect())
    })
}
#[pyfunction]
pub fn _family_replay_merge_outcomes(
    base: &super::outcome::PyOutcome,
    corruption: &super::outcome::PyOutcome,
    budget: &PyFamilyReplayBudget,
) -> PyResult<super::outcome::PyOutcome> {
    scoped(budget, |b| {
        let mut charge = OutcomeCopyCharge::default();
        charge.outcome(&base.inner)?;
        charge.outcome(&corruption.inner)?;
        // Vec::extend may grow the initially cloned ledger allocations; reserve
        // a conservative second copy before invoking the ordinary merge path.
        charge.charge(b, 2)?;
        super::lineage::merge_lineage_corruption(base, corruption).map_err(|e| e.to_string())
    })
}
fn projection_charge(
    outcome: &crate::pass::Outcome,
    budget: &mut FamilyReplayBudget,
) -> Result<(), String> {
    let state = outcome.final_simulation();
    budget.charge_bases(state.pool.len() as u64)?;
    budget.charge_history(state)?;
    let mut work = OutcomeCopyCharge::default();
    work.outcome(outcome)?;
    budget.charge_work(work.work)?;
    Ok(())
}
#[pyfunction]
pub fn _family_replay_airr_record<'py>(
    py: Python<'py>,
    outcome: &super::outcome::PyOutcome,
    refdata: &super::refdata::PyRefDataConfig,
    sequence_id: &str,
    budget: &PyFamilyReplayBudget,
) -> PyResult<Bound<'py, PyDict>> {
    scoped(budget, |b| {
        projection_charge(&outcome.inner, b)?;
        outcome
            .airr_record(py, refdata, sequence_id)
            .map_err(|e| e.to_string())
    })
}
#[pyfunction]
pub fn _family_replay_validate_record<'py>(
    py: Python<'py>,
    outcome: &super::outcome::PyOutcome,
    refdata: &super::refdata::PyRefDataConfig,
    sequence_id: &str,
    budget: &PyFamilyReplayBudget,
) -> PyResult<Vec<Bound<'py, PyDict>>> {
    scoped(budget, |b| {
        projection_charge(&outcome.inner, b)?;
        outcome
            .validate_record(py, refdata, sequence_id)
            .map_err(|e| e.to_string())
    })
}
#[pyfunction]
pub fn _family_replay_final_simulation(
    outcome: &super::outcome::PyOutcome,
    budget: &PyFamilyReplayBudget,
) -> PyResult<super::simulation::PySimulation> {
    scoped(budget, |b| {
        b.charge_work(1)?;
        b.charge_history_copy(1, std::mem::size_of::<crate::ir::Simulation>() as u64)?;
        Ok(super::simulation::PySimulation::new(
            outcome.inner.final_simulation().clone(),
        ))
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_family_replay_has_observed, m)?)?;
    m.add_function(wrap_pyfunction!(_family_replay_observed_outcomes, m)?)?;
    m.add_function(wrap_pyfunction!(_family_replay_node_outcomes, m)?)?;
    m.add_function(wrap_pyfunction!(_family_replay_merge_outcomes, m)?)?;
    m.add_function(wrap_pyfunction!(_family_replay_airr_record, m)?)?;
    m.add_function(wrap_pyfunction!(_family_replay_validate_record, m)?)?;
    m.add_function(wrap_pyfunction!(_family_replay_final_simulation, m)?)?;
    m.add_class::<PyFamilyReplayBudget>()?;
    m.add_function(wrap_pyfunction!(_family_replay_pack, m)?)?;
    m.add_function(wrap_pyfunction!(_family_replay_unpack, m)?)?;
    m.add_function(wrap_pyfunction!(_family_replay_digest, m)?)?;
    m.add_function(wrap_pyfunction!(_family_replay_parse_json, m)?)?;
    m.add_function(wrap_pyfunction!(_family_replay_numeric_profile, m)?)?;
    Ok(())
}

#[cfg(test)]
mod output_copy_tests {
    use super::*;
    #[test]
    fn copied_event_and_trace_containers_are_charged_separately_from_shared_states() {
        let mut trace = crate::trace::Trace::new();
        trace.record("site", crate::trace::ChoiceValue::Int(1));
        trace.record("bases", crate::trace::ChoiceValue::Bases(b"ACGT".to_vec()));
        let outcome = crate::pass::Outcome {
            revisions: vec![crate::ir::Simulation::new()],
            pass_names: vec!["pass".into()],
            trace,
            events: vec![],
        };
        let mut charge = OutcomeCopyCharge::default();
        charge.outcome(&outcome).unwrap();
        assert_eq!(charge.entries, 4); // one revision, one name, two choices
        assert!(charge.bytes > std::mem::size_of::<crate::ir::Simulation>() as u64);
        let mut limits = FamilyReplayLimits::default();
        limits.max_history_entries = 3;
        assert!(charge
            .charge(&mut FamilyReplayBudget::new(limits), 1)
            .unwrap_err()
            .contains("resource limit"));
    }
}
