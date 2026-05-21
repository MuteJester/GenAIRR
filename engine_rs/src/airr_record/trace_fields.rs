use crate::address;
use crate::trace::{ChoiceValue, Trace};

pub(super) fn trace_int(trace: &Trace, address: &str) -> i64 {
    match trace.find(address) {
        Some(rec) => match rec.value {
            ChoiceValue::Int(v) => v,
            _ => 0,
        },
        None => 0,
    }
}

pub(super) fn trace_bool(trace: &Trace, address: &str) -> bool {
    match trace.find(address) {
        Some(rec) => match rec.value {
            ChoiceValue::Bool(v) => v,
            _ => false,
        },
        None => false,
    }
}

pub(super) fn mutation_count(trace: &Trace) -> i64 {
    if let Some(r) = trace.find(address::MUTATE_S5F_COUNT) {
        if let ChoiceValue::Int(v) = r.value {
            return v;
        }
    }
    if let Some(r) = trace.find(address::MUTATE_UNIFORM_COUNT) {
        if let ChoiceValue::Int(v) = r.value {
            return v;
        }
    }
    0
}
