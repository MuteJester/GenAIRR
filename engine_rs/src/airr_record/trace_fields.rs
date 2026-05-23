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
