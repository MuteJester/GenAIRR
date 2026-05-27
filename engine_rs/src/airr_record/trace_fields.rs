use crate::address::ChoiceAddress;
use crate::trace::{ChoiceValue, Trace};

pub(super) fn trace_int_choice(trace: &Trace, address: ChoiceAddress) -> i64 {
    match trace.find_choice(address) {
        Some(rec) => match rec.value {
            ChoiceValue::Int(v) => v,
            _ => 0,
        },
        None => 0,
    }
}

pub(super) fn trace_bool_choice(trace: &Trace, address: ChoiceAddress) -> bool {
    match trace.find_choice(address) {
        Some(rec) => match rec.value {
            ChoiceValue::Bool(v) => v,
            _ => false,
        },
        None => false,
    }
}
