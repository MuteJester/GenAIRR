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

/// Optional variant of [`trace_int_choice`]. Returns `None` when
/// the address is absent from the trace; returns `None` (not the
/// `0` sentinel) when the recorded value is the wrong kind so the
/// caller can distinguish "not opted in" from "recorded as 0". The
/// paired-end AIRR projection path needs this distinction — three
/// absent addresses means "no layout requested," not "every
/// length zero."
pub(super) fn trace_int_choice_opt(trace: &Trace, address: ChoiceAddress) -> Option<i64> {
    let rec = trace.find_choice(address)?;
    match rec.value {
        ChoiceValue::Int(v) => Some(v),
        _ => None,
    }
}
