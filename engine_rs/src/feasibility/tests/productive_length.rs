use super::*;

#[test]
fn productive_length_in_frame_with_no_forced_stops_is_feasible() {
    let v_tail = b"TGT";
    let j_head = b"TGG";
    assert!(productive_length_is_feasible(v_tail, 6, j_head));
}

#[test]
fn productive_length_out_of_frame_is_not_feasible() {
    let v_tail = b"TGT";
    let j_head = b"TGG";
    assert!(!productive_length_is_feasible(v_tail, 5, j_head));
}

#[test]
fn productive_length_with_forced_stop_in_known_region_is_not_feasible() {
    let v_tail = b"TAA";
    let j_head = b"TGG";
    assert!(!productive_length_is_feasible(v_tail, 0, j_head));
}

#[test]
fn no_known_stop_treats_np_bases_as_wildcards() {
    let v_tail = b"T";
    let j_head = b"GG";
    assert!(no_known_stop_for_length(v_tail, 2, j_head));
}

#[test]
fn no_known_stop_detects_stop_in_fully_known_v_tail() {
    let v_tail = b"TAA";
    let j_head = b"GGG";
    assert!(!no_known_stop_for_length(v_tail, 0, j_head));
}

#[test]
fn no_known_stop_detects_stop_in_fully_known_j_head() {
    let v_tail = b"TGT";
    let j_head = b"TAG";
    assert!(!no_known_stop_for_length(v_tail, 0, j_head));
}
