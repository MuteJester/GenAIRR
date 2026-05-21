use super::*;

#[test]
fn anchor_tail_no_trim_returns_anchor_onwards() {
    let seq = b"AAAAACGTGCN";
    let out = anchor_tail(seq, 7, 0, 0).expect("anchor fits, no trim");
    assert_eq!(out, &seq[7..]);
}

#[test]
fn anchor_tail_with_5prime_trim_under_anchor_keeps_anchor_onwards() {
    let seq = b"AAAAACGTGCN";
    let out = anchor_tail(seq, 7, 4, 0).expect("5' trim doesn't touch anchor");
    assert_eq!(out, &seq[7..]);
}

#[test]
fn anchor_tail_5prime_trim_past_anchor_is_none() {
    let seq = b"AAAAACGTGCN";
    assert!(anchor_tail(seq, 7, 8, 0).is_none());
}

#[test]
fn anchor_tail_3prime_trim_into_anchor_codon_is_none() {
    let seq = b"AAAAACGTGCN";
    assert!(anchor_tail(seq, 7, 0, 2).is_none());
}

#[test]
fn anchor_tail_trim_3_larger_than_seq_is_none_not_panic() {
    let seq = b"AAAAACGTGCN";
    assert!(anchor_tail(seq, 7, 0, 999).is_none());
}

#[test]
fn anchor_head_no_trim_returns_start_through_anchor_codon_end() {
    let seq = b"AAATGGNNN";
    let out = anchor_head(seq, 3, 0, 0).expect("anchor fits, no trim");
    assert_eq!(out, &seq[0..6]);
}

#[test]
fn anchor_head_with_5prime_trim_under_anchor_starts_at_trim_5() {
    let seq = b"AAATGGNNN";
    let out = anchor_head(seq, 3, 2, 0).expect("trim_5 stays under anchor");
    assert_eq!(out, &seq[2..6]);
}

#[test]
fn anchor_head_5prime_trim_past_anchor_is_none() {
    let seq = b"AAATGGNNN";
    assert!(anchor_head(seq, 3, 5, 0).is_none());
}

#[test]
fn anchor_head_3prime_trim_into_anchor_codon_is_none() {
    let seq = b"AAATGGNNN";
    assert!(anchor_head(seq, 3, 0, 4).is_none());
}
