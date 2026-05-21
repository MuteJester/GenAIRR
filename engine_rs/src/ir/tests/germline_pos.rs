use super::*;

// ── GermlinePos newtype invariants ───────────────────────────────

#[test]
#[should_panic(expected = "u16::MAX is reserved for NONE")]
fn germline_pos_pos_rejects_max() {
    // u16::MAX is reserved for the NONE sentinel — constructing a
    // positioned GermlinePos with that value must panic, otherwise
    // a caller could silently create a value that compares equal
    // to NONE under derive(PartialEq).
    let _ = GermlinePos::pos(u16::MAX);
}

#[test]
fn germline_pos_none_projects_to_option_none() {
    assert_eq!(GermlinePos::NONE.get(), None);
    assert!(GermlinePos::NONE.is_none());
    assert!(!GermlinePos::NONE.is_some());
}

#[test]
fn germline_pos_pos_round_trips() {
    let p = GermlinePos::pos(42);
    assert_eq!(p.get(), Some(42));
    assert!(p.is_some());
    assert!(!p.is_none());
}

#[test]
fn nucleotide_size_unchanged() {
    // Pin the layout guarantee from `#[repr(transparent)]` on
    // GermlinePos: Nucleotide stays 6 bytes after the migration.
    // If this fails after a future change, the newtype lost its
    // niche / transparency and the migration regressed memory
    // cost on the hottest data structure in the engine.
    // (Plan's spec said 8 bytes; actual baseline was 6 -- both
    // pre- and post-migration measured 6 with the current field
    // set: base/germline/germline_pos/segment/flags. The point
    // of the assertion is unchangedness, not the magic number.)
    assert_eq!(std::mem::size_of::<Nucleotide>(), 6);
}
