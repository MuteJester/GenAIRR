//! v3.0 review fix: pin the `admits_fixed_base_at` override against
//! drift. The default trait impl returns `true` for every site/byte;
//! the override walks junction codons with the byte substituted at
//! `site` and rejects iff a resulting codon translates to a stop.

use super::*;

/// V anchor codon = `TAC` (Y), J anchor codon = `GGG` (G). Junction
/// is `[V_anchor_pool=6, J_anchor_pool+3=12)`. Position 8 is the
/// third base of the V anchor codon — under junction codon walk
/// from start 6, the codon is `TA·`, so `A` at site 8 creates a
/// TAA stop.
fn fixture() -> (RefDataConfig, Simulation) {
    let cfg = make_vj_with_anchor_codons(b"TAC", b"GGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    (cfg, sim)
}

#[test]
fn admits_fixed_base_at_uppercase_stop_creating_base_is_rejected() {
    let (cfg, sim) = fixture();
    let c = NoStopCodonInJunction::new();
    assert!(!c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'A'));
}

#[test]
fn admits_fixed_base_at_lowercase_stop_creating_base_is_rejected() {
    // `translate_codon` is case-insensitive — a lowercase `a`
    // candidate at the same junction site must yield the same
    // verdict as uppercase.
    let (cfg, sim) = fixture();
    let c = NoStopCodonInJunction::new();
    assert!(!c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'a'));
}

#[test]
fn admits_fixed_base_at_safe_base_is_admitted() {
    let (cfg, sim) = fixture();
    let c = NoStopCodonInJunction::new();
    assert!(c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'C'));
    assert!(c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'c'));
}

#[test]
fn admits_fixed_base_at_non_canonical_byte_does_not_produce_stop() {
    // Non-canonical bytes (`N`, IUPAC ambiguities) translate to `X`,
    // not a stop. The contract has no opinion → admitted. (Any
    // additional restriction on `N` writes is a different
    // contract's concern, not this one.)
    let (cfg, sim) = fixture();
    let c = NoStopCodonInJunction::new();
    assert!(c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(8), b'N'));
}

#[test]
fn admits_fixed_base_at_out_of_junction_is_unconstrained() {
    let (cfg, sim) = fixture();
    let c = NoStopCodonInJunction::new();
    // Site 0 is V[0], deep outside the junction window [6, 12).
    assert!(c.admits_fixed_base_at(&sim, Some(&cfg), NucHandle::new(0), b'A'));
}

#[test]
fn admits_fixed_base_at_no_refdata_is_unconstrained() {
    let (_cfg, sim) = fixture();
    let c = NoStopCodonInJunction::new();
    assert!(c.admits_fixed_base_at(&sim, None, NucHandle::new(8), b'A'));
}

#[test]
fn admits_fixed_base_at_agrees_with_mask_for_canonical_bases() {
    // Pin the invariant the v3.0 design depends on: at every
    // canonical base, `admits_fixed_base_at` and
    // `admissible_bases_at(site).admits(b)` must agree. This is
    // the equivalence that lets the substitution helper trust the
    // mask for canonical writes and only consult
    // `admits_fixed_base_at` for transforms that change the byte.
    let (cfg, sim) = fixture();
    let c = NoStopCodonInJunction::new();
    for site in 0..(sim.pool.len() as u32) {
        let handle = NucHandle::new(site);
        let mask = c.admissible_bases_at(&sim, Some(&cfg), handle);
        for &b in &[b'A', b'C', b'G', b'T'] {
            assert_eq!(
                mask.admits(b),
                c.admits_fixed_base_at(&sim, Some(&cfg), handle, b),
                "mismatch at site {} base {}",
                site,
                b as char
            );
        }
    }
}
