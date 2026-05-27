//! v3.0 acceptance test: `admissible_bases_at` mask MUST agree with
//! `admits_with_context(TargetedBaseSubstitution)` for every
//! canonical base at every relevant site.
//!
//! This is the regression guard for the architectural invariant:
//!
//! > **Support comes from contracts. Probabilities come from the
//! > natural pass distribution restricted to that support.**
//!
//! If `mask.admits(b)` ever disagrees with the existing predicate
//! API, the constrain-before-propose sampler will pick a base the
//! reject-after-propose path would reject (or vice versa), and the
//! invariant is broken. We don't enforce this at runtime (the cost
//! would defeat the whole point of the cache); we enforce it here.

use crate::contract::test_support::{make_assembled_sim_from_refdata, make_vj_with_anchor_codons};
use crate::contract::{
    AnchorPreserved, BaseMask, ChoiceContext, Contract, ContractSet, NoStopCodonInJunction,
};
use crate::ir::{NucHandle, Segment};
use crate::trace::ChoiceValue;

const CANONICAL_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const NO_STOP_ADDRESS: &str = "mutate.uniform.base[0]";

/// Assert that for every canonical base, the mask's verdict matches
/// the predicate's verdict at the given `(sim, site)`. Builds a
/// `TargetedBaseSubstitution` choice context so the predicate
/// dispatches to the per-position substitution logic.
fn assert_mask_matches_predicate<C: Contract>(
    contract: &C,
    sim: &crate::ir::Simulation,
    refdata: &crate::refdata::RefDataConfig,
    site: NucHandle,
    address: &str,
) {
    let mask = contract.admissible_bases_at(sim, Some(refdata), site);
    for &base in &CANONICAL_BASES {
        let context = ChoiceContext::targeted_base_substitution(0, 1, site)
            .with_address_if_missing(crate::address::ChoiceAddress::parse(address));
        let predicate_ok = contract
            .admits_typed(sim, Some(refdata), context, &ChoiceValue::Base(base))
            .is_ok();
        let mask_ok = mask.admits(base);
        assert_eq!(
            predicate_ok,
            mask_ok,
            "mask/predicate mismatch for {} at site {}: predicate={} mask={} (mask bits 0b{:04b})",
            base as char,
            site.index(),
            predicate_ok,
            mask_ok,
            mask.0,
        );
    }
}

// ──────────────────────────────────────────────────────────────────
// NoStopCodonInJunction
// ──────────────────────────────────────────────────────────────────

#[test]
fn no_stop_codon_mask_matches_predicate_at_every_junction_site() {
    // V anchor codon = "TGT" (Cys), J anchor codon = "TGG" (Trp).
    // Junction is `[6, 12)` — 6 bases, in-frame. The pool around
    // the junction has many possible stop-creating substitutions.
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let contract = NoStopCodonInJunction::new();

    // Junction is pool positions [6, 12). Test every position
    // inside it.
    for pos in 6..12 {
        assert_mask_matches_predicate(&contract, &sim, &cfg, NucHandle::new(pos), NO_STOP_ADDRESS);
    }
}

#[test]
fn no_stop_codon_mask_is_unconstrained_outside_junction() {
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let contract = NoStopCodonInJunction::new();

    // Pool positions outside [6, 12) — V framework + first 0..6
    // bases, J framework after the anchor codon — should all be
    // unconstrained (no opinion).
    for pos in 0..6 {
        let mask = contract.admissible_bases_at(&sim, Some(&cfg), NucHandle::new(pos));
        assert_eq!(
            mask,
            BaseMask::UNCONSTRAINED,
            "site {} (V framework, out of junction) should be unconstrained, got 0b{:04b}",
            pos,
            mask.0
        );
        // And predicate agrees: every canonical base admitted.
        assert_mask_matches_predicate(&contract, &sim, &cfg, NucHandle::new(pos), NO_STOP_ADDRESS);
    }
}

#[test]
fn no_stop_codon_mask_correctly_excludes_stop_creators() {
    // Construct a junction where some single-base substitutions
    // create stop codons. Anchor codons TGT (Cys) + TGG (Trp)
    // mean junction = "TGT" + "TGG" = "TGTTGG", which translates
    // to Cys-Trp — no stops. Substituting at position 7 (the G
    // of TGT) to A would give "TAT" (Tyr, no stop). To create a
    // stop, we'd need e.g. position 6 substitution to keep T but
    // make the next codon stop-shaped. Specifically: at junction
    // codon 1 (positions 9..12 = "TGG", Trp), substituting
    // position 10 (the first G) to A would yield "TAG" — a stop.
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let contract = NoStopCodonInJunction::new();

    let mask = contract.admissible_bases_at(&sim, Some(&cfg), NucHandle::new(10));
    // 'A' at position 10 would make codon "TAG" → stop. Mask must
    // not contain 'A'. Other bases (C, G, T) should be admissible
    // unless they also create stops at this codon.
    assert!(
        !mask.admits(b'A'),
        "site 10 substitute A would create TAG stop; mask must exclude it (got 0b{:04b})",
        mask.0
    );
}

#[test]
fn no_stop_codon_mask_unconstrained_when_no_refdata() {
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let contract = NoStopCodonInJunction::new();

    // No refdata → contract has no opinion → unconstrained.
    let mask = contract.admissible_bases_at(&sim, None, NucHandle::new(7));
    assert_eq!(mask, BaseMask::UNCONSTRAINED);
}

// ──────────────────────────────────────────────────────────────────
// AnchorPreserved
// ──────────────────────────────────────────────────────────────────

#[test]
fn anchor_preserved_mask_matches_predicate_at_anchor_codon() {
    // V anchor = "TGT" (Cys) at V_pool[6..9]. AnchorPreserved.V
    // guards changes that would alter the anchor amino acid.
    // The codon table tells us synonymous Cys codons are TGT and
    // TGC. So substituting position 8 (T → C) preserves Cys.
    // Other substitutions at positions 6/7/8 may change it.
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let v_anchor = AnchorPreserved::new(Segment::V);

    // V anchor codon lives at pool positions [6, 9).
    for pos in 6..9 {
        assert_mask_matches_predicate(
            &v_anchor,
            &sim,
            &cfg,
            NucHandle::new(pos),
            "mutate.uniform.base[0]",
        );
    }
}

#[test]
fn anchor_preserved_mask_unconstrained_outside_anchor() {
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let v_anchor = AnchorPreserved::new(Segment::V);

    // V anchor is at pool [6, 9). Sites outside that codon
    // should be unconstrained.
    for pos in [0, 5, 9, 11] {
        let mask = v_anchor.admissible_bases_at(&sim, Some(&cfg), NucHandle::new(pos));
        assert_eq!(
            mask,
            BaseMask::UNCONSTRAINED,
            "site {} (outside V anchor codon) should be unconstrained, got 0b{:04b}",
            pos,
            mask.0
        );
    }
}

#[test]
fn anchor_preserved_v_admits_synonymous_substitutions() {
    // Anchor codon TGT (Cys) is synonymous with TGC. So at
    // position 8 (the third base of the codon), bases C and T
    // both preserve Cys; A and G change the amino acid.
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let v_anchor = AnchorPreserved::new(Segment::V);

    let mask = v_anchor.admissible_bases_at(&sim, Some(&cfg), NucHandle::new(8));
    assert!(
        mask.admits(b'T'),
        "T preserves TGT (Cys), expected admitted"
    );
    assert!(
        mask.admits(b'C'),
        "C makes TGC (still Cys), expected admitted"
    );
    assert!(!mask.admits(b'A'), "A makes TGA (stop) — not Cys");
    assert!(!mask.admits(b'G'), "G makes TGG (Trp) — not Cys");
}

#[test]
fn anchor_preserved_mask_unconstrained_when_no_refdata() {
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);
    let v_anchor = AnchorPreserved::new(Segment::V);

    let mask = v_anchor.admissible_bases_at(&sim, None, NucHandle::new(7));
    assert_eq!(mask, BaseMask::UNCONSTRAINED);
}

// ──────────────────────────────────────────────────────────────────
// ContractSet composition: bundle = intersection of masks
// ──────────────────────────────────────────────────────────────────

#[test]
fn contract_set_intersects_no_stop_and_anchor_masks_at_anchor_codon_inside_junction() {
    // The V anchor codon at pool [6, 9) sits inside the junction
    // [6, 12). So at site 8 (third base of anchor codon = junction
    // codon 0's third base), BOTH contracts have an opinion:
    //   - AnchorPreserved.V: only T or C admissible (Cys synonyms)
    //   - NoStopCodonInJunction: any base that doesn't make a
    //     junction stop
    // The bundle's admissible mask must be the intersection
    // (bitwise AND).
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);

    let no_stop = NoStopCodonInJunction::new();
    let anchor_v = AnchorPreserved::new(Segment::V);

    let mut bundle = ContractSet::new();
    bundle.add(Box::new(NoStopCodonInJunction::new()));
    bundle.add(Box::new(AnchorPreserved::new(Segment::V)));

    let site = NucHandle::new(8);
    let m_no_stop = no_stop.admissible_bases_at(&sim, Some(&cfg), site);
    let m_anchor = anchor_v.admissible_bases_at(&sim, Some(&cfg), site);
    let m_bundle = bundle.admissible_bases_at(&sim, Some(&cfg), site);
    assert_eq!(
        m_bundle,
        BaseMask(m_no_stop.0 & m_anchor.0),
        "bundle mask should be AND of contract masks"
    );
}

#[test]
fn contract_set_bundle_matches_predicate_everywhere_in_junction() {
    // The composability invariant: bundle.admissible_bases_at(site).admits(b)
    // must equal bundle.admits_with_context(..., TargetedBaseSubstitution, ...).is_ok()
    // for every site in the relevant range and every canonical base.
    let cfg = make_vj_with_anchor_codons(b"TGT", b"TGG");
    let sim = make_assembled_sim_from_refdata(&cfg);

    let mut bundle = ContractSet::new();
    bundle.add(Box::new(NoStopCodonInJunction::new()));
    bundle.add(Box::new(AnchorPreserved::new(Segment::V)));
    bundle.add(Box::new(AnchorPreserved::new(Segment::J)));

    for pos in 0..12 {
        let site = NucHandle::new(pos);
        let mask = bundle.admissible_bases_at(&sim, Some(&cfg), site);
        for &base in &CANONICAL_BASES {
            let context = ChoiceContext::targeted_base_substitution(0, 1, site)
                .with_address_if_missing(crate::address::ChoiceAddress::parse(NO_STOP_ADDRESS));
            let predicate_ok = bundle
                .admits_typed(&sim, Some(&cfg), context, &ChoiceValue::Base(base))
                .is_ok();
            let mask_ok = mask.admits(base);
            assert_eq!(
                predicate_ok, mask_ok,
                "bundle mask/predicate mismatch at site {} base {}: predicate={} mask={}",
                pos, base as char, predicate_ok, mask_ok
            );
        }
    }
}
