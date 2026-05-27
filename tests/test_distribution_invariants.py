"""Distribution-invariant golden tests.

Companion to [docs/distribution_invariant_audit.md](../docs/distribution_invariant_audit.md).
Pins the central rule across all sampling mechanisms:

    Pr[draw = x] ∝ natural_weight(x) * 1[contracts.admit(x)]

i.e., natural distribution restricted to admissible support, with
implicit renormalization via inverse-CDF.

All tests use Monte Carlo with N=4000 deterministic seeds and a 5σ
tolerance, matching the Rust-side oracle tests'
(`indel_pass_count_1_under_productive_matches_exact_enumeration`,
`indel_pass_count_2_matches_exact_continuation_weighting`)
precedent.

5σ for a Bernoulli proportion at p, N:
    σ = sqrt(p * (1-p) / N)
    5σ = 5 * sqrt(p * (1-p) / N)

For p ≈ 0.5, N = 4000:  5σ ≈ ±0.0395 (≈ ±4 percentage points).
For p ≈ 0.55, N = 4000: 5σ ≈ ±0.0393.
"""
from __future__ import annotations

import math
from collections import Counter

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Fixtures and helpers
# ──────────────────────────────────────────────────────────────────


N_SEEDS = 4000


def _vj_basic() -> "ge.RefDataConfig":
    """V_anchor=6 in 9bp V (V_anchor_to_end=3, anchor codon TGT);
    J_anchor=0 in 6bp J (J_anchor_offset=0, anchor codon TGG).
    Junction length under NP1=L: 3 + L + 3 = L + 6.
    Frame-divisible iff L ≡ 0 (mod 3)."""
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TGGAAA", anchor=0)
    return cfg


def _within_5sigma(observed: float, expected: float, n: int) -> bool:
    """Bernoulli 5σ check: |obs - exp| ≤ 5 * sqrt(exp * (1-exp) / n)."""
    if expected <= 0.0 or expected >= 1.0:
        # For degenerate p, fall back to absolute tolerance.
        return abs(observed - expected) < 0.01
    sigma = math.sqrt(expected * (1.0 - expected) / n)
    return abs(observed - expected) <= 5.0 * sigma


def _collect_choices(exp: "ga.Experiment", n: int, address_filter):
    """Run `n` seeds, yield every trace choice whose address passes
    `address_filter(addr)`."""
    compiled = exp.compile()
    for seed in range(n):
        outcome = compiled.simulator.run(seed=seed)
        for choice in outcome.trace().choices():
            if address_filter(choice.address):
                yield choice


def _int_value(choice_value) -> int:
    """Extract the int from a ChoiceValue::Int(N)-formatted string."""
    return int(str(choice_value).strip("Int()"))


def _bool_value(choice_value) -> bool:
    return "true" in str(choice_value).lower()


# ──────────────────────────────────────────────────────────────────
# D1: NP length narrowing under productive (uniform natural)
# ──────────────────────────────────────────────────────────────────


def test_np_length_narrowing_under_productive_is_uniform_on_admissible() -> None:
    """Audit §5 D1: with NP1 length distribution uniform over
    {0..5} and productive_only active, the admissible subset is
    {0, 3} (the lengths that make junction divisible by 3).
    Expected: Pr[0] ≈ Pr[3] ≈ 0.5 — natural weights are equal, so
    renormalization preserves equality."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(L, 1.0) for L in range(6)])
        .trim(enabled=False)
        .productive_only()
    )
    counts = Counter()
    for seed in range(N_SEEDS):
        rec = exp.run_records(n=1, seed=seed).records[0]
        counts[rec["np1_length"]] += 1

    # Only admissible values should appear.
    assert set(counts.keys()) == {0, 3}, (
        f"non-admissible NP1 lengths leaked: {set(counts.keys())}"
    )

    # Empirical proportions within 5σ of 0.5.
    p0 = counts[0] / N_SEEDS
    p3 = counts[3] / N_SEEDS
    assert _within_5sigma(p0, 0.5, N_SEEDS), (
        f"Pr[NP1=0] = {p0:.4f}, expected 0.5 ± 5σ"
    )
    assert _within_5sigma(p3, 0.5, N_SEEDS), (
        f"Pr[NP1=3] = {p3:.4f}, expected 0.5 ± 5σ"
    )


# ──────────────────────────────────────────────────────────────────
# D2: NP length renormalization (non-uniform natural)
#
# Negative control: if the engine were doing "filter then sample
# uniformly from survivors," this test would observe ≈ 50/50 for the
# {0, 3} subset. The correct natural-weight-preserving renormalization
# yields ≈ 67/33.
# ──────────────────────────────────────────────────────────────────


def test_np_length_renormalization_preserves_natural_weights() -> None:
    """Audit §5 D2: natural distribution {0: weight=2.0, 3: weight=1.0}.
    Both lengths are admissible under productive_only. Renormalization
    must preserve the 2:1 ratio → Pr[0] = 2/3, Pr[3] = 1/3.

    This is the audit's primary **negative control**: a buggy
    implementation that filters-then-samples-uniformly would yield
    Pr[0] = Pr[3] = 0.5, distinguishable from 0.667 at high
    confidence with N=4000."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(0, 2.0), (3, 1.0)])
        .trim(enabled=False)
        .productive_only()
    )
    counts = Counter()
    for seed in range(N_SEEDS):
        rec = exp.run_records(n=1, seed=seed).records[0]
        counts[rec["np1_length"]] += 1

    assert set(counts.keys()) == {0, 3}
    p0 = counts[0] / N_SEEDS
    p3 = counts[3] / N_SEEDS

    # Expected: 2/3 vs 1/3.
    assert _within_5sigma(p0, 2.0 / 3.0, N_SEEDS), (
        f"Pr[NP1=0] = {p0:.4f}, expected 0.6667 ± 5σ "
        f"(natural-weight-preserving renormalization)"
    )
    assert _within_5sigma(p3, 1.0 / 3.0, N_SEEDS), (
        f"Pr[NP1=3] = {p3:.4f}, expected 0.3333 ± 5σ"
    )

    # Negative-control bound: a buggy "filter then uniform"
    # implementation would put p0 near 0.5. Confirm we're NOT there.
    assert abs(p0 - 0.5) > 0.05, (
        f"Pr[NP1=0] = {p0:.4f} is suspiciously close to the buggy 0.5; "
        f"renormalization may be missing"
    )


# ──────────────────────────────────────────────────────────────────
# D3: PCR site distribution is uniform
# ──────────────────────────────────────────────────────────────────


def test_pcr_site_distribution_is_uniform() -> None:
    """Audit §5 D3: `pcr_amplify` samples sites uniformly from
    `[0, pool_len)`. Without active contracts narrowing the per-site
    mask, every position should be drawn with equal probability.

    Pool length is 18 (baseline VJ); expected per-site count ≈
    total / 18."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .pcr_amplify(rate=1.0 / 18.0)
    )
    site_counts = Counter()
    for choice in _collect_choices(
        exp, N_SEEDS, lambda a: "pcr.error_site" in a
    ):
        site_counts[_int_value(choice.value)] += 1

    total = sum(site_counts.values())
    assert total >= 1000, f"too few PCR errors ({total}) for statistical power"

    # Every position 0..17 should appear with frequency ≈ 1/18.
    p_per_site = 1.0 / 18.0
    for site in range(18):
        observed = site_counts.get(site, 0) / total
        assert _within_5sigma(observed, p_per_site, total), (
            f"site {site}: observed Pr = {observed:.4f}, "
            f"expected {p_per_site:.4f} ± 5σ (N={total})"
        )


# ──────────────────────────────────────────────────────────────────
# D4: Indel count=1 kind distribution under productive
#
# The Rust oracle test
# `indel_pass_count_1_under_productive_matches_exact_enumeration`
# uses the same fixture and asserts the same theoretical ratio.
# This Python reproduction confirms the Rust test's invariant
# from the user-facing API surface.
# ──────────────────────────────────────────────────────────────────


def test_indel_count_1_kind_distribution_matches_rust_oracle() -> None:
    """Audit §5 D4: under productive_only + count=1 +
    insertion_prob=0.5 on the baseline VJ fixture, the empirical
    Pr[kind = insertion] should match the Rust oracle theoretical
    value of ≈ 0.5517 within 5σ.

    The oracle value is derived in the Rust test's comment from
    exact enumeration of admissible (kind, site, base) tuples
    weighted by per-event natural mass; the asymmetry from 0.5
    reflects that insertion has one more admissible site than
    deletion (end-of-pool insertion exists; end-of-pool deletion
    doesn't)."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .polymerase_indels(count=1, insertion_prob=0.5)
        .productive_only()
    )
    insertions = 0
    deletions = 0
    for choice in _collect_choices(
        exp, N_SEEDS, lambda a: a == "corrupt.indel.kind[0]"
    ):
        if _bool_value(choice.value):
            insertions += 1
        else:
            deletions += 1

    total = insertions + deletions
    assert total == N_SEEDS, f"expected {N_SEEDS} kind records, got {total}"

    p_ins = insertions / total
    expected = 0.5517  # Rust oracle
    assert _within_5sigma(p_ins, expected, total), (
        f"Pr[insertion] = {p_ins:.4f}, expected {expected:.4f} ± 5σ "
        f"(Rust oracle from indel_pass_count_1_under_productive_matches_exact_enumeration)"
    )


# ──────────────────────────────────────────────────────────────────
# D5: Indel count=2 continuation weighting (the gold-standard
# negative control)
#
# Without continuation weighting (the DP multiplier), the empirical
# Pr[step_0_delta = +1] would converge to ≈ 0.2222 instead of the
# correct ≈ 0.1353. The 5σ tolerance at N=4000 is tight enough that
# breaking the DP multiplier WOULD fail this test.
# ──────────────────────────────────────────────────────────────────


def test_indel_count_2_continuation_weighting_matches_rust_oracle() -> None:
    """Audit §5 D5: under productive_only + count=2 +
    insertion_prob=0.5, the empirical Pr[step-0 event is an
    insertion in V/NP1/D/NP2 (FrameDelta(+1) site)] should match
    the Rust oracle theoretical value of ≈ 0.1353 within 5σ.

    **Negative control**: without continuation weighting (DP
    multiplier in `sample_admissible_tuple`), the per-step
    distribution would be unconstrained-natural ≈ 0.2222 — a
    distinguishing 4σ gap at N=4000. If the DP is ever broken,
    this test fails before the Rust oracle does.

    See ``engine_rs/src/passes/corrupt/indel/tests/constraints.rs``
    for the theoretical derivation."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .polymerase_indels(count=2, insertion_prob=0.5)
        .productive_only()
    )

    # Step 0 of each 2-event tuple: classify by (kind, site).
    # FrameDelta(+1) = insertion at site in V/NP region (pool position
    # < J_region.start = 12 on baseline fixture).
    plus_one_count = 0
    total = 0
    compiled = exp.compile()
    for seed in range(N_SEEDS):
        outcome = compiled.simulator.run(seed=seed)
        choices = outcome.trace().choices()
        kind0 = None
        site0 = None
        for choice in choices:
            if choice.address == "corrupt.indel.kind[0]":
                kind0 = _bool_value(choice.value)
            elif choice.address == "corrupt.indel.site[0]":
                site0 = _int_value(choice.value)
            if kind0 is not None and site0 is not None:
                break
        if kind0 is None or site0 is None:
            continue
        total += 1
        # FrameDelta(+1): insertion at site in [V_region.start=0, J_region.start=12).
        if kind0 and 0 <= site0 < 12:
            plus_one_count += 1

    assert total >= N_SEEDS - 50, f"too many missing step-0 records: {total}"
    p_plus_one = plus_one_count / total
    expected = 0.1353  # Rust oracle (continuation-weighted)

    assert _within_5sigma(p_plus_one, expected, total), (
        f"Pr[step_0_delta=+1] = {p_plus_one:.4f}, "
        f"expected {expected:.4f} ± 5σ "
        f"(continuation-weighted; without DP would be ≈ 0.2222)"
    )

    # Negative-control bound: if the DP were missing, p_plus_one
    # would be near 0.2222 (the unconstrained natural marginal).
    # Confirm we're NOT there.
    assert abs(p_plus_one - 0.2222) > 0.03, (
        f"Pr[step_0_delta=+1] = {p_plus_one:.4f} is suspiciously close "
        f"to the unconstrained-natural 0.2222; the DP continuation "
        f"weighting may be missing"
    )


# ──────────────────────────────────────────────────────────────────
# D6: SHM under productive narrows V_anchor_codon substitutions
# ──────────────────────────────────────────────────────────────────


def test_shm_under_productive_narrows_v_anchor_codon_substitutions() -> None:
    """Audit §5 D6: substitutions inside the V anchor codon
    (pool[6:9] on baseline) must preserve the Cys amino acid under
    productive_only. The per-site mask narrows to synonymous codons
    only, so the empirical mutation rate at sites 6/7/8 under
    productive is **measurably lower** than under unconstrained
    (where any substitution is allowed).

    Qualitative test — exact ratio depends on S5F mutability at
    each site."""
    # IMPORTANT: build two separate Experiment instances. The DSL's
    # `productive_only()` returns `self` (in-place mutation), so
    # sharing a base would entangle the two query branches.
    def _build_exp(productive: bool) -> "ga.Experiment":
        exp = (
            ga.Experiment.on(_vj_basic())
            .recombine(np1_lengths=[(3, 1.0)])
            .trim(enabled=False)
            .mutate(rate=1.0 / 18.0)
        )
        if productive:
            exp = exp.productive_only()
        return exp

    exp_unconstrained = _build_exp(productive=False)
    exp_productive = _build_exp(productive=True)

    def count_anchor_site_hits(exp):
        anchor_hits = 0
        compiled = exp.compile()
        for seed in range(N_SEEDS):
            outcome = compiled.simulator.run(seed=seed)
            for choice in outcome.trace().choices():
                if choice.address.startswith("mutate.s5f.site"):
                    site = _int_value(choice.value)
                    if 6 <= site <= 8:  # V anchor codon positions
                        anchor_hits += 1
        return anchor_hits

    uncon_hits = count_anchor_site_hits(exp_unconstrained)
    prod_hits = count_anchor_site_hits(exp_productive)

    # Productive should hit V_anchor_codon strictly less often than
    # unconstrained — the mask narrows admissible substitutions, so
    # the integrated mass at these sites drops.
    assert prod_hits < uncon_hits, (
        f"productive_only didn't narrow V_anchor_codon substitutions: "
        f"prod_hits={prod_hits} >= uncon_hits={uncon_hits}"
    )


# ──────────────────────────────────────────────────────────────────
# D7: PCR base distribution is uniform
# ──────────────────────────────────────────────────────────────────


def test_pcr_substitution_base_distribution_is_uniform() -> None:
    """Audit §5 D7: PCR error base sampling is uniform over A/C/G/T
    when no contracts narrow the per-site mask. With pool length 18
    and rate=1/18 → ~1 error per record, accumulated over 4000
    seeds, each base should appear with frequency ≈ 0.25."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .pcr_amplify(rate=1.0 / 18.0)
    )
    base_counts = Counter()
    for choice in _collect_choices(
        exp, N_SEEDS, lambda a: "pcr.error_base" in a
    ):
        # ChoiceValue::Base(N) stringifies as "Base(N)"; extract the
        # byte character. The repr varies — accept any single-char
        # canonical base inside the string.
        s = str(choice.value)
        for b in "ACGT":
            if b in s:
                base_counts[b] += 1
                break

    total = sum(base_counts.values())
    assert total >= 1000, f"too few PCR errors ({total})"
    assert set(base_counts.keys()) == set("ACGT")

    for b in "ACGT":
        observed = base_counts[b] / total
        assert _within_5sigma(observed, 0.25, total), (
            f"base {b}: observed Pr = {observed:.4f}, "
            f"expected 0.25 ± 5σ (N={total})"
        )
