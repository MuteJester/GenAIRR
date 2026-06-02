"""Contract pins for the Experiment DSL / Plan Signature
Completeness audit.

Companion to
[`docs/plan_signature_completeness_audit.md`](../docs/plan_signature_completeness_audit.md).

The audit's central claim is that **every parameterized
public DSL surface that changes proposal support or output
folds into `pass_plan_signature`**, with three documented
soft gaps (SampleAllele distribution narrowing, receptor-
revision replacement V distribution, S5F kernel payload),
each backed by an independent safety mechanism.

Pin set:

- ``pin_change_*`` — per-parameter-family change-detection:
  two experiments differing only by the named parameter
  produce DIFFERENT signatures.
- ``pin_default_*`` — equivalent-default normalization:
  no-kwarg vs explicit canonical-default produce IDENTICAL
  signatures.
- ``pin_canonical_*`` — dict-order insensitivity for the
  canonicalizers (`segment_rates`, `v_subregion_rates`,
  Markov rows, empirical distribution pairs).
- ``pin_replay_gate_*`` — replay across a mismatched
  signature fails the gate BEFORE consuming any trace.
- ``pin_soft_gap_*`` — the three documented soft gaps +
  their backstops are pinned so a regression that silently
  changes either layer surfaces here.
- ``pin_legacy_*`` — v1/v2 trace compatibility surfaces
  remain stable.
- ``pin_unfolded_*`` — every production Pass with non-default
  configurable parameters implements a non-empty
  ``parameter_signature`` when those parameters are
  configured.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "audit-docs" / "plan_signature_completeness_audit.md"


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────


def _sig(exp: "ga.Experiment", *, seed: int = 0) -> str:
    """Compile, run one record, and return the resulting
    plan signature. Cheap proxy for "what signature does this
    experiment compile to" — the signature is a function of
    the plan, not the sampled state, so seed only affects
    which choices are recorded (irrelevant to the test)."""
    compiled = exp.compile()
    outcome = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(outcome, seed=seed)
    return json.loads(tf.to_json())["pass_plan_signature"]


def _cfg_with(np_lengths=None, np_bases=None, p_lengths=None, trims=None):
    """Deep-copy the bundled HUMAN_IGH_OGRDB and attach the
    requested typed reference-models plane(s)."""
    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    existing = getattr(cfg, "reference_models", None)
    kwargs = {}
    if isinstance(existing, ReferenceEmpiricalModels):
        kwargs.update(
            np_lengths=existing.np_lengths,
            trims=existing.trims,
            np_bases=existing.np_bases,
            p_nucleotide_lengths=existing.p_nucleotide_lengths,
        )
    if np_lengths is not None:
        kwargs["np_lengths"] = np_lengths
    if np_bases is not None:
        kwargs["np_bases"] = np_bases
    if p_lengths is not None:
        kwargs["p_nucleotide_lengths"] = p_lengths
    if trims is not None:
        kwargs["trims"] = trims
    cfg.reference_models = ReferenceEmpiricalModels(**kwargs)
    return cfg


# ──────────────────────────────────────────────────────────────────
# 1. pin_change_* — per-parameter-family change detection
# ──────────────────────────────────────────────────────────────────


def test_pin_change_np1_length_distribution_changes_signature() -> None:
    cfg_a = _cfg_with(np_lengths={"NP1": EmpiricalDistributionSpec([(0, 1.0)])})
    cfg_b = _cfg_with(np_lengths={"NP1": EmpiricalDistributionSpec([(1, 1.0)])})
    assert _sig(ga.Experiment.on(cfg_a).recombine()) != _sig(
        ga.Experiment.on(cfg_b).recombine()
    )


def test_pin_change_p_v3_length_distribution_changes_signature() -> None:
    cfg_a = _cfg_with(
        p_lengths={"V_3": EmpiricalDistributionSpec([(0, 1.0), (1, 1.0)])}
    )
    cfg_b = _cfg_with(
        p_lengths={"V_3": EmpiricalDistributionSpec([(0, 1.0), (2, 1.0)])}
    )
    assert _sig(ga.Experiment.on(cfg_a).recombine()) != _sig(
        ga.Experiment.on(cfg_b).recombine()
    )


def test_pin_change_trim_v_3_distribution_changes_signature() -> None:
    cfg_a = _cfg_with(trims={"V_3": EmpiricalDistributionSpec([(0, 1.0)])})
    cfg_b = _cfg_with(trims={"V_3": EmpiricalDistributionSpec([(0, 1.0), (1, 1.0)])})
    assert _sig(ga.Experiment.on(cfg_a).recombine()) != _sig(
        ga.Experiment.on(cfg_b).recombine()
    )


def test_pin_change_np_base_model_kind_changes_signature() -> None:
    """Uniform → empirical_first_base → markov are three
    different generator paths; each must produce a distinct
    plan signature."""
    cfg_uniform = _cfg_with(np_bases={"NP1": NpBaseModelSpec(kind="uniform")})
    cfg_emp = _cfg_with(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="empirical_first_base",
                first_base={"A": 4.0, "C": 1.0, "G": 1.0, "T": 1.0},
            )
        }
    )
    cfg_markov = _cfg_with(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                transitions={
                    b: {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0} for b in "ACGT"
                },
            )
        }
    )
    sigs = {
        _sig(ga.Experiment.on(c).recombine())
        for c in (cfg_uniform, cfg_emp, cfg_markov)
    }
    assert len(sigs) == 3, "three np_base model kinds collapsed to fewer signatures"


def test_pin_change_markov_transition_row_changes_signature() -> None:
    """Two cartridges with identical first_base but differing
    A-row transitions must produce different plan signatures
    (the Markov v1 boundary)."""
    base = {b: {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0} for b in "ACGT"}
    cfg_a = _cfg_with(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                transitions=base,
            )
        }
    )
    perturbed = copy.deepcopy(base)
    perturbed["A"]["T"] = 99.0  # only A's row diverges
    cfg_b = _cfg_with(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                transitions=perturbed,
            )
        }
    )
    assert _sig(ga.Experiment.on(cfg_a).recombine()) != _sig(
        ga.Experiment.on(cfg_b).recombine()
    )


def test_pin_change_mutation_rate_vs_count_form_changes_signature() -> None:
    """`mutate(rate=...)` lowers to `CountSource::Rate` and
    `mutate(count_pairs=...)` lowers to
    `CountSource::Distribution`. The two forms must produce
    different signature substrings (`count=rate:...` vs
    `count=dist:...`)."""
    a = ga.Experiment.on("human_igh").recombine().mutate(rate=0.01)
    b = ga.Experiment.on("human_igh").recombine().mutate(count=[(1, 1.0)])
    assert _sig(a) != _sig(b)


def test_pin_change_segment_rates_non_default_changes_signature() -> None:
    base = ga.Experiment.on("human_igh").recombine().mutate(rate=0.01)
    perturbed = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(rate=0.01, segment_rates={"V": 2.0, "D": 1.0, "J": 1.0, "NP": 1.0})
    )
    assert _sig(base) != _sig(perturbed)


def test_pin_change_v_subregion_rates_non_default_changes_signature() -> None:
    base = ga.Experiment.on("human_igh").recombine().mutate(rate=0.01)
    perturbed = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            rate=0.01,
            v_subregion_rates={"FWR1": 2.0, "CDR1": 1.0, "FWR2": 1.0, "CDR2": 1.0, "FWR3": 1.0},
        )
    )
    assert _sig(base) != _sig(perturbed)


def test_pin_change_polymerase_indels_insertion_prob_changes_signature() -> None:
    a = (
        ga.Experiment.on("human_igh")
        .recombine()
        .polymerase_indels(count=1, insertion_prob=0.1)
    )
    b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .polymerase_indels(count=1, insertion_prob=0.9)
    )
    assert _sig(a) != _sig(b)


def test_pin_change_end_loss_length_changes_signature() -> None:
    a = ga.Experiment.on("human_igh").recombine().end_loss_5prime(length=(1, 3))
    b = ga.Experiment.on("human_igh").recombine().end_loss_5prime(length=(5, 7))
    assert _sig(a) != _sig(b)


def test_pin_change_contaminate_apply_prob_changes_signature() -> None:
    a = ga.Experiment.on("human_igh").recombine().contaminate(prob=0.1)
    b = ga.Experiment.on("human_igh").recombine().contaminate(prob=0.5)
    assert _sig(a) != _sig(b)


def test_pin_change_random_strand_orientation_prob_changes_signature() -> None:
    a = ga.Experiment.on("human_igh").recombine().random_strand_orientation(prob=0.1)
    b = ga.Experiment.on("human_igh").recombine().random_strand_orientation(prob=0.5)
    assert _sig(a) != _sig(b)


def test_pin_change_invert_d_prob_changes_signature() -> None:
    a = ga.Experiment.on("human_igh").recombine().invert_d(prob=0.1)
    b = ga.Experiment.on("human_igh").recombine().invert_d(prob=0.5)
    assert _sig(a) != _sig(b)


def test_pin_change_receptor_revision_prob_changes_signature() -> None:
    a = ga.Experiment.on("human_igh").recombine().receptor_revision(prob=0.1)
    b = ga.Experiment.on("human_igh").recombine().receptor_revision(prob=0.5)
    assert _sig(a) != _sig(b)


def test_pin_change_paired_end_r1_length_changes_signature() -> None:
    a = ga.Experiment.on("human_igh").recombine().paired_end(r1_length=80, insert_size=200)
    b = ga.Experiment.on("human_igh").recombine().paired_end(r1_length=120, insert_size=200)
    assert _sig(a) != _sig(b)


def test_pin_change_paired_end_insert_size_changes_signature() -> None:
    a = ga.Experiment.on("human_igh").recombine().paired_end(r1_length=80, insert_size=200)
    b = ga.Experiment.on("human_igh").recombine().paired_end(r1_length=80, insert_size=400)
    assert _sig(a) != _sig(b)


# ──────────────────────────────────────────────────────────────────
# 2. pin_default_* — equivalent-default normalization
# ──────────────────────────────────────────────────────────────────


def test_pin_default_segment_rates_default_collides_with_explicit_all_ones() -> None:
    a = ga.Experiment.on("human_igh").recombine().mutate(rate=0.01)
    b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(rate=0.01, segment_rates={"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0})
    )
    assert _sig(a) == _sig(b)


def test_pin_default_v_subregion_rates_default_collides_with_explicit_all_ones() -> None:
    a = ga.Experiment.on("human_igh").recombine().mutate(rate=0.01)
    b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            rate=0.01,
            v_subregion_rates={"FWR1": 1.0, "CDR1": 1.0, "FWR2": 1.0, "CDR2": 1.0, "FWR3": 1.0},
        )
    )
    assert _sig(a) == _sig(b)


def test_pin_default_empty_p_nucleotide_plane_is_byte_identical_to_no_plane() -> None:
    """Empty `p_nucleotide_lengths={}` lowers identically to
    omitting the field entirely (the lowering layer omits
    every `push_p_addition` call when the dict is empty)."""
    cfg_a = copy.deepcopy(ga.HUMAN_IGH_OGRDB)  # no reference_models touch
    cfg_b = _cfg_with(p_lengths={})  # empty plane
    assert _sig(ga.Experiment.on(cfg_a).recombine()) == _sig(
        ga.Experiment.on(cfg_b).recombine()
    )


# ──────────────────────────────────────────────────────────────────
# 3. pin_canonical_* — dict-order insensitivity
# ──────────────────────────────────────────────────────────────────


def test_pin_canonical_segment_rates_dict_order_insensitive() -> None:
    a = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(rate=0.01, segment_rates={"V": 2.0, "D": 0.5, "J": 1.0, "NP": 1.0})
    )
    b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(rate=0.01, segment_rates={"NP": 1.0, "D": 0.5, "V": 2.0, "J": 1.0})
    )
    assert _sig(a) == _sig(b)


def test_pin_canonical_v_subregion_rates_dict_order_insensitive() -> None:
    a = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            rate=0.01,
            v_subregion_rates={"FWR1": 2.0, "CDR1": 0.5, "FWR2": 1.0, "CDR2": 1.0, "FWR3": 1.0},
        )
    )
    b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            rate=0.01,
            v_subregion_rates={"CDR2": 1.0, "FWR3": 1.0, "FWR1": 2.0, "FWR2": 1.0, "CDR1": 0.5},
        )
    )
    assert _sig(a) == _sig(b)


def test_pin_canonical_markov_dict_order_insensitive() -> None:
    """Two cartridges with semantically identical Markov
    matrices but different Python dict insertion order must
    produce equal plan signatures. The canonical A/C/G/T
    walk in `MarkovBaseGenerator::signature()` enforces this."""
    transitions_a = {
        "A": {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
        "C": {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
        "G": {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
        "T": {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
    }
    transitions_b = {
        "T": {"T": 1.0, "G": 1.0, "C": 1.0, "A": 1.0},
        "G": {"T": 1.0, "G": 1.0, "C": 1.0, "A": 1.0},
        "C": {"T": 1.0, "G": 1.0, "C": 1.0, "A": 1.0},
        "A": {"T": 1.0, "G": 1.0, "C": 1.0, "A": 1.0},
    }
    cfg_a = _cfg_with(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                transitions=transitions_a,
            )
        }
    )
    cfg_b = _cfg_with(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                transitions=transitions_b,
            )
        }
    )
    assert _sig(ga.Experiment.on(cfg_a).recombine()) == _sig(
        ga.Experiment.on(cfg_b).recombine()
    )


def test_pin_canonical_empirical_distribution_pair_insertion_order_insensitive() -> None:
    """Distinct `(value, weight)` pairs lower through
    `_np_lengths_from_models` (or the trim / P-length
    equivalents), which sorts by value. Two cartridges with
    the same multiset of distinct pairs but different
    Python-list insertion orders MUST produce equal plan
    signatures.

    **Note:** `EmpiricalLengthDist::from_pairs` does NOT
    deduplicate — `[(1, 0.5), (1, 0.5)]` is kept as two
    entries and produces a DIFFERENT signature from
    `[(1, 1.0)]`. The bridge sort handles input-order
    insensitivity but not multiplicity collapsing. See
    `docs/plan_signature_completeness_audit.md` §4 for the
    boundary."""
    cfg_a = _cfg_with(np_lengths={"NP1": EmpiricalDistributionSpec([(2, 0.3), (1, 0.7)])})
    cfg_b = _cfg_with(np_lengths={"NP1": EmpiricalDistributionSpec([(1, 0.7), (2, 0.3)])})
    assert _sig(ga.Experiment.on(cfg_a).recombine()) == _sig(
        ga.Experiment.on(cfg_b).recombine()
    )


# ──────────────────────────────────────────────────────────────────
# 4. pin_replay_gate_* — mismatched signature fails the gate
# ──────────────────────────────────────────────────────────────────


def _trace_file(exp: "ga.Experiment", *, seed: int) -> "ga._engine.TraceFile":
    compiled = exp.compile()
    return compiled.simulator.trace_file_from(
        compiled.simulator.run(seed=seed), seed=seed
    )


def test_pin_replay_gate_mismatched_np_length_fails_before_consuming_trace() -> None:
    cfg_a = _cfg_with(np_lengths={"NP1": EmpiricalDistributionSpec([(0, 1.0)])})
    cfg_b = _cfg_with(np_lengths={"NP1": EmpiricalDistributionSpec([(1, 1.0)])})
    tf = _trace_file(ga.Experiment.on(cfg_a).recombine(), seed=0)
    cb = ga.Experiment.on(cfg_b).recombine().compile()
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        cb.simulator.replay_from_trace_file(tf)


def test_pin_replay_gate_mismatched_markov_matrix_fails_before_consuming_trace() -> None:
    base = {b: {"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0} for b in "ACGT"}
    perturbed = copy.deepcopy(base)
    perturbed["A"]["T"] = 99.0
    cfg_a = _cfg_with(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                transitions=base,
            )
        }
    )
    cfg_b = _cfg_with(
        np_bases={
            "NP1": NpBaseModelSpec(
                kind="markov",
                first_base={"A": 1.0, "C": 1.0, "G": 1.0, "T": 1.0},
                transitions=perturbed,
            )
        }
    )
    tf = _trace_file(ga.Experiment.on(cfg_a).recombine(), seed=0)
    cb = ga.Experiment.on(cfg_b).recombine().compile()
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        cb.simulator.replay_from_trace_file(tf)


def test_pin_replay_gate_mismatched_p_v3_length_fails_before_consuming_trace() -> None:
    cfg_a = _cfg_with(p_lengths={"V_3": EmpiricalDistributionSpec([(1, 1.0)])})
    cfg_b = _cfg_with(p_lengths={"V_3": EmpiricalDistributionSpec([(2, 1.0)])})
    tf = _trace_file(ga.Experiment.on(cfg_a).recombine(), seed=0)
    cb = ga.Experiment.on(cfg_b).recombine().compile()
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        cb.simulator.replay_from_trace_file(tf)


def test_pin_replay_gate_mismatched_segment_rates_fails_before_consuming_trace() -> None:
    exp_a = ga.Experiment.on("human_igh").recombine().mutate(rate=0.01)
    exp_b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(rate=0.01, segment_rates={"V": 2.0, "D": 1.0, "J": 1.0, "NP": 1.0})
    )
    tf = _trace_file(exp_a, seed=0)
    cb = exp_b.compile()
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        cb.simulator.replay_from_trace_file(tf)


def test_pin_replay_gate_mismatched_paired_end_insert_fails_before_consuming_trace() -> None:
    exp_a = (
        ga.Experiment.on("human_igh")
        .recombine()
        .paired_end(r1_length=80, insert_size=200)
    )
    exp_b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .paired_end(r1_length=80, insert_size=400)
    )
    tf = _trace_file(exp_a, seed=0)
    cb = exp_b.compile()
    with pytest.raises(ValueError, match="pass plan signature mismatch"):
        cb.simulator.replay_from_trace_file(tf)


# ──────────────────────────────────────────────────────────────────
# 5. pin_soft_gap_* — documented gaps + their backstops
# ──────────────────────────────────────────────────────────────────


def test_pin_soft_gap_restrict_alleles_does_not_change_plan_signature() -> None:
    """Documented soft gap (1) — `restrict_alleles` narrows
    the SampleAllele distribution but is NOT folded into the
    plan signature. Pinned so a future "fix" surfaces here
    explicitly (the fix would also need to update the audit
    doc + tighten the backstop tests below)."""
    a = ga.Experiment.on("human_igh").recombine()
    b = ga.Experiment.on("human_igh").restrict_alleles(v="IGHVF1-G1*01").recombine()
    assert _sig(a) == _sig(b), (
        "restrict_alleles now changes plan signature — the documented "
        "soft gap (1) has been tightened. Update "
        "docs/plan_signature_completeness_audit.md §8 and §10 "
        "accordingly."
    )


def test_pin_soft_gap_allele_weights_does_not_change_plan_signature() -> None:
    """Documented soft gap (1) extension — `v_allele_weights`
    reweights the SampleAllele distribution but is NOT
    folded."""
    a = ga.Experiment.on("human_igh").recombine()
    b = (
        ga.Experiment.on("human_igh")
        .recombine(v_allele_weights={"IGHVF1-G1*01": 100.0})
    )
    assert _sig(a) == _sig(b), (
        "v_allele_weights now changes plan signature — the documented "
        "soft gap has been tightened. Update audit doc."
    )


def test_pin_soft_gap_sample_allele_strict_backstop_rejects_out_of_support() -> None:
    """The backstop for soft gap (1): when the recorded
    allele ID is outside the narrowed distribution's support,
    the strict sampler rejects rather than silently substituting.

    Construct two single-allele restrictions to DIFFERENT V
    alleles. The trace recorded under restriction-to-A
    carries A's allele ID; replaying it against
    restriction-to-B MUST fail at the SampleAllele sampler
    because A is outside B's support. Plan signatures still
    collide because the distribution narrowing isn't folded
    (soft gap 1)."""
    # Pick two genuinely-different V alleles (from different
    # gene families so their IDs are guaranteed to differ).
    v_alleles_sorted = sorted(
        a.name for gene in ga.HUMAN_IGH_OGRDB.v_alleles.values() for a in gene
    )
    allele_a = v_alleles_sorted[0]
    allele_b = v_alleles_sorted[len(v_alleles_sorted) // 2]
    assert allele_a != allele_b
    exp_a = ga.Experiment.on("human_igh").restrict_alleles(v=allele_a).recombine()
    exp_b = ga.Experiment.on("human_igh").restrict_alleles(v=allele_b).recombine()
    # Plan signatures collide (soft gap 1).
    assert _sig(exp_a) == _sig(exp_b)
    # But strict replay against the OTHER restriction fails
    # at the sampler — the recorded allele ID isn't in B's
    # support of size 1.
    ca = exp_a.compile()
    tf = ca.simulator.trace_file_from(ca.simulator.run(seed=0), seed=0)
    cb = exp_b.compile()
    with pytest.raises(
        Exception, match="sample_allele|distribution_support|StrictSamplingError"
    ):
        cb.simulator.replay_from_trace_file(tf)


def test_pin_soft_gap_s5f_kernel_digest_in_cartridge_manifest() -> None:
    """The backstop for soft gap (3): the cartridge manifest
    surfaces `s5f_kernel_digest` independently of the plan
    signature. A regression that removes the manifest digest
    would surface here."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    shm = m["models"]["shm"]
    assert "s5f_kernel_digest" in shm, (
        "manifest no longer exposes s5f_kernel_digest — soft gap (3)'s "
        "backstop is broken. Restore the digest OR fold the kernel "
        "payload into the plan signature in lockstep."
    )
    assert shm["in_content_hash"] is False, (
        "shm.in_content_hash flipped to True — verify the kernel "
        "choice is now folded into refdata_content_hash AND the "
        "audit doc soft gap (3) is updated to reflect the tightening"
    )


# ──────────────────────────────────────────────────────────────────
# 6. pin_legacy_* — v1/v2 trace compatibility
# ──────────────────────────────────────────────────────────────────


def test_pin_legacy_address_schema_version_is_one() -> None:
    """The on-disk trace-address vocabulary is versioned by
    `ADDRESS_SCHEMA_VERSION` in
    `engine_rs/src/address.rs`. The slice-by-slice additive
    discipline (new variants only, no rename / drop) means
    the version stays at 1 across the recent biology slices
    (P-nucleotide, Markov, V-subregion rates / counters,
    paired-end, etc.). A regression that bumps it without
    updating the v1 / v2 fixture set would surface here.
    """
    src = (_REPO_ROOT / "engine_rs" / "src" / "address.rs").read_text(encoding="utf-8")
    assert "pub const ADDRESS_SCHEMA_VERSION: u32 = 1;" in src, (
        "address schema version bumped — refresh the v1/v2 trace "
        "compat fixtures and the address-schema documentation in "
        "lockstep"
    )


def test_pin_legacy_frozen_address_spellings_test_exists() -> None:
    """The `frozen_address_spellings_for_choice_address_schema_v1`
    Rust unit test in `address.rs` pins one representative
    of every typed variant. The Python contract pin here
    asserts the Rust test fixture is present (the Rust
    test itself is the actual enforcement)."""
    src = (_REPO_ROOT / "engine_rs" / "src" / "address.rs").read_text(encoding="utf-8")
    assert "fn frozen_address_spellings_for_choice_address_schema_v1" in src


# ──────────────────────────────────────────────────────────────────
# 7. pin_unfolded_* — every production pass with non-default
# configurable parameters folds them into parameter_signature
# ──────────────────────────────────────────────────────────────────


def test_pin_unfolded_production_passes_define_parameter_signature() -> None:
    """Every production `Pass` implementation (excluding the
    documented test-only `EchoPass` / `SampleBasePass`) with
    at least one non-discriminator constructor parameter
    overrides `parameter_signature`. The default trait impl
    returns `""` — passes that have parameters MUST override
    or their signatures silently regress.

    Pinned by scanning the engine source. New production
    passes with parameters that fall through to the default
    surface here."""
    passes_dir = _REPO_ROOT / "engine_rs" / "src" / "passes"
    # Files that don't need a parameter_signature override:
    # - `assemble_segment.rs` — only takes a `segment`, which
    #   is in `name()`.
    # - `echo.rs` / `sample_base.rs` — test-only passes used
    #   exclusively by `run_smoke_plan`.
    test_only = {"echo.rs", "sample_base.rs"}
    no_params = {"assemble_segment.rs"}
    missing: list[str] = []
    for rs in passes_dir.rglob("*.rs"):
        # Skip submodule organisational files (`mod.rs`,
        # nested `tests/`, etc.) and the test-only / no-
        # params passes.
        if rs.name in test_only or rs.name in no_params:
            continue
        if rs.name == "mod.rs":
            continue
        if "tests" in rs.parts:
            continue
        if rs.parent != passes_dir and rs.parent.name in {
            "generate_np",
            "mutate",
            "corrupt",
            "mutation_transaction",
            "assemble_segment",
        }:
            # These are the per-pass submodule files. Some
            # of them (`generate_np/execution.rs`,
            # `mutate/uniform.rs`, etc.) DO define
            # parameter_signature for the containing pass.
            # Allow only the parent module's file to be the
            # declaration site.
            pass
        src = rs.read_text(encoding="utf-8")
        if "impl Pass for" not in src:
            continue
        if "fn parameter_signature" not in src:
            missing.append(str(rs.relative_to(_REPO_ROOT)))
    assert not missing, (
        f"production Pass implementations missing parameter_signature: "
        f"{missing}. Add a fold (or document the exception in the "
        f"audit doc + this test's allow-list)."
    )


# ──────────────────────────────────────────────────────────────────
# 8. Doc anchor
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_plan_signature_completeness_contract.py" in doc, (
        "audit doc no longer references the contract file"
    )
    for marker in (
        "## 1. Existing scaffolding",
        "## 2. Q1",
        "## 8. The three documented soft gaps",
        "## 12. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
