"""Contract pins for the V-Subregion SHM Rate audit.

Companion to
[`docs/v_subregion_shm_rate_design.md`](../docs/v_subregion_shm_rate_design.md).
The audit is pre-implementation: it freezes today's surfaces
(``pin_scaffold_*``) and the gaps a future implementation slice
would close (``pin_absence_*``).

This audit follows two prior shipped slices:

- **Targeted SHM** (`shm_segment_rate_design.md`) — per-segment
  rate scalars on the V / D / J / NP buckets.
- **V-Subregion Cartridge Annotation Surface**
  (`v_region_substructure_audit.md` Slice 1) — per-V-allele
  IMGT FWR1 / CDR1 / FWR2 / CDR2 / FWR3 intervals as a cartridge
  property.

The next slice connects the two: a
``v_subregion_rates={"CDR1": 2.0, "FWR": 0.5}`` kwarg on
``Experiment.mutate`` that layers per-V-subregion scalars on top
of the existing per-segment rates.

Split:

- ``pin_scaffold_*`` tests freeze the pre-existing surfaces the
  audit's Slice B builds on: the segment-rate classifier
  signature, the S5F profile-build call site, the V-subregion
  cartridge annotations from Slice 1, today's mutate trace
  address vocabulary, and the (currently broken)
  ``pass_plan_signature`` shape that omits compile-time pass
  parameters.
- ``pin_absence_*`` tests freeze the gaps Slice B closes: the
  ``v_subregion_rates`` kwarg, the DSL validator, the Rust rate
  struct, the manifest's rate-support block, and the
  per-region AIRR counter fields (the latter deferred to a
  separate counters slice).

A finding the audit surfaces in §6: the segment-rates slice's
audit promised ``pass_plan_signature`` would fold compile-time
pass parameters into the signature so a mismatched rate vector
fails the replay-safety check. That hasn't shipped — pinned
here as a ``pin_scaffold_*`` so the gap-closure flips it in
lockstep.
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path

import pytest

import GenAIRR as ga


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "audit-docs" / "v_subregion_shm_rate_design.md"


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — segment_at_position signature
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_segment_at_position_takes_only_sequence_and_pos() -> None:
    """Audit §2: the existing ``segment_at_position`` classifier
    takes only ``(sequence, pos)`` and does NOT receive any
    allele identity. Pinned at source so the slice's sibling
    helper (`v_subregion_at_position`) is added with a
    deliberately wider signature: ``(pos, sequence, assignments,
    refdata)`` to reach the assigned V allele's subregions."""
    src = (
        _REPO_ROOT / "engine_rs" / "src" / "passes" / "mutate" / "segment_rates.rs"
    ).read_text(encoding="utf-8")
    import re

    match = re.search(
        r"pub fn segment_at_position\s*\((.*?)\)\s*->\s*Option<Segment>",
        src,
        re.DOTALL,
    )
    assert match is not None, (
        "segment_at_position signature not found; audit §2 reference "
        "drifted"
    )
    params = match.group(1)
    # Today's parameters: sequence + pos. No assignments, no refdata.
    assert "sequence:" in params
    assert "pos:" in params
    for forbidden in ("assignments:", "refdata:", "allele_id:"):
        assert forbidden not in params, (
            f"segment_at_position now takes {forbidden!r}; the slice's "
            "v_subregion_at_position sibling assumption changed — flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — S5F profile-build applies segment-rate factor
# ──────────────────────────────────────────────────────────────────


def test_pin_present_s5f_build_profile_applies_combined_rate_factor() -> None:
    """Audit §2 — Slice B landed. ``S5FMutationPass::build_profile``
    is the single per-position weight walk where both segment
    rates and V-subregion rates apply. The slice unified the two
    multipliers into a single ``combined_site_factor`` helper
    (in ``mutation_transaction::substitution``) so the
    build_profile body has one call site that handles both. Pin
    the unified call structure."""
    src = (
        _REPO_ROOT
        / "engine_rs"
        / "src"
        / "passes"
        / "mutate"
        / "s5f"
        / "sampling.rs"
    ).read_text(encoding="utf-8")
    assert "fn build_profile" in src, (
        "S5F build_profile no longer exists; sampling-loop "
        "architecture regressed"
    )
    # The combined helper (segment × V-subregion factor) replaces
    # the old direct segment_at_position call.
    assert "combined_site_factor" in src, (
        "build_profile no longer calls combined_site_factor; "
        "Slice B unification regressed"
    )
    # Both rate vectors must still flow into the build_profile body.
    assert "segment_rates" in src
    assert "v_subregion_rates" in src


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — Simulation.assignments carries AlleleInstance
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_simulation_assignments_carries_allele_identity() -> None:
    """Audit §2: ``Simulation.assignments: AlleleAssignments`` is
    populated before the mutate pass runs and exposes the V
    allele's ``allele_id`` + ``trim_5``. Pin the struct shape so
    a refactor that hides allele identity from the mutate pass
    surfaces."""
    sim_src = (
        _REPO_ROOT / "engine_rs" / "src" / "ir" / "simulation.rs"
    ).read_text(encoding="utf-8")
    assignment_src = (
        _REPO_ROOT / "engine_rs" / "src" / "assignment.rs"
    ).read_text(encoding="utf-8")
    # Simulation carries the assignments field (the path-prefixed
    # form ``crate::assignment::AlleleAssignments`` is current).
    assert "pub assignments:" in sim_src and "AlleleAssignments" in sim_src, (
        "Simulation.assignments field missing or renamed; "
        "v_subregion_rates feasibility regressed"
    )
    # AlleleInstance exposes allele_id + trim_5.
    assert "pub allele_id: AlleleId" in assignment_src
    assert "pub trim_5: u16" in assignment_src


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — bundled human cartridges have 100% V-subregion coverage
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "name",
    ["HUMAN_IGH_OGRDB", "HUMAN_IGK_OGRDB", "HUMAN_IGL_OGRDB"],
)
def test_pin_scaffold_bundled_v_has_full_subregion_coverage(name: str) -> None:
    """Slice 1 prerequisite. Pinned here too because the
    v_subregion_rates slice is unsatisfiable on cartridges with
    zero annotated V alleles. The bundled human Ig cartridges
    must continue to carry full coverage so the rate kwarg has
    something to weight."""
    cfg = getattr(ga, name)
    sup = cfg.cartridge_manifest()["models"]["shm"]["v_subregion_support"]
    assert sup["available"] is True
    assert sup["annotated_v_count"] == sup["total_v_count"] > 0


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — V alleles carry the five canonical labels
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_v_alleles_carry_five_canonical_labels() -> None:
    """Every V allele in ``HUMAN_IGH_OGRDB`` exposes exactly the
    five canonical IMGT labels through the bridged refdata. The
    rate kwarg's dict accepts the same vocabulary; pinning the
    cartridge surface here keeps the two in lockstep."""
    from GenAIRR._refdata_resolver import dataconfig_to_refdata

    canonical = {"FWR1", "CDR1", "FWR2", "CDR2", "FWR3"}
    rd = dataconfig_to_refdata(ga.HUMAN_IGH_OGRDB)
    sampled = 0
    for v_id in range(rd.v_pool_size()):
        subs = rd.v_allele(v_id).subregions
        if not subs:
            continue
        labels = {label for label, _, _ in subs}
        assert labels == canonical, (
            f"V allele {v_id} has labels {labels}, expected {canonical}; "
            "Slice 1 surface regressed"
        )
        sampled += 1
    assert sampled > 0, "no V allele had subregions; Slice 1 regressed"


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — existing segment_rates still works post-Slice-1
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_segment_rates_still_validates_on_annotated_cartridge() -> None:
    """No regression check: the existing per-segment rate kwarg
    still produces clean simulations when the cartridge carries
    V-subregion annotations. Confirms Slice 1 did not perturb the
    targeted-SHM path the v_subregion_rates slice composes with."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .productive_only()
        .mutate(model="s5f", rate=0.03, segment_rates={"V": 2.0, "NP": 0.0})
    )
    refdata = exp.refdata
    result = exp.run_records(n=20, seed=4242)
    report = result.validate_records(refdata)
    assert report, (
        f"segment_rates regression on annotated cartridge: "
        f"{len(report.failures)}/{report.count} records failed validation"
    )


# ──────────────────────────────────────────────────────────────────
# 7. Scaffold — today's mutate trace address vocabulary
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_mutate_trace_address_vocabulary() -> None:
    """Audit §6: today the mutate pass emits exactly three
    address variants per model — ``Count``, ``Site(u32)``,
    ``Base(u32)``. The v_subregion_rates slice promised no new
    addresses; pin the existing vocabulary at source so a slice
    that adds a new variant surfaces."""
    address_src = (
        _REPO_ROOT / "engine_rs" / "src" / "address.rs"
    ).read_text(encoding="utf-8")
    for expected in (
        "MutateUniformCount",
        "MutateUniformSite(u32)",
        "MutateUniformBase(u32)",
        "MutateS5fCount",
        "MutateS5fSite(u32)",
        "MutateS5fBase(u32)",
    ):
        assert expected in address_src, (
            f"mutate trace address {expected!r} missing; audit §6 "
            "address-surface assumption regressed"
        )
    # And the slice's hypothetical address names must NOT exist.
    for forbidden in (
        "MutateS5fSubregion",
        "MutateUniformSubregion",
        "MutateSubregionWeight",
    ):
        assert forbidden not in address_src, (
            f"mutate trace address {forbidden!r} now exists; slice has "
            "landed — flip pin and update §6"
        )


# ──────────────────────────────────────────────────────────────────
# 8. Scaffold — pass_plan_signature omits compile-time params today
# ──────────────────────────────────────────────────────────────────


def test_pin_present_pass_plan_signature_folds_compile_time_params() -> None:
    """Audit §6 — **Slice A landed.** ``pass_plan_signature`` now
    folds each pass's compile-time parameters into the signature
    via the ``Pass::parameter_signature`` trait method, and the
    on-disk schema version is bumped to v3. A trace recorded
    under one ``segment_rates`` vector now refuses to replay
    against a plan with a different vector (signature gate fires
    before any choice is consumed).

    Positive surface this pin asserts:

      1. The Rust trait carries a ``parameter_signature`` method.
      2. The new format emits ``name(params)`` per pass — the
         signature contains ``(`` characters.
      3. A v1/v2-shape fallback (``pass_plan_signature_names_only``)
         exists for legacy fixture replay.
      4. Two plans differing only by ``segment_rates`` produce
         different plan signatures end-to-end through the Python
         DSL.
      5. Default ``segment_rates`` and an explicit all-ones dict
         produce the **same** plan signature (behaviourally
         equivalent inputs → equal signatures).
    """
    import json

    trace_file_src = (
        _REPO_ROOT / "engine_rs" / "src" / "trace_file.rs"
    ).read_text(encoding="utf-8")
    traits_src = (
        _REPO_ROOT / "engine_rs" / "src" / "pass" / "traits.rs"
    ).read_text(encoding="utf-8")

    # (1) The Pass trait carries the parameter_signature method.
    assert "fn parameter_signature" in traits_src, (
        "Pass trait no longer declares parameter_signature; Slice A "
        "surface regressed"
    )

    # (2) New format emits `name(params)` per pass — the
    # `pass_plan_signature` body uses the `{}({})` envelope.
    import re

    body_match = re.search(
        r"pub fn pass_plan_signature\(plan: &PassPlan\) -> String \{(.*?)\n\}",
        trace_file_src,
        re.DOTALL,
    )
    assert body_match is not None
    body = body_match.group(1)
    assert "parameter_signature" in body, (
        "pass_plan_signature no longer folds parameter_signature; "
        "Slice A surface regressed"
    )

    # (3) v1/v2 fallback function exists.
    assert "pub fn pass_plan_signature_names_only" in trace_file_src, (
        "pass_plan_signature_names_only fallback missing; legacy "
        "fixture replay would break"
    )

    # (4) Two plans differing only by segment_rates produce
    # different signatures end-to-end.
    e_a = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03)
    )
    e_b = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03, segment_rates={"V": 2.0})
    )
    c_a, c_b = e_a.compile(), e_b.compile()
    tf_a = c_a.simulator.trace_file_from(
        c_a.simulator.run(seed=42), seed=42
    )
    tf_b = c_b.simulator.trace_file_from(
        c_b.simulator.run(seed=42), seed=42
    )
    sig_a = json.loads(tf_a.to_json())["pass_plan_signature"]
    sig_b = json.loads(tf_b.to_json())["pass_plan_signature"]
    assert sig_a != sig_b, (
        "plans differing in segment_rates produced equal plan "
        "signatures; Slice A surface regressed"
    )
    # And the v3 envelope `name(...)` is observable.
    assert "(" in sig_a, (
        "plan signature does not use the `name(params)` envelope; "
        "Slice A schema regressed"
    )

    # (5) Default vs explicit all-ones produce equal signatures.
    e_default = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03)
    )
    e_explicit = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(
            model="s5f",
            rate=0.03,
            segment_rates={"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0},
        )
    )
    c_d, c_e = e_default.compile(), e_explicit.compile()
    tf_d = c_d.simulator.trace_file_from(
        c_d.simulator.run(seed=42), seed=42
    )
    tf_e = c_e.simulator.trace_file_from(
        c_e.simulator.run(seed=42), seed=42
    )
    sig_d = json.loads(tf_d.to_json())["pass_plan_signature"]
    sig_e = json.loads(tf_e.to_json())["pass_plan_signature"]
    assert sig_d == sig_e, (
        "default segment_rates and explicit all-ones produced "
        "different signatures; behavioural-equivalence guarantee "
        "regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 9. Absence — no v_subregion_rates kwarg on mutate
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_rates_kwarg_on_mutate() -> None:
    """Audit §1 / §3 — Slice B landed. The kwarg exists with
    default ``None`` (flat-default fast path). The other names
    listed below remain forbidden — only ``v_subregion_rates``
    is the canonical surface."""
    sig = inspect.signature(ga.Experiment.mutate)
    assert "v_subregion_rates" in sig.parameters, (
        "Experiment.mutate is missing the v_subregion_rates kwarg; "
        "Slice B surface regressed"
    )
    assert sig.parameters["v_subregion_rates"].default is None
    for forbidden in (
        "subregion_rates",
        "v_region_rates",
        "cdr_rates",
        "fr_rates",
        "fwr_rates",
    ):
        assert forbidden not in sig.parameters, (
            f"Experiment.mutate now accepts {forbidden!r}; only "
            "``v_subregion_rates`` is the canonical Slice B surface — "
            "additional aliases would require a manifest update"
        )


# ──────────────────────────────────────────────────────────────────
# 10. Absence — no _validate_v_subregion_rates in experiment.py
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_rates_validator_in_experiment() -> None:
    """Audit §3 / §4 — Slice B landed. The DSL-side validator
    ``_validate_v_subregion_rates`` exists and exposes the
    canonical label vocabulary plus the FWR / CDR alias-expansion
    rule. The other (legacy) names listed below remain
    forbidden — only ``_validate_v_subregion_rates`` is the
    canonical surface."""
    from GenAIRR import experiment as exp_module

    assert hasattr(exp_module, "_validate_v_subregion_rates"), (
        "_validate_v_subregion_rates missing from experiment.py; "
        "Slice B DSL boundary regressed"
    )
    # The default tuple (all 1.0) is exposed so other code (manifest,
    # compile-time checks) can identify the flat-default case.
    assert hasattr(exp_module, "_DEFAULT_V_SUBREGION_RATES")
    assert exp_module._DEFAULT_V_SUBREGION_RATES == (1.0, 1.0, 1.0, 1.0, 1.0)
    # Alias and canonical label constants are exposed for the
    # manifest's ``v_subregion_rate_support`` block.
    assert exp_module._V_SUBREGION_RATE_LABELS == (
        "FWR1",
        "CDR1",
        "FWR2",
        "CDR2",
        "FWR3",
    )
    assert exp_module._V_SUBREGION_RATE_ALIASES == {
        "FWR": ("FWR1", "FWR2", "FWR3"),
        "CDR": ("CDR1", "CDR2"),
    }
    for forbidden in (
        "validate_v_subregion_rates",  # without leading underscore
        "_expand_v_subregion_aliases",
        "_normalize_v_subregion_rates",
    ):
        assert not hasattr(exp_module, forbidden), (
            f"experiment.py now exposes {forbidden!r}; only "
            "``_validate_v_subregion_rates`` is the canonical Slice B "
            "boundary"
        )


# ──────────────────────────────────────────────────────────────────
# 11. Absence — no VSubregionRateWeights struct in Rust
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_rate_weights_struct_in_rust() -> None:
    """Audit §11 Slice B — landed. The Rust crate now carries
    ``VSubregionRateWeights`` (sibling of ``SegmentRateWeights``)
    in ``engine_rs/src/passes/mutate/v_subregion_rates.rs``. The
    five-field shape matches the canonical label order."""
    import subprocess

    result = subprocess.run(
        ["grep", "-rn", "struct VSubregionRateWeights", "engine_rs/src/"],
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
    )
    assert result.stdout.strip(), (
        "Rust no longer declares struct VSubregionRateWeights; "
        "Slice B surface regressed"
    )


def test_pin_present_v_subregion_at_position_helper_in_rust() -> None:
    """Audit §2 — landed. ``v_subregion_at_position`` is the
    sibling of ``segment_at_position`` introduced by Slice B; the
    mutate passes call it once per V-segment site (when the rate
    vector is non-default) to look up the assigned V allele's
    IMGT subregion. Pin at source."""
    import subprocess

    result = subprocess.run(
        ["grep", "-rn", "fn v_subregion_at_position", "engine_rs/src/"],
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
    )
    assert result.stdout.strip(), (
        "Rust no longer declares fn v_subregion_at_position; "
        "Slice B helper regressed"
    )


# ──────────────────────────────────────────────────────────────────
# 12. Absence — no v_subregion_rate_support block in manifest
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_rate_support_in_manifest() -> None:
    """Audit §7 — Slice B landed. The manifest now carries the
    ``v_subregion_rate_support`` block describing the user-facing
    rate kwarg (canonical labels + ``FWR`` / ``CDR`` aliases +
    flat default + the ``in_plan_signature=True`` /
    ``in_content_hash=False`` boundary). The Slice-1
    ``v_subregion_support`` annotation block stays present."""
    m = ga.HUMAN_IGH_OGRDB.cartridge_manifest()
    shm = m["models"]["shm"]
    # Slice 1 annotation block.
    assert "v_subregion_support" in shm
    # Slice B rate-capability block.
    assert "v_subregion_rate_support" in shm, (
        "manifest no longer carries v_subregion_rate_support; "
        "Slice B regressed"
    )
    sup = shm["v_subregion_rate_support"]
    assert sup["available"] is True
    assert sup["labels"] == ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"]
    assert sup["aliases"] == {
        "FWR": ["FWR1", "FWR2", "FWR3"],
        "CDR": ["CDR1", "CDR2"],
    }
    assert sup["default"] == {
        "FWR1": 1.0,
        "CDR1": 1.0,
        "FWR2": 1.0,
        "CDR2": 1.0,
        "FWR3": 1.0,
    }
    assert sup["in_plan_signature"] is True
    assert sup["in_content_hash"] is False


# ──────────────────────────────────────────────────────────────────
# 13. Present — V-subregion mutation counters now partition n_v
# ──────────────────────────────────────────────────────────────────


def test_pin_present_v_subregion_mutation_counters() -> None:
    """Audit §8 — the counters slice landed. ``n_fwr1_mutations``
    / ``n_cdr1_mutations`` / ``n_fwr2_mutations`` /
    ``n_cdr2_mutations`` / ``n_fwr3_mutations`` /
    ``n_v_unannotated_mutations`` now partition ``n_v_mutations``
    on every AIRR record. Two-bucket aggregates remain absent."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(model="s5f", rate=0.03)
        .run_records(n=1, seed=0)
    )
    rec = result.records[0]
    for required in (
        "n_cdr1_mutations",
        "n_cdr2_mutations",
        "n_fwr1_mutations",
        "n_fwr2_mutations",
        "n_fwr3_mutations",
        "n_v_unannotated_mutations",
    ):
        assert required in rec, (
            f"AIRR record missing {required!r}; V-subregion counters "
            "slice regressed"
        )
    for forbidden in ("n_cdr_mutations", "n_fwr_mutations"):
        assert forbidden not in rec, (
            f"two-bucket aggregate {forbidden!r} now exists; counters "
            "audit §1 recommended five-label fields only"
        )


# ──────────────────────────────────────────────────────────────────
# 14. Absence — no compile-time rejection for unsatisfiable rates
# ──────────────────────────────────────────────────────────────────


def test_pin_present_compile_time_rejection_for_unsatisfiable_rates() -> None:
    """Audit §4 — Slice B landed. The DSL now rejects a
    non-default ``v_subregion_rates`` configuration when the
    bound cartridge has zero annotated V alleles. Verifies the
    user-facing error: a cartridge stripped of subregions raises
    ``ValueError`` before the step is appended."""
    import copy

    cfg = copy.deepcopy(ga.HUMAN_IGH_OGRDB)
    for alleles in cfg.v_alleles.values():
        for a in alleles:
            a.subregions = None
            a.gapped_seq = ""
    exp = ga.Experiment.on(cfg).recombine()
    with pytest.raises(ValueError, match="no V-subregion annotations"):
        exp.mutate(
            model="s5f",
            rate=0.03,
            v_subregion_rates={"CDR1": 2.0},
        )
    # Default rates against an unannotated cartridge stay accepted
    # (they're a no-op — the existing fast path runs).
    exp.mutate(model="s5f", rate=0.03)

    # Legacy-name absence: no DSL-side helper that checks
    # subregion availability under a future alternate name.
    from GenAIRR import experiment as exp_module

    for forbidden in (
        "_check_v_subregion_rates_satisfiable",
        "_check_subregion_rate_coverage",
        "_v_subregion_rate_satisfiability",
    ):
        assert not hasattr(exp_module, forbidden), (
            f"experiment.py now exposes {forbidden!r}; Slice B has "
            "landed — flip pin and add a positive test for "
            "unsatisfiable-rate rejection"
        )


# ──────────────────────────────────────────────────────────────────
# 15. Doc anchor — audit doc exists and references contract file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; 14-section structure stays intact."""
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_v_subregion_shm_rate_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted"
    )
    for marker in (
        "## 1. Q1",
        "## 2. Q2",
        "## 3. Q3",
        "## 11. Implementation order",
        "## 14. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
