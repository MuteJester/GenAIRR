"""Contract pins for the clonal plan-split architecture.

Companion to
[`docs/clonal_plan_split_design.md`](../docs/clonal_plan_split_design.md).
The audit's predecessor passes (the Bug-A-through-F slices) shipped
the user-visible changes; this contract suite freezes the resulting
architecture as positive scaffold pins plus deferred-architecture
absence pins.

The split:

- ``pin_scaffold_*`` tests freeze the corrected behaviour: the
  ancestor/descendant phase classification, the DSL ordering
  guards (both directions), parent/child trace contents, IR-vs-
  trace projection sourcing, parent/record counts, and validator
  posture.
- ``pin_absence_*`` tests freeze the deferred-architecture
  surfaces (``ClonalFamily`` aggregate, clonal trace-file bundle,
  pre-SHM junction validator, mutation-distance aggregator) so a
  future slice flips them in lockstep.

When a contributor changes how clonal lowering works, the audit
doc + this file are the change-control surface: update both, then
any other layer that needs to follow. Behavioural regression in
the underlying slices surfaces here as well as in the per-slice
spec files.
"""
from __future__ import annotations

import inspect
import re
from pathlib import Path

import pytest

import GenAIRR as ga
from GenAIRR import _compiled, _pipeline_ir
from GenAIRR.experiment import _descendant_phase_step_classifier
from GenAIRR.result import SimulationResult


_REPO_ROOT = Path(__file__).resolve().parent.parent


# ──────────────────────────────────────────────────────────────────
# Shared canonical fixtures
# ──────────────────────────────────────────────────────────────────


def _canonical_clonal(*, n_clones: int = 2, per_clone: int = 3):
    """The user-spec canonical full clonal pipeline. Every
    ancestor-phase + descendant-phase guarded method appears in its
    canonical position relative to ``expand_clones``."""
    return (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .receptor_revision(prob=1.0)
        .expand_clones(n_clones=n_clones, per_clone=per_clone)
        .mutate(count=5)
        .end_loss_5prime(length=10)
        .end_loss_3prime(length=10)
        .random_strand_orientation(prob=0.5)
        .paired_end(r1_length=150, insert_size=300)
    )


def _refdata():
    return ga.Experiment.on("human_igh").compile().refdata


# ──────────────────────────────────────────────────────────────────
# 1. Scaffold — fork-marker structural shape
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_clonal_fork_step_is_partition_marker() -> None:
    """``_ClonalForkStep`` is the structural partition marker the
    audit's plan-split sits on. The doc's §1 describes the
    compile-time split as keyed off this exact type."""
    step = _pipeline_ir._ClonalForkStep(n_clones=3, size=5)
    assert step.n_clones == 3
    assert step.size == 5
    # Frozen dataclass — mutation raises.
    with pytest.raises(Exception):
        step.n_clones = 7  # type: ignore[misc]


def test_pin_scaffold_compile_returns_clonal_when_fork_present() -> None:
    """``Experiment.compile()`` returns ``CompiledClonalExperiment``
    iff a fork is appended. Audit §1 hinges on this branching."""
    plain = ga.Experiment.on("human_igh").recombine().compile()
    clonal = _canonical_clonal(n_clones=1, per_clone=1).compile()
    assert isinstance(plain, _compiled.CompiledExperiment)
    assert not isinstance(plain, _compiled.CompiledClonalExperiment)
    assert isinstance(clonal, _compiled.CompiledClonalExperiment)


# ──────────────────────────────────────────────────────────────────
# 2. Scaffold — ancestor / descendant phase classification
# ──────────────────────────────────────────────────────────────────


_ANCESTOR_PHASE_BUILDERS = {
    "recombine": lambda e: e.recombine(),
    "invert_d": lambda e: e.invert_d(prob=0.5),
    "receptor_revision": lambda e: e.receptor_revision(prob=0.5),
}

_DESCENDANT_PHASE_BUILDERS = {
    "mutate": lambda e: e.mutate(count=5),
    "pcr_amplify": lambda e: e.pcr_amplify(count=3),
    "ambiguous_base_calls": lambda e: e.ambiguous_base_calls(count=4),
    "sequencing_errors": lambda e: e.sequencing_errors(count=3),
    "polymerase_indels": lambda e: e.polymerase_indels(count=2),
    "end_loss_5prime": lambda e: e.end_loss_5prime(length=10),
    "end_loss_3prime": lambda e: e.end_loss_3prime(length=10),
    "random_strand_orientation": lambda e: e.random_strand_orientation(prob=1.0),
    "paired_end": lambda e: e.paired_end(r1_length=150, insert_size=300),
}


@pytest.mark.parametrize("method", sorted(_DESCENDANT_PHASE_BUILDERS))
def test_pin_scaffold_descendant_phase_classifier_recognizes_step(method) -> None:
    """Every guarded descendant-phase DSL method appends a step
    that the classifier recognizes as belonging to that method.
    The classifier is the single source of truth for the unified
    guard."""
    builder = _DESCENDANT_PHASE_BUILDERS[method]
    # Empty Experiment so the appended step is the only thing in
    # ``_steps``; the classifier scans element-by-element.
    base = ga.Experiment.on("human_igh")
    exp = builder(base)
    classifications = [
        _descendant_phase_step_classifier(s) for s in exp._steps
    ]
    assert method in classifications, (
        f"_descendant_phase_step_classifier failed to recognize a "
        f"step from {method!r}; the unified guard is broken."
    )


@pytest.mark.parametrize("method", sorted(_ANCESTOR_PHASE_BUILDERS))
def test_pin_scaffold_ancestor_phase_steps_not_classified_as_descendant(method) -> None:
    """Ancestor-phase steps must NOT be misclassified as
    descendant-phase — that would make the guard fire on
    ``recombine`` / ``invert_d`` / ``receptor_revision`` at
    ``expand_clones`` time, which would break every clonal
    pipeline."""
    builder = _ANCESTOR_PHASE_BUILDERS[method]
    base = ga.Experiment.on("human_igh")
    exp = builder(base)
    for step in exp._steps:
        assert _descendant_phase_step_classifier(step) is None, (
            f"ancestor-phase step from {method!r} misclassified as "
            "descendant-phase."
        )


def test_pin_scaffold_contaminate_unguarded_pin() -> None:
    """``contaminate(prob=...)`` is deliberately NOT in the
    descendant-phase classifier table. Audit §3 documents the
    omission as a deliberate decision; pin it so a follow-up
    that adds it is a tracked change."""
    base = ga.Experiment.on("human_igh").contaminate(prob=0.5)
    for step in base._steps:
        assert _descendant_phase_step_classifier(step) is None, (
            "contaminate now classifies as descendant-phase; audit §3 "
            "decision changed — flip this pin in lockstep."
        )


# ──────────────────────────────────────────────────────────────────
# 3. Scaffold — ordering guards (both directions)
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "method",
    ["invert_d", "receptor_revision"],
)
def test_pin_scaffold_ancestor_phase_post_fork_rejects(method) -> None:
    """The ancestor-phase guards (method-level, opposite direction
    from the unified table) still fire for ``invert_d`` and
    ``receptor_revision`` appended after ``expand_clones``."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=1, per_clone=2)
    )
    with pytest.raises(ValueError, match="before expand_clones"):
        _ANCESTOR_PHASE_BUILDERS[method](exp)


@pytest.mark.parametrize("method", sorted(_DESCENDANT_PHASE_BUILDERS))
def test_pin_scaffold_descendant_phase_pre_fork_rejects(method) -> None:
    """The unified guard at ``expand_clones`` rejects every
    descendant-phase method placed pre-fork."""
    builder = _DESCENDANT_PHASE_BUILDERS[method]
    exp = builder(ga.Experiment.on("human_igh").recombine())
    with pytest.raises(ValueError, match=f"{method} must be called after"):
        exp.expand_clones(n_clones=1, per_clone=2)


@pytest.mark.parametrize("method", sorted(_DESCENDANT_PHASE_BUILDERS))
def test_pin_scaffold_descendant_phase_non_clonal_still_works(method) -> None:
    """Non-clonal pipelines remain unguarded — each descendant-
    phase method can be appended freely when no ``expand_clones``
    is in play. Audit §3 calls this out explicitly."""
    builder = _DESCENDANT_PHASE_BUILDERS[method]
    exp = builder(ga.Experiment.on("human_igh").recombine())
    # Compiles + runs.
    result = exp.run_records(n=2, seed=0)
    assert len(result) == 2


# ──────────────────────────────────────────────────────────────────
# 4. Scaffold — parent trace contents
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_parent_trace_has_recombination_and_ancestor_passes() -> None:
    """For a canonical pipeline with ``invert_d`` and
    ``receptor_revision`` configured, the parent ``Outcome``'s
    ``pass_names()`` lists every ancestor-phase pass. Audit §4
    enumerates the canonical 14-pass set."""
    result = _canonical_clonal(n_clones=1, per_clone=1).run_records(seed=0)
    assert result.parents is not None
    parent = result.parents[0]
    names = parent.pass_names()
    # Recombination machinery.
    for required in (
        "sample_allele.v",
        "sample_allele.d",
        "sample_allele.j",
        "assemble.v",
        "assemble.d",
        "assemble.j",
        "generate_np.np1",
        "generate_np.np2",
    ):
        assert required in names, (
            f"parent missing required ancestor-phase pass {required!r}; "
            f"got {names}"
        )
    # Configured ancestor-phase decisions.
    assert "invert_d" in names
    assert "receptor_revision" in names
    # Non-empty trace + event ledger.
    assert parent.event_count() == len(names)
    assert len(parent.trace()) > 0


# ──────────────────────────────────────────────────────────────────
# 5. Scaffold — descendant trace contents
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_descendant_trace_has_only_descendant_phase_passes() -> None:
    """Descendants' ``pass_names()`` contains exactly the post-fork
    passes from the canonical configuration, in order. None of the
    parent's pass names leak through."""
    result = _canonical_clonal(n_clones=1, per_clone=2).run_records(seed=0)
    assert result.outcomes is not None
    for desc in result.outcomes:
        names = desc.pass_names()
        assert names == [
            "mutate.s5f",
            "corrupt.end_loss.5",
            "corrupt.end_loss.3",
            "corrupt.rev_comp",
            "paired_end",
        ], f"unexpected descendant pass_names: {names}"
        # Specifically: no ancestor-phase passes.
        for forbidden in (
            "sample_allele.v",
            "assemble.v",
            "invert_d",
            "receptor_revision",
        ):
            assert forbidden not in names, (
                f"descendant carries ancestor-phase pass {forbidden!r}; "
                "boundary handoff regressed."
            )


def test_pin_scaffold_descendant_trace_does_not_duplicate_parent() -> None:
    """The descendant's ``trace()`` length is strictly less than a
    non-clonal run of ``recombine + same post-fork plan`` (the
    difference equals the recombine-pass choices the parent ran
    but didn't pass on). Audit §4 calls this out as a load-bearing
    invariant."""
    clonal = _canonical_clonal(n_clones=1, per_clone=1).run_records(seed=0)
    assert clonal.outcomes is not None
    desc_trace_len = len(clonal.outcomes[0].trace())

    # Same post-fork tail as the canonical, layered on a non-clonal
    # recombine pipeline. The ancestor-phase calls go through their
    # method-level guards (they require pre-fork placement and a
    # canonical chain), so we replay only the descendant-phase tail
    # from a fresh recombine.
    non_clonal = (
        ga.Experiment.on("human_igh")
        .recombine()
        .invert_d(prob=1.0)
        .receptor_revision(prob=1.0)
        .mutate(count=5)
        .end_loss_5prime(length=10)
        .end_loss_3prime(length=10)
        .random_strand_orientation(prob=0.5)
        .paired_end(r1_length=150, insert_size=300)
        .compile()
        .run(n=1, seed=0)[0]
    )
    non_clonal_trace_len = len(non_clonal.trace())
    assert desc_trace_len < non_clonal_trace_len, (
        f"descendant trace_len ({desc_trace_len}) >= non-clonal "
        f"({non_clonal_trace_len}); parent trace may have leaked "
        "onto descendants."
    )


# ──────────────────────────────────────────────────────────────────
# 6. Scaffold — parent count + record count
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_parent_count_equals_clone_count_not_record_count() -> None:
    """``len(parents) == n_clones`` and
    ``len(outcomes) == n_clones * per_clone``. Audit §4 § "Parent
    count vs. record count" pin."""
    result = _canonical_clonal(n_clones=3, per_clone=4).run_records(seed=0)
    assert result.parents is not None
    assert len(result.parents) == 3
    assert result.outcomes is not None
    assert len(result.outcomes) == 12


# ──────────────────────────────────────────────────────────────────
# 7. Scaffold — projection inheritance via IR
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_descendants_inherit_ir_sourced_provenance() -> None:
    """Audit §5: ``d_inverted``, ``receptor_revision_applied``,
    ``original_v_call`` survive the parent→descendant boundary via
    the IR (assignments). Verify the descendants' AIRR records
    carry the parent's values in the canonical full clonal
    pipeline."""
    result = _canonical_clonal(n_clones=2, per_clone=3).run_records(seed=0)
    for r in result:
        assert r["d_inverted"] is True, (
            "descendant lost d_inverted; IR-sourced inheritance broken."
        )
        assert r["receptor_revision_applied"] is True, (
            "descendant lost receptor_revision_applied; Bug D fix "
            "regressed."
        )
        assert r["original_v_call"] != "", (
            "descendant lost original_v_call; Bug D fix regressed."
        )
        assert r["original_v_call"].startswith("IGHV")


def test_pin_scaffold_parent_projection_has_no_observation_fields() -> None:
    """Audit §5: projecting the parent ``Outcome`` directly produces
    a record with default values for observation fields
    (``rev_comp=False``, ``end_loss_*_length=0``,
    ``r1_sequence==""``, ``r2_sequence==""``). The parent ran no
    observation-stage passes, so the projection is clean defaults."""
    from GenAIRR._airr_record import outcome_to_airr_record

    result = _canonical_clonal(n_clones=1, per_clone=1).run_records(seed=0)
    refdata = _refdata()
    assert result.parents is not None
    parent_proj = outcome_to_airr_record(
        result.parents[0], refdata, sequence_id="p0"
    )
    # Observation fields ABSENT on parent projection.
    assert parent_proj["rev_comp"] is False
    assert parent_proj["end_loss_5_length"] == 0
    assert parent_proj["end_loss_3_length"] == 0
    assert parent_proj["r1_sequence"] == ""
    assert parent_proj["r2_sequence"] == ""
    # Ancestor IR fields POPULATED on parent projection.
    assert parent_proj["d_inverted"] is True
    assert parent_proj["receptor_revision_applied"] is True
    assert parent_proj["original_v_call"] != ""


# ──────────────────────────────────────────────────────────────────
# 8. Scaffold — descendant-phase divergence is per-descendant
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_random_strand_diverges_across_siblings() -> None:
    """Audit §5: ``rev_comp`` is descendant-trace-sourced and
    sampled independently per descendant. At ``prob=0.5`` with
    multiple siblings, both values should appear in the batch —
    proves per-descendant divergence."""
    exp = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=2, per_clone=4)
        .random_strand_orientation(prob=0.5)
    )
    result = exp.run_records(seed=0)
    rev_comp_values = {r["rev_comp"] for r in result}
    assert rev_comp_values == {True, False}, (
        f"rev_comp not diverging across siblings: only got "
        f"{rev_comp_values}. Per-descendant independence broken."
    )


# ──────────────────────────────────────────────────────────────────
# 9. Scaffold — validator posture on the canonical stack
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_canonical_clonal_passes_all_three_validators() -> None:
    """Audit §6: the canonical full clonal pipeline produces records
    that pass the per-record postcondition validator, the field-
    only family validator, AND the parent-aware family validator
    with refdata. All three reports are ok."""
    result = _canonical_clonal(n_clones=2, per_clone=3).run_records(
        seed=0, expose_provenance=True
    )
    refdata = _refdata()
    record_report = result.validate_records(refdata)
    assert record_report.ok, f"per-record failures: {record_report.failures[:3]}"

    family_report = result.validate_families()
    assert family_report.ok, f"family failures: {family_report.failures[:3]}"
    assert family_report.family_count == 2

    parent_report = result.validate_families_with_parents(refdata)
    assert parent_report.ok, (
        f"parent-aware failures: {parent_report.failures[:3]}"
    )
    assert parent_report.family_count == 2


def test_pin_scaffold_validate_records_true_does_not_invoke_parent_aware() -> None:
    """Audit §6: ``validate_records=True`` on a clonal
    ``run_records`` runs per-record + field-only family, NOT
    parent-aware. Pinned via source inspection of the clonal
    ``run_records`` body so a refactor that auto-wires the
    parent-aware validator surfaces here."""
    src = inspect.getsource(_compiled.CompiledClonalExperiment.run_records)
    # The field-only family raiser must still be called.
    assert "_raise_on_family_validation_failure" in src
    # The parent-aware validator must NOT be invoked automatically.
    assert "validate_families_with_parents" not in src, (
        "validate_records=True now auto-invokes the parent-aware "
        "validator; audit §6 architectural decision changed — flip "
        "this pin in lockstep with the slice that made the change."
    )


# ──────────────────────────────────────────────────────────────────
# 10. Scaffold — `_descendant_phase_step_classifier` is the SoT
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_classifier_is_module_level_helper() -> None:
    """The classifier lives at module scope so the unified guard at
    ``expand_clones`` and any future contributor can reuse it.
    Pinned so a refactor that buries it inside ``expand_clones``
    surfaces here."""
    from GenAIRR import experiment as exp_module

    assert hasattr(exp_module, "_descendant_phase_step_classifier")
    assert callable(exp_module._descendant_phase_step_classifier)


# ──────────────────────────────────────────────────────────────────
# 11. Absence — deferred architecture
# ──────────────────────────────────────────────────────────────────


def test_pin_absence_no_clonalfamily_aggregate() -> None:
    """No ``ClonalFamily`` / ``FamilyOutcome`` / ``FamilyRecord``
    aggregate type yet. The audit's §13 out-of-scope list keeps
    this deferred until lineage biology slices arrive."""
    public_names = set(getattr(ga, "__all__", []))
    forbidden = {"ClonalFamily", "FamilyOutcome", "FamilyRecord", "ClonalFamilyRecord"}
    assert not (public_names & forbidden), (
        f"public package now exports family aggregate type(s) "
        f"{public_names & forbidden}; flip pin in lockstep."
    )


def test_pin_absence_no_clonal_trace_file_bundle() -> None:
    """No clonal trace-file bundling helper on
    ``CompiledClonalExperiment``. Per-outcome ``trace_file_from``
    still lives on the underlying simulators; a clonal-aware
    bundler is Slice 5+ scope."""
    clonal = _canonical_clonal(n_clones=1, per_clone=1).compile()
    for forbidden in (
        "trace_file_from_family",
        "clonal_trace_file_from",
        "family_trace_file",
        "trace_file_from_clone",
    ):
        assert not hasattr(clonal, forbidden), (
            f"CompiledClonalExperiment.{forbidden} now exists; flip "
            "pin."
        )


def test_pin_absence_no_pre_shm_junction_validator() -> None:
    """Pre-SHM junction invariance validator still deferred (needs
    either a parent-derived projected field or a future
    ``FamilyRecord``). Audit §13."""
    for forbidden in (
        "validate_pre_shm_junction_invariance",
        "validate_families_pre_shm_junction",
        "validate_pre_shm_junction",
    ):
        assert not hasattr(SimulationResult, forbidden), (
            f"SimulationResult.{forbidden} now exists; flip pin."
        )


def test_pin_absence_no_mutation_distance_aggregator() -> None:
    """Mutation-distance distribution aggregator still deferred.
    Audit §13."""
    result = _canonical_clonal(n_clones=1, per_clone=1).run_records(seed=0)
    for forbidden in (
        "mutation_distances_per_clone",
        "family_distance_distribution",
        "clone_distance_matrix",
        "per_family_shm_distance",
        "mutations_per_clone",
    ):
        assert not hasattr(result, forbidden), (
            f"SimulationResult.{forbidden} now exists; flip pin."
        )


# ──────────────────────────────────────────────────────────────────
# 12. Doc anchor — audit doc exists + references this file
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc must continue to exist and reference this
    contract file; the canonical 14-section structure is intact.
    A regression here means the change-control surface drifted."""
    docs_dir = Path(__file__).resolve().parent.parent / "docs"
    if not docs_dir.is_dir():
        import pytest
        pytest.skip("docs/ is contributor-only; not present in this checkout")
    doc_path = _REPO_ROOT / "docs" / "clonal_plan_split_design.md"
    assert doc_path.exists(), "clonal_plan_split_design.md missing"
    doc = doc_path.read_text(encoding="utf-8")
    assert "test_clonal_plan_split_contract.py" in doc, (
        "audit doc no longer references the contract file; lockstep "
        "convention drifted."
    )
    for marker in (
        "## 1. Current state",
        "## 5. Projection inheritance",
        "## 12. Test surface",
        "## 14. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
