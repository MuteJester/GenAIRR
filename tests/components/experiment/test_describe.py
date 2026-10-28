"""describe(): the biology-style narration of a configured experiment.
"""
from __future__ import annotations
from GenAIRR import CompiledExperiment, Experiment
import GenAIRR as ga  # noqa: E402
from tests.support.fixtures.experiment_dsl import _vdj_refdata, _vj_refdata

import pytest

pytestmark = pytest.mark.behavior

def test_describe_empty_experiment_says_no_steps():
    out = Experiment.on(_vj_refdata()).allow_curatable_refdata().describe()
    assert "Experiment on" in out
    assert "no steps" in out.lower()


def test_describe_recombine_only_vj_lists_segments_and_np1():
    out = (
        Experiment.on(_vj_refdata()).allow_curatable_refdata()
        .recombine(np1_lengths=[(0, 1.0), (3, 1.0)])
        .trim(enabled=False)
        .describe()
    )
    # Header line + numbered step
    assert "1. V(D)J recombination" in out
    assert "V/J alleles" in out  # VJ, not VDJ
    assert "NP1" in out
    assert "NP2" not in out  # VJ has no NP2


def test_describe_recombine_vdj_lists_both_np_regions():
    out = (
        Experiment.on(_vdj_refdata()).allow_curatable_refdata()
        .recombine(np1_lengths=[(0, 1.0), (3, 1.0)], np2_lengths=[(2, 1.0)])
        .trim(enabled=False)
        .describe()
    )
    assert "V/D/J alleles" in out
    assert "NP1" in out
    assert "NP2" in out


def test_describe_mutate_s5f_names_kernel():
    out = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=(5, 15))
        .describe()
    )
    assert "S5F" in out
    assert "5–15" in out  # range, en-dash
    # readable kernel label, not the bare s5f_model name
    assert "HH_S5F" in out or "hh_s5f" in out


def test_describe_corrupt_pcr_reads_as_pcr_errors():
    out = (
        Experiment.on(_vj_refdata()).allow_curatable_refdata()
        .recombine(np1_lengths=[(2, 1.0)]).trim(enabled=False)
        .pcr_amplify(count=(0, 3))
        .describe()
    )
    assert "PCR substitution errors" in out
    assert "0–3" in out


def test_describe_clonal_pipeline_shows_visible_fork_divider():
    out = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=50, per_clone=20)
        .mutate(count=8)
        .describe()
    )
    assert "Clonal expansion: 50 lineages × 20 reads/clone" in out
    assert "once per clone" in out
    assert "once per descendant" in out
    # numbering continues across the fork
    assert "1. V(D)J recombination" in out
    assert "2. Somatic hypermutation" in out


def test_describe_includes_metadata_when_stamped():
    out = (
        ga.Experiment.on("human_igh")
        .recombine()
        .with_metadata(donor="P1", tissue="blood")
        .describe()
    )
    assert "donor='P1'" in out
    assert "tissue='blood'" in out


def test_describe_compiled_experiment_renders_productive_constraint_in_biology_terms():
    compiled = (
        ga.Experiment.on("human_igh")
        .recombine()
        .mutate(count=(5, 15))
        .productive_only().compile(allow_curatable_refdata=True)
    )
    out = compiled.describe()
    # Internal contract names collapse to a single biology label.
    assert "productive sequences only" in out
    # And the technical names should NOT leak.
    assert "productive_junction_frame" not in out
    assert "anchor_preserved" not in out


def test_describe_compiled_clonal_experiment_shows_fork_and_numbers_continue():
    compiled = (
        ga.Experiment.on("human_igh")
        .recombine()
        .expand_clones(n_clones=3, per_clone=4)
        .mutate(count=8)
        .pcr_amplify(count=(0, 2))
        .compile(allow_curatable_refdata=True)
    )
    out = compiled.describe()
    assert "Clonal expansion: 3 lineages × 4 reads/clone" in out
    # 1. recombine, 2. mutate, 3. pcr — numbering continues through the fork
    assert "1. V(D)J recombination" in out
    assert "2. Somatic hypermutation" in out
    assert "3. PCR substitution errors" in out
