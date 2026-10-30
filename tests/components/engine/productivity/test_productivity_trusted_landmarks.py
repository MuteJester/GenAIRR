"""Verified coordinates identify noncanonical landmarks without inventing phase."""

import hashlib
import json

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from tests.biology.catalog import verifies

pytestmark = pytest.mark.behavior


def _reference(segment, *, status="verified", phase=0):
    v = b"AAACCCGGG" if segment == "V" else b"AAACCCTGT"
    j = b"GGGAAA" if segment == "J" else b"TTTAAA"
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("V*01", "V", v, anchor=6)
    cfg.add_j_allele("J*01", "J", j, anchor=0)
    if status is not None:
        sequence = v if segment == "V" else j
        cfg.set_allele_evidence(
            segment,
            0,
            json.dumps(
                {
                    "schema_version": 1,
                    "sequence_sha256": "sha256:" + hashlib.sha256(sequence).hexdigest(),
                    "coding_phase": phase,
                    "anchor_status": status,
                    "boundary_5": "unknown",
                    "boundary_3": "unknown",
                    "source": "independently positioned test landmark",
                    "method": "explicit_coordinate_fixture",
                    "notes": [],
                }
            ),
        )
    return cfg


def _experiment(cfg):
    return (
        ga.Experiment.on(cfg)
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(0, 1)])
        .trim(enabled=False)
    )


@verifies("expression.productivity",
          claim="productivity_reasons names which assessment actually failed, "
                "not merely that one did: a rejected anchor residue reports "
                "anchor_residue_rejected and specifically not "
                "reference_anchor_unresolved, so a consumer can tell a "
                "noncanonical landmark from an unresolvable one",
          fields=("productive", "vj_in_frame", "stop_codon", "productivity_reasons"),
          basis="reference")
@pytest.mark.parametrize("segment", ["V", "J"])
def test_verified_noncanonical_landmark_has_defined_junction_and_rejected_productivity(segment):
    cfg = _reference(segment)
    result = _experiment(cfg).run_records(n=1, seed=3)
    rec = result[0]
    assert rec["productive"] is False
    assert rec["vj_in_frame"] is True
    assert rec["stop_codon"] is False
    assert rec["junction_start"] == 6
    assert rec["junction_end"] == 12
    assert rec["junction"] == rec["sequence"][6:12]
    assert "anchor_residue_rejected" in rec["productivity_reasons"]
    assert "reference_anchor_unresolved" not in rec["productivity_reasons"]
    assert result.validate_records(cfg).ok

    rec["productive"] = True
    report = result.validate_records(cfg)
    assert not report.ok
    assert any(
        issue.get("field") == "productive"
        for failure in report.failures
        for issue in failure["issues"]
    )


def test_verified_noncanonical_v_coordinate_does_not_supply_implicit_coding_phase():
    cfg = _reference("V", phase=None)
    result = _experiment(cfg).run_records(n=1, seed=3)
    rec = result[0]
    assert rec["productive"] is False
    assert rec["junction_start"] == 6 and rec["junction_end"] == 12
    assert rec["vj_in_frame"] is None and rec["stop_codon"] is None
    assert "coding_phase_unresolved" in rec["productivity_reasons"]
    assert "anchor_residue_rejected" in rec["productivity_reasons"]
    assert result.validate_records(cfg).ok


def test_verified_noncanonical_j_coordinate_does_not_supply_implicit_j_tail_phase():
    cfg = _reference("J", phase=None)
    result = _experiment(cfg).run_records(n=1, seed=3)
    rec = result[0]
    assert rec["productive"] is False
    assert rec["junction_start"] == 6 and rec["junction_end"] == 12
    # The independent V phase still defines the observed ORF and junction frame.
    assert rec["vj_in_frame"] is True and rec["stop_codon"] is False
    assert "j_tail_unresolved" in rec["productivity_reasons"]
    assert result.validate_records(cfg).ok


@pytest.mark.parametrize("status", [None, "inferred", "unknown", "invalid"])
def test_unverified_noncanonical_reference_does_not_identify_a_landmark(status):
    cfg = _reference("V", status=status, phase=None)
    result = _experiment(cfg).run_records(n=1, seed=3)
    rec = result[0]
    assert rec["productive"] is None
    assert rec["junction_start"] is None and rec["junction_end"] is None
    assert rec["junction"] == ""
    assert rec["vj_in_frame"] is None and rec["stop_codon"] is None
    assert "reference_anchor_unresolved" in rec["productivity_reasons"]
    assert result.validate_records(cfg).ok


@pytest.mark.parametrize("segment", ["V", "J"])
@pytest.mark.parametrize("strict", [False, True])
def test_verified_position_does_not_admit_noncanonical_source_to_productive_sampling(segment, strict):
    cfg = _reference(segment)
    with pytest.raises(ValueError):
        _experiment(cfg).productive_only().run_records(n=1, seed=3, strict=strict)
