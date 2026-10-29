"""Nullable duplex records through the ordinary Result and export boundary."""
import csv
import json

import pytest

from GenAIRR import Experiment, SimulationResult, _engine as ge
from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin
from tests.support.fixtures.duplex_observation import _fixture, _compile

pytestmark = pytest.mark.composition


def _result(*, clipped=False):
    refs, _, seed = _fixture()
    plan = ge.PassPlan()
    if clipped:
        plan.push_corrupt_5prime_loss([(4, 1.0)])
    outcome = _compile(plan, refs).run_from(seed, 0)
    return refs, SimulationResult.from_outcomes([outcome], refs)


@pytest.mark.parametrize("clipped", [False, True])
def test_nullable_result_keeps_full_provenance_and_validates_actual_records(clipped):
    refs, result = _result(clipped=clipped)
    rec = result[0]
    assert rec["sequence_id"] == "seq0"
    assert rec["p_v_3_length"] is None
    assert rec["productive"] is (None if clipped else True)
    assert rec["duplex_repair"]["product"]["sequence"] == "AAATGTCCCTTTAAA"
    assert rec["duplex_repair"]["current_paths"]["v"]["reference_start"] == (4 if clipped else 0)
    assert result.validate_records(refs).ok


@pytest.mark.parametrize("replacement", [False, 0.0, "0"])
def test_nested_provenance_type_tampering_is_detected(replacement):
    refs, result = _result()
    point = result[0]["duplex_repair"]["current_paths"]["v"]["body"][0]
    assert type(point["reference_position"]) is int and point["reference_position"] == 0
    point["reference_position"] = replacement
    report = result.validate_records(refs)
    assert not report.ok
    assert any(i["field"] == "duplex_repair" for i in report.failures[0]["issues"])


def test_nullable_dataframe_and_json_ancestry_cell_preserve_unknown_values(tmp_path):
    _, result = _result(clipped=True)
    pd = pytest.importorskip("pandas")
    frame = result.to_dataframe(airr_strict=True)
    assert str(frame["productive"].dtype) == "boolean"
    assert str(frame["p_v_3_length"].dtype) == "Int64"
    assert str(frame["v_sequence_start"].dtype) == "Int64"
    assert pd.isna(frame.loc[0, "productive"]) and pd.isna(frame.loc[0, "p_v_3_length"])
    assert frame.loc[0, "v_germline_start"] == 5
    assert frame.loc[0, "duplex_repair"]["current_paths"]["v"]["reference_start"] == 4
    target = tmp_path / "duplex.tsv"
    result.to_tsv(target)
    with target.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["productive"] == "" and row["p_v_3_length"] == ""
    assert row["rev_comp"] == "F"
    assert json.loads(row["duplex_repair"]) == result[0]["duplex_repair"]


def test_strict_result_cigar_uses_oriented_inverted_d_offset(tmp_path):
    refs = ge.RefDataConfig.vdj()
    refs.add_v_allele("V*01", "V", b"ACG", anchor=None)
    refs.add_d_allele("D*01", "D", b"AATG")
    refs.add_j_allele("J*01", "J", b"CCC", anchor=None)
    opened = lambda name: DuplexJoin(name).open("left", offset=0).open("right", offset=0)
    repair = CoupledDuplexRepair("ACG", "AATG", "CCC", d_inverted=True).join(
        opened("vd").synapse(right_start=3).ligate()
    ).join(
        opened("dj").resect("left", "top", count=1).resect("left", "bottom", count=1)
        .synapse(right_start=3).ligate()
    )
    seed = repair.seed(refs, v_call="V*01", d_call="D*01", j_call="J*01",
                       reference_ends_are_coding_ends=True)
    outcome = _compile(ge.PassPlan(), refs).run_from(seed, 0)
    result = SimulationResult.from_outcomes([outcome], refs)
    assert result[0]["d_germline_start"] == 1
    assert result[0]["duplex_repair"]["current_paths"]["d"]["oriented_reference_start"] == 0
    target = tmp_path / "inverted.tsv"
    result.to_tsv(target)
    with target.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["d_cigar"] == "3S3M"
    assert row["d_germline_start"] == "2"
    assert result.validate_records(refs).ok


def test_duplex_provenance_cannot_be_overwritten_by_user_metadata():
    with pytest.raises(ValueError, match="reserved"):
        Experiment.on("human_igh").with_metadata(duplex_repair={"forged": True})
