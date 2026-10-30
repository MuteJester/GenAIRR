"""The SimulationResult wrapper and its AIRR export.
"""
from __future__ import annotations
import pytest
import GenAIRR as ga  # noqa: E402
from tests.support.fixtures.experiment_dsl import _human_igh_exp, _human_igh_records

pytestmark = pytest.mark.behavior

def test_run_records_returns_simulation_result():
    from GenAIRR.result import SimulationResult

    result = _human_igh_records(n=3)
    assert isinstance(result, SimulationResult)
    assert len(result) == 3


def test_simulation_result_is_listlike():
    result = _human_igh_records(n=4)

    # __len__
    assert len(result) == 4

    # __getitem__ by int
    rec0 = result[0]
    assert isinstance(rec0, dict)

    # __getitem__ by slice
    pair = result[1:3]
    assert isinstance(pair, list)
    assert len(pair) == 2

    # __iter__
    seen = list(result)
    assert len(seen) == 4
    assert all(isinstance(r, dict) for r in seen)


def test_simulation_result_repr():
    result = _human_igh_records(n=2)
    assert repr(result) == "<SimulationResult n=2>"


def test_simulation_result_exposes_outcomes_and_records():
    result = _human_igh_records(n=2)
    assert result.outcomes is not None
    assert len(result.outcomes) == 2
    # ``records`` is the underlying list view.
    assert result.records is not None
    assert len(result.records) == 2


def test_simulation_result_from_records_only_has_no_outcomes():
    from GenAIRR.result import SimulationResult

    rec = {"sequence": "ACGT", "v_call": "x", "junction": ""}
    result = SimulationResult([rec])
    assert len(result) == 1
    assert result.outcomes is None
    assert result[0]["sequence"] == "ACGT"


def test_record_core_fields_populated():
    result = _human_igh_records(n=1)
    rec = result[0]
    # Core sequence fields
    assert isinstance(rec["sequence"], str)
    assert rec["sequence_length"] == len(rec["sequence"])
    assert rec["sequence_length"] > 0


def test_record_v_call_resolves_to_imgt_name():
    result = _human_igh_records(n=3)
    for rec in result:
        v_call = rec["v_call"]
        assert isinstance(v_call, str)
        # Human IGH V allele names follow the IGHV...*NN pattern.
        assert v_call.startswith("IGH"), v_call
        assert "*" in v_call, v_call


def test_record_d_and_j_calls_resolve_for_vdj():
    result = _human_igh_records(n=3)
    for rec in result:
        assert rec["d_call"].startswith("IGH"), rec["d_call"]
        assert rec["j_call"].startswith("IGH"), rec["j_call"]
        assert "*" in rec["j_call"]


def test_record_d_call_empty_for_vj_chain():
    result = (
        ga.Experiment.on("human_igk").recombine().run_records(n=2, seed=7)
    )
    for rec in result:
        assert rec["d_call"] == ""
        assert rec["v_call"].startswith("IGK")
        assert rec["j_call"].startswith("IGK")


def test_record_v_d_j_coordinates_are_within_sequence():
    result = _human_igh_records(n=3)
    for rec in result:
        slen = rec["sequence_length"]
        for prefix in ("v", "d", "j"):
            start = rec[f"{prefix}_sequence_start"]
            end = rec[f"{prefix}_sequence_end"]
            assert start is not None and end is not None
            assert 0 <= start <= end <= slen


def test_record_junction_window_is_consistent():
    # Junction is derived from V-anchor → J-anchor + 3. Verify the
    # slice matches the recorded start/end and the length is correct.
    result = _human_igh_records(n=5)
    for rec in result:
        if rec["junction_start"] is None:
            continue
        seq = rec["sequence"]
        start = rec["junction_start"]
        end = rec["junction_end"]
        assert seq[start:end] == rec["junction"]
        assert rec["junction_length"] == end - start


def test_record_productive_flag_is_consistent_with_frame_and_stop():
    result = _human_igh_exp().productive_only().run_records(n=10, seed=42)
    for rec in result:
        assert rec["productive"] is True
        assert rec["vj_in_frame"] is True
        assert rec["stop_codon"] is False
        # Junction must be in-frame and translate cleanly.
        assert rec["junction_length"] is not None
        assert rec["junction_length"] % 3 == 0
        assert "*" not in rec["junction_aa"]


def test_record_n_mutations_matches_trace():
    # When .mutate(count=N) is fixed, the record's ``n_mutations``
    # should equal the trace count.
    result = _human_igh_exp().mutate(count=12).run_records(n=2, seed=0)
    for rec in result:
        assert rec["n_mutations"] == 12
        assert rec["mutation_rate"] == pytest.approx(
            12 / rec["sequence_length"]
        )


def test_record_corruption_counters_default_to_zero():
    result = _human_igh_records(n=2)
    for rec in result:
        assert rec["n_pcr_errors"] == 0
        assert rec["n_quality_errors"] == 0
        assert rec["n_indels"] == 0
        assert rec["is_contaminant"] is False


def test_record_corruption_counters_track_pcr_errors():
    # `count=4` is four attempts. A replacement is drawn from A/C/G/T
    # including the base already present, so the counter reports the subset
    # that changed a base. Bounded above by the attempts and asserted
    # non-zero in aggregate, which a counter reading the drawn count would
    # fail by hitting the attempt total exactly.
    result = (
        _human_igh_exp().pcr_amplify(count=4).run_records(n=25, seed=0)
    )
    total = 0
    for rec in result:
        assert rec["n_pcr_errors"] <= 4
        total += rec["n_pcr_errors"]
    assert 0 < total < 100, (
        f"{total} realized errors across 25 records of 4 attempts is either "
        f"nothing or every attempt")


def test_compiled_run_records_matches_experiment_run_records():
    # The two paths (Experiment.run_records vs CompiledExperiment.run_records)
    # should produce identical record sequences for the same seed.
    exp = _human_igh_exp().mutate(count=5)
    a = exp.run_records(n=3, seed=99)
    b = exp.compile(allow_curatable_refdata=True).run_records(n=3, seed=99)
    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert x["sequence"] == y["sequence"]
        assert x["v_call"] == y["v_call"]
        assert x["junction"] == y["junction"]


def test_to_tsv_roundtrip(tmp_path):
    import csv as _csv

    result = _human_igh_records(n=3)
    path = tmp_path / "batch.tsv"
    result.to_tsv(str(path))

    with open(path, encoding="utf-8") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))

    assert len(rows) == 3
    # Spot-check a few canonical columns.
    for row, rec in zip(rows, result):
        assert row["sequence"] == rec["sequence"]
        assert row["v_call"] == rec["v_call"]
        assert int(row["sequence_length"]) == rec["sequence_length"]


def test_to_csv_uses_comma_delimiter(tmp_path):
    result = _human_igh_records(n=2)
    path = tmp_path / "batch.csv"
    result.to_csv(str(path))

    text = path.read_text(encoding="utf-8")
    header = text.splitlines()[0]
    assert "," in header
    assert "\t" not in header
    assert "sequence" in header
    assert "v_call" in header


def test_to_tsv_writes_canonical_column_order(tmp_path):
    from GenAIRR.result import _DEFAULT_COLUMN_ORDER

    result = _human_igh_records(n=1)
    path = tmp_path / "order.tsv"
    result.to_tsv(str(path))

    header = path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    # Header must start with the canonical column order.
    assert header[: len(_DEFAULT_COLUMN_ORDER)] == _DEFAULT_COLUMN_ORDER


def test_to_tsv_serializes_none_as_empty(tmp_path):
    from GenAIRR.result import SimulationResult

    rec = {
        "sequence": "ACGT",
        "sequence_length": 4,
        "v_call": "v1",
        "d_call": "",
        "junction": "",
        "junction_start": None,
        "junction_end": None,
    }
    path = tmp_path / "none.tsv"
    SimulationResult([rec]).to_tsv(str(path))
    text = path.read_text(encoding="utf-8")
    # ``None`` columns should render as empty fields, not literal "None".
    assert "None" not in text


def test_to_fasta_emits_header_and_sequence(tmp_path):
    result = _human_igh_records(n=2)
    path = tmp_path / "batch.fa"
    result.to_fasta(str(path))

    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 4  # 2 records × (header + seq)
    assert lines[0].startswith(">seq0|v_call=")
    assert "|j_call=" in lines[0]
    assert lines[1] == result[0]["sequence"]
    assert lines[2].startswith(">seq1|v_call=")
    assert lines[3] == result[1]["sequence"]


def test_to_fasta_respects_custom_prefix(tmp_path):
    result = _human_igh_records(n=1)
    path = tmp_path / "custom.fa"
    result.to_fasta(str(path), prefix="read_")
    first = path.read_text(encoding="utf-8").splitlines()[0]
    assert first.startswith(">read_0|")
