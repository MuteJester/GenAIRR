"""Regressions from the independent biological metadata audit."""

import csv

import pytest

import GenAIRR as ga

pytestmark = pytest.mark.behavior


@pytest.fixture
def clean_batch():
    compiled = ga.Experiment.on("human_igh").recombine().compile()
    return compiled.run_records(n=2, seed=42), compiled.refdata


@pytest.mark.parametrize("field,value", [
    ("sequence", "AAA"), ("productive", "wrong"),
    ("junction_start", 999999), ("rev_comp", 0),
])
def test_validator_checks_returned_field_values_and_types(clean_batch, field, value):
    result, refdata = clean_batch
    result[0][field] = value
    report = result.validate_records(refdata)
    assert not report.ok
    assert any(
        issue["kind"] == "RecordFieldMismatch" and issue["field"] == field
        for failure in report.failures for issue in failure["issues"]
    )


def test_validator_detects_missing_fields(clean_batch):
    result, refdata = clean_batch
    del result[0]["junction"]
    report = result.validate_records(refdata)
    assert report.summary()["RecordFieldMissing"] == 1


def test_validator_detects_missing_outcomes(clean_batch):
    result, refdata = clean_batch
    result.outcomes.pop()
    report = result.validate_records(refdata)
    assert report.summary()["RecordOutcomeCountMismatch"] == 1
    assert report.count == 2


@pytest.mark.parametrize("field", [
    "sequence", "productive", "v_call", "p_v_3_length", "truth_v_call", "clone_id",
    "recalled_v_call", "recalled_d_call", "recalled_j_call",
])
def test_metadata_cannot_overwrite_generated_fields(field):
    with pytest.raises(ValueError, match="generated fields"):
        ga.Experiment.on("human_igh").with_metadata(**{field: "wrong"})


def test_custom_metadata_remains_allowed_and_validates(clean_batch):
    result, refdata = clean_batch
    result[0]["sample_id"] = "sample-1"
    assert result.validate_records(refdata).ok


def test_metadata_cannot_replace_genotype_provenance():
    genotype = (
        ga.Genotype.from_dataconfig(ga.HUMAN_IGH_OGRDB)
        .complete_from_reference().with_subject("subject-1")
    )
    experiment = (
        ga.Experiment.on("human_igh").with_genotype(genotype).recombine()
        .with_metadata(subject_id="wrong-subject")
    )
    with pytest.raises(ValueError, match="reserved generated field"):
        experiment.run_records(n=1, seed=42)


def test_default_tsv_is_airr_and_keeps_python_records_unchanged(clean_batch, tmp_path):
    result, _ = clean_batch
    original = dict(result[0])
    path = tmp_path / "rearrangements.tsv"
    result.to_tsv(str(path))
    with path.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["v_sequence_start"] == "1"
    assert row["rev_comp"] == "F"
    assert row["productive"] in ("T", "F", "")
    for segment in ("d", "j"):
        start = original[f"{segment}_sequence_start"]
        ref_start = original[f"{segment}_germline_start"]
        prefix = (f"{start}S" if start else "") + (f"{ref_start}N" if ref_start else "")
        assert row[f"{segment}_cigar"] == prefix + original[f"{segment}_cigar"]
    assert result[0] == original


def test_legacy_tsv_coordinate_mode_is_explicit(clean_batch, tmp_path):
    result, _ = clean_batch
    path = tmp_path / "python.tsv"
    result.to_tsv(str(path), airr_strict=False)
    with path.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["v_sequence_start"] == "0"
    assert row["rev_comp"] == "False"


def test_reverse_observation_preserves_coding_annotations():
    forward = ga.Experiment.on("human_igh").recombine().run_records(n=4, seed=42)
    experiment = (
        ga.Experiment.on("human_igh").recombine().random_strand_orientation(prob=1)
    )
    reverse = experiment.run_records(n=4, seed=42)
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    for fwd, rev in zip(forward, reverse):
        assert rev["sequence"] == fwd["sequence"].translate(complement)[::-1]
        assert rev["rev_comp"] is True
        for field in (
            "junction", "junction_aa", "junction_start", "junction_end", "sequence_aa",
            "np1", "np2", "np1_aa", "np2_aa", "productive", "stop_codon",
            "v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end",
            "sequence_alignment", "germline_alignment",
        ):
            assert rev[field] == fwd[field], field
    assert reverse.validate_records(experiment.refdata).ok


def test_complete_read_loss_keeps_realized_provenance():
    reference = ga.Experiment.on("human_igh").recombine().run_records(n=1, seed=42)[0]
    experiment = ga.Experiment.on("human_igh").recombine().end_loss_5prime(length=9999)
    result = experiment.run_records(n=1, seed=42)
    assert result[0]["sequence"] == ""
    assert result[0]["end_loss_5_length"] == len(reference["sequence"])
    assert result.validate_records(experiment.refdata).ok
