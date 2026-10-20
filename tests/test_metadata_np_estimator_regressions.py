"""Keep AIRR combined N/P observations out of N-only generation models."""
import copy
import io

import pytest

import GenAIRR as ga
from GenAIRR.reference_models import EmpiricalDistributionSpec, ReferenceEmpiricalModels


_P_COLUMNS = ("p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length")
_ESTIMATORS = ("estimate_np_length_distributions", "estimate_np_base_model")


def _builder():
    return ga.ReferenceCartridgeBuilder.from_fasta(
        v_fasta=">IGKV1-MOCK*01\nAAACCCTGT\n",
        j_fasta=">IGKJ1-MOCK*01\nTGGAAA\n",
        chain_type="BCR_LIGHT_KAPPA",
    )


def _estimate(builder, estimator, records):
    options = {"kind": "empirical_first_base"} if estimator.endswith("base_model") else {}
    return getattr(builder, estimator)(records, **options)


@pytest.mark.parametrize("estimator", _ESTIMATORS)
@pytest.mark.parametrize("column", _P_COLUMNS)
@pytest.mark.parametrize("value", [2, "2"])
def test_nonzero_p_provenance_rejects_before_changing_model_or_report(
    estimator, column, value,
):
    builder = _builder()
    _estimate(builder, estimator, [{"np1": "CAT", "np1_length": 3}])
    original_models = builder._reference_models
    original_stages = copy.deepcopy(builder.report().stages)
    records = [
        {"np1": "AAA", "np1_length": 3},
        {"np1": "TGCATT", "np1_length": 6, column: value},
    ]
    with pytest.raises(ValueError, match=rf"row 1 has nonzero {column}.*combine N and P"):
        _estimate(builder, estimator, records)
    assert builder._reference_models is original_models
    assert builder.report().stages == original_stages


@pytest.mark.parametrize("estimator", _ESTIMATORS)
@pytest.mark.parametrize("provenance", ["missing", "zero", "tsv_zero"])
def test_missing_and_zero_p_provenance_remain_accepted(estimator, provenance):
    builder = _builder()
    row = {"np1": "CAT", "np1_length": 3}
    if provenance == "zero":
        row.update({column: 0 for column in _P_COLUMNS})
    if provenance == "tsv_zero":
        records = io.StringIO(
            "np1\tnp1_length\t" + "\t".join(_P_COLUMNS) + "\n"
            "CAT\t3\t0\t0\t0\t0\n"
        )
    else:
        records = [row]
    _estimate(builder, estimator, records)
    models = builder.build().reference_models
    if estimator.endswith("base_model"):
        assert dict(models.np_bases["NP1"].first_base) == {
            "A": 1 / 3, "C": 1 / 3, "T": 1 / 3,
        }
    else:
        assert dict(models.np_lengths["NP1"].values) == {3: 1.0}


@pytest.mark.parametrize("estimator", _ESTIMATORS)
def test_corrected_generated_np_record_cannot_double_count_p_during_estimation(estimator):
    cfg = copy.deepcopy(ga.HUMAN_IGK_OGRDB)
    cfg.reference_models = ReferenceEmpiricalModels(
        p_nucleotide_lengths={
            "V_3": EmpiricalDistributionSpec([(2, 1)]),
            "J_5": EmpiricalDistributionSpec([(1, 1)]),
        },
    )
    record = (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(3, 1)])
        .trim(enabled=False).run_records(n=1, seed=0)[0]
    )
    assert record["np1_length"] == len(record["np1"]) == 6
    assert record["np1"][2:-1] == "CAT"
    assert (record["p_v_3_length"], record["p_j_5_length"]) == (2, 1)

    builder = _builder()
    builder.estimate_p_nucleotide_lengths([record])
    with pytest.raises(ValueError, match="AIRR NP fields combine N and P"):
        _estimate(builder, estimator, [record])
    models = builder.build().reference_models
    assert dict(models.p_nucleotide_lengths["V_3"].values) == {2: 1.0}
    assert not models.np_lengths
    assert not models.np_bases
