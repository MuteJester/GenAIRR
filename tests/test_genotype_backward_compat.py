"""The genotype machinery preserves unconstrained simulation output.

``MASTER_DIGEST`` was captured on ``master`` (pre-genotype) for
``Experiment.on(HUMAN_IGH_OGRDB).recombine().run_records(n=100, seed=12345)``
The biological metadata fixes deliberately change calls, alignments, and
translations. The sampling digest below pins the unchanged molecules,
generating alleles, trims, and event counts; derived metadata is checked
by the dedicated biological regression suites.
"""
import hashlib

import GenAIRR as ga
import GenAIRR.data as gdata

MASTER_DIGEST = "3be8e5ea124e1dfff256f93b5ddbb925fdb738d4d6f2eeb5763a81e5b6213460"
SAMPLING_DIGEST = "72332c92093bf1d4dc2f21ae6020ad19458bfdba2e9d16b0cbff021bb898c24f"
SAMPLING_FIELDS = (
    "sequence_id", "sequence", "truth_v_call", "truth_d_call", "truth_j_call",
    "v_trim_5", "v_trim_3", "d_trim_5", "d_trim_3", "j_trim_5", "j_trim_3",
    "n_mutations", "n_indels", "n_pcr_errors", "n_quality_errors",
    "p_v_3_length", "p_d_5_length", "p_d_3_length", "p_j_5_length",
)


def _digest(records):
    h = hashlib.sha256()
    for rec in records:
        h.update(repr(sorted(rec.items())).encode())
    return h.hexdigest()


def test_no_genotype_preserves_sampling_while_correcting_metadata():
    res = ga.Experiment.on(gdata.HUMAN_IGH_OGRDB).recombine().run_records(
        n=100, seed=12345, expose_provenance=True
    )
    sampled_records = []
    for rec in res:
        coding_stops = any(
            rec["sequence"][i:i + 3].upper() in {"TAA", "TAG", "TGA"}
            for i in range(0, len(rec["sequence"]) - 2, 3)
        )
        assert rec["stop_codon"] is coding_stops
        for seg in ("v", "d", "j"):
            assert rec[f"{seg}_call"] == rec[f"truth_{seg}_call"]
        sampled_records.append({k: rec[k] for k in SAMPLING_FIELDS})
    assert _digest(sampled_records) == SAMPLING_DIGEST


def test_no_genotype_output_is_deterministic():
    a = ga.Experiment.on(gdata.HUMAN_IGH_OGRDB).recombine().run_records(n=100, seed=12345)
    b = ga.Experiment.on(gdata.HUMAN_IGH_OGRDB).recombine().run_records(n=100, seed=12345)
    assert _digest(a) == _digest(b)
