"""The genotype machinery preserves unconstrained simulation output.

``MASTER_DIGEST`` was captured on ``master`` (pre-genotype) for
``Experiment.on(HUMAN_IGH_OGRDB).recombine().run_records(n=100, seed=12345)``
The biological metadata fixes deliberately change calls, alignments, and
translations. The sampling digest below pins the unchanged molecules,
generating alleles, trims, and event counts; derived metadata is checked
by the dedicated biological regression suites.

**What this digest can and cannot be asked.** It answers "did a change to the
genotype machinery move the no-genotype sampling path", which is the question
it was written for, and it answers it well: anything that perturbs an
unconstrained run shows up here immediately. It is not a cross-version
promise, because the engine's sampling model is versioned separately and a
bump moves every seeded draw by design. ``SAMPLING_DIGEST`` was recaptured
once, at ``DETERMINISM_VERSION`` 2, when per-pass RNG substreams landed; a
later bump will move it again, and the move is correct exactly when
``DETERMINISM_VERSION`` moved with it. ``MASTER_DIGEST`` is older still and
describes the pre-genotype model under version 1.
"""
import hashlib

import GenAIRR as ga
import GenAIRR.data as gdata

import pytest

pytestmark = pytest.mark.behavior

MASTER_DIGEST = "3be8e5ea124e1dfff256f93b5ddbb925fdb738d4d6f2eeb5763a81e5b6213460"
#: Recaptured at DETERMINISM_VERSION 2 (per-pass RNG substreams). See the
#: module docstring for when moving this is correct and when it is a bug.
SAMPLING_DIGEST = "0c29c80527193f4f12d7ed52d378da7f29ed9f753e7726051b3d2cc1582b1905"
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
    assert _digest(sampled_records) == SAMPLING_DIGEST, (
        f"the unconstrained sampling path moved. If that was deliberate, "
        f"DETERMINISM_VERSION (now {ga.DETERMINISM_VERSION}) has to move with "
        "it and this digest is recaptured; if it was not, a change to the "
        "genotype machinery has reached a run that uses no genotype.")


def test_no_genotype_output_is_deterministic():
    a = ga.Experiment.on(gdata.HUMAN_IGH_OGRDB).recombine().run_records(n=100, seed=12345)
    b = ga.Experiment.on(gdata.HUMAN_IGH_OGRDB).recombine().run_records(n=100, seed=12345)
    assert _digest(a) == _digest(b)
