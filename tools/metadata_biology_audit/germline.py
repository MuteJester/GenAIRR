"""Independent, deterministic germline/recombination metadata audit probes.

These report reproduced discrepancies rather than changing simulator behavior.
Run through the audit entry point, or import ``run_checks``.
"""
from __future__ import annotations

import copy

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR.reference_models import EmpiricalDistributionSpec, ReferenceEmpiricalModels


def _run(exp, seed=0):
    result = exp.run_records(n=1, seed=seed, expose_provenance=True)
    report = result.validate_records(exp.refdata)
    issues = [
        issue
        for failure in report.failures
        for issue in failure["issues"]
    ]
    return result.records[0], issues


def _revision_experiment():
    ref = ge.RefDataConfig.vdj()
    ref.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    ref.add_v_allele("v2*01", "v2", b"GGGTTTAAACCC", anchor=6)
    ref.add_d_allele("d1*01", "d1", b"ACGTTA")
    ref.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return (
        ga.Experiment.on(ref).allow_curatable_refdata()
        .recombine(np1_lengths=[(0, 1)], np2_lengths=[(0, 1)])
        .trim(enabled=False).receptor_revision(prob=1)
    )


def _p_experiment(trim):
    cfg = copy.deepcopy(ga.HUMAN_IGK_OGRDB)
    cfg.reference_models = ReferenceEmpiricalModels(
        p_nucleotide_lengths={
            "V_3": EmpiricalDistributionSpec([(2, 1)]),
            "J_5": EmpiricalDistributionSpec([(1, 1)]),
        }
    )
    return (
        ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)])
        .trim(v_3=[(trim, 1)], j_5=[(trim, 1)])
    )


def run_checks() -> list[dict]:
    checks = []
    rec, issues = _run(_revision_experiment(), seed=0)
    checks.append({
        "id": "revision_stale_v_trim",
        "title": "V replacement exports the original trim and invents deletion columns",
        "expected": {
            "v_call": "v2*01", "v_trim_3": 3,
            "v_germline_end": 9, "v_cigar": "9M",
        },
        "observed": {
            k: rec[k] for k in (
                "v_call", "original_v_call", "v_trim_3",
                "v_germline_end", "v_cigar", "sequence_alignment",
                "germline_alignment",
            )
        },
        "reproduced": rec["v_trim_3"] != 3 and rec["v_cigar"] == "9M3D",
        "validator_issues": issues,
        "source": [
            "engine_rs/src/airr_record/builder.rs:80",
            "engine_rs/src/passes/receptor_revision.rs:642",
            "engine_rs/src/passes/receptor_revision.rs:525",
        ],
    })

    rec, issues = _run(_revision_experiment(), seed=1)
    checks.append({
        "id": "revision_same_allele",
        "title": "Revision reports replacement when it selects the existing V allele",
        "expected": "A replacement event selects a different V allele or reports not applied.",
        "observed": {
            k: rec[k] for k in (
                "original_v_call", "truth_v_call", "v_call",
                "receptor_revision_applied",
            )
        },
        "reproduced": (
            rec["receptor_revision_applied"]
            and rec["original_v_call"] == rec["truth_v_call"]
        ),
        "validator_issues": issues,
        "source": ["engine_rs/src/passes/receptor_revision.rs:427"],
    })

    rec, issues = _run(_p_experiment(trim=0))
    between = rec["sequence"][rec["v_sequence_end"]:rec["j_sequence_start"]]
    checks.append({
        "id": "p_missing_from_np",
        "title": "AIRR combined N/P metadata omits P nucleotides",
        "expected": {"np1": between, "np1_length": len(between)},
        "observed": {
            k: rec[k] for k in (
                "np1", "np1_length", "p_v_3_length", "p_j_5_length",
                "v_sequence_end", "j_sequence_start",
            )
        },
        "reproduced": rec["np1"] != between or rec["np1_length"] != len(between),
        "validator_issues": issues,
        "source": [
            "engine_rs/src/airr_record/builder.rs:191",
            "engine_rs/src/passes/p_addition.rs:291",
        ],
        "references": [
            "https://docs.airr-community.org/en/stable/datarep/rearrangements.html",
        ],
    })

    rec, issues = _run(_p_experiment(trim=1))
    checks.append({
        "id": "p_after_coding_end_deletion",
        "title": "P additions are modeled from exonuclease-trimmed coding ends",
        "classification": "biological model limitation requiring a modeling decision",
        "expected": (
            "Hairpin-derived P bases come from intact original coding ends. "
            "Do not label complements of newly trimmed germline ends as such P bases."
        ),
        "observed": {
            k: rec[k] for k in (
                "v_trim_3", "j_trim_5", "p_v_3_length", "p_j_5_length",
            )
        },
        "reproduced": (
            rec["v_trim_3"] > 0 and rec["p_v_3_length"] > 0
            and rec["j_trim_5"] > 0 and rec["p_j_5_length"] > 0
        ),
        "validator_issues": issues,
        "source": ["engine_rs/src/passes/p_addition.rs:134"],
        "references": ["https://www.imgt.org/textes/FAQ/"],
    })

    ref = ge.RefDataConfig.vj()
    ref.add_v_allele("vLong*01", "vLong", b"CCCCCCCCC", anchor=6)
    ref.add_v_allele("vShort*01", "vShort", b"ACCCCC", anchor=3)
    ref.add_j_allele("j1*01", "j1", b"TTTGGG", anchor=0)
    exp = (
        ga.Experiment.on(ref).allow_curatable_refdata()
        .recombine(np1_lengths=[(0, 1)]).trim(enabled=False)
        .restrict_alleles(v="vLong*01").mutate(model="uniform", rate=1)
    )
    rec, issues = _run(exp, seed=19)
    called_length = {"vLong*01": 9, "vShort*01": 6}[rec["v_call"].split(",")[0]]
    checks.append({
        "id": "shorter_called_allele_hybrid_germline",
        "title": "An alternate shorter allele produces out-of-bounds germline metadata",
        "expected": (
            "Coordinates and germline alignment must refer to the reported allele; "
            "vShort*01 has only 6 bases, so its germline end cannot exceed 6."
        ),
        "observed": {
            **{k: rec[k] for k in (
                "v_call", "recalled_v_call", "truth_v_call", "v_germline_start",
                "v_germline_end", "v_cigar", "germline_alignment",
            )},
            "called_reference_length": called_length,
        },
        "reproduced": (
            rec["v_call"] == "vShort*01" and rec["v_germline_end"] > 6
        ),
        "validator_issues": issues,
        "source": [
            "engine_rs/src/airr_record/walk/mod.rs:246",
            "engine_rs/src/airr_record/segment_projection.rs:212",
        ],
        "scope": (
            "Synthetic high-mutation fixture. Separate built-in human_igh "
            "1000-record sweeps at rates 0.03 and 0.1, seed 42, found no "
            "out-of-bounds coordinates or germline-substring mismatches."
        ),
    })
    return checks
