"""Reproduce the two formerly retry-exhausted human IGH indel precursors.

Run from the repository with an up-to-date native extension:
    .venv/bin/python tools/metadata_biology_audit/indel_sampler_benchmark.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import statistics
import time

import GenAIRR as ga
from GenAIRR import _engine


def precursor_experiment(original_end: bool):
    experiment = ga.Experiment.on("human_igh").recombine(
        np1_lengths=[(n, 1) for n in range(7)],
        np2_lengths=[(n, 1) for n in range(7)],
    )
    if original_end:
        experiment.process_coding_ends(
            reference_ends_are_coding_ends=True,
            v_3=[(-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1)],
            j_5=[(-1, 1), (0, 1), (1, 1)],
            d_5=[(-1, 1), (0, 1), (1, 1)],
            d_3=[(-2, 1), (0, 1), (2, 1)],
        )
    return experiment.invert_d(prob=0.5).productive_only().mutate(model="uniform", rate=0.12)


def run(samples: int = 40):
    reports = []
    for original, seed in ((True, 20261028), (False, 20261027)):
        compiled = precursor_experiment(original).compile()
        precursor = compiled.simulator.run(seed=seed, strict=True).final_simulation()
        plan = _engine.PassPlan()
        plan.push_corrupt_indel([(2, 1.0)])
        sampler = plan.compile(refdata=compiled.refdata, respect=_engine.productive(), strict=True)
        durations, diagnostics = [], []
        for continuation_seed in range(samples):
            start = time.perf_counter()
            outcome = sampler.run_from(precursor, seed=continuation_seed, strict=True)
            durations.append(time.perf_counter() - start)
            assert outcome.airr_record(compiled.refdata)["productive"] is True
            report = outcome.trace().diagnostics()[0]
            assert report["requested_count"] == report["applied_count"] == 2
            diagnostics.append(report)
        reports.append({
            "sampler_model": "natural_sequential_v1",
            "reference_hash": compiled.refdata.content_hash(),
            "extension_sha256": hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
            "original_end": original,
            "precursor_seed": seed,
            "sequence_length": len(precursor),
            "continuation_seeds": [0, samples - 1],
            "successes": samples,
            "median_ms": statistics.median(durations) * 1000,
            "maximum_ms": max(durations) * 1000,
            "strategies": sorted({d["strategy"] for d in diagnostics}),
            "attempts": [d["attempts"] for d in diagnostics],
        })
    return reports


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=40)
    args = parser.parse_args()
    if args.samples < 1:
        parser.error("--samples must be positive")
    print(json.dumps(run(args.samples), indent=2))
