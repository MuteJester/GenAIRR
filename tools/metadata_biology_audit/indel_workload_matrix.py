"""Deterministic count/length/probability matrix for the native indel sampler.

Latency includes trace construction and final validation, excludes AIRR
projection. RSS is the process high-water mark (not a per-record allocation).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import statistics
import time

import GenAIRR as ga
from GenAIRR import _engine


def precursor(length):
    reference = ga.RefDataConfig.vdj()
    v = b"ATG" + b"AAA" * ((length - 21) // 3) + b"TGT"
    reference.add_v_allele("v*01", "v", v, anchor=len(v)-3)
    reference.add_d_allele("d*01", "d", b"AAA")
    reference.add_j_allele("j*01", "j", b"TTTAAAAAAAAA", anchor=0)
    compiled = (ga.Experiment.on(reference).recombine(np1_lengths=[(0, 1)], np2_lengths=[(0, 1)])
                .trim(enabled=False).compile())
    outcome = compiled.simulator.run(seed=20261020, strict=True)
    assert outcome.airr_record(compiled.refdata)["productive"] is True
    return compiled.refdata, outcome.final_simulation()


def run(samples):
    profiles = []
    for length in (369, 501, 651):
        reference, state = precursor(length)
        for constrained, probability in ((False, .5), (True, .5), (True, .95)):
            for count in (1, 2, 4, 6):
                plan = _engine.PassPlan()
                plan.push_corrupt_indel([(count, 1.)], insertion_prob=probability)
                simulator = plan.compile(refdata=reference,
                                         respect=_engine.productive() if constrained else None,
                                         strict=True)
                latencies, reports, failures = [], [], []
                cpu_start = time.process_time()
                for seed in range(samples):
                    start = time.perf_counter()
                    try:
                        outcome = simulator.run_from(state, seed=seed, strict=True)
                        reports.extend(outcome.trace().diagnostics())
                    except _engine.StrictSamplingError as error:
                        failures.append({"seed": seed, "args": list(error.args),
                                         "explored_states": getattr(error, "explored_states", None)})
                    latencies.append(time.perf_counter() - start)
                cpu = time.process_time() - cpu_start
                ordered = sorted(latencies)
                profiles.append({
                    "sequence_length": len(state), "reference_hash": reference.content_hash(),
                    "count": count, "insertion_probability": probability, "productive": constrained,
                    "seeds": [0, samples-1], "attempted_records": samples,
                    "completed_records": len(reports), "failures": failures,
                    "throughput_completed_per_second": len(reports) / sum(latencies),
                    "median_ms": statistics.median(latencies)*1000,
                    "p95_ms": ordered[math.ceil(.95*samples)-1]*1000,
                    "maximum_ms": max(latencies)*1000, "process_cpu_seconds": cpu,
                    "process_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    "strategies": {name: sum(d["strategy"] == name for d in reports)
                                   for name in sorted({d["strategy"] for d in reports})},
                    "fallback_completions": sum(d["strategy"] not in
                                                ("independent_tuple_rejection", "unconstrained_sequential")
                                                for d in reports),
                    "maximum_attempts": max((d["attempts"] for d in reports), default=0),
                    "maximum_explored_states": max((d.get("explored_states", 0) for d in reports), default=0),
                })
    return {
        "schema": "indel_workload_matrix_v1", "model": "natural_sequential_v1",
        "extension_sha256": hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
        "samples_per_profile": samples, "profiles": profiles,
        "limitations": ["literal valid references, not an empirical repertoire sample",
                        "single-process latency under ambient machine load",
                        "RSS is cumulative process high-water mark on Linux",
                        "finite seeds do not establish a universal success guarantee"],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=30)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.samples < 1:
        parser.error("--samples must be positive")
    args.output.write_text(json.dumps(run(args.samples), indent=2, allow_nan=False) + "\n")
