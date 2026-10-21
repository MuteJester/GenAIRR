"""Isolated full-length synthetic durable-family timing and process RSS.

Uses only GenAIRR and the Python standard library. Each mode/case/repeat is
run in a fresh process. RSS is whole-process resident memory/high water, not
the incremental allocation of a family or proof of native Arc sharing.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
import time

from tools.metadata_biology_audit.duplex_family_benchmark import _identity, _rss


ROOT = Path(__file__).resolve().parents[2]
CASES = {
    "shallow": {"max_generations": 2, "n_max": 8, "n_sample": 2, "lambda_base": 2.4},
    "deeper": {"max_generations": 6, "n_max": 8, "n_sample": 2, "lambda_base": 2.4},
    "survival_retries": {"max_generations": 3, "n_max": 8, "n_sample": 2, "lambda_base": 0.7},
}
COMMON = {"n_clones": 1, "rate": 0.01, "s5f_model": "hh_s5f", "selection_strength": 0.0,
          "strict": True, "validate_records": True, "allow_extinction": False}


def _digest(value):
    from GenAIRR import _engine
    return hashlib.sha256(_engine._family_replay_pack(value).encode()).hexdigest()


def _truth(result):
    return {
        "records_sha256": _digest(result.records),
        "trees_sha256": _digest([tree.to_node_table_tsv() for tree in result.lineage_trees]),
        "observed_records": len(result),
        "accepted_nodes": sum(len(tree.nodes()) for tree in result.lineage_trees),
        "accepted_max_generation": max(node.generation for tree in result.lineage_trees for node in tree.nodes()),
        "sample_draws": sum(record["duplicate_count"] for record in result),
        "observed_sequence_lengths": sorted({len(record["sequence"]) for record in result}),
    }


def _experiment(case):
    import GenAIRR as ga
    from GenAIRR.experimental import CoupledDuplexRepair, DuplexHistory, DuplexJoin, DuplexRepairModel

    # Coding length and landmarks are synthetic controls, not fitted germlines.
    codons = ("GCT", "GAC", "CAG", "GGT", "AAG", "TCT", "GAA", "ACT", "GTG", "TAC")
    v = "".join(codons[i % len(codons)] for i in range(97)) + "TGTGCTGGT"
    d = "GGTGACTACGGTCAGGCT"
    j = "TGG" + "".join(codons[i % len(codons)] for i in range(16))
    assert (len(v), len(d), len(j)) == (300, 18, 51)
    ref = ga.RefDataConfig.vdj()
    ref.add_v_allele("IGHV-SYN*01", "IGHV-SYN", v.encode(), anchor=291)
    ref.add_d_allele("IGHD-SYN*01", "IGHD-SYN", d.encode())
    ref.add_j_allele("IGHJ-SYN*01", "IGHJ-SYN", j.encode(), anchor=0)
    repair = CoupledDuplexRepair(v, d, j)
    repair = repair.join(
        DuplexJoin("vd").open("left", offset=0).open("right", offset=0)
        .resect("left", "top", count=1).synapse(right_start=len(v)).fill("top").ligate()
    ).join(
        DuplexJoin("dj").open("left", offset=0).open("right", offset=0)
        .synapse(right_start=len(d)).ligate()
    )
    model = DuplexRepairModel.bind(
        ref, [DuplexHistory(name, repair, "IGHV-SYN*01", "IGHD-SYN*01", "IGHJ-SYN*01")
              for name in ("restored-a", "restored-b")],
        reference_ends_are_coding_ends=True,
    )
    exp = ga.Experiment.on(ref).recombine(duplex_repair=model)
    exp.clonal_lineage(**CASES[case], **{k: COMMON[k] for k in
                      ("n_clones", "rate", "s5f_model", "selection_strength", "allow_extinction")})
    exp.pcr_amplify(count=1).sequencing_errors(count=1).productive_only()
    return ref, exp


def _worker(mode, case, seed, artifact):
    from GenAIRR.experimental import FamilyReplayBundle
    result = {"mode": mode, "case": case, "seed": seed, "rss_before_setup": _rss(),
              "cpu_count": os.cpu_count(), "load_average_before": list(os.getloadavg())}
    if mode in ("plain", "capture"):
        start = time.perf_counter()
        ref, exp = _experiment(case)
        result["setup_seconds"] = time.perf_counter() - start
        start = time.perf_counter()
        compiled = exp.compile()
        result["compile_seconds"] = time.perf_counter() - start
        gc.collect()
        result["rss_before_run"] = _rss()
        start = time.perf_counter()
        simulated = compiled.run_records(seed=seed, strict=COMMON["strict"],
                                         validate_records=COMMON["validate_records"],
                                         capture_family_replay=mode == "capture")
        result["run_seconds"] = time.perf_counter() - start
        result["rss_after_run"] = _rss()
        result.update(_truth(simulated))
        assert all(record["productive"] is True for record in simulated)
        if mode == "capture":
            payload = json.loads(simulated._family_replay_seal)
            attempts = [attempt for clone in payload["clones"] for attempt in clone["attempts"]]
            result.update({
                "attempts": len(attempts),
                "rejected_attempts": len(attempts) - COMMON["n_clones"],
                "all_attempt_nodes": sum(len(attempt["nodes"]) for attempt in attempts),
                "all_attempt_edges": sum(len(attempt["edges"]) for attempt in attempts),
                "growth_proposals": sum(len(attempt["growth"]) for attempt in attempts),
                "ordered_sample_choices": sum(len(attempt["samples"]) for attempt in attempts),
                "attempt_terminations": [attempt["termination"] for attempt in attempts],
            })
            del payload, attempts
            gc.collect()
            result["rss_before_save"] = _rss()
            start = time.perf_counter()
            simulated.save_family_replay(artifact)
            result["save_seconds"] = time.perf_counter() - start
            result["artifact_bytes"] = artifact.stat().st_size
            result["rss_after_save"] = _rss()
    else:
        gc.collect()
        result["rss_before_parse"] = _rss()
        start = time.perf_counter()
        bundle = FamilyReplayBundle.read_from(artifact)
        result["parse_seconds"] = time.perf_counter() - start
        result["rss_after_parse"] = _rss()
        result["artifact_bytes"] = artifact.stat().st_size
        if mode == "replay":
            result["rss_before_replay"] = _rss()
            start = time.perf_counter()
            restored = bundle.replay_records()
            result["replay_seconds"] = time.perf_counter() - start
            result["rss_after_replay"] = _rss()
            result.update(_truth(restored))
            assert restored.validate_records(bundle.refdata).ok
            assert restored.validate_families_with_parents(bundle.refdata).ok
    result["identity"] = _identity()
    result["load_average_after"] = list(os.getloadavg())
    result["benchmark_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result


def run(repeats, seed, selected_cases):
    from GenAIRR.experimental import FamilyReplayLimits
    identity = _identity()
    benchmark_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    measurements = []
    with tempfile.TemporaryDirectory(prefix="genairr-durable-benchmark-") as directory:
        for case in selected_cases:
            for repeat in range(repeats):
                artifact = Path(directory) / f"{case}-{repeat}.json"
                group = {}
                for mode in ("plain", "capture", "parse", "replay"):
                    command = [sys.executable, "-m", "tools.metadata_biology_audit.durable_family_benchmark",
                               "--worker", "--mode", mode, "--cases", case, "--seed", str(seed + repeat),
                               "--artifact", str(artifact)]
                    completed = subprocess.run(command, cwd=ROOT, env=env, text=True, capture_output=True,
                                               timeout=180, check=True)
                    row = json.loads(completed.stdout)
                    row["repeat"] = repeat
                    for field in ("source_sha256", "extension_sha256"):
                        if row["identity"][field] != identity[field]:
                            raise RuntimeError(f"benchmark inputs changed during measurements: {field}")
                    if row["benchmark_sha256"] != benchmark_hash:
                        raise RuntimeError("benchmark script changed during measurements")
                    group[mode] = row
                    measurements.append(row)
                for field in ("records_sha256", "trees_sha256", "accepted_nodes", "observed_records", "sample_draws"):
                    assert group["plain"][field] == group["capture"][field] == group["replay"][field], field
    if "survival_retries" in selected_cases:
        assert any(row.get("rejected_attempts", 0) > 0 for row in measurements if row["case"] == "survival_retries"), "retry control needs a seed with extinction"
    summaries = []
    for case in selected_cases:
        rows = [row for row in measurements if row["case"] == case]
        summary = {"case": case, "parameters": {**COMMON, **CASES[case]}, "repeats": repeats}
        for mode, fields in {"plain": ("compile_seconds", "run_seconds"),
                             "capture": ("compile_seconds", "run_seconds", "save_seconds", "artifact_bytes", "attempts", "rejected_attempts", "all_attempt_nodes", "accepted_nodes"),
                             "parse": ("parse_seconds",), "replay": ("replay_seconds",)}.items():
            subset = [row for row in rows if row["mode"] == mode]
            for field in fields:
                summary[f"{mode}_median_{field}"] = statistics.median(row[field] for row in subset)
            peaks = [snapshot["peak_kib"] for row in subset for name, snapshot in row.items() if name.startswith("rss_")]
            summary[f"{mode}_max_process_peak_kib"] = max(peaks)
        summary["capture_run_over_plain_ratio"] = summary["capture_median_run_seconds"] / summary["plain_median_run_seconds"]
        summaries.append(summary)
    return {"schema": "durable_family_benchmark_v1", "identity": identity,
            "benchmark_sha256": benchmark_hash, "synthetic_lengths_nt": {"v": 300, "d": 18, "j": 51, "molecule": 369},
            "resource_limits": json.loads(FamilyReplayLimits()._json()),
            "repair_support": {"authored_histories": 2, "weights": [1.0, 1.0], "physical_resection_restored": "one top-strand V nucleotide"},
            "checks": ["fresh process per mode/repeat", "identical seeded plain/captured records and trees", "artifact-only fresh-process semantic replay", "replayed records/tree/abundance identity", "all rejected attempts retained"],
            "scope": "Three small synthetic full-length workloads; no empirical repair calibration, production-throughput claim, hard RSS bound, or incremental native-allocation claim. Capture includes sealing; save includes encoding and atomic fsync publication; parse is bounded structural/integrity loading; replay includes reconstructing plans, all failed attempts, accepted trees, and observed reads.",
            "summaries": summaries, "measurements": measurements}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20261022)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--mode", choices=("plain", "capture", "parse", "replay"), help=argparse.SUPPRESS)
    parser.add_argument("--artifact", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.repeats <= 0:
        parser.error("repeats must be positive")
    if args.worker:
        print(json.dumps(_worker(args.mode, args.cases[0], args.seed, args.artifact), allow_nan=False))
    else:
        report = run(args.repeats, args.seed, args.cases)
        text = json.dumps(report, indent=2, allow_nan=False) + "\n"
        if args.output:
            args.output.write_text(text)
        else:
            print(text, end="")


if __name__ == "__main__":
    main()
