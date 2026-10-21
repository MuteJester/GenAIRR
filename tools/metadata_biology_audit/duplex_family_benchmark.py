"""Isolated duplex lineage timing, retained-memory and observed-node replay controls.

Run only after building the intended release extension and making the CPU quiet:
  .venv/bin/python -m tools.metadata_biology_audit.duplex_family_benchmark \
    --output audit-docs/duplex_family_benchmark_2026_10_20.json

The harness uses only the standard library beyond GenAIRR itself. Each depth and
repeat gets a separate process. RSS is process memory, not evidence that any
particular native founder allocation is shared through Arc.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gc
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
PARAMETERS = {
    "n_clones": 2,
    "n_max": 32,
    "n_sample": 16,
    "lambda_base": 3.0,
    "rate": 0.15,
    "selection_strength": 0.0,
    "productive_only": True,
    "strict": True,
    "validate_records": True,
}


def _rss():
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if platform.system() == "Darwin":
        peak /= 1024
    current = None
    if platform.system() == "Linux":
        resident_pages = int(Path("/proc/self/statm").read_text().split()[1])
        current = resident_pages * os.sysconf("SC_PAGE_SIZE") / 1024
    return {"current_kib": current, "peak_kib": peak}


def _identity():
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()

    paths = sorted([
        *ROOT.joinpath("engine_rs/src").rglob("*.rs"),
        *ROOT.joinpath("src/GenAIRR").rglob("*.py"),
        ROOT / "engine_rs/Cargo.toml", ROOT / "engine_rs/Cargo.lock",
    ])
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.relative_to(ROOT)).encode() + b"\0")
        digest.update(path.read_bytes() + b"\0")
    from GenAIRR import _engine
    extension = Path(_engine.__file__).resolve()
    return {
        "commit": git("rev-parse", "HEAD"),
        "branch": git("branch", "--show-current"),
        "dirty": bool(git("status", "--porcelain")),
        "source_sha256": digest.hexdigest(),
        "source_hash_scope": "engine_rs/src/**/*.rs,src/GenAIRR/**/*.py,engine_rs/Cargo.toml,Cargo.lock",
        "benchmark_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "extension_path": str(extension),
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
        "python": sys.version,
        "platform": platform.platform(),
    }


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _worker(depth, seed):
    process_baseline = _rss()
    import GenAIRR as ga
    from GenAIRR.experimental import CoupledDuplexRepair, DuplexHistory, DuplexJoin, DuplexRepairModel

    refdata = ga.RefDataConfig.vdj()
    refdata.add_v_allele("V*01", "V", b"AAATGT", anchor=3)
    refdata.add_d_allele("D*01", "D", b"CCC")
    refdata.add_j_allele("J*01", "J", b"TTTAAA", anchor=0)
    repair = CoupledDuplexRepair("AAATGT", "CCC", "TTTAAA")
    repair = repair.join(
        DuplexJoin("vd").open("left", offset=0).open("right", offset=0)
        .resect("left", "top", count=1).synapse(right_start=6).fill("top").ligate()
    ).join(
        DuplexJoin("dj").open("left", offset=0).open("right", offset=0)
        .synapse(right_start=3).ligate()
    )
    histories = [
        DuplexHistory(name, repair, "V*01", "D*01", "J*01")
        for name in ("restored-a", "restored-b")
    ]
    authored = DuplexRepairModel.bind(
        refdata, histories, reference_ends_are_coding_ends=True,
    )
    experiment = (
        ga.Experiment.on(refdata).recombine(duplex_repair=authored)
        .clonal_lineage(
            n_clones=PARAMETERS["n_clones"], max_generations=depth,
            n_max=PARAMETERS["n_max"], n_sample=PARAMETERS["n_sample"],
            lambda_base=PARAMETERS["lambda_base"], rate=PARAMETERS["rate"],
            selection_strength=PARAMETERS["selection_strength"],
        ).productive_only()
    )
    gc.collect()
    before_compile = _rss()
    start = time.perf_counter()
    compiled = experiment.compile()
    compile_seconds = time.perf_counter() - start
    after_compile = _rss()

    start = time.perf_counter()
    result = compiled.run_records(
        seed=seed, strict=PARAMETERS["strict"],
        validate_records=PARAMETERS["validate_records"],
    )
    run_records_seconds = time.perf_counter() - start
    after_run = _rss()
    assert len(result.parents) == len(result.families) == len(result.lineage_trees) == 2
    assert all(record["productive"] is True for record in result)
    assert all(record["p_v_3_length"] is None for record in result)
    assert sum(record["duplicate_count"] for record in result) == 32

    native_replay_seconds = 0.0
    replay_check_seconds = 0.0
    replayed_records = 0
    start = time.perf_counter()
    for record in result:
        family = result.families[record["parent_id"]]
        replay_start = time.perf_counter()
        replayed = family.replay_node(record["lineage_node_id"], refdata)
        native_replay_seconds += time.perf_counter() - replay_start
        check_start = time.perf_counter()
        assert replayed.validate_record(refdata) == []
        replay_record = replayed.airr_record(refdata, sequence_id=record["sequence_id"])
        # The ordinary wrapper also stamps clone/node abundance metadata.
        # Compare every native generated field against its actual emitted value.
        emitted_native = {field: record[field] for field in replay_record}
        assert _canonical(replay_record) == _canonical(emitted_native)
        replay_check_seconds += time.perf_counter() - check_start
        replayed_records += 1
    replay_seconds = time.perf_counter() - start
    after_replay = _rss()
    nodes_by_family = [len(tree.nodes()) for tree in result.lineage_trees]
    observed_by_family = [sum(node.observed for node in tree.nodes())
                          for tree in result.lineage_trees]
    selected = Counter(parent.airr_record(refdata)["duplex_repair"]["selection"]["history_id"]
                       for parent in result.parents)
    assert sum(observed_by_family) == len(result) == replayed_records
    return {
        "max_generations": depth,
        "seed": seed,
        "compile_seconds": compile_seconds,
        "run_records_seconds": run_records_seconds,
        "replay_seconds": replay_seconds,
        "native_replay_seconds": native_replay_seconds,
        "replay_check_seconds": replay_check_seconds,
        "process_baseline_rss": process_baseline,
        "before_compile_rss": before_compile,
        "after_compile_rss": after_compile,
        "after_run_records_rss": after_run,
        "after_replay_rss": after_replay,
        "compile_peak_increase_kib": max(0, after_compile["peak_kib"] - before_compile["peak_kib"]),
        "run_records_peak_increase_kib": max(0, after_run["peak_kib"] - after_compile["peak_kib"]),
        "replay_peak_increase_kib": max(0, after_replay["peak_kib"] - after_run["peak_kib"]),
        "actual_total_nodes": sum(nodes_by_family),
        "nodes_by_family": nodes_by_family,
        "observed_records": len(result),
        "observed_records_by_family": observed_by_family,
        "sampled_cell_draws": sum(record["duplicate_count"] for record in result),
        "selected_histories": dict(selected),
        "replayed_records": replayed_records,
        "replay_errors": 0,
        "serialized_record_bytes": sum(len(_canonical(record).encode()) for record in result),
        "identity": _identity(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depths", nargs="+", type=int, default=[2, 4, 6])
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20261020)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.repeats <= 0 or any(depth <= 0 or depth > 1000 for depth in args.depths):
        parser.error("repeats must be positive and depths must be between 1 and 1000")
    if args.worker:
        print(json.dumps(_worker(args.depths[0], args.seed)))
        return

    identity = _identity()
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    measurements = []
    for depth in args.depths:
        for repeat in range(args.repeats):
            command = [
                sys.executable, str(Path(__file__).resolve()), "--worker",
                "--depths", str(depth), "--seed", str(args.seed + repeat),
            ]
            completed = subprocess.run(
                command, cwd=ROOT, env=env, text=True, capture_output=True,
                check=True, timeout=120,
            )
            measurements.append(dict(json.loads(completed.stdout), repeat=repeat))
    final_identity = _identity()
    for entry in [final_identity, *(measurement["identity"] for measurement in measurements)]:
        for field in ("source_sha256", "extension_sha256", "benchmark_sha256"):
            if entry[field] != identity[field]:
                raise RuntimeError(f"benchmark inputs changed during measurements: {field}")

    summaries = []
    for depth in args.depths:
        group = [measurement for measurement in measurements
                 if measurement["max_generations"] == depth]
        summary = {"max_generations": depth, "repeats": len(group)}
        for field in ("compile_seconds", "run_records_seconds", "replay_seconds",
                      "native_replay_seconds", "replay_check_seconds"):
            summary[f"median_{field}"] = statistics.median(measurement[field] for measurement in group)
        for field in ("compile_peak_increase_kib", "run_records_peak_increase_kib",
                      "replay_peak_increase_kib"):
            summary[f"max_{field}"] = max(measurement[field] for measurement in group)
        for field in ("actual_total_nodes", "observed_records"):
            summary[f"min_{field}"] = min(measurement[field] for measurement in group)
            summary[f"max_{field}"] = max(measurement[field] for measurement in group)
        summaries.append(summary)
    result = {
        "model": "duplex_lineage_isolated_benchmark_v1",
        "identity": identity,
        "parameters": PARAMETERS,
        "source": {
            "sequence": "AAATGTCCCTTTAAA",
            "history_ids": ["restored-a", "restored-b"],
            "history_weights": [1.0, 1.0],
            "restored_v_top_bases": 1,
            "reference_bound": True,
        },
        "measurement": {
            "process_isolation": "one sequential subprocess per depth and repeat",
            "compile": "normal Experiment.compile, after model/DSL construction",
            "run_records": "normal compiled lineage run_records, full metadata export and all record/family validation",
            "replay": "each emitted observed node; native replay, validation, metadata export and exact canonical JSON comparison",
            "memory": "current resident RSS from Linux procfs plus monotonic process RSS high-water marks; KiB",
        },
        "limitations": [
            "Short synthetic 15-base founders do not establish performance for realistic full-length references or larger histories.",
            "Equal authored history weights are unfitted controls, not measured biochemical rates.",
            "Neutral offspring and mutation settings are declared benchmark controls, not an empirical family model fit.",
            "All emitted metadata and native family replay state remain retained during measurements; this is not a bounded-stream benchmark.",
            "Process RSS, allocator retention and high-water marks cannot prove individual immutable founder Arc sharing; no allocation introspection is performed.",
            "Replay timing includes validation and metadata equality work; native replay time is also reported separately.",
            "Two repeats at each depth are a small descriptive sample, not a statistical latency guarantee.",
        ],
        "summaries": summaries,
        "measurements": measurements,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
        print(json.dumps({"output": str(args.output), "summaries": summaries}, indent=2))
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
