"""Isolated finite-duplex compile and bounded-stream benchmark (stdlib only).

Run with the built extension:
  .venv/bin/python -m tools.metadata_biology_audit.duplex_compile_benchmark \
    --output audit-docs/duplex_experiment_compile_benchmark_2026_10_20.json

Each measurement has its own process. Histories differ by stable authored ID;
the physical blunt-join recipe is intentionally shared. RSS is a process high
water mark, not an assertion about individual native Arc allocations.
"""
from __future__ import annotations

import argparse
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


def _rss_kib():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value / 1024 if platform.system() == "Darwin" else value


def _worker(rows, contexts, records, seed):
    process_baseline = _rss_kib()
    import GenAIRR as ga
    from GenAIRR import _engine
    from GenAIRR.experimental import CoupledDuplexRepair, DuplexHistory, DuplexJoin, DuplexRepairModel

    refdata = ga.RefDataConfig.vdj()
    for i in range(contexts):
        refdata.add_v_allele(f"V{i+1}*01", f"V{i+1}", b"AAATGT", anchor=3)
    refdata.add_d_allele("D*01", "D", b"CCC")
    refdata.add_j_allele("J*01", "J", b"TTTAAA", anchor=0)
    script = CoupledDuplexRepair("AAATGT", "CCC", "TTTAAA", max_operations=8, max_bases=64)
    for join, shift in (("vd", 6), ("dj", 3)):
        script = script.join(DuplexJoin(join).open("left", offset=0).open("right", offset=0)
                             .synapse(right_start=shift).ligate())
    histories = [
        DuplexHistory(f"history-{i:06}", script, f"V{i % contexts + 1}*01", "D*01", "J*01", 1.0)
        for i in range(rows)
    ]
    model = DuplexRepairModel.bind(
        refdata, histories, reference_ends_are_coding_ends=True,
        # One natural and one distinct admitted mask per flat V/D/J source.
        max_preflight_operations=contexts + rows * 8 + 2 * (contexts + 2),
        max_preflight_bases=rows * script.max_bases,
    )
    experiment = ga.Experiment.on(refdata).recombine(duplex_repair=model).productive_only()
    gc.collect()
    baseline_rss = _rss_kib()
    start = time.perf_counter()
    compiled = experiment.compile()
    compile_seconds = time.perf_counter() - start
    compile_rss = _rss_kib()
    counts = {}
    first = None
    start = time.perf_counter()
    for outcome in compiled.stream(n=records, seed=seed, strict=True):
        index = outcome.trace().find("duplex_repair.v1.history").value
        counts[index] = counts.get(index, 0) + 1
        if first is None:
            record = outcome.airr_record(refdata)
            assert record["sequence"] == "AAATGTCCCTTTAAA"
            assert record["productive"] is True
            selection = record["duplex_repair"]["selection"]
            assert selection["history_id"] == histories[index].id
            first = {
                "history_id": selection["history_id"],
                "model_digest": selection["model_digest"],
                "sequence": record["sequence"],
                "physical_ledger_bases": len(record["duplex_repair"]["product"]["ledger"]),
            }
    stream_seconds = time.perf_counter() - start
    stream_rss = _rss_kib()
    return {
        "rows": rows, "natural_contexts": contexts, "records": records, "seed": seed,
        "compile_seconds": compile_seconds, "stream_seconds": stream_seconds,
        "process_baseline_rss_kib": process_baseline,
        "before_compile_rss_kib": baseline_rss, "after_compile_rss_kib": compile_rss,
        "after_stream_rss_kib": stream_rss,
        "compile_peak_increase_kib": max(0, compile_rss - baseline_rss),
        "stream_peak_increase_kib": max(0, stream_rss - compile_rss),
        "distinct_histories_sampled": len(counts), "sampled_records": sum(counts.values()),
        "max_preflight_operations": model.max_preflight_operations,
        "max_preflight_bases": model.max_preflight_bases,
        "extension_path": str(Path(_engine.__file__).resolve()), "first_record": first,
    }


def _identity():
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()

    paths = sorted(
        [*ROOT.joinpath("engine_rs/src").rglob("*.rs"),
         *ROOT.joinpath("src/GenAIRR").rglob("*.py"),
         ROOT / "engine_rs/Cargo.toml", ROOT / "engine_rs/Cargo.lock"]
    )
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.relative_to(ROOT)).encode() + b"\0")
        digest.update(path.read_bytes() + b"\0")
    from GenAIRR import _engine
    extension = Path(_engine.__file__).resolve()
    return {
        "commit": git("rev-parse", "HEAD"), "branch": git("branch", "--show-current"),
        "dirty": bool(git("status", "--porcelain")),
        "source_sha256": digest.hexdigest(),
        "source_hash_scope": "engine_rs/src/**/*.rs,src/GenAIRR/**/*.py,engine_rs/Cargo.toml,Cargo.lock",
        "extension_path": str(extension),
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
        "python": sys.version, "platform": platform.platform(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", nargs="+", type=int, default=[8, 32, 128])
    parser.add_argument("--contexts", nargs="+", type=int, default=[1, 8])
    parser.add_argument("--records", type=int, default=1000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20261020)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if any(value <= 0 for value in [*args.rows, *args.contexts, args.records, args.repeats]):
        parser.error("rows, contexts, records, and repeats must be positive")
    if args.worker:
        print(json.dumps(_worker(args.rows[0], args.contexts[0], args.records, args.seed)))
        return
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    measurements = []
    for contexts in args.contexts:
        for rows in args.rows:
            if rows < contexts:
                continue
            for repeat in range(args.repeats):
                command = [
                    sys.executable, str(Path(__file__).resolve()), "--worker",
                    "--rows", str(rows), "--contexts", str(contexts),
                    "--records", str(args.records), "--seed", str(args.seed + repeat),
                ]
                completed = subprocess.run(command, cwd=ROOT, env=env, text=True,
                                           capture_output=True, check=True, timeout=120)
                measurements.append(dict(json.loads(completed.stdout), repeat=repeat))
    summaries = []
    for contexts in args.contexts:
        for rows in args.rows:
            group = [m for m in measurements if m["rows"] == rows and m["natural_contexts"] == contexts]
            if not group:
                continue
            summaries.append({
                "rows": rows, "natural_contexts": contexts, "repeats": len(group),
                "median_compile_seconds": statistics.median(m["compile_seconds"] for m in group),
                "median_stream_seconds": statistics.median(m["stream_seconds"] for m in group),
                "max_compile_peak_increase_kib": max(m["compile_peak_increase_kib"] for m in group),
                "max_stream_peak_increase_kib": max(m["stream_peak_increase_kib"] for m in group),
            })
    result = {
        "model": "finite_authored_duplex_history_v1", "identity": _identity(),
        "measurement": "isolated processes; monotonic process RSS high-water marks; bounded streaming",
        "limitations": [
            "Synthetic 15-base molecules and duplicate physical recipes test vocabulary scaling only.",
            "Reference hashing/binding and physical replay remain per active authored row.",
            "RSS does not prove native Arc sharing; no entire outcome batch is retained.",
            "Timings include Python/native execution and one first-record metadata export.",
        ],
        "summaries": summaries, "measurements": measurements,
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
