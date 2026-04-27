#!/usr/bin/env python3
"""
bench_e2e.py — Python benchmark equivalent for comparison with C.

Usage: python3 benchmarks/bench_e2e.py [N]
  N = number of sequences to generate (default 10000)

Measures wall-clock time for the same pipeline as bench_e2e.c:
  Phase 1: Load DataConfig + compile
  Phase 2: N × simulate (rearrange + S5F mutate)
  Phase 3: Convert to AIRR dicts
"""

import sys
import time
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from GenAIRR import Experiment


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10000

    print("GenAIRR Python End-to-End Benchmark")
    print("====================================")
    print(f"Generating {N} sequences (HUMAN_IGH + S5F)\n")

    # Phase 1: Setup
    t0 = time.perf_counter()
    from GenAIRR.ops import rate, model
    sim = (Experiment.on("human_igh")
           .mutate(rate(0.05, 0.15), model("s5f"))
           .compile(seed=42))
    t_setup = time.perf_counter() - t0
    print(f"Phase 1 — Compile:          {t_setup*1000:8.3f} ms")

    # Phase 2: Generate N sequences
    t1 = time.perf_counter()
    result = sim.simulate(n=N)
    records = list(result)
    t_generate = time.perf_counter() - t1

    print(f"Phase 2 — Generate:         {t_generate*1000:8.3f} ms  ({len(records)} sequences)")
    print(f"  Per sequence:             {t_generate/len(records)*1e6:8.3f} µs")
    print(f"  Throughput:               {len(records)/t_generate:8.0f} seq/sec")

    # Compute stats
    lengths = [len(r.get('sequence', '')) for r in records]
    mut_rates = [r.get('mutation_rate', 0) for r in records]
    productive = sum(1 for r in records if r.get('productive', False))

    avg_len = sum(lengths) / len(lengths)
    avg_rate = sum(mut_rates) / len(mut_rates)

    print(f"  Avg length:               {avg_len:8.1f} bp")
    print(f"  Avg mutation rate:        {avg_rate:8.4f}")
    print(f"  Productive:               {productive:8d} / {len(records)} ({100*productive/len(records):.1f}%)")

    # Phase 3: Write TSV
    t2 = time.perf_counter()
    outpath = "/tmp/genairr_bench_py.tsv"
    import csv
    if records:
        keys = records[0].keys()
        with open(outpath, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
            w.writeheader()
            w.writerows(records)
    t_write = time.perf_counter() - t2
    print(f"\nPhase 3 — Write TSV:        {t_write*1000:8.3f} ms")

    t_total = t_setup + t_generate + t_write
    print(f"\nTotal wall time:            {t_total*1000:8.3f} ms")
    print(f"  Compile:                  {100*t_setup/t_total:8.1f}%")
    print(f"  Generate:                 {100*t_generate/t_total:8.1f}%")
    print(f"  Write:                    {100*t_write/t_total:8.1f}%")

    fsize = os.path.getsize(outpath) if os.path.exists(outpath) else 0
    print(f"\nOutput: {outpath} ({fsize/1e6:.1f} MB, {len(records)} rows)")


if __name__ == "__main__":
    main()
