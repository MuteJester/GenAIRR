# Durable family replay: full-length synthetic benchmark

All six seeded controls passed: uncaptured, captured and artifact-only replay
produced identical typed records and tree exports. The 24 isolated subprocesses
covered 16 attempted families, including 10 rejected survival attempts, with
101 attempted nodes and 79 nodes in accepted families. Default resource limits
passed without overrides. [Raw measurements](durable_family_benchmark_2026_10_21.json)
include every stage, seed, node count, counter and identity.

The synthetic coding molecule contains a 300 nt V, 18 nt D and 51 nt J. Two
equally weighted authored histories restore one resected top-strand V nucleotide.
Runs use neutral growth, the existing `hh_s5f` kernel at rate 0.01, productive
contracts, strict execution, validation and one PCR/one sequencing error per
observed cell. Each family samples two cells with replacement; repeated draws
can produce one emitted record. These controls do not fit biochemical rates or
represent an empirical repertoire.

Each case used seeds 20261022 and 20261023 in separate processes for uncaptured
run, capture/publication, parse and replay. The shallow/deeper cases use carrying
capacity 8, offspring mean 2.4 and depths 2/6; the retry case uses mean 0.7 and
depth 3. The retry seeds required 7/2 attempts, including 6/1 rejected attempts.

| Case | Plain run, s | Captured run, s | Save, s | Parse, s | Replay, s | File, MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Shallow | 0.258 | 1.703 | 3.116 | 0.420 | 2.795 | 6.30 |
| Deeper | 0.272 | 1.967 | 2.509 | 0.443 | 3.211 | 6.37 |
| Survival retries | 0.190 | 1.219 | 2.029 | 0.421 | 2.126 | 6.27 |

Values are medians of two realizations. Compilation took approximately 4–5 ms.
Captured run time includes sealing; save includes checking, content encoding,
commitments and atomic publication with file/directory fsync. Parse includes
bounded loading and integrity checks. Replay time excludes that initial parse
but includes fresh plan/reference reconstruction, all rejected attempts, accepted
trees and observed read stages. Python startup/import time is excluded.

Maximum whole-process high-water RSS measured immediately after each stage:

| Case | Plain run, MiB | Captured run, MiB | After save, MiB | Parse, MiB | After replay, MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| Shallow | 65.4 | 201.0 | 249.2 | 169.6 | 284.0 |
| Deeper | 65.8 | 200.3 | 246.2 | 171.0 | 285.2 |
| Survival retries | 56.0 | 133.9 | 225.1 | 168.9 | 195.5 |

These are process peaks, not incremental native allocations, retained family
size or a hard RSS guarantee. Save peaks include preceding benchmark checks;
replay peaks include parsing. Baseline process peaks were approximately 32 MiB.
The current JSON codec has substantial time and memory cost: capture took
6.4–7.2 times the uncaptured run time, with publication costing another 2–3 s.
This supports a future codec/copy optimization milestone, without changing the
ordinary simulation path or weakening replay evidence.

Measurements ran on a shared 16-CPU development host after GenAIRR verification
loads finished; recorded one-minute loads ranged approximately 1.42–1.59. OS
caches were shared. Two runs are descriptive observations, not throughput,
scaling, latency or hardware-independent guarantees.

The release extension SHA256 is
`a31be980a932d7333714f2ddf51c78ff9296943450d2d105313638bcff482244`.
The complete engine/Python source scope recorded in the JSON hashes to
`bb602141a97971330c88db8e6457df8328ffe573a13cd9960946ff06ea7aa8f3`;
the benchmark script hashes to
`09b3df1de490fbcb9bdd14850a7f75707db870ca4b4a196644fe218941dc828f`.
The worktree was uncommitted during measurement, so these content identities
bind the measured implementation more precisely than the preceding commit.

Reproduce with:

```bash
.venv/bin/python -m tools.metadata_biology_audit.durable_family_benchmark \
  --repeats 2 --output audit-docs/durable_family_benchmark_2026_10_21.json
```
