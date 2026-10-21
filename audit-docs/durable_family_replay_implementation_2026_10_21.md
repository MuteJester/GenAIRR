# Durable duplex lineage replay

The [approved design](durable_family_replay_design_2026_10_21.md) is implemented for source-backed BCR VDJ `clonal_lineage` experiments, including multiple clones. Capture, persistence and reconstruction run entirely inside GenAIRR. No simulator, executable archive, pickle loader or new dependency is required.

## Public workflow

Enable `capture_family_replay=True` in `Experiment.run_records` or the corresponding compiled lineage entry point. Optional `family_replay_limits=FamilyReplayLimits(...)` supplies computational ceilings. Save with `result.save_family_replay(path)`. In a fresh process, `FamilyReplayBundle.read_from(path).replay_records()` reconstructs ordinary native outcomes, parents, families, lineage trees and nullable AIRR dictionaries. The verified bundle exposes its reconstructed `refdata`; `replay_node(clone_id=..., node_id=...)` reconstructs retained cellular states, including roots and unsampled nodes, before read corruption.

The [executable guide](../site_docs/guides/experimental-duplex-repair.md#saving-and-replaying-a-complete-lineage) contains the complete example. The APIs are exported from `GenAIRR.experimental` and the artifact declares `genairr_family_replay`, schema 1, workflow `duplex_lineage_v1`.

## Scientific evidence and replay

Each clone has exactly one source-backed repaired founder. Every growth attempt is captured, including rejected extinct attempts and its own automatically generated target. Receipts preserve explicit/automatic/neutral target policy, offspring proposals even when capacity suppresses children, actual incoming-edge mutation choices, final live populations, ordered sample draws, abundance and accepted-attempt addressing. Actual mutation transitions include repeated substitutions and reversions. Accepted extinct clones retain their founder and empty observed family.

Observed cells execute the recorded read plan once in cell arena order; repeated samples of that cell contribute abundance. Unmerged read receipts retain PCR and quality substitutions, ambiguous bases, indels, terminal losses and reversal. Reconstruction reuses native constructors, the physical repair/projection contracts, division replay and checked continuation composition. Native record, family and parent validators run before results are returned. The final requested node reconstruction shares the same caller budget as verification of the complete bundle.

Replay consumes recorded choices under the recorded execution policy with a native fence forbidding fresh RNG words. It validates positive effective support, control flow, event/state commitments, topology, counts, affinity, ancestry and exact output. This establishes consistency with the declared generating process. It does not establish seed-stream provenance, authorship or empirical biological frequencies.

## Durable execution context

Closed versioned descriptors preserve the effective ordered reference catalogue and allele IDs, sequence-bound evidence, reference rules and identity, contracts, productivity policy, finite repair vocabulary, exact native distribution CDFs and numeric mutation kernels. Reconstruction uses ordinary checked native plan compilation. Unsupported passes and schedules that change when explicit ordering edges are omitted fail before captured sampling. External cartridges and built-in mutation-kernel loaders are unnecessary on replay.

A new length-framed, type-tagged SHA-256 identity prevents delimiter aliases in names and metadata. Existing trace v4 bytes, hashes, signatures and ordinary source/continuation gates retain their meanings. Legacy traces embedded inside the new format receive additional closed-field validation without changing their historical schemas.

Binary64 parameters and dynamic output values use exact bits. A qualified architecture/OS/numeric-probe profile rejects incompatible execution environments explicitly. Recorded package/engine versions and extension SHA-256 are informative build provenance, rather than a requirement to load the identical binary. JSON output tagging preserves nulls, Boolean/integer/float distinctions, signed zero, tiny finite floats and literal hexadecimal strings; unbounded integers and nonfinite values fail.

## Sealing, resources and publication

Capture seals the actual output after the selected entry point applies its metadata. Direct compiled runs and Experiment runs retain their existing metadata behavior. Nested annotation snapshots are independent of user containers. Saving rejects record changes, reordering, removed rows and replacement native evidence. Earlier results generated without capture cannot manufacture missing attempt evidence.

The strict reader rejects duplicate keys, integer overflow, unknown fields/variants, unresolved or cyclic component references, surplus components and incompatible capabilities. A content DAG deduplicates structured values. Iterative inspection computes expanded depth/value/byte costs before materialization, including repeated references. Native accounting aggregates source construction, preflight, attempts, choices, candidate scans, pool copies, retained branch histories and commitments; output rows have cumulative work/byte fences. Exhaustion returns a structured error without changing the sampled process. These are software-accounted limits, not an absolute process-RSS guarantee.

Captured growth currently requires `lambda_base * (1 + selection_strength) <= 64`; this conservative computational restriction does not alter uncaptured sampling or represent a fitted biological prior. Extreme nonfinite S5F normalization is rejected during bounded capture rather than producing an unreplayable no-op. Permissive legal no-ops retain their original policy during replay, and impossible counts or zero-weight choices are rejected.

Publication uses a same-directory temporary file, file flush/fsync, race-safe no-clobber linking by default or explicit atomic replacement, then directory fsync on supported systems. Prepublication failures preserve the previous target and remove temporary files. If directory sync fails after publication, the complete new file may already be visible and the error reports uncertain crash persistence; no partial bundle is published.

## Verification

The [source/binary-pinned verification artifact](durable_family_replay_verification_2026_10_21.json) records commands, complete suite results and changed-source hashes. The [independent durable sweep](durable_family_sweep_2026_10_21.json) compares captured and uncaptured generation, saved/reloaded native truth, literal ancestry/edit/counter reconstruction, trees and all retained cells. Existing family, projection and observation sweeps are rerun against the final extension as separate dated artifacts; earlier milestone evidence remains historical.

All 3,138 nonperformance Python tests pass, including 143 durable replay tests.
The final full run uses two disjoint file groups, checking runtime source and
extension hashes before and after execution. The 135 focused compatibility
tests also pass. Native debug and release suites pass 1,663 and 1,662 tests
respectively, with three inherited ignored cases. All eight existing
performance checks pass under unchanged thresholds in a separate run. The
strict documentation build and all seven actual Python guide examples pass.

The durable sweep covers 52 plans and 104 founders, including 18 rejected
attempts and four accepted extinct families. It checks 298 observed records and
1,354 cellular states with zero errors. The compatibility sweeps also report
zero errors across 96 existing family plans, 32,768 projection cases and
18,432 observation continuations from 2,048 founders.

The [isolated full-length benchmark](durable_family_benchmark_2026_10_21.md), with its [raw measurements](durable_family_benchmark_2026_10_21.json), reports capture, publication, parsing and replay separately, including file sizes, process peak RSS and retry/deep-history controls. The JSON/evidence-heavy first implementation can incur substantial memory and runtime overhead; these measurements describe the declared workloads and do not promise streaming or arbitrary family-size performance. Exact measurements and test totals are recorded in the verification artifact.

Its 24 isolated processes cover six paired realizations of a 369-nucleotide
synthetic molecule, depths two and six, and a low-growth retry control. All
plain, captured and replayed record/tree hashes agree under the unchanged
default limits. The captured runs retain 16 attempts, including ten rejected
attempts, and 101 total nodes. Median capture takes 1.22–1.97 seconds, publication
2.03–3.12 seconds and replay 2.13–3.21 seconds. Artifact sizes are approximately
6.27–6.37 MiB; maximum whole-process peak RSS is approximately 285 MiB. These
are descriptive two-repeat measurements on a shared development host.

## Boundaries

V1 excludes fixed expansion/repertoire bundles, legacy lineage, arbitrary external initial snapshots, VJ, phased genotype families, revision, contaminants and paired-end adapters. Existing in-memory APIs for other supported workflows remain available. Additional workflow adapters and persistence migrations require explicit versioned contracts.

Equal authored repair-history weights remain unfitted controls over the supplied vocabulary. Stochastic generation of new repair histories, general flap/mismatch processing, reopening completed joins, enzyme kinetics and empirical biochemical calibration remain separate milestones. Measured data can later fit the native model; this implementation adds no invented biochemical rates.
