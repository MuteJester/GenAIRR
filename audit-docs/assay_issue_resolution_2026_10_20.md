# Assay issue resolution

Branch: `fix/biological-metadata-correctness`. Reviewed 2026-10-20.

Subsequent approved work adds an opt-in original-coding-end outcome model;
see [its implementation and verification](original_coding_end_processing_2026_10_20.md).
The assay results below describe the preceding fixes. The legacy post-trim P
mode remains available; the new mode uses original coding-tip templates and
joint net processing outcomes with persistent base ancestry.

This follows the twelve findings in the AlignAIR assay report for baseline
commit `1ededa6`. The report was found in the adjacent AlignAIR worktree;
that worktree's source was left unchanged. The baseline observations remain
valid historical evidence. This document describes the changed behavior.

| Assay issue | Resolution |
| --- | --- |
| 1. End loss becomes alignment deletion | Explicit clipping now updates observed reference coverage. Missing terminal coverage creates no `D` or identity penalty. Structural deletions inside coverage remain recorded. Fully unobserved segments have no germline bounds, including an already-deleted segment crossed by clipping. |
| 2. Inverted D uses forward walker | D alignments traverse the reverse complement with forward-allele coordinate bounds. Strings, identity, CIGAR and real deletions share the same orientation. Streaming and fresh re-call scoring also share orientation-aware trim limits. |
| 3. V alignment overrun | Projection uses the generating reference and its retained interval. Random NP bases cannot extend the truth alignment beyond its allele. |
| 4. Post-SHM re-call mixes frames | As requested, standard `v_call`, `d_call`, and `j_call` now identify the final generating alleles. `recalled_v_call`, `recalled_d_call`, and `recalled_j_call` preserve post-mutation inference separately. All truth alignment fields use the generating reference. Donor re-calls consider the union of alleles carried across both chromosomes. |
| 5. Short re-call loses source deletion | A re-call's length no longer limits the truth CIGAR. Source deletions remain represented, including terminal deletions beyond a shorter inferred allele. |
| 6. NP length one short | Generating N/P boundaries remain intact when an inserted base happens to match neighboring trimmed germline. The combined NP interval includes intervening P bases. |
| 7. Missing nonproductive junction AA | `junction_aa` translates complete codons from the physical junction's local start even when it is nonproductive or out of frame. Trailing one or two nucleotides are omitted. Productivity remains a separate annotation. |
| 8. Read-length mutation denominator | `mutation_rate` sums each SHM stage's realized events divided by its input molecule length. Read loss cannot inflate previous contributions. Repeated hits and reversions remain events; this is not an observed-read mismatch fraction. Lineage keeps its net-difference summary divided by pre-observation node length. |
| 9. D overtrim panic | Statically impossible finite trim support raises a compile error; an impossible realized trim returns a structured execution error in either strict or permissive runs. No silent clamping. |
| 10. Optional pandas required by kernels | The two affected bundled artifacts now contain plain dictionaries. All 5,120 numeric kernel values per model are bit-identical to the originals. A subprocess blocks both pandas and numpy while loading and simulating all four models. Serialized artifact hashes intentionally changed. |
| 11. No P model/control | Authorable per-end P-length distributions and productive-only nonzero-P support are implemented. Defaults remain zero without an empirical P plane. No calibrated human-IGH frequencies were invented. The existing post-trim palindrome mechanism remains a modeling limitation, described below. |
| 12. Description/version | D inversion has a proper description. The rebuilt editable package and source both report `2.3.1`; stale installed version metadata required reinstalling rather than a source version change. |

The call-field change is intentional and affects consumers of the old inference
surface. Standard calls are single generating alleles. Re-called fields may
contain comma-separated ties. `expose_provenance=True` retains `truth_*_call`
as compatibility aliases of the standard calls. Receptor revision uses the
replacement allele, with `original_v_call` preserving its predecessor.
Without observed inference, re-called fields are empty rather than copying the
generating identity into an evidence field.

Clonal continuation outcomes now carry ancestral events, revisions and trace
once, so inherited mutation counters remain independently checkable. Lineage
outcomes retain their synthesized net-SHM summary rather than a full chronological
mutation history. Appended lineage corruption events now point to their actual
choice spans in the combined trace.

Verification includes source-position reconstruction for 800 high-SHM records
with indels and both D orientations, 1,000 exact-NP-length records, explicit
clipping/export cases, donor-restricted inference and true replay, model loading
without optional dependencies, and mutation burden across changing molecule
lengths. The native tests additionally cover source deletions beyond a shorter
re-call and clipping across an already-empty segment.

P additions still use the reverse complement of a retained, trimmed coding
flank. This does not claim to implement the biochemical sequence of hairpin
opening and exonuclease processing. Default cartridges do not supply calibrated
P frequencies. Combined NP strings also remain unsuitable for N-only cartridge
estimation when explicit P provenance is nonzero; those estimators reject such
inputs until a trustworthy N-only representation is available.

Passing these checks demonstrates the stated invariants and implemented-model
consistency. It is not a guarantee of universal biological realism.

## Verification

- Full Python correctness suite: **2,579 passed**, eight performance tests
  deselected, 168 warnings, 477.31 seconds. Log:
  `/tmp/genairr_assay_python_final.log`.
- Full Rust suites: **1,445 debug tests and 1,444 release tests passed**.
  Each profile retains two ignored integration tests and two ignored doctests;
  one panic regression is debug-only.
- The [external assay rerun](assay_oracle_verification_2026_10_20.json) rebuilt
  **all 13 dev suites at 2% scale** against consistent simulator and assay
  sources: **638/638 oracle cases passed**, including donor-aware scoring,
  with no exclusions, warning cases or unscorable segments. This is an oracle
  consistency check, not an external aligner accuracy benchmark.
- The [coordinate sweep](biological_metadata_fix_coordinate_sweep_2026_10_20.json)
  checked **1,200 IGH/IGK records with zero errors** across twelve pipeline
  profiles.
- Repeating the original high-mutation comparison (1,000 IGH records, seed 11,
  SHM 0.15, two indels) gave **zero standard/source call disagreements**.
  Re-calls still excluded the source in V/D/J on **71/53/2 records**, exactly
  matching the pre-change comparison. All records passed validation.
- The [original diagnostic probes](biological_metadata_fix_verification_2026_10_20.json)
  now reproduce two of nineteen observations: the post-trim P modeling question
  and segment-local native CIGAR strings. The native-format probe still flags
  missing AIRR placement prefixes in Python dictionaries; strict TSV export
  supplies those prefixes and passes its separate check. Internal dictionaries
  retain their documented zero-based, segment-local representation.
- A final validator-only adjustment accepts absent inference in direct native
  partial outcomes and rejects invented re-calls there. The full native suites
  and **53 Python validator/assay regression tests** passed after rebuilding it.
  This adjustment followed collection of the full Python run; that run and the
  subsequent focused verification are reported separately rather than claiming
  one full run against a later binary.
- `git diff --check` passed.
