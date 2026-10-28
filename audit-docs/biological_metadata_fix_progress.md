# Biological metadata repair progress

Branch: `fix/biological-metadata-correctness`. Review date: 2026-10-20.

The approved [original coding-end processing follow-up](original_coding_end_processing_2026_10_20.md)
adds an opt-in finite original-tip outcome model with persistent ancestry and
shared productive feasibility. The assay counts below precede that addition.

The latest follow-up addresses the twelve AlignAIR assay findings; see the
[issue-by-issue resolution and current verification](assay_issue_resolution_2026_10_20.md).
Verification counts farther below record earlier stages, before this follow-up.

The [baseline audit](biological_metadata_audit_2026_10_20.md) and its
[recorded observations](biological_metadata_audit_2026_10_20.json) describe
commit `1ededa676e8904430f7ea0bdcf03092de4888a41`. They remain historical
evidence; their reproduction counts are not the status of this branch.

## Latest assay follow-up verification

- Standard calls and truth alignments now use generating alleles. Inferred
  calls are separate `recalled_*_call` fields and honor donor genotype support.
- End loss shortens observed reference coverage. Genuine structural deletions
  and generated NP boundaries remain represented correctly.
- Mutation rates preserve mutation-stage lengths, nonproductive junctions
  retain AA translation, explicit overtrims return errors, and all bundled
  S5F kernels load without pandas or numpy.
- **2,579 Python correctness tests passed**, eight performance tests excluded.
  A subsequent native-only missing-inference validator adjustment passed
  **53 focused Python checks** against its rebuilt extension.
- Full native suites after that adjustment: **1,445 debug / 1,444 release**
  tests passed; two integration tests and two doctests remain ignored per profile.
- **638/638 external oracle cases passed** across all thirteen dev suites at
  2% scale, with no exclusions, warnings or unscorable segments.
- The current coordinate sweep has **1,200 records and zero errors**. The
  nineteen original diagnostic probes now reproduce only the post-trim P
  modeling question and the documented segment-local native CIGAR format;
  strict TSV placement passes.

See the [twelve-issue resolution](assay_issue_resolution_2026_10_20.md) and
[external oracle evidence](assay_oracle_verification_2026_10_20.json).

The implementation now addresses these concrete failures:

| Area | Changed behavior |
| --- | --- |
| Coding interpretation | Translation follows the retained V coding origin. Stop reporting covers the retained coding sequence, including stops outside the junction. Productivity checks V/J anchor phase and the J coding endpoint, so internal indels and terminal structural deletions cannot masquerade as missing read coverage. Explicit observation-loss provenance preserves the phase of clipped reads. |
| Productive generation | The productive contract includes coding-sequence integrity. Reference/trim feasibility excludes retained coding stops, substitutions preserve the coding invariant, and structural tuples receive a complete post-state check. |
| Alignment order and orientation | P bases appear at their actual positions. Inverted D uses reverse-complement reference traversal consistently for alignment strings, CIGAR operations, identity and bounds. Reverse-strand observations preserve coding-oriented biological annotations. |
| Generating references | Standard V/D/J calls, alignment bases, bounds, and deletion fills use the final generating allele. Post-mutation inference is retained separately in `recalled_*_call`; its allele length cannot alter the truth alignment. |
| Junction annotations | Combined NP strings and lengths retain generating boundaries, including P bases; matches to adjacent trimmed germline do not absorb sampled junction bases. Physical junction AA is available for nonproductive records too. |
| Persistent provenance | Final assignment trims and P-addition totals survive clonal forks. V replacement exports its resulting trims and selects a different eligible V allele. Complete read loss retains the realized removal count. |
| Live-call cache | No-op mutation passes leave no pending calls, and staging/refresh share boundary overlap and D-orientation rules. Inferred calls in genotype experiments use only donor-carried alleles. |
| Validation and custom metadata | Validation compares the actual returned fields, types and presence against the outcome projection and checks record/outcome counts. Custom metadata cannot overwrite generated annotations. Native invariant checks remain in addition to that comparison. |
| AIRR export | Default TSV export converts starts to one-based closed coordinates, writes `T`/`F` booleans, and adds leading query/reference CIGAR placement offsets without changing the Python records. |

The relevant regressions include
[native biological alignments](../engine_rs/src/airr_record/tests/biological_alignment.rs),
[Python metadata repairs](../tests/test_metadata_biology_fixes.py),
[germline provenance](../tests/test_metadata_germline_regressions.py), and
[indel distribution checks](../tests/components/engine/test_distribution_invariants.py).
The indel tests enumerate complete productive support independently; they
separate retry-exhaustion NoOps from deletions and retain 5σ comparisons for
nondegenerate distributions. In the one-event fixture, valid edits are now
insertion-only, which is checked exactly.

Compatibility and remaining limitations:

- **Original-tip processing is opt-in.** The new model uses original coding-tip
  templates and coordinated net outcomes. Existing configurations retain their
  post-trim P semantics. Explicit duplex repair and empirical calibration remain
  modeling limitations; see the follow-up linked above.
- **End loss now changes observed coverage.** Terminal reference bases outside
  the read are not deletion columns. Actual structural gaps inside the
  observed interval remain `D`, including inverted D. Unobserved segments have
  no germline coordinates, even when earlier structural edits emptied them.
- **Standard call semantics changed deliberately.** `v_call`, `d_call`, and
  `j_call` name generating alleles. Inference consumers must use `recalled_*_call`.
  Optional `truth_*_call` fields remain compatibility aliases.
- **Mutation rate is event burden.** It sums each SHM stage's realized events
  divided by that stage's input molecule length; later read loss cannot inflate
  it. It is not an observed mismatch fraction. Lineage retains its net-difference
  summary using the pre-observation node length.
- **Internal CIGARs remain segment-local.** Python dictionaries and default
  DataFrames retain zero-based half-open coordinates and core CIGAR runs.
  Raw native CIGARs therefore still reproduce the baseline finding about
  missing AIRR placement offsets.
  Strict export adds placement before converting coordinates.
  `to_tsv(airr_strict=False)` preserves the previous TSV representation.
  Strict export performs these conversions; it is not a general validator of
  arbitrary custom metadata against every AIRR content restriction.
- **Reverse orientation is an intentional compatibility change.**
  V/D/J and junction coordinates index the coding orientation, while paired
  read windows continue to index the emitted observation. Inverted D retains
  forward-reference bounds; consumers of its reverse-complement alignment
  must honor `d_inverted`.
- **Productive P support uses local completion conditioning.** Fixed and
  mixed nonzero P lengths participate in a completion check over deterministic
  P bytes, allele/trim choices, N support, and D orientation. Each draw
  renormalizes its original local weights over admissible candidates; this
  does not define a globally conditioned joint probability law. Impossible
  support errors even in permissive mode instead of dropping requested
  additions. The P-after-trim model remains unchanged.
- **NP estimators need N-only evidence when P is present.** Combined NP
  annotations cannot safely be learned as TdT-only lengths or bases and then
  combined with a separate P model. The repair uses an explicit unsupported
  input guard for records carrying nonzero P counts; recovering N-only
  training evidence remains future work.
- **Seeded results and validation expectations can change.** Productive
  support, revision choices and repaired annotations differ from the
  baseline. A passing validator establishes agreement with the implemented
  model, not universal biological realism or complete receptor function.

The review checked the alignment/export interactions: both alignment
rendering and segment bookkeeping account for intervening P columns; D
orientation is shared by both; strict CIGAR prefixes use the original
zero-based starts; and reverse observations leave the coding annotations
consistent with their alignments.

## Follow-up: NP estimation with P-bearing data

Follow-up testing after correcting `np1`/`np2` projection exposed an estimator compatibility problem. Human IGK, seed 0, with N length three and P lengths V3=2/J5=1 produces `np1="tgCATt"` and `np1_length=6`: three N bases (`CAT`) and three P bases. The NP length estimator learned six N bases, while the separate P estimator retained lengths two and one. Regenerating with both models produced a nine-base combined interval. The NP base estimator also included the P bases, changing the N-only A/C/T weights from one third each to A/C/G weights of one sixth and T weight of one half. Neither estimator originally warned.

The cause is that [np_lengths.py](../src/GenAIRR/_cartridge_estimators/np_lengths.py) and [np_base_model.py](../src/GenAIRR/_cartridge_estimators/np_base_model.py) consumed combined AIRR measurements directly as models for the engine's separate N-addition sampler, while ignoring explicit P provenance.

The repair now rejects explicitly nonzero P-length provenance in either NP estimator before changing models or reports. Missing and zero P fields, including TSV string zeros, retain their previous behavior. This prevents the demonstrated double counting; it does **not** implement inference of N-only distributions from P-bearing data. Simply subtracting historical P counts or clipping that many characters is unsafe after read loss, revision, or alignment reassignment. Inputs without P provenance can still contain unrecognized P sequence. A future estimator needs an explicit, trustworthy N-only representation or a joint model.

## Follow-up: productive sampling with nonzero P additions

The temporary blanket rejection of nonzero P support has been replaced with
joint completion checks. Fixed and mixed P lengths must remain within the
authored support, and every deterministic P base must participate in frame and
stop checks, including codons spanning N/P boundaries. Future D inversion and
first-order Markov N-base support are included in those checks.
Phased genotypes restrict future allele support to the selected chromosome,
including positive gene usage and allele-copy weights. This prevents an
earlier choice from depending on a completing J allele on another chromosome.

The Python recombination pipeline uses the supported coding-order schedule.
Custom Rust plans require enumerable length supports, canonical first-order
N-base support, one N-generation pass per gap, and one P-length pass per end
in that order. Unsupported custom declarations produce a compilation error.

P-addition metadata still records emitted bases through persistent assignments
and `PRegionAdded` events. Valid trace replay must consume the recorded lengths
without substituting zero or drawing replacement randomness. Nonzero P remains
derived from retained oriented coding flanks; this work does not change the
biochemical model discussed above.

Focused genotype verification passes: nine native sampler tests and 17 Python
tests covering phased P completion, exact trace replay, and existing genotype
behavior.

## Verification of productive P support

The extension was rebuilt after the P-aware solver, runtime sampling, phased
genotype restrictions, and pool-relative codon phases were implemented.

- Full Rust suites: **1,428 passed in debug**, **1,427 in release**; the same
  two ignored integration tests and two ignored doctests remain per profile.
- Positive P Python suites: **44 passed**, including both strict/permissive
  policies, fixed/mixed lengths, boundary stops, Markov lookahead, trimming,
  inversion, replay, and clonal projection. Seven additional phased-genotype
  cases pass, including true trace replay.
- Six native exhaustive oracles independently concatenate complete molecules
  and compare solver decisions for allele/trim/P/N/orientation choices and
  N prefixes. All pass. Custom contracts remain active during sampling.
- Full Python regression run: **2,534 passed**, eight performance tests
  deselected, and one legacy documentation assertion failed. That assertion
  required replay documentation to say admissibility was never re-evaluated.
  It was corrected to describe proposal/support/completion validation, then
  all **14 tests** in `test_productive_failure_modes.py` passed. Thus every
  one of the 2,535 selected tests has a passing result against the final
  implementation, using the full run plus that focused rerun.
- [Additional chain verification](productive_p_support_verification_2026_10_20.json):
  **550/550 records passed** across 11 IGL/TCRA/TCRB/TCRG/TCRD configurations,
  including 102 D inversions and an IGL SHM configuration. Fixed P lengths
  remained exactly one at every applicable end; 519 records had positive P.
  TCRA/B/G used explicit curatable-reference opt-in; the default reference
  gate failures were independently reproduced without P as well.
- The original 19 diagnostic probes and 1,200-record independent coordinate
  sweep were rerun with the rebuilt extension. Their results remain as
  recorded below: zero coordinate failures and the same four explicitly
  retained model/format discrepancies.
- `git diff --check` passes.

The following historical totals describe the preceding repair pass, before
the productive-P follow-up.

## Verification before the productive-P follow-up

The Rust extension was rebuilt from this branch before verification.

- Full Python correctness suite: **2,476 passed**, eight performance tests
  deselected (`pytest tests -q -m 'not performance' --disable-warnings`).
- NP estimator suites rerun after the additional Python guard: **112 passed**,
  including 24 new guard regressions.
- Full Rust suite: **1,414 passed in debug**, **1,413 passed in release**.
  Each profile retains two ignored integration tests and two ignored doctests;
  the count difference is the debug-only panic regression.
- [Independent coordinate sweep](biological_metadata_fix_coordinate_sweep_2026_10_20.json):
  **1,200 records, zero failures**, across 12 IGH/IGK pipeline profiles.
- [Original diagnostic probes rerun](biological_metadata_fix_verification_2026_10_20.json):
  19 probes, with four discrepancies retained. These are the two partial-loss
  alignment conventions, the P-after-trim model, and native segment-local
  CIGAR placement. Strict TSV placement offsets pass the added check.
- The original fixed-seed batches now contain **zero falsely productive
  records with coding stops in each 1,000-record batch**, compared with 23
  clean and 292 SHM records at baseline. This is evidence under those exact
  settings, not a guarantee for every model and feature combination.
- `git diff --check` passes.

These results support the repaired behaviors above. The unresolved model
decisions and supported-format boundaries remain part of the review.
