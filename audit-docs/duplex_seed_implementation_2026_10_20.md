# Persistent duplex seeds and downstream observations

Implementation checkpoint for slice B of the
[projection and integration design](duplex_airr_projection_and_integration_design_2026_10_20.md),
on `fix/biological-metadata-correctness`. Slice A remains recorded separately
at commit `a805688`; its historical verification artifacts are unchanged.

## Implemented contract

`CoupledDuplexRepair.seed(...)` installs an exactly bound, completed repair
product into the native persistent simulation IR. It checks reference identity,
replayed repair evidence, unique source ownership, V-D-J order, native coordinate
limits and the full query partition before returning a seed. Valid overlaps
still have standalone partial projection; they cannot enter this single-owner
execution representation. Conversion never resamples a history.

`project_duplex_simulation(simulation, refdata)` produces the nullable AIRR
projection for that simulation's current observation. The generating calls,
net repair trims and immutable repair snapshot retain their founder meanings.
Current query coordinates, alignments, certified body claims, productivity and
separately named re-calls describe the observation. Both join orders and inverted
D retain canonical reference coordinates and oriented alignment traversal.

The persistent sidecar has shared immutable founder evidence and independent
copy-on-write query links, edit events, removed links and sampled PCR/quality
stages. The current query and ancestry must replay exactly, including old/new
bytes and removal causes. A founder-column link represents chemical/information
ancestry after copying, rather than identical physical atoms in descendants.
Observation insertion references are explicitly branch local.

Native continuations accept only concrete audited implementations: uniform/S5F
mutation, PCR/quality substitution, N masking, structural indels, terminal read
loss and reverse-complement presentation. The whole plan is checked before
execution, including an unsupported suffix. Matching a trusted pass's name or
effect declaration cannot authorize a custom implementation. Bindings,
assignments, coverage and structural layout are checked on re-entry.

This is an experimental native continuation API. Normal `Experiment` finite
history sampling, family integration and a durable replay format bound to the
initial duplex snapshot remain future slices. Legacy AIRR and family exporters,
legacy record validation, durable trace-file export and preservation-only
`AnchorPreserved` fail explicitly on duplex states. The checked Rust legacy
builder returns an error; its existing infallible wrapper fails closed on
misuse. The standalone repair script itself remains deterministically replayable.
No external tool or dependency was added.

## Biological and metadata safeguards

- Body placement can be supported by either physical strand. Synthesized
  restoration retains its chemical birth and template ancestry.
- Structural deletions consume reference positions. Terminal clipping reduces
  coverage without changing founder repair trims or becoming a structural indel.
  A wholly deleted gene retains a deletion-only path when coverage supports it;
  its query interval is missing and a zero-width path boundary locates the CIGAR.
- Duplex productivity uses certified landmarks and actual current codons under
  `duplex_observed_vdj_v1`. The canonical legacy policy is unchanged. Structural
  phase shifts and read clipping are distinct; a coincidental matching motif
  cannot restore a removed landmark. Productive indel sampling uses authoritative
  materialized views rather than a cache that drops duplex evidence.
- NP amino-acid strings include only complete coding-frame codons contained in
  the NP interval. Missing phase produces a missing translation.
- Mutation rate retains each mutation stage's input molecule denominator.
  Cumulative per-mechanism counts survive successive continuations and cloning.
  PCR/quality scalar counts represent sampled attempts, including empty-query
  stages. Applied substitutions are separately named in the extension.
- Reverse-complement output changes `sequence` and its presentation flag;
  alignment and provenance coordinates stay coding oriented.
- Unknown legacy P-emission lengths stay missing. Raw edits without an attributed
  mechanism do not receive fabricated mutation or corruption counters.

## Boundary bugs found while implementing this slice

The observation audit identified additional errors before completion. A final
query-boundary insertion was excluded from the final region under ordinary
range-shift rules; duplex ownership now follows the explicitly chosen source
region at every insertion boundary, including an empty query. NP translation
was incorrectly restarted at each tract rather than using the coding frame.
A stale observer candidate could retain the old J endpoint after a terminal
insertion; duplex structural re-calls now refresh against the actual final IR.
The J-tail phase check now uses observed reference coverage, so clipping an
inserted suffix or clipping past an earlier structural deletion cannot confuse
query loss with reference loss. Very short clipped queries safely yield empty
complete-codon translations.
PCR/quality substitutions alone could not recover sampled attempt counts,
particularly after complete clipping; persistent stage summaries now preserve
them separately. Integrity checking now includes current byte replay and
independent observation-coverage/partition validation.

## Validation and remaining scope

Final build identity, regression results and literal sweeps are recorded in
[the verification artifact](duplex_seed_verification_2026_10_20.json).
The release extension SHA-256 is
`f29f1b3cf99847cd435436b997c27141bf2203b5bc12c78ff5ed76154fc901b3`.

- Full Python regression: 2,891 passed; eight performance cases excluded.
- Performance: all eight passed with unchanged thresholds.
- Native debug: 1,583 passed; native release: 1,582 passed. Both retain three
  inherited ignored tests; release excludes one debug-only test.
- The 97 affected Python tests pass independently.
- Projection sweep: 32,768 controls, 17,408 accepted and 15,360 correctly
  rejected as incompatible, with no errors.
- Observation sweep: 2,048 founders and 18,432 native continuations, with no
  errors.
- Strict documentation build passed. All four guide examples executed in order
  in a temporary directory and produced their expected exports.
- Staged diff check and source-manifest verification passed.

The permanent observation sweep replays emitted edits independently and verifies founder-column
conservation, current body remapping, alignment/query conservation and immutable
repair counts across both join orders and D orientations.

These are authored mechanistic and synthetic controls. They do not fit measured
enzyme rates, establish repertoire frequencies or validate unrestricted repair
kinetics. General flap/mismatch processing, reopening completed joins and measured
biochemical calibration remain explicit model extensions. The next integration
slice is a finite, reference-bound authored history source in `Experiment`, with
support-wide capability validation and no rejection resampling.
