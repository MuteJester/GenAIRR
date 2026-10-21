# Productivity, indels and coding-end repair implementation

This checkpoint was subsequently committed as `f066b44`. The
[coupled D repair follow-up](coupled_duplex_implementation_2026_10_20.md) extends
the experimental local kernel; the verification below describes the preceding
checkpoint.

Implemented on `fix/biological-metadata-correctness`, following the
[approved plan](productivity_indel_and_repair_plan_2026_10_20.md). The baseline
commit is `1ededa676e8904430f7ea0bdcf03092de4888a41`; the changes remain in the
working tree. All simulation, annotation resolution and estimation run inside
GenAIRR. No external simulator, inference tool or new runtime dependency was
added.

## Productivity and reference evidence

The shared `canonical_observed_vdj_v1` assessment evaluates original landmark
provenance, live V Cys/J Phe-or-Trp residues, coding frame, observed coverage and
possible coding stops. Its scope is `observed_coding_span`. Unknown reference
landmarks, missing required coverage or unresolved stop possibilities can yield
`None`; a reference-admission override cannot turn that uncertainty into
productivity. IUPAC codons are evaluated as sets of possible translations:
`GCN` is definitely sense, `TAR` definitely stop and `NNN` unresolved.

Ordinary and original-end productive recombination use complete continuation
support. In particular, zero-P VDJ generation now checks future D and junction
stops before selecting support. Mutation, indels, lineage divisions, final
emission and independent record validation use the same policy semantics.
Both permissive and strict productive runs check the final assessment. The
record includes policy, scope and reason codes; a known junction keeps its
actual translation when mutation changes an anchor residue.

The final guard exposed an additional receptor-revision defect: permissive V
replacement sampled length-compatible alleles without checking biological
admissibility. Both ordinary and genotype-aware revision now filter hypothetical
replacements before drawing, preserve the eligible weights and emit events only
for the committed replacement. Replay enforces those contracts in both modes.
Known failing human IGH seeds are permanent regressions.

The complete solver now iterates reachable nucleotide-automaton states and uses
direct IUPAC stop certainty instead of repeatedly enumerating translations.
An unchanged 100-seed full-record comparison checks that these optimizations do
not change generated choices or metadata.

Native reference authoring now resolves supported IMGT V landmarks and unique
J motifs explicitly. Optional sequence-bound evidence carries landmark status,
independent coding phase and coding-boundary completeness. Source functionality
and unknown species are preserved. A new portable authoring-checksum scheme
survives serialization round trips without changing historical archive hashes.

The reproducible archive audit covers **106** bundled cartridges: **80** have
strict validation findings and are marked quarantined; the other **26** pass
structural validation but still lack independent source/phase/boundary
certification. No archived allele was renamed, removed or corrected by a motif
guess. Quarantine is an explicit unresolved result, not a repaired reference
release. Source-backed repair still requires pinned authoritative source data.

## Natural-law indels and replay

The versioned `natural_sequential_v1` sampler conditions actual sequential edit
histories on the final contracts while preserving the already drawn count.
Recoverable temporary frameshifts, stops and inserted interruptions are allowed
inside a batch. Counts zero through two have exact finite fallback support;
short-tuple evaluation uses indexed edit views and is checked against fully
materialized histories and independent probability oracles.

Higher counts use bounded complete-history proposals and weighted exact search.
The higher-count fallback interleaves independent full-history probes with
deterministically scheduled search, charging both to the existing work budget.
A successful independent probe has the same conditional history law; an
unfinished weighted search is never sampled. Setting `rejection_attempts=0`
retains forced exhaustive mode for oracle checks.
Unfinished search reports `SamplingBudgetExceeded`; only a complete enumeration
or sound proof can establish `EmptyAdmissibleSupport`. No partial reservoir is
returned as an exact sample. Diagnostics record requested/applied counts,
strategy and work independently of biological events. Arbitrary large counts
and custom constraints are not promised cheap or guaranteed completion.

Replay checks typed coordinates, the configured count, positive kind/base
support, evolving bounds and final constraints. Explicit legacy compatibility
is restricted to historical choice replay; a seed rerun cannot silently use an
old sampler law. Trace schema 4 binds reference identity, address schema,
productivity policy and active contracts. Seed reruns and choice replay apply
the same context gates.

## Retained-outcome calibration

`CodingEndOutcomeModel` carries finite signed V/J outcomes, joint D rows,
compatible reference identities and training provenance. Direct `d_joint`
support preserves D-end correlation in simulation. Supervised fitting reads
completed generating histories and requires explicit independent-event and
unconditioned-sampling declarations. Known lineage/clone/collapsed observations
are rejected. Counts, marginal Wilson intervals and held-out table scores are
reported. These are retained generating-outcome frequencies, not physical
resection or enzyme-rate estimates.

The finite compatible-history fitter sums caller-supplied complete explanations
with their forward likelihoods. Its report includes convergence and residual
history entropy, explicitly leaves identifiability unestablished, and supplies
no invented biochemical confidence intervals. Generic AIRR annotations are not
automatically converted into unique processing histories.

As requested, calibration validation uses synthetic recovery now. The larger
native experiment uses 2,000 independent training events and 500 held-out
events, one feasible allele per segment, correlated D outcomes and D inversion.
Measured biochemical fitting and donor/clone validation await measured data.

## Experimental duplex repair

`GenAIRR.experimental.DuplexRepair` exposes a bounded, explicitly authored local
two-ended repair kernel. Both physical strands retain separate nucleotide
identities and original coordinates. Opening transfers original bases across
the hairpin turn; resection retains a removed-base ledger; annealed product
columns retain both contributors. Template synthesis records its template and
requires a paired 3′ primer and continuous opposing strand. Ligation requires
two resolved nicks and a fully paired product. Exact scripts replay through the
same kernel.

This first experimental slice handles one local join. It does **not** yet feed
the standard AIRR pipeline, model successive coupled D joins or provide general
flap/mismatch processing, kinetic rates or biochemical calibration. Unsupported
substrates fail explicitly. The existing `original_coding_end_v1` AIRR model
remains the supported retained-outcome simulator; the experimental kernel does
not add a second set of N/P bases to it.

Verification includes independent short-substrate examples, paired-primer
rejections and an exhaustive dinucleotide opening matrix: 320 accepted products
and 192 incompatible overhang pairs, with strand symmetry and nucleotide
conservation checked for every accepted product.

## Verification results

The final release extension SHA-256 is
`03b21c9769157a78eb4797477be6acfb02265acdd85d58ef8b362950ac8d0e78`.
The [verification artifact](productivity_indel_repair_verification_2026_10_20.json)
records commands, logs, reference/build identities, artifact checksums and the
modified package source hashes.

| Check | Result |
| --- | --- |
| Full Python suite excluding performance | **2,765 passed**, 8 performance tests deselected |
| Existing performance suite | **8 passed**, budgets unchanged |
| Full native debug suite | **1,521 passed**, 0 failures |
| Full native release suite | **1,520 passed**, 0 failures |
| Strict MkDocs build | Passed |
| Cross-locus original-end sweep with strict productive two-indel batches | **660 records**, 22 completed profiles, 0 record errors; replay identical |
| Explicit reference/support guards in that sweep | **6 rejected profiles**, recorded separately |
| Fixed former 370/376-nt failing precursors | **80/80** continuation seeds completed |
| Indel count/length/probability matrix | **1,080/1,080** completed across 36 profiles |
| Native synthetic retained-outcome recovery | 2,000 training / 500 held-out events; maximum probability error **0.011** |

The debug/release difference is an existing debug-only test. Each ordinary
native suite leaves three opt-in/example tests ignored: two documentation
examples and the exhaustive count-two timing probe. That timing probe was also
run separately and passed at all three receptor lengths. The previously ignored
mutation/deletion metadata tests were repaired to exercise real native passes;
all **86** native metadata tests now run and pass. They check stable generating
calls/references separately from re-calls and preserve event-driven refresh.

Current build workload timings are **0.985 s** for 1,000 productive full-stack
records (W3, unchanged 2 s budget) and **0.792 s** for 1,000 productive records
with zero through five indels (W6, unchanged 1 s budget). Timings depend on the
machine and ambient load. The forced exact two-indel fallback visits 3.43,
6.31 and 10.63 million states at 369, 501 and 651 bases, taking 2.035, 2.896 and
4.237 seconds here. Those are fallback enumeration costs, not normal throughput.
The matrix reports cumulative process peak RSS explicitly rather than claiming
per-record native allocation measurements.

Reproducible evidence:

- [Cross-locus metadata and replay sweep](productivity_indel_cross_locus_2026_10_20.json)
- [Synthetic calibration recovery](coding_end_synthetic_recovery_2026_10_20.json)
- [Indel workload matrix](indel_workload_matrix_2026_10_20.json)
- [Former failing fixed precursors](indel_fixed_precursors_2026_10_20.json)
- [Forced exact fallback timings](indel_exact_fallback_2026_10_20.json)
- [Current productive workload timings](productivity_performance_2026_10_20.json)

These checks support the declared sequence/provenance policies and bounded
models. They do not certify every reference annotation, expression or receptor
function, identify unobserved microscopic histories, or establish measured
biochemical rates. Coupled D duplex joins, general flap/mismatch repair and
AIRR integration of the experimental duplex ancestry remain further model
milestones. The earlier limitation report and proposed plan are historical
evidence rather than the current correctness status.
