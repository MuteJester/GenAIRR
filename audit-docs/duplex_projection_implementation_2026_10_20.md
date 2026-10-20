# Standalone duplex AIRR projection

Implemented on `fix/biological-metadata-correctness` after design checkpoint
`c480374`. This completes the standalone projection portion of slice A in the
[integration design](duplex_airr_projection_and_integration_design_2026_10_20.md).
Persistent simulation installation and normal `Experiment` integration remain
subsequent slices.

## Resulting behavior

`CoupledDuplexRepair.project(reference, v_call=..., d_call=..., j_call=...)`
returns a bound AIRR record, independent gene paths, availability reasons and
the complete paired repair ancestry. `export_duplex_projection` writes the
full JSON envelope or a strict AIRR TSV row with the ancestry as a JSON cell.
The [guide](../site_docs/guides/experimental-duplex-repair.md) includes a runnable
example.

Bindings require exact generating allele names and whole canonical allele
substrates. They retain reference content identity, sequence checksums,
source-instance identity, D orientation and the reference evidence used for
landmark/phase readiness. Unknown junction-facing boundaries require an
explicit coding-end assertion; known partial boundaries cannot be overridden.
No reference annotations are silently curated.

The existing native transition executor now accepts an optional evidence
observer. Opening, resection, synthesis, component translation and ligation
emit deterministic operation evidence. The old local and coupled `run()`
products retain their existing schema and physical behavior. There is no
second implementation of the repair chemistry.

Physical ancestry, template ancestry and body placement remain separate.
Resecting a D-top base and restoring it from D-bottom records the physical
removal and synthesis while retaining the certified D alignment position.
Transferring that synthesis around a later hairpin revokes its previous body
placement without losing its ID or template ancestry. Historical template
links remain resolvable after template removal.

The query is coding-oriented top. A retained body certificate on either
paired strand can support a column's reference alignment. Agreeing claims are
deduplicated; genuine V-body/D-body overlap retains both independent gene
paths. Unsupported combined alignments and their relative coordinates remain
missing, with explicit reasons. The adapter never chooses one overlapping
gene or uses a best sequence match to invent an owner.

Canonical D coordinates remain forward-reference coordinates. Its alignment
traverses the oriented reference, and strict export uses that oriented CIGAR
offset. Operation tips retain both oriented and canonical ends under inversion.
Stable record IDs distinguish different bindings/histories even when their
final sequences agree.

## Productivity and measurement contracts

`duplex_observed_vdj_v1` accepts a synthesized landmark only with certified
in-place placement. It reuses existing canonical residue, phase, stop and
coverage semantics without changing `canonical_observed_vdj_v1`.
Verified positional evidence can locate a junction without coding phase;
a known rejected residue remains a policy failure. A defined junction keeps
its actual translation on nonproductive products.

One regression produces the same `AAATGTCCCTTTAAA` query by restoring the V
landmark, or by removing that landmark and adding untemplated `TGT`. The
first can be productive; the second retains its N ancestry, lost V landmark
and determinate rejection. Identical sequence does not authorize identical
causal metadata.

Net trims describe missing terminal body positions after repair, not physical
resection totals. NP strings describe query intervals between body paths;
two physical N additions paired in one column contribute one query base.
Repair events do not increment SHM, PCR or downstream indel counters.

Existing P scalar fields remain missing because the coupled kernel has no
equivalent legacy P-emission contract. Separate counters describe physical
operations and top-query placement roles after both joins. Missing values
remain distinct from known empty tracts through JSON and TSV export.

The legacy P, NP, trim and retained-outcome fitters now explicitly reject
duplex provenance. Previously an NP fitter could accept the missing P scalars
and train an N-only model on a combined tract. Rejection occurs before model
or report mutation, for both dictionary and TSV inputs. Learning a distribution
over coupled histories requires a separate fitting contract.

## Verification

Final extension SHA-256:
`3b601a4d6df473e8c45bc3b31da1be8148a67899d594f248e823e49e40f9950c`.
Source hashes, commands and logs are in the
[verification artifact](duplex_projection_verification_2026_10_20.json).

| Check | Result and scope |
| --- | --- |
| Final affected Python suites | **334 passed**, including 29 new projection/export/fitting controls; the short blunt test also enumerates 256 inner controls |
| Final native debug suite | **1,541 passed**, 0 failed, 3 inherited ignored tests/examples |
| Final native release suite | **1,540 passed**, 0 failed, 3 inherited ignored tests/examples |
| Final exhaustive projection sweep | **32,768 controls**, 0 mismatches: 17,408 accepted, 15,360 expected incompatible rejections |
| Existing performance suite | **8 passed**, budgets unchanged |
| Strict documentation build and formatting of changed Rust files | Passed |
| Initial full Python regression | **2,813 passed**, 8 performance tests deselected; run before final record-ID, canonical-tip and fitter-guard additions |

After the initial full Python run, the final release build was checked with
the affected suites, both full native suites, exhaustive sweep and performance
budgets above. The full Python run was not repeated after those final
additions. The first performance invocation exceeded W6's one-second limit
at 1.335 seconds; the isolated W6 check and subsequent complete eight-test
suite passed without changing code or thresholds.

The sweep independently predicts complementarity, literal sequence,
canonical D paths, empty NP intervals and zero net trims across all
dinucleotide V/D/J triples, both orientations, both orders and symmetric
opening offsets zero/one. It checks placement and physical conservation;
it does not measure biochemical frequencies. Separate literal fixtures cover
restoration, true body overlap, mixed physical/body ancestry, removed-template
chains, nonproductive junctions and asymmetric inverted-D export. Corrupt
native snapshots/certificates fail deterministic replay rather than becoming
partial biological annotations.

## Remaining integration work

This adapter is experimental and standalone. `capabilities.simulation_seed`
reports compatibility with the planned single-segment layout; it does not
install a simulation. `experiment_integration` remains false. The current
normal `AirrRecord` and assignment scalar migrations, persistent observation
maps, mutation/indel/clipping updates, finite weighted history sampling and
clonal continuation remain to implement and validate before normal
`Experiment` support is exposed.

Whole-allele resolved A/C/G/T substrates remain required. This slice adds no
VJ adapter, receptor revision, general overlap execution IR, flap/mismatch
repair, reopening of protected ligations or measured biochemical rates. All
implemented repair, projection, export and validation remain inside GenAIRR,
without an external simulator or inference dependency.
