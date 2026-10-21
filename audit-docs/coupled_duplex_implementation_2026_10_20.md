# Coupled D duplex repair

The completed metadata/productivity/indel work was checkpointed as
`f066b44` on `fix/biological-metadata-correctness` before this slice began.
This slice implements successive V-D and D-J repair on one shared D substrate
inside GenAIRR. It adds no external simulator, inference service or runtime
dependency.

## Resulting behavior

`GenAIRR.experimental.CoupledDuplexRepair` accepts canonical V, D and J
sequences, an explicit D orientation and two authored `DuplexJoin` scripts.
Exactly one `vd` and one `dj` join must complete; either order is supported.
The existing `DuplexRepair` API and result schema remain available.

The implementation allocates all original physical nucleotides once in V/D/J
order. D inversion reverses placement and exchanges physical strands without
changing an original nucleotide's identity, chemical base, canonical strand or
canonical reference coordinate. Joining consumes existing molecular components.
The first product is never reinitialized as a new germline sequence.

Both APIs use one native transition executor for opening, resection, N addition,
synapsis, template-dependent synthesis and ligation. Opening moves actual
surviving material, including synthesis from a preceding join. Synthesis requires
a paired 3′ primer and a continuous opposing template. Each local join must
resolve its duplex; final output requires one connected, fully paired VDJ
product.

`right_start` uses stable original oriented origins: D relative to V for `vd`,
J relative to oriented D for `dj`. Resection never redefines an origin. Component
translation moves both physical strands, all source-origin anchors and previous
ligation boundaries together. The result exposes final `origin_positions` and
both physical strand boundaries in `ligation_seams`.

## Coupling and ancestry

The original D length alone cannot determine support for the second opening.
After the first join, pairing supplied by V/J or earlier synthesis can support
the remaining D hairpin. Permanent regressions include:

- V `AAG`, D `CA`, J `TAA`, with opening offset one at all four tips:
  `AAGCATAA`, both join orders and D orientations. After both openings no
  original-original D pair remains, while the physical product stays paired.
- A second opening of offset two on the original two-base D, supported by the
  actual larger component and bounded by the previous ligation bond.
- V `AAG`, D `C`, J `GAA`: a first-join N addition and its synthesized mate
  support the remaining single-base-D opening, producing `AAGCCGAA`.
- Synthesis from the first join moving onto the opposite physical strand during
  the second opening, retaining its identity and template link.
- A template being removed after synthesis, with a later synthesized base
  copying that first-generation synthesis. Both historical template edges stay
  resolvable.

The global `ledger` retains every original and synthesized nucleotide once.
Each allocated ID occurs exactly once among final columns and removed bases.
The two contributors to each paired column remain separate. An original record
contains `segment`, canonical `strand` and `reference_position`; synthesized
bases retain their `kind` and optional `template_id`. Sequence matches never
turn synthesized material into germline provenance.

Budgets cover both joins together. The base budget includes both initial
strands of V/D/J and every new nucleotide. Resection gives no refund. Invalid,
incomplete or over-budget scripts return errors without a completed product.

## Verification

Final extension SHA-256:
`a59585ffc7e20db37e858d3ec1be13484b12d4367453afcd04cc72421187a302`.
Commands, source/build hashes and scope are in the
[verification artifact](coupled_duplex_verification_2026_10_20.json).

| Check on the rebuilt extension/source | Result |
| --- | --- |
| Full native debug suite | **1,530 passed**, 0 failures |
| Full native release suite | **1,529 passed**, 0 failures |
| Focused Python coupled/local/calibration suites | **67 passed**: 29 coupled, 19 existing local, 19 calibration |
| Existing performance suite | **8 passed**, budgets unchanged |
| Strict documentation build | Passed |
| Permanent exhaustive coupled sweep | **32,768 cases**, 0 mismatches |
| Independent implementation of the exhaustive audit | Same expected accept/reject results, no defects found |

The exhaustive sweep covers all 4,096 dinucleotide V/D/J triples, both D
orientations, both join orders and symmetric opening offsets zero/one. Literal
complementarity predicts **17,408 accepted products** and **15,360 incompatible
overhang rejections**. Checks include exact product sequence, canonical base
identity, conservation, stable IDs under inversion, a single D allocation and
commutation for these independent compatible fixtures. These are the same
32,768 cases checked by two implementations, not 65,536 distinct histories.
The [sweep artifact](coupled_duplex_sweep_2026_10_20.json) records its exact build.

The Python suite also checks 1,280 short blunt histories with hand-derived
source IDs, exact JSON replay, cumulative limits, removed-template ancestry and
rejection at prior ligation bonds. Existing local-duplex expectations still pass
after the shared-kernel refactor. The native debug/release difference and three
ignored opt-in/example tests are inherited from the checkpoint. The full
2,765-test Python run belongs to that preceding checkpoint; this slice reruns
the affected Python suites, full native suites and existing performance budgets.

## Explicit model boundaries

The input contract supplies three coding duplexes and four junction-facing
hairpins. It does not claim RSS recognition, RAG activation or a detailed
cleavage chronology. Both join orders are authored mechanistic controls under
that contract, without a fitted biological prevalence.

Previously ligated bonds and their adjacent physical nucleotides are protected.
Opening or erosion that would consume/cross those bonds is explicitly outside
this bounded model. Such errors do not assert universal biological impossibility
under other repair mechanisms. General flap/mismatch repair, reopening a
completed junction and biochemical kinetics remain extensions.

The coupled model is still experimental and separate from standard AIRR
simulation. Projecting its many-origin paired ancestry into AIRR fields is a
further milestone. Its authored choices are not measured enzyme rates; measured
calibration awaits the later data discussed with the user. Existing retained-
outcome fitting and synthetic recovery remain available inside GenAIRR.

The [experimental guide](../site_docs/guides/experimental-duplex-repair.md)
documents the API, substrate contract and primary biochemical sources. The
short-substrate examples above are independent GenAIRR geometry controls, not
experimental measurements from those studies.
