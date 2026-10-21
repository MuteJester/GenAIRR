# Duplex family continuation and unfitted defaults

This is slice D of the
[approved integration design](duplex_airr_projection_and_integration_design_2026_10_20.md),
following finite-source checkpoint `a75db8c` on
`fix/biological-metadata-correctness`. The user authorized implementation and
explicit logical defaults while measured repair data are unavailable. All
mechanisms remain inside GenAIRR.

## Default probability contract

Omitting `DuplexHistory.weight` means `1.0`. For a realized source context,
equal positive weights select uniformly among its eligible authored history
IDs. Distinct histories with identical DNA retain separate mass. The existing
local conditioning law and complete physical/capability preflight still apply;
uniform weights never admit an invalid operation or unsupported layout.

This is an unfitted distribution over a declared finite support. It supplies
neither a uniform distribution over all possible chemical histories nor a
measured enzyme rate. Later fitted row weights use the same immutable model
API. Existing cartridge usage is retained when available, with the existing
uniform allele fallback when absent. Clone sizes, division rates, SHM and
selection retain their existing explicit GenAIRR parameters; neutral selection
remains the default. No outside simulator or new data dependency is introduced.

## Founder and branch truth

Each clone samples and repairs one founder. Descendants inherit its generating
alleles, history selection, full paired repair ledger, template ancestry,
certified founder placements and physical operation counts. Immutable repair
evidence is shared. Observation maps, removed links and event histories diverge
independently after the fork. Repair IDs describe founder ancestry after
replication, rather than the same physical atoms in different cells.

The new mode supports fixed clonal expansion, clonal repertoires and BCR lineage
trees. Existing pre/post-fork ordering guards remain. Phased genotype plus
families, VJ repair, receptor revision, contaminants and paired-end adapters
retain their existing explicit restrictions.

Every node outcome must retain actual division events and cumulative mechanism
counters. Reconstructing a synthetic net difference from germline loses
back-mutation history and is unsuitable for duplex nodes. The legacy node
projection remains unchanged; the duplex path composes founder and actual
incoming-edge outcomes along the ancestry path. Downstream errors, clipping
and structural edits preserve their existing distinct semantics.

Productive constraints apply to the founder, each division and the final
observed read under `duplex_observed_vdj_v1`. Custom/preservation-only contracts
remain gated. Direct native family entry validates the bound reference and
repair/observation evidence. The direct simulation-only entry point requires
explicit reference data for duplex founders.

Affinity scoring translates a contiguous, phase-certified observed coding
span. Per-region translation can drop codons crossing V/NP/D/J boundaries and
cannot define duplex affinity. Neutral growth without a target needs no
affinity translation. Requested scoring with unresolved coding evidence fails
explicitly rather than assigning a fabricated zero affinity.

## Observation and abundance

Identical DNA does not establish identical ancestry or counters. Duplex tree
sampling therefore aggregates repeated draws of the same cell ID; different
sampled cell IDs remain separate observations even when their DNA matches.
Their abundance counts sample occurrences of that cell. Legacy genotype
collapse remains unchanged.

Duplex repertoire collapse keys include all generated record truth except the
presentation sequence ID. Records with equal DNA but different repair/observed
histories remain separate. User metadata is applied after truth-based collapse.
Fixed clonal expansion continues to emit each descendant separately.

All three result paths preserve nullable AIRR output and canonical JSON
ancestry cells in delimited exports. Family validation compares immutable
founder evidence and detects dropped sidecars, including against retained
parent outcomes. Current paths and observation changes may differ. Parent
indices are distinct from clone IDs; extinct accepted lineages retain their
founder/tree positions. Branch-local insertion identities require the enclosing
clone/node/read context when compared across records.

## Exact replay boundary

`FamilyOutcome.replay_node(node_id, refdata)` reconstructs any retained node,
including the root and unobserved cells, from its founder and actual recorded
division choices. It checks reference
identity, contracts and the stored mutator configuration, makes no replacement
draws, and compares sequence, counters and full ancestry with the retained
node. It returns the cellular state before downstream read observation; native
family objects remain accessible from the lineage result.

This is in-memory branch replay. A family is not one ordinary compiled pass
plan. Generic durable plan traces continue to reject synthesized/merged family
outcomes. A durable full-family package and durable arbitrary initial-snapshot
format remain separate work; no persisted-family replay is implied here.

## Verification gates

- Equal default and explicit-one weights produce identical models, source
  choices, records and source replay; upstream contexts retain their natural
  law despite unequal history counts.
- Shared founder evidence survives flat/tree branches, restoration, either D
  orientation and join order, mutation, structural edits, clipping and reversal.
- Independent literal mutation histories include repeated and back mutations;
  node counters and replay retain every actual event.
- Equal-DNA observations with different ancestry are not silently collapsed;
  sampled cell abundances conserve the draw count.
- Nullable JSON/DataFrame/TSV output, actual-record and family/parent validators,
  wrong-reference replay, unsupported contracts and mutated evidence fail as
  specified.
- Existing family outputs, seeded behavior and trace signatures remain intact;
  run affected and full native/Python suites, unchanged performance budgets and
  strict documentation before committing the checkpoint.
