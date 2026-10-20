# Experimental duplex repair

`GenAIRR.experimental.DuplexRepair` is an explicitly authored local experiment
for two complete coding-end hairpins. The Rust kernel represents both strands
and preserves their separate origins. `CoupledDuplexRepair` extends this to two
successive joins sharing one D molecule. These experimental APIs do not yet feed
the standard AIRR experiment pipeline.

```python
from GenAIRR.experimental import DuplexRepair

script = (
    DuplexRepair("AAG", "CAA")
    .open("left", offset=1)
    .open("right", offset=1)
    .synapse(right_start=3)
    .ligate()
)
product = script.run()
assert product["sequence"] == "AAGCAA"
replayed = DuplexRepair.from_dict(script.to_dict()).run()
assert replayed == product
```

Top strands follow the coding coordinate axis from 5′ to 3′. Bottom strands
follow the same axis from 3′ to 5′. `right_start` positions the original right
coding origin on the left coding axis. Opening at offset `h` transfers original
bases across the hairpin turn and produces a `2*h` 3′ overhang; opening itself
does not synthesize nucleotides. The paired annealed columns can retain origins
from both input ends.

## Supported transitions

- `open(side, offset=...)` opens each original hairpin once while retaining a
  paired stem.
- `resect(side, strand, count=...)` removes bases from the exposed end of the
  selected strand and preserves the removed origins.
- `add_n(side, sequence)` extends the exposed 3′ terminus. Its sequence follows
  polymerization order, 5′ to 3′, for either side.
- `synapse(right_start=...)` aligns the ends. Opposing bases must complement;
  overlapping same-direction strands requiring flap processing are rejected.
- `fill(strand)` extends a paired 3′ primer across a gap using a continuous
  opposing template. Each new base records its template identifier.
- `ligate()` requires two resolved 3′-OH/5′-phosphate nicks and a fully paired
  product.

Opening, resection and N extension precede synapsis in this bounded model.
Ligation completes the history. Unsupported transitions, ambiguous input,
unresolved mismatches, gaps and work-budget violations raise explicit errors.
An isolated 3′ overhang cannot be filled from a recessed 5′ end. Discontinuous
or unpaired-primer polymerase mechanisms require additional explicit operations.

Physical resection and final net loss are distinct. For example, resecting the
last top base of a blunt left end and filling it from the retained bottom strand
can restore the same final coding sequence while recording the resection and
new synthesis separately.

The result contains both original references, both final strand sequences,
paired columns, nucleotide identities, original strand coordinates, template
links, removed bases and the operation history. Only completed duplexes produce
a sequence. `max_operations` and `max_bases` bound the experiment.

## Model scope

The operations are specified mechanisms with authored choices. They carry no
fitted enzyme rates, universal opening probabilities or inferred unique repair
history. `DuplexRepair` handles a single local join, compatible annealing
and template completion; the coupled API below shares material across two joins.
Flap processing, mismatch resolution, general repair
kinetics and experimental calibration remain explicit extensions.

Verification covers nucleotide conservation during opening, paired ancestry,
opposing N additions sharing a product column, resection/restoration, strand
symmetry, exact script replay and rejection of invalid primer/template geometry.

## Coupled V-D and D-J repair

`CoupledDuplexRepair` initializes the two strands of V, D and J once. Each
`DuplexJoin` acts on the actual connected components left by preceding joins.
It preserves the original nucleotide identities and any earlier synthesis.
The following inverted-D example performs D-J first:

```python
from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin

def blunt(name, right_start):
    return (
        DuplexJoin(name)
        .open("left", offset=0)
        .open("right", offset=0)
        .synapse(right_start=right_start)
        .ligate()
    )

script = (
    CoupledDuplexRepair("ACG", "AAT", "TCC", d_inverted=True)
    .join(blunt("dj", 3))
    .join(blunt("vd", 3))
)
product = script.run()
assert product["sequence"] == "ACGATTTCC"
assert CoupledDuplexRepair.from_dict(script.to_dict()).run() == product
```

Exactly one `vd` and one `dj` join must complete. Either order is supported.
The order is an authored part of the history, without a fitted probability.
Each join supports the same opening, resection, N extension, synapsis, fill and
ligation operations as a local script. Joining consumes existing components;
it never constructs a new germline substrate from an intermediate sequence.

### Ends, coordinates and inversion

For `vd`, `left` means the exposed V tip and `right` means the exposed left tip
of oriented D. For `dj`, `left` means the exposed right tip of oriented D and
`right` means the exposed J tip. D's canonical ends exchange assembled roles
when inverted. `top`/`bottom` always describe the current physical coding axis;
an original base separately retains its canonical source strand.

For `vd`, `right_start` locates D's original oriented origin relative to V's
origin. For `dj`, it locates J's original origin relative to oriented D's origin.
These origin anchors survive resection and move with the whole component at
synapsis. They are not recomputed from the first surviving base. The result's
`origin_positions` gives their final common-axis positions, and
`ligation_seams` records both strand junction boundaries.

Canonical V/D/J references are supplied as input. Inverting D exchanges its
physical strand maps and reverses their placement. It preserves each original
base's identity, canonical strand, canonical reference coordinate and character.
It does not create an independent reverse-complement copy of D.

### Shared substrate and ancestry

Both D ends share finite physical material. The first join can alter the
substrate available to the second, including supplying paired material or
template-dependent synthesis. Opening transfers the bases actually present;
their identity and original/template ancestry survive movement between strands.

A sum of the two opening offsets is insufficient to decide support. For
example, V `AAG`, D `CA` and J `TAA`, with offset one at all four facing tips,
can produce `AAGCATAA` successively. The first partner supplies pairing retained
through the second opening; after both openings no original-original D pair
remains.
Transition checks use the current component geometry and the protected prior
ligation bonds.

Each `columns` entry contains both physical strand contributors. Its bases use
canonical `original` records with `segment`, `strand` and `reference_position`.
The global `ledger` contains every original or synthesized nucleotide once.
Every allocated ID is present either in final columns or in `removed`, and
every template ID resolves in the ledger even if that template was later
removed. New synthesis is never relabeled as germline by a subsequent join.

`max_operations` is shared across both scripts. `max_bases` counts the two
original strands of all three references and every newly synthesized base;
removal refunds none. Partial, incompatible, over-budget or unresolved products
raise errors and export no completed sequence.

### Bounded substrate contract

The input supplies three coding duplexes and four junction-facing hairpins.
The model does not simulate RSS recognition, RAG cleavage activation or a
detailed cleavage chronology. Both join orders are mechanistic controls under
this supplied-substrate contract, rather than a repertoire-wide claim about
their prevalence.

Completed ligation bonds and their two adjacent physical nucleotides are
protected. Processing that would consume or cross those bonds is outside this
bounded model and raises an explicit error. It is not classified as biologically
impossible in every repair mechanism. General flap/mismatch repair, reopening
completed joins, stochastic enzyme kinetics and biochemical calibration remain
extensions. Fully paired completion is required for the final VDJ product.

Hairpin opening and end processing are supported by primary biochemical work,
including sequence-dependent opening in
[Lu, Schwarz and Lieber (2007)](https://pubmed.ncbi.nlm.nih.gov/17932067/)
and coding-joint reconstitution in
[Lu and colleagues (2008)](https://pubmed.ncbi.nlm.nih.gov/18722175/).
The coupled short-substrate controls above are GenAIRR geometry tests, not
measurements from those studies. No universal enzyme rates or empirical
frequencies are assigned to the authored scripts.
