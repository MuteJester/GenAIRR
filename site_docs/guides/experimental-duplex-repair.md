# Experimental duplex repair

`GenAIRR.experimental.DuplexRepair` is an explicitly authored local experiment
for two complete coding-end hairpins. The Rust kernel represents both strands
and preserves their separate origins. This experimental API does not yet feed
the standard AIRR experiment pipeline or simulate coupled successive D joins.

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
history. The kernel currently handles a single local join, compatible annealing
and template completion. Flap processing, mismatch resolution, general repair
kinetics and experimental calibration remain explicit extensions.

Verification covers nucleotide conservation during opening, paired ancestry,
opposing N additions sharing a product column, resection/restoration, strand
symmetry, exact script replay and rejection of invalid primer/template geometry.
