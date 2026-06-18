# Sample only productive sequences

<p class="lead">`productive_only()` constrains V(D)J sampling so the engine only ever
produces <strong>productive</strong> rearrangements - in-frame junction, no stop
codon in the junction, and the conserved V (Cys) and J (Trp/Phe) anchors intact.
It is a chainable, order-independent step on the <code>Experiment</code>.</p>

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")
      .recombine()
      .productive_only()          # constrain to productive rearrangements
      .run_records(n=1000, seed=7)
)
# every record: rec["productive"] is True
```

## What "productive" means here

A rearrangement is productive when its assembled coding sequence could yield a
functional receptor. GenAIRR enforces the **productive triad**:

- the junction is **in-frame** (the V(D)J join preserves the reading frame),
- there is **no stop codon in the junction**, and
- the conserved anchors are intact - the **V cysteine** and the **J
  tryptophan/phenylalanine** that bound the CDR3.

## How it works - constrain before propose

`productive_only()` is **not** a reject-and-resample loop that throws away bad
draws after the fact. GenAIRR follows a **constrain-before-propose** discipline:
the constraint narrows the *candidate support at each draw*, so out-of-frame or
stop-bearing candidates are never proposed in the first place. The engine samples
only from the admissible set.

This is why it is efficient and why the productive fraction jumps from the
unconstrained baseline (heavy chains are productive only a minority of the time
by chance - about [18.5% in the demo](../demo.md)) to 100%. It also composes
with the other sampling constraints (e.g.
[`restrict_alleles`](recombination-junction.md#controlling-the-allele-universe)) -
the engine intersects all active constraints into one admissible support, so the
choice of where in the chain you call `productive_only()` doesn't matter.

## Interaction with SHM

The same discipline extends to somatic hypermutation. When `productive_only()` is
in the pipeline, the [mutation pass](shm-targeting.md#interaction-with-productive_only)
never proposes a substitution that would break the productive triad - the
constraint masks the per-site proposal support, so a mutation that would
introduce a junction stop or destroy an anchor simply can't be drawn. `productive`
stays `True`, and `n_mutations` counts only the substitutions that landed and
survived the mask.

## What you get on the record

Every record carries the AIRR `productive` field; with `productive_only()` it is
`True` for all of them. The output is productive *by construction* - you do not
need to filter post-hoc, and [`validate_records`](../validation/validate-records.md)
will confirm the productive invariants hold for every row.

## When to use it

- **Realistic functional repertoires** - most expressed immunoglobulin (BCR)
  sequences are productive; `productive_only()` matches that without discarding
  samples.
- **Benchmarks that assume productive input** - aligners/annotators are usually
  evaluated on productive reads.

Leave it off when you specifically want the unconstrained mix (including
non-productive rearrangements) - e.g. to study out-of-frame rates or to
stress-test a tool on non-productive input.

!!! warning "Immunoglobulin (BCR) only"
    The productive contract's anchor checks assume BCR semantics, so it is an
    immunoglobulin feature. TCR refdata accepts the call but raises a
    `ValueError` at compile time - don't add `productive_only()` to a TCR
    pipeline.

The constraint is implemented as a bundle of engine contracts (anchor
preservation, junction frame, no-junction-stop); the contributor-facing details
live in the [architecture guide](../architecture/index.md).
