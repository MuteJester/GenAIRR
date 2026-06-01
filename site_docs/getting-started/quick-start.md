# Quick start

<p class="lead">A complete first simulation in 20 lines: bind to a
bundled cartridge, append a recombination pass with the productive
constraint, run 1,000 seeded records, and inspect the output.</p>

!!! tip "Your learning path"
    You're at step 2 of the **"I want to simulate sequences"**
    path. Next: [Your first AIRR record](first-airr-record.md)
    explains the fields in the output; then
    [The Experiment builder](../guides/experiment-builder.md)
    covers the full DSL surface.
    [See all paths →](../learn.md)

## The whole thing

```python
import GenAIRR as ga

result = (
    ga.Experiment.on("human_igh")     # bind to bundled human IGH cartridge
      .recombine()                    # append V(D)J recombination pass
      .productive_only()              # constraint-aware: only productive sequences
      .run_records(n=1000, seed=42)   # compile + run + project to AIRR records
)

print(len(result))               # 1000
rec = result[0]
print(rec["v_call"],
      rec["d_call"],
      rec["j_call"])              # 'IGHVF10-G38*04 IGHD2-15*01 IGHJ2*01'
print(rec["junction_aa"])         # 'CVKDDGNRGYCSGGSCYGRCCALDYWYFDLW'
print(rec["productive"])          # True
```

If this runs and prints, you're done with the first simulation.

## The mental model

Three nouns carry GenAIRR's surface. Once you have these mapped, the
rest of the API follows from them.

| Noun | Role | Example |
|---|---|---|
| **`Experiment`** | The pipeline. A fluent builder of biology + library + sequencing passes. | `ga.Experiment.on("human_igh").recombine().mutate(rate=0.05)` |
| **`DataConfig`** | The reference cartridge. Identity, allele catalogue, rules, empirical models — everything the engine needs to know about the biology it's simulating. | `ga.HUMAN_IGH_OGRDB`, or load via the string shortcut `"human_igh"` |
| **`SimulationResult`** | The output. A list-like wrapper around a batch of AIRR record dicts plus the underlying `Outcome` objects for advanced introspection. | `result = exp.run_records(n=100)`; `result[0]["v_call"]` |

Every method call on `Experiment` returns the same `Experiment`,
extended with one more step. The pipeline only actually runs when
you call `.run_records(...)` (or `.compile().run(...)` if you want
the compiled plan first). This means the DSL composes cleanly:
`recombine().mutate().pcr_amplify()` builds the plan, then
`run_records(n=..., seed=...)` executes it.

## What each step does

| Step | Effect |
|---|---|
| `Experiment.on("human_igh")` | Bind to the bundled human IGH cartridge. Other shortcuts: `"human_igk"`, `"human_igl"`, `"mouse_igh"`, `"human_tcrb"`. Pass a `DataConfig` instead of a string for a custom cartridge. |
| `.recombine()` | Append a V(D)J recombination pass — sample alleles, trim, fill NP1/NP2, assemble. Defaults to the cartridge's empirical models. |
| `.productive_only()` | Constraint-aware: the engine samples only choices that produce a productive sequence (in-frame junction, no stop codons, anchors preserved). No retry loops. |
| `.run_records(n=1000, seed=42)` | Compile the plan, run 1,000 seeded draws, project each into an AIRR-format record. Same seed → byte-identical output across runs and platforms. |

## What you can do next

The simulation's done; now you can:

- **[Inspect your first AIRR record](first-airr-record.md)** — what
  each field means and where it came from.
- **[Export the results](export-results.md)** — AIRR TSV, FASTA,
  FASTQ, paired-end FASTQ, pandas DataFrame.
- **[Validate the output](../validation/validate-records.md)** — run
  the postcondition validator to confirm every field is internally
  consistent.
- **[Build your own cartridge](../concepts/reference-cartridge.md)**
  — the four-plane cartridge model and how to author a custom
  reference from FASTA.

---

## Next step

→ [Your first AIRR record](first-airr-record.md) — read the output
field by field.
