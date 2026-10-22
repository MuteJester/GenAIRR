# Reference

<p class="lead">A portable, hash-identified local germline reference built from
your FASTA files and an annotation sidecar. Build it with
<code>Reference.from_fasta</code>, store it with <code>write</code>, reload it
offline with <code>Reference.read</code>, and run it with
<code>Experiment.on(reference)</code>.</p>

See [Portable local references](../guides/portable-references.md) for the
workflow, the sidecar schema and the capability rules. Errors raised are
`GenAIRR.ReferenceArtifactError` (with a `kind`) and
`GenAIRR.ReferenceCapabilityError`.

## Common members

| Member | Purpose |
|---|---|
| `Reference.annotation_template(...)` | Build an annotation sidecar skeleton from FASTA files |
| `Reference.from_fasta(...)` | Build a reference from FASTA files and a sidecar |
| `Reference.read(path)` / `.write(path, replace=False)` | Load / store the portable file |
| `.report()` | Records, annotations, findings, defaults and capabilities |
| `.capabilities()` | Which features the evidence supports (`.available`, `.explain`, `.phase_support`) |
| `.require(*names, evidence="assertion")` | Raise `ReferenceCapabilityError` unless the capabilities hold |
| `.alleles`, `.records`, `.entities`, `.findings` | Inspect the catalogue |
| `.catalogue_hash`, `.execution_hash`, `.artifact_hash` | Content identity |
| `.to_json()` | The stored file text |

## Class reference

::: GenAIRR._reference.artifact.Reference
    options:
      show_source: false
      show_root_heading: true
      heading_level: 3
      members:
        - annotation_template
        - from_fasta
        - read
        - write
        - report
        - capabilities
        - require
        - to_json
