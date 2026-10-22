# Portable local references

<p class="lead">A <code>GenAIRR.Reference</code> turns your own germline FASTA files and
a small annotation sidecar into a single portable file. The file records which
alleles you gave, what you claimed about them, and how strongly each claim is
backed. It can be moved to another machine and run offline, and GenAIRR
reports which simulation features your evidence supports.</p>

Use it when the bundled cartridges do not cover your species, locus or allele
set. For the cartridge route (with fitted models and the full builder), see
[Build a reference cartridge](build-reference-cartridge.md).

## The workflow

1. Generate an annotation **template** from your FASTA files.
2. Fill in what you know; leave the rest as `null`.
3. **Build** the reference and read its `report()`.
4. **Write** the portable file.
5. **Read** it back (offline, anywhere) and **run** it.
6. **Validate** the records against the reference.

The examples below assume four files in the working directory: `v.fasta`,
`d.fasta`, `j.fasta` and `annotations.json`. Omit `d` for a VJ chain
(light chain, TRA, TRG).

### 1. Generate a template

```python
from GenAIRR import Reference

template = Reference.annotation_template(
    v="v.fasta", d="d.fasta", j="j.fasta",
    reference={
        "name": "my_igh", "version": "1", "locus": "IGH", "chain_type": "vdj",
        "taxon": {"scientific_name": "Homo sapiens", "taxonomy_id": 9606},
        "build_date": "2026-10-21",
    },
)
print([(a["segment"], a["label"]) for a in template["alleles"]][:3])
print(sorted(template["alleles"][0]))
```

The template lists every allele with each annotation set to `null`
(unknown), already carries the `sequence_sha256` of each sequence, and adds
`suggestions` that are candidates only (for example a unique J motif
position). A suggestion is never applied for you: copy it only if your evidence
confirms it. If your V sequences are IMGT-gapped, pass
`v_alignment="imgt_gapped"` (see [Conventions](#conventions)).

Edit the JSON into your `annotations.json` (the template is a plain `dict`;
write it with `json.dump(..., indent=2)` and edit it).

### 2. Build and inspect

```python
ref = Reference.from_fasta(
    v="v.fasta", d="d.fasta", j="j.fasta", annotations="annotations.json"
)
print(ref.report())
print(ref.capabilities().available())
```

`print(ref.report())` shows the record accounting (accepted, duplicate,
rejected, with reasons), the per-allele annotations with their support, named
entities, findings, the defaults profile and the capabilities. `report()` also
has `to_dict()` for JSON. `ref.records`, `ref.alleles`, `ref.entities`
and `ref.findings` expose the same information as objects.

### 3. Write and read

```python
ref.write("my_igh.genairr-ref.json")
# ... move the file anywhere, even to a machine without network access ...
again = Reference.read("my_igh.genairr-ref.json")
assert again.execution_hash == ref.execution_hash
```

`write` refuses to overwrite an existing file unless you pass
`replace=True`. The file contains no local paths and does not depend on
where it was built. Written files follow your umask (for example mode 0644).
`Reference.read` verifies the file's integrity and raises
`GenAIRR.ReferenceArtifactError` for a malformed or modified file.

Three hashes identify a reference: `catalogue_hash` (the alleles and the
claims about them), `execution_hash` (what the simulator actually runs) and
`artifact_hash` (the stored file). `ref.to_json()` returns the file text.

### 4. Simulate and validate

```python
from GenAIRR import Experiment

result = Experiment.on(again).recombine().run_records(n=200, seed=7, validate_records=True)
print(len(result), result.reference_content_hash == again.execution_hash)
print(result.validate_records(again).ok)
```

A path is never read implicitly. `Experiment.on("my_igh.genairr-ref.json")`
(any string ending in `.json`) or `Experiment.on(Path(...))` raises `TypeError`,
and any other unknown name raises `ValueError`; both messages point you to
`GenAIRR.Reference.read(path)`.

## The annotation sidecar

The sidecar is a JSON file with `"schema": "genairr_reference_annotations"`
and `"schema_version": 1`. It has a `reference` header (name, version,
locus, chain type, taxon, build date), a `conventions` block and an `alleles`
list:

The citation below is illustrative. `support: "source"` requires a real, checkable source.

```json
{
  "segment": "V",
  "label": "IGHV1-2*02",
  "sequence_sha256": "sha256:e6619303a938c2c78fe5586e643b9abadd7e4b0512f0b2bfc8bfef5659e3d157",
  "gene": null,
  "functionality": {"value": "functional", "support": "source",
                    "source": "IMGT/GENE-DB release 202543-2", "method": "IMGT functionality field"},
  "anchor": {"value": 285, "support": "assertion", "source": null,
             "method": "conserved cysteine located by hand"},
  "coding_phase": {"value": 0, "support": "assertion", "source": null,
                   "method": "derived from the anchor", "derivation": "anchor_derived"},
  "boundary_5": null,
  "boundary_3": {"value": "complete", "support": "assertion", "source": null,
                 "method": "sequence ends at the germline end"},
  "notes": []
}
```

### Support levels

Every annotation is either `null` (unknown) or an object with `value`,
`support`, `source` and `method`.

| `support` | Meaning | Requirements |
|---|---|---|
| `"source"` | The value comes from a citable source (database release, paper, assay). | `source` must be a non-blank citation. |
| `"assertion"` | You state it, with no external source. | `source` may be `null`. |

An **assertion is not source support.** GenAIRR keeps the difference all the
way through to `capabilities()` and to the report, and never promotes an
assertion to source. Use `method` to say how you arrived at the value.

### Derivations

`coding_phase` additionally carries a `derivation`: `"independent"` when
the phase is supported on its own, or `"anchor_derived"` when you computed
it from the anchor. An **anchor-derived phase is not independent evidence**:
it can only be as strong as the anchor it came from, so it does not upgrade
any capability. `ref.capabilities().phase_support` reports which kind the
reference has (`"none"`, `"anchor_derived"`, `"independent_assertion"` or `"independent_source"`). Anchors
derived from a declared IMGT alignment are recorded with the derivation
`"imgt_alignment"`.

### The sequence digest

`sequence_sha256` binds each annotation to one sequence. It is the SHA-256
of the sequence after normalization: uppercase, gap characters removed, only
A/C/G/T, formatted as `sha256:<hex>`. Case and IMGT gap dots therefore do not
change the digest, but any change to the bases does, so an annotation can never
silently attach to a different sequence. The template fills the digest in for
you.

### Labels

An allele label is taken from the FASTA header and kept exactly as written; it
is never renamed. A label must be 1 to 256 characters with no whitespace,
control or format characters and no comma (commas separate multiple calls in
AIRR fields). A record whose label breaks the rule is rejected with a reason in
the report. Identical sequences under two labels, or the same label twice, are
reported as findings or duplicates rather than merged silently.

### Conventions

`conventions.coordinates` is `"zero_based_half_open"`. IMGT alignment is
declared only by `conventions.v_alignment: "imgt_gapped"` (the default is
`"none"`). When declared, GenAIRR derives V anchors and the FWR/CDR
subregions from the gapped V sequences; GenAIRR never guesses an alignment from
the sequence.

## Capabilities

`ref.capabilities()` answers which features your evidence supports:

| Capability | Needs |
|---|---|
| `catalogue_inspection` | Always available. |
| `observed_span_simulation` | At least one accepted V and J allele (and D for a VDJ chain). |
| `source_annotated_junction` | An accepted anchor for every V and J allele. |
| `canonical_productive` | Accepted anchors, with the V anchor codon encoding C and the J anchor codon encoding the locus residue (W for IGH, F otherwise). |
| `original_coding_ends` | `boundary_3: "complete"` on V and `boundary_5: "complete"` on J (and complete D ends). |
| `complete_reference_span` | Complete `boundary_5` and `boundary_3` on every V and J allele. |
| `genotype_constrained` | Not yet supported for portable references. |
| `duplex_family_replay` | Not yet supported for portable references. |

Each available capability carries an evidence level, `"source"` or
`"assertion"`, equal to the weakest support among the claims it relies on.
`available(evidence="source")` lists only source-backed capabilities. Use
`require` to fail early:

```python
ref.require("observed_span_simulation", "canonical_productive")  # raises if unavailable
ref.require("observed_span_simulation", evidence="source")       # needs source-backed evidence
```

`require` raises `GenAIRR.ReferenceCapabilityError` with the full
explanation when a capability is missing or too weakly supported.

### Worked example: a light chain with partial evidence

This IGK example has two V and two J alleles, but only
`kappa_v_a` and `kappa_j_1` have anchors.

```python
igk = Reference.from_fasta(v="igk/v.fasta", j="igk/j.fasta", annotations="igk/annotations.json")
print(igk.capabilities().explain())
```

```text
available   catalogue_inspection
available   observed_span_simulation
unavailable source_annotated_junction: 2 requirement(s) unmet
- V:kappa_v_b: accepted anchor — add an "anchor" annotation (support "source" with a citation, or "assertion")
- J:kappa_j_2: accepted anchor — add an "anchor" annotation (support "source" with a citation, or "assertion"); candidate at 8 (unique J motif) — copy it only if your evidence confirms it
unavailable canonical_productive: 2 requirement(s) unmet
- V:kappa_v_b: accepted anchor — add an "anchor" annotation (support "source" with a citation, or "assertion")
- J:kappa_j_2: accepted anchor — add an "anchor" annotation (support "source" with a citation, or "assertion"); candidate at 8 (unique J motif) — copy it only if your evidence confirms it
unavailable original_coding_ends: 4 requirement(s) unmet
- V:kappa_v_a: accepted boundary_3 == "complete" — add "boundary_3": {"value": "complete", ...} with its support
- V:kappa_v_b: accepted boundary_3 == "complete" — add "boundary_3": {"value": "complete", ...} with its support
- J:kappa_j_1: accepted boundary_5 == "complete" — add "boundary_5": {"value": "complete", ...} with its support
- J:kappa_j_2: accepted boundary_5 == "complete" — add "boundary_5": {"value": "complete", ...} with its support
unavailable complete_reference_span: 8 requirement(s) unmet
- V:kappa_v_a: accepted boundary_5 == "complete" — add "boundary_5": {"value": "complete", ...} with its support
- V:kappa_v_a: accepted boundary_3 == "complete" — add "boundary_3": {"value": "complete", ...} with its support
- V:kappa_v_b: accepted boundary_5 == "complete" — add "boundary_5": {"value": "complete", ...} with its support
- V:kappa_v_b: accepted boundary_3 == "complete" — add "boundary_3": {"value": "complete", ...} with its support
- J:kappa_j_1: accepted boundary_5 == "complete" — add "boundary_5": {"value": "complete", ...} with its support
- J:kappa_j_1: accepted boundary_3 == "complete" — add "boundary_3": {"value": "complete", ...} with its support
- J:kappa_j_2: accepted boundary_5 == "complete" — add "boundary_5": {"value": "complete", ...} with its support
- J:kappa_j_2: accepted boundary_3 == "complete" — add "boundary_3": {"value": "complete", ...} with its support
unavailable genotype_constrained: not yet supported for portable references
unavailable duplex_family_replay: not yet supported for portable references
phase support: none
```

Each unmet requirement names the allele, the missing annotation and what to
add. The `candidate at 8` hint is a suggestion only. The reference still
simulates (`observed_span_simulation`), but a run on it cannot decide
productivity for reads that use `kappa_v_b`: such records report
`productive` as `None` with `coding_phase_unresolved` in
`productivity_reasons`.

## Running a portable reference

`Experiment.on(reference)` accepts a `Reference` anywhere a config name or
`DataConfig` is accepted. Because a portable reference carries no fitted
models, the run uses the explicit **unfitted defaults profile**
`genairr.unfitted_defaults.v1`:

- uniform allele usage over the reference's alleles;
- NP1 and NP2 lengths uniform over 0 to 6;
- uniform NP bases (A, C, G, T);
- exonuclease trimming disabled;
- P-nucleotides disabled.

These are placeholders, not fitted to data and not a model of any real
repertoire. `exp.describe()` names the portable reference and the profile, and
`ref.report()` lists the defaults. Every result records
`reference_content_hash`, the content hash of the reference it ran against.

Override any default per run with the usual arguments:

```python
exp = (
    Experiment.on(ref)
    .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(i, 1.0) for i in range(7)])
    .trim(v_3=[(2, 1.0)])
)
print(exp.describe())
records = exp.run_records(n=20, seed=1)
```

`result.validate_records(reference)` checks the records against the selected
reference. Passing a different reference raises `ValueError`. This applies to
bundled references too: every result now records the hash of the reference it
ran against. A genotype run with novel alleles must be validated against
`compiled.refdata` (the effective reference).

## Not yet supported

| Use | Behaviour |
|---|---|
| `with_genotype(...)` or `run_cohort(...)` on a `Reference` | Raises `GenAIRR.ReferenceCapabilityError` (a subclass of `ReferenceArtifactError`, `kind == "capability"`): "genotype- and cohort-constrained runs are not yet supported for portable references". |
| Family-replay capture (`capture_family_replay=True`) | Rejected before sampling with `GenAIRR.experimental.FamilyReplayError` (category `"unsupported"`): "family-replay capture is not yet supported for portable references". |

## Errors

`GenAIRR.ReferenceArtifactError` has a `kind`:

| `kind` | Raised when |
|---|---|
| `malformed` | The sidecar or file is not valid JSON of the expected shape (the error names the path inside the document). |
| `unsupported_schema` | The schema name or version is not understood. |
| `unsupported_convention` | A declared convention (coordinates, alignment) is not supported. |
| `conflict` | An annotation contradicts the catalogue, for example a sidecar entry whose digest matches no accepted sequence. |
| `integrity` | A value is inconsistent with its sequence (for example an anchor out of bounds) or a stored file fails its checksum. |
| `capability` | A requested use is not supported (raised as the subclass `ReferenceCapabilityError`). |
| `empty_selection` | No allele of a segment was accepted, so nothing can run. |
| `io` | A file cannot be read or written (including refusing to overwrite an existing file). |

See the [Reference API](../reference/reference.md) for every method.
