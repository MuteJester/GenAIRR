# Reference evidence and catalogue quarantine

A reference cartridge can pass structural checks without having independently
verified coding boundaries or landmark positions. GenAIRR keeps those claims
separate. The bundled catalogue audit is a reproducible inspection of the
archived data, not a certification of its biological completeness.

## Authoring new references

`ReferenceCartridgeBuilder` and direct `VAllele` / `JAllele` construction use
the GenAIRR native resolver. No external annotation program is required.

- A declared IMGT-gapped V alignment maps conserved residue 104 to an ungapped
  coordinate and validates the cysteine codon. The result is **inferred** from
  that alignment declaration.
- A unique locus-compatible J F/W-G-X-G motif is an **inferred** candidate.
  Multiple candidates stay unresolved; the resolver does not select one
  arbitrarily. The J expectation is W for IGH and F for the other supported
  immunoglobulin/TCR loci.
- A bare ungapped V sequence does not supply enough positional evidence to
  identify the conserved cysteine. Its anchor stays unknown unless the author
  supplies a source-backed coordinate.
- IMGT FASTA header accession, species and F/ORF/P status are retained on each
  allele. Missing functionality becomes unknown. Cartridge-wide species stays
  unknown until explicitly declared with `infer_identity(species=...)`.
- Missing native functionality raises an error. Authoring no longer silently
  drops anchors when an obsolete resolver cannot be imported.

An anchor-derived modulo-three value is not stored as independent coding-phase
evidence. Bare FASTA input also does not establish complete genomic coding ends.
The new `reference_evidence` document records these limitations explicitly.

For curated annotations, use `builder.with_reference_evidence(...)`, keyed by
allele name. An optional `anchors={name: position}` supplies source-backed
zero-based landmark positions. Each evidence document has these fields:

```python
{
    "schema_version": 1,
    "sequence_sha256": "sha256:<SHA-256 of exact ungapped sequence bytes>",
    "coding_phase": 0,  # 0, 1, 2, or None when not independently established
    "anchor_status": "verified",  # verified / inferred / unknown / invalid
    "boundary_5": "complete",     # complete / partial / unknown
    "boundary_3": "complete",
    "source": "<accession/release or a versioned local evidence identifier>",
    "method": "<how this annotation was established>",
    "notes": ["<evidence and any limitations>"],
}
```

Complete/verified values are author claims that require evidence; do not copy
these example values onto unreviewed sequences. The bridge checks the schema,
sequence checksum, phase range and agreement with a usable anchor. Verified
landmark claims require a source identifier. It does not establish whether the
cited source proves the claim. A wrong checksum is an error, never a fallback to
legacy interpretation. Original-end processing excludes alleles whose relevant
coding boundary is explicitly partial, even when its Boolean acknowledgement
is true. Unknown boundary status can still be covered by that explicit caller
assertion; it does not become source-verified evidence. Evidence is available from
`refdata.allele_evidence_json(segment, allele_id)` and the cartridge manifest.

A verified landmark coordinate identifies its position even when the residue is
noncanonical: the junction remains defined and canonical productivity is false.
That position alone does not establish coding phase for a noncanonical source
residue; independent phase evidence is needed. Legacy or inferred noncanonical
landmarks remain unresolved. Productive sampling still requires canonical
landmark residues.

The native registry uses segment and allele name, so curated pool reindexing
cannot attach an annotation to another allele. Same-length synthetic novel
alleles retain coordinate evidence with explicit derivative provenance and a
new sequence checksum. They do not inherit a claim of source verification.

## Bundled audit and unresolved data

```python
from GenAIRR import data
status = data.audit_config("MOUSE_TCRA_IMGT")
print(status["status"], status["issue_counts"])
```

The registry covers all 106 archived cartridges. At this audit revision, 80
have structural or anchor findings and are labelled `quarantined`; 26 pass
those checks and are labelled `structurally_valid_evidence_unresolved`.
A finding can reflect an intentionally included pseudogene/ORF as well as a
bad annotation: quarantine is a statement about default strict admission,
not a claim that every reported allele is biologically impossible.

Default strict compilation rejects catalogue findings. The explicit
`allow_curatable_refdata()` diagnostic override permits curatable findings;
fatal structural problems remain errors. Inspection/loading of archived
cartridges remains available. Unknown or invalid reference landmarks cannot
become verified biological measurements simply by selecting the override.

The registry binds each result to both the archived file SHA-256 and the
DataConfig checksum. `audit_config` rejects stale file evidence. A manifest
includes a bundled audit only when its cartridge content matches the audited
checksum. Authoritative source snapshots supporting repairs to the known mouse
TCRA/TCRG anchor errors and mouse IGK/IGL structural errors are not bundled.
Their status remains unresolved; no replacement anchors or sequences were
invented, and no archived pickle was rewritten.

Reproduce the registry from the local GenAIRR installation:

```bash
python tools/audit_reference_cartridges.py --output reference_audit_v1.json
```

This does not estimate empirical recombination distributions. Structural
cartridges and fitted empirical models remain separate data products.

## Compatibility and trace identity

Legacy cartridges without evidence keep their original reference content-hash
encoding. Explicit evidence uses a versioned hash extension and participates in
reference identity. Python pickles with no new per-allele evidence remain
unchanged. New builder artifacts explicitly use `pickle4_roundtrip_v2` for the
Python DataConfig checksum: a local serialization roundtrip normalizes first-load
attribute-key interning before hashing. Archived `pickle4_v1` checksums retain
exactly their existing interpretation. Adding, changing or correcting evidence creates a different
reference artifact; old traces must still be interpreted against the original
reference and semantics, never silently retargeted to a curated replacement.
