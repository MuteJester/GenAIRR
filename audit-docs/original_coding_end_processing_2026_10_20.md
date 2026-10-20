# Original coding-end processing implementation

Branch: `fix/biological-metadata-correctness`. Date: 2026-10-20.

This implements the approved first stage: a finite model of retained original
coding-end outcomes, entirely inside GenAIRR. It adds no simulator, inference
service, subprocess dependency, or package dependency. Explicit strand repair,
hairpin-opening kinetics, and calibration of natural outcome frequencies remain
separate work.

## Configuration and model contract

```python
import GenAIRR as ga

experiment = (
    ga.Experiment.on("human_igh")
    .recombine(np1_lengths=[(n, 1.0) for n in range(7)],
               np2_lengths=[(n, 1.0) for n in range(7)])
    .process_coding_ends(
        v_3=[(-1, 1.0), (0, 2.0), (1, 1.0)],
        d_5=[(0, 1.0)], d_3=[(-2, 1.0), (0, 1.0)],
        j_5=[(0, 1.0)],
        reference_ends_are_coding_ends=True,
    )
)
```

Weights in this example are authored demonstrations. A positive outcome removes
reference bases, zero preserves the end, and a negative outcome retains an
original-tip palindromic extension. The boundary declaration is an explicit
input contract: the junction-facing reference ends must be the original complete
biological coding ends. GenAIRR does not infer this from a cropped sequence.
Copied P templates must have resolved A/C/G/T bases.

Omitted end distributions are fixed at zero. Choosing this mode replaces
automatic trim/P settings. Explicit trim settings conflict in either order.
Existing configurations retain the legacy model. New mode parameters and
version are part of plan/replay signatures, with new outcome-choice addresses.

The Python API builds a joint D table from products of canonical-end weights.
The native constructor also accepts authored joint rows directly, permitting
future correlated learned distributions. Opposite deletions may meet at an
empty D interval; crossing deletions are rejected. Extensions use original
template positions even when the opposite deletion removes those positions
from the retained reference body. D inversion reverses the complete processed
segment and swaps its physical P ends while preserving canonical coordinates.

Physical feasibility filters initial allele choices, even without productivity
contracts. Productive support and runtime emission share the pure materializer.
Constrained draws retain authored weights among choices admitting a completion.
This is local restriction to feasible support, without downstream completion-mass
weighting. Strict failures remain explicit; no zero-P fallback is fabricated.

Hand-built native plans must respect coding-order assembly/emission and finish
recombination before base or observation edits. Mixing legacy P additions into
the new model is rejected, including additions on another segment.

## Provenance and metadata

Standard germline calls and alignments keep the generating-reference semantics
established by the earlier assay fixes. AIRR trim fields remain nonnegative
canonical reference losses. P template ancestry does not make a P base a retained
germline alignment base. Combined NP fields include intervening N and P.

An optional `coding_end_processing` dictionary records model/version, the coding
boundary contract, canonical outcome choices, original allele identities,
retained P counts at joining, current observed P counts, and per-base origins.
Origins retain template position, source orientation, canonical and assembled
ends, original base, current base and current coding-oriented query position.
Removed bases keep their origins with null current positions. Mutation, masking,
indels, observation loss, forks, lineages and V replacement preserve this history.

`recombination_complete` distinguishes full assembly from partial native output.
On partial output, historical P counts are provisional emitted-to-date counts.
Completion stays true when later clipping removes all observed bases. No nick
offset, physical strand-resection amount or polymerase history is invented.

CSV/TSV export stores the dictionary as a JSON cell. Its explicit zero-based,
coding-oriented coordinate convention survives strict AIRR conversion of the
standard columns. The generated metadata field is reserved against user stamps.

The native validator independently checks geometry, source templates, original
bases, per-choice template multiplicities, historical counts and current P flags.
It detects missing emissions even when origin records, observed pool bases and
legacy counters are removed together. Public validation additionally compares
the returned metadata against its outcome.

## Calibration boundary

Finite weighted outcomes are usable now; learning them is future GenAIRR-native
work. Existing empirical trims are never reinterpreted as physical erosion.
Existing empirical N inputs are not automatically refitted: historical NP
estimates can contain P unless their training provenance establishes otherwise.
Current N estimators reject explicit nonzero P training provenance.

The model restricts same-end net deletion and retained original-tip P to separate
outcomes. This restriction does not assert that recessed palindromes, strand
restoration or other repair pathways are impossible. A future paired-duplex
model can express those pathways through additional explicit mechanisms.

## Verification

- Complete Rust debug suite: **1,479 passed**.
- Complete Rust release suite: **1,478 passed**. The difference is one existing
  debug-only panic regression. Each profile retains two ignored integration
  tests and two ignored doctests.
- Cross-locus stress: **660 records across 22 profiles** passed native
  validation, independent template/count/position reconstruction, and true replay.
  Six additional profiles were rejected explicitly: four profile combinations
  have fatal mouse IGK/IGL reference integrity problems, and two productive mouse
  TCRA/TCRG profiles have invalid anchors. The artifact preserves every rejection.
- Complete Python correctness suite: **2,636 passed**, eight performance-budget
  cases deselected, 169 warnings, 461.40 seconds. The initial run exposed two
  whitespace-sensitive source assertions after rustfmt; both assertions now
  accept equivalent Rust formatting, and the final full rerun passes.
- `git diff --check` passed.

The [verification artifact](original_coding_end_verification_2026_10_20.json)
records exact counts, logs and extension identity. Its extension SHA-256 matches
the independently recorded cross-locus stress artifact.

Independent small-sequence oracles cover VJ/VDJ, uniform/Markov N support,
orientation, joint D outcomes and committed outcome prefixes. Additional tests
cover filtered weights, invalid schedules, linked genotype, true trace replay,
all-base masking, complete read loss, receptor revision, clones and JSON export.

The reproducible cross-locus sweep is
`python -m tools.metadata_biology_audit.original_end_stress`; its artifact is
`original_coding_end_stress_2026_10_20.json`. It tests authored distributions and
metadata/replay consistency, independently reconstructing source template bytes.
It is not an empirical validation of natural hairpin-processing probabilities.

## Existing issues found during the sweep

The [separate evidence note](original_end_existing_limits_2026_10_20.md)
records two problems outside the new processing kernel. Under an explicit
invalid-reference override, shared productivity labeling can accept bad annotated
anchors; this also affects unconstrained new-mode output. Default validation
rejects these cartridges, and new-mode productive support rejects the bad anchors
even after the override. Passing ancestry/coordinate checks does not establish
biological correctness of those invalid annotations.

The strict productive-indel sampler can also exhaust its bounded retry budget
despite nonempty valid support. This pre-existing behavior was reproduced in both
end models. The sweep keeps indels in unconstrained profiles and uses SHM only
for productive IG profiles; TCR profiles use no SHM. This validation does not
claim to resolve that independent indel sampler limitation.
