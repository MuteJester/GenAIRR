# Biological metadata correctness audit

Audit date: 2026-10-20. Source commit: `1ededa676e8904430f7ea0bdcf03092de4888a41`.

This records the baseline audit at that commit. Implementation work on
`fix/biological-metadata-correctness` changes the behavior described below;
the original JSON remains baseline evidence.

**GenAIRR's current metadata cannot be treated as universally correct biological ground truth.** Reproducible failures affect functionality, translation, germline alignments, recombination provenance and AIRR interchange. Some occur in ordinary built-in human IGH pipelines. Several others require optional features or deliberately small custom references that isolate a particular calculation.

The audit rebuilt the Rust extension against this checkout and used the editable Python package. No simulator implementation was changed. The accompanying probes and JSON preserve concrete examples for fixes and regression tests.

## Evidence and scope

- [Runnable audit](../tools/metadata_biology_audit/run_audit.py): 19 diagnostic probes covering the findings below, including the separately identified modeling and interpretation concerns.
- [Recorded results](biological_metadata_audit_2026_10_20.json): exact observations, expected behavior, source pointers and current validator responses.
- [Independent alignment sweep](../tools/metadata_biology_audit/coordinate_sweep.py) and [results](biological_metadata_coordinate_sweep_2026_10_20.json): 1,200 human IGH/IGK records across six pipeline profiles.

Run from the repository root after building the current extension:

```bash
VIRTUAL_ENV="$PWD/.venv" .venv/bin/maturin develop --release -m engine_rs/Cargo.toml
.venv/bin/python tools/metadata_biology_audit/run_audit.py --output /tmp/metadata-audit.json
.venv/bin/python tools/metadata_biology_audit/coordinate_sweep.py --output /tmp/metadata-sweep.json
```

The diagnostic command's exit status indicates successful execution. A result with `reproduced: true` means the discrepancy remains present; it does not mean the record is correct. Reproduction counts are cases, not independent root causes.

Coordinates below use the Python API's zero-based, half-open convention unless explicitly discussing exported AIRR TSV. Small synthetic references are mechanical test fixtures, not claims about complete natural receptors.

P1 means a high-priority error in downstream biological interpretation or alignment data. P2 means a narrower correctness or interoperability defect. Neither priority is a claim about prevalence across every species or parameter setting.

## Confirmed findings

### 1. P1 — Internal stops outside the junction are ignored

```python
ga.Experiment.on("human_igh").recombine().productive_only() \
    .mutate(rate=0.03).run_records(n=1, seed=42)[0]
```

The first record contains an in-frame V-region `TAA` at nucleotide 243, but reports `productive=True` and `stop_codon=False`. Its junction contains no stop. This is not merely a disagreement with the returned translation: an independent genetic-code table detects the stop.

The problem also occurs without mutation. Row 119 of `recombine().productive_only().run_records(n=120, seed=42)` retains the reference allele `IGHVF9-G32*02`'s `TAG` at `[282,285)`, immediately before the Cys anchor. The reference has no assigned functional status. It is still labeled productive.

In two fixed-seed batches of 1,000 records, the independent checks found 23 such clean records and 292 records with the 3% SHM setting. These are demonstrations under these exact parameters, not estimated population rates.

An associated branch forces `stop_codon=False` whenever the junction is out of frame. Plain recombination, seed 42, row 1 contains V-frame stops but follows this branch.

**Cause:** [builder.rs](../engine_rs/src/airr_record/builder.rs), lines 461–485, checks only the junction and suppresses stop reporting in the out-of-frame branch. The validator repeats that interpretation.

**Biological consequence:** nonfunctional coding sequences can enter productive-only benchmarks. AIRR defines `stop_codon` over the aligned sequence and recommends productivity requiring an open reading frame without internal stops. [AIRR schema](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)

**Fix direction:** derive stop status from the retained coding sequence in its correct biological frame; make productivity, constraints and validation use that definition. Retain pseudogene/nonfunctional references for unrestricted simulation, but do not equate an intact junction with a functional rearrangement.

### 2. P1 — Upstream frameshifts are hidden by changing the translation origin

A one-base insertion before the V Cys anchor shifts the entire `sequence_aa` translation start to keep the displaced anchor in frame.

Minimal example:

| Field | Value |
|---|---|
| Observed sequence | `ATCGAAATGTAAAAAATGGGGG` |
| V CIGAR | `2M1I10M` |
| Reported protein | `SKCKKWG` |
| Translation from the retained V coding origin | `IEM*KMG` |
| Reported productivity | `True` |

A built-in human IGH example also reproduces: `recombine().polymerase_indels(count=1, insertion_prob=1).run_records(n=20, seed=42)[19]`, with V CIGAR `136M1I163M`.

**Cause:** [builder.rs](../engine_rs/src/airr_record/builder.rs), lines 501–503, uses `junction_start % 3` for the entire translation. Junction length alone also determines `vj_in_frame`.

**Consequence:** the reported protein is wrong, and internal frameshift damage can be labeled productive. IMGT identifies frameshift mutations as an unproductive feature. [IMGT functionality](https://imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html)

**Fix direction:** track coding phase from the retained reference origin, distinguishing missing terminal coverage from internal indels. If corrected translations are offered, expose the correction explicitly. The existing curated test at `engine_rs/src/compiled/tests/curated.rs:821` currently endorses productivity after an upstream insertion and needs its biological expectation reviewed.

### 3. P1 — P nucleotides move to the end of the alignment

The fixture assembles V, three V-end P bases, three N bases, two J-end P bases and J:

```text
sequence:           AAACCCTGTACACATAATTTGGG
sequence_alignment: AAACCCTGTCATTTTGGGACAAA
```

All five P bases are emitted after J in the alignment. There are no query deletions here: removing alignment gaps should reproduce the observed sequence in order.

**Cause:** P bases are absent from `sim.sequence.regions`; [walk/mod.rs](../engine_rs/src/airr_record/walk/mod.rs), lines 115–121, appends every uncovered base after traversing all regions, assuming uncovered bases are trailing indel insertions.

**Consequence:** the alignment is not an alignment of the reported molecule.

**Fix direction:** render alignment columns in pool order and classify P bases in their actual junction positions.

### 4. P1 — Combined N/P annotations omit P bases

With built-in human IGK references, forced V-end P length 2, J-end P length 1, no N addition and no trimming, the three bases between V and J exist in the molecule but `np1=""` and `np1_length=0`.

**Cause:** [builder.rs](../engine_rs/src/airr_record/builder.rs), around line 191, reads only the structural NP region. P additions are outside it.

**Consequence:** junction decomposition and insertion lengths disagree with the actual sequence. AIRR's NP fields describe combined N/P sequence between the assigned segments. [AIRR schema](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)

**Fix direction:** derive these fields from the resolved inter-segment interval, including P bases. Keep mechanistic N-only and per-end P provenance in separately named fields.

### 5. P1 — D inversion produces an impossible germline alignment

An untrimmed, unmutated D allele `AACGTA` is correctly assembled as `TACGTT` under inversion. Its metadata nevertheless reports:

```text
d_cigar = 5D6M5D
d_identity = 0.0
n_indels = 0
```

The CIGAR consumes 16 reference positions even though the allele and reported germline span contain only six.

**Cause:** [segment_projection.rs](../engine_rs/src/airr_record/segment_projection.rs), around line 176, and [walk/mod.rs](../engine_rs/src/airr_record/walk/mod.rs), around line 217, assume increasing forward reference positions and use forward allele bases.

**Fix direction:** make the alignment walker, coordinates, identity and CIGAR consistently honor the D orientation. A provenance flag alone does not repair the alignment.

### 6. P1 — Reverse-strand output changes the receptor's biological annotations

Otherwise identical forward and reverse observations report junction proteins `CTF` and `KCT`. V query coordinates change from `[0,9)` to `[9,18)`, while the sequence/germline alignment stays forward.

**Cause:** [sequence.rs](../engine_rs/src/airr_record/sequence.rs), lines 18–44, reverse-complements junction/NP strings, flips query coordinates and translates the antisense strings, while leaving other annotation fields forward.

AIRR requires annotations for `rev_comp=True` to use the reverse complement of `sequence`. [AIRR schema](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)

**Fix direction:** keep biological annotations in the coding orientation and represent the observed read's strand consistently. Test a forward/reverse pair for invariant junction protein, calls and coding-orientation alignment coordinates.

### 7. P1 — Clonal descendants lose founder trims and acquire fictitious deletions

Compare these two pipelines with seed 42:

```python
ga.Experiment.on("human_igh").recombine() \
    .clonal_repertoire(n_clones=1, max_size=1)

ga.Experiment.on("human_igh").recombine() \
    .clonal_repertoire(n_clones=1, max_size=1).sequencing_errors(count=0)
```

The returned molecules are identical, but:

| Field | No post-fork step | Zero-error post-fork step |
|---|---:|---:|
| `d_trim_3` | 4 | 0 |
| `j_trim_5` | 1 | 0 |
| `d_cigar` | `27M` | `27M4D` |
| `j_cigar` | `52M` | `1D52M` |
| `d_identity` | 1.0 | 0.8709677 |
| `j_identity` | 1.0 | 0.9811321 |

**Cause:** [repertoire.py](../src/GenAIRR/_compiled/repertoire.py), lines 123–130, projects a descendant whose trace lacks the founder's recombination choices. [builder.rs](../engine_rs/src/airr_record/builder.rs), lines 80–111, reads trims from that incomplete trace and defaults to zero. The supported, deprecated `expand_clones` path exhibits the same issue. This reproduction does not establish the same defect in `clonal_lineage`, which has separate provenance merging.

**Fix direction:** read final allele-instance trims from persistent state or merge the required founder provenance before projection.

### 8. P1 — V replacement exports the old trim

Replacing a retained nine-base V with the first nine bases of a twelve-base V correctly sets the new internal trim to three. The record still exports `v_trim_3=0`, `v_germline_end=12` and `v_cigar="9M3D"`.

**Cause:** revision records the replacement trim separately and updates the allele instance, while [builder.rs](../engine_rs/src/airr_record/builder.rs), around line 83, reads the original recombination choice.

**Consequence:** intentional recombination trimming is mislabeled as deletion damage; reference span and identity are affected.

**Fix direction:** project the final allele-instance assignment while preserving the original decisions in explicitly historical provenance.

### 9. P2 — A shorter called allele can produce hybrid germline sequence

A valid custom-reference fixture starts from a nine-base V and mutates toward a six-base alternative. It reports `v_call="vShort*01"`, but `v_germline_end=9` and `v_cigar="9M"`. The final three germline bases are copied from the originally sampled allele, creating a reference sequence that is not the reported allele.

**Cause:** [walk/mod.rs](../engine_rs/src/airr_record/walk/mod.rs), lines 246–249, falls back to `nuc.germline` when the called allele is too short. The projection kernel does the same.

**Scope:** confirmed in a synthetic high-mutation fixture. Two additional built-in human IGH batches, 1,000 records each at mutation rates 0.03 and 0.1, found no reference-bound or reference-substring failures despite some allele-call switches.

**Fix direction:** require a valid alignment to each candidate's own reference; never fill an out-of-bounds called-reference position from a different allele.

### 10. P2 — Complete read loss discards its realized count

Removing an 18-base molecule reports `end_loss_5_length=0`, even though 18 bases were removed. [builder.rs](../engine_rs/src/airr_record/builder.rs), lines 75–77, returns before filling the counter. Unlike most findings here, the existing validator detects this empty-read discrepancy.

**Fix direction:** preserve realized event counts even when the observation becomes empty. Partial end-loss alignment conventions raise a separate interpretation concern described below.

### 11. P2 — CIGARs omit required leading offsets

For plain human IGH recombination, seed 42:

```text
d_sequence_start = 309; d_cigar = 27M
j_sequence_start = 336; j_germline_start = 1; j_cigar = 52M
```

The respective AIRR prefixes should be `309S` and `336S1N`. AIRR requires leading query clipping and reference skipping when those starts are nonzero. [AIRR CIGAR specification](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#cigar-specification)

**Cause:** [segment_projection.rs](../engine_rs/src/airr_record/segment_projection.rs), around line 338, serializes only the segment's internal operation runs.

**Consequence:** consumers reconstructing an alignment from CIGAR alone receive the wrong placement.

### 12. P2 — TSV export is not AIRR compliant, including strict export

`to_tsv()` defaults to zero-based starts. `to_tsv(airr_strict=True)` adjusts starts, but still writes booleans as `True`/`False`, instead of AIRR's `T`/`F`. [AIRR data values](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#data-values)

**Cause:** [_result_export.py](../src/GenAIRR/_result_export.py), lines 235 and 490–508, makes coordinate conversion optional and passes native booleans to `csv.DictWriter`.

**Scope:** zero-based coordinates are a documented, reasonable Python API convention. The defect is presenting the default TSV as AIRR output, and the incomplete behavior of strict interchange. Strict export also does not repair the CIGAR/orientation issues above.

### 13. P1 — Validation does not inspect the emitted record

After generating a record, replace its sequence with `AAA`, set its `productive` value to a non-boolean string and set `junction_start=999999`. `result.validate_records(refdata)` still returns `ok=True`.

The public generation path also accepts:

```python
ga.Experiment.on("human_igh").recombine() \
    .with_metadata(sequence="AAA", productive="wrong") \
    .run_records(n=1, seed=42, validate_records=True)
```

**Cause:** [_result_validation.py](../src/GenAIRR/_result_validation.py), lines 242–244, passes only the sequence identifier to `outcome.validate_record`; it validates a newly projected record, not the supplied dictionary. Metadata stamps are applied after validation and can overwrite biological fields.

**Consequence:** the validation guarantee does not cover what consumers actually receive. This is separate from validators sharing the same incorrect biological definitions as projection.

**Fix direction:** validate the actual emitted fields against independently derived expectations, check record/outcome cardinality, reserve generated field names against metadata collisions, and validate after all allowed transformations.

### 14. P2 — Revision can report a replacement with the same V allele

The custom two-V fixture, seed 1, returns `receptor_revision_applied=True` while the original and final truth V calls are identical.

**Cause:** the non-genotype replacement sampler in [receptor_revision.rs](../engine_rs/src/passes/receptor_revision.rs), lines 427–438, filters by length but does not exclude the current allele.

**Consequence:** this disagrees with the documented operation of replacing V with a different germline allele and can overcount actual replacements.

**Fix direction:** exclude the existing assignment, or explicitly distinguish an attempted revision draw from an effective replacement.

## Modeling and interpretation concerns requiring decisions

### P nucleotides are reconstructed after deleting coding-end bases

The simulator permits positive P length and positive germline trimming at the same end, and constructs the palindrome from the newly trimmed sequence. A built-in IGK-reference fixture yields `v_trim_3=1`, `j_trim_5=1`, `p_v_3_length=2` and `p_j_5_length=1`.

[PAdditionPass](../engine_rs/src/passes/p_addition.rs), around line 134, deliberately obtains post-trim coding bytes. IMGT identifies P nucleotides at intact coding ends, resulting from asymmetric hairpin opening without subsequent exonuclease removal of those ends. [IMGT P-nucleotide explanation](https://www.imgt.org/textes/FAQ/)

This is a mechanistic modeling issue in addition to the hard P-metadata defects above. Decide whether the simulator should model hairpin opening followed by end processing, or enforce intact-end eligibility for retained P additions. A palindrome synthesized after germline deletion should not silently be labeled the retained product of the original hairpin.

### Missing terminal coverage is encoded as deletion columns

Removing three bases from the 5′ end of an otherwise exact nine-base V gives `v_cigar="3D6M"`, reference interval `[0,9)` and identity 2/3, although the observed bases exactly match `[3,9)`. The symmetric J 3′ case also reproduces.

This behavior is intentional in `builder.rs:574–582`: the walker aligns missing terminal read positions as deletions against the retained germline. An alignment confined to observed coverage would instead start at reference position three and retain identity 1.0.

The discrepancy is confirmed, but terminal deletion columns alone do not prove an invalid alignment, and AIRR allows tool-dependent identity definitions. [AIRR definition clarifications](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#definition-clarifications)

Decide and document whether these fields describe a global alignment to the retained precursor or only observed reference coverage. Downstream users must be able to distinguish missing coverage from structural deletion damage before using CIGAR or identity to quantify biological mutation. The two partial-loss probes describe this contract choice; they are not counted as independently confirmed correctness defects.

## What the passing checks establish

The independent 1,200-record sweep covered human IGH and IGK, 100 records each under baseline, substitutions/N/PCR/quality corruption, indels, end loss, paired-end projection, and a combined forward pipeline. It checked ungapped query order, alignment lengths, segment slices, CIGAR consumption, recomputed identity and paired-read windows. All passed. These checks deliberately do not equate coherent arithmetic with correct biological semantics; terminal deletion columns can be internally consistent yet conflate missing coverage with sequence deletion for downstream consumers.

An additional 200 ordinary records each for human IGH, IGK and IGL had consistent junction slices, lengths and independent translations in the declared frame. Alternative-allele reference checks passed on the two 1,000-record built-in SHM batches described in finding 9.

The full Python run (`pytest tests -q -m 'not performance' --disable-warnings`) passed **2,448 tests**, with eight performance tests deselected, in 483.97 seconds. The existing Rust debug suite passed **1,387 tests**, with two ignored integration tests and two ignored doctests. The full release-mode Rust invocation hit an existing test that unconditionally expects a debug-only panic (`pool_range_after_segment_replacement_overlap_panics_in_debug`); the same test passes under the debug profile. The release metadata/property integration suites passed 104 tests, with two ignored.

Passing existing tests does not refute the findings: several validators recompute the same assumptions, and some tests explicitly encode them. Most diagnostic examples produce no validator issues.

This audit does not certify every bundled species/reference, every genotype or lineage combination, the empirical realism of allele/trim/NP distributions, or full receptor function. The simulated V(D)J sequence alone does not establish promoter, splice-site or constant-region integrity. Historical source alleles and evidence-based `v_call`/`d_call`/`j_call` are also distinct concepts: a call differing from `truth_*_call` is not automatically a bug.

## Recommended repair order

1. Correct whole-coding-region functionality and translation, then update constraints and independently specified expectations.
2. Repair alignment projection for P additions, D orientation and final allele assignments, and settle the observed-coverage convention. These fields underpin downstream mutation and recombination analysis.
3. Preserve founder/revision provenance, then make reverse-strand projection and AIRR export consistent.
4. Make validation inspect actual outputs and check reference bounds, query order, coding phase and schema semantics independently.
5. Resolve the P-end mechanism and same-allele revision contract, then expand the matrix across species, clonal/genotype paths and feature combinations.

Each repair should turn its diagnostic into a positive regression test. Merely changing both the builder and its mirrored validator to agree would not establish biological correctness.
