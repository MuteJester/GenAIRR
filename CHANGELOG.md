# Changelog

All notable changes to GenAIRR are documented here.

## [2.4.0] - 2026-11-05

**Read this first if you have stored data.** `DETERMINISM_VERSION` moved from
the unnumbered model that shipped through 2.3.1 to **4**. A given seed does not
produce the 2.3.1 output any more — not approximately, not statistically. Every
exported record and every trace file now carries the number as
`determinism_version`, and `GenAIRR.DETERMINISM_VERSION` reads it without
running anything. A dataset with no such column was drawn by the pre-2.4.0
model. See [MIGRATING.md](MIGRATING.md) for what to do with data you already
have.

### Breaking and behaviour-changing

- **`DETERMINISM_VERSION` is now 4** — three separate increments, each of which
  moves what a seed produces:

  | version | what changed |
  |---|---|
  | 2 | one random substream per pass and per clonal-family role, keyed rather than offset or XOR-salted; a substitution draws a base the site does not already hold |
  | 3 | clonal lineage regulation is Beverton-Holt; `n_max` is the population a lineage settles at |
  | 4 | affinity selection scores the protein the cell actually expresses |

  Each is described in full below. Version 1 is the unnumbered model that
  shipped up to and including 2.3.1.

- **Clonal lineage regulation obeys a different law, and `n_max` means
  something different.** `n_max` is now the live population a lineage
  **settles at**, not the pivot of a law that settled lower and not a hard cap.
  **The same `n_max` therefore yields about 3x more cells at the default
  `lambda_base=1.5`**, with proportional memory and runtime; lower `n_max` to
  recover the old population, whose plateau was `n_max * (1 - 1/lambda_base)`.
  `lambda_base` is documented as what it is: the Poisson mean for a
  founder-baseline cell in an empty germinal center, not the mean a cell draws
  against once capacity damping or selection is active.

- **`Experiment.expand_clones()` is deprecated** and raises a
  `DeprecationWarning`. It keeps working, and it is still the only surface that
  delivers an exact `n_clones x per_clone` count. Use
  `clonal_lineage()` for BCR affinity-maturation trees with a ground-truth tree
  per clone, or `clonal_repertoire()` for TCR and flat-BCR abundance
  repertoires with clone-size distributions and `duplicate_count`. Neither is a
  drop-in replacement for a fixed `per_clone`: a lineage run returns at most
  `n_sample` records per clone after genotype collapse, and a repertoire run
  collapses identical copies into one record carrying the count.

- **Most bundled catalogues are quarantined under strict compilation.** A
  reproducible reference audit marks every bundled catalogue carrying
  unresolved anchor or pseudogene/ORF findings, and strict compilation — the
  default — refuses them. 80 of the 106 bundled configs are affected, including
  every mouse catalogue and `human_tcrb`; the three human OGRDB cartridges,
  `HUMAN_IGH_EXTENDED` and 22 IMGT catalogues still compile strictly. Add
  `.allow_curatable_refdata()` to sample from a quarantined catalogue anyway,
  or `.curate_refdata(policy)` to drop the offending alleles. The error lists
  the findings allele by allele.

- **`DETERMINISM_VERSION` is now 4.** Affinity selection scores the protein the
  cell actually expresses. Seed-moving wherever selection is active;
  `selection_strength=0` runs with no target are unaffected.

  The cell's amino acids came from concatenating a per-region translation, and a
  codon was kept only if it fitted wholly inside one region. Every region
  boundary in a rearrangement is a junction boundary, so the residues selection
  could not see were exactly the ones it exists to select on: measured over 200
  `human_igh` records, 2.68 codons dropped per record (range 0 to 4) and **536
  of 536 of them straddling V/NP1, NP1/D, D/NP2 or NP2/J**. The frame did not
  restart, so every scored codon was real - the dropped ones simply vanished and
  shifted every amino-acid position after them left. For one record the scored
  protein began `VYYCARGYCSSTSC` where the real one is `VYYCARDRGYCSSTSC*KK`:
  the loss starts at the conserved CDR3 anchor.

  Measured through the public API: a cell grown with `rate=0`, and therefore
  identical to its founder, scored against the founder's own protein reported
  an affinity of **3.57e-67 instead of 1.0**. The model could not recognise a
  cell as identical to itself.

  Both paths now translate the coding span continuously, in the frame resolved
  from the V allele's native reference phase - the same protein a record reports
  as `sequence_aa`. Because the frame comes from the reference, an affinity
  model now travels with its reference data, and `simulate_lineage` requires
  `refdata` when selection is active.

- **A stop codon no longer scores as a perfect match.** `substitution_cost` is
  `S(a,a) - S(a,b)`, and for a stop or an unresolved `X` both terms fell through
  to the same `-4`, so the cost cancelled to exactly zero against every target
  residue and `affinity` came out at 1.0 - the maximum. A `KPGT*` cell reported
  affinity 1.0 against `KPGTA`, `KPGTW` and `KPGTAAAAAAAA` alike. Non-residues
  are now charged the maximal substitution cost, so a receptor that cannot be
  translated is scored as far from the target rather than close to it.

  It is deliberately *not* treated as lethal. SHM here is not constrained to
  keep a receptor productive, so stops accumulate fast -- **96% of sampled
  cells carry one at the default `rate=0.05`**, 89% at `rate=0.01`, and 50% even
  at `rate=0`, inherited from founders, 54% of which carry one already. Giving a
  stop-bearing cell no offspring would therefore extinguish essentially every
  lineage under selection rather than modelling a germinal centre. That the
  mutation model produces mostly untranslatable receptors is a separate matter
  from how affinity scores them.

- **A `target_aa` of the wrong length is refused.** The distance compared
  `min(len, len)` positions and ignored the rest, so a truncated or over-long
  target scored silently and reported a number that was not a distance to it.
  A supplied target must now match the receptor protein's length; the error
  names both.


- **`DETERMINISM_VERSION` is now 3.** Clonal lineage population regulation is
  Beverton-Holt, `lambda_base / (1 + (lambda_base - 1) * P / n_max)`, in place
  of the discrete logistic map `lambda_base * max(0, 1 - P / n_max)` and the
  hard cap that truncated any generation exceeding `n_max`. A seed does not
  carry across the two: every `clonal_lineage()` and duplex family run moves.
  Datasets generated before this are datasets from the old law and carry
  `determinism_version` 2 to say so.

  The old law was wrong in three ways. Reaching capacity was fatal: the cap
  truncated a generation to exactly `n_max`, where `1 - P/n_max` is exactly 0,
  so the whole next generation drew zero offspring and the lineage died with
  certainty. Measured across 180 runs, 67 reached capacity and all 67 went
  extinct in the following generation; over 40 generations at `lambda_base=4.0`
  every cap hit was fatal at every capacity tested. With the default
  `allow_extinction=False` such a family is silently discarded and re-grown
  under a fresh seed, so the delivered data quietly excluded the clones that
  expanded most; when every attempt hit the cap the run raised an error whose
  advice -- raise `lambda_base` -- made the failure more likely, because a
  faster-growing population reaches capacity sooner.

  `n_max` was also not the capacity it is named for. The logistic law settles
  at `n_max * (1 - 1/lambda_base)`, so a requested 1000 delivered 333 cells at
  the default `lambda_base=1.5`. And above `lambda_base=3` the map is past its
  period-doubling threshold, so there was no settled population at all: it
  cycled between 383 and 875 at 3.5 and moved chaotically between 218 and 942
  at 4.0, while the documentation promised a plateau.

  Beverton-Holt is monotone in `P`, so it settles at `n_max` for every
  `lambda_base` and cannot oscillate. It is exactly 1 at `P = n_max`, so a
  population at capacity replaces itself. The hard cap is gone rather than
  repaired: a generation may overshoot capacity -- about 6% in practice -- and
  the next one is regulated back below replacement. Truncating has to choose
  which cells die, and the previous choice let vector position decide, keeping
  the broods of early-indexed parents and dropping every later parent's brood
  entire. **`n_max` now means what it says, so the same `n_max` yields more
  cells than before** -- about 3x at the default `lambda_base=1.5`, with
  proportional memory and runtime. Lower `n_max` to recover the old
  population: the old plateau was `n_max * (1 - 1/lambda_base)`.

### Added

- **Two clonal models beside the deprecated star.**
  `Experiment.clonal_lineage(...)` grows a real BCR genealogy per clone: a
  generation-synchronous birth-death process with per-division S5F
  hypermutation, optional affinity selection, live-cell sampling and
  genotype collapse. Each clone comes back with an AIRR record per observed
  cell -- `clone_id`, `lineage_node_id`, `lineage_parent_id`,
  `lineage_generation`, `lineage_abundance`, `lineage_affinity`,
  `duplicate_count` -- and with the ground-truth tree itself on
  `result.lineage_trees`, exportable as Newick, FASTA or a node table. It is
  BCR-only and refuses a TCR locus, because T cells do not hypermutate.
  `Experiment.clonal_repertoire(...)` covers what a TCR clone actually is: one
  rearrangement proliferated to a heavy-tailed clone size (power-law or
  lognormal, with an `unexpanded_fraction` of never-expanded singletons), read
  through library prep and collapsed to `duplicate_count`. Neither has a
  per-clone size you set directly; see the deprecation note above.

- **Paired-end read layout and FASTQ export.**
  `Experiment.paired_end(r1_length=..., r2_length=..., insert_size=...)`
  projects every record as an Illumina-style read pair: `r1_sequence`,
  `r2_sequence` (already reverse-complemented), `r1_start` / `r1_end` /
  `r2_start` / `r2_end`, `insert_size` and `read_layout`. It is a projection
  and mutates nothing. `result.to_fastq(path, quality=...)` writes single-end
  FASTQ and `result.to_paired_fastq(r1, r2)` writes two synchronised files
  with portable `/1` and `/2` header suffixes; quality strings come from a
  pluggable model (`"illumina"` or `"constant"`). A requested insert longer
  than the surviving molecule leaves the windows empty rather than clipping,
  which `validate_records` reports as `ReadWindowOutOfBounds`.

- **Each NP region says where its bases came from.** `np1_n_tdt`, `np1_n_p`
  and `np1_n_inserted` (and the `np2_*` triple) partition that region's length
  into TdT N additions, templated palindromic P bases and bases an indel
  inserted. An AIRR `np1` string has always combined all three, so a model
  fitted against it as a pure N addition was fitting something else; the
  partition makes the question answerable from the record. The validator
  checks it through `NpPartitionMismatch` and `NpPartitionOriginAmbiguous`.

- **The reading frame and every predicate of `productive`.**
  `coding_frame_start` / `coding_frame_end` say which span the verdict was
  computed over, and `v_anchor_accepted`, `j_anchor_accepted`, `vj_in_frame`,
  `stop_codon` and `j_tail_in_frame` are the five predicates it is the
  conjunction of, each three-valued, beside `productivity_policy`,
  `productivity_scope` and `productivity_reasons`. A record that disagrees with
  a tool can now be localised to the predicate instead of the verdict -- and
  `j_tail_in_frame` is the one a record can fail on its own, which the other
  four cannot show.

- `n_ambiguous_bases`: the positions `ambiguous_base_calls()` actually masked.
  It is not `sequence.count("N")`, which also counts any `N` the reference
  allele carries, and it moves neither error counter -- the documentation used
  to attribute the masking count to `n_quality_errors`, which
  `sequencing_errors()` owns.

- `result.to_genotype_tsv(path)` and `CohortResult.to_genotype_tsv(path)`:
  a run writes the genotype it was given beside the records it produced, one
  row per carried allele, with `subject_id`, `locus`, `segment`, `gene`,
  `allele`, `chromosome`, `novel` and `gene_state`. Until now the planted
  truth and the records left the library by separate doors -- the records from
  the result, the genotype from the `Genotype` object the caller still had to
  be holding -- so a dataset handed to someone else arrived without the thing
  it was meant to be scored against. A gene deleted on a chromosome keeps its
  row with the literal allele `(none)`, because an empty TSV cell cannot be
  told apart from a missing value, and a result with no genotype raises rather
  than writing a zero-row file that looks like success. Whether a given read
  drew a novel allele is a join of `truth_v_call` against this table on
  `(subject_id, allele)`, which the genotype guide now shows; it is
  deliberately not a record column, since the record already carries both
  halves.

  The writer refuses more than an absent genotype. It refuses a permissive
  genotype that still leaves reference genes unspecified, because the engine
  samples those genes over the whole reference and their alleles would reach
  the records but no row of the file: measured on `HUMAN_IGH_OGRDB`, a
  permissive genotype with one gene pinned writes two rows while the records
  draw 119 alleles that are not in it, and a benchmark scoring a caller
  against that file marks correct calls wrong. It also refuses a cohort where
  only some subjects carry a genotype, and two subjects sharing one
  `subject_id`, since both produce a table that looks complete and is not.
  Per-chromosome expression weight (`chromosome_weights`) is *not* in the
  sidecar and is documented as such.

- `determinism_version` on a trace file as well as on every exported record.
  A trace recorded a seed and a plan signature, which together say what to
  draw and from where, but not which model does the drawing: replaying a file
  from another engine produced a different molecule and nothing in the
  artifact disagreed. Trace-file schema version 5 carries the number, and a
  replay whose number disagrees with the installed engine now fails by name.
  Files written before the stamp report `0`, meaning "unknown", and keep the
  treatment they have always had.
- `determinism_version` on every exported record, and
  `GenAIRR.DETERMINISM_VERSION` to read the same number without running
  anything. The engine's reproducibility promise is that one seed and one
  determinism version give one answer, on every supported platform; two
  different versions promise each other nothing at all. The number increments
  whenever a released change moves what a seed produces, and it travels on the
  record rather than in a sidecar because a dataset leaves this library as a
  file of rows and is then passed on by itself. It ships at `4`; a dataset with
  no such column was produced by the unnumbered model that preceded the stamp,
  which is why the numbering does not start at 1. The table at the top of this
  release says what separates each version from the one before it.

  The reproducibility goldens assert the promise in CI on Linux, macOS and
  Windows across Python 3.10-3.13. Linux aarch64 ships a wheel but is not
  covered by them.

- Six per-V-subregion recoverable mutation counters —
  `n_fwr1_observable_mutations`, `n_cdr1_observable_mutations`,
  `n_fwr2_observable_mutations`, `n_cdr2_observable_mutations`,
  `n_fwr3_observable_mutations`, `n_v_unannotated_observable_mutations` — the
  regional counterpart of `n_observable_mutations`. **Use these for a CDR/FWR
  selection-pressure ratio computed from a read.** The existing per-subregion
  counters are a correct partition of `n_v_mutations`, but they count
  hypermutation *events* over the molecule: measured against an independent
  oracle that buckets the visible differences by the same IMGT table, they
  agree with what a read shows on 60 of 60 records at one mutation, 32 of 60
  at twenty, and 17 of 60 once twelve nucleotides are lost from the 5' end,
  where the record counts mutations in bases the read no longer contains. The
  six new fields agree on 60 of 60 in every one of those configurations. They
  sum to the V part of `n_observable_mutations` and each is at most its event
  counterpart; `validate_records` checks both through
  `VSubregionObservableExceedsEvents` and
  `VSubregionObservableSumExceedsTotal`. The committed goldens gain six
  columns with no existing value changed — though every record in them
  carries zero mutations, so that stability evidences nothing about these
  fields.

- `n_observable_mutations`: the somatic-hypermutation substitutions a mutation
  caller can recover from `sequence` — the positions a mutation pass wrote a
  changed base at that are still in the molecule and still differ from the
  germline base recorded there. **Score a mutation caller against this field,
  not against `n_mutations`.** It is below `n_mutations` by three intended
  effects: a mutation on an untemplated junction base changes a base but has
  no germline base to differ from (about 5% of records that carry a mutation),
  two mutations can hit one position with the later one restoring the germline
  base, and a position lost to end loss or a deletion is no longer in the
  molecule — so unlike `mutation_rate`, which is a burden preserved through
  observation loss, this count shrinks with the molecule.

  Attributable to hypermutation and to nothing else: a PCR misincorporation, a
  sequencing miscall, an `N` mask and a background contaminant all leave a
  position differing from its germline, and none of them moves this count.
  Those are reported on `n_pcr_errors` and `n_quality_errors`, which is what a
  caller's false positives should be measured against. Two scoping notes: the
  count is over `sequence`, which stays the whole molecule under
  `paired_end()` while `r1_sequence` / `r2_sequence` cover less; and on a
  contaminant record the molecule is a background sequence, so no mutation
  field on it describes this rearrangement's hypermutation.

  Produced on the native, duplex-observed and duplex-projected planes, and
  checked by `validate_records` through `NObservableMutationsMismatch` and the
  cross-field `ObservableMutationsExceedCount`. It appears in the default TSV
  column order beside `n_mutations`, so the committed goldens gain one column
  and no existing value changes.

- `GenAIRR.Reference`: portable local germline references built from V/D/J
  FASTA and a versioned annotation sidecar (`Reference.annotation_template`,
  `from_fasta`, `report`, `capabilities`, `require`, `write`, `read`). Every
  annotation is unknown or carries its support (`source` with a citation, or
  `assertion`) and, for coding phase, whether it is anchor-derived. Capabilities
  report their evidence level and never upgrade an assertion to source.
  `Experiment.on(reference)` runs offline with the explicit
  `genairr.unfitted_defaults.v1` profile; genotype- and cohort-constrained runs
  and family-replay capture are not yet supported for portable references.
  New errors `ReferenceArtifactError` and `ReferenceCapabilityError`. See the
  "Portable local references" guide.
- Experimental `CoupledDuplexRepair` and `DuplexJoin` scripts for successive
  V-D/D-J joins on one shared two-strand D substrate, including inversion,
  canonical source ancestry, a global nucleotide ledger and cumulative budgets.
- `CodingEndOutcomeModel` with joint D outcomes, compatible-reference identities,
  generating-history estimation, held-out scoring and finite compatible-history
  likelihood fitting. Calibration records its source and sampling assumptions.
- Experimental `GenAIRR.experimental.DuplexRepair` for authored two-strand
  coding-end opening, resection, N extension, compatible annealing, paired-primer
  template synthesis and ligation, with both strand ancestries and exact replay.
- Trace schema 4 records productivity policy and active contracts. Replay and
  seed reruns share reference-content and address-version validation.

- Opt-in `Experiment.process_coding_ends()` with explicit original coding-end
  boundary confirmation and finite signed net outcomes. Positive values remove
  reference bases; negative values retain original-tip P extensions. The mode
  replaces automatic trim/P settings, preserves per-base origin provenance,
  and supports productive completion constraints. It does not claim calibrated
  microscopic hairpin-opening probabilities.

### Changed

- **Every pass now draws from its own random substream.** The engine held one
  random stream for the whole pipeline, so each pass's draws began wherever the
  passes in front of it had left off. Adding a pass, removing one, or
  configuring one to do nothing therefore moved every draw after it: measured
  over 200 human IGH records with a `mutate` behind it,
  `receptor_revision(prob=0.0)` left the rearrangement identical on 200 of 200
  records and the final sequence identical on **0 of 200**. A seeded ablation —
  the same repertoire with and without one step, which is how anyone measures
  what a step does — was therefore a valid comparison about which alleles were
  drawn and an invalid one about the molecule.

  Each pass is now handed a stream derived from the run seed and the pass's own
  identity: its name and which occurrence of that name it is. Two
  `pcr_amplify` steps in one pipeline stay independent, and `mutate` draws the
  same words whether or not an `invert_d` sits in front of it. A pass
  configured off is now indistinguishable from a pass that is not there, the
  molecule included.

  The same derivation now covers a clonal family's three roles — growth,
  sampling and affinity-target selection. Those were separated by salting the
  family's seed with XOR, which separates nothing: XOR is invertible, so a
  salted stream is the plain stream of another seed, and for the two salts in
  use there was an offset (`0x27203f3729310000`) at which one family's
  sampling draws were another family's target draws, word for word. Folding a
  role's name into the seed has no such offset.

  **This changes the output of every seeded run.** Not only runs that use a
  particular pass: the stream every pass draws from has moved. A dataset
  produced before this change cannot be reproduced after it, which is what
  `determinism_version` below exists to make visible.

- **A clonal run's clone, retry and read seeds are derived, not offset.** Each
  sub-simulation of a clonal run took its seed by adding to the run seed: clone
  *c* used `seed + c × 1_000_000`, a survival retry used
  `clone_seed + attempt × 7919`, a read or an observed node used
  `clone_seed + 1 + index`. Addition is a translation, so every run seed shared
  one lattice and a different seed reused streams instead of drawing new ones.
  Clone 1 of a run at seed 0 was the same simulation as clone 0 of a run at
  seed 1,000,000 — identical founder, identical `v_call`/`d_call`/`j_call`,
  identical junction, and with one read per clone the identical corrupted read;
  on the lineage path a whole family repeated, mutations included. **Replicate
  runs at nearby seeds were therefore not independent replicates.** Two further
  collisions came from the same arithmetic: a node's observation seed met a
  retry's seed inside one clone (node 7918 and attempt 1), and nothing bounded
  a clone's node offsets from crossing into the next clone's block.

  Every axis is now a keyed substream of the run seed with the axis as its
  role, through the same primitive the passes use. The run seed is folded in
  rather than translated, `read` and `attempt` cannot meet, and no axis has a
  block to overflow. This changes the output of every seeded clonal run.

  A family-replay bundle records which scheme drew its seeds, and a bundle
  written under the previous one still replays: the recorded label picks the
  scheme instead of being compared against the current one.

- **A cohort subject's seed and a genotype attempt's seed come from the same
  primitive as everything else.** Both were drawn from
  `random.Random(run_seed).getrandbits(63)`. That mixes, so neither had the
  collision described above — subject 1 of one cohort was not subject 0 of
  another. What it did have was a second generator in the derivation chain:
  CPython's Mersenne Twister decided engine seeds, somewhere the engine's
  determinism version does not reach and neither seed-discipline guard could
  see, so "how does GenAIRR turn a seed into a sub-simulation's seed" had two
  answers. `run_cohort(..., seed=S)` now seeds subject *i* with a substream of
  *S* named `subject`, and `Genotype.sample(..., seed=S)` seeds viability
  attempt *a* with one named `genotype_attempt`; a failed attempt still cannot
  land on what a direct draw at `S + 1` would get, which is what the base
  generator was there for. **This changes which genotype a given seed draws
  and which records a given cohort seed produces.** A package-wide guard now
  fails if a seed is assigned from anything but that primitive.

- **`count` is now a count of events, not of attempts, for `mutate(model="uniform")`,
  `pcr_amplify` and `sequencing_errors`.** Each of the three drew its replacement
  base from all four bases including the one already at the position, so about a
  quarter of draws wrote back what was there. Such a write changes nothing and was
  correctly not counted, which left three of the five corruption budgets delivering
  roughly 76% of their nominal value with nothing in the record or the documentation
  saying so: `pcr_amplify(count=8)` produced about six errors. The replacement is now
  drawn from the three bases the position does not hold, so `n_pcr_errors`,
  `n_quality_errors` and `n_mutations` equal the configured count. Measured
  realized/requested on human_igh, 300 records: `pcr_amplify` 0.757 → 1.000,
  `sequencing_errors` 0.757 → 1.000, `mutate(model="uniform")` 0.780 → 1.000.
  `mutate(model="s5f")` and `ambiguous_base_calls` were already at 1.000 and are
  unchanged — the S5F kernel's substitution row carries no mass on the base at the
  site, and `ambiguous_base_calls` pins its write to `N` rather than sampling.

  **This changes the output of every seeded run that uses one of the three.** The
  sequences differ, not only the counts. A dataset produced before this change
  cannot be reproduced after it.

  A requested error can still go undelivered, in one case: under a contract such as
  `productive_only()` a position that admits no replacement is left alone, and the
  counters report what was applied rather than what was asked for.

- Every simulated `SimulationResult` (bundled references included) now records
  the content hash of the reference it ran against. `validate_records` and
  `validate_families_with_parents` against a different reference now raise
  `ValueError` instead of misreporting. A genotype run with novel alleles must
  be validated against `compiled.refdata` (the effective reference).
- Artifact files written by GenAIRR (portable references and family-replay
  bundles) now follow the user's umask (for example 0644) instead of owner-only
  0600. An existing file is still never overwritten without `replace=True`, and
  the destination appears only once complete. On a filesystem without hard links
  (FAT/exFAT, some SMB shares) a write without `replace=True` publishes by rename
  on Windows; elsewhere it fails with an `io` error and leaves the destination
  untouched, because such a filesystem has no atomic no-overwrite step. Write to
  a local disk and copy the file, or pass `replace=True` to overwrite.
- `Experiment.on("<unknown name>")` errors now point to
  `GenAIRR.Reference.read(path)`; a path string is never read implicitly.
  `Experiment.on("x.json")` (any string ending in `.json`) now raises
  `TypeError` where it previously raised `ValueError`; a path object such as
  `Experiment.on(Path("x.json"))` still raises `TypeError`, now with the same
  pointer.

### Fixed

- Lineage and clone metadata is now validated on the legacy `clonal_lineage`
  path. The check that resolves a record back to its tree node ran only in
  duplex mode, so on the legacy path **every one of the eight clonal fields
  could be edited and `validate_records` still returned ok** — `clone_id`,
  `sequence_id`, `lineage_node_id`, `lineage_parent_id`,
  `lineage_generation`, `lineage_affinity`, `duplicate_count` and
  `lineage_abundance`. The immutable (tree, clone, node) origin addressing is
  now collected on both paths, and `sequence_id` is checked against the
  `clone{clone_id}_node{lineage_node_id}` composition it is built from —
  which is what the addressing exists for, since two cells can carry
  identical DNA and history and an edited namespace tag is invisible to every
  other check. `parent_id` and `lineage_mutations_from_parent` remain
  duplex-only and are required there rather than checked wherever they appear.

- `n_pcr_errors` and `n_quality_errors` now count the substitutions each
  corruption pass **applied that changed a base**, not the number it was asked
  to introduce. A replacement is drawn from A/C/G/T including the base already
  present, so some draws write a byte and introduce no error: measured at 500
  reported against 371 realized over 100 records at `count=5`, for each pass.
  This is the same defect `n_mutations` carried, in the counters beside it, and
  the record validator read the drawn count too — so it agreed with the record
  because it reproduced the record's error. Fixed on the native and
  duplex-observed planes; the duplex ledger already tracked the applied counts
  separately, and its attempt counts are retained under
  `duplex_repair.sampled_error_attempts`.

  The documentation attributed `n_quality_errors` to `ambiguous_base_calls()`.
  It is populated by `sequencing_errors()`. `ambiguous_base_calls()` masks
  positions with `N` and reports on neither error counter; its count is
  available from the trace as `corrupt.ns.count`.

- `end_loss_5_length` and `end_loss_3_length` now report the total removed when
  more than one end-loss pass is configured at the same end. Both passes record
  at one trace address (`corrupt.end_loss.5`), and the record builder read it
  with a first-match lookup, so `end_loss_5prime(2).end_loss_5prime(7)` removed
  nine bases and reported two, in either order. The record validator read it the
  same way, so it agreed with the record and certified it. `primer_trim_5prime`
  is an alias for `end_loss_5prime`, so the same applied to any mix of the two.
  A single end-loss pass per end — every packaged preset and all committed
  reference output — is unaffected, so no golden or packaged artifact changes.

- `n_mutations` and the per-segment mutation counters now count only
  substitutions that changed a base. A uniform replacement is drawn from
  A/C/G/T including the base already present, and those draws were counted, so
  `mutate(model="uniform", count=7)` reported 7 where roughly a quarter of the
  draws had changed nothing — 23.3% of draws measured at `count=1`, 24.4% at
  `count=20`, against 25% in theory. **`count=N` is the number of substitution
  attempts; `n_mutations` is how many of them changed something, so for the
  uniform model it is now usually below `N`.** `mutation_rate`,
  `mutations_from_parent` and the per-segment and per-subregion counters follow.

  No single `model="s5f"` write is ever a no-op — that kernel does not draw the
  base already present — so the per-pass S5F counters do not move. One S5F
  number does: a lineage node's count is a net difference against the founder's
  germline, and the kernel excludes the base currently present rather than the
  germline base, so a second hit at a mutated position can draw the germline
  base back. Such a position no longer differs from its germline and is no
  longer counted. Measured on a two-generation, rate 0.05 lineage: five
  restored positions across five nodes, moving their counts from
  `[32, 28, 38, 31, 33]` to `[31, 27, 36, 30, 33]`.

  The draw itself is untouched: the RNG consumption, the recorded trace and
  every simulated sequence are byte-for-byte as before, and all committed
  goldens are unchanged. Traces from earlier versions replay identically. One
  exported channel beyond the counters does move: a mutation pass whose every
  draw was a no-op now has a zero delta, and a zero-delta commit emits no
  `MutationCountChanged` event, so such a pass reports one simulation event
  where it previously reported two. Its `BaseChanged` events are unchanged, so
  every write is still recorded; the duplex edit ledger likewise records every
  write including the no-ops, because both are physical records of what the
  engine did rather than counters.

  One gap remains and is intended: a mutation landing on an untemplated
  junction base changes a base and is counted, but has no germline base to
  differ from, so no mutation caller can recover it — measured at 5.2% of
  uniform and 4.5% of S5F records that carry a mutation. See
  `n_observable_mutations` above for the recoverable count.

- The junction of a simulated record now runs from the first surviving base of
  the V Cys codon to the last surviving base of the J Trp/Phe codon, read from
  each base's recorded germline position. Exonuclease trimming eats these codons
  only from the inside, so their outer limits stay put and N/P bases filling
  trimmed codon positions count as junction. Previously a landmark codon partly
  removed by trimming left the junction undefined (about 0.5% of human IGH and
  1.6% of human IGL records under default trimming), as did an indel deleting a
  codon's outer base, and an indel inside the J codon was ignored (the end was
  the codon's first base plus 3; about 0.5% of polymerase indel events). The
  junction is now undefined only when a whole conserved codon is gone, a codon's
  outer base is trimmed from the segment's outer end, or a kept codon base lies
  outside the read. Affected records gain correct `junction`, `junction_aa`,
  `junction_start`, `junction_end`, `junction_length` and `vj_in_frame` values,
  with the matching frame reason in `productivity_reasons`. Record validation
  derives the same bounds independently. `productive` is unchanged (these
  records already lacked an intact landmark), and generated sequences are
  unchanged.
- The DataConfig bridge no longer assigns IMGT FWR/CDR subregions to V alleles
  whose `gapped_seq` declares no IMGT gaps (e.g. ungapped FASTA). Explicit
  subregions and IMGT-gapped sequences are unchanged; all 106 bundled cartridges
  keep identical native content hashes. Older user DataConfigs built from
  ungapped V FASTA now compile without fabricated subregions, so their native
  content hash changes. `compute_v_region_boundaries()` raises `ValueError` for
  ungapped input.
- `ReferenceCartridgeBuilder.infer_identity()` fills per-allele source, species
  and locus only where the loader left them unset, preserving FASTA-header
  provenance, and reports preserved species labels. New `with_build_date()`
  fixes `metadata.last_updated` for reproducible checksums with or without
  `infer_identity()`.
- The native reference bridge rejects V/D/J sequences longer than 65,535 nt
  with `ValueError` before pool mutation, instead of narrowing lengths to `u16`.
- `ReferenceCartridgeBuilder.from_fasta()` decides whether a V sequence is an
  IMGT alignment from its gap markers, not its header. Header-less IMGT-gapped
  FASTA (e.g. OGRDB `gapped`) keeps its validated Cys104 anchor instead of
  losing it. An IMGT-header ungapped V no longer reads codon 309 as IMGT residue
  104 and fabricates an anchor; its anchor stays unknown.
- Reference authoring uses the shared native landmark resolver, retains FASTA
  functionality/source declarations, and keeps missing species unknown.
  Optional checksum-bound evidence transports independent coding phase and
  coding-end confidence without changing legacy reference hashes. A reproducible
  audit marks bundled catalogues with findings as quarantined for default strict
  compilation; unresolved source evidence is explicit and archived data is not
  silently corrected.
- Replace approximate constrained-indel continuation weights with versioned
  natural sequential history sampling and exact weighted fallback. Report work
  exhaustion separately from proven empty support; validate replay count, kind,
  base weights and checked coordinates. Preserve sampling diagnostics separately
  from biological events.
- Use complete productive continuation support for ordinary zero-P VDJ
  recombination as well as P-aware processing, excluding future D/junction stops
  before choosing an allele or length.

- Use versioned observed-span productivity with nullable evidence and stable
  diagnostics. Reference admission overrides no longer legitimize invalid
  landmarks or frames; ambiguous possible stops remain unknown. Default
  anchors require V Cys and J Phe/Trp classes with original provenance.
  Productive-only runs assert accepted final metadata in both modes, and the
  record validator derives an independent assessment. Lineage divisions
  inherit the same constraints instead of silently mutating outside them.
- Filter ordinary and genotype-aware receptor-revision candidates by active
  biological contracts before drawing, preserving their original weights.
  Permissive execution and replay now enforce the same replacement checks.

- Keep standard V/D/J calls and alignment references tied to the generating
  allele; retain post-mutation inference separately as `recalled_*_call`.
  Ground-truth NP boundaries no longer move when random junction bases match
  a neighboring trimmed allele.
- Represent end loss as missing reference coverage, preserving actual
  structural deletions inside the observed interval.

- Preserve mutation-stage molecule lengths in `mutation_rate`, so later read
  loss does not inflate the reported cumulative SHM event burden. Clonal
  descendants retain inherited mutation counts and their event provenance.
- Translate the physical nucleotide junction even when it is out of frame
  or nonproductive; complete codons remain in `junction_aa`.
- Report impossible explicit recombination trims as structured errors instead
  of panicking, including in permissive runs.
- Load every bundled S5F model without pandas or NumPy. The two migrated
  artifacts retain exactly the same mutation-model values.
- Describe D inversion and its configured probability in `Experiment.describe()`.
- Calculate productivity and translation from the retained V coding origin,
  including internal stops and structural frameshifts outside the junction.
- Preserve coding-oriented annotations when emitting reverse-strand observations.
- Preserve P-nucleotide order in alignments and include P bases in combined NP
  strings; align inverted D references in their assembled orientation.
- Keep germline alignments within the called allele, including when its length
  differs from the sampled allele.
- Preserve recombination trims and P-addition provenance through clonal forks
  and receptor revision. Revision selects a different eligible V allele.
- Keep realized loss counters when a read is removed completely.
- Support nonzero P additions with productive-only sampling by retaining
  choices that admit a complete stop-free, in-frame sequence. Fixed P
  lengths are preserved; impossible support fails explicitly.
- Refresh live-call coordinates after no-op and boundary-edit passes so later
  end loss cannot reuse stale D/J positions.
- Use the correct canonical trim-end budgets for cached inverted-D re-calls,
  keeping NP-extension bounds and allele ties consistent with fresh scoring.
- Validate returned record fields and record/outcome counts; prevent custom
  metadata from overwriting generated sequence annotations.

### Documentation

- **Documented that a seed names a record, not a run, on the plain path.** `n`
  records are drawn at `seed` .. `seed + n - 1`, which is what makes a dataset
  splittable and resumable, and the guides already said to offset the starting
  seed to stitch two batches. What none of them said is what happens if you
  offset by less: the batches overlap. Measured on human IGH at `n=200`,
  `seed=0` and `seed=1` share 199 of 199 comparable records, while `seed=0` and
  `seed=200` share none. So runs at seeds 0, 1, 2 are *not* three independent
  replicates of an experiment — they are one dataset read from three starting
  points, and variance measured across them is almost entirely shared. The
  behaviour is unchanged and deliberate; `run()`'s docstring,
  [the builder guide](site_docs/guides/experiment-builder.md),
  [the validation summary](site_docs/validation/index.md) and a new section of
  [the trace/replay guide](site_docs/guides/trace-replay.md) now say to step the
  seed by at least `n` for replicates. A clonal run's own sub-simulations do not
  work this way and never did after the change above: each takes a keyed
  substream, so adjacent clonal seeds share nothing.

- Every code example in the README and on the documentation site was re-run
  against this release and the values it prints were re-measured, so the
  allele names, counts and field values on the page are the ones the engine
  produces at `determinism_version` 4 rather than the ones the pre-stamp model
  produced. The reproducibility promise now states what is actually checked
  (same seed and same determinism version, on Linux, macOS and Windows across
  Python 3.10-3.13), the quarantined-catalogue set is documented with the
  `allow_curatable_refdata()` opt-in, and the pages that said the per-region
  mutation counters, the cartridge model estimators or the NP origin fields
  were deferred now describe them, because they ship.

- `MIGRATING.md` gained a 2.3.1 to 2.4.0 section; the 1.x to 2.0 guide is
  unchanged beneath it. A stale find/replace row that told 1.x users to rename
  `expand_clones(n_clones=...)` to `n=...` was removed -- the keyword is and
  always was `n_clones`, and the rename it described would raise `TypeError`.

### Compatibility

- Standard `v_call`, `d_call`, and `j_call` now name single generating alleles.
  Consumers of the previous evidence-based calls should read `recalled_v_call`,
  `recalled_d_call`, and `recalled_j_call`; optional `truth_*_call` columns remain
  aliases of the standard calls. Alignments now consistently use the generating
  reference, including when the inferred allele differs in length.
- Donor-genotype re-calls consider only alleles carried across that donor's
  two chromosomes, while retaining their original reference-pool IDs.

- `mutation_rate` sums each mutation stage's realized events divided by that
  stage's input molecule length. It is an event burden, not the mismatch
  fraction in an observed read; lineage records retain net-difference semantics.
- Clonal descendant outcomes now include ancestral trace, revisions, and events
  once, so inherited mutation counters have complete provenance.
- The raw-file checksums of `hh_s5f_60` and `hh_s5f_opposite` changed when their
  pandas objects were replaced by portable dictionaries; kernel values are
  bit-for-bit unchanged.
- `to_tsv()` now defaults to AIRR coordinates, CIGAR placement offsets and
  `T`/`F` boolean encoding. Pass `airr_strict=False` for the previous zero-based
  TSV representation. Python record dictionaries and default DataFrames keep
  zero-based half-open coordinates and segment-local CIGAR operation runs.
- AIRR `np1`/`np2` strings and lengths now include P bases between resolved
  segment alignments; they must not be interpreted as pure TdT N additions.
- Reverse-strand `*_sequence_start/end`, junction and amino-acid annotations
  now use the coding orientation, as AIRR requires. Read-window coordinates
  continue to index the emitted observation.
- Productive-only sampling excludes retained coding stops and protects the
  complete observed coding sequence. Seeded productive results can therefore
  differ from previous versions.
- Productive P-addition sampling restricts each draw's original weights to
  candidates with a productive completion across deterministic P bytes,
  N-addition support, and D orientation. Impossible support raises an error
  even in permissive mode instead of substituting zero P or sentinel N bases.
- NP length/base estimators reject records with explicitly nonzero P-length
  provenance. AIRR NP fields combine N and P; fitting them as N-only models
  while separately fitting P lengths would count P additions twice.
- `pandas` is not part of any extra and never was. `result.to_dataframe()`
  needs it; install it yourself. The install docs and the README said it
  shipped under `[all]`, which installs numpy, scipy, graphviz, tqdm and
  fastmcp.
- The pytest matrix runs on Python 3.10, 3.11, 3.12 and 3.13 across Linux,
  macOS and Windows. Python 3.9 remains supported and its wheels are
  smoke-tested on all three platforms, but the suite is not run against it.

## [2.3.1] - 2026-07-05

Maintenance release. No public API or simulation-behaviour changes.

### Added
- `tools/build_imgt_configs.py` — a first-class maintainer tool that builds
  structural cartridges from IMGT V-QUEST germline FASTA via
  `ReferenceCartridgeBuilder`.
- `__all__` on the documented public authoring modules (`reference_models`,
  `reference_rules`, `cartridge_builder`, `genotype_priors`, `cohort`).

### Changed
- Large behaviour-preserving internal reorganization: several god-files split
  into focused packages/modules (`experiment.py`, `_cartridge_estimators`,
  `_compiled`, `result`, `genotype`, `utilities/visualize`, and the Rust
  `airr_record/validate.rs`), plus dead-code and duplicate-module removal. No
  user-visible import paths changed.

### Documentation
- Clarify that only the human IGH/IGK/IGL and TCRB cartridges carry
  data-derived distributions; every other bundled cartridge uses
  uniform/placeholder parameters (with guidance on fitting real ones).

## [2.3.0] - 2026-06-18

### Added - Per-individual diploid genotypes

A new `Genotype` layer makes V(D)J recombination haplotype-phased: attach a
genotype to an `Experiment` and the V, D and J of each rearrangement are drawn
from a single chromosome, honouring allele presence/absence, zygosity, gene
deletion, copy number, and novel/private alleles.

- **`Genotype` builder** - `from_dataconfig` / `permissive`, `homozygous` /
  `heterozygous` / `delete_gene` / `duplicate_gene`, `add_novel_allele` (private
  alleles injected into an effective reference), `complete_from_reference`,
  `chromosome_weights`, `with_subject`, `to_table` / `to_tsv`.
- **`Experiment.with_genotype(g)`** - phased recombination, with per-record
  `subject_id` / `haplotype` provenance and `result.genotypes`.
- **`Genotype.sample(...)`** - draw a genotype from population priors
  (independent per-gene Hardy-Weinberg: allele frequencies + per-haplotype
  deletion).
- **Cartridge genotype-priors plane** - persist a donor-population germline
  prior on a cartridge (`DataConfig.genotype_priors` / `PopulationGenotypeModel`);
  author or estimate it (`PopulationGenotypeModel.from_genotypes`, builder
  `set_genotype_priors` / `estimate_genotype_priors`); `Genotype.sample(cfg)`
  consumes it with per-input provenance.
- **`Experiment.run_cohort(...)`** - simulate many subjects (each with its own
  genotype) in one call, returning a `CohortResult` with per-subject access and
  combined export.
- **Genotype-aware receptor revision** -
  `receptor_revision(..., same_haplotype=True)` restricts the replacement V to
  the carried alleles on the drawn rearrangement chromosome.

### Documentation

- New genotype guide cluster: overview, sampling & population priors, cohorts,
  and benchmarking genotype inference (with a worked TIgGER / IgDiscover example
  recovering a planted genotype).

### Compatibility

- Purely additive. An experiment with no genotype attached produces
  byte-identical output to 2.2.0 (pinned by a checksum test).

## [2.0.0] — 2026-05-25

### Breaking Changes — Experiment DSL refactor for biological readability

The Experiment DSL was reworked around a single principle: a fluent chain
should read like a wet-lab protocol, not a configuration file. The pre-v2
DSL leaked engine mechanics (`corrupt_*` as an umbrella verb, `respect=`
as a runtime kwarg, `count=` as the unit for somatic hypermutation) and
hid biological assumptions inside `recombine(trim=True)`. v2.0 fixes both
ends. **Old scripts will not run unmodified** — see
[MIGRATING.md](MIGRATING.md) for the full map.

#### New methods

- **`.productive_only()`** — declares the productive-sequence contract
  bundle as a chainable method instead of a `respect=ga.productive()`
  kwarg. Order-independent in the chain.
- **`.restrict_alleles(v=..., d=..., j=...)`** — was `.using(...)`.
  Renamed for semantic precision (narrows sampling support, doesn't
  pin one allele).
- **`.trim(enabled=False)` / `.trim(v_3=..., d_5=..., d_3=..., j_5=...)`** —
  exonuclease-trim override method. Default-on inside `.recombine()`;
  call `.trim()` only to disable or supply custom distributions.
  Chain-position-guarded (raises if called after a mutation/corruption
  step). Was previously hidden behind `.recombine(trim=...)`.
- **`mutate(rate=0.03)`** — per-base SHM rate; engine samples count
  per record from `Poisson(rate × pool_len)`. The canonical
  biology default. `mutate(count=...)` remains as the explicit
  benchmark-friendly form; passing both raises.
- **`exp.describe()`** — renders the chain as a biology-style narrative
  per-step. Use it as a quick sanity check before `compile()`.

#### Renamed methods (full map)

| v1 | v2 |
|---|---|
| `respect=ga.productive()` (kwarg) | `.productive_only()` (method) |
| `.using(v=...)` | `.restrict_alleles(v=...)` |
| `.recombine(trim=True/False)` | `.recombine()` + `.trim(enabled=...)` |
| `mutate(count=15)` | `mutate(rate=0.03)` *(canonical)* or `mutate(count=15)` *(retained)* |
| `corrupt_pcr` | `pcr_amplify` |
| `corrupt_quality` | `sequencing_errors` |
| `corrupt_contaminants` | `contaminate` |
| `corrupt_indels` | `polymerase_indels` |
| `corrupt_ns` | `ambiguous_base_calls` |
| `corrupt_reverse_complement` | `random_strand_orientation` |
| `corrupt_5prime_loss` | `primer_trim_5prime` |
| `corrupt_3prime_loss` | `primer_trim_3prime` |
| `with_clonal_structure(n_clones=N, size=K)` | `expand_clones(n_clones=N, per_clone=K)` |

Both `pcr_amplify` and `sequencing_errors` accept `rate=` (per-base
error probability, runtime Poisson) in addition to the legacy
`count=`. Exactly one must be provided. `mutate(rate=...)` shares
the same shape.

#### Stricter validation

- `recombine(np2_lengths=...)` on a VJ chain is now a `ValueError`
  (was silent ignore).
- Binding a raw `RefDataConfig` (no `DataConfig` backing) emits a
  `UserWarning` when synthetic NP-length defaults or no-op trim are
  injected. Pass explicit `np1_lengths=` / `np2_lengths=` and call
  `.trim(enabled=False)` to silence.

#### Engine

- New `MutationCountSource` enum in `engine_rs/src/passes/mutate/`
  with `Distribution` and `Rate(f64)` variants. The rate variant
  samples count via a Knuth Poisson sampler at execute time
  against the current pool length. `S5FMutationPass::new_rate(...)`
  and `UniformMutationPass::new_rate(...)` constructors expose the
  new mode; the existing count-based constructors are unchanged.

#### Tests + benchmarks

- 631 pytest, 637 cargo lib, 128 cargo integration/property — all green.
- Throughput unchanged on the count path (no engine-side cost when
  rate-mode is unused). Rate-mode adds ~1 Poisson sample per record.

---

## [1.0.0] — 2025-03-10

### Breaking Changes

This is a complete rewrite. The Python simulation engine, pipeline system, and
public API have all been replaced. Code written for 0.6.x will not work without
changes.

- **New Experiment DSL** — Single entry point replaces the old
  `AugmentationPipeline` + step-based system. All simulation is now configured
  through `Experiment.on("config").mutate(...).observe(...).run(n=1000)`.
- **C simulation engine** — All rearrangement, mutation, corruption, and AIRR
  serialization now run in compiled C. Typical throughput is 50–100k
  sequences/second on a single core.
- **Removed packages** — `simulation/`, `sequence/`, `mutation/`, `pipeline/`,
  `steps/`, `validation/`, `container/`, `unbiased/`, `TCR/` have all been
  removed. Their functionality is now in the C backend.
- **No runtime dependencies** — The core package has zero mandatory
  dependencies. Optional extras (`numpy`, `scipy`, `graphviz`, `fastmcp`) are
  available via `pip install GenAIRR[all]`.

### Added

- 106+ built-in DataConfigs covering 23 species (human, mouse, rat, rabbit,
  dog, cat, cow, sheep, pig, horse, gorilla, rhesus, cynomolgus, ferret,
  chicken, salmon, trout, zebrafish, platypus, alpaca, dromedary, goat).
- GDC binary format for fast DataConfig serialization to/from the C engine.
- Selection pressure simulation (CDR/FWR replacement acceptance rates).
- Class switch recombination (CSR) with isotype-specific rates.
- D-gene inversion and receptor revision simulation.
- Pipeline introspection hooks — snapshot the internal sequence state at any
  pipeline stage for debugging and analysis.
- Execution tracing — human-readable log of every internal decision.
- MCP server with 21 tools for AI-assisted sequence analysis.
- Cross-platform support: Linux, macOS, and Windows (pre-built wheels).
- Allele locking — constrain simulation to specific V/D/J alleles.
- Streaming interface via `CompiledSimulator.stream()`.

### Fixed

- S5F silent mutations (substitution could draw the same base).
- TCRB anchorless alleles (anchor=0 treated as valid position).
- Indel mutation string truncation (null bytes in annotation strings).
- Junction bounds with 3' corruption (junction_end could exceed sequence length).

## [0.6.3] — 2024-12-01

Last release of the Python-only engine. See the
[0.6.x documentation](https://genairr.readthedocs.io/) for details.
