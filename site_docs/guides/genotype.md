# Genotypes: per-individual diploid germline

<p class="lead">A <strong>genotype</strong> in GenAIRR is one person's diploid germline
complement - which V/D/J alleles they carry, on which chromosome, in what copy
number. Attach a <code>Genotype</code> to an <code>Experiment</code> and V(D)J
recombination becomes <strong>haplotype-phased</strong>: the V, D and J of each
rearrangement are drawn from a single chromosome, honouring allele
presence/absence, zygosity, and gene deletion. With no genotype attached the
engine is byte-for-byte unchanged. This page explains exactly what a genotype is,
how the engine samples from it, and how to build one - nothing here is a black
box. Three companion pages cover the rest:
<a href="genotype-priors/">sampling &amp; population priors</a>,
<a href="genotype-cohorts/">cohorts</a> (many subjects at once), and
<a href="genotype-benchmarking/">benchmarking genotype inference</a>.</p>

## What a genotype is (and why it matters)

Every person inherits two copies of the immunoglobulin heavy-chain locus - one on
each homologous chromosome (one **haplotype** from each parent). Across the
population the locus is extraordinarily polymorphic: a reference set may list
dozens of alleles per gene, but a *single individual* carries only a handful -
typically **one or two alleles per gene** - and may be **missing entire genes**
(deletion) or carry **extra copies** (duplication). That per-individual set is the
**genotype**.

GenAIRR models four things a genotype encodes:

| Concept | Meaning | In GenAIRR |
|---|---|---|
| **Allele presence/absence** | only carried alleles can rearrange | alleles not in the genotype are never sampled |
| **Diploid zygosity** | per gene: 1 allele (homozygous) or 2 (heterozygous) | `homozygous` / `heterozygous` |
| **Gene deletion / copy number** | a gene can be absent on one or both chromosomes, or duplicated | `delete_gene` / `duplicate_gene` |
| **Haplotype phasing** | V, D, J of one rearrangement come from one chromosome | drawn automatically, recorded per record |

Phasing is what makes a genotype more than "a list of alleles to allow". The IGH
locus is physically on a chromosome, so a single recombination event splices a V,
a D and a J **from the same chromosome**. That linkage is exactly the signal
haplotype-inference methods exploit (e.g. the IGHJ6-anchor approach), and GenAIRR
reproduces it.

!!! note "Supported loci and chains"
    Genotypes work on any GenAIRR **reference cartridge** - a packaged germline
    reference set (V/D/J alleles plus the empirical models for one locus; see
    [Reference cartridge](../concepts/reference-cartridge.md)). They cover BCR
    **and** TCR, heavy **and** light/α/β chains. On **VDJ** loci (IGH, TRB, TRD)
    the genotype spans V, D and J and each rearrangement draws all three from one
    chromosome. On **VJ** loci (IGK, IGL, TRA, TRG) there is no D segment:
    genotype V and J, D rows are simply not required and are ignored. The
    examples below use the human IGH cartridge, but the same API applies to every
    locus; just use that cartridge's gene/allele names.

!!! note "A note on the gene/allele names in these examples"
    The examples use the bundled human IGH cartridge (`HUMAN_IGH_OGRDB`), whose
    `metadata.reference_set` reads `"OGRDB V8"` - it is derived from
    [OGRDB](https://ogrdb.airr-community.org/), the AIRR Community's *Open
    Germline Receptor Database*, which provides curated immunoglobulin and T-cell
    receptor germline gene sequences. In this cartridge the genes and alleles are
    labelled like **`IGHVF1-G1*01`**, not the IMGT positional names
    (**`IGHV1-69*01`**) you may be used to - that is simply this cartridge's
    labelling. GenAIRR uses each cartridge's own gene/allele names **verbatim**:
    they are exactly what appears in `v_call` / `truth_v_call` and in the
    genotype tables, so they line up with the germline database you score
    against. Substitute your own cartridge's names in the examples - many bundled
    cartridges (for example `HUMAN_TCRB_IMGT` and the per-species `*_IMGT` sets)
    use IMGT naming instead.

## Quick start

```python
import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype

cfg = gdata.HUMAN_IGH_OGRDB

# Build a diploid genotype: start from the reference, then edit specific genes.
g = (
    Genotype.from_dataconfig(cfg)
      .complete_from_reference("homozygous_first_reference")  # fill the rest
      .heterozygous("IGHVF1-G1", "IGHVF1-G1*01", "IGHVF1-G1*02")  # two alleles
      .homozygous("IGHVF2-G4", "IGHVF2-G4*01")                    # one allele
      .delete_gene("IGHVF3-G7", haplotype="both")                 # gene absent
      .with_subject("DONOR01")
)

result = (
    ga.Experiment.on(cfg)
      .with_genotype(g)        # recombination is now haplotype-phased
      .recombine()
      .run_records(n=1000, seed=7)
)

result[0]["subject_id"]          # 'DONOR01'   - provenance on every record
result[0]["haplotype"]           # 0 or 1      - which chromosome this read used
result.genotypes[0].to_table()   # ground-truth genotype (per gene, per haplotype)
```

## How recombination samples from a genotype

When a genotype is attached, recombination runs a single phased sampling pass per
rearrangement. The steps, in order:

1. **Draw a chromosome.** One of the two haplotypes is chosen, weighted by the
   `chromosome_weights` (default `[0.5, 0.5]`). This choice is made **once** and
   shared by V, D and J - that is the phasing.
2. **Per segment, draw a gene then an allele.** Among the genes *present on the
   chosen chromosome*, a gene is sampled (weighted by usage - see below), then the
   allele follows from that chromosome's slot for the gene. A deleted gene is
   simply not offered on the chromosome that lacks it.
3. **Assign and continue.** The chosen V/D/J alleles are assigned and the rest of
   the pipeline (trimming, NP, assembly, SHM, corruption) runs unchanged.

Every random choice (chromosome, gene, within-slot allele) is recorded to the
trace, so seeded runs are byte-stable and fully replayable.

### Viability and `productive_only`

A chromosome is only drawn if it is **viable** - it must carry at least one usable
allele for every required segment (V and J, plus D on heavy chains). This matters
with deletions: if one haplotype lacks a J gene entirely, only the other
chromosome is ever drawn. If **neither** chromosome can produce a rearrangement,
the genotype is rejected at compile time with a clear error (rather than failing
at run time).

Under [`productive_only`](productive.md), viability also accounts for
productive-junction feasibility, and the V chosen earlier in the pass constrains
the J drawn later (the phased choices are evaluated together, not independently).

### Strict vs permissive

`Genotype.from_dataconfig(cfg)` is **strict**: any gene that could be used during
recombination but was never specified is an error when the experiment is compiled
(`compile()` / `run_records()`) - you must define
the whole genotype (use `complete_from_reference` to fill the genes you don't care
about). This guarantees a genuine diploid complement, which is what you want for a
ground-truth benchmark.

`Genotype.permissive(cfg)` is a separate, explicitly **non-diploid** fallback:
unspecified genes are left to sample over *all* their reference alleles without
phasing. It exists for the "I only want to constrain a few genes" case - it is
**not** a biological genotype, and it is labelled as such in `to_table()` and
`repr`.

In both modes, feasibility (e.g. `productive_only`) is applied the same way the
non-genotype path applies it: candidates are filtered to the feasible set, and the
unfiltered set is used only as a last resort when nothing is feasible - so a
genotype run never silently samples alleles a normal run would have avoided.

## Building a genotype

The `Genotype` builder is a fluent, validated editor over a `DataConfig`'s
reference alleles. Every method returns `self`, so calls chain.

```python
g = Genotype.from_dataconfig(cfg)                       # strict (recommended)

g.homozygous("IGHVF2-G4", "IGHVF2-G4*01")               # 1 allele on both chromosomes
g.heterozygous("IGHVF1-G1", "IGHVF1-G1*01", "IGHVF1-G1*02")  # different allele per chromosome
g.delete_gene("IGHVF3-G7", haplotype="both")            # gene absent entirely (homozygous deletion)
g.homozygous("IGHVF3-G8", "IGHVF3-G8*01")               # carried on both chromosomes...
g.delete_gene("IGHVF3-G8", haplotype=1)                 # ...then removed on chromosome 1 (hemizygous)
g.duplicate_gene("IGHVF1-G2", ["IGHVF1-G2*01", "IGHVF1-G2*02"], haplotype=0)  # >1 copy on one chromosome
g.chromosome_weights(0.6, 0.4)                          # allelic-expression imbalance
g.with_subject("DONOR01")                               # provenance label

g.complete_from_reference("homozygous_first_reference") # fill every unspecified gene
```

Every editing method takes a `segment` argument (`"V"` default, or `"D"` / `"J"`),
so genotype the D and J loci too - important since J anchors and D/J usage drive
haplotype-inference methods:

| Method | Signature | Notes |
|---|---|---|
| `homozygous` | `(gene, allele, segment="V")` | one allele on both chromosomes |
| `heterozygous` | `(gene, allele0, allele1, segment="V")` | one allele per chromosome |
| `delete_gene` | `(gene, haplotype="both"\|0\|1, segment="V")` | whole-gene or one-chromosome (hemizygous) deletion |
| `duplicate_gene` | `(gene, alleles=[...], haplotype=0\|1, segment="V")` | >1 copy on one chromosome |
| `add_novel_allele` | `(name, *, base, mutations\|sequence, segment="V", allow_nonfunctional=False)` | define a private allele (see below) |
| `chromosome_weights` | `(w0, w1)` | allelic-expression imbalance (default 0.5/0.5) |
| `with_subject` | `(sid)` | provenance label stamped on every record |
| `complete_from_reference` | `(policy="homozygous_first_reference"\|"heterozygous_first_two")` | fill unspecified genes |

```python
# Genotype the J locus too - e.g. heterozygous IGHJ6 + a homozygous IGHJ4:
g.heterozygous("IGHJ6", "IGHJ6*02", "IGHJ6*03", segment="J")
g.homozygous("IGHJ4", "IGHJ4*02", segment="J")
```

Notes and guard-rails:

- **`delete_gene(..., haplotype=0|1)`** (one chromosome) requires the gene to be
  specified first - deleting a single haplotype of an *unspecified* gene would
  silently delete both, so it raises instead.
- **`complete_from_reference(policy=...)`** fills only genes you haven't touched.
  `"homozygous_first_reference"` (default) makes each unspecified gene homozygous
  for its first cartridge allele - note this is the *first listed* allele, **not**
  a population-frequency-common one (this policy consults no frequency prior; for
  frequency-driven genotypes use
  [`Genotype.sample`](genotype-priors.md)). `"heterozygous_first_two"` uses
  the first two alleles.
- Unknown gene/allele names, NaN/inf chromosome weights, and segments left with no
  usable allele are all rejected at build/attach with clear messages.
- `with_genotype` is **mutually exclusive** with
  [`restrict_alleles`](../reference/experiment.md) and the
  `recombine(*_allele_weights=...)` kwargs - the genotype owns allele presence and
  expression. `receptor_revision` **is** supported (see
  [Receptor revision with a genotype](#receptor-revision-with-a-genotype)); the
  clonal forks (`expand_clones` / `clonal_lineage` / `clonal_repertoire`) are still
  rejected with a genotype in this release (see
  [Limitations](#limitations-this-release)).

### Gene usage

Within a chromosome, which *gene* is used is weighted by the cartridge's typed
allele-usage model (`reference_models.allele_usage`), aggregated to the gene
level, and scaled by **copy-number dosage** (a duplicated gene recombines
proportionally more often). Cartridges that don't author a typed `allele_usage`
fall back to uniform-over-present-genes (× dosage). See
[Allele usage](v-usage.md) and [Estimate models from data](estimate-cartridge-models.md)
for authoring usage.

### More genotype recipes

**A richer diploid genotype** - several heterozygous genes, a homozygous gene,
a whole-gene (homozygous) deletion, a hemizygous deletion, and allelic-expression
imbalance, with everything else filled from the reference:

```python
g = (
    Genotype.from_dataconfig(cfg)
      .heterozygous("IGHVF1-G1", "IGHVF1-G1*01", "IGHVF1-G1*02")
      .heterozygous("IGHVF1-G2", "IGHVF1-G2*01", "IGHVF1-G2*02")
      .homozygous("IGHVF2-G4", "IGHVF2-G4*01")
      .delete_gene("IGHVF3-G7", haplotype="both")     # absent on both chromosomes
      .homozygous("IGHVF3-G8", "IGHVF3-G8*01")        # carried on both...
      .delete_gene("IGHVF3-G8", haplotype=1)          # ...then removed on chr 1 (hemizygous)
      .chromosome_weights(0.65, 0.35)                 # chromosome 0 expressed more
      .complete_from_reference()                      # the remaining genes
      .with_subject("DONOR_A")
)
```

**Gene duplication** - one chromosome carries two alleles of the same gene
(specify the gene on both chromosomes first, then add the extra copy to one):

```python
g = (
    Genotype.from_dataconfig(cfg)
      .homozygous("IGHVF1-G3", "IGHVF1-G3*01")                          # both chromosomes carry *01
      .duplicate_gene("IGHVF1-G3", ["IGHVF1-G3*01", "IGHVF1-G3*02"], haplotype=0)  # chr 0 now carries two copies
      .complete_from_reference()
      .with_subject("DONOR_DUP")
)
# chromosome 0 carries {*01, *02}, chromosome 1 carries {*01};
# the extra copy raises this gene's recombination share (copy-number dosage).
```

**Build a fully-specified strict genotype programmatically** - drive the builder
from a per-gene plan (the natural shape if you load a genotype from a table or
generate many subjects):

```python
plan = {
    "IGHVF1-G1": ("IGHVF1-G1*01", "IGHVF1-G1*02"),  # 2 alleles  -> heterozygous
    "IGHVF1-G2": ("IGHVF1-G2*01",),                 # 1 allele   -> homozygous
    "IGHVF3-G7": (),                                 # 0 alleles  -> deleted
    # ... one entry per gene you want to pin
}

g = Genotype.from_dataconfig(cfg)
for gene, alleles in plan.items():
    if not alleles:
        g.delete_gene(gene, haplotype="both")
    elif len(alleles) == 1:
        g.homozygous(gene, alleles[0])
    else:
        g.heterozygous(gene, alleles[0], alleles[1])
g.complete_from_reference().with_subject("DONOR_B")
```

**Inspect the non-trivial genes** of any genotype:

```python
for row in g.to_table():
    if row["zygosity"] != "homozygous":
        print(row["gene"], row["zygosity"], row["haplotype_0"], row["haplotype_1"])
```

## Receptor revision with a genotype

[Receptor revision](../reference/experiment.md) models a post-recombination V
replacement. With a genotype attached, the replacement V is drawn from the
**carried** V alleles on the **drawn rearrangement chromosome** (the haplotype the
original V came from), excluding the current V - so the revised receptor stays
consistent with the individual's germline:

```python
g = Genotype.sample(cfg, seed=0, subject_id="donor")
res = (ga.Experiment.on(cfg).with_genotype(g)
       .recombine().receptor_revision(prob=0.2)        # same_haplotype=True by default
       .run_records(n=500, seed=1, expose_provenance=True))
# revised records: original_v_call = pre-revision V; v_call / truth_v_call = the
# carried replacement; receptor_revision_applied = True
```

`same_haplotype=False` is a **synthetic control** that draws the replacement from
either chromosome's carried V alleles - useful for ablation studies, but not a
realistic model of secondary V rearrangement (which is a *cis*, same-chromosome
event). Either way the record's `haplotype` provenance keeps naming the original
rearrangement chromosome.

This is **haplotype-aware V replacement**: it guarantees the replacement is an
allele the individual carries, but it does not model genomic V order, RSS
constraints, upstream-V availability, or deletion of intervening loci. Carried
**novel** alleles on the drawn chromosome are valid replacement targets. Receptor
revision works the same way inside [`run_cohort`](genotype-cohorts.md) (per
subject) and is still not combined with the clonal forks.

## Novel / private alleles

Individuals carry germline alleles that aren't in any reference - *private* or
*novel* alleles. Discovering them is a central task for IgDiscover, partis, and
TIgGER's `findNovelAlleles`. GenAIRR can plant them as ground truth.

`add_novel_allele` derives a private allele from a reference **base** allele by
applying point `mutations` (or supplying an explicit `sequence` of the same
length), inheriting the base's gene, anchor, functional status and V sub-regions.
The novel allele is then placed like any allele, and at `compile()` time it is
injected into an **effective reference** (base catalogue + your private alleles)
so it flows through alignment and AIRR output as a genuine allele:

```python
g = (
    Genotype.from_dataconfig(cfg)
      .add_novel_allele("IGHVF1-G1*i01", base="IGHVF1-G1*01",
                        mutations=[(38, "C"), (41, "A")])   # two point variants
      .complete_from_reference()
      .heterozygous("IGHVF1-G1", "IGHVF1-G1*01", "IGHVF1-G1*i01")  # one reference + one private
      .with_subject("DONOR_N")
)

result = (
    ga.Experiment.on(cfg).with_genotype(g).recombine()
      .run_records(n=500, seed=3, expose_provenance=True)
)
# The private allele is sampled, assembled and reported like any allele -
# its name appears in v_call / truth_v_call and the reads carry its variants.
```

The novel allele's **gene is taken from its name** and must match the base
allele's gene; it must be a same-length (substitution-only) variant. The
synthesized coding sequence is **validated** - for V/J the conserved anchor codon
must still encode the conserved residue (Cys for V, Trp/Phe for J) and the coding
frame must be stop-free. A variant that breaks either is rejected unless you pass
`allow_nonfunctional=True` (then it is kept and marked non-functional). Novel
alleles are flagged in the ground truth: each `to_table()`/`to_tsv()` row carries
a `novel` list of the private alleles carried at that gene.

**Benchmarking novel-allele discovery.** Plant a novel allele, simulate, then run
the discovery tool against the **base** germline (the cartridge *without* your
private alleles) so the tool must rediscover it from the reads - and score its
output against the planted novel sequence. (Write the base germline FASTA from
`cfg.v_alleles`; write the truth from `genotype.to_table()`.)

## Ground truth and provenance

A genotype experiment emits, by construction, everything an evaluation needs:

- **Per-record fields:** `subject_id` and `haplotype` are **always** stamped on
  every record; the truth/provenance columns (`truth_v_call` / `truth_d_call` /
  `truth_j_call`, and `original_v_call` / `receptor_revision_applied` when
  revision is used) require `expose_provenance=True`.
- **AIRR-standard vs GenAIRR extensions:** `sequence_id`, `sequence`,
  `sequence_alignment`, `v_call` / `d_call` / `j_call` etc. are standard
  [AIRR Rearrangement](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)
  fields. `haplotype` (a `0`/`1` chromosome index - **not** the AIRR
  `*_germline_alignment`/haplotype-set sense), `subject_id`, `truth_*`,
  `original_v_call`, and `receptor_revision_applied` are **GenAIRR extension
  columns** added for ground-truth benchmarking; you won't find them in the AIRR
  schema.
- **`result.genotypes`:** the list of attached `Genotype` objects (one per subject).
- **`Genotype.to_table()` / `to_tsv(path)`:** the ground-truth genotype as a table
  one row per (segment, gene) with `zygosity`
  (`homozygous` / `heterozygous` / `hemizygous` / `deleted`), the carried alleles
  per haplotype, and per-haplotype `allele:copies:weight` detail. This is the
  reference a genotype-inference benchmark compares against.

```python
for row in result.genotypes[0].to_table():
    if row["zygosity"] != "homozygous":      # show the interesting genes
        print(row["gene"], row["zygosity"], row["haplotype_0"], row["haplotype_1"])
```

## More genotype topics

- **[Genotype sampling & population priors](genotype-priors.md)** - draw a
  plausible genotype with `Genotype.sample`, or author a donor-population prior
  on the cartridge so sampling is data-driven.
- **[Genotype cohorts](genotype-cohorts.md)** - `run_cohort` to simulate many
  subjects, each with their own genotype, in one call.
- **[Benchmarking genotype inference](genotype-benchmarking.md)** - the
  end-to-end recipe + a worked TIgGER / IgDiscover example recovering a planted
  genotype.

## Limitations (this release)

**Model assumptions** (what the genotype machinery does *not* try to capture, so
you can judge whether it fits your study):

- **Sampling is independent per gene.** `Genotype.sample` draws each gene
  independently under a per-chromosome Hardy–Weinberg model - **no** linkage
  disequilibrium, gene co-deletion blocks, ancestry, or donor-specific haplotype
  structure (see the caveat box in
  [Sampling & population priors](genotype-priors.md)). Hand-built genotypes are
  exactly what you specify.
- **Deletion, not duplication, in sampling.** Sampling models gene presence/
  deletion only; gene duplication / copy-number > 1 must be built explicitly with
  `duplicate_gene`.
- **Receptor revision is haplotype-aware, not mechanistic.** It restricts the
  replacement V to carried alleles on the drawn chromosome, but does not model
  genomic V order, RSS constraints, upstream-V availability, or deletion of
  intervening loci.

**Deferred features:**

- **External loaders** - importing genotypes from VDJbase / TIgGER / IgDiscover /
  partis output. (You can build the equivalent `Genotype` by hand today.)

## Backward compatibility

The genotype machinery is purely additive. An experiment with **no** genotype
attached produces byte-identical output to previous releases (pinned by a
checksum test). Attaching a genotype is the only thing that switches recombination
onto the phased path.
