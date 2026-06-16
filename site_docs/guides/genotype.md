# Genotypes: per-individual diploid germline

<p class="lead">A <strong>genotype</strong> in GenAIRR is one person's diploid germline
complement — which V/D/J alleles they carry, on which chromosome, in what copy
number. Attach a <code>Genotype</code> to an <code>Experiment</code> and V(D)J
recombination becomes <strong>haplotype-phased</strong>: the V, D and J of each
rearrangement are drawn from a single chromosome, honouring allele
presence/absence, zygosity, and gene deletion. With no genotype attached the
engine is byte-for-byte unchanged. This page explains exactly what a genotype is,
how the engine samples from it, how to build one, and how to use it to benchmark
genotype-inference tools — nothing here is a black box.</p>

## What a genotype is (and why it matters)

Every person inherits two copies of the immunoglobulin heavy-chain locus — one on
each homologous chromosome (one **haplotype** from each parent). Across the
population the locus is extraordinarily polymorphic: a reference set may list
dozens of alleles per gene, but a *single individual* carries only a handful —
typically **one or two alleles per gene** — and may be **missing entire genes**
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

result[0]["subject_id"]          # 'DONOR01'   — provenance on every record
result[0]["haplotype"]           # 0 or 1      — which chromosome this read used
result.genotypes[0].to_table()   # ground-truth genotype (per gene, per haplotype)
```

## How recombination samples from a genotype

When a genotype is attached, recombination runs a single phased sampling pass per
rearrangement. The steps, in order:

1. **Draw a chromosome.** One of the two haplotypes is chosen, weighted by the
   `chromosome_weights` (default `[0.5, 0.5]`). This choice is made **once** and
   shared by V, D and J — that is the phasing.
2. **Per segment, draw a gene then an allele.** Among the genes *present on the
   chosen chromosome*, a gene is sampled (weighted by usage — see below), then the
   allele follows from that chromosome's slot for the gene. A deleted gene is
   simply not offered on the chromosome that lacks it.
3. **Assign and continue.** The chosen V/D/J alleles are assigned and the rest of
   the pipeline (trimming, NP, assembly, SHM, corruption) runs unchanged.

Every random choice (chromosome, gene, within-slot allele) is recorded to the
trace, so seeded runs are byte-stable and fully replayable.

### Viability and `productive_only`

A chromosome is only drawn if it is **viable** — it must carry at least one usable
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
recombination but was never specified is an error at attach time — you must define
the whole genotype (use `complete_from_reference` to fill the genes you don't care
about). This guarantees a genuine diploid complement, which is what you want for a
ground-truth benchmark.

`Genotype.permissive(cfg)` is a separate, explicitly **non-diploid** fallback:
unspecified genes are left to sample over *all* their reference alleles without
phasing. It exists for the "I only want to constrain a few genes" case — it is
**not** a biological genotype, and it is labelled as such in `to_table()` and
`repr`.

In both modes, feasibility (e.g. `productive_only`) is applied the same way the
non-genotype path applies it: candidates are filtered to the feasible set, and the
unfiltered set is used only as a last resort when nothing is feasible — so a
genotype run never silently samples alleles a normal run would have avoided.

## Building a genotype

The `Genotype` builder is a fluent, validated editor over a `DataConfig`'s
reference alleles. Every method returns `self`, so calls chain.

```python
g = Genotype.from_dataconfig(cfg)                       # strict (recommended)

g.homozygous("IGHVF2-G4", "IGHVF2-G4*01")               # 1 allele on both chromosomes
g.heterozygous("IGHVF1-G1", "IGHVF1-G1*01", "IGHVF1-G1*02")  # different allele per chromosome
g.delete_gene("IGHVF3-G7", haplotype="both")            # gene absent entirely (homozygous deletion)
g.delete_gene("IGHVF3-G8", haplotype=1)                 # absent on chromosome 1 only (hemizygous)
g.duplicate_gene("IGHVF1-G2", ["IGHVF1-G2*01", "IGHVF1-G2*03"], haplotype=0)  # >1 copy on one chromosome
g.chromosome_weights(0.6, 0.4)                          # allelic-expression imbalance
g.with_subject("DONOR01")                               # provenance label

g.complete_from_reference("homozygous_first_reference") # fill every unspecified gene
```

Notes and guard-rails:

- **`delete_gene(..., haplotype=0|1)`** (one chromosome) requires the gene to be
  specified first — deleting a single haplotype of an *unspecified* gene would
  silently delete both, so it raises instead.
- **`complete_from_reference(policy=...)`** fills only genes you haven't touched.
  `"homozygous_first_reference"` (default) makes each unspecified gene homozygous
  for its first cartridge allele — note this is the *first listed* allele, **not**
  a population-frequency-common one (GenAIRR has no frequency prior in this
  release; the name says exactly what it does). `"heterozygous_first_two"` uses
  the first two alleles.
- Unknown gene/allele names, NaN/inf chromosome weights, and segments left with no
  usable allele are all rejected at build/attach with clear messages.
- `with_genotype` is **mutually exclusive** with
  [`restrict_alleles`](../reference/experiment.md) and the
  `recombine(*_allele_weights=...)` kwargs — the genotype owns allele presence and
  expression. It is also rejected together with `receptor_revision` and with the
  clonal forks (`expand_clones` / `clonal_lineage` / `clonal_repertoire`) in this
  release (see [Limitations](#limitations-this-release)).

### Gene usage

Within a chromosome, which *gene* is used is weighted by the cartridge's typed
allele-usage model (`reference_models.allele_usage`), aggregated to the gene
level, and scaled by **copy-number dosage** (a duplicated gene recombines
proportionally more often). Cartridges that don't author a typed `allele_usage`
fall back to uniform-over-present-genes (× dosage). See
[Allele usage](v-usage.md) and [Estimate models from data](estimate-cartridge-models.md)
for authoring usage.

## Ground truth and provenance

A genotype experiment emits, by construction, everything an evaluation needs:

- **Per-record fields:** `subject_id` and `haplotype` (`0`/`1`, the chromosome the
  rearrangement used) are stamped on every AIRR record. Standard truth columns
  (`truth_v_call`, …) are available with `expose_provenance=True`.
- **`result.genotypes`:** the list of attached `Genotype` objects (one per subject).
- **`Genotype.to_table()` / `to_tsv(path)`:** the ground-truth genotype as a table
  — one row per (segment, gene) with `zygosity`
  (`homozygous` / `heterozygous` / `hemizygous` / `deleted`), the carried alleles
  per haplotype, and per-haplotype `allele:copies:weight` detail. This is the
  reference a genotype-inference benchmark compares against.

```python
for row in result.genotypes[0].to_table():
    if row["zygosity"] != "homozygous":      # show the interesting genes
        print(row["gene"], row["zygosity"], row["haplotype_0"], row["haplotype_1"])
```

## Research workflow: benchmarking genotype inference

The point of simulating from a *known* genotype is that you can run a
genotype-inference tool on the resulting repertoire and score it against the
planted truth — with no real-data uncertainty about what the right answer is.

The recipe is the same for any tool:

1. Build a `Genotype`, simulate a repertoire, write the AIRR table
   (`result.to_tsv(...)`) and/or reads FASTA, and the ground truth
   (`genotype.to_tsv(...)`).
2. Run the inference tool to recover the per-individual allele set.
3. Compare recovered vs planted: presence/absence, zygosity, and (for
   discovery tools) any novel alleles.

### Worked example: TIgGER and IgDiscover recover a planted genotype

To show this end to end we planted a diploid IGH genotype in `human_igh` —
**3 heterozygous** V genes (two alleles each), **3 homozygous** (one allele),
and **3 fully deleted** genes — and filled the rest from the reference. We
simulated 4,000 reads with light SHM, then ran two independent AIRR
genotype-inference tools on the result:
[**TIgGER**](https://tigger.readthedocs.io) (Immcantation; consumes the AIRR
table) and [**IgDiscover**](https://igdiscover.se) (germline discovery from the
raw reads, with its own IgBLAST).

Because GenAIRR already emits AIRR records with `v_call` **and**
`sequence_alignment`, TIgGER's `inferGenotype` consumes the rearrangement table
**directly — no separate IgBLAST step is needed**:

```r
library(tigger); library(airr)
rep    <- read_rearrangement("repertoire.tsv")     # GenAIRR's AIRR output
germ_v <- readIgFasta("germline_V.fasta")          # cartridge V germline (names match v_call)
geno   <- inferGenotype(rep, germline_db = germ_v, find_unmutated = TRUE)
plotGenotype(geno)
```

TIgGER recovered the planted genotype **exactly**: every heterozygous gene → two
alleles, every homozygous gene → one, every deleted gene → **absent**. Across all
52 V genes, allele-presence **precision = 1.00**, **recall = 1.00**, and the
per-gene allele count matched the truth for **52/52** genes.

**IgDiscover**, run on the raw reads with the cartridge as its starting database,
independently agreed: **precision = 1.00** (zero false-positive alleles),
**recall = 0.96** (50/52 carried alleles), with **all three deletions correct**
and **all heterozygous genes fully resolved** (both alleles recovered). The two
missed alleles were low-expression single-copy genes below IgDiscover's default
expression threshold — a tool-tuning matter, not a simulation artefact.

![GenAIRR-simulated genotype recovered by TIgGER and IgDiscover: planted vs inferred allele counts agree for every gene](../assets/genotype-tigger-recovery.png)

*(A) The nine study genes: both tools' inferred allele counts match the planted
zygosity for each (heterozygous → 2, homozygous → 1, deleted → 0). (B) All 52 V
genes fall on the agreement diagonal; presence precision = 1.00 for both tools,
recall 1.00 (TIgGER) / 0.96 (IgDiscover), zero false-positive alleles, all
deletions correct.*

### Reproduce it

```python
import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype

cfg = gdata.HUMAN_IGH_OGRDB
g = (
    Genotype.from_dataconfig(cfg)
      .complete_from_reference("homozygous_first_reference")
      .heterozygous("IGHVF1-G1", "IGHVF1-G1*01", "IGHVF1-G1*02")
      .homozygous("IGHVF2-G4", "IGHVF2-G4*01")
      .delete_gene("IGHVF3-G7", haplotype="both")
      .with_subject("DONOR01")
)
res = (
    ga.Experiment.on(cfg).with_genotype(g).recombine()
      .mutate(rate=0.004)                       # light SHM, as in real data
      .run_records(n=4000, seed=7)
)
res.to_tsv("repertoire.tsv")                    # AIRR table → TIgGER
g.to_tsv("truth_genotype.tsv")                  # ground truth to score against

# export the cartridge V germline (names match v_call) for TIgGER's germline_db
with open("germline_V.fasta", "w") as fh:
    for gene, alleles in cfg.v_alleles.items():
        for a in alleles:
            fh.write(f">{a.name}\n{a.ungapped_seq.upper()}\n")
```

Then run the R snippet above and compare `geno` against `truth_genotype.tsv`.

### Running other tools on the same data

The only difference between tools is whether they consume the **AIRR table**
(TIgGER) or the **raw reads** (`reads.fasta`, which you can write from `result`),
running their own aligner:

- **[IgDiscover](https://igdiscover.se)** — germline *discovery* from reads (its
  own IgBLAST + iterative filtering). Initialise with the cartridge germline as
  the starting database and the simulated reads, then run the pipeline; the
  `final/database/V.fasta` expressed-allele set is the recovered genotype:

  ```bash
  igdiscover init --database db/ --single-reads reads.fasta project/
  cd project && igdiscover run
  ```

- **[partis](https://github.com/psathyrella/partis)** — HMM annotation with
  per-sample germline inference (`partis cache-parameters --infname reads.fa
  --initial-germline-dir db/`). partis also reports per-sample allele support and
  novel alleles, scored the same way.

Because the genotype is planted, every tool is scored identically: recovered
allele set vs `genotype.to_table()` — presence precision/recall, zygosity, and
deletion calls.

## Limitations (this release)

The genotype foundation is deliberately scoped. Deferred to later work:

- **Novel / private alleles** — per-individual germline variants not in the
  reference (synthesis + provenance). Today a genotype draws from reference
  alleles only.
- **Cohorts** — many subjects, each with their own genotype, in one run
  (`with_genotype` is single-subject; `result.genotypes` is a one-element list).
- **Population priors** — sampling a plausible diploid genotype from
  allele/deletion frequencies (today genotypes are specified explicitly).
- **External loaders** — importing genotypes from VDJbase / TIgGER / IgDiscover /
  partis output.
- **Cartridge genotype plane** — persisting a population genotype model in a
  cartridge.
- **Same-haplotype receptor revision** — `receptor_revision` with a genotype is
  rejected for now.

## Backward compatibility

The genotype machinery is purely additive. An experiment with **no** genotype
attached produces byte-identical output to previous releases (pinned by a
checksum test). Attaching a genotype is the only thing that switches recombination
onto the phased path.
