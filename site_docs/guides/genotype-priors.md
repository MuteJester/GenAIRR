# Genotype sampling & population priors

Two related ways to get a per-individual genotype without pinning every gene by
hand: **sample** one from population priors, or **author a population prior on
the cartridge** so sampling is data-driven. Both produce an ordinary `Genotype`
you attach with `with_genotype` (see the [Genotypes overview](genotype.md) for
the core model, and [Genotype cohorts](genotype-cohorts.md) to draw many at
once). Examples assume the imports and `cfg` from the overview.

## Sampling from population priors

Instead of specifying every gene, draw a plausible diploid genotype with
`Genotype.sample`. It uses an **independent per-gene, per-chromosome
Hardy-Weinberg model**: each gene on each chromosome is independently deleted
(`haplotype_deletion_prob`) or assigned an allele from that gene's frequencies, so
homozygous / heterozygous / hemizygous / deleted states emerge at the expected
rates.

```python
g = Genotype.sample(
    cfg,
    seed=7,
    allele_frequencies={                 # {segment: {gene: {allele: weight}}}
        "V": {"IGHVF1-G1": {"IGHVF1-G1*01": 8, "IGHVF1-G1*02": 2}},
    },                                   # unspecified genes -> uniform within gene
    haplotype_deletion_prob=0.05,        # 5% per-haplotype gene-absence
    subject_id="DONOR_R1",
)
res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(n=500, seed=1)
```

`Genotype.sample` always returns a **fully-specified, runnable** genotype. By
default `ensure_viable=True` re-draws (up to `max_resamples`, deterministically)
until at least one **expressible** (positive-weight) chromosome carries every
required segment, raising a clear error only if your deletion settings make that
impossible; pass `ensure_viable=False` to allow infeasible draws.

!!! note "Default sampling is HW *conditioned on viability*"
    With the default `ensure_viable=True`, draws that leave no complete usable
    haplotype are rejected, so per-gene deletion/zygosity rates are Hardy-Weinberg
    **conditioned on at least one viable chromosome** — not the unconditional HW
    rates. This only matters at high `haplotype_deletion_prob` (e.g. a single-J
    cartridge with a high J-deletion prob). Use `ensure_viable=False` for the raw,
    unconditional draw (which may be infeasible and rejected at compile).

Frequency priors accept the segment-aware
`{segment: {gene: {allele: weight}}}` shape (or a flat `{gene: …}` when the gene is
unambiguous); a supplied gene's listed alleles define its distribution (weight `0`
excludes an allele), and unspecified genes fall back to uniform.
`haplotype_deletion_prob` is a float or a per-gene/per-segment dict.

!!! warning "What this model is — and isn't"
    This is an **independent per-gene** sampler. It does **not** model linkage
    disequilibrium, gene co-deletion blocks, ancestry, or donor-specific haplotype
    structure, and it samples **catalogue alleles only** (no novel alleles) and
    **deletion only** (no duplication). Supply explicit `allele_frequencies` from a
    population source (e.g. VDJbase) for realistic per-gene frequencies; the default
    is uniform within each gene. `allele_frequencies="usage_as_prior"` is an
    explicit opt-in that reuses the cartridge's recombination `allele_usage` as a
    frequency proxy — convenient but biologically approximate.

## Population genotype models on a cartridge

The `allele_frequencies` / `haplotype_deletion_prob` you pass to
`Genotype.sample` can instead be **authored once on the cartridge** as a
*population genotype model* — a donor-population germline prior. It is a distinct
plane from `reference_models.allele_usage`: `allele_usage` weights how often each
allele is *expressed* during recombination, whereas a genotype prior describes
which alleles a *donor population carries* (frequencies, gene-deletion rates, and
population novel alleles). It lives on `DataConfig.genotype_priors`.

### Authoring and attaching a model

```python
from GenAIRR.genotype_priors import PopulationGenotypeModel, PopulationNovelAllele

model = PopulationGenotypeModel(
    model_id="IGH-toy-1", source="hand-authored",     # identity is required
    allele_frequencies={"V": {"IGHVF1-G1": {"IGHVF1-G1*01": 3.0, "IGHVF1-G1*02": 1.0}}},
    haplotype_deletion_prob={"V": {"IGHVF1-G1": 0.1}},
)
# Attach via the cartridge builder (validated against the chain type + catalogue):
#   cfg = builder.set_genotype_priors(model).build()
```

`set_genotype_priors` validates the model against the cartridge: unknown
genes/alleles raise, and any population novel allele goes through the same
functional validation as `add_novel_allele` (conserved anchor, stop-free frame,
base-allele match). A non-`None` plane becomes part of cartridge identity (it
folds into `compute_checksum()` and the manifest).

### Estimating a model from observed genotypes

Given a list of observed `Genotype` objects (e.g. one per donor), estimate a
prior directly:

```python
model = PopulationGenotypeModel.from_genotypes(
    [g_donor1, g_donor2, g_donor3],
    cfg=cfg, model_id="cohort-est", source="my-cohort",
)
# or, attaching to a builder in one chained step:
#   builder.estimate_genotype_priors([g_donor1, g_donor2, g_donor3], source="my-cohort")
```

Estimator conventions: allele frequencies are counted **per carried chromosome**
(homozygous contributes 2, hemizygous 1, deleted 0); gene deletion probability is
`deleted_haplotypes / (2 × n_subjects)`. `pseudocount` smooths the per-gene
catalogue-allele counts only (deletion gets none). Genotypes carrying a
duplicated gene are rejected — the plane is deletion-only.

### Drawing from the cartridge plane

When a cartridge carries a plane, `Genotype.sample(cfg)` **auto-uses it** and
records where every input came from:

```python
g = Genotype.sample(cfg, seed=7)        # plane supplies freqs / deletion / weights
print(g.prior_provenance)
# {'allele_frequencies': 'cartridge', 'haplotype_deletion_prob': 'cartridge',
#  'chromosome_weights': 'cartridge', 'novel_alleles': 'none',  # 'cartridge' if the model has novels
#  'model_id': 'IGH-toy-1', 'model_checksum': '…'}
print(g.to_metadata())                  # subject_id + provenance + source/effective refdata hashes
```

Pass `use_cartridge_priors=False` for a clean uniform, catalogue-only draw (all
plane consumption — including novels — disabled). Each input is sourced
**independently**, so you can mix explicit and cartridge values:

```python
g = Genotype.sample(
    cfg, seed=7,
    allele_frequencies={"V": {"IGHVF1-G1": {"IGHVF1-G1*01": 1.0}}},  # explicit
    # haplotype_deletion_prob and chromosome_weights left to the plane
)
print(g.prior_provenance["allele_frequencies"])      # 'explicit'
print(g.prior_provenance["haplotype_deletion_prob"]) # 'cartridge'
print(g.prior_provenance["chromosome_weights"])      # 'cartridge'
```

**Population novel alleles** on the plane are *candidate* alleles for the draw: a
novel that gets sampled is carried (and flows into `v_call` / reads / truth like
any allele); a novel that isn't drawn never pollutes the output reference. By
default (`include_cartridge_novel_alleles="auto"`) novels are injected only when
the allele frequencies are cartridge- or uniform-sourced — supplying an
**explicit** frequency table keeps it explicit. Pass
`include_cartridge_novel_alleles=True` to inject them anyway, or `False` to never.

### Auditing the plane

`cfg.cartridge_manifest()["models"]["genotype_priors"]` reports availability,
`model_id` / `source` / `version`, the plane's `model_checksum`, per-segment gene
counts, the novel-allele count, and `source_field` (`"DataConfig.genotype_priors"`).

