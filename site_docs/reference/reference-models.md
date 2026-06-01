# Reference models and rules

<p class="lead">The typed dataclasses that describe a cartridge's
empirical-model plane (allele usage, NP-base model, trim and
length distributions) and its rules plane (anchor expectations,
alphabet, severity policy). These specs are the schemas the
<a href="reference-cartridge-builder.md"><code>ReferenceCartridgeBuilder</code></a>
estimators write into; you can also hand-author them when
estimating from records isn't possible.</p>

For the conceptual model — four typed planes (identity,
catalogue, rules, empirical models) — see
[the reference-cartridge concept page](../concepts/reference-cartridge.md).

## Empirical models

### `ReferenceEmpiricalModels`

The bundle that carries the four empirical planes. One of these
attaches to `DataConfig.reference_models` after `.build()`.

::: GenAIRR.reference_models.ReferenceEmpiricalModels
    options:
      show_source: false
      show_root_heading: true
      heading_level: 4
      members: false

### `EmpiricalDistributionSpec`

The `[(value, weight), ...]` shape used for trim and NP-length
distributions. Weights renormalise; zero-weight entries drop out
of the proposal support.

::: GenAIRR.reference_models.EmpiricalDistributionSpec
    options:
      show_source: false
      show_root_heading: true
      heading_level: 4
      members: false

### `NpBaseModelSpec`

The NP-base sampling model. Two kinds: empirical-first-base
(per-base marginal) and Markov (first-base marginal + 4×4
transition matrix). The Markov form is the default produced by
`estimate_np_base_model(records, kind="markov")`.

::: GenAIRR.reference_models.NpBaseModelSpec
    options:
      show_source: false
      show_root_heading: true
      heading_level: 4
      members: false

### `AlleleUsageSpec`

Per-segment allele weights for V, D, and J. Weights renormalise
within each segment.

::: GenAIRR.reference_models.AlleleUsageSpec
    options:
      show_source: false
      show_root_heading: true
      heading_level: 4
      members: false

## Rules

### `ReferenceRulesSpec`

Anchor expectations + alphabet + severity policy. Hand-authored
or attached via `.with_rules(...)` on the builder.

::: GenAIRR.reference_rules.ReferenceRulesSpec
    options:
      show_source: false
      show_root_heading: true
      heading_level: 4
      members: false

### `AnchorRuleSpec`

Per-anchor expected amino-acid plus a `required` flag. Members of
`ReferenceRulesSpec.anchors` are `AnchorRuleSpec` instances.

::: GenAIRR.reference_rules.AnchorRuleSpec
    options:
      show_source: false
      show_root_heading: true
      heading_level: 4
      members: false
