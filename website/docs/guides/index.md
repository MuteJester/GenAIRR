---
title: Guides
sidebar_label: Overview
---

# Guides

Practical recipes for common GenAIRR tasks. Each guide is self-contained — pick the one that matches what you need to do.

| Section | Guide | Description |
|---------|-------|-------------|
| **Simulation Basics** | [The Experiment DSL](/docs/guides/basics/experiment-dsl) | How to build simulation pipelines using the fluent `Experiment` API. |
| | [Chain Types](/docs/guides/basics/chain-types) | Understanding the structural differences between VDJ and VJ rearrangements. |
| | [Export Formats](/docs/guides/basics/export) | Saving results as AIRR TSV, FASTA, or pandas DataFrames, and memory-efficient streaming. |
| **Biological Variation** | [Somatic Hypermutation](/docs/guides/options/shm) | Deep dive into context-dependent S5F and uniform mutation models. |
| | [Biological Events](/docs/guides/options/biology) | Controlling biological sampling via Allele Locking and Weighting. |
| | [Clonal Families](/docs/guides/options/clonal-families) | Simulating hierarchical lineages and shared mutation history. |
| **Technical Realism** | [Sequencing Artifacts](/docs/guides/options/artifacts) | Modeling technical noise like primer trimming, PCR errors, and indels. |
| | [Reproducibility](/docs/guides/options/reproducibility) | Ensuring bit-for-bit identity across platforms using seeded PRNGs and the **Addressed Trace**. |
| **Safety & Constraints** | [Productive Sequences](/docs/guides/options/productive) | Using **Contracts** to ensure functional, in-frame rearrangements by construction. |
