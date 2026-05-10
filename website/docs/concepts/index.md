---
title: Concepts
sidebar_label: Overview
---

# Concepts

This section explains **how GenAIRR works internally** — the data structures, algorithms, and design decisions that make every output field surgically accurate.

You don't need to read this to use GenAIRR. But if you want to understand why the metadata is always correct, how mutations interact with the reading frame, or how coordinates survive corruption — this is where you'll find the answers.

| Page | What it covers |
|------|---------------|
| [Persistent IR Architecture](/docs/concepts/persistent-ir) | The heart of the engine — a persistent, immutable data structure that tracks the nucleotide pool, segment regions, and allele assignments across a complete simulation history. |
| [The Simulation Pipeline](/docs/concepts/simulation-pipeline) | How the engine executes a modular **PassPlan** — from recombination and somatic hypermutation to sequencing artifacts — ensuring every step is tracked and reproducible. |
| [How Sampling Works](/docs/concepts/sampling-mechanics) | Transparency into the engine's internal sampling logic: empirical weights, Markov-chain trimming, and TdT transition matrices. |
| [Metadata Accuracy](/docs/concepts/metadata-accuracy) | How the engine derives ~70 ground-truth fields (coordinates, CIGARs, alignments) through its single-pass builder, maintaining perfect consistency even under extreme mutation and corruption. |
