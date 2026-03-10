---
title: Concepts
sidebar_label: Overview
---

# Concepts

This section explains **how GenAIRR works internally** — the data structures, algorithms, and design decisions that make every output field surgically accurate.

You don't need to read this to use GenAIRR. But if you want to understand why the metadata is always correct, how mutations interact with the reading frame, or how coordinates survive corruption — this is where you'll find the answers.

| Page | What it covers |
|------|---------------|
| [The ASeq Linked List](/docs/concepts/aseq-linked-list) | The core data structure — each nucleotide is a node carrying its own segment identity, germline origin, mutation flags, and reading frame position |
| [The Simulation Pipeline](/docs/concepts/simulation-pipeline) | How a sequence flows through assembly → functionality → mutation → corruption → serialization, and what each stage does to the linked list |
| [Metadata Accuracy](/docs/concepts/metadata-accuracy) | How coordinates, germline alignment, mutation strings, and allele calls remain correct through every transformation — including boundary ambiguity resolution |
