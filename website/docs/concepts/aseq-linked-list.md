---
title: The ASeq Linked List
sidebar_label: ASeq Linked List
---

# The ASeq Linked List

At the heart of GenAIRR's C engine is a data structure called **ASeq** (Annotated Sequence). It's a doubly-linked list where each node represents a single nucleotide, and each node carries all the metadata about that nucleotide's origin, identity, and history.

<div className="callout-card">
  <div className="cc-title">The core design principle</div>
  <div className="cc-body">
    Metadata lives <strong>in the nodes</strong>, not in external dictionaries. When you mutate, insert, delete, or corrupt a nucleotide, its metadata moves with it. There are no external arrays to keep in sync, no coordinate offsets to recalculate. The annotations <em>are</em> the sequence.
  </div>
</div>

## The Nuc node

Every nucleotide in a simulated sequence is a `Nuc` struct — 44 bytes carrying everything the engine needs to know about that position.

<div className="node-diagram">
  <div className="nd-title">Nuc — Annotated Nucleotide Node</div>
  <div className="nd-field"><span className="nd-name">current</span><span className="nd-type">char</span><span className="nd-desc">The base right now (A / C / G / T / N)</span></div>
  <div className="nd-field"><span className="nd-name">germline</span><span className="nd-type">char</span><span className="nd-desc">Original germline base ('\0' for NP-region nodes)</span></div>
  <div className="nd-field"><span className="nd-name">segment</span><span className="nd-type">Segment</span><span className="nd-desc">V, NP1, D, NP2, J, C, UMI, or ADAPTER</span></div>
  <div className="nd-field"><span className="nd-name">germline_pos</span><span className="nd-type">uint16</span><span className="nd-desc">Position within the germline allele sequence</span></div>
  <div className="nd-field"><span className="nd-name">flags</span><span className="nd-type">uint16</span><span className="nd-desc">Bitmask of events (mutated, error, anchor, ...)</span></div>
  <div style={{borderTop: '1px dashed rgba(0,0,0,0.1)', margin: '0.5rem 0'}} />
  <div className="nd-field"><span className="nd-name">prev</span><span className="nd-type">*Nuc</span><span className="nd-desc">← previous node in the doubly-linked list</span></div>
  <div className="nd-field"><span className="nd-name">next</span><span className="nd-type">*Nuc</span><span className="nd-desc">→ next node in the doubly-linked list</span></div>
  <div style={{borderTop: '1px dashed rgba(0,0,0,0.1)', margin: '0.5rem 0'}} />
  <div className="nd-field"><span className="nd-name">frame_phase</span><span className="nd-type">uint8</span><span className="nd-desc">0, 1, or 2 — position within reading frame</span></div>
  <div className="nd-field"><span className="nd-name">amino_acid</span><span className="nd-type">char</span><span className="nd-desc">Translated amino acid (phase-0 nodes only)</span></div>
  <div className="nd-field"><span className="nd-name">codon_next</span><span className="nd-type">*Nuc</span><span className="nd-desc">Skip pointer → next phase-0 node (codon rail)</span></div>
  <div className="nd-field"><span className="nd-name">productive</span><span className="nd-type">bool</span><span className="nd-desc">False if part of a stop codon or broken anchor</span></div>
</div>

### Primary fields

**`current`** — the base as it exists right now. This starts as the germline base but can be changed by SHM, PCR errors, sequencing errors, or N-insertion. This is what ends up in the output `sequence` field.

**`germline`** — the original germline base, frozen at assembly time. This field is **never modified** after initial assignment (with one narrow exception: D-inversion and receptor revision update it because they change the germline reference itself). For NP-region nodes, germline is `'\0'` because these bases have no germline origin. This is what ends up in the `germline_alignment` field.

**`segment`** — which segment this nucleotide belongs to. The segment tag is permanent — it never changes. Even after 5' corruption removes the first 20 V nodes, the remaining V nodes still know they're V.

<div style={{display: 'flex', flexWrap: 'wrap', gap: '0.5rem', margin: '1rem 0'}}>
  <span className="seg-chip seg-v">SEG_V</span>
  <span className="seg-chip seg-np">SEG_NP1</span>
  <span className="seg-chip seg-d">SEG_D</span>
  <span className="seg-chip seg-np">SEG_NP2</span>
  <span className="seg-chip seg-j">SEG_J</span>
  <span className="seg-chip seg-v">SEG_C</span>
  <span className="seg-chip seg-umi">SEG_UMI</span>
  <span className="seg-chip seg-adp">SEG_ADAPTER</span>
</div>

**`germline_pos`** — the position of this base within the germline allele sequence. For a V allele like `IGHVF10-G50*04`, germline_pos tells you exactly which base of that allele this node represents. At serialization, the engine reads the germline_pos of the first and last V nodes to compute `v_germline_start` and `v_germline_end` — no counting needed.

### Event flags

The `flags` field is a bitmask tracking everything that has happened to this nucleotide. Flags accumulate — a node can be both `MUTATED` and `SEQ_ERROR` if an SHM mutation was later hit by a sequencing error.

<div style={{display: 'flex', flexWrap: 'wrap', gap: '0.3rem', margin: '1rem 0'}}>
  <span className="flag-badge">MUTATED</span>
  <span className="flag-badge">SEQ_ERROR</span>
  <span className="flag-badge">PCR_ERROR</span>
  <span className="flag-badge">P_NUCLEOTIDE</span>
  <span className="flag-badge">N_NUCLEOTIDE</span>
  <span className="flag-badge">IS_N</span>
  <span className="flag-badge">ANCHOR</span>
  <span className="flag-badge">INDEL_INS</span>
</div>

| Flag | Bit | Set by |
|------|-----|--------|
| `MUTATED` | 0 | S5F somatic hypermutation |
| `SEQ_ERROR` | 1 | Position-dependent sequencing errors |
| `PCR_ERROR` | 2 | PCR polymerase errors |
| `P_NUCLEOTIDE` | 3 | Palindromic nucleotide at junction |
| `N_NUCLEOTIDE` | 4 | TdT-generated N-nucleotide (NP region) |
| `IS_N` | 5 | Base corrupted to ambiguous N |
| `ANCHOR` | 6 | Conserved anchor (V-Cys or J-Trp/Phe) |
| `INDEL_INS` | 7 | Inserted by indel operation |

The serialization code uses these flags to separate mutations, PCR errors, and sequencing errors into their respective output fields — a clean separation that requires no post-hoc heuristics.

## A sequence as a linked list

Here's what a real assembled sequence looks like as an ASeq. Each box is a Nuc node, colored by segment:

<div className="seq-vis">
  <div>
    <span className="sv-label">Segment</span>
    <span className="sv-v">V V V V V V V V V</span>
    <span> </span>
    <span className="sv-np">N N N N</span>
    <span> </span>
    <span className="sv-d">D D D D D</span>
    <span> </span>
    <span className="sv-np">N N N</span>
    <span> </span>
    <span className="sv-j">J J J J J J</span>
  </div>
  <div>
    <span className="sv-label">Current</span>
    <span className="sv-v">g a g g t g c a g</span>
    <span> </span>
    <span className="sv-np">T G T C</span>
    <span> </span>
    <span className="sv-d">a t a t t</span>
    <span> </span>
    <span className="sv-np">T T A</span>
    <span> </span>
    <span className="sv-j">a c t a c t</span>
  </div>
  <div>
    <span className="sv-label">Germline</span>
    <span className="sv-v">g a g g t g c a g</span>
    <span> </span>
    <span className="sv-n">∅ ∅ ∅ ∅</span>
    <span> </span>
    <span className="sv-d">a t a t t</span>
    <span> </span>
    <span className="sv-n">∅ ∅ ∅</span>
    <span> </span>
    <span className="sv-j">a c t a c t</span>
  </div>
  <div>
    <span className="sv-label">Germline Pos</span>
    <span className="sv-v">0 1 2 3 4 5 6 7 8</span>
    <span> </span>
    <span className="sv-n">- - - -</span>
    <span> </span>
    <span className="sv-d">3 4 5 6 7</span>
    <span> </span>
    <span className="sv-n">- - -</span>
    <span> </span>
    <span className="sv-j">9 . . . . .</span>
  </div>
</div>

NP nodes have `germline = '\0'` (shown as ∅) and no germline position — they were created by random nucleotide generation, not copied from an allele.

### After SHM mutation

When S5F mutates position 6 (V node, `c` → `T`), only that one node changes:

<div className="seq-vis">
  <div>
    <span className="sv-label">Current</span>
    <span className="sv-v">g a g g t g </span><span className="sv-mut">T</span><span className="sv-v"> a g</span>
    <span> </span>
    <span className="sv-np">T G T C</span>
    <span> ...</span>
  </div>
  <div>
    <span className="sv-label">Germline</span>
    <span className="sv-v">g a g g t g c a g</span>
    <span> </span>
    <span className="sv-n">∅ ∅ ∅ ∅</span>
    <span> ...</span>
  </div>
  <div>
    <span className="sv-label">Flags</span>
    <span className="sv-v">· · · · · · </span><span className="sv-mut">M</span><span className="sv-v"> · ·</span>
    <span> </span>
    <span className="sv-np">N N N N</span>
    <span> ...</span>
  </div>
</div>

The `germline` field still shows `c`. The `flags` field now has `MUTATED` set. At serialization, this produces the mutation string entry `6:c>T`.

## The ASeq container

<div className="node-diagram">
  <div className="nd-title">ASeq — Annotated Sequence Container</div>
  <div className="nd-field"><span className="nd-name">pool[1024]</span><span className="nd-type">Nuc[]</span><span className="nd-desc">Pre-allocated arena for all nodes (no malloc)</span></div>
  <div className="nd-field"><span className="nd-name">pool_used</span><span className="nd-type">int</span><span className="nd-desc">Next free slot in the pool</span></div>
  <div style={{borderTop: '1px dashed rgba(0,0,0,0.1)', margin: '0.5rem 0'}} />
  <div className="nd-field"><span className="nd-name">head / tail</span><span className="nd-type">*Nuc</span><span className="nd-desc">First and last nodes in the linked list</span></div>
  <div className="nd-field"><span className="nd-name">length</span><span className="nd-type">int</span><span className="nd-desc">Number of active (linked) nodes</span></div>
  <div style={{borderTop: '1px dashed rgba(0,0,0,0.1)', margin: '0.5rem 0'}} />
  <div className="nd-field"><span className="nd-name">seg_first[8]</span><span className="nd-type">*Nuc[]</span><span className="nd-desc">First node of each segment — O(1) boundary access</span></div>
  <div className="nd-field"><span className="nd-name">seg_last[8]</span><span className="nd-type">*Nuc[]</span><span className="nd-desc">Last node of each segment — O(1) boundary access</span></div>
  <div style={{borderTop: '1px dashed rgba(0,0,0,0.1)', margin: '0.5rem 0'}} />
  <div className="nd-field"><span className="nd-name">codon_rail_valid</span><span className="nd-type">bool</span><span className="nd-desc">Is the codon rail current?</span></div>
  <div className="nd-field"><span className="nd-name">n_stop_codons</span><span className="nd-type">int</span><span className="nd-desc">Live count — updated on every mutation</span></div>
  <div className="nd-field"><span className="nd-name">v_anchor_node</span><span className="nd-type">*Nuc</span><span className="nd-desc">Cached conserved Cys (junction start)</span></div>
  <div className="nd-field"><span className="nd-name">j_anchor_node</span><span className="nd-type">*Nuc</span><span className="nd-desc">Cached conserved W/F (junction end)</span></div>
</div>

### Arena allocation

All Nuc nodes come from a contiguous `pool[1024]` array — no `malloc` per node. This gives cache-friendly allocation and zero-cost cleanup (just reset `pool_used` to 0). The maximum sequence length of 1024 nucleotides is more than enough for any immunoglobulin or TCR sequence.

### Segment boundary cache

`seg_first[8]` and `seg_last[8]` cache the first and last node of each segment type, giving O(1) access to any segment boundary. These are automatically maintained by all insert/delete operations. When serializing, the engine reads `seg_first[SEG_V]` to find where V starts — no scanning needed.

## The codon rail

The codon rail is a secondary structure overlaid on the nucleotide list. It tracks the reading frame and provides **O(1) amino acid updates** on point mutations.

```mermaid
graph LR
  subgraph codon_rail["Codon Rail"]
    direction LR
    N0(["G  phase 0  Glu"])
    N1["a  ph 1"]
    N2["g  ph 2"]
    N3(["G  phase 0  Val"])
    N4["t  ph 1"]
    N5["g  ph 2"]
    N6(["C  phase 0  Gln"])
    N7["a  ph 1"]
    N8["g  ph 2"]
  end

  N0 --> N1 --> N2 --> N3 --> N4 --> N5 --> N6 --> N7 --> N8
  N0 -.->|codon_next| N3
  N3 -.->|codon_next| N6

  style N0 fill:#2563eb,stroke:#1d4ed8,color:#fff,stroke-width:2px
  style N3 fill:#2563eb,stroke:#1d4ed8,color:#fff,stroke-width:2px
  style N6 fill:#2563eb,stroke:#1d4ed8,color:#fff,stroke-width:2px
```

After assembly, `aseq_build_codon_rail()` walks the entire list, assigning each node a `frame_phase` (0, 1, 2, repeating). Phase-0 nodes are **codon heads** — they hold the translated amino acid and link to the next phase-0 node via `codon_next` skip pointers.

### O(1) mutation updates

When `aseq_mutate()` changes a base, it doesn't retranslate the entire sequence. It:

1. Finds the codon head for the mutated node (walk back 0-2 nodes using `frame_phase`)
2. Retranslates just that one codon (3 bases → 1 amino acid)
3. Updates `n_stop_codons` if the amino acid changed to/from `*`

This means checking productivity after a mutation is effectively free.

### Frame shift propagation

Insertions and deletions shift the reading frame for all downstream nodes. The engine walks forward from the insertion/deletion point, reassigning `frame_phase` and retranslating each affected codon. For batch operations (5'/3' corruption), the codon rail is rebuilt once from scratch instead of propagating individually.

## Operations

Every pipeline step manipulates the ASeq through a small set of primitives:

| Operation | What it does | Metadata effect | Cost |
|-----------|-------------|----------------|------|
| `aseq_mutate()` | Change `current`, set flag | `germline` untouched, 1 codon retranslated | O(1) |
| `aseq_revert()` | Restore `current` to `germline` | Clear `MUTATED` flag, retranslate codon | O(1) |
| `aseq_insert_after()` | Splice new node into list | New node gets `germline='\0'`, propagate frame shift | O(n) |
| `aseq_delete()` | Unlink node from list | Node vanishes with its metadata, propagate frame | O(n) |
| `aseq_delete_head_n()` | Batch remove from 5' end | Rebuild codon rail once | O(n) |
| `aseq_delete_tail_n()` | Batch remove from 3' end | Rebuild codon rail once | O(n) |
| `aseq_reverse_complement()` | Reverse list + complement bases | Flags preserved, codon rail invalidated | O(n) |

<div className="callout-card">
  <div className="cc-title">Why this design works</div>
  <div className="cc-body">
    In many simulators, the sequence is a flat string and annotations live in parallel arrays. Every mutation requires updating these external structures — a constant source of bugs. In GenAIRR, there's nothing to synchronize. When you delete a node, its metadata disappears with it. When you mutate a node, the germline field is untouched. The engine can apply 10+ pipeline steps and still produce metadata that matches the sequence at every position — because the metadata was never separated from the sequence in the first place.
  </div>
</div>
