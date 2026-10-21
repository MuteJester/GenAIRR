# Existing limitations encountered during original-end verification

Historical findings: follow-up fixes and the remaining model boundaries are
tracked in the [implementation report](productivity_indel_repair_implementation_2026_10_20.md).
The examples below describe the pre-fix behavior.

Verified on 2026-10-20 on `fix/biological-metadata-correctness`, using the
current source extension. These checks identify limitations outside the new
original-end outcome model. No cartridge or legacy sampler changes were made.

## Curatable reference anchors and legacy productivity

The bundled `mouse_tcra` (`MOUSE_TCRA_IMGT`) and `mouse_tcrg`
(`MOUSE_TCRG_IMGT`) cartridges contain incorrect anchor annotations:

| Cartridge | Reference evidence | Default compilation |
| --- | --- | --- |
| `mouse_tcra` | All 197 V alleles have `anchor=0`; none encodes C at that position. Counts: G 99, Q 45, A 31, S 7, T 6, K 4, E 2, D 2, L 1. | Rejected: 205 curatable validation issues. |
| `mouse_tcrg` | All 4 J alleles have `anchor=0`: three encode I and one encodes S. All 12 V anchors encode C; 2 V alleles also contain coding stops. | Rejected: 4 curatable validation issues. |

The explicit `allow_curatable_refdata()` override permits legacy compilation.
For each cartridge, this exact legacy profile ran seeds
`20261020..20261049` successfully with `strict=True`:

```python
exp = (
    ga.Experiment.on(name).allow_curatable_refdata()
    .recombine(np1_lengths=[(n, 1) for n in range(7)])
    .trim(enabled=False).productive_only()
)
compiled = exp.compile()
record = compiled.simulator.run(
    seed=20261020, strict=True
).airr_record(compiled.refdata)
```

Both profiles reported `productive=True` for all 30 records. All 30
`mouse_tcra` junction translations lacked the initial C; one also lacked
terminal F/W. All 30 `mouse_tcrg` junction translations lacked terminal F/W.
Concrete seed `20261020`:

- `mouse_tcra`: `TRAV14D-3/DV8*03` / `TRAJ19*01`; reported junction
  translation begins `QQQVRQSP` and ends `ASVVI`. The junction incorrectly
  spans almost the entire V segment.
- `mouse_tcrg`: `TRGV2*01` / `TRGJ4*01`; reported junction translation
  `CAVWMSRS`, ending in S.

The legacy anchor contract preserves the reference anchor amino acid,
without requiring that the reference amino acid itself be appropriate
([retention.rs](../engine_rs/src/contract/anchor_preserved/retention.rs)).
Record productivity likewise uses preservation of the annotated reference
anchor ([builder.rs](../engine_rs/src/airr_record/builder.rs)).
Consequently, metadata correctness under this override depends on valid
reference annotations; successful legacy compilation is insufficient evidence.

This shared record-classification limitation also affects **unconstrained
new-model output** under the same explicit override. Replacing
`.trim(enabled=False).productive_only()` in the example above with
`.process_coding_ends(reference_ends_are_coding_ends=True,
v_3=[(0, 1)], j_5=[(0, 1)])`, without `productive_only()`, gives:

| Cartridge | Seeds | Records labeled productive | Labeled productive with noncanonical junction endpoints |
| --- | --- | --- | --- |
| `mouse_tcra` | `20261020..20261049` | 9/30 | 9/9 lack initial C. |
| `mouse_tcrg` | `20261020..20261049` | 1/30 | 1/1 lacks terminal F/W. |

Examples are `mouse_tcra` seed `20261022`, `TRAV10D*03` / `TRAJ34*01`,
whose reported junction translation starts `GEKVEQHE`, and `mouse_tcrg`
seed `20261024`, `TRGV7*01` / `TRGJ4*01`, whose reported junction translation
is `CASWAGLS`. Both have `productive=True`. Successful original-end ancestry,
coordinate and replay checks do not establish complete biological correctness
when invalid reference annotations are explicitly admitted.

The new original-end productive solver correctly rejects both cartridges
with “no productive completion,” including with all-zero outcomes. The sweep's
V outcomes `[-2,-1,0,1,2]` and J outcomes `[-1,0,1]` cannot repair those anchor
annotations. As a diagnostic control only, an in-memory copy restricted to
`TRAV1*01` / `TRAJ11*01`, with the V anchor changed from 0 to its terminal
in-frame C codon at 261 and J anchor unchanged at 25, produced 30/30 productive
new-model records with C…F junction translations. This is not a curated
replacement cartridge.

Temporary reproducers: `/tmp/genairr_mouse_tcra_diagnostic.py`,
`/tmp/genairr_mouse_tcrg_diagnostic.py`,
`/tmp/genairr_mouse_tcra_control.py`. Reference annotations need independent
curation; legacy acceptance under the explicit override remains a separate
known limitation.

## Strict productive indel sampling can exhaust its retry budget

For `human_igh`, use N1/N2 lengths `0..6` with equal weights, D inversion
probability 0.5, `productive_only()`, uniform SHM rate 0.12, and
`polymerase_indels(count=2)`. In original-end mode, use equal-weight outcomes:
V3 `[-2,-1,0,1,2]`, D5 `[-1,0,1]`, D3 `[-2,0,2]`, J5 `[-1,0,1]`.
Strict runs fail at seeds `20261028`, `20261030`, and `20261031` with:

```text
('corrupt.indel', 'corrupt.indel.site[1]', 'empty_admissible_support')
```

Legacy end processing with the same N, inversion, mutation and indel settings
also fails at seed `20261027`. This is not specific to original-end processing.

The tuple sampler retries at most 16 times, then reports the above error;
the retry limit and failure mapping already exist at the baseline commit.
Its DP filters junction-frame deltas; the final contract check additionally
tests complete coding validity. Exhaustion does **not** establish empty
biological support
([tuple_sampler.rs](../engine_rs/src/passes/corrupt/indel/tuple_sampler.rs)).

On the exact 370-nt precursor from failing original-end seed `20261028`,
running only the two-indel pass with the same productive contracts succeeded
for 38/40 independent seeds. Seed 0 inserts C at query position 113, then
deletes position 175 in the updated sequence. The failing legacy precursor
likewise succeeded for 35/40 independent seeds. Both comparisons establish
nonempty support. Permissive reruns of the original failing seeds retain the
sequence and record two `site=-1` no-ops.

Temporary reproducers: `/tmp/genairr_indel_diagnostic.py`,
`/tmp/genairr_indel_fixed_precursor.py`. A productive stress run that omits
indels does not validate this sampler limitation.
