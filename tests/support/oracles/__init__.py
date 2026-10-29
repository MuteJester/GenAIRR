"""Derive expected truth independently of the code under test.

Three rules. The first two are enforced mechanically by tests in support/meta;
the third is enforced as far as it can be, and the limit is stated below.

- An oracle must not import GenAIRR. It is handed a recorded result and reads
  it; it never calls the library whose output it is checking. A junction oracle
  that calls the production junction builder proves nothing.
- Truth comes from the engine's recorded per-base provenance, never from
  re-aligning sequences. Re-alignment introduces a second, unvalidated
  inference step and is not a source of truth.
- **Every identity must be probed for circularity before it is claimed.**
  Perturb *both* sides of the comparison together -- the coordinate and the
  byte, the string and the index it is read at -- and measure what the oracle
  reports. If it reports nothing, the identity is a tautology about one
  producer and the claim must say so.

  This rule exists because the same over-claim was made twice. The
  reconstruction oracle compares a recorded germline base to the reference at
  the coordinate the record names, and the engine lays byte and coordinate
  down from the same two numbers, so shifting both agrees: zero mismatches
  over 366 positions. That was found, written down, and then the projection
  oracle asserted the opposite of it -- 226 columns carrying a different
  germline base produced one complaint. A review caught it both times, after
  the claims were already committed.

  Each reporting oracle therefore declares a module-level ``CIRCULARITY``
  mapping: what was perturbed, and what the oracle said. **Its limit:** the
  keys are the author's own names for the identities and nothing checks that
  they cover every identity the module can emit. What the meta test enforces
  is that the mapping exists, that every note carries a measurement, and that
  a module which reports disagreements has one at all. The judgement that the
  list is complete is still the author's, which is exactly the judgement that
  was wrong twice -- so a reviewer should check it rather than trust it.

Oracles are deliberately slow and obvious. They are the thing that has to be
right, so clarity beats cleverness.
"""
