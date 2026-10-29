"""Derive expected truth independently of the code under test.

Two rules, both enforced by tests in support/meta:

- An oracle must not import GenAIRR. It is handed a recorded result and reads
  it; it never calls the library whose output it is checking. A junction oracle
  that calls the production junction builder proves nothing.
- Truth comes from the engine's recorded per-base provenance, never from
  re-aligning sequences. Re-alignment introduces a second, unvalidated
  inference step and is not a source of truth.

Oracles are deliberately slow and obvious. They are the thing that has to be
right, so clarity beats cleverness.
"""
