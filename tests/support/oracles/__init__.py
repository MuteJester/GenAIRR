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
  mapping: what was perturbed, and what the oracle said. Coverage is enforced
  rather than trusted -- every module also declares ``IDENTITIES``, the
  mapping must cover exactly those keys, and four of the five validate the
  identity at construction so a report cannot name one nobody declared. The
  fifth, ``reconstruction``, has no identity-carrying result type and
  attributes every public function instead. The old limit here read "nothing
  checks that they cover every identity the module can emit", and under it two
  modules were found emitting identities absent from their mapping and a third
  declaring none at all.

  **The limit that remains:** a note is prose, and prose does not re-measure
  itself. tests/biology/test_oracle_circularity.py executes the probes for the
  identities with no sensitivity coverage otherwise and pins the rest as
  prose-only with what each would take, so the unexecuted count can only go
  down. Executing the first five corrected one note whose numbers came from a
  different fixture than its own sentence claimed.

Oracles are deliberately slow and obvious. They are the thing that has to be
right, so clarity beats cleverness.
"""


def reporting_modules() -> dict[str, "object"]:
    """``name -> module`` for every oracle here that reports disagreements.

    Derived by scanning for a result type rather than listed, so a new oracle
    is subject to the rules above the moment it is written. Lives in the
    package rather than in the meta test that first needed it, because the
    circularity probes need the same scan and a test module may not import
    another test module.
    """
    import importlib

    from tests.support import paths

    found = {}
    root = paths.tests_root() / "support" / "oracles"
    for module in sorted(root.glob("*.py")):
        if module.name == "__init__.py":
            continue
        source = module.read_text(encoding="utf-8")
        if "class Disagreement" not in source and "class Mismatch" not in source:
            continue
        found[module.stem] = importlib.import_module(f"{__name__}.{module.stem}")
    return found
