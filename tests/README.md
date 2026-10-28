# The GenAIRR test suite

GenAIRR is used to benchmark other tools, so everything downstream treats its
metadata as a perfect oracle: when it says a read's `np1` is five untemplated
bases, every comparison built on it assumes that is what actually happened. The
suite's job is therefore not "does the code run" but "is the recorded truth
actually true."

## What to run

```bash
make test-fast        # the inner loop: drops the slowest modules, runs parallel
make test-parallel    # everything except the wall-time budgets
make test-budgets     # the budgets, serially and alone
make test-lint        # ruff over tests/
```

`test-fast` is *parallel minus the heavy tail*, not a serial fast subset. The
`slow` marker alone accounts for about half the runtime; most of the speed comes
from `pytest-xdist`. Serially the same selection takes several minutes.

The budgets run on their own deliberately. Under `-n auto` the other workers
saturate the cores, so the measured time reflects contention rather than the
engine, and the budgets fail intermittently.

Install what the suite needs with `pip install -e ".[test,mcp]"`. Without the
`mcp` extra the MCP server tests are skipped; without `airr` the schema
compliance tests are skipped.

## Where a new test goes

```
tests/
  components/     by what part of GenAIRR it covers
  biology/        the mechanism catalog, and tests of biological truth
  acceptance/     end-to-end journeys a user actually takes
  performance/    wall-time budgets
  support/        shared machinery; not a test package
  data/ golden/   fixture data and reproducibility references
```

`components/` is named for functional areas rather than import paths, because
most of GenAIRR lives in top-level modules and a literal package mirror would
have nowhere to put a large share of these tests. Put a test next to the area it
covers; if two areas could claim it, prefer the one whose code you would edit to
make it fail.

A test of **biological truth** — "is this field what actually happened?" — goes
under `biology/`, organised by developmental stage rather than by code, because
a biological claim cuts across modules: the recorded junction depends on
trimming, P addition, N addition and coding-end processing at once.

## Declaring a tier

Every test module outside `support/` declares one — `components/`, `biology/`,
`acceptance/` and `performance/` alike:

```python
pytestmark = pytest.mark.behavior
```

- `pin` — pins a live surface a feature reuses; fails when that surface shifts
- `behavior` — exercises a feature's own behaviour
- `composition` — runs the feature against the full release stack

`-m pin` then selects across the whole suite. Module-level `skip` or `xfail` is
refused; a conditional `skipif` must be recorded in
`support/meta/test_structure.py` with its reason, so an area that never runs is
visible rather than invisible.

## Registering a biological mechanism

Coverage is claimed by registration, not by location, so a test counts wherever
it lives:

```python
@verifies("junction.n_addition",
          claim="np1 is exactly the bases carrying no germline coordinate",
          fields=("np1", "np1_length"))
def test_np1_matches_the_provenance():
    ...
```

`fields` is required and validated. It names the truth fields this test actually
checks, because the report stands on it: a test checking one of a mechanism's
six fields is a partial result, and defaulting to "all of them" would let one
test report the mechanism verified. An unknown mechanism id or a field the
mechanism does not own fails at collection.

```bash
pytest --biology-coverage        # verified / partial / unverified per mechanism
```

Status comes from registered Python tests only. A registration on an `xfail`ed
or skipped test does not count — crediting a mechanism from a test that cannot
pass is no better than not testing it.

Mechanisms are declared in `biology/catalog/`. Adding an engine pass without
classifying it there fails `test_every_engine_pass_is_classified_exactly_once`.

## Writing an oracle

An oracle computes what the answer *should* be without using the code path it
checks. Two rules, both enforced by tests:

- Oracles import nothing from `GenAIRR`. They are handed a recorded result and
  read it. A junction oracle that called the production junction builder would
  prove nothing.
- Truth comes from recorded per-base provenance or from the reference tables,
  never from re-aligning sequences. Reading the allele a record *names* and
  looking up the base at the coordinate it *claims* is the independent step.

Two things to know before working on junction truth. A **P nucleotide is
templated** — the reverse complement of its own segment's post-trim tip — so
"no germline coordinate" does not mean "untemplated". And **P versus N survives
only in letter case**, lowercase for P: calling `.upper()` destroys the
distinction.

## Shared machinery

```
support/
  fixtures/     pytest fixtures and session-scoped builds
  builders/     construct inputs; imports no pytest
  oracles/      derive expected truth independently
  asserts/      domain assertions with legible failure messages
  cases/        shared parametrization tables
  strategies/   Hypothesis strategies
  golden.py     the reproducibility fixture definition
  paths.py      resolve tests/data and tests/golden
  meta/         the suite's own structural tests
```

Two rules, both enforced:

- **No test module imports another test module.** Shared setup goes in
  `support/`. Borrowing a fixture from a neighbour couples their collection
  order and makes the lender impossible to move.
- **Do not build a path from `__file__`.** Use `support.paths`, which is anchored
  to a fixed location. File-relative arithmetic is what makes a module unsafe to
  relocate, and it broke silently when this tree was reorganised — in one case
  as an empty `parametrize` set that collected zero tests and reported nothing.

  Three uses are legitimate and are listed, with their reasons, in
  `FILE_ANCHOR_ALLOWED` in `support/meta/test_structure.py`: `paths.py` itself,
  the two guards that audit the support layer and so cannot be positioned by
  it, and a pin that reads its own source to check its own docstring. Anything
  else fails. The list is the rule — a new exception is argued for there rather
  than added by precedent.

## Reproducibility fixtures

`tests/golden/seed42_*.tsv` are the cross-platform references. Regenerate them
only when the output format changes on purpose:

```bash
python tools/regenerate_golden_files.py
```

Then commit them: they are what every CI platform is checked against.
