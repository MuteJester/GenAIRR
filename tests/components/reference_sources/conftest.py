"""Fixtures for the reference-source tests.

These live in a conftest rather than being imported by name, and the reason is
concrete rather than stylistic. The original design preferred explicit imports
so a reader could see where a fixture came from, but a fixture imported by name
has two costs:

- It is syntactically indistinguishable from an unused import, so ruff's F401
  cannot be enforced on any module that uses one. Running `ruff --fix` over the
  suite removed these imports and broke 177 tests.
- A fixture's own dependencies must be imported too. ``inventory`` takes
  ``inv``, and ``report_d`` takes both; importing ``inventory`` alone fails at
  run time with "fixture 'inv' not found", which is how this surfaced.

pytest resolves the whole chain from here. The implementations stay in
tests/support/fixtures/ so there is still one definitional home; this file only
makes them visible to the directory that uses them.
"""
from __future__ import annotations

from tests.support.fixtures.reference_source_models import (  # noqa: F401
    inv,
    inventory,
    report,
    report_d,
)
from tests.support.fixtures.source_api import offline  # noqa: F401
