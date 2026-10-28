"""Hypothesis strategies for the suite.

Oracle verification is naturally property-based: the claim is about any
simulated sequence, not three hand-picked ones. Strategies live here so a
property test states what it is asserting rather than how to build inputs.
"""

import pytest

#: Hypothesis is a test-only extra. Skipping here rather than raising keeps a
#: checkout without the extra collectable: conftest guards its own use of
#: hypothesis for the same reason, and an unguarded import here would have
#: made the whole run abort instead of skipping these tests.
pytest.importorskip("hypothesis")
