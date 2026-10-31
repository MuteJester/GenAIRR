"""Fixtures shared by several test modules."""
from __future__ import annotations

from dataclasses import dataclass
import importlib.util
import json
import os
from types import MappingProxyType
from typing import Any, Mapping

import pytest


@dataclass(frozen=True)
class StarterSession:
    """The human IGH starter built once for the session.  Nothing in it can be changed: the
    attributes are frozen, ``blobs`` is a read-only mapping of bytes, every model object is
    immutable, and ``document`` returns a fresh copy of the recipe on every access."""

    source: Any                 # OgrdbSource
    lock: Any                   # SourceLock
    blobs: Mapping[str, bytes]
    inventory: Any              # SourceInventory, interpreted
    projection: Any             # Projection
    build: Any                  # ReferenceBuild of that projection
    recipe_text: str            # the recipe document as canonical JSON

    @property
    def document(self) -> dict:
        """A fresh copy of the recipe document."""
        return json.loads(self.recipe_text)


@pytest.fixture(scope="session")
def igh_starter_session(tmp_path_factory) -> StarterSession:
    """The human IGH release 10 starter, built once per test session: the pinned exports recorded
    as local files (``OgrdbSource.from_files``) and the selection fixture's recipe projected over
    them once."""
    from GenAIRR._artifact_io import canonical_json
    from GenAIRR._reference_sources import model, projection
    from GenAIRR._reference_sources.build import ReferenceBuild
    from tests.support.fixtures.source_api import ALL, IGH, fixture_bytes, from_files
    from tests.support.fixtures.source_projection import interpreted, starter_recipe
    source = from_files(tmp_path_factory.mktemp("igh-starter"), IGH)
    blobs = MappingProxyType({fmt: fixture_bytes(IGH, fmt) for fmt in ALL})
    inventory = interpreted(blobs, source.lock)
    document = starter_recipe(source.lock, inventory)
    recipe = model.Recipe.from_dict(document)
    result = projection.project(source.lock, blobs, recipe)
    build = ReferenceBuild._create(lock=source.lock, recipe=recipe, assets=blobs, inputs=result.inputs,
                                   reference=result.reference, report=result.report)
    return StarterSession(source=source, lock=source.lock, blobs=blobs, inventory=inventory, projection=result,
                          build=build, recipe_text=canonical_json(document))


@pytest.fixture(scope="session")
def igl_two_allele_session(tmp_path_factory) -> tuple:
    """``(source, build)``: two IGL alleles built once per session over the complete pinned IGL
    release 4 exports (recorded as local files; no scope patch).  Both are immutable."""
    from tests.support.fixtures.source_api import from_files, full_recipe
    root = tmp_path_factory.mktemp("igl-two-alleles")
    source = from_files(root)
    return source, source.build(full_recipe(root, source))


# ── Hypothesis profiles ──────────────────────────────────────────────────────
# Registered centrally so property-based tests get workable limits without
# repeating settings at each call site. Guarded on the import: the suite must
# still collect when the optional `test` extra is not installed.
try:
    from hypothesis import HealthCheck, settings
except ImportError:  # the extra is optional
    pass
else:
    _COMMON = {
        # Building a simulated repertoire is orders of magnitude slower than
        # Hypothesis's default deadline, which would flake rather than inform.
        "deadline": None,
        "suppress_health_check": [HealthCheck.too_slow],
    }
    settings.register_profile("dev", max_examples=25, **_COMMON)
    # derandomize so a CI failure reproduces from the commit alone, with no
    # seed to hunt for in the log. `database=None` because CI checks out a
    # fresh tree and so never has one: with the database on, running this
    # profile locally replays stored failures from earlier runs and fails a
    # clean tree on a defect already fixed, which happened and cost an hour.
    settings.register_profile("ci", max_examples=100, derandomize=True,
                              database=None, **_COMMON)
    # `deep` is the discovery campaign's profile: many examples, and no
    # example database, so a run carries no state from the last one. NOT
    # derandomized, unlike `ci` -- a derandomized campaign explores the same
    # examples for ever, which is the opposite of what a campaign is for. The
    # seed comes from the caller instead (CI derives it from the commit), so
    # every commit covers new ground and every failure reproduces from the
    # seed in the log.
    # The discovery campaign's profiles. NOT derandomized, unlike `ci`: a
    # derandomized campaign explores the same examples for ever, which is the
    # opposite of what a campaign is for. The seed comes from the caller
    # instead (CI derives it from the commit), so every commit covers new
    # ground and every failure reproduces from the seed in the log. No example
    # database either, so a run carries no state from the last one.
    #
    # `explain` is off in both. It re-runs a failing example dozens of times
    # to attribute it to lines of code, which cost 88 seconds of a 92-second
    # run and tells us nothing: the campaign prints the pipeline and the seed
    # itself. Shrinking is capped the same way -- on a property where every
    # example builds a repertoire, an unbounded shrink spent five minutes and
    # gave up, which is 16 minutes of a pull request spent on presentation.
    from hypothesis import Phase

    _CAMPAIGN_PHASES = (Phase.explicit, Phase.reuse, Phase.generate,
                        Phase.target, Phase.shrink)
    settings.register_profile(
        "campaign", max_examples=200, database=None, phases=_CAMPAIGN_PHASES,
        report_multiple_bugs=False, **_COMMON)
    settings.register_profile(
        "deep", max_examples=600, database=None, phases=_CAMPAIGN_PHASES,
        report_multiple_bugs=False, **_COMMON)
    # `scout` is for working a finding by hand: generate only, so the first
    # failing example arrives in seconds instead of after a shrink.
    settings.register_profile(
        "scout", max_examples=300, database=None,
        phases=(Phase.explicit, Phase.reuse, Phase.generate, Phase.target),
        report_multiple_bugs=False, **_COMMON)
    settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "dev"))


# ── biology coverage ─────────────────────────────────────────────────────────
# Hooks imported rather than defined here: the implementation lives with the
# catalog it reports on.
from tests.biology.catalog._plugin import (  # noqa: E402,F401
    pytest_addoption,
    pytest_collection_modifyitems,
    pytest_runtest_logreport,
    pytest_sessionfinish,
    pytest_terminal_summary,
)

# pytest_testnodedown is declared by pytest-xdist, so implementing it without
# xdist installed is an unknown hook and a plugin validation error. It is only
# needed when the run is distributed, which is exactly when xdist is present.
if importlib.util.find_spec("xdist") is not None:  # noqa: E402
    from tests.biology.catalog._plugin import pytest_testnodedown  # noqa: E402,F401
