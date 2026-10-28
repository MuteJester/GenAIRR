"""Fixtures shared by several test modules."""
from __future__ import annotations

from dataclasses import dataclass
import json
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
    from tests.test_reference_source_api import ALL, IGH, fixture_bytes, from_files
    from tests.test_reference_source_projection import interpreted, starter_recipe
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
    from tests.test_reference_source_api import from_files, full_recipe
    root = tmp_path_factory.mktemp("igl-two-alleles")
    source = from_files(root)
    return source, source.build(full_recipe(root, source))
