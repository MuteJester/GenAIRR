"""Fixtures shared by several test modules."""
from __future__ import annotations

from types import SimpleNamespace

import pytest


@pytest.fixture(scope="session")
def igh_starter_session(tmp_path_factory):
    """The human IGH release 10 starter, built once per test session.

    The pinned exports are recorded as local files (``OgrdbSource.from_files``) and the selection
    fixture's recipe is projected over them once.  Every part is immutable, so tests share it
    without affecting each other.  Fields: ``source``, ``lock``, ``blobs`` (format -> bytes),
    ``inventory`` (interpreted), ``document`` (the recipe as a dict), ``projection`` (the
    ``Projection``) and ``build`` (the ``ReferenceBuild`` of that projection).
    """
    from GenAIRR._reference_sources import model, projection
    from GenAIRR._reference_sources.build import ReferenceBuild
    from tests.test_reference_source_api import ALL, IGH, fixture_bytes, from_files
    from tests.test_reference_source_projection import interpreted, starter_recipe
    source = from_files(tmp_path_factory.mktemp("igh-starter"), IGH)
    blobs = {fmt: fixture_bytes(IGH, fmt) for fmt in ALL}
    inventory = interpreted(blobs, source.lock)
    document = starter_recipe(source.lock, inventory)
    recipe = model.Recipe.from_dict(document)
    result = projection.project(source.lock, blobs, recipe)
    build = ReferenceBuild._create(lock=source.lock, recipe=recipe, assets=blobs, inputs=result.inputs,
                                   reference=result.reference, report=result.report)
    return SimpleNamespace(source=source, lock=source.lock, blobs=blobs, inventory=inventory, document=document,
                           projection=result, build=build)


@pytest.fixture(scope="session")
def igl_two_allele_session(tmp_path_factory):
    """``(source, build)``: two IGL alleles built once per session over the complete pinned IGL
    release 4 exports (recorded as local files; no scope patch).  Both are immutable."""
    from tests.test_reference_source_api import from_files, full_recipe
    root = tmp_path_factory.mktemp("igl-two-alleles")
    source = from_files(root)
    return source, source.build(full_recipe(root, source))
