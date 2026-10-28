"""The engine's own list of the passes it defines.

A mechanism catalog that claims to cover the simulation needs to know what the
simulation consists of. ``Outcome.pass_names()`` reports only the passes a
given compiled plan scheduled, so a pass nobody schedules is invisible -- which
is exactly the drift a completeness check exists to catch.
"""
from __future__ import annotations

import pytest

from GenAIRR._compiled import pass_registry

pytestmark = pytest.mark.pin


def test_registry_lists_the_passes_biology_depends_on():
    names = set(pass_registry())
    expected = {"TrimPass", "GenerateNPPass", "PAdditionPass", "AssembleSegmentPass",
                "InvertDPass", "ReceptorRevisionPass", "SampleAllelePass",
                "ProcessCodingEndsPass", "EmitCodingEndPPass"}
    assert expected <= names, f"missing from the registry: {sorted(expected - names)}"


def test_registry_lists_the_assay_passes():
    names = set(pass_registry())
    expected = {"PCRErrorPass", "QualityErrorPass", "NCorruptionPass", "ContaminantPass",
                "RevCompPass", "EndLossPass", "IndelPass", "PairedEndSamplingPass"}
    assert expected <= names, f"missing from the registry: {sorted(expected - names)}"


def test_registry_excludes_test_doubles():
    """The doubles live behind cfg(test) and must not reach a released build."""
    doubles = {"FooPass", "StubPass", "SpoofedPass", "TestEchoPass", "TestRandomBasePass",
               "NamedPass", "TraceProbePass"}
    assert not (doubles & set(pass_registry()))


def test_registry_is_sorted_and_has_no_duplicates():
    names = pass_registry()
    assert list(names) == sorted(set(names))


def test_registry_is_a_tuple_of_str():
    names = pass_registry()
    assert isinstance(names, tuple) and all(isinstance(n, str) for n in names)


def _declared_passes():
    """Every pass type the engine defines, read from its own source.

    ``impl Pass for X`` at column zero under ``engine_rs/src/passes`` is the
    definition of a production pass. Test doubles are either in a ``tests``
    module, in ``test_support.rs``, or indented inside a ``#[cfg(test)]``
    block, and all three are excluded.
    """
    import re

    from tests.support import paths

    root = paths.repo_root() / "engine_rs" / "src" / "passes"
    pattern = re.compile(r"^impl Pass for ([A-Za-z0-9_]+)", re.MULTILINE)
    found = {}
    for source in sorted(root.rglob("*.rs")):
        parts = source.parts
        if source.name in ("test_support.rs", "tests.rs") or "tests" in parts:
            continue
        for name in pattern.findall(source.read_text(encoding="utf-8")):
            found[name] = str(source.relative_to(paths.repo_root()))
    return found


def test_the_registry_names_every_pass_the_engine_defines():
    """The registry is a hand-maintained list, so the dangerous direction is the
    one it cannot check itself.

    ``registry_names_types_the_crate_exports`` on the Rust side proves every
    name listed is a real type. Nothing proved the converse: a pass added to
    the crate and forgotten here is invisible to the catalog's completeness
    test, which would stay green while a whole mechanism went unclassified --
    the precise drift that test exists to catch.
    """
    declared = _declared_passes()
    registered = set(pass_registry())
    assert len(declared) >= 20, (
        f"only {len(declared)} pass implementations found; the scan is not "
        "reading the engine source")
    missing = {name: where for name, where in declared.items() if name not in registered}
    assert missing == {}, (
        "passes the engine defines but the registry does not name:\n"
        + "\n".join(f"  {name} ({where})" for name, where in sorted(missing.items()))
        + "\n\nAdd them to PASS_NAMES in engine_rs/src/python/registry.rs and "
          "classify them in the mechanism catalog.")
    assert registered - set(declared) == set(), (
        f"the registry names passes the engine does not define: "
        f"{sorted(registered - set(declared))}")
