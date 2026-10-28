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
