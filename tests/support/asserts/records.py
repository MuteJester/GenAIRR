"""Assertions about AIRR records."""
from __future__ import annotations

from typing import Any, Mapping


def assert_field(record: Mapping[str, Any], field: str, expected: Any, *, because: str) -> None:
    """Assert one record field, reporting the claim that field is supposed to satisfy.

    ``because`` states what the field means, not what the test does: a failure
    should say which guarantee broke, not merely that two values differed.
    """
    actual = record.get(field)
    assert actual == expected, (
        f"record field {field!r} is {actual!r}, expected {expected!r}\n"
        f"  because: {because}")
