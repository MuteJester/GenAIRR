"""Errors raised by the portable reference artifact code."""
from __future__ import annotations

KINDS = frozenset({
    "malformed", "unsupported_convention", "unsupported_schema", "conflict",
    "integrity", "capability", "empty_selection", "io",
})


class ReferenceArtifactError(ValueError):
    """A reference artifact failure; ``kind`` classifies it and ``path`` locates it."""

    def __init__(self, kind: str, path: str, message: str):
        if kind not in KINDS:
            raise ValueError(f"unknown reference error kind {kind!r}")
        self.kind = kind
        self.path = path
        self.message = message
        super().__init__(f"{path}: {message}" if path else message)


class ReferenceCapabilityError(ReferenceArtifactError):
    """A requested capability is not supported by the reference's claims."""
