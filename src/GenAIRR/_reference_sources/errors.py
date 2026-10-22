"""Typed errors raised while importing published germline sources."""
from __future__ import annotations

KINDS = (
    "io", "network", "resource_limit", "malformed", "unsupported_schema", "unsupported_profile",
    "identity_mismatch", "digest_mismatch", "missing_asset", "incomplete_dossier",
    "unresolved_selection", "stale_decision", "projection", "capability",
)


class ReferenceSourceError(ValueError):
    """A published-source failure.

    ``kind`` is one of :data:`KINDS`; ``path`` is a JSON-pointer-like location inside the
    document being read ("" for the whole document); ``message`` names the affected
    asset, record or allele and the expected and observed values where they apply.
    """

    def __init__(self, kind: str, path: str, message: str):
        if kind not in KINDS:
            raise ValueError(f"unknown reference source error kind {kind!r}")
        self.kind = kind
        self.path = path
        self.message = message
        super().__init__(f"{path}: {message}" if path else message)

    def under(self, prefix: str) -> "ReferenceSourceError":
        """The same error located below ``prefix`` in an enclosing document."""
        return ReferenceSourceError(self.kind, prefix + self.path, self.message)
