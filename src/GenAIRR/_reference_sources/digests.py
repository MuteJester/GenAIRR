"""The typed digest graph of a published-source build.

Raw assets are identified by ordinary SHA-256 over their bytes.  Every canonical model
digest uses its own domain, and the dependency graph is acyclic::

    raw asset hashes -> lock hash -> inventory hash
    lock hash + recipe -> recipe hash
    lock + inventory + recipe + adapter/projection ids -> build-input hash
    build-input hash -> generated FASTA/sidecar -> runtime reference hashes
    inventory + decisions + projection + runtime hashes -> report hash
    member hashes + report hash + runtime hashes + actual builder identity -> manifest hash

Acquisition receipts are audit events and enter no digest.  Runtime citations may name
the build-input hash, never the report or manifest hash.
"""
from __future__ import annotations

from GenAIRR._reference_sources._fields import digest
from GenAIRR._reference_sources.lock import DIGEST_DOMAIN as LOCK_DOMAIN, SourceLock

DOMAINS = (
    LOCK_DOMAIN,
    "genairr_reference_source_inventory_v1",
    "genairr_reference_source_recipe_v1",
    "genairr_reference_source_build_inputs_v1",
    "genairr_reference_source_build_report_v1",
    "genairr_reference_source_build_manifest_v1",
)


def _require(value, cls, name):
    if not isinstance(value, cls):
        raise TypeError(f"{name} must be a {cls.__name__}, got {type(value).__name__}")


def lock_hash(lock: SourceLock) -> str:
    """``genairr_reference_source_lock_v1`` over every lock field but the commitment (= the commitment)."""
    _require(lock, SourceLock, "lock")
    return digest(lock.body(), LOCK_DOMAIN)
