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

from GenAIRR._reference_sources._fields import digest, digest_of_canonical, require_type
from GenAIRR._reference_sources.consistency import check_build_inputs
from GenAIRR._reference_sources.inventory import SourceInventory
from GenAIRR._reference_sources.lock import DIGEST_DOMAIN as LOCK_DOMAIN, SourceLock
from GenAIRR._reference_sources.recipe import Recipe

INVENTORY_DOMAIN = "genairr_reference_source_inventory_v1"
RECIPE_DOMAIN = "genairr_reference_source_recipe_v1"
BUILD_INPUTS_DOMAIN = "genairr_reference_source_build_inputs_v1"
REPORT_DOMAIN = "genairr_reference_source_build_report_v1"
MANIFEST_DOMAIN = "genairr_reference_source_build_manifest_v1"

DOMAINS = (
    LOCK_DOMAIN,
    INVENTORY_DOMAIN,
    RECIPE_DOMAIN,
    BUILD_INPUTS_DOMAIN,
    REPORT_DOMAIN,
    MANIFEST_DOMAIN,
)


def lock_hash(lock: SourceLock) -> str:
    """``genairr_reference_source_lock_v1`` over every lock field but the commitment (= the commitment)."""
    require_type(lock, SourceLock, "lock")
    return digest(lock.body(), LOCK_DOMAIN)


def recipe_hash(recipe: Recipe) -> str:
    """``genairr_reference_source_recipe_v1`` over the complete recipe (which holds the lock hash)."""
    require_type(recipe, Recipe, "recipe")
    return digest(recipe.to_dict(), RECIPE_DOMAIN)


def _memoised(document, name: str, compute) -> str:
    """A digest of a frozen document, computed once per instance (the document never changes)."""
    cached = document.__dict__.get(name)
    if cached is None:
        cached = compute()
        object.__setattr__(document, name, cached)
    return cached


def inventory_hash(inventory: SourceInventory) -> str:
    """``genairr_reference_source_inventory_v1`` over the complete inventory (which holds the lock hash)."""
    require_type(inventory, SourceInventory, "inventory")
    return _memoised(inventory, "_inventory_hash", lambda: digest(inventory.to_dict(), INVENTORY_DOMAIN))


def build_inputs(lock: SourceLock, inventory: SourceInventory, recipe: Recipe) -> dict:
    """The build-input record: lock, inventory and recipe hashes plus adapter and projection ids.

    Runs every cross-document check first, so inconsistent inputs never get an identity.
    """
    require_type(lock, SourceLock, "lock")
    require_type(inventory, SourceInventory, "inventory")
    require_type(recipe, Recipe, "recipe")
    check_build_inputs(lock, inventory, recipe)
    return unchecked_build_inputs(lock, inventory, recipe)


def unchecked_build_inputs(lock: SourceLock, inventory: SourceInventory, recipe: Recipe) -> dict:
    """The build-input record without re-running the cross-document checks (callers ran them)."""
    return {"source_lock": lock_hash(lock), "inventory": inventory_hash(inventory), "recipe": recipe_hash(recipe),
            "adapter": recipe.adapter.to_dict(), "projection": recipe.projection.to_dict()}


def build_input_hash(lock: SourceLock, inventory: SourceInventory, recipe: Recipe) -> str:
    """``genairr_reference_source_build_inputs_v1`` over :func:`build_inputs`.

    Exists before projection, so runtime citations may name it; it cannot depend on a
    report, manifest or receipt.
    """
    return digest(build_inputs(lock, inventory, recipe), BUILD_INPUTS_DOMAIN)


def report_hash(report) -> str:
    """``genairr_reference_source_build_report_v1`` over the complete report (never stored in it)."""
    from GenAIRR._reference_sources.report import BuildReport
    require_type(report, BuildReport, "report")
    return _memoised(report, "_report_hash", lambda: digest(report.to_dict(), REPORT_DOMAIN))


def remember_report_hash(report, member: bytes) -> str:
    """:func:`report_hash` from the report's member bytes (its canonical JSON and one LF), without
    serialising the report again, and remembered for ``report``.

    The caller guarantees that ``member`` is exactly that serialisation of ``report``; a dossier
    check has just compared the two byte for byte.
    """
    from GenAIRR._reference_sources.report import BuildReport
    require_type(report, BuildReport, "report")
    if type(member) is not bytes or not member.endswith(b"\n"):
        raise ValueError("expected the report's member bytes: canonical JSON followed by one LF")
    text = member[:-1].decode("utf-8")
    return _memoised(report, "_report_hash", lambda: digest_of_canonical(text, REPORT_DOMAIN))


def manifest_hash(manifest) -> str:
    """``genairr_reference_source_build_manifest_v1`` over the complete manifest (never stored in it)."""
    from GenAIRR._reference_sources.manifest import BuildManifest
    require_type(manifest, BuildManifest, "manifest")
    return digest(manifest.to_dict(), MANIFEST_DOMAIN)
