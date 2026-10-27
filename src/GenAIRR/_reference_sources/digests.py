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


def _require(value, cls, name):
    if not isinstance(value, cls):
        raise TypeError(f"{name} must be a {cls.__name__}, got {type(value).__name__}")


def lock_hash(lock: SourceLock) -> str:
    """``genairr_reference_source_lock_v1`` over every lock field but the commitment (= the commitment)."""
    _require(lock, SourceLock, "lock")
    return digest(lock.body(), LOCK_DOMAIN)


def recipe_hash(recipe: Recipe) -> str:
    """``genairr_reference_source_recipe_v1`` over the complete recipe (which holds the lock hash)."""
    _require(recipe, Recipe, "recipe")
    return digest(recipe.to_dict(), RECIPE_DOMAIN)


def inventory_hash(inventory: SourceInventory) -> str:
    """``genairr_reference_source_inventory_v1`` over the complete inventory (which holds the lock hash)."""
    _require(inventory, SourceInventory, "inventory")
    return digest(inventory.to_dict(), INVENTORY_DOMAIN)


def build_inputs(lock: SourceLock, inventory: SourceInventory, recipe: Recipe) -> dict:
    """The build-input record: lock, inventory and recipe hashes plus adapter and projection ids.

    Runs every cross-document check first, so inconsistent inputs never get an identity.
    """
    _require(lock, SourceLock, "lock")
    _require(inventory, SourceInventory, "inventory")
    _require(recipe, Recipe, "recipe")
    check_build_inputs(lock, inventory, recipe)
    return {"source_lock": lock_hash(lock), "inventory": inventory_hash(inventory), "recipe": recipe_hash(recipe),
            "adapter": recipe.adapter.to_dict(), "projection": recipe.projection.to_dict()}


def build_input_hash(lock: SourceLock, inventory: SourceInventory, recipe: Recipe) -> str:
    """``genairr_reference_source_build_inputs_v1`` over :func:`build_inputs`.

    Exists before projection, so runtime citations may name it; it cannot depend on a
    report, manifest or receipt.
    """
    return digest(build_inputs(lock, inventory, recipe), BUILD_INPUTS_DOMAIN)
