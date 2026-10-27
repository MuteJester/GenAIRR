"""Cross-document consistency checks between a lock, recipe and source inventory.

Each document validates itself on construction; these checks make sure documents that
reference one another agree, so a rehashed but inconsistent document is rejected.
Paths locate the offending field inside the first document argument.
"""
from __future__ import annotations

from GenAIRR._reference_sources._fields import bad, ptr
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.recipe import Recipe


def _require(value, cls, name):
    if not isinstance(value, cls):
        raise TypeError(f"{name} must be a {cls.__name__}, got {type(value).__name__}")


def check_recipe(recipe: Recipe, lock: SourceLock) -> None:
    """The recipe builds from exactly this lock and cannot relabel its identity.

    ``digest_mismatch``: ``source_lock_hash`` is not this lock's commitment.
    ``identity_mismatch``: adapter, taxon, locus or chain type differ from the lock, or a
    selected record / decision evidence names an asset the lock does not hold (selected
    records must come from a record export).
    """
    _require(recipe, Recipe, "recipe")
    _require(lock, SourceLock, "lock")
    if recipe.source_lock_hash != lock.commitment:
        raise bad("/source_lock_hash", f"recipe builds from lock {recipe.source_lock_hash}, but this lock is "
                                       f"{lock.commitment} ({lock.set_id} release {lock.release})", "digest_mismatch")
    if recipe.adapter != lock.adapter:
        raise bad("/adapter", f"recipe adapter {recipe.adapter.to_dict()} differs from the lock adapter "
                              f"{lock.adapter.to_dict()}", "identity_mismatch")
    reference, identity = recipe.reference, lock.identity
    if reference.taxon != identity.taxon:
        raise bad("/reference/taxon", f"recipe taxon {reference.taxon.to_dict()} differs from the source taxon "
                                      f"{identity.taxon.to_dict()}", "identity_mismatch")
    if reference.locus != identity.locus:
        raise bad("/reference/locus", f"recipe locus {reference.locus!r} differs from the source locus "
                                      f"{identity.locus!r}", "identity_mismatch")
    if reference.chain_type != identity.chain_type:
        raise bad("/reference/chain_type", f"recipe chain type {reference.chain_type!r} differs from the source "
                                           f"chain type {identity.chain_type!r}", "identity_mismatch")
    roles = {a.asset_id: a.role for a in lock.assets}
    for i, entry in enumerate(recipe.selection):
        path = ptr(ptr(ptr("/selection", i), "record"), "asset")
        if entry.record.asset not in roles:
            raise bad(path, f"{entry.label}: asset {entry.record.asset!r} is not in the lock "
                            f"(locked: {list(roles)})", "identity_mismatch")
        if roles[entry.record.asset] != "record_export":
            raise bad(path, f"{entry.label}: asset {entry.record.asset!r} is a {roles[entry.record.asset]} "
                            "asset, not a record export", "identity_mismatch")
    for i, decision in enumerate(recipe.decisions):
        for j, locator in enumerate(decision.evidence):
            if locator.asset not in roles:
                raise bad(ptr(ptr(ptr(ptr("/decisions", i), "evidence"), j), "asset"),
                          f"evidence asset {locator.asset!r} is not in the lock (locked: {list(roles)})",
                          "identity_mismatch")
