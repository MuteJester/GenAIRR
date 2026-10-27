"""Cross-document consistency checks between a lock, recipe and source inventory.

Each document validates itself on construction; these checks make sure documents that
reference one another agree, so a rehashed but inconsistent document is rejected.
Paths locate the offending field inside the first document argument.
"""
from __future__ import annotations

from GenAIRR._reference_sources._fields import bad, ptr
from GenAIRR._reference_sources.common import PROFILE_SCOPES
from GenAIRR._reference_sources.inventory import COORDINATE_PREDICATES, RUNTIME_PREDICATES, SourceInventory
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.recipe import Recipe


def _require(value, cls, name):
    if not isinstance(value, cls):
        raise TypeError(f"{name} must be a {cls.__name__}, got {type(value).__name__}")


def profile_scope_problem(lock: SourceLock):
    """``None`` if the lock lies inside its interpretation profile's validated scope, else why not."""
    profile = lock.adapter.coordinate_profile
    scope = PROFILE_SCOPES.get((profile.id, profile.version), ())
    key = (lock.identity.taxon.taxonomy_id, lock.set_id, lock.release)
    formats = next((entry[3] for entry in scope if entry[:3] == key), None)
    if formats is None:
        validated = [f"{s} release {r} (taxon {t})" for t, s, r, _ in scope]
        return (f"profile {profile.id} version {profile.version} has been validated only for {validated}; "
                f"this lock is {lock.set_id} release {lock.release} (taxon {key[0]})")
    unvalidated = [a.format for a in lock.assets if a.role in ("record_export", "gapped_alignment")
                   and a.format not in formats]
    if unvalidated:
        return f"profile {profile.id} version {profile.version} has not been validated for formats {unvalidated}"
    return None


def check_profile_scope(lock: SourceLock) -> None:
    """``unsupported_profile`` unless the lock lies inside its profile's validated scope."""
    _require(lock, SourceLock, "lock")
    problem = profile_scope_problem(lock)
    if problem:
        raise bad("", problem, "unsupported_profile")


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


def check_inventory(inventory: SourceInventory, lock: SourceLock) -> None:
    """The inventory was generated from exactly this lock and refers only to its assets.

    ``digest_mismatch``: ``source_lock_hash`` is not this lock's commitment.
    ``identity_mismatch``: a different adapter, or a locator into an asset the lock lacks.
    """
    _require(inventory, SourceInventory, "inventory")
    _require(lock, SourceLock, "lock")
    if inventory.source_lock_hash != lock.commitment:
        raise bad("/source_lock_hash", f"inventory was generated from lock {inventory.source_lock_hash}, but this "
                                       f"lock is {lock.commitment}", "digest_mismatch")
    if inventory.adapter != lock.adapter:
        raise bad("/adapter", f"inventory adapter {inventory.adapter.to_dict()} differs from the lock adapter "
                              f"{lock.adapter.to_dict()}", "identity_mismatch")
    locked = set(lock.asset_ids)
    for path, asset_id in inventory.asset_references():
        if asset_id not in locked:
            raise bad(path, f"asset {asset_id!r} is not in the lock (locked: {sorted(locked)})", "identity_mismatch")
    problem = profile_scope_problem(lock)
    if problem:
        for i, claim in enumerate(inventory.claims):
            if claim.status == "accepted" and claim.predicate in COORDINATE_PREDICATES:
                raise bad(ptr(ptr("/claims", i), "status"), f"accepted {claim.predicate} claim about "
                                                            f"{claim.subject.id}: {problem}", "unsupported_profile")


def check_selection(recipe: Recipe, inventory: SourceInventory) -> None:
    """Every selection entry resolves to one exact record and one of its candidates.

    ``unresolved_selection``: no such record, a malformed/unsupported record, or no
    candidate of that record with the selected digest (no first-name matching).
    ``identity_mismatch``: the record's segment, label, taxon or locus differ from the
    entry or the recipe's reference (names are never changed).
    """
    _require(recipe, Recipe, "recipe")
    _require(inventory, SourceInventory, "inventory")
    if recipe.source_lock_hash != inventory.source_lock_hash:
        raise bad("/source_lock_hash", f"recipe builds from lock {recipe.source_lock_hash}, but the inventory is of "
                                       f"lock {inventory.source_lock_hash}", "digest_mismatch")
    reference = recipe.reference
    for i, entry in enumerate(recipe.selection):
        path = ptr("/selection", i)
        try:
            record = inventory.record(entry.record)
        except KeyError:
            raise bad(ptr(path, "record"), f"{entry.label}: no source record at {entry.record.id}",
                      "unresolved_selection") from None
        if record.status != "parsed":
            raise bad(ptr(path, "record"), f"{entry.label}: record {record.id} is {record.status} "
                                           f"({record.problem}) and cannot be selected", "unresolved_selection")
        if record.segment != entry.segment:
            raise bad(ptr(path, "segment"), f"{entry.label}: record {record.id} declares segment {record.segment}, "
                                            f"not {entry.segment}", "identity_mismatch")
        if record.label != entry.label:
            raise bad(ptr(path, "label"), f"record {record.id} is labelled {record.label!r}, not {entry.label!r}; "
                                          "labels are never changed", "identity_mismatch")
        if record.taxonomy_id != reference.taxon.taxonomy_id or record.locus != reference.locus:
            raise bad(ptr(path, "record"), f"{entry.label}: record {record.id} declares taxon {record.taxonomy_id} "
                                           f"locus {record.locus}, not {reference.taxon.taxonomy_id} "
                                           f"{reference.locus}", "identity_mismatch")
        if not any(c.sha256 == entry.sequence_sha256 for c in inventory.candidates_of(entry.record)):
            raise bad(ptr(path, "sequence_sha256"), f"{entry.label}: record {record.id} supplies no candidate with "
                                                    f"digest {entry.sequence_sha256}", "unresolved_selection")


def selected_subjects(recipe: Recipe, inventory: SourceInventory) -> set:
    """Ids of the selected records and of the one selected candidate of each."""
    found = set()
    for entry in recipe.selection:
        found.add(entry.record.id)
        found.update(c.id for c in inventory.candidates_of(entry.record) if c.sha256 == entry.sequence_sha256)
    return found


def _subject_digests(inventory: SourceInventory, claim) -> set:
    try:
        return {inventory.candidate(claim.subject).sha256}
    except KeyError:
        return {c.sha256 for c in inventory.candidates_of(claim.subject)}


def check_decisions(recipe: Recipe, inventory: SourceInventory) -> None:
    """Every decision still names current findings, claims, evidence and subject digest.

    ``stale_decision``: an unknown finding/claim id, an evidence locator that is not an
    observation, a subject digest that changed, a claim about an unselected record or
    candidate, or a ``select_alternative`` whose digest is not a selected candidate among
    its findings' subjects.  ``projection``: an
    ``accept_claim`` of a claim the runtime reference cannot represent or that has no value.
    """
    _require(recipe, Recipe, "recipe")
    _require(inventory, SourceInventory, "inventory")
    selected = selected_subjects(recipe, inventory)
    for i, decision in enumerate(recipe.decisions):
        path = ptr("/decisions", i)
        findings, claims = [], []
        for j, finding_id in enumerate(decision.findings):
            try:
                findings.append(inventory.finding(finding_id))
            except KeyError:
                raise bad(ptr(ptr(path, "findings"), j), f"finding {finding_id} no longer exists", "stale_decision") \
                    from None
        for j, claim_id in enumerate(decision.claims):
            try:
                claims.append(inventory.claim(claim_id))
            except KeyError:
                raise bad(ptr(ptr(path, "claims"), j), f"claim {claim_id} no longer exists", "stale_decision") from None
        for j, locator in enumerate(decision.evidence):
            try:
                inventory.observation(locator)
            except KeyError:
                raise bad(ptr(ptr(path, "evidence"), j), f"evidence {locator.id} is not an observed source field",
                          "stale_decision") from None
        for claim in claims:
            if decision.sequence_sha256 not in _subject_digests(inventory, claim):
                raise bad(ptr(path, "sequence_sha256"), f"claim {claim.id} is about {claim.subject.id}, whose "
                                                        f"sequence is not {decision.sequence_sha256}", "stale_decision")
        if decision.action == "accept_claim":
            for j, claim in enumerate(claims):
                if claim.predicate not in RUNTIME_PREDICATES:
                    raise bad(ptr(ptr(path, "claims"), j), f"{claim.predicate} claims are not representable in a "
                                                           f"runtime reference ({list(RUNTIME_PREDICATES)})",
                              "projection")
                if claim.status == "unresolved":
                    raise bad(ptr(ptr(path, "claims"), j), f"claim {claim.id} is unresolved and has no value to "
                                                           "accept", "projection")
        for j, claim in enumerate(claims):
            if claim.subject.id not in selected:
                raise bad(ptr(ptr(path, "claims"), j), f"claim {claim.id} is about {claim.subject.id}, which is not "
                                                       "a selected record or candidate", "stale_decision")
        if decision.action == "select_alternative":
            subjects = {s.id for f in findings for s in f.subjects}
            chosen = [e for e in recipe.selection if e.sequence_sha256 == decision.sequence_sha256 and (
                e.record.id in subjects or any(c.id in subjects and c.sha256 == e.sequence_sha256
                                               for c in inventory.candidates_of(e.record)))]
            if not chosen:
                raise bad(ptr(path, "sequence_sha256"), f"no selected record with digest {decision.sequence_sha256} "
                                                        "is a subject of the decision's findings", "stale_decision")


def check_build_inputs(lock: SourceLock, inventory: SourceInventory, recipe: Recipe) -> None:
    """All cross-document checks a build performs before any projection or native work."""
    check_profile_scope(lock)
    check_recipe(recipe, lock)
    check_inventory(inventory, lock)
    check_selection(recipe, inventory)
    check_decisions(recipe, inventory)
