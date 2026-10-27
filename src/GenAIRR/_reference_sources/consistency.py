"""Cross-document consistency checks between a lock, recipe and source inventory.

Each document validates itself on construction; these checks make sure documents that
reference one another agree, so a rehashed but inconsistent document is rejected.
Paths locate the offending field inside the first document argument.
"""
from __future__ import annotations

from GenAIRR._reference_sources._fields import bad, ptr, require_type
from GenAIRR._reference_sources.common import PROFILE_SCOPES
from GenAIRR._reference_sources.inventory import COORDINATE_PREDICATES, RUNTIME_PREDICATES, SourceInventory
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.recipe import Recipe


def profile_scope_problem(lock: SourceLock):
    """``None`` if the lock lies inside its interpretation profile's validated scope, else why not."""
    profile = lock.adapter.coordinate_profile
    scope = PROFILE_SCOPES.get((profile.id, profile.version), ())
    key = (lock.identity.taxon.taxonomy_id, lock.set_id, lock.release)
    pinned = next((entry[3] for entry in scope if entry[:3] == key), None)
    if pinned is None:
        validated = [f"{s} release {r} (taxon {t})" for t, s, r, _ in scope]
        return (f"profile {profile.id} version {profile.version} has been validated only for {validated}; "
                f"this lock is {lock.set_id} release {lock.release} (taxon {key[0]})")
    exact = {fmt: (sha256, size) for fmt, sha256, size in pinned}
    for asset in lock.assets:
        if asset.role not in ("record_export", "gapped_alignment"):
            continue
        if asset.format not in exact:
            return (f"profile {profile.id} version {profile.version} has not been validated for format "
                    f"{asset.format!r} of {lock.set_id} release {lock.release}")
        if (asset.sha256, asset.bytes) != exact[asset.format]:
            sha256, size = exact[asset.format]
            return (f"{asset.asset_id} of {lock.set_id} release {lock.release} is {asset.sha256} ({asset.bytes} "
                    f"bytes), not the validated {sha256} ({size} bytes); profile {profile.id} version "
                    f"{profile.version} has not been validated for these bytes")
    return None


def check_profile_scope(lock: SourceLock) -> None:
    """``unsupported_profile`` unless the lock lies inside its profile's validated scope."""
    require_type(lock, SourceLock, "lock")
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
    require_type(recipe, Recipe, "recipe")
    require_type(lock, SourceLock, "lock")
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
    ``unsupported_profile``: a lock outside its profile's validated scope with an accepted or
    conflicting coordinate claim (only candidates and unresolved claims are allowed there).
    """
    require_type(inventory, SourceInventory, "inventory")
    require_type(lock, SourceLock, "lock")
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
            if claim.status in ("accepted", "conflicting") and claim.predicate in COORDINATE_PREDICATES:
                raise bad(ptr(ptr("/claims", i), "status"), f"{claim.status} {claim.predicate} claim about "
                                                            f"{claim.subject.id}: outside the validated scope a "
                                                            f"coordinate claim is only a candidate or unresolved; "
                                                            f"{problem}", "unsupported_profile")


def check_selection(recipe: Recipe, inventory: SourceInventory) -> None:
    """Every selection entry resolves to one exact record and one of its candidates.

    ``unresolved_selection``: no such record, a malformed/unsupported record, or no
    candidate of that record with the selected digest (no first-name matching).
    ``identity_mismatch``: the record's segment, label, taxon or locus differ from the
    entry or the recipe's reference (names are never changed).
    """
    require_type(recipe, Recipe, "recipe")
    require_type(inventory, SourceInventory, "inventory")
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


def _selected_candidates(recipe: Recipe, inventory: SourceInventory) -> dict:
    """``{record id: the selected candidate of that record}`` for every selection entry."""
    chosen = {}
    for entry in recipe.selection:
        for candidate in inventory.candidates_of(entry.record):
            if candidate.sha256 == entry.sequence_sha256:
                chosen[entry.record.id] = candidate
    return chosen


def selected_subjects(recipe: Recipe, inventory: SourceInventory) -> set:
    """Ids of the selected records and of the one selected candidate of each."""
    return {e.record.id for e in recipe.selection} | {c.id for c in _selected_candidates(recipe, inventory).values()}


def check_decisions(recipe: Recipe, inventory: SourceInventory) -> None:
    """Every decision still names current findings, claims, evidence and subject digest.

    Checked per decision, in order:

    * ``stale_decision``: an unknown finding or claim id, or evidence that is not an
      observed source field;
    * ``projection``: ``withhold_claim``/``accept_claim`` of a claim the runtime reference
      cannot represent, or ``accept_claim`` of an unresolved claim (no value);
    * ``stale_decision``: a claim about an unselected record or candidate, or a
      ``sequence_sha256`` that is not the subject's selected sequence (a candidate's own
      digest; for a record, its selection entry's digest);
    * ``identity_mismatch``: ``accept_claim`` evidence not bound to every decided claim's
      subject (the claim evidence-root rule);
    * ``stale_decision``: a ``select_alternative`` whose digest is not a selected candidate
      among its findings' subjects.
    """
    require_type(recipe, Recipe, "recipe")
    require_type(inventory, SourceInventory, "inventory")
    entries = {e.record.id: e for e in recipe.selection}
    chosen = _selected_candidates(recipe, inventory)
    selected = set(entries) | {c.id for c in chosen.values()}
    by_digest: dict = {}
    for entry in recipe.selection:
        by_digest.setdefault(entry.sequence_sha256, []).append(entry)
    for i, decision in enumerate(recipe.decisions):
        path = ptr("/decisions", i)
        findings, claims = [], []
        for j, finding_id in enumerate(decision.findings):
            try:
                findings.append(inventory.finding(finding_id))
            except KeyError:
                raise bad(ptr(ptr(path, "findings"), j), f"finding {finding_id} no longer exists",
                          "stale_decision") from None
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
        for j, claim in enumerate(claims):
            where = ptr(ptr(path, "claims"), j)
            if claim.predicate not in RUNTIME_PREDICATES:
                raise bad(where, f"{decision.action}: {claim.predicate} claims are not representable in a runtime "
                                 f"reference ({list(RUNTIME_PREDICATES)})", "projection")
            if decision.action == "accept_claim" and claim.status == "unresolved":
                raise bad(where, f"claim {claim.id} is unresolved and has no value to accept", "projection")
        for j, claim in enumerate(claims):
            if claim.subject.id not in selected:
                raise bad(ptr(ptr(path, "claims"), j), f"claim {claim.id} is about {claim.subject.id}, which is not "
                                                       "a selected record or candidate", "stale_decision")
        for claim in claims:
            record = inventory.subject_record(claim)
            expected = entries[record.id].sequence_sha256
            if decision.sequence_sha256 != expected:
                raise bad(ptr(path, "sequence_sha256"), f"claim {claim.id} is about {claim.subject.id}, whose "
                                                        f"selected sequence is {expected}, not "
                                                        f"{decision.sequence_sha256}", "stale_decision")
        if decision.action == "accept_claim":
            for j, locator in enumerate(decision.evidence):
                for claim in claims:
                    problem = inventory.evidence_binding_problem(claim, locator)
                    if problem:
                        raise bad(ptr(ptr(path, "evidence"), j), f"claim {claim.id}: {problem}",
                                  "identity_mismatch")
        if decision.action == "select_alternative":
            subjects = {s.id for f in findings for s in f.subjects}
            if not any(e.record.id in subjects or (e.record.id in chosen and chosen[e.record.id].id in subjects)
                       for e in by_digest.get(decision.sequence_sha256, ())):
                raise bad(ptr(path, "sequence_sha256"), f"no selected record with digest {decision.sequence_sha256} "
                                                        "is a subject of the decision's findings", "stale_decision")


def check_build_inputs(lock: SourceLock, inventory: SourceInventory, recipe: Recipe) -> None:
    """All cross-document checks a build performs before any projection or native work."""
    check_profile_scope(lock)
    check_recipe(recipe, lock)
    check_inventory(inventory, lock)
    check_selection(recipe, inventory)
    check_decisions(recipe, inventory)
