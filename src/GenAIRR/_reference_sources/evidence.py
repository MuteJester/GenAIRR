"""Evidence classification: claims, mappings and findings from an OGRDB source inventory.

:func:`interpret` reads the uninterpreted inventory of a verified lock (the inventory
:func:`~GenAIRR._reference_sources.ogrdb.build_inventory` gives for exactly these bytes)
under the interpretation profile the lock names, and returns a new inventory that adds
claims, mappings and findings.  Records, candidates, mentioned identities, observations,
support references and the earlier findings are kept unchanged.

Every claim says how strong it is:

* ``support``: ``source`` (a provider field states the value), ``assertion`` (GenAIRR
  derives it from source values) or ``none``;
* ``derivation``: ``stated``, ``anchor_derived`` (a phase from the anchor modulo 3) or
  ``variant_comparison`` (from comparing the two record exports);
* ``status``: ``accepted`` only for a value whose meaning the profile establishes on
  validated bytes and whose per-record checks hold; ``candidate`` for such a value outside
  the profile's validated scope or one that disagrees with a related value; otherwise
  ``unresolved`` with no value.

Coordinate values are never repaired, relocated, extended or chosen by vote: a value that
fails a check stays unresolved (or a candidate) with a finding that says why.  Each claim
cites only the fields of its own record and the gapped record its own exact mapping binds,
so identities with identical DNA never share evidence.
"""
from __future__ import annotations

from collections import Counter
from typing import Any, Mapping, Optional

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference_sources import consistency, coordinates as c, profiles
from GenAIRR._reference_sources._fields import MAX_REFS, bad, ptr, require_type
from GenAIRR._reference_sources.common import GAPPED_VARIANT, RECORD_EXPORT_FORMATS, Locator
from GenAIRR._reference_sources.digests import lock_hash
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.inventory import (
    IMGT_GAP, Block, Claim, Finding, Mapping as SourceMapping, SourceInventory, alignment_blocks, text_digest,
)
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.ogrdb import EXPORT_ANCHOR, build_inventory

_JUNCTION = ("source_annotated_junction", "canonical_productive")

# code -> (severity, category, affected capabilities, meaning).
FINDING_CODES = {
    "coordinate_profile_out_of_scope": (
        "warning", "unsupported_interpretation", _JUNCTION,
        "The lock is outside the interpretation profile's validated scope (another set or release, or bytes "
        "other than the validated ones); its coordinate claims are only candidates or unresolved."),
    "coding_interval_unusable": (
        "warning", "invalid_source_syntax", (),
        "gene_start, gene_end or the record sequence give no interval of the record sequence; the coding "
        "interval is unresolved."),
    "coding_interval_not_reproduced": (
        "warning", "contradictory_evidence", ("original_coding_ends",),
        "sequence[gene_start - 1 : gene_end] is not the supplied coding sequence; how they relate is recorded "
        "and no interval is applied beyond what reproduces."),
    "v_landmark_field_unusable": (
        "warning", "invalid_source_syntax", _JUNCTION,
        "A V record has no single IMGT delineation, or its cdr3_start is absent, null, 0 (served for a missing "
        "value), a Boolean, not an integer or out of the coding sequence's range; the anchor is unresolved."),
    "v_landmark_alignment_unavailable": (
        "warning", "unavailable_evidence", _JUNCTION,
        "No single gapped record of the same export variant carries exactly the V record's label (or the "
        "gapped export is not locked, or the label is not unique in its export), so cdr3_start cannot be "
        "bound to the record and the anchor is unresolved."),
    "v_landmark_binding_failed": (
        "warning", "contradictory_evidence", _JUNCTION,
        "The V record's gapped record does not ungap to its coding sequence, has fewer than 312 columns, holds "
        "a gap in the landmark codon at columns 310-312, or gives a cdr3_start other than the record's; the "
        "anchor is unresolved and nothing is relocated."),
    "j_landmark_field_unusable": (
        "warning", "invalid_source_syntax", _JUNCTION,
        "A J record's j_cdr3_end or j_codon_frame is absent, null, a Boolean, not an integer or out of range; "
        "that value is unresolved."),
    "phase_landmark_projection_unsupported": (
        "warning", "contradictory_evidence", _JUNCTION,
        "The stated J phase is not the stated landmark position modulo 3; neither value is shifted and both "
        "stay candidates."),
    "noncanonical_landmark_codon": (
        "warning", "informational", ("canonical_productive",),
        "The codon at an anchor does not encode the conserved residue (Cys for V; Trp for IGH J; Phe for IGK "
        "and IGL J); the anchor stays where the source puts it and no other codon is used."),
}


def _locator(record, *tokens) -> Locator:
    pointer = record.locator.pointer
    for token in tokens:
        pointer = ptr(pointer, token)
    return Locator(record.locator.asset, pointer)


def _canonical(locators) -> tuple:
    return tuple(sorted(set(locators), key=Locator.sort_key))


class _Findings:
    """Findings keyed by code and observations: equal observations (for example the same
    allele in both record exports) form one finding naming every subject."""

    def __init__(self, lock_hash_: str):
        self.lock_hash = lock_hash_
        self.groups: dict = {}

    def add(self, code: str, *, subjects=(), locators=(), observations: dict) -> None:
        entry = self.groups.setdefault((code, canonical_json(observations)), (set(), set(), observations))
        entry[0].update(subjects)
        entry[1].update(locators)

    def build(self) -> list:
        built = []
        for (code, _), (subjects, locators, observations) in self.groups.items():
            severity, category, capabilities, _ = FINDING_CODES[code]
            subjects, locators = _canonical(subjects), _canonical(locators)
            parts = max(1, -(-len(subjects) // MAX_REFS), -(-len(locators) // MAX_REFS))
            for part in range(parts):
                details = dict(observations)
                if parts > 1:
                    details.update(part=part + 1, parts=parts)
                window = slice(part * MAX_REFS, (part + 1) * MAX_REFS)
                built.append(Finding(source_lock_hash=self.lock_hash, code=code, severity=severity,
                                     category=category, subjects=subjects[window], locators=locators[window],
                                     observations=details, capabilities=capabilities))
        return built


class _Interpretation:
    """The claims, mappings and findings of one inventory under one profile."""

    def __init__(self, lock: SourceLock, inventory: SourceInventory, profile):
        self.lock = lock
        self.inventory = inventory
        self.profile = profile
        self.lock_hash = inventory.source_lock_hash
        self.scope_problem = consistency.profile_scope_problem(lock)
        self.locus = lock.identity.locus
        self.claims: list = []
        self.mappings: list = []
        self.findings = _Findings(self.lock_hash)
        self.exports = [a.asset_id for a in lock.assets if a.asset_id in RECORD_EXPORT_FORMATS]
        self.parsed = {asset: [r for r in inventory.records if r.locator.asset == asset and r.status == "parsed"]
                       for asset in self.exports}
        self.gapped: dict = {}
        for observation in inventory.observations:
            if observation.locator.role == "gapped_alignment":
                value = observation.value.get()
                self.gapped.setdefault(observation.locator.asset, {}).setdefault(value["header"], []).append(
                    (observation.locator, value["sequence"]))

    # -- reading ---------------------------------------------------------

    def value(self, locator: Locator) -> Any:
        observation = self.inventory.observation(locator)
        return c.ABSENT if observation.presence == "absent" else observation.value.get()

    def field(self, record, *tokens) -> tuple:
        locator = _locator(record, *tokens)
        return locator, self.value(locator)

    # -- building --------------------------------------------------------

    def status(self, resolved: bool, *, support: str = "source", candidate: bool = False) -> tuple:
        """``(status, support)`` of a value: unresolved without one, a candidate outside the
        validated scope or when a related value disagrees, accepted otherwise."""
        if not resolved:
            return "unresolved", "none"
        if candidate or self.scope_problem is not None:
            return "candidate", "none"
        return "accepted", support

    def claim(self, subject: Locator, predicate: str, value: Any, system: Optional[str], derivation: str,
              graded: tuple, evidence, mappings=()) -> Claim:
        status, support = graded
        claim = Claim(source_lock_hash=self.lock_hash, subject=subject, predicate=predicate,
                      value=None if status == "unresolved" else value, coordinate_system=system, support=support,
                      derivation=derivation, status=status, evidence=_canonical(evidence),
                      mappings=tuple(sorted(set(mappings))))
        self.claims.append(claim)
        return claim

    def mapping(self, **fields) -> SourceMapping:
        mapping = SourceMapping(source_lock_hash=self.lock_hash, **fields)
        self.mappings.append(mapping)
        return mapping

    def finding(self, code: str, record, candidate=None, *, locators=(), **observations) -> None:
        subjects = [record.locator] + ([candidate.locator] if candidate is not None else [])
        self.findings.add(code, subjects=subjects, locators=locators,
                          observations={"label": record.label, "segment": record.segment, **observations})

    # -- the interpretation ----------------------------------------------

    def run(self) -> None:
        if self.scope_problem is not None:
            self.findings.add("coordinate_profile_out_of_scope",
                              locators=[Locator(a, EXPORT_ANCHOR) for a in self.exports],
                              observations={"profile": self.profile.id, "version": self.profile.version,
                                            "problem": self.scope_problem})
        for asset in self.exports:
            labels = Counter(r.label for r in self.parsed[asset])
            for record in self.parsed[asset]:
                candidates = self.inventory.candidates_of(record.locator)
                if not candidates:
                    continue
                (candidate,) = candidates
                self.coding_interval(record, candidate)
                if record.segment == "V":
                    self.v_landmark(record, candidate, labels[record.label])
                elif record.segment == "J":
                    self.j_landmark(record, candidate)

    def coding_interval(self, record, candidate) -> c.CodingInterval:
        (seq_at, sequence), (start_at, start), (end_at, end) = (self.field(record, name) for name in
                                                                ("sequence", "gene_start", "gene_end"))
        found = c.coding_interval(sequence, start, end, candidate.sequence)
        mappings = []
        if found.start is not None:
            blocks = (Block((found.start, found.end), (0, found.end - found.start)),) if found.reproduces else ()
            mappings.append(self.mapping(
                kind="coding_interval", input=seq_at, output=candidate.locator, input_sha256=text_digest(sequence),
                output_sha256=candidate.sha256, input_space="record_sequence", output_space="coding_sequence",
                blocks=blocks, transform="extract_interval", validation="exact" if found.reproduces else "mismatch").id)
        self.claim(candidate.locator, "coding_interval", [found.start, found.end], "record_sequence", "stated",
                   self.status(found.reproduces is True), (start_at, end_at, seq_at), mappings)
        if found.problem is not None:
            self.finding("coding_interval_unusable", record, candidate, locators=(start_at, end_at, seq_at),
                         reason=found.problem)
        elif not found.reproduces:
            self.finding("coding_interval_not_reproduced", record, candidate, locators=(start_at, end_at, seq_at),
                         gene_start=start, gene_end=end, declared_length=found.end - found.start,
                         coding_length=len(candidate.sequence), relation=found.relation, extra=found.extra)
        return found

    def gapped_record(self, record, carriers: int) -> tuple:
        """``(locator, gapped sequence, None)`` of the one record of the same variant's gapped
        export labelled exactly as ``record``, or ``(None, None, why not)``."""
        gapped_asset = GAPPED_VARIANT[record.locator.asset]
        if gapped_asset not in self.lock.asset_ids:
            return None, None, f"the {gapped_asset} export is not locked"
        if carriers != 1:
            return None, None, (f"{carriers} records of {record.locator.asset} carry the label {record.label}; "
                                "a gapped record cannot be bound to one of them")
        found = self.gapped.get(gapped_asset, {}).get(record.label, [])
        if not found:
            return None, None, f"no {gapped_asset} record is labelled exactly {record.label}"
        if len(found) > 1:
            return None, None, f"{len(found)} {gapped_asset} records are labelled {record.label}"
        (locator, sequence), = found
        return locator, sequence, None

    def gapped_mapping(self, candidate, locator: Locator, gapped: str) -> Optional[SourceMapping]:
        exact = gapped.replace(IMGT_GAP, "") == candidate.sequence
        blocks = alignment_blocks(gapped) if exact else ()
        if len(blocks) > MAX_REFS:
            return None
        return self.mapping(kind="gapped_alignment", input=candidate.locator, output=locator,
                            input_sha256=candidate.sha256, output_sha256=text_digest(gapped),
                            input_space="coding_sequence", output_space="imgt_gapped", blocks=blocks,
                            transform="remove_imgt_gaps", validation="exact" if exact else "mismatch")

    def v_landmark(self, record, candidate, carriers: int) -> None:
        delineations_at, delineations = self.field(record, "v_gene_delineations")
        scheme = self.profile.v_delineation_scheme
        imgt = [j for j, d in enumerate(delineations if isinstance(delineations, list) else ())
                if isinstance(d, dict) and d.get("delineation_scheme") == scheme]
        evidence = [delineations_at]
        if len(imgt) == 1:
            evidence = [_locator(record, "v_gene_delineations", imgt[0], name)
                        for name in ("cdr3_start", "delineation_scheme")]
        gapped_at, gapped, unavailable = self.gapped_record(record, carriers)
        mapping = None
        if gapped_at is not None:
            mapping = self.gapped_mapping(candidate, gapped_at, gapped)
            if mapping is None:
                gapped_at, gapped = None, None
                unavailable = f"the gapped record of {record.label} has more than {MAX_REFS} aligned blocks"
            else:
                evidence.append(gapped_at)
        if len(imgt) != 1:
            found = c.VLandmark(None, f"{len(imgt) or 'no'} {scheme} delineation{'s' * (len(imgt) > 1)} in "
                                      "v_gene_delineations", "field", None, None)
        else:
            found = c.v_landmark(candidate.sequence, self.value(evidence[0]), gapped)
        resolved = found.anchor is not None
        mappings = [mapping.id] if mapping is not None else []
        self.claim(candidate.locator, "anchor", found.anchor, "coding_sequence", "stated", self.status(resolved),
                   evidence, mappings)
        phase = found.anchor % 3 if resolved else None
        graded = self.status(resolved, support="assertion")
        self.claim(candidate.locator, "coding_phase", phase, "coding_sequence", "anchor_derived", graded, evidence,
                   mappings)
        self.claim(candidate.locator, "phase_independence", "anchor_derived", None, "anchor_derived", graded,
                   evidence, mappings)
        if not resolved:
            code = {"field": "v_landmark_field_unusable", "alignment": "v_landmark_alignment_unavailable",
                    "binding": "v_landmark_binding_failed"}[found.kind]
            reason = unavailable if found.kind == "alignment" else found.problem
            self.finding(code, record, candidate, locators=evidence, reason=reason)
        else:
            self.landmark_codon(record, candidate, found.anchor, evidence)

    def j_landmark(self, record, candidate) -> None:
        (end_at, end), (frame_at, frame) = self.field(record, "j_cdr3_end"), self.field(record, "j_codon_frame")
        found = c.j_landmark(candidate.sequence, end, frame)
        consistent = found.consistent
        self.claim(candidate.locator, "anchor", found.anchor, "coding_sequence", "stated",
                   self.status(found.anchor is not None, candidate=consistent is False), (end_at,))
        self.claim(candidate.locator, "coding_phase", found.phase, "coding_sequence", "stated",
                   self.status(found.phase is not None, candidate=consistent is not True), (frame_at,))
        self.claim(candidate.locator, "phase_independence", None, None, "stated", self.status(False),
                   (end_at, frame_at))
        for name, at, problem in (("j_cdr3_end", end_at, found.anchor_problem),
                                  ("j_codon_frame", frame_at, found.phase_problem)):
            if problem is not None:
                self.finding("j_landmark_field_unusable", record, candidate, locators=(at,), field=name,
                             reason=problem)
        if consistent is False:
            self.finding("phase_landmark_projection_unsupported", record, candidate, locators=(end_at, frame_at),
                         anchor=found.anchor, phase=found.phase, anchor_modulo_3=found.anchor % 3)
        if found.anchor is not None:
            self.landmark_codon(record, candidate, found.anchor, (end_at,))

    def landmark_codon(self, record, candidate, anchor: int, evidence) -> None:
        expected = self.profile.landmark_residues[record.segment].get(self.locus, ())
        codon = candidate.sequence[anchor:anchor + 3]
        residue = c.residue(codon)
        if expected and residue not in expected:
            self.finding("noncanonical_landmark_codon", record, candidate, locators=evidence, anchor=anchor,
                         codon=codon, residue=residue, expected_residues=list(expected),
                         coordinate_system="coding_sequence")


def interpret(lock: SourceLock, assets: Mapping[str, bytes], inventory: SourceInventory) -> SourceInventory:
    """The inventory of ``lock`` with its claims, mappings and findings.

    ``assets`` are the verified locked bytes (asset id -> bytes) and ``inventory`` the
    uninterpreted inventory built from exactly them.  The interpretation profile is the one
    the lock's adapter names (``unsupported_profile`` if it is not implemented); it is never
    chosen from the data.  Raises ``digest_mismatch`` at ``/source_lock_hash`` for an
    inventory of another lock, the byte checks of ``build_inventory`` for bytes that are not
    the locked ones, and ``identity_mismatch`` at ``/inventory`` for an inventory that is not
    the one those bytes give (an edited or already interpreted inventory).
    """
    require_type(lock, SourceLock, "lock")
    require_type(inventory, SourceInventory, "inventory")
    if inventory.source_lock_hash != lock_hash(lock):
        raise bad("/source_lock_hash", f"the inventory was built under lock {inventory.source_lock_hash}, not this "
                                       f"lock {lock_hash(lock)}", "digest_mismatch")
    profile = profiles.profile_for(lock.adapter.coordinate_profile)
    if build_inventory(lock, assets) != inventory:
        raise bad("/inventory", "the inventory is not the uninterpreted inventory of these locked bytes",
                  "identity_mismatch")
    work = _Interpretation(lock, inventory, profile)
    work.run()
    fields = {name: getattr(inventory, name) for name in ("records", "candidates", "mentioned_identities",
                                                         "observations", "support")}
    try:
        result = SourceInventory.assemble(source_lock_hash=inventory.source_lock_hash, adapter=inventory.adapter,
                                          claims=work.claims, mappings=work.mappings,
                                          findings=inventory.findings + tuple(work.findings.build()), **fields)
    except ReferenceSourceError as exc:   # an inventory bound such as the number of claims
        raise exc.under("/inventory") from None
    consistency.check_inventory(result, lock)
    return result
