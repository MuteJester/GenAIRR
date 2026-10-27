"""Evidence classification: claims, mappings and findings from an OGRDB source inventory.

:func:`interpret` reads the uninterpreted inventory of a verified lock (the inventory
:func:`~GenAIRR._reference_sources.ogrdb.build_inventory` gives for exactly these bytes)
under the interpretation profile the lock names, and returns a new inventory that adds
claims, mappings and findings.  Records, candidates, mentioned identities, observations,
support references and the earlier findings are kept unchanged.

Every claim says how strong it is:

* ``support``: ``source`` (a provider field states the value), ``assertion`` (GenAIRR
  derives it from source values: a phase from the anchor, an extension from comparing the
  two record exports) or ``none``;
* ``derivation``: ``stated``, ``anchor_derived`` (a phase from the anchor modulo 3) or
  ``variant_comparison`` (from comparing the two record exports);
* ``status``: ``accepted`` only for a value whose meaning the profile establishes on
  validated bytes and whose per-record checks hold; ``candidate`` for such a value outside
  the profile's validated scope or one that disagrees with a related value; otherwise
  ``unresolved`` with no value.

What is claimed:

* J: the anchor and the reported phase value are source-stated.  Whether that phase was
  derived independently of the landmark is unresolved (``phase_independence``), so a J
  phase is not independent phase evidence even though its derivation is ``stated``.
* V: the anchor is source-stated and accepted only when the record's own gapped record
  binds it; the phase is at most anchor-derived (``assertion``).
* Coding intervals are accepted only where they reproduce the coding sequence;
  provider extensions are recorded by comparing the two record exports, and a record
  context that is only the record's own end is marked as such, never read as support.
* Both coding boundaries and the provider functionality are always unresolved: no tag,
  Boolean, RSS adjacency or exact slice establishes a complete end or a functionality
  class.

Coordinate values are never repaired, relocated, extended or chosen by vote: a value that
fails a check stays unresolved (or a candidate) with a finding that says why.  Each claim
cites the fields of its own record and the gapped record its own mapping binds; an
extension claim also cites the same identity's record in the other record export, which
it compares.  Identities with identical DNA never share evidence.  All exports of one provider set
are one evidence lineage (:func:`evidence_lineage`): their agreement is consistency, not
independent confirmation, and embedded support entries are never retrieved or verified.
"""
from __future__ import annotations

from collections import Counter
import re
from typing import Any, Mapping, Optional

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference_sources import consistency, coordinates as c, profiles
from GenAIRR._reference_sources._fields import MAX_REFS, bad, ptr, require_type, show
from GenAIRR._reference_sources.common import GAPPED_VARIANT, RECORD_EXPORT_FORMATS, Locator
from GenAIRR._reference_sources.digests import lock_hash
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.inventory import (
    IMGT_GAP, Block, Claim, Finding, Mapping as SourceMapping, SourceInventory, alignment_blocks, text_digest,
)
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.ogrdb import EXPORT_ANCHOR, build_inventory

_JUNCTION = ("source_annotated_junction", "canonical_productive")
_ENDS = ("original_coding_ends",)
# The endpoint comparison's relation for one export's coding sequence extending the other's.
_EXTENDS = re.compile(r"(airr|airr_ex)_extends_([35])prime")
# The first bases of a heptamer-like V RSS start (CACAGTG); only reported, never used to place one.
_HEPTAMER_CORE = "CAC"
_EXTENSION_CODES = {"matches": "coding_extension_matches_record_context",
                    "without_record_context": "coding_extension_without_record_context",
                    "disagrees": "coding_extension_disagrees_with_record_context"}
# The record fields a boundary claim considers (none of them establishes a complete end).
_BOUNDARY_FIELDS = {
    ("V", "boundary_5"): ("gene_start", "curational_tags"),
    ("V", "boundary_3"): ("gene_end", "v_rs_start", "v_rs_end", "curational_tags"),
    ("D", "boundary_5"): ("gene_start", "d_rs_5_prime_start", "d_rs_5_prime_end", "curational_tags"),
    ("D", "boundary_3"): ("gene_end", "d_rs_3_prime_start", "d_rs_3_prime_end", "curational_tags"),
    ("J", "boundary_5"): ("gene_start", "j_rs_start", "j_rs_end", "curational_tags"),
    ("J", "boundary_3"): ("gene_end", "j_donor_splice", "curational_tags"),
}

# The profile conventions this module implements.  A profile declaring any other value is refused
# (``unsupported_profile``), so a profile field and the code that applies it cannot drift apart.
IMPLEMENTED_CONVENTIONS = {
    "rules": c.RULES,
    "landmark_columns": (c.LANDMARK_FIRST_COLUMN, c.LANDMARK_LAST_COLUMN),
    "missing_delineation_value": c.MISSING_DELINEATION_VALUE,
    "v_delineation_scheme": "IMGT",
    "j_phase_independence": "unresolved",
    "v_phase_derivation": "anchor_derived",
    "functionality_classes": {},
}

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
    "coding_extension_matches_record_context": (
        "info", "informational", _ENDS,
        "One export's coding sequence extends the other's, and the extension equals the record bases just "
        "beyond the declared coding interval.  The extension is a provider inference: neither end is adopted "
        "as the physical coding end, and where those bases are the record's own end "
        "(context_is_sequence_terminus) they are no genomic flank."),
    "coding_extension_without_record_context": (
        "info", "unavailable_evidence", _ENDS,
        "One export's coding sequence extends the other's, and the record holds no bases (or not all of them) "
        "at the extension's positions, or the declared interval does not prove those positions; the extension "
        "is not compared with the record."),
    "coding_extension_disagrees_with_record_context": (
        "warning", "contradictory_evidence", _ENDS,
        "One export's coding sequence extends the other's, and the extension differs from the record bases at "
        "the same positions; the same bases elsewhere in the record are no positional support."),
    "coding_end_before_declared_v_rss": (
        "warning", "contradictory_evidence", _ENDS,
        "Record bases lie between the declared coding end (gene_end) and the declared V RSS start; the shorter "
        "declared end may describe an affirmed span and a longer coding sequence is a candidate physical end. "
        "Neither is corrected or adopted; RSS coordinates are read as 1-based record positions only here."),
    "v_rss_start_inside_coding_interval": (
        "warning", "contradictory_evidence", (),
        "The declared V RSS start lies inside the declared coding interval, contradicting the record's own "
        "coding interval; the RSS coordinates are kept as observations and no coordinate is moved."),
    "v_rss_start_before_coding_interval": (
        "warning", "contradictory_evidence", (),
        "The declared V RSS start lies before the declared coding interval; the RSS coordinates are kept as "
        "observations and no coordinate is moved."),
    "v_rss_gap_after_coding_end": (
        "info", "informational", (),
        "Per record export: the gap between the declared coding end and the declared V RSS start (how many V "
        "records have none, the labels at every other gap), the labels without a declared or usable RSS start, "
        "those whose declared RSS start does not begin with a heptamer-like CAC, and those whose start has fewer "
        "than three record bases, so the core cannot be judged."),
    "declared_v_rss_start_without_cac": (
        "info", "informational", (),
        "The record bases at the declared V RSS start do not begin with CAC, the core of the heptamer CACAGTG; "
        "the declared position is reported as given and no other position is searched."),
    "v_rss_coordinate_beyond_record_sequence": (
        "warning", "invalid_source_syntax", (),
        "A declared V RSS start or end lies beyond the record sequence; a start there gives no position and "
        "no record bases are read past the sequence."),
    "affirmed_span_mapping_unresolved": (
        "info", "unavailable_evidence", _ENDS,
        "A curation note limits the allele's affirmed span (in gapped numbering) or marks likely further "
        "bases; without a reviewed interpretation, the exact note text, the sequence it was read against and an "
        "exact gapped alignment, no affirmed interval is derived.  Where the record holds the note text the "
        "profile read and the sequence it was read against, and the note states a column, the supplied span is "
        "compared with that column through the exact gapped alignment (a span running past it holds bases that "
        "were not affirmed).  For a record whose note text changed, only the digest of the text read is named."),
    "affirmed_span_mapped": (
        "info", "informational", _ENDS,
        "A reviewed historical note, bound to its exact text and sequence, maps through the record's exact "
        "gapped alignment to an affirmed interval of the coding sequence (not a physical coding end)."),
    "provider_functionality_unclassified": (
        "info", "unsupported_interpretation", (),
        "Per record export: the provider's functional values, kept raw.  The profile establishes no "
        "classification from them (a Boolean says neither functional, ORF nor pseudogene), so every "
        "functionality claim is unresolved."),
    "source_evidence_lineage": (
        "info", "informational", (),
        "All locked exports of one provider set are one evidence lineage: agreement between them (a record and "
        "its gapped alignment, or the two record exports) is consistency within that lineage, not independent "
        "confirmation.  Embedded support entries are kept as given, never retrieved or verified."),
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
        self.intervals: dict = {}        # record id -> CodingInterval of its candidate
        self.alignments: dict = {}       # candidate id -> (gapped-alignment mapping, gapped sequence)
        self.extensions: dict = {}       # (segment, label) -> extension between the record exports
        self.notes: dict = {}            # record id -> historical note outcome
        self.support: dict = {}          # record id -> [SupportReference]
        for reference in inventory.support:
            self.support.setdefault(reference.record.id, []).append(reference)
        self.lineage_id = _lineage(lock)
        self.observations = {o.id: o for o in inventory.observations}
        self.gapped: dict = {}
        for observation in inventory.observations:
            if observation.locator.role == "gapped_alignment":
                value = observation.value.get()
                self.gapped.setdefault(observation.locator.asset, {}).setdefault(value["header"], []).append(
                    (observation.locator, value["sequence"]))

    # -- reading ---------------------------------------------------------

    def value(self, locator: Locator) -> Any:
        observation = self.observations[locator.id]
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
        self.lineage()
        if self.scope_problem is not None:
            self.findings.add("coordinate_profile_out_of_scope",
                              locators=[Locator(a, EXPORT_ANCHOR) for a in self.exports],
                              observations={"profile": self.profile.id, "version": self.profile.version,
                                            "problem": self.scope_problem})
        for asset in self.exports:
            labels = Counter(r.label for r in self.parsed[asset])
            self.functionality(asset)
            for record in self.parsed[asset]:
                candidates = self.inventory.candidates_of(record.locator)
                if not candidates:
                    continue
                (candidate,) = candidates
                self.intervals[record.id] = self.coding_interval(record, candidate)
                self.boundaries(record, candidate)
                if record.segment == "V":
                    self.v_landmark(record, candidate, labels[record.label])
                elif record.segment == "J":
                    self.j_landmark(record, candidate)
        self.compare_extensions()
        self.historical_notes()
        for asset in self.exports:
            self.v_rss(asset)

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

    def boundaries(self, record, candidate) -> None:
        """Both coding boundaries, unresolved: no provider field or tag establishes a complete end."""
        for predicate in ("boundary_5", "boundary_3"):
            self.claim(candidate.locator, predicate, None, None, "stated", self.status(False),
                       [_locator(record, name) for name in _BOUNDARY_FIELDS[(record.segment, predicate)]])

    def functionality(self, asset: str) -> None:
        """An unresolved functionality claim per parsed record, and the raw values per export."""
        values = {"true": 0, "false": 0, "null": 0, "absent": 0, "other": 0}
        for record in self.inventory.records:
            if record.locator.asset != asset:
                continue
            at = _locator(record, "functional")
            observation = self.observations.get(at.id)    # none for a record that is not an object
            value = c.ABSENT if observation is None or observation.presence == "absent" else observation.value.get()
            values["absent" if value is c.ABSENT else "null" if value is None else
                   str(value).lower() if type(value) is bool else "other"] += 1
            if record.status == "parsed":
                self.claim(record.locator, "functionality", None, None, "stated", self.status(False), (at,))
        self.findings.add("provider_functionality_unclassified", locators=[Locator(asset, EXPORT_ANCHOR)],
                          observations={"asset": asset, "values": values, "classified": 0,
                                        "reason": dict(self.profile.unresolved)["functional_false_meaning"]})

    def lineage(self) -> None:
        assets = [a.asset_id for a in self.lock.assets if a.role != "data_terms"]
        self.findings.add("source_evidence_lineage", locators=[Locator(a, EXPORT_ANCHOR) for a in self.exports],
                          observations={"lineage": self.lineage_id, "provider": self.lock.provider,
                                        "assets": assets, "independent_lineages": 1,
                                        "support_references": len(self.inventory.support),
                                        "support_status": "embedded_unretrieved"})

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
                self.alignments[candidate.id] = (mapping, gapped)
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

    # -- endpoints: extensions, notes and declared RSS ----------------------

    def compare_extensions(self) -> None:
        """An extension claim and a context finding for every identity whose coding sequence in
        one record export extends the other's (from the endpoint comparison's findings)."""
        records = {r.id: r for r in self.inventory.records}
        for finding in self.inventory.findings:
            if finding.code != "export_coding_sequence_differs":
                continue
            observations = finding.observations.get()
            match = _EXTENDS.fullmatch(observations.get("relation") or "")
            if match is None:
                continue
            sides: dict = {}
            for subject in finding.subjects:
                if subject.id in records:
                    sides.setdefault(subject.asset, []).append(records[subject.id])
            extending = match.group(1)
            (other,) = [a for a in RECORD_EXPORT_FORMATS if a != extending]
            (longer,), (shorter,) = sides[extending], sides[other]
            self.extension(longer, shorter, match.group(2), observations["extra"])

    def extension(self, longer, shorter, side: str, extra: str) -> None:
        (long_candidate,) = self.inventory.candidates_of(longer.locator)
        (short_candidate,) = self.inventory.candidates_of(shorter.locator)
        (seq_at, sequence), (start_at, start), (end_at, end) = (self.field(longer, name) for name in
                                                                ("sequence", "gene_start", "gene_end"))
        interval = self.intervals[longer.id]
        bases = reason = None
        if interval.start is None:
            status, reason = "without_record_context", f"the declared coding interval is unusable: {interval.problem}"
        elif sequence[interval.start:interval.end] != short_candidate.sequence:
            status = "without_record_context"
            reason = (f"the declared coding interval of the {longer.locator.asset} record does not reproduce the "
                      f"{shorter.locator.asset} coding sequence, so no record position is proven for the extension")
        else:
            context = c.extension_context(sequence, (interval.start, interval.end), side, extra)
            status, bases = context.status, context.bases
        # A context ending (or starting) exactly where the record sequence does is only the record's
        # own terminal bases, not genomic flank around the extension.
        terminus = None
        if status != "without_record_context":     # judged on the record bases actually present
            terminus = (interval.end + len(bases) == len(sequence)) if side == "3" else (interval.start == len(bases))
        evidence = [long_candidate.locator, short_candidate.locator, seq_at, start_at, end_at]
        rs_at, rs = self.field(longer, "v_rs_start")
        versus_rss = None
        if longer.segment == "V":
            evidence.append(rs_at)
            if side == "3" and interval.start is not None:
                versus_rss = _extension_end(interval.end, extra, c.stated_position(rs).value)
        self.claim(long_candidate.locator, "extension",
                   {"side": side, "length": len(extra), "matches_record_context": status == "matches"}, None,
                   "variant_comparison", self.status(True, support="assertion"), evidence[:5])
        self.findings.add(_EXTENSION_CODES[status], locators=evidence, subjects=(
            longer.locator, long_candidate.locator, shorter.locator, short_candidate.locator), observations={
            "label": longer.label, "segment": longer.segment, "extending_export": longer.locator.asset,
            "side": side, "extension": extra, "record_context": bases, "reason": reason, "gene_start": _raw(start),
            "gene_end": _raw(end), "v_rs_start": _raw(rs) if longer.segment == "V" else None,
            "extension_end_vs_v_rss": versus_rss, "context_is_sequence_terminus": terminus,
            "lineage": self.lineage_id})
        self.extensions[(longer.segment, longer.label)] = {"export": longer.locator.asset, "side": side,
                                                           "bases": extra,
                                                           "matches_record_context": status == "matches"}

    def historical_notes(self) -> None:
        """One finding per profile note of this set, naming each record of its allele: an
        affirmed interval only with a reviewed interpretation, the exact note text and
        sequence, and an exact gapped alignment; otherwise unresolved with the reason.  Where
        the note states a column, each record's supplied span is compared with it through the
        record's exact gapped mapping."""
        for note in self.profile.historical_notes:
            if (note.set_id, note.release) != (self.lock.set_id, self.lock.release):
                continue
            outcomes: dict = {}      # (finding code, text matches) -> (record entries, subjects, locators)
            label, segment = note.label, None
            for asset in self.exports:
                for record in self.parsed[asset]:
                    candidates = self.inventory.candidates_of(record.locator)
                    if record.label != note.label or not candidates:
                        continue
                    (candidate,) = candidates
                    segment = record.segment
                    entry, locators, code = self.note_record(note, asset, record, candidate)
                    found = outcomes.setdefault((code, entry["text_matches"]), ([], [], []))
                    found[0].append(entry)
                    found[1].extend((record.locator, candidate.locator))
                    found[2].extend(locators)
            for (code, matches), (entries, subjects, locators) in outcomes.items():
                # The profile's reading (quotes, stated column) is the record's note only where the
                # record still holds the reviewed text; otherwise only that text's digest is named.
                observations = {"label": label, "segment": segment, "field": note.field, "numbering": note.numbering,
                                "text_matches": matches, "reviewed_text_sha256": note.text_sha256,
                                "quotes": list(note.quotes) if matches else [],
                                "stated_column": note.stated_column if matches else None}
                if code == "affirmed_span_mapped":
                    observations.update(first_affirmed_column=note.interpretation.first_affirmed_column,
                                        last_affirmed_column=note.interpretation.last_affirmed_column,
                                        coordinate_system="coding_sequence")
                observations["records"] = entries
                self.findings.add(code, subjects=subjects, locators=locators, observations=observations)

    def note_record(self, note, asset: str, record, candidate) -> tuple:
        """``(record entry, locators, finding code)`` of one record carrying a profile note."""
        text_at, text = self.field(record, note.field)
        matches = type(text) is str and text_digest(text) == note.text_sha256
        bound = dict(note.sequences).get(asset) == candidate.sha256
        mapping, gapped = self.alignments.get(candidate.id, (None, None))
        span = c.affirmed_span(note, asset, candidate.sequence, gapped)
        interval, reason = span.interval, span.reason
        if not matches:
            interval, reason = None, (f"the record's {note.field} text is not the note text the profile lists, so "
                                      "the note is not read and the profile's stated column is not applied")
        elif not bound:
            interval, reason = None, (f"the {asset} coding sequence is not the sequence this note was read against, "
                                      "so the stated column is not applied to it")
        versus = None
        if (matches and bound and note.stated_column is not None and mapping is not None
                and mapping.validation == "exact"):
            versus = c.span_against_column(mapping.blocks, note.stated_column)
        entry = {"asset": asset, "text_matches": matches, "sequence_bound": bound,
                 "gapped_mapping": None if mapping is None else mapping.validation,
                 "supplied_last_column": None if versus is None else versus.last_column,
                 "span_vs_stated_column": None if versus is None else versus.relation,
                 "bases_beyond_stated_column": None if versus is None else versus.bases_beyond,
                 "columns_short_of_stated_column": None if versus is None else versus.columns_short,
                 "affirmed_interval": None if interval is None else list(interval), "reason": reason}
        if matches:
            self.notes[record.id] = (text_at, {"field": note.field, "quotes": list(note.quotes),
                                               "affirmed_span": "unresolved" if interval is None else "mapped"})
        locators = [text_at] + ([mapping.output] if mapping is not None else [])
        return entry, locators, "affirmed_span_mapping_unresolved" if interval is None else "affirmed_span_mapped"

    def v_rss(self, asset: str) -> None:
        """Declared V RSS coordinates against the declared coding interval of one record export:
        the gap between the coding end and the RSS start, the bases the RSS starts with, and
        coordinates beyond the record sequence.  Nothing here is projected."""
        counts = {"adjacent": 0, "gaps": {}, "no_v_rs_start": [], "unusable": [], "starts_without_cac": [],
                  "starts_too_close_to_sequence_end": []}
        v_records = 0
        for record in self.parsed[asset]:
            if record.segment != "V":
                continue
            v_records += 1
            (rs_at, rs), (rs_end_at, rs_end) = self.field(record, "v_rs_start"), self.field(record, "v_rs_end")
            seq_at, sequence = self.field(record, "sequence")
            candidates = self.inventory.candidates_of(record.locator)
            if type(sequence) is str and candidates:
                for name, at, value in (("v_rs_start", rs_at, rs), ("v_rs_end", rs_end_at, rs_end)):
                    if type(value) is int and value > len(sequence):
                        self.finding("v_rss_coordinate_beyond_record_sequence", record, candidates[0],
                                     locators=(at, seq_at), field=name, value=value, sequence_length=len(sequence))
            interval = self.intervals.get(record.id)
            if interval is None or interval.start is None:
                counts["unusable"].append(record.label)
                continue
            position = c.v_rss_position(interval.start + 1, interval.end, rs)
            if position.gap is None or rs > len(sequence):
                counts["no_v_rs_start" if rs is None or rs is c.ABSENT else "unusable"].append(record.label)
                continue
            if position.gap == 0:
                counts["adjacent"] += 1
            else:
                counts["gaps"].setdefault(str(position.gap), []).append(record.label)
            (candidate,) = candidates
            at_start = sequence[rs - 1:rs + 6]
            if len(at_start) < len(_HEPTAMER_CORE):        # too few record bases to judge the core
                counts["starts_too_close_to_sequence_end"].append(record.label)
            elif not at_start.startswith(_HEPTAMER_CORE):
                counts["starts_without_cac"].append(record.label)
                self.finding("declared_v_rss_start_without_cac", record, candidate, locators=(rs_at, seq_at),
                             v_rs_start=rs, gap=position.gap, record_bases_at_v_rs_start=at_start)
            if position.position == "after_coding_end" and position.gap > 0:
                self.coding_end_before_rss(record, candidate, interval, rs, position.gap)
            elif position.position != "after_coding_end":
                self.finding(f"v_rss_start_{position.position}", record, candidate, locators=(
                    _locator(record, "gene_start"), _locator(record, "gene_end"), rs_at, rs_end_at, seq_at),
                    gene_start=interval.start + 1, gene_end=interval.end, v_rs_start=rs, v_rs_end=_raw(rs_end),
                    gap=position.gap, record_bases_at_v_rs_start=at_start,
                    record_bases_after_gene_end=sequence[interval.end:interval.end + 7])
        self.findings.add("v_rss_gap_after_coding_end", locators=[Locator(asset, EXPORT_ANCHOR)], observations={
            "asset": asset, "v_records": v_records, "adjacent": counts["adjacent"],
            "gaps": {gap: sorted(labels) for gap, labels in sorted(counts["gaps"].items(), key=lambda x: int(x[0]))},
            "no_v_rs_start": sorted(counts["no_v_rs_start"]), "unusable": sorted(counts["unusable"]),
            "starts_without_cac": sorted(counts["starts_without_cac"]),
            "starts_too_close_to_sequence_end": sorted(counts["starts_too_close_to_sequence_end"])})

    def coding_end_before_rss(self, record, candidate, interval, rs: int, gap: int) -> None:
        (seq_at, sequence), (end_at, _), (rs_at, _) = (self.field(record, name) for name in
                                                       ("sequence", "gene_end", "v_rs_start"))
        found = self.extensions.get((record.segment, record.label))
        extension = None
        if found is not None and found["side"] == "3":
            extension = {"export": found["export"], "bases": found["bases"],
                         "matches_record_context": found["matches_record_context"],
                         "reaches_rss_minus_1": _extension_end(interval.end, found["bases"], rs)
                         == "reaches_v_rss_minus_1"}
        support = [s for s in self.support.get(record.id, ()) if s.kind == "unrearranged"]
        note_at, note = self.notes.get(record.id, (None, None))
        locators = [end_at, rs_at, seq_at] + [s.locator for s in support] + ([note_at] if note_at else [])
        self.finding("coding_end_before_declared_v_rss", record, candidate, locators=locators,
                     gene_end=interval.end, v_rs_start=rs, gap=gap,
                     bases_through_rss_minus_1=sequence[interval.end:rs - 1], extension=extension,
                     historical_note=note, genomic_support=[
                         {"repository": s.repository, "reference": s.reference, "status": s.status}
                         for s in support], lineage=self.lineage_id)

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
    for name, implemented in IMPLEMENTED_CONVENTIONS.items():
        declared = getattr(profile, name)
        if declared != implemented:
            raise bad("/coordinate_profile", f"profile {profile.id} version {profile.version} declares {name} "
                                             f"{show(dict(declared) if isinstance(declared, Mapping) else declared)}; "
                                             f"this interpretation implements {show(implemented)}",
                      "unsupported_profile")
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


def _lineage(lock: SourceLock) -> str:
    return f"{lock.provider}:{lock.set_id}/{lock.release}"


def evidence_lineage(lock: SourceLock, locator: Locator) -> str:
    """The evidence lineage of a locator into ``lock``'s assets: every asset of one provider set
    is one lineage (``<provider>:<set id>/<release>``), so evidence roots in several of its
    exports never count as independent confirmations."""
    require_type(lock, SourceLock, "lock")
    require_type(locator, Locator, "locator")
    if locator.asset not in lock.asset_ids:
        raise bad("/asset", f"asset {locator.asset!r} is not in the lock (locked: {list(lock.asset_ids)})",
                  "identity_mismatch")
    return _lineage(lock)


def _raw(value: Any) -> Any:
    """A field value for a finding's observations (an absent field as null)."""
    return None if value is c.ABSENT else value


def _extension_end(gene_end: int, extension: str, v_rs_start: Optional[int]) -> Optional[str]:
    """Where a 3' extension ends relative to a declared V RSS start (1-based record positions)."""
    if v_rs_start is None:
        return None
    last = gene_end + len(extension)
    if last == v_rs_start - 1:
        return "reaches_v_rss_minus_1"
    return "before_v_rss" if last < v_rs_start - 1 else "overlaps_v_rss"
