"""Projection of an explicit recipe over a locked source into an unchanged runtime reference.

:func:`project` takes a source lock, its exact locked bytes and a parsed :class:`Recipe`.
It regenerates the interpreted source inventory from the bytes itself, reading them once (an
inventory supplied by a caller is never trusted), runs every cross-document check, resolves the ordered
selection against the inventory, and writes deterministic inputs for the existing local
importer:

* ``v.fasta``, ``d.fasta`` (VDJ only) and ``j.fasta``: for each selected allele, ``>`` and its
  literal label, then its supplied coding sequence (uppercase A/C/G/T, never repaired) in
  80-column lines; LF line ends, a trailing LF, UTF-8;
* ``annotations.json``: an annotation sidecar version 1, canonical JSON plus LF, with
  ``v_alignment`` ``none`` so that no alignment-derived assertion competes with a
  source-stated anchor.

The files are built into a runtime :class:`~GenAIRR._reference.artifact.Reference` with the
public ``Reference.from_fasta`` in a private temporary directory, whose name enters no output.
The complete :class:`~GenAIRR._reference_sources.report.BuildReport` gives one outcome per
explicit source record and per claim.

What reaches the runtime reference, per selected allele:

* an anchor only from an accepted coding-sequence anchor claim, at that claim's support (a
  source-supported one cites its exact source fields and the raw digest of their assets);
* a coding phase only from an accepted claim whose value is the projected anchor modulo 3,
  always as ``anchor_derived`` at the claim's own support, never as independent phase
  evidence (whether a source phase is independent of the landmark is not established);
* functionality and coding boundaries only from accepted claims (none under the current
  interpretation profile: the provider's Boolean ``functional`` is never mapped);
* a candidate claim only through an ``accept_claim`` decision, as an assertion citing that
  decision and its evidence; a conflicting claim never (it stays null and can only be
  withheld); ``withhold_claim`` withholds a claim;
* representable alias and paralog relationships, citing the fields that state them; a name
  that is not a selected allele stays a non-executable entity without sequence or weight;
* notes citing the source record (asset, raw digest, pointer), the build-input hash and the
  projection.

A selected allele carrying ``phase_landmark_projection_unsupported`` (a stated phase that is
not its stated landmark modulo 3) needs a ``withhold_claim`` decision covering its anchor and
coding-phase claims; neither value is ever shifted or replaced.

Everything the runtime reference cannot hold stays in the build report:
:data:`PROJECTION_LOSSES` and the claim and record outcomes say what and why.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass
import json
from pathlib import Path
import re
from types import MappingProxyType
from typing import Mapping, Optional

from GenAIRR._artifact_io import canonical_json
from GenAIRR._reference.artifact import Reference
from GenAIRR._reference.errors import ReferenceArtifactError
from GenAIRR._reference.importer import FASTA_MAX_BYTES, MAX_SEQUENCE as RUNTIME_MAX_SEQUENCE, _label_of
from GenAIRR._reference.model import check_label as runtime_label_problem
from GenAIRR._reference.sidecar import (
    MAX_BYTES as SIDECAR_MAX_BYTES, SCHEMA_NAME as SIDECAR_SCHEMA, SCHEMA_VERSION as SIDECAR_VERSION,
)
from GenAIRR._reference_sources import evidence, ogrdb
from GenAIRR._reference_sources._fields import (
    MAX_TEXT, bad, digest, document_bytes, private_directory, ptr, require_type,
)
from GenAIRR._reference_sources.common import RECORD_EXPORT_FORMATS, Generator, Locator, RuntimeHashes
from GenAIRR._reference_sources.consistency import (
    check_decisions, check_inventory, check_profile_scope, check_recipe, check_selection,
)
from GenAIRR._reference_sources.coordinates import RULES
from GenAIRR._reference_sources.digests import BUILD_INPUTS_DOMAIN, unchecked_build_inputs
from GenAIRR._reference_sources.errors import ReferenceSourceError
from GenAIRR._reference_sources.inventory import RUNTIME_PREDICATES, SourceInventory
from GenAIRR._reference_sources.lock import SourceLock
from GenAIRR._reference_sources.recipe import Recipe
from GenAIRR._reference_sources.report import (
    BuildHashes, BuildReport, ClaimOutcome, Counts, ExportCounts, OutcomeCounts, RecordOutcome, RequirementOutcome,
    RuntimeCapabilities, SegmentCounts, SelectedAllele,
)

# What this projection never puts into a runtime reference; all of it stays in the build report.
PROJECTION_LOSSES = (
    ("v_subregions", "The runtime reference has no V FWR/CDR subregion annotations: its V alignment is 'none', "
                     "so no subregion of a source gapped alignment survives projection. Simulated records therefore "
                     "carry no per-region V mutation counts (n_fwr1_mutations to n_fwr3_mutations stay 0) and "
                     "count every V mutation in n_v_unannotated_mutations."),
    ("gapped_sequences", "Source gapped sequences and their alignments to each coding sequence stay in the build "
                         "report; the runtime reference holds ungapped coding sequences only."),
    ("affirmed_spans", "Curation notes that limit an allele's affirmed span stay in the build report; no affirmed "
                       "interval is projected."),
    ("genomic_context", "Record sequences around the coding interval (flanks, RSS and leader coordinates, and the "
                        "coding extensions of the other record export) stay in the build report; the runtime "
                        "reference holds the selected coding sequence only."),
)

FASTA_MEMBERS = MappingProxyType({"V": "v.fasta", "D": "d.fasta", "J": "j.fasta"})
ANNOTATIONS_MEMBER = "annotations.json"
FASTA_LINE_WIDTH = 80
SIDECAR_CONVENTIONS = MappingProxyType({"coordinates": "zero_based_half_open", "v_alignment": "none"})

_SIDECAR_ANNOTATIONS = ("functionality", "anchor", "coding_phase", "boundary_5", "boundary_3")
_BASES = frozenset("ACGT")
_WITHHELD = {
    "unresolved": "unresolved: no usable value is established for it; the record's findings say why",
    "conflicting": "conflicting source claims; the competing values stay in the build report",
    "candidate": "a candidate only (support none); an accept_claim decision with evidence would project it as "
                 "an assertion",
}
# Unresolved claims whose predicate the interpretation profile never resolves, and why.
_UNRESOLVED = {
    "functionality": "unresolved: the provider states a functional value, but the interpretation profile cannot "
                     "classify it (true includes ORF-like alleles, and false is not established as pseudogene), so "
                     "the runtime functionality stays unknown",
    "boundary_5": "unresolved: no provider field or tag establishes a complete or partial 5' coding end under the "
                  "interpretation profile",
    "boundary_3": "unresolved: no provider field or tag establishes a complete or partial 3' coding end under the "
                  "interpretation profile",
}
_NOT_REPRESENTABLE = {
    "phase_independence": "the runtime reference holds a coding phase and its derivation, not a separate claim "
                          "about its independence of the landmark",
    "coding_interval": "the runtime reference holds the selected coding sequence, not its interval in the record "
                       "sequence",
    "extension": "the runtime reference holds the selected coding sequence, not the coding extensions of the other "
                 "record export",
}
_PHASE_AS_ANCHOR_DERIVED = ("projected with derivation anchor_derived: the {derivation} value equals the anchor "
                            "modulo 3, but whether it was derived independently of the landmark is unresolved, so it "
                            "is never projected as independent phase evidence")
_CONVERSIONS = {
    ("anchor", "V", "stated"): f"anchor = {RULES['v_anchor_coding0']} in the coding sequence, the cdr3_start bound "
                               f"to the record's exact gapped record ({RULES['v_cdr3_start_from_gapped']})",
    ("anchor", "J", "stated"): f"anchor = {RULES['j_anchor_coding0']} in the coding sequence",
    ("coding_phase", "V", "anchor_derived"): "coding phase = anchor % 3 (anchor-derived)",
    ("coding_phase", "J", "stated"): f"coding phase = {RULES['j_reported_phase']}, equal to anchor % 3; projected "
                                     "as anchor-derived because its independence of the landmark is unresolved",
}
_BOUNDARY_NEEDS = {"V": ("boundary_3",), "D": ("boundary_5", "boundary_3"), "J": ("boundary_5",)}
# A relationship name stated by a record: ``/aliases/<n>`` or ``/paralogs/<n>`` below the record.
_RELATION_FIELD = re.compile(r"/(aliases|paralogs)/(?:0|[1-9][0-9]*)")
_RELATION_KINDS = {"aliases": "alias", "paralogs": "paralog"}


@dataclass(frozen=True)
class Projection:
    """A projected build: the exact input files, the runtime reference and the build report.

    ``fasta`` is ``((segment, bytes), ...)`` in V, (D,) J order and ``annotations`` the sidecar
    bytes; ``build_input_hash`` is the identity the runtime citations name.
    """

    build_input_hash: str
    fasta: tuple
    annotations: bytes
    reference: Reference
    report: BuildReport

    def fasta_of(self, segment: str) -> bytes:
        for name, data in self.fasta:
            if name == segment:
                return data
        raise KeyError(segment)

    @property
    def inputs(self) -> dict:
        """A fresh ``{member name: exact bytes}`` of every projected input file."""
        files = {FASTA_MEMBERS[segment]: data for segment, data in self.fasta}
        files[ANNOTATIONS_MEMBER] = self.annotations
        return files


def _under(prefix: str, check, *args) -> None:
    try:
        check(*args)
    except ReferenceSourceError as exc:
        raise exc.under(prefix) from None


def _bounded(parts) -> Optional[str]:
    """Detail items joined with ``"; "``, or ``None`` for none.  A report detail holds at most
    ``MAX_TEXT`` characters: beyond that, whole leading items are kept and an explicit marker says
    how many were omitted.  Nothing is ever cut silently."""
    parts = list(parts)
    if not parts:
        return None
    text = "; ".join(parts)
    keep = len(parts)
    while len(text) > MAX_TEXT:
        keep -= 1
        marker = f"... and {len(parts) - keep} more not shown (a detail holds at most {MAX_TEXT} characters)"
        text = "; ".join(parts[:keep] + [marker])
    return text


def _require_record_exports(lock: SourceLock) -> None:
    """A build compares both record exports, so both must be locked."""
    absent = [f for f in RECORD_EXPORT_FORMATS if f not in lock.asset_ids]
    if absent:
        recorded = {a.format: a for a in lock.missing_assets}
        why = "; ".join(f"{f}: {recorded[f].reason}" + (f" ({recorded[f].detail})" if recorded[f].detail else "")
                        for f in absent if f in recorded)
        raise bad("/source_lock/assets", f"a build compares both record exports {list(RECORD_EXPORT_FORMATS)}; the "
                                         f"lock of {lock.set_id} release {lock.release} lacks {absent}"
                                         + (f" ({why})" if why else ""), "missing_asset")


def _sequence_problem(value) -> Optional[str]:
    """Why a supplied coding sequence cannot be a runtime allele's sequence, or ``None``."""
    if type(value) is not str:
        return "its coding_sequence is not text"
    if not 1 <= len(value) <= RUNTIME_MAX_SEQUENCE:
        return f"its coding sequence has {len(value)} bases; a runtime allele holds 1 to {RUNTIME_MAX_SEQUENCE}"
    for position, char in enumerate(value):
        if char not in _BASES:
            return (f"its coding sequence holds {char!r} at position {position}; a runtime allele holds only "
                    "uppercase A, C, G and T")
    return None


def _label_problem(label) -> Optional[str]:
    """Why a label cannot be a runtime allele's label, written as a FASTA header, or ``None``."""
    problem = runtime_label_problem(label)
    if problem is not None:
        return f"label {label!r} {problem}"
    read = _label_of(label)
    if read != label:
        return f"a FASTA header >{label} is read back as the label {read!r}, so the label would not be kept"
    return None


def _supplied_sequence_problem(inventory: SourceInventory, record) -> str:
    """Why a parsed record supplies no candidate sequence, from its observed coding_sequence."""
    try:
        observation = inventory.observation(Locator(record.locator.asset,
                                                    ptr(record.locator.pointer, "coding_sequence")))
    except KeyError:
        return "the record has no coding_sequence field"
    if observation.presence == "absent":
        return "the record has no coding_sequence"
    return _sequence_problem(observation.value.get()) or "the record supplies no coding sequence candidate"


def _check_supplied_sequences(recipe: Recipe, inventory: SourceInventory) -> None:
    """A selected record whose coding sequence the runtime cannot hold is refused in Python,
    before any native work, naming the allele (the inventory forms no candidate for it)."""
    reference = recipe.reference
    for i, entry in enumerate(recipe.selection):
        try:
            record = inventory.record(entry.record)
        except KeyError:
            continue
        if (record.status != "parsed" or (record.segment, record.label) != (entry.segment, entry.label)
                or (record.taxonomy_id, record.locus) != (reference.taxon.taxonomy_id, reference.locus)
                or inventory.candidates_of(record.locator)):
            continue                # check_selection decides these
        raise bad(ptr("/selection", i), f"{entry.allele_id}: record {record.id} cannot be projected: "
                                        f"{_supplied_sequence_problem(inventory, record)}", "projection")


def _noncanonical(segment: str, observed: dict) -> str:
    """A non-canonical anchor codon, from its ``noncanonical_landmark_codon`` finding."""
    position = "IMGT 104" if segment == "V" else "IMGT 118"
    return (f"its {position} codon {observed.get('codon')} encodes {observed.get('residue')}, not "
            f"{' or '.join(observed.get('expected_residues') or ())} (at coding position {observed.get('anchor')}); "
            "the anchor is present, but its codon is not canonical")


def _relation(observed: Optional[dict]) -> str:
    """How the record exports' coding sequences relate, from an ``export_coding_sequence_differs`` finding."""
    relation = (observed or {}).get("relation")
    extra = (observed or {}).get("extra") or ""
    parts = relation.split("_extends_") if relation else ()
    if len(parts) == 2:
        other = "airr" if parts[0] == "airr_ex" else "airr_ex"
        return (f"the {parts[0]} coding sequence extends the {other} one at the {parts[1][0]}' end by {len(extra)} "
                f"nt ({extra})")
    if relation == "different":
        return "the airr and airr_ex coding sequences differ, neither extending the other"
    if relation == "not_comparable":
        return "the record exports cannot be compared one to one (several records or no candidate on a side)"
    return "alternatives within one record export"


def _require_segments(recipe: Recipe) -> None:
    present = {entry.segment for entry in recipe.selection}
    absent = [s for s in recipe.reference.identity.segments if s not in present]
    if absent:
        raise bad("/selection", f"the recipe selects no {' allele and no '.join(absent)} allele; a "
                                f"{recipe.reference.chain_type} reference needs at least one allele of each of "
                                f"{list(recipe.reference.identity.segments)}", "unresolved_selection")


# ── resolution of the recipe against the inventory ───────────────────


class _Build:
    """The selection resolved against the inventory, and everything projected from it."""

    def __init__(self, lock: SourceLock, inventory: SourceInventory, recipe: Recipe, build_input_hash: str):
        self.lock, self.inventory, self.recipe = lock, inventory, recipe
        self.build_input_hash = build_input_hash
        self.segments = lock.identity.segments
        self.assets = {a.asset_id: a for a in lock.assets}
        self.entries = recipe.selection
        self.chosen = []
        for entry in self.entries:
            (candidate,) = [c for c in inventory.candidates_of(entry.record) if c.sha256 == entry.sequence_sha256]
            self.chosen.append(candidate)
        # Runtime order: V, then D, then J, each in recipe order (the order of the FASTA files).
        self.order = [i for segment in self.segments for i, e in enumerate(self.entries) if e.segment == segment]
        self.by_subject = {}
        for i, entry in enumerate(self.entries):
            self.by_subject[entry.record.id] = i
            self.by_subject[self.chosen[i].id] = i
        self.identities = {(e.segment, e.label): i for i, e in enumerate(self.entries)}
        self.decided = {claim: (i, d) for i, d in enumerate(recipe.decisions) for claim in d.claims}
        self.claims_of: dict = {}
        for claim in inventory.claims:
            self.claims_of.setdefault(claim.subject.id, []).append(claim)

    def entry_claims(self, i: int) -> list:
        """Every claim about selection entry ``i``: about its record or its chosen candidate."""
        return self.claims_of.get(self.entries[i].record.id, []) + self.claims_of.get(self.chosen[i].id, [])

    # -- citations ------------------------------------------------------

    def cite(self, locators) -> str:
        """``<asset> <raw sha256> <pointer>, <pointer>; ...`` in canonical locator order."""
        groups: dict = {}
        for locator in sorted(set(locators), key=Locator.sort_key):
            groups.setdefault(locator.asset, []).append(locator.pointer)
        return "; ".join(f"{asset} {self.assets[asset].sha256} {', '.join(pointers)}"
                         for asset, pointers in groups.items())

    def method(self, claim, segment: str, accepted_by: Optional[int]) -> str:
        adapter, projection = self.lock.adapter, self.recipe.projection
        profile = adapter.coordinate_profile
        rule = _CONVERSIONS.get((claim.predicate, segment, claim.derivation))
        if rule is None and claim.predicate == "coding_phase":
            rule = (f"coding phase ({claim.derivation}) equal to anchor % 3, projected as anchor_derived because its "
                    "independence of the landmark is unresolved")
        elif rule is None:
            rule = f"{claim.predicate} ({claim.derivation})"
        text = f"{adapter.id} v{adapter.version} with {profile.id} v{profile.version}: {rule}"
        if accepted_by is not None:
            text += f"; accepted by recipe decision {accepted_by}"
        return f"{text}; projection {projection.id} v{projection.version}"

    def notes(self, i: int) -> list:
        record, projection = self.entries[i].record, self.recipe.projection
        return [f"source record {record.asset} {self.assets[record.asset].sha256} {record.pointer}",
                f"build inputs {self.build_input_hash}, projection {projection.id} version {projection.version}"]

    # -- checks -----------------------------------------------------------

    def check_projectable(self) -> None:
        """Every selected label survives a FASTA header and every sequence is one the runtime holds,
        checked here in Python so that nothing unrepresentable reaches the native layer."""
        for i in self.order:
            entry = self.entries[i]
            problem = _label_problem(entry.label) or _sequence_problem(self.chosen[i].sequence)
            if problem is not None:
                raise bad(ptr("/selection", i), f"{entry.allele_id} cannot be projected: {problem}", "projection")

    def check_frameshifts(self) -> None:
        """A selected allele whose stated phase is not its stated landmark modulo 3 needs both
        observations withheld: the runtime holds one phase, the anchor modulo 3."""
        withheld = {c for d in self.recipe.decisions if d.action == "withhold_claim" for c in d.claims}
        for finding in self.inventory.findings:
            if finding.code != "phase_landmark_projection_unsupported":
                continue
            observed = finding.observations.get()
            for i in sorted({self.by_subject[s.id] for s in finding.subjects if s.id in self.by_subject}):
                claims = [c for c in self.entry_claims(i) if c.predicate in ("anchor", "coding_phase")]
                open_ = [c.id for c in claims if c.id not in withheld]
                if open_:
                    raise bad(ptr("/selection", i), (
                        f"{self.entries[i].allele_id} carries finding {finding.id} (phase_landmark_projection_"
                        f"unsupported: the stated phase {observed.get('phase')} is not the stated anchor "
                        f"{observed.get('anchor')} modulo 3, {observed.get('anchor_modulo_3')}); the runtime "
                        "reference holds a phase only as its anchor modulo 3, and neither source value is shifted "
                        "or replaced. Exclude the allele from the selection, or withhold its anchor and coding_phase "
                        f"claims with a withhold_claim decision (not yet withheld: {', '.join(open_)})"),
                        "projection")

    # -- claims -----------------------------------------------------------

    def claims(self) -> tuple:
        """One outcome per inventory claim (inventory order); fills ``self.annotations``."""
        outcomes, slots, anchors, phases = {}, {}, {}, []
        for claim in self.inventory.claims:
            if claim.subject.id not in self.by_subject:
                outcomes[claim.id] = ClaimOutcome(claim.id, "not_selected", None, None)
            elif claim.predicate not in RUNTIME_PREDICATES:
                outcomes[claim.id] = ClaimOutcome(claim.id, "not_representable", None,
                                                  _NOT_REPRESENTABLE.get(claim.predicate))
            elif claim.predicate == "coding_phase":
                phases.append(claim)            # judged against the projected anchor
            else:
                outcomes[claim.id] = self._runtime_claim(claim, slots, anchors)
        for claim in phases:
            outcomes[claim.id] = self._runtime_claim(claim, slots, anchors)
        self.annotations, weakened = self._annotations(slots)
        for claim_id, why in weakened.items():
            outcome = outcomes[claim_id]
            reason = why if outcome.reason is None else f"{outcome.reason}; {why}"
            outcomes[claim_id] = ClaimOutcome(claim_id, outcome.outcome, outcome.decision, reason)
        return tuple(outcomes[c.id] for c in self.inventory.claims)

    def _runtime_claim(self, claim, slots: dict, anchors: dict) -> ClaimOutcome:
        i = self.by_subject[claim.subject.id]
        entry = self.entries[i]
        index, decision = self.decided.get(claim.id, (None, None))
        if decision is not None and decision.action == "withhold_claim":
            return ClaimOutcome(claim.id, "withheld", index, f"withheld by recipe decision {index}")
        accepted_by = index if decision is not None and decision.action == "accept_claim" else None
        if accepted_by is not None and claim.status == "conflicting":
            values = sorted({canonical_json(c.value.get()) for c in self.entry_claims(i)
                             if (c.predicate, c.coordinate_system, c.status)
                             == (claim.predicate, claim.coordinate_system, "conflicting")})
            raise bad(ptr("/decisions", accepted_by), (
                f"{entry.allele_id}: accept_claim names {claim.predicate} claim {claim.id}, which is conflicting: "
                f"the source claims about this {claim.predicate} disagree (values {', '.join(values)}), and a "
                "conflicting annotation can only be withheld (withhold_claim), never projected"), "projection")
        if claim.status != "accepted" and accepted_by is None:
            reason = _UNRESOLVED.get(claim.predicate) if claim.status == "unresolved" else None
            return ClaimOutcome(claim.id, "withheld", None, reason or _WITHHELD[claim.status])
        value = claim.value.get()
        problem = self._gate(claim, value, anchors.get(i))
        if problem is not None:
            if accepted_by is not None:
                raise bad(ptr("/decisions", accepted_by), f"{entry.allele_id}: accept_claim names {claim.predicate} "
                                                          f"claim {claim.id}, which cannot be projected: {problem}",
                          "projection")
            return ClaimOutcome(claim.id, "withheld", None, problem)
        # An anchor-derived phase also rests on the fields of the anchor it is derived from.
        resting_on = anchors[i][1] if claim.predicate == "coding_phase" else ()
        if claim.status == "accepted":
            support, source = claim.support, self.cite(claim.evidence + resting_on)
        else:                               # a candidate named by an accept_claim decision
            support = "assertion"
            source = f"recipe decision {accepted_by}: {self.cite(decision.evidence + resting_on)}"
        if claim.predicate == "anchor":
            anchors[i] = (value, claim.evidence)
        slots.setdefault((i, claim.predicate), []).append(
            (claim, value, support, source, self.method(claim, entry.segment, accepted_by)))
        reason = None
        if claim.predicate == "coding_phase" and claim.derivation != "anchor_derived":
            reason = _PHASE_AS_ANCHOR_DERIVED.format(derivation=claim.derivation)
        return ClaimOutcome(claim.id, "projected", accepted_by, reason)

    @staticmethod
    def _gate(claim, value, anchor: Optional[tuple]) -> Optional[str]:
        """Why a claim with a value cannot become a runtime annotation, or ``None``; ``anchor`` is the
        ``(value, evidence)`` of the anchor projected for the same allele, if any."""
        if claim.predicate == "v_subregions":
            return "the runtime projection holds no V subregions (its V alignment is 'none')"
        if claim.predicate in ("anchor", "coding_phase") and claim.coordinate_system != "coding_sequence":
            return (f"stated in {claim.coordinate_system} coordinates; the runtime reference takes coding-sequence "
                    "coordinates")
        if claim.predicate == "coding_phase":
            if anchor is None:
                return "no anchor of this sequence is projected, and the runtime reference holds a phase only as " \
                       "the anchor modulo 3"
            if value != anchor[0] % 3:
                return f"phase {value} is not the projected anchor {anchor[0]} modulo 3 ({anchor[0] % 3})"
        return None

    def _annotations(self, slots: dict) -> tuple:
        """``({entry index: {predicate: sidecar annotation}}, {claim id: why its support was lowered})``
        from the projected claims."""
        annotations: dict = {}
        weakened: dict = {}
        for (i, predicate), items in slots.items():
            values = {canonical_json(value) for _, value, _, _, _ in items}
            if len(values) > 1:
                raise bad(ptr("/selection", i), f"{self.entries[i].allele_id}: claims "
                                                f"{[c.id for c, *_ in items]} would project different {predicate} "
                                                f"values {sorted(values)}", "projection")
            supports = {support for _, _, support, _, _ in items}     # the weaker one when they differ
            annotation = {"value": items[0][1],
                          "support": supports.pop() if len(supports) == 1 else "assertion",
                          "source": "; ".join(dict.fromkeys(source for *_, source, _ in items)),
                          "method": "; ".join(dict.fromkeys(method for *_, method in items))}
            if predicate == "coding_phase":
                annotation["derivation"] = "anchor_derived"
            annotations.setdefault(i, {})[predicate] = annotation
            weaker = [c.id for c, _, support, _, _ in items if support == annotation["support"]]
            for claim, _, support, _, _ in items:
                if support != annotation["support"]:
                    weakened[claim.id] = (f"projected at support {annotation['support']}, weaker than its own "
                                          f"{support} support, because {', '.join(weaker)} projected into the same "
                                          f"{predicate} annotation at support {annotation['support']}")
        return annotations, weakened

    # -- relationships ----------------------------------------------------

    def relationships(self) -> list:
        """Sidecar relationships for the aliases and paralogs the selected records state.

        A name that is a selected allele of the same segment is that allele (both executable);
        any other name is a non-executable entity, with no sequence, support or usage weight.
        A relationship the runtime reference cannot hold (a name equal to the allele's own label,
        or an alias of an executable label of the same segment) is not projected; the selected
        record's outcome says so.
        """
        stated: dict = {}
        for observation in self.inventory.observations:
            record = observation.record
            if record is None or record.id not in self.by_subject:
                continue
            match = _RELATION_FIELD.fullmatch(observation.locator.pointer[len(record.pointer):])
            if match is None or (observation.presence, observation.interpretation) != ("present", "interpreted"):
                continue
            name = observation.value.get()
            if runtime_label_problem(name) is None:
                stated.setdefault(record.id, {}).setdefault((_RELATION_KINDS[match.group(1)], name), []).append(
                    observation.locator)
        relationships, self.unprojected = [], {}
        for i in self.order:
            entry = self.entries[i]
            for (kind, name), locators in stated.get(entry.record.id, {}).items():
                if name == entry.label:
                    problem = f"{kind} {name} names the allele's own label"
                elif kind == "alias" and (entry.segment, name) in self.identities:
                    problem = (f"alias {name} names an executable {entry.segment} allele of this reference, and "
                               "an alias is another name of the same allele")
                else:
                    relationships.append({"kind": kind, "allele": entry.allele_id, "entity": name,
                                          "source": self.cite(locators)})
                    continue
                self.unprojected.setdefault(entry.record.id, []).append(
                    f"relationship not projected: {problem} ({self.cite(locators)})")
        return relationships

    # -- records ----------------------------------------------------------

    def records(self) -> tuple:
        """One outcome per inventory record (inventory order)."""
        inventory = self.inventory
        sequences: dict = {}
        for record in inventory.records:
            if record.status == "parsed":
                sequences.setdefault((record.segment, record.label), set()).update(
                    c.sha256 for c in inventory.candidates_of(record.locator))
        self.noncanonical = {s.id: f.observations.get() for f in inventory.findings
                             if f.code == "noncanonical_landmark_codon" for s in f.subjects}
        self.relations = {(o["segment"], o["label"]): o for o in (f.observations.get() for f in inventory.findings
                                                                  if f.code == "export_coding_sequence_differs")}
        outcomes = []
        for record in inventory.records:
            if record.status != "parsed":
                outcomes.append(RecordOutcome(record.locator, "rejected", None, None, None, record.problem))
            elif record.id in self.by_subject:
                entry = self.entries[self.by_subject[record.id]]
                detail = _bounded(self.unprojected.get(record.id, ()))
                outcomes.append(RecordOutcome(record.locator, "selected", entry.allele_id, None, None, detail))
            else:
                exclusion, related, details = self._exclusion(record, sequences)
                outcomes.append(RecordOutcome(record.locator, "excluded", None, exclusion, related,
                                              _bounded(details)))
        return tuple(outcomes)

    def _exclusion(self, record, sequences: dict) -> tuple:
        """``(reason, related selected record, detail items)`` for an unselected parsed record; first
        rule wins."""
        candidates = self.inventory.candidates_of(record.locator)
        key = (record.segment, record.label)
        if key in self.identities:
            entry = self.entries[self.identities[key]]
            if any(c.sha256 == entry.sequence_sha256 for c in candidates):
                return ("duplicate_occurrence", entry.record,
                        [f"another occurrence of {entry.allele_id}, selected at {entry.record.id}; it supplies the "
                         "same coding sequence"])
            supplied = "a different coding sequence" if candidates else "no coding sequence"
            return ("alternative_for_selected_identity", entry.record,
                    [f"an alternative record of {entry.allele_id}, selected at {entry.record.id}; it supplies "
                     f"{supplied}"])
        problem = self._structural_problem(record, candidates)
        if problem is not None:
            return "structural_incompatibility", None, [problem]
        if len(sequences.get(key, ())) > 1:
            return ("unresolved_conflicting_sequence", None,
                    [f"{record.label} has {len(sequences[key])} different coding sequences across its records "
                     f"({_relation(self.relations.get(key))}), and the recipe selects none of them"])
        unmet = self._missing_evidence(record, candidates)
        if unmet:
            return "missing_required_evidence", None, unmet
        return "outside_explicit_selection", None, ()

    def _structural_problem(self, record, candidates) -> Optional[str]:
        """Why the record could not be projected at all, or ``None``."""
        if record.segment not in self.segments:
            return f"segment {record.segment} is not part of a {self.lock.identity.chain_type} reference"
        problem = _label_problem(record.label)
        if problem is None and not candidates:
            problem = _supplied_sequence_problem(self.inventory, record)
        return problem

    def _missing_evidence(self, record, candidates) -> list:
        """Each recipe requirement this record's own claims cannot meet, and why."""
        (candidate,) = candidates
        claims = self.claims_of.get(candidate.id, ())

        def accepted(predicate, value=None):
            return [c for c in claims if c.predicate == predicate and c.status == "accepted"
                    and (value is None or c.value.get() == value)]

        anchors = [c for c in accepted("anchor") if c.coordinate_system == "coding_sequence"]
        problems = []
        for requirement in self.recipe.requires:
            capability, level, problem = requirement.capability, requirement.evidence, None
            if capability in ("source_annotated_junction", "canonical_productive") and record.segment != "D":
                if not anchors:
                    problem = "no accepted anchor"
                elif capability == "source_annotated_junction" and any(c.support != "source" for c in anchors):
                    problem = "its anchor is not source-supported"
                elif capability == "canonical_productive" and candidate.id in self.noncanonical:
                    problem = _noncanonical(record.segment, self.noncanonical[candidate.id])
                elif capability == "canonical_productive" and level == "source":
                    problem = ("a source-certified canonical_productive needs an independently source-supported "
                               "coding phase, and this projection projects phases as anchor-derived only")
            elif capability in ("original_coding_ends", "complete_reference_span"):
                needs = _BOUNDARY_NEEDS[record.segment] if capability == "original_coding_ends" else \
                    ("boundary_5", "boundary_3")
                found = [accepted(p, "complete") for p in needs]
                if not all(found):
                    problem = "no accepted complete " + " and ".join(p for p, f in zip(needs, found) if not f)
                elif level == "source" and any(c.support != "source" for f in found for c in f):
                    problem = "its complete boundaries are not source-supported"
            if problem is not None:
                problems.append(f"{capability} at evidence {level}: {problem}")
        return problems

    # -- files ------------------------------------------------------------

    def fasta(self) -> tuple:
        files = []
        for segment in self.segments:
            lines = []
            for i in self.order:
                if self.entries[i].segment != segment:
                    continue
                sequence = self.chosen[i].sequence
                lines.append(">" + self.entries[i].label)
                lines.extend(sequence[k:k + FASTA_LINE_WIDTH] for k in range(0, len(sequence), FASTA_LINE_WIDTH))
            files.append((segment, ("\n".join(lines) + "\n").encode("utf-8")))
        return tuple(files)

    def sidecar(self, relationships: list) -> bytes:
        alleles = []
        for i in self.order:
            entry, annotations = self.entries[i], self.annotations.get(i, {})
            allele = {"segment": entry.segment, "label": entry.label, "sequence_sha256": entry.sequence_sha256,
                      "gene": None}
            for name in _SIDECAR_ANNOTATIONS:
                allele[name] = annotations.get(name)
            allele["notes"] = self.notes(i)
            alleles.append(allele)
        document = {"schema": SIDECAR_SCHEMA, "schema_version": SIDECAR_VERSION,
                    "reference": self.recipe.reference.to_dict(), "conventions": dict(SIDECAR_CONVENTIONS),
                    "alleles": alleles, "relationships": relationships}
        return document_bytes(document)


# ── the runtime reference ────────────────────────────────────────────


def _check_sizes(fasta: tuple, annotations: bytes) -> None:
    """The projected files fit the runtime importer's byte limits."""
    for segment, data in fasta:
        if len(data) > FASTA_MAX_BYTES:
            raise bad("/selection", f"{FASTA_MEMBERS[segment]} would hold {len(data)} bytes; the runtime importer "
                                    f"reads at most {FASTA_MAX_BYTES}", "resource_limit")
    if len(annotations) > SIDECAR_MAX_BYTES:
        raise bad("/selection", f"{ANNOTATIONS_MEMBER} would hold {len(annotations)} bytes; the runtime importer "
                                f"reads at most {SIDECAR_MAX_BYTES}", "resource_limit")


def _runtime(fasta: tuple, annotations: bytes) -> Reference:
    """Build through the public ``Reference.from_fasta`` from files in a private temporary directory."""
    with private_directory("genairr_projection_", "the projected files") as root:
        paths = {}
        try:
            for segment, data in fasta:
                paths[segment] = Path(root, FASTA_MEMBERS[segment])
                paths[segment].write_bytes(data)
            sidecar = Path(root, ANNOTATIONS_MEMBER)
            sidecar.write_bytes(annotations)
        except OSError as exc:
            raise bad("", f"cannot write the projected files into a private temporary directory: {exc}",
                      "io") from exc
        try:
            return Reference.from_fasta(v=paths["V"], j=paths["J"], d=paths.get("D"), annotations=sidecar)
        except ReferenceArtifactError as exc:
            raise bad("", f"the runtime reference rejected the projected files ({exc.kind} at "
                          f"{exc.path or 'the whole input'}): {exc.message}", "projection") from exc


def _check_runtime(reference: Reference, build: _Build, document: dict) -> None:
    """The runtime alleles are exactly the selection, and the defaults are the recipe's."""
    expected = [build.entries[i] for i in build.order]
    alleles = reference.alleles
    for n, entry in enumerate(expected):
        allele = alleles[n] if n < len(alleles) else None
        found = None if allele is None else (allele.segment, allele.label, allele.sequence_sha256)
        if found != (entry.segment, entry.label, entry.sequence_sha256):
            raise bad(ptr("/selection", build.order[n]), f"runtime allele {n} is {found}, not {entry.allele_id} "
                                                         f"with {entry.sequence_sha256}", "projection")
    if len(alleles) != len(expected):
        raise bad("/selection", f"the runtime reference holds {len(alleles)} alleles for {len(expected)} selected",
                  "projection")
    profile = document["models"]["profile"]
    if profile != build.recipe.defaults_profile:
        raise bad("/defaults_profile", f"the runtime reference uses defaults profile {profile!r}, not the recipe's "
                                       f"{build.recipe.defaults_profile!r}", "projection")


_MAX_CAPABILITY_MESSAGE = MAX_TEXT


def _check_requirements(recipe: Recipe, capabilities) -> None:
    """Each recipe requirement is satisfied by the built reference's capabilities: available at an
    evidence level at least the required one.  One error names every unmet requirement (the report
    is not returned on failure); its path is the first unmet one.

    Every unmet requirement and every group of its blockers (those sharing a requirement and a
    hint) is always stated, with the group's count; there are at most as many requirements as
    capabilities.  Only the alleles of each group are bounded: the groups share what is left of
    the message's budget, and each group ends with its own "... and N more" marker when it names
    fewer alleles than it counts.
    """
    unmet = []                                  # (index, header, [(group header, alleles)])
    for i, requirement in enumerate(recipe.requires):
        if requirement.capability in capabilities.available(evidence=requirement.evidence):
            continue
        status = capabilities.status(requirement.capability)
        header = (f"{ptr('/requires', i)}: required capability {requirement.capability} at evidence "
                  f"{requirement.evidence} is not available on the built reference")
        if status.available:
            header += f"; it is available only at evidence {status.evidence}"
        if status.note:
            header += f"; {status.note}"
        groups: dict = {}
        for blocker in status.unmet:
            groups.setdefault((blocker.requirement, blocker.hint), []).append(blocker.allele)
        unmet.append((i, header, [(f"unmet for {len(alleles)} alleles: {needed} (resolved by: {hint})", alleles)
                                  for (needed, hint), alleles in groups.items()]))
    if not unmet:
        return
    lead = f"{len(unmet)} of {len(recipe.requires)} required capabilities are not met: "
    marker = " ... and {} more not shown (the message is bounded)"
    fixed = len(lead) + sum(len(header) + 3 + sum(len(group) + 2 + len(marker) + 8 for group, _ in groups)
                            for _, header, groups in unmet)
    count = sum(len(groups) for _, _, groups in unmet)
    budget = max(0, (_MAX_CAPABILITY_MESSAGE - fixed) // count) if count else 0
    parts = []
    for _, header, groups in unmet:
        text = header
        for group, alleles in groups:
            named, used = [], 0
            for allele in alleles:
                if used + len(allele) + 2 > budget:
                    break
                named.append(allele)
                used += len(allele) + 2
            text += f"; {group}: {', '.join(named)}"
            if len(named) < len(alleles):
                text += marker.format(len(alleles) - len(named))
        parts.append(text)
    raise bad(ptr("/requires", unmet[0][0]), lead + " | ".join(parts), "capability")


# ── the entry point ──────────────────────────────────────────────────


def project(lock: SourceLock, assets: Mapping[str, bytes], recipe: Recipe) -> Projection:
    """Build ``recipe`` over ``lock`` and its exact locked bytes into a runtime reference.

    ``assets`` maps each locked asset id to its bytes.  Errors are
    :class:`~GenAIRR._reference_sources.errors.ReferenceSourceError`; their paths are
    ``/source_lock/...`` for the lock, ``/assets/<asset id>...`` for the bytes,
    ``/inventory/...`` for the regenerated inventory, recipe paths (``/selection/<i>``,
    ``/decisions/<i>``, ``/requires/<i>``, ``/defaults_profile``) for everything the recipe
    decides, and ``/report/...`` for the build report's own re-checks.  Arguments of the wrong type
    raise ``TypeError``.
    """
    require_type(lock, SourceLock, "lock")
    require_type(recipe, Recipe, "recipe")
    _require_record_exports(lock)
    _under("/source_lock", check_profile_scope, lock)
    check_recipe(recipe, lock)
    inventory = evidence._interpret_locked(lock, assets)
    _under("/inventory", check_inventory, inventory, lock)
    _check_supplied_sequences(recipe, inventory)
    check_selection(recipe, inventory)
    check_decisions(recipe, inventory)
    _require_segments(recipe)
    inputs = unchecked_build_inputs(lock, inventory, recipe)
    build_input_hash = digest(inputs, BUILD_INPUTS_DOMAIN)

    build = _Build(lock, inventory, recipe, build_input_hash)
    build.check_projectable()
    build.check_frameshifts()
    claims = build.claims()
    relationships = build.relationships()
    records = build.records()
    fasta = build.fasta()
    annotations = build.sidecar(relationships)
    _check_sizes(fasta, annotations)

    reference = _runtime(fasta, annotations)
    document = json.loads(reference.to_json())
    _check_runtime(reference, build, document)
    _check_requirements(recipe, reference.capabilities())
    try:
        report = _report(build, inputs, records, claims, reference, Generator.from_dict(document["generator"]))
    except ReferenceSourceError as exc:          # the report's own re-checks, located in the report
        raise exc.under("/report") from None
    return Projection(build_input_hash, fasta, annotations, reference, report)


def _runtime_accounting(recipe: Recipe, reference: Reference) -> dict:
    """The report fields that describe the runtime reference: its hashes, its capabilities and
    whether each recipe requirement is satisfied by them."""
    capabilities = reference.capabilities()
    available = {level: capabilities.available(evidence=level) for level in ("assertion", "source")}
    return {"runtime": RuntimeHashes(reference.catalogue_hash, reference.models_hash, reference.execution_hash,
                                     reference.artifact_hash),
            "capabilities": RuntimeCapabilities.from_dict(capabilities.to_dict()),
            "requirements": tuple(RequirementOutcome(r.capability, r.evidence, r.capability in available[r.evidence])
                                  for r in recipe.requires)}


def _report(build: _Build, inputs: dict, records: tuple, claims: tuple, reference: Reference,
            generator: Generator) -> BuildReport:
    """The complete build report; its construction re-checks every document, hash and count."""
    lock, recipe, inventory = build.lock, build.recipe, build.inventory
    comparison = ogrdb.comparison(lock, inventory)
    segments = {s: 0 for s in ("V", "D", "J", "unassigned")}
    for entry in recipe.selection:
        segments[entry.segment] += 1
    counts = Counts(
        exports=tuple(ExportCounts(asset, **{k: SegmentCounts(**v) for k, v in inventory.export_counts(asset).items()})
                      for asset in comparison.exports),
        selected_alleles=SegmentCounts(**segments),
        record_outcomes=OutcomeCounts(**{k: sum(r.outcome == k for r in records)
                                         for k in ("selected", "excluded", "rejected")}))
    return BuildReport(
        generator=generator, source_lock=lock, recipe=recipe, inventory=inventory,
        hashes=BuildHashes(source_lock=inputs["source_lock"], inventory=inputs["inventory"], recipe=inputs["recipe"],
                           build_inputs=build.build_input_hash),
        comparison=comparison, counts=counts,
        selection=tuple(SelectedAllele(e.allele_id, c.locator) for e, c in zip(build.entries, build.chosen)),
        records=records, claims=claims, **_runtime_accounting(recipe, reference))


def report_for(report: BuildReport, reference: Reference, generator: Generator) -> BuildReport:
    """The build report of an already projected build, stated for ``reference`` and ``generator``.

    ``reference`` is a runtime artifact with the same runtime content as the projected one, for
    example the artifact a dossier stores, written by the build that ``generator`` names.  The
    source, recipe, inventory, selection and every outcome are ``report``'s; the runtime hashes,
    capabilities and requirement outcomes are ``reference``'s.  A dossier check compares the
    result with the stored report, so an earlier build is never judged by the artifact identity
    the checking version would give it.
    """
    require_type(report, BuildReport, "report")
    require_type(reference, Reference, "reference")
    require_type(generator, Generator, "generator")
    try:
        return dataclasses.replace(report, generator=generator, **_runtime_accounting(report.recipe, reference))
    except ReferenceSourceError as exc:
        raise exc.under("/report") from None
