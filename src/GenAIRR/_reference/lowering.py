"""Lower a validated catalogue view into the engine's reference configuration.

Lowering is lossy in one direction only: claims that are not accepted, or
that disagree with each other, stay unknown natively.  Curatable native
findings that the observed-span profile can tolerate are admitted and
recorded; anything else is refused.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from GenAIRR import _engine
from GenAIRR._refdata_resolver import _j_anchor_residues_for_locus
from GenAIRR._reference.capabilities import _usable
from GenAIRR._reference.errors import ReferenceArtifactError
from GenAIRR._reference.model import SUBREGION_LABELS, Catalogue

PROFILE_ID = "observed_span_simulation"
ADMITTED_CURATABLE_KINDS = frozenset(
    {"MissingAnchor", "VAnchorNotCys", "JAnchorUnexpectedAa", "ReferenceEvidenceProblem"})
_TAG = "genairr_reference_v1"
_SEGMENT_NAMES = {"V": "v", "D": "d", "J": "j"}


@dataclass(frozen=True)
class Lowered:
    refdata: "_engine.RefDataConfig"
    spec: dict
    admission: dict


def default_view(catalogue: Catalogue) -> Dict[str, List[str]]:
    """Executable allele ids per segment, in accepted-record order."""
    segments = ("V", "J") if catalogue.header.chain_type == "vj" else ("V", "D", "J")
    view: Dict[str, List[str]] = {s: [] for s in segments}
    by_id = {a.id: a for a in catalogue.alleles}
    for record in catalogue.records:
        if record.status != "accepted" or record.allele is None:
            continue
        allele = by_id[record.allele]
        if allele.segment in view and record.allele not in view[allele.segment]:
            view[allele.segment].append(record.allele)
    return view


def _single(catalogue, allele_id, predicate):
    claims = _usable(catalogue, allele_id, predicate)
    return claims[0] if claims else None


def _evidence(catalogue, allele, phase) -> dict:
    anchor_claim = _single(catalogue, allele.id, "anchor") if allele.segment != "D" else None
    anchor_claims = _usable(catalogue, allele.id, "anchor") if allele.segment != "D" else ()
    boundaries = {}
    notes = []
    for predicate in ("boundary_5", "boundary_3"):
        c = _single(catalogue, allele.id, predicate)
        boundaries[predicate] = c.value if c else "unknown"
        if c:
            notes.append(f"{predicate}={c.support}:{c.derivation}")
    phase_claim = _single(catalogue, allele.id, "coding_phase")
    if phase_claim:
        notes.append(f"coding_phase={phase_claim.support}:{phase_claim.derivation}")
    func = _single(catalogue, allele.id, "functionality")
    if func:
        notes.append(f"functionality={func.support}:{func.derivation}")
    if not anchor_claims:
        status, method, source = "unknown", f"{_TAG}:none", None
    else:
        # Evidence is only as strong as the weakest accepted anchor claim.
        weakest = next((c for c in anchor_claims if c.support != "source"), anchor_claims[0])
        status = "verified" if weakest.support == "source" else "inferred"
        # When several claims share the weakest support, the first in claim order supplies the method.
        method = f"{_TAG}:{weakest.support}:{weakest.derivation}"
        # Source: first non-null source among source-supported claims, else among all.
        source = next((c.source for c in anchor_claims if c.support == "source" and c.source), None)
        if source is None:
            source = next((c.source for c in anchor_claims if c.source), None)
    return {
        "schema_version": 1,
        "sequence_sha256": allele.sequence_sha256,
        "coding_phase": phase,
        "anchor_status": status,
        "boundary_5": boundaries["boundary_5"],
        "boundary_3": boundaries["boundary_3"],
        "source": source,
        "method": method,
        "notes": sorted(notes),
    }


def _check_view(catalogue: Catalogue, view: Dict[str, List[str]]) -> None:
    valid = ("V", "J") if catalogue.header.chain_type == "vj" else ("V", "D", "J")
    executable = {a for ids in default_view(catalogue).values() for a in ids}
    for segment, ids in view.items():
        if segment not in valid:
            raise ReferenceArtifactError(
                "integrity", f"view.{segment}", f"segment not valid for chain type {catalogue.header.chain_type}")
        for allele_id in ids:
            if allele_id not in executable or catalogue.allele(allele_id).segment != segment:
                raise ReferenceArtifactError(
                    "integrity", f"view.{segment}", f"{allele_id!r} is not an executable {segment} allele")


def lower(catalogue: Catalogue, view: Dict[str, List[str]]) -> Lowered:
    header = catalogue.header
    _check_view(catalogue, view)
    ref = _engine.RefDataConfig(header.chain_type)
    ref.set_allowed_bases(["A", "C", "G", "T"])
    ref.set_v_anchor_rule(expected_aa=["C"])
    ref.set_j_anchor_rule(
        expected_aa=_j_anchor_residues_for_locus(header.locus) or ["F", "W"])
    ref.set_identity(
        species=header.taxon.scientific_name, locus=header.locus,
        reference_set=header.name if header.version is None else f"{header.name}@{header.version}",
        name=header.name, source=_TAG)

    adders = {"V": ref.add_v_allele, "D": getattr(ref, "add_d_allele", None), "J": ref.add_j_allele}
    segments = ("V", "J") if header.chain_type == "vj" else ("V", "D", "J")
    for segment in segments:
        for index, allele_id in enumerate(view.get(segment, [])):
            allele = catalogue.allele(allele_id)
            anchor_claim = _single(catalogue, allele_id, "anchor") if segment != "D" else None
            anchor = anchor_claim.value if anchor_claim else None
            func = _single(catalogue, allele_id, "functionality")
            phase_claim = _single(catalogue, allele_id, "coding_phase")
            phase = phase_claim.value if phase_claim and phase_claim.derivation == "stated" else None
            if phase is not None and anchor is not None and phase != anchor % 3:
                raise ReferenceArtifactError(
                    "conflict", allele_id,
                    f"stated coding phase {phase} disagrees with anchor {anchor}")
            kwargs: Dict[str, Any] = {"functional_status": func.value if func else None}
            if segment != "D":
                kwargs["anchor"] = anchor
            if segment == "V":
                sub = _single(catalogue, allele_id, "v_subregions")
                # IMGT order, never the claim value's key order (which canonical JSON sorts).
                kwargs["subregions"] = (
                    [(label, *sub.value[label]) for label in SUBREGION_LABELS if label in sub.value]
                    if sub else None)
            adders[segment](allele.label, allele.gene, allele.sequence.encode("ascii"), **kwargs)
            ref.set_allele_evidence(
                segment, index, json.dumps(_evidence(catalogue, allele, phase)))

    findings = ref.validate()
    fatal = [f for f in findings if f.get("severity") == "fatal"]
    if fatal:
        raise ReferenceArtifactError(
            "integrity", "", "native validation failed: " + "; ".join(json.dumps(f, sort_keys=True) for f in fatal))
    admitted = []
    for f in findings:
        if f["kind"] not in ADMITTED_CURATABLE_KINDS:
            raise ReferenceArtifactError(
                "capability", "", f"native finding {f['kind']} is not admitted by profile {PROFILE_ID}")
        segment = f.get("segment")
        index = f.get("allele_id")
        allele = view[segment][index] if segment in view and index is not None else None
        admitted.append({"kind": f["kind"], "segment": segment, "allele": allele})
    admitted.sort(key=lambda a: (a["kind"], a["segment"] or "", a["allele"] or ""))
    return Lowered(ref, json.loads(ref.reference_spec_json()),
                   {"profile": PROFILE_ID, "accepted_curatable_findings": admitted})


def seal(lowered: Lowered, context: dict) -> "_engine.RefDataConfig":
    lowered.refdata.set_reference_context(json.dumps(context))
    return lowered.refdata
