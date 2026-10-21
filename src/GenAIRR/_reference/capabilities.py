"""Capability evaluation for a reference catalogue, with explicit evidence levels."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

from GenAIRR._airr_record import _translate_codon
from GenAIRR._reference.model import Catalogue

CAPABILITIES = (
    "catalogue_inspection",
    "observed_span_simulation",
    "source_annotated_junction",
    "canonical_productive",
    "original_coding_ends",
    "complete_reference_span",
    "genotype_constrained",
    "duplex_family_replay",
)

_EVIDENCE_RANK = {"assertion": 1, "source": 2}
_UNSUPPORTED_NOTE = "not yet supported for portable references"
_ANCHOR_HINT = 'add an "anchor" annotation (support "source" with a citation, or "assertion")'
_BOUNDARY_HINT = 'add "{p}": {{"value": "complete", ...}} with its support'
_J_RESIDUE = {"IGH": {"W"}, "IGK": {"F"}, "IGL": {"F"}, "TRA": {"F"}, "TRB": {"F"},
              "TRG": {"F"}, "TRD": {"F"}}


@dataclass(frozen=True)
class Unmet:
    allele: str
    requirement: str
    hint: str

    def to_dict(self):
        return {"allele": self.allele, "requirement": self.requirement, "hint": self.hint}


@dataclass(frozen=True)
class CapabilityStatus:
    name: str
    available: bool
    evidence: Optional[str]
    unmet: tuple
    note: Optional[str]

    def to_dict(self):
        return {"name": self.name, "available": self.available, "evidence": self.evidence,
                "unmet": [u.to_dict() for u in self.unmet], "note": self.note}


@dataclass(frozen=True)
class CapabilityReport:
    statuses: tuple
    phase_support: str

    def status(self, name: str) -> CapabilityStatus:
        for s in self.statuses:
            if s.name == name:
                return s
        raise KeyError(name)

    def available(self, *, evidence: str = "assertion") -> tuple:
        if evidence not in _EVIDENCE_RANK:
            raise ValueError(f"evidence must be one of {sorted(_EVIDENCE_RANK)}")
        need = _EVIDENCE_RANK[evidence]
        return tuple(s.name for s in self.statuses
                     if s.available and (s.evidence is None or _EVIDENCE_RANK[s.evidence] >= need))

    def explain(self) -> str:
        lines = []
        for s in self.statuses:
            if s.available:
                tag = f" [evidence: {s.evidence}]" if s.evidence else ""
                lines.append(f"available   {s.name}{tag}")
            elif s.unmet:
                lines.append(f"unavailable {s.name}: {len(s.unmet)} requirement(s) unmet")
            else:
                lines.append(f"unavailable {s.name}: {s.note}")
            for u in s.unmet:
                lines.append(f"- {u.allele}: {u.requirement} — {u.hint}")
        lines.append(f"phase support: {self.phase_support}")
        return "\n".join(lines)

    def to_dict(self):
        return {"statuses": [s.to_dict() for s in self.statuses], "phase_support": self.phase_support}


def _accepted(catalogue, allele_id, predicate):
    """Accepted claims, or () if any conflict claim exists for the predicate."""
    claims = catalogue.claims_for(allele_id, predicate)
    if any(c.status == "conflict" for c in claims):
        return ()
    return tuple(c for c in claims if c.status == "accepted")


def _usable(catalogue, allele_id, predicate, want=None):
    """Accepted claims when they agree on one value (and match ``want``), else ()."""
    claims = _accepted(catalogue, allele_id, predicate)
    values = {repr(c.value) for c in claims}
    if len(values) != 1:
        return ()
    if want is not None and claims[0].value != want:
        return ()
    return claims


def _weakest(claims):
    if not claims:
        return None
    return "source" if all(c.support == "source" for c in claims) else "assertion"


def _phase_support(catalogue, ids):
    if not ids:
        return "none"
    claims = [_accepted(catalogue, i, "coding_phase") for i in ids]
    if not all(claims):
        return "none"
    flat = [c for cs in claims for c in cs]
    if any(c.derivation == "anchor_derived" for c in flat):
        return "anchor_derived"
    if all(c.derivation == "stated" and c.support in ("source", "assertion") for c in flat):
        return "independent_source" if all(c.support == "source" for c in flat) else "independent_assertion"
    return "none"


def _anchor_hint(catalogue, allele_id):
    for c in catalogue.claims_for(allele_id, "anchor"):
        if c.status == "candidate" and c.derivation == "unique_motif":
            return (f"{_ANCHOR_HINT}; candidate at {c.value} (unique J motif) — "
                    "copy it only if your evidence confirms it")
    return _ANCHOR_HINT


def _codon(catalogue, allele_id, position):
    return _translate_codon(catalogue.allele(allele_id).sequence[position:position + 3])


def evaluate(catalogue: Catalogue, view: dict, *, native_fatal: Sequence[dict] = ()) -> CapabilityReport:
    vs, ds, js = (list(view.get(k, ())) for k in ("V", "D", "J"))
    vdj = catalogue.header.chain_type == "vdj"
    statuses = []

    def add(name, unmet=(), evidence=None, note=None, force_unavailable=False):
        ok = not unmet and not force_unavailable
        statuses.append(CapabilityStatus(name, ok, evidence if ok else None, tuple(unmet), note))

    add("catalogue_inspection")

    unmet = []
    if not vs:
        unmet.append(Unmet("(view)", "at least one V allele", "select at least one V allele in the view"))
    if not js:
        unmet.append(Unmet("(view)", "at least one J allele", "select at least one J allele in the view"))
    if vdj and not ds:
        unmet.append(Unmet("(view)", "at least one D allele", "select at least one D allele in the view"))
    for f in native_fatal:
        kind = f.get("kind", "finding") if isinstance(f, dict) else str(f)
        detail = f.get("detail", "") if isinstance(f, dict) else ""
        unmet.append(Unmet("(view)", f"no fatal engine finding ({kind})",
                           detail or "resolve the reported engine finding"))
    add("observed_span_simulation", unmet)

    view_unmet = [u for u in unmet if u.requirement.startswith("at least one") and
                  not u.requirement.endswith("D allele")]

    # Junction anchors.
    junction_unmet, anchor_claims, anchors = list(view_unmet), [], {}
    for aid in vs + js:
        claims = _usable(catalogue, aid, "anchor")
        if claims:
            anchor_claims.extend(claims)
            anchors[aid] = claims[0].value
        else:
            junction_unmet.append(Unmet(aid, "accepted anchor", _anchor_hint(catalogue, aid)))
    add("source_annotated_junction", junction_unmet, _weakest(anchor_claims))

    # Canonical productive.
    expected = _J_RESIDUE.get(catalogue.header.locus, {"F", "W"})
    canon = list(junction_unmet)
    for aid in vs:
        if aid in anchors and _codon(catalogue, aid, anchors[aid]) != "C":
            canon.append(Unmet(aid, "V anchor codon encodes C",
                               "check the anchor position or leave this allele out of the view"))
    for aid in js:
        if aid in anchors and _codon(catalogue, aid, anchors[aid]) not in expected:
            canon.append(Unmet(aid, f"J anchor codon encodes {' or '.join(sorted(expected))}",
                               "check the anchor position or leave this allele out of the view"))
    phase = _phase_support(catalogue, vs + js)
    add("canonical_productive", canon, _weakest(anchor_claims), note=f"phase support: {phase}")

    # Boundaries.
    def boundaries(required):
        unmet, claims = list(view_unmet), []
        for aid, preds in required:
            for p in preds:
                found = _usable(catalogue, aid, p, "complete")
                if found:
                    claims.extend(found)
                else:
                    unmet.append(Unmet(aid, f'accepted {p} == "complete"', _BOUNDARY_HINT.format(p=p)))
        return unmet, _weakest(claims)

    req = [(a, ("boundary_3",)) for a in vs] + [(a, ("boundary_5",)) for a in js] + \
          [(a, ("boundary_5", "boundary_3")) for a in ds]
    unmet, ev = boundaries(req)
    add("original_coding_ends", unmet, ev)
    req = [(a, ("boundary_5", "boundary_3")) for a in vs + ds + js]
    unmet, ev = boundaries(req)
    add("complete_reference_span", unmet, ev)

    for name in ("genotype_constrained", "duplex_family_replay"):
        add(name, note=_UNSUPPORTED_NOTE, force_unavailable=True)
    if tuple(s.name for s in statuses) != CAPABILITIES:
        raise RuntimeError("capability statuses out of order")
    return CapabilityReport(tuple(statuses), phase)
