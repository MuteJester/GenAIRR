"""Builders for the reference-source model objects.

Typed models with a lot of required structure; these assemble valid instances
so a test can state only the field it is interested in.
"""
from __future__ import annotations
import copy
import json
from pathlib import Path
import pytest
from GenAIRR._artifact_io import canonical_json, strict_loads, typed_digest
from GenAIRR._reference_sources.errors import KINDS, ReferenceSourceError
from GenAIRR._reference_sources import model as m
from tests.support import paths

ROOT = paths.data("reference_sources")
IGH = "9606.IGH_VDJ"
ROLE = {"airr": "record_export", "airr_ex": "record_export", "gapped": "gapped_alignment",
        "gapped_ex": "gapped_alignment"}
ADAPTER = {"id": "ogrdb_germline_set", "version": 1,
           "coordinate_profile": {"id": "ogrdb_api_v2_human_ig_observed_v1", "version": 1}}
def _load_json(path: Path, max_bytes: int = 1 << 23):
    return strict_loads(path.read_bytes().decode("utf-8"), max_bytes=max_bytes)
SOURCES = _load_json(ROOT / "sources.json")
def _source_set(set_id):
    (found,) = [s for s in SOURCES["sets"] if s["set_id"] == set_id]
    return found
def _raises(kind, path=None):
    """Context manager: a ReferenceSourceError of ``kind`` (and ``path`` if given)."""
    class _Check:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            assert exc_type is not None, f"expected ReferenceSourceError({kind!r})"
            assert issubclass(exc_type, ReferenceSourceError), f"unexpected {exc_type.__name__}: {exc}"
            assert exc.kind == kind, f"expected kind {kind!r}, got {exc.kind!r}: {exc}"
            if path is not None:
                assert exc.path == path, f"expected path {path!r}, got {exc.path!r}: {exc}"
            self.error = exc
            return True
    return _Check()
def _reseal(lock: dict) -> dict:
    """Recompute a lock dict's commitment after an edit (a rehashed document)."""
    body = {k: v for k, v in lock.items() if k != "commitment"}
    lock["commitment"] = typed_digest(body, "genairr_reference_source_lock_v1")
    return lock
def lock_dict(set_id=IGH, formats=("airr", "airr_ex", "gapped", "gapped_ex"),
              origin="provider_download") -> dict:
    source = _source_set(set_id)
    assets = []
    for asset in source["assets"]:
        if asset["format"] in formats:
            fmt = asset["format"]
            assets.append({"asset_id": fmt, "role": ROLE[fmt], "format": fmt,
                           "sha256": "sha256:" + asset["sha256"], "bytes": asset["bytes"],
                           "declared_url": asset["url"] if origin == "provider_download" else None,
                           "origin_kind": origin})
    missing = [{"format": f, "role": ROLE[f], "reason": "not_requested", "detail": None}
               for f in ("airr", "airr_ex", "gapped", "gapped_ex") if f not in formats]
    locus = source["locus"]
    lock = {
        "artifact": "genairr_reference_source_lock", "schema_version": 1, "provider": "ogrdb",
        "set_id": source["set_id"], "release": source["release"],
        "identity": {"taxon": {"scientific_name": source["species"]["label"],
                               "taxonomy_id": int(source["species"]["id"].split(":")[1])},
                     "locus": locus, "chain_type": "vdj" if locus == "IGH" else "vj"},
        "adapter": copy.deepcopy(ADAPTER), "assets": assets, "missing_assets": missing,
        "terms": {"status": "unresolved", "assets": [], "licence": None, "citation": None,
                  "reason": "data terms statement not captured with these exports"},
    }
    return _reseal(lock)
def receipt_dict(fmt="airr", origin="provider_download") -> dict:
    (asset,) = [a for a in _source_set(IGH)["assets"] if a["format"] == fmt]
    provider = origin == "provider_download"
    return {
        "artifact": "genairr_reference_source_receipt", "schema_version": 1,
        "event_time_utc": asset["retrieved_at_utc"], "origin": origin, "format": fmt,
        "sha256": "sha256:" + asset["sha256"], "bytes": asset["bytes"],
        "requested_url": asset["url"] if provider else None,
        "final_url": asset["final_url"] if provider else None,
        "redirects": [],
        "http_status": asset["http_status"] if provider else None,
        "headers": ({"content-type": asset["content_type"], "content-disposition": asset["content_disposition"]}
                    if provider else {}),
    }
_RECORDS = {}
def _records(fmt, set_id=IGH):
    key = (set_id, fmt)
    if key not in _RECORDS:
        source = _source_set(set_id)
        document = _load_json(ROOT / "ogrdb" / set_id / str(source["release"]) / f"{fmt}.json")
        (germline_set,) = document["GermlineSet"]
        _RECORDS[key] = germline_set["allele_descriptions"]
    return _RECORDS[key]
def _pointer(index):
    return f"/GermlineSet/0/allele_descriptions/{index}"
def _find(label, fmt="airr", set_id=IGH):
    (found,) = [(i, r) for i, r in enumerate(_records(fmt, set_id)) if r["label"] == label]
    return found
def selection_entry(label, fmt="airr", reason="Use this supplied observed sequence for the declared workflow"):
    index, record = _find(label, fmt)
    return {"segment": record["sequence_type"], "label": label,
            "record": {"asset": fmt, "pointer": _pointer(index)},
            "sequence_sha256": m.sequence_digest(record["coding_sequence"]), "reason": reason}
SELECTED = ("IGHV1-58*03", "IGHV1-69*01", "IGHV1-69D*01", "IGHD1-1*01", "IGHJ4*02")
def recipe_dict(lock=None, selected=SELECTED) -> dict:
    lock = lock or lock_dict()
    return {
        "artifact": "genairr_reference_recipe", "schema_version": 1,
        "source_lock_hash": lock["commitment"],
        "adapter": copy.deepcopy(ADAPTER), "projection": {"id": "observed_span", "version": 1},
        "reference": {"name": "human_igh_ogrdb_r10", "version": "1", "taxon": copy.deepcopy(lock["identity"]["taxon"]),
                      "locus": lock["identity"]["locus"], "chain_type": lock["identity"]["chain_type"],
                      "build_date": "2026-10-22"},
        "selection": [selection_entry(label) for label in selected],
        "decisions": [],
        "requires": [{"capability": "observed_span_simulation", "evidence": "assertion"}],
        "defaults_profile": "genairr.unfitted_defaults.v1",
    }
CLAIM_ID = "claim:" + "a" * 64
FINDING_ID = "finding:" + "b" * 64
def decision(action="select_alternative", **overrides):
    entry = selection_entry("IGHV1-58*03")
    d = {"action": action, "findings": [FINDING_ID] if action == "select_alternative" else [],
         "claims": [] if action == "select_alternative" else [CLAIM_ID],
         "sequence_sha256": entry["sequence_sha256"],
         "reason": "The airr record ends at gene_end like every other selected V; the extension is not adopted",
         "evidence": [{"asset": "airr", "pointer": entry["record"]["pointer"] + "/coding_sequence"}]}
    d.update(overrides)
    return d
L = m.Locator
def _fasta(fmt, set_id=IGH):
    source = _source_set(set_id)
    text = (ROOT / "ogrdb" / set_id / str(source["release"]) / f"{fmt}.fasta").read_bytes().decode("ascii")
    records = []
    for line in text.split("\n"):
        if line.startswith(">"):
            records.append([line[1:], ""])
        elif line:
            records[-1][1] += line
    return records
def _gapped_index(label, fmt):
    (found,) = [(i, seq) for i, (name, seq) in enumerate(_fasta(fmt)) if name == label]
    return found
def _blocks(gapped):
    """Exact interval map from ungapped positions to gapped columns."""
    blocks, position = [], 0
    for column, char in enumerate(gapped):
        if char == ".":
            continue
        if blocks and blocks[-1][1][1] == column:
            blocks[-1] = ((blocks[-1][0][0], position + 1), (blocks[-1][1][0], column + 1))
        else:
            blocks.append(((position, position + 1), (column, column + 1)))
        position += 1
    return tuple(m.Block(input=i, output=o) for i, o in blocks)
def _sorted(locators):
    return tuple(sorted(locators, key=L.sort_key))
INV_RECORDS = (("airr", "IGHV1-58*03"), ("airr", "IGHV1-69*01"), ("airr", "IGHV1-69D*01"), ("airr", "IGHD1-1*01"),
               ("airr", "IGHJ4*02"), ("airr_ex", "IGHV1-58*03"), ("airr_ex", "IGHV1-69*01"),
               ("airr_ex", "IGHD1-1*01"), ("airr_ex", "IGHJ4*02"))
class Inv:
    """A small realistic inventory over pinned IGH release 10 records, built field by field."""

    def __init__(self, lock=None):
        self.lock = lock or m.SourceLock.from_dict(lock_dict())
        self.records, self.candidates, self.observations, self.support = [], [], [], []
        self.mentioned, self.mappings, self.claims, self.findings = [], [], [], []
        self.rec, self.cand, self.obs = {}, {}, {}
        for fmt, label in INV_RECORDS:
            index, raw = _find(label, fmt)
            loc = L(fmt, _pointer(index))
            self.rec[fmt, label] = loc
            self.records.append(m.SourceRecord(
                locator=loc, status="parsed", problem=None, label=raw["label"],
                provider_id=raw["allele_description_id"], provider_ref=raw["allele_description_ref"],
                segment=raw["sequence_type"], taxonomy_id=int(raw["species"]["id"].split(":")[1]),
                locus=raw["locus"]))
            fields = ["label", "coding_sequence", "functional", "sequence_type"]
            if raw["sequence_type"] == "V":
                fields.append("v_gene_delineations/0/cdr3_start")
            if raw["sequence_type"] == "J":
                fields += ["j_codon_frame", "j_cdr3_end"]
            if raw["paralogs"]:
                fields.append("paralogs/0")
            for field in fields:
                value = raw
                for token in field.split("/"):
                    value = value[int(token)] if isinstance(value, list) else value[token]
                self._observe(L(fmt, f"{loc.pointer}/{field}"), loc, value)
            cand = L(fmt, loc.pointer + "/coding_sequence")
            self.cand[fmt, label] = cand
            self.candidates.append(m.SequenceCandidate(locator=cand, record=loc, sequence=raw["coding_sequence"],
                                                       sha256=m.sequence_digest(raw["coding_sequence"])))
            self.claims.append(self.claim(loc, "functionality", "functional", [self.obs_at(loc, "functional")]))
        # Paralog-only identity in the extended export, mentioned by IGHV1-69*01.
        ex69 = self.rec["airr_ex", "IGHV1-69*01"]
        self.mentioned.append(m.MentionedIdentity(asset="airr_ex", segment="V", label="IGHV1-69D*01", mentions=(
            m.Mention(record=ex69, locator=L("airr_ex", ex69.pointer + "/paralogs/0"), relationship="paralog"),)))
        # Embedded, unretrieved repository support of one record.
        rec5803 = self.rec["airr", "IGHV1-58*03"]
        _, raw5803 = _find("IGHV1-58*03")
        entry = raw5803["unrearranged_support"][0]
        self.support.append(m.SupportReference(
            locator=L("airr", rec5803.pointer + "/unrearranged_support/0"), record=rec5803, kind="unrearranged",
            repository=entry["repository_name"], reference=entry["repository_ref"], status="embedded_unretrieved"))
        # Set-level observation: the API version is not the schema version or the data release.
        info = L("airr", "/Info/version")
        self._observe(info, None, _load_json(ROOT / "ogrdb" / IGH / "10" / "airr.json")["Info"]["version"],
                      interpretation="preserved")
        # Variant-specific gapped alignments, each mapping one record occurrence to its own variant.
        for fmt, gapped_fmt in (("airr", "gapped"), ("airr_ex", "gapped_ex")):
            gindex, gapped = _gapped_index("IGHV1-58*03", gapped_fmt)
            gloc = L(gapped_fmt, f"/{gindex}")
            self._observe(gloc, None, {"header": "IGHV1-58*03", "sequence": gapped})
            cand = self.cand[fmt, "IGHV1-58*03"]
            sequence = gapped.replace(".", "")
            mapping = m.Mapping(source_lock_hash=self.lock.commitment, kind="gapped_alignment", input=cand,
                                output=gloc,
                                input_sha256=m.sequence_digest(sequence),
                                output_sha256=m.text_digest(gapped), input_space="coding_sequence",
                                output_space="imgt_gapped", blocks=_blocks(gapped), transform="remove_imgt_gaps",
                                validation="exact")
            self.mappings.append(mapping)
        # V junction anchor: source-stated, and re-derived from the variant's own gapped record.
        cand, rec = self.cand["airr", "IGHV1-58*03"], self.rec["airr", "IGHV1-58*03"]
        cdr3 = self.obs_at(rec, "v_gene_delineations/0/cdr3_start")
        anchor = _find("IGHV1-58*03")[1]["v_gene_delineations"][0]["cdr3_start"] - 4
        gapped = _gapped_index("IGHV1-58*03", "gapped")[1]
        assert len(gapped[:312].replace(".", "")) - 3 == anchor
        self.anchor = anchor
        self.claims.append(self.claim(cand, "anchor", anchor, [cdr3], system="coding_sequence"))
        self.claims.append(self.claim(cand, "anchor", anchor, [L("gapped", f"/{_gapped_index('IGHV1-58*03', 'gapped')[0]}")],
                                      system="coding_sequence", derivation="imgt_alignment",
                                      mappings=[self.mappings[0].id]))
        self.claims.append(self.claim(cand, "coding_phase", anchor % 3, [cdr3], system="coding_sequence",
                                      derivation="anchor_derived"))
        self.claims.append(self.claim(cand, "phase_independence", "anchor_derived", [cdr3],
                                      derivation="anchor_derived"))
        # J: the reported phase value is accepted; whether it is independent of the landmark is unresolved.
        jcand, jrec = self.cand["airr", "IGHJ4*02"], self.rec["airr", "IGHJ4*02"]
        _, jraw = _find("IGHJ4*02")
        frame = self.obs_at(jrec, "j_codon_frame")
        self.claims.append(self.claim(jcand, "anchor", jraw["j_cdr3_end"] - 1, [self.obs_at(jrec, "j_cdr3_end")],
                                      system="coding_sequence"))
        self.claims.append(self.claim(jcand, "coding_phase", jraw["j_codon_frame"] - 1, [frame],
                                      system="coding_sequence"))
        self.claims.append(self.claim(jcand, "phase_independence", None, [frame], status="unresolved",
                                      support="none"))
        # Extended-export alternative: provider 3' extension beyond the observed span.
        ex = self.cand["airr_ex", "IGHV1-58*03"]
        self.claims.append(self.claim(ex, "extension", {"side": "3", "length": 1, "matches_record_context": True},
                                      [L("airr_ex", ex.pointer), L("airr", cand.pointer)],
                                      derivation="variant_comparison"))
        # D: no universal coding phase; its physical 3' end stays unresolved.
        dcand = self.cand["airr", "IGHD1-1*01"]
        self.claims.append(self.claim(dcand, "boundary_3", None, [L("airr", dcand.pointer)], status="unresolved",
                                      support="none"))
        self.findings.append(m.Finding(
            source_lock_hash=self.lock.commitment, code="variant_coding_sequence_differs", severity="warning", category="contradictory_evidence",
            subjects=_sorted([cand, ex]), locators=_sorted([cand, ex]),
            observations={"label": "IGHV1-58*03", "extension_3": "A"},
            capabilities=("observed_span_simulation", "canonical_productive")))
        self.findings.append(m.Finding(
            source_lock_hash=self.lock.commitment, code="paralog_only_identity", severity="info", category="informational",
            subjects=(self.rec["airr_ex", "IGHV1-69*01"],), locators=(L("airr_ex", ex69.pointer + "/paralogs/0"),),
            observations={"label": "IGHV1-69D*01"}, capabilities=()))
        self.findings.append(m.Finding(
            source_lock_hash=self.lock.commitment, code="api_version_is_not_release", severity="info", category="unsupported_interpretation",
            subjects=(), locators=(info,), observations={"value": "2.0.0"}, capabilities=()))

    def _observe(self, locator, record, value, interpretation="interpreted"):
        self.obs[locator.id] = locator
        self.observations.append(m.Observation(locator=locator, record=record, presence="present", value=value,
                                               interpretation=interpretation))

    def obs_at(self, record, field):
        return self.obs[L(record.asset, f"{record.pointer}/{field}").id]

    def claim(self, subject, predicate, value, evidence, *, system=None, derivation="stated", status="accepted",
              support="source", mappings=()):
        return m.Claim(source_lock_hash=self.lock.commitment, subject=subject, predicate=predicate, value=value, coordinate_system=system, support=support,
                       derivation=derivation, status=status, evidence=_sorted(evidence), mappings=tuple(mappings))

    def build(self, **overrides):
        fields = dict(source_lock_hash=self.lock.commitment, adapter=self.lock.adapter, records=self.records,
                      candidates=self.candidates, mentioned_identities=self.mentioned,
                      observations=self.observations, support=self.support, claims=self.claims,
                      mappings=self.mappings, findings=self.findings)
        fields.update(overrides)
        return m.SourceInventory.assemble(**fields)
@pytest.fixture(scope="module")
def inv():
    return Inv()
@pytest.fixture(scope="module")
def inventory(inv):
    return inv.build()
def _inventory_dict(inventory):
    return copy.deepcopy(inventory.to_dict())
def _nested(depth):
    value = 0
    for _ in range(depth):
        value = [value]
    return value
def _recipe(inventory, *, selected=SELECTED, decisions=(), edit=None):
    d = recipe_dict(selected=selected)
    d["decisions"] = list(decisions)
    if edit:
        edit(d)
    return m.Recipe.from_dict(d)
def _finding(inventory, code):
    (found,) = [f for f in inventory.findings if f.code == code]
    return found
def _claim(inventory, predicate, subject):
    (found,) = [c for c in inventory.claims if c.predicate == predicate and c.subject == subject]
    return found
RUNTIME = {"catalogue_hash": "sha256:" + "1" * 64, "models_hash": "sha256:" + "2" * 64,
           "execution_hash": "sha256:" + "3" * 64, "artifact_hash": "sha256:" + "4" * 64}
GENERATOR = {"name": "GenAIRR", "version": "2.3.1"}
def _capabilities():
    def status(name, available, evidence=None, unmet=(), note=None):
        return {"name": name, "available": available, "evidence": evidence,
                "unmet": [{"allele": a, "requirement": r, "hint": h} for a, r, h in unmet], "note": note}
    blocker = ("V:IGHV1-69*01", 'accepted boundary_3 == "complete"', 'add "boundary_3": {"value": "complete", ...}')
    return {"statuses": [
        status("catalogue_inspection", True), status("observed_span_simulation", True),
        status("source_annotated_junction", True, "source"),
        status("canonical_productive", True, "assertion", note="phase support: anchor_derived"),
        status("original_coding_ends", False, unmet=[blocker]),
        status("complete_reference_span", False, unmet=[blocker]),
        status("genotype_constrained", False, note="not yet supported for portable references"),
        status("duplex_family_replay", False, note="not yet supported for portable references")],
        "phase_support": "anchor_derived"}
def report_decisions(inventory, inv):
    jcand, dcand = inv.cand["airr", "IGHJ4*02"], inv.cand["airr", "IGHD1-1*01"]
    return [decision(findings=[_finding(inventory, "variant_coding_sequence_differs").id]),
            decision("withhold_claim", claims=[_claim(inventory, "boundary_3", dcand).id],
                     sequence_sha256=inventory.candidate(dcand).sha256,
                     reason="The physical 3' end of this D is not established by the source"),
            decision("accept_claim", claims=[_claim(inventory, "coding_phase", jcand).id],
                     sequence_sha256=inventory.candidate(jcand).sha256,
                     evidence=[inv.obs_at(inv.rec["airr", "IGHJ4*02"], "j_codon_frame").to_dict()],
                     reason="j_codon_frame - 1 is the reported phase under the pinned producer profile")]
def report_dict(inventory, inv, *, decisions=None) -> dict:
    lock = m.SourceLock.from_dict(lock_dict())
    recipe = _recipe(inventory, decisions=report_decisions(inventory, inv) if decisions is None else decisions)
    selected = {e.record.id: e for e in recipe.selection}
    chosen = {}
    for e in recipe.selection:
        (c,) = [c for c in inventory.candidates_of(e.record) if c.sha256 == e.sequence_sha256]
        chosen[e.record.id] = c
    by_label = {(e.segment, e.label): e for e in recipe.selection}
    records = []
    for r in inventory.records:
        if r.id in selected:
            records.append({"record": r.locator.to_dict(), "outcome": "selected", "allele": selected[r.id].allele_id,
                            "exclusion": None, "related": None, "detail": None})
            continue
        if r.status != "parsed":
            records.append({"record": r.locator.to_dict(), "outcome": "rejected", "allele": None,
                            "exclusion": None, "related": None, "detail": r.problem})
            continue
        entry = by_label.get((r.segment, r.label))
        if entry is None:
            records.append({"record": r.locator.to_dict(), "outcome": "excluded", "allele": None,
                            "exclusion": "outside_explicit_selection", "related": None, "detail": None})
            continue
        same = any(c.sha256 == entry.sequence_sha256 for c in inventory.candidates_of(r.locator))
        records.append({"record": r.locator.to_dict(), "outcome": "excluded", "allele": None,
                        "exclusion": "duplicate_occurrence" if same else "alternative_for_selected_identity",
                        "related": entry.record.to_dict(), "detail": None})
    decided = {cid: i for i, d in enumerate(recipe.decisions) for cid in d.claims}
    selected_subjects = set(selected) | {c.id for c in chosen.values()}
    claims = []
    for c in inventory.claims:
        if c.subject.id not in selected_subjects:
            outcome, reason = "not_selected", None
        elif c.predicate not in m.RUNTIME_PREDICATES:
            outcome, reason = "not_representable", "the runtime reference has no such claim"
        elif c.id in decided and recipe.decisions[decided[c.id]].action == "withhold_claim":
            outcome, reason = "withheld", recipe.decisions[decided[c.id]].reason
        elif c.status == "unresolved":
            outcome, reason = "withheld", "unresolved"
        else:
            outcome, reason = "projected", None
        claims.append({"claim": c.id, "outcome": outcome, "decision": decided.get(c.id), "reason": reason})
    exports = [{"asset": a, **inventory.export_counts(a)} for a in ("airr", "airr_ex")]
    seg = {"V": 0, "D": 0, "J": 0, "unassigned": 0}
    for e in recipe.selection:
        seg[e.segment] += 1
    outcomes = {k: sum(r["outcome"] == k for r in records) for k in ("selected", "excluded", "rejected")}
    return {
        "artifact": "genairr_reference_source_build_report", "schema_version": 1,
        "generator": dict(GENERATOR), "source_lock": lock.to_dict(), "recipe": recipe.to_dict(),
        "inventory": inventory.to_dict(),
        "hashes": {"source_lock": m.lock_hash(lock), "inventory": m.inventory_hash(inventory),
                   "recipe": m.recipe_hash(recipe), "build_inputs": m.build_input_hash(lock, inventory, recipe)},
        "runtime": dict(RUNTIME),
        "comparison": {"complete": True, "exports": ["airr", "airr_ex"], "missing": [],
                       "findings": [_finding(inventory, "variant_coding_sequence_differs").id]},
        "counts": {"exports": exports, "selected_alleles": seg, "record_outcomes": outcomes},
        "selection": [{"allele": e.allele_id, "candidate": chosen[e.record.id].locator.to_dict()}
                      for e in recipe.selection],
        "records": records, "claims": claims, "capabilities": _capabilities(),
        "requirements": [{"capability": "observed_span_simulation", "evidence": "assertion", "satisfied": True}],
    }
@pytest.fixture(scope="module")
def report_d(inventory, inv):
    return report_dict(inventory, inv)
@pytest.fixture(scope="module")
def report(report_d):
    return m.BuildReport.from_dict(report_d)
FIXED = ("inputs/annotations.json", "inputs/d.fasta", "inputs/j.fasta", "inputs/v.fasta", "recipe.json",
         "reference.genairr-ref.json", "report.json", "source.lock.json")
def manifest_dict(report, *, d_fasta=True) -> dict:
    members = [{"name": name, "sha256": "sha256:" + format(i, "064x"), "bytes": 100 + i}
               for i, name in enumerate(FIXED) if d_fasta or name != "inputs/d.fasta"]
    for asset in report.source_lock.assets:
        members.append({"name": "cache/blobs/sha256/" + asset.blob_name, "sha256": asset.sha256, "bytes": asset.bytes})
    members.sort(key=lambda x: x["name"])
    hashes = report.to_dict()["hashes"]
    return {"artifact": "genairr_reference_source_build_manifest", "schema_version": 1,
            "builder": dict(GENERATOR),
            "versions": {"source_lock": 1, "recipe": 1, "report": 1, "reference": 1, "annotations": 1},
            "adapter": copy.deepcopy(ADAPTER), "projection": {"id": "observed_span", "version": 1},
            "hashes": {**hashes, "report": m.report_hash(report)}, "runtime": dict(RUNTIME), "members": members}
