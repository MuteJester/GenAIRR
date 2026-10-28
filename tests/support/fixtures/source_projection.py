"""Recipes and projections over the pinned exports.

The starter recipe and the small two-allele recipe, plus the helpers that read a
projection's entries back.
"""
from __future__ import annotations

import hashlib
import json
from GenAIRR._reference_sources import evidence, ogrdb
from tests.support.fixtures.ogrdb_coordinates import V_LABEL
from tests.support.fixtures.ogrdb_inventory import (
    IGH, RECORDS, SELECTION, raw_records,
)

J1, J7 = "IGLJ1*01", "IGLJ7*02"

MINI = (V_LABEL, J1)

CHOICES = {row["label"]: row for row in SELECTION["record_choices"]}

STARTER_REQUIRES = (("observed_span_simulation", "assertion"), ("source_annotated_junction", "source"),
                    ("canonical_productive", "assertion"))

def sha(sequence: str) -> str:
    return "sha256:" + hashlib.sha256(sequence.encode("ascii")).hexdigest()

def records_in(blobs, fmt) -> list:
    """The records of one record export, read with the standard library."""
    return json.loads(blobs[fmt].decode("utf-8"))["GermlineSet"][0]["allele_descriptions"]

def entry(blobs, fmt, label, reason="selected for this test") -> dict:
    records = records_in(blobs, fmt)
    (index,) = [i for i, r in enumerate(records) if r["label"] == label]
    return {"segment": records[index]["sequence_type"], "label": label,
            "record": {"asset": fmt, "pointer": f"{RECORDS}/{index}"},
            "sequence_sha256": sha(records[index]["coding_sequence"]), "reason": reason}

def recipe(lock, selection, *, requires=(("observed_span_simulation", "assertion"),), decisions=(),
           name="test reference", **reference) -> dict:
    identity = lock.identity
    head = {"name": name, "version": "1", "taxon": identity.taxon.to_dict(), "locus": identity.locus,
            "chain_type": identity.chain_type, "build_date": "2026-10-27"}
    head.update(reference)
    return {"artifact": "genairr_reference_recipe", "schema_version": 1, "source_lock_hash": lock.commitment,
            "adapter": lock.adapter.to_dict(), "projection": {"id": "observed_span", "version": 1},
            "reference": head, "selection": list(selection), "decisions": list(decisions),
            "requires": [{"capability": c, "evidence": e} for c, e in requires],
            "defaults_profile": "genairr.unfitted_defaults.v1"}

def interpreted(blobs, lock):
    return evidence.interpret(lock, blobs, ogrdb.build_inventory(lock, blobs))

def _airr_index() -> dict:
    return {r["label"]: i for i, r in enumerate(raw_records(IGH, "airr"))}

def starter_recipe(lock, inventory) -> dict:
    """The selection fixture as a recipe: every selected allele from its own 'airr' record, in
    fixture order, and one select_alternative decision per recorded record choice."""
    index, records = _airr_index(), raw_records(IGH, "airr")
    selection = []
    for segment in ("V", "D", "J"):
        for label in SELECTION["selected"][segment]:
            reason = SELECTION["record_choice_reason"] if label in CHOICES else "selected by the starter rules"
            selection.append({"segment": segment, "label": label,
                              "record": {"asset": "airr", "pointer": f"{RECORDS}/{index[label]}"},
                              "sequence_sha256": sha(records[index[label]]["coding_sequence"]), "reason": reason})
    differs = {f.observations.get()["label"]: f.id for f in inventory.findings
               if f.code == "export_coding_sequence_differs"}
    decisions = [{"action": "select_alternative", "findings": [differs[label]], "claims": [],
                  "sequence_sha256": sha(records[index[label]]["coding_sequence"]),
                  "reason": SELECTION["record_choice_reason"], "evidence": []} for label in CHOICES]
    return recipe(lock, selection, requires=STARTER_REQUIRES, decisions=decisions,
                  name="human IGH starter (OGRDB IGH_VDJ release 10)")
