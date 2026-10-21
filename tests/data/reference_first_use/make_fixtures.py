"""Write the first-use acceptance fixtures and the golden human IGH reference artifact.

One-off generator, committed with its outputs.  Run it from anywhere:

    python tests/data/reference_first_use/make_fixtures.py

Sequences, anchors and coding phases are copied from GenAIRR's bundled ``HUMAN_IGH_IMGT``
and ``HUMAN_IGK_IMGT`` cartridges.  These are software fixtures, not independent evidence:
every annotation is an assertion (``"support": "assertion"``) whose method names the
cartridge it was copied from, coding phases are ``anchor_derived``, and nothing claims
source support.  Each sidecar starts from ``Reference.annotation_template`` and is then
filled in, which is the workflow a user follows.  Sequence digests in the sidecars are over
the uppercase, ungapped sequence.

Outputs (all relative to ``tests/``):

* ``data/reference_first_use/human_igh/{v,d,j}.fasta`` and ``annotations.json`` -- vdj IGH,
  three V, two D and two J alleles.
* ``data/reference_first_use/human_igh_gapped/v.fasta`` and ``annotations.json`` -- the same
  V alleles, IMGT-gapped as distributed (lower case), under ``v_alignment: "imgt_gapped"``
  with no stated V anchors; D and J come from ``human_igh``.
* ``data/reference_first_use/igk_local/{v,j}.fasta`` and ``annotations.json`` -- vj IGK with
  local allele names, an exact duplicate record, a record with an unsupported base, a record
  whose label cannot be represented, a J without an anchor, and an alias-only name.
* ``golden/reference/human_igh.genairr-ref.json`` -- the artifact built from ``human_igh``.
"""
from __future__ import annotations

import json
from pathlib import Path

from GenAIRR import Reference, data

HERE = Path(__file__).resolve().parent
GOLDEN = HERE.parents[1] / "golden" / "reference" / "human_igh.genairr-ref.json"
BUILD_DATE = "2026-10-21"
HUMAN = {"scientific_name": "Homo sapiens", "taxonomy_id": 9606}

IGH_V = {"IGHV1-2*02": 285, "IGHV3-23*01": 285, "IGHV4-34*01": 282}
IGH_D = ("IGHD3-10*01", "IGHD2-2*01")
IGH_J = {"IGHJ4*02": 14, "IGHJ6*02": 29}
IGK_V_A, IGK_V_B = "IGKV1-39*01", "IGKV3-20*01"
IGK_J_1, IGK_J_2 = "IGKJ1*01", "IGKJ2*01"
KAPPA_V_A_ANCHOR, KAPPA_J_1_ANCHOR, KAPPA_J_2_CANDIDATE = 261, 7, 8
BAD_BASE_POSITION = 100


def _method(config: str) -> str:
    return f"software fixture copied from GenAIRR {config}; not independent evidence"


def _assertion(value, config: str) -> dict:
    return {"value": value, "support": "assertion", "source": None, "method": _method(config)}


def _anchor_phase(anchor: int, config: str) -> dict:
    return {**_assertion(anchor % 3, config), "derivation": "anchor_derived"}


def _bundled(config, segment: str, name: str):
    pool = {"V": config.v_alleles, "D": config.d_alleles, "J": config.j_alleles}[segment]
    for alleles in pool.values():
        for allele in alleles:
            if allele.name == name:
                return allele
    raise KeyError(f"{segment} allele {name} not found")


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _write_fasta(path: Path, records) -> None:
    _write_text(path, "".join(f">{header}\n{sequence}\n" for header, sequence in records))


def _write_json(path: Path, document: dict) -> None:
    _write_text(path, json.dumps(document, indent=2, ensure_ascii=False) + "\n")


def _header(name, version, locus, chain_type) -> dict:
    return {"name": name, "version": version, "locus": locus, "chain_type": chain_type,
            "taxon": dict(HUMAN), "build_date": BUILD_DATE}


def _fill(template: dict, annotations: dict, relationships=()) -> dict:
    """Fill a template's unknown annotations from ``{(segment, label): {name: annotation}}``."""
    seen = set()
    for entry in template["alleles"]:
        key = (entry["segment"], entry["label"])
        for name, value in annotations.get(key, {}).items():
            entry[name] = value
        seen.add(key)
    missing = set(annotations) - seen
    if missing:
        raise RuntimeError(f"annotations for alleles the template does not list: {sorted(missing)}")
    template["relationships"] = list(relationships)
    return template


def _igh_dj_annotations(igh) -> dict:
    config = "HUMAN_IGH_IMGT"
    out = {}
    for name in IGH_D:
        out[("D", name)] = {"functionality": _assertion("functional", config),
                            "boundary_5": _assertion("complete", config),
                            "boundary_3": _assertion("complete", config)}
    for name, anchor in IGH_J.items():
        if _bundled(igh, "J", name).anchor != anchor:
            raise RuntimeError(f"{name}: bundled anchor changed")
        out[("J", name)] = {"functionality": _assertion("functional", config),
                            "anchor": _assertion(anchor, config),
                            "coding_phase": _anchor_phase(anchor, config),
                            "boundary_5": _assertion("complete", config)}
    return out


def make_human_igh(igh) -> None:
    config = "HUMAN_IGH_IMGT"
    directory = HERE / "human_igh"
    _write_fasta(directory / "v.fasta", [(n, _bundled(igh, "V", n).ungapped_seq.upper()) for n in IGH_V])
    _write_fasta(directory / "d.fasta", [(n, _bundled(igh, "D", n).ungapped_seq.upper()) for n in IGH_D])
    _write_fasta(directory / "j.fasta", [(n, _bundled(igh, "J", n).ungapped_seq.upper()) for n in IGH_J])
    annotations = _igh_dj_annotations(igh)
    for name, anchor in IGH_V.items():
        if _bundled(igh, "V", name).anchor != anchor:
            raise RuntimeError(f"{name}: bundled anchor changed")
        annotations[("V", name)] = {"functionality": _assertion("functional", config),
                                    "anchor": _assertion(anchor, config),
                                    "coding_phase": _anchor_phase(anchor, config),
                                    "boundary_3": _assertion("complete", config)}
    template = Reference.annotation_template(
        v=directory / "v.fasta", d=directory / "d.fasta", j=directory / "j.fasta",
        reference=_header("human_igh_fixture", "1", "IGH", "vdj"))
    _write_json(directory / "annotations.json", _fill(template, annotations))


def make_human_igh_gapped(igh) -> None:
    config = "HUMAN_IGH_IMGT"
    directory = HERE / "human_igh_gapped"
    igh_dir = HERE / "human_igh"
    _write_fasta(directory / "v.fasta", [(n, _bundled(igh, "V", n).gapped_seq) for n in IGH_V])
    annotations = _igh_dj_annotations(igh)
    for name, anchor in IGH_V.items():
        # No stated anchor: the declared IMGT alignment supplies it.
        annotations[("V", name)] = {"functionality": _assertion("functional", config),
                                    "coding_phase": _anchor_phase(anchor, config),
                                    "boundary_3": _assertion("complete", config)}
    template = Reference.annotation_template(
        v=directory / "v.fasta", d=igh_dir / "d.fasta", j=igh_dir / "j.fasta",
        reference=_header("human_igh_gapped_fixture", "1", "IGH", "vdj"), v_alignment="imgt_gapped")
    _write_json(directory / "annotations.json", _fill(template, annotations))


def make_igk_local(igk) -> None:
    config = "HUMAN_IGK_IMGT"
    directory = HERE / "igk_local"
    v_a = _bundled(igk, "V", IGK_V_A)
    v_b = _bundled(igk, "V", IGK_V_B)
    j_1 = _bundled(igk, "J", IGK_J_1)
    j_2 = _bundled(igk, "J", IGK_J_2)
    if (v_a.anchor, j_1.anchor, j_2.anchor) != (KAPPA_V_A_ANCHOR, KAPPA_J_1_ANCHOR, KAPPA_J_2_CANDIDATE):
        raise RuntimeError("bundled IGK anchors changed")
    seq_a, seq_b = v_a.ungapped_seq.upper(), v_b.ungapped_seq.upper()
    bad = seq_a[:BAD_BASE_POSITION] + "N" + seq_a[BAD_BASE_POSITION + 1:]
    _write_fasta(directory / "v.fasta", [
        ("kappa_v_a", seq_a),
        ("kappa_v_b", seq_b),
        ("kappa_v_a", seq_a),       # an exact duplicate record
        ("kappa_v_bad", bad),       # unsupported base
        ("kappa,v", seq_b),         # label cannot be represented (contains ',')
    ])
    _write_fasta(directory / "j.fasta", [("kappa_j_1", j_1.ungapped_seq.upper()),
                                         ("kappa_j_2", j_2.ungapped_seq.upper())])
    annotations = {
        ("V", "kappa_v_a"): {"anchor": _assertion(KAPPA_V_A_ANCHOR, config),
                             "coding_phase": _anchor_phase(KAPPA_V_A_ANCHOR, config)},
        ("J", "kappa_j_1"): {"anchor": _assertion(KAPPA_J_1_ANCHOR, config),
                             "coding_phase": _anchor_phase(KAPPA_J_1_ANCHOR, config)},
    }
    relationships = [{"kind": "alias", "allele": "V:kappa_v_a", "entity": "lab_alias_kva", "source": None}]
    template = Reference.annotation_template(
        v=directory / "v.fasta", j=directory / "j.fasta",
        reference=_header("igk_local_fixture", None, "IGK", "vj"))
    candidates = {(s["allele"], s["value"]) for s in template["suggestions"]}
    if ("J:kappa_j_2", KAPPA_J_2_CANDIDATE) not in candidates:
        raise RuntimeError(f"expected a motif candidate at {KAPPA_J_2_CANDIDATE} for kappa_j_2: {candidates}")
    _write_json(directory / "annotations.json", _fill(template, annotations, relationships))


def make_golden() -> None:
    igh_dir = HERE / "human_igh"
    reference = Reference.from_fasta(v=igh_dir / "v.fasta", d=igh_dir / "d.fasta", j=igh_dir / "j.fasta",
                                     annotations=igh_dir / "annotations.json")
    GOLDEN.parent.mkdir(parents=True, exist_ok=True)
    reference.write(GOLDEN, replace=True)
    print(f"catalogue_hash = {reference.catalogue_hash}")
    print(f"execution_hash = {reference.execution_hash}")


def main() -> None:
    igh, igk = data.HUMAN_IGH_IMGT, data.HUMAN_IGK_IMGT
    make_human_igh(igh)
    make_human_igh_gapped(igh)
    make_igk_local(igk)
    make_golden()


if __name__ == "__main__":
    main()
