"""Native reference annotation entry points shared by authoring consumers.

A matching codon is not proof of a conserved landmark or a complete coding end.
The optional evidence document records exactly which claim the author supplied.
"""
import hashlib
import json


def sequence_sha256(sequence):
    data = sequence.encode("ascii") if isinstance(sequence, str) else bytes(sequence)
    return "sha256:" + hashlib.sha256(data).hexdigest()


def resolve_allele_anchor(allele, segment, *, imgt_gapped=False):
    from . import _engine
    resolver = getattr(_engine, "resolve_reference_anchor", None)
    if resolver is None:
        raise ImportError("The installed GenAIRR native engine lacks the reference resolver; rebuild or reinstall the matching GenAIRR version")
    result = json.loads(resolver(
        allele.ungapped_seq.encode("ascii"), segment,
        locus=allele.name[:3].upper(),
        gapped_sequence=allele.gapped_seq.encode("ascii") if imgt_gapped else None,
        imgt_gapped=imgt_gapped,
    ))
    allele.anchor_meta = result
    allele.anchor = result["position"]
    allele.reference_evidence = {
        "schema_version": 1,
        "sequence_sha256": sequence_sha256(allele.ungapped_seq),
        # A landmark-derived modulo-three value is not independent frame evidence.
        "coding_phase": None,
        "anchor_status": result["status"],
        "boundary_5": "unknown", "boundary_3": "unknown",
        "source": None, "method": result["method"],
        "notes": [result["reason"]],
    }
    return result


def annotate_fasta_header(allele, raw_header, segment):
    """Retain IMGT header declarations without upgrading missing evidence.

    An IMGT-shaped header identifies its coordinate convention; a bare FASTA
    sequence does not. Partial/completeness flags are preserved in the raw
    header but not guessed into per-end coding-boundary claims.
    """
    fields = raw_header.lstrip(">").strip().split("|")
    is_imgt = len(fields) >= 5 and fields[4].strip().upper() == f"{segment}-REGION"
    if is_imgt:
        allele.species = fields[2].strip() or None
        declared = fields[3].strip().strip("()[]").upper()
        allele.functional_status = {"F": "functional", "ORF": "orf", "P": "pseudogene"}.get(declared, "unknown")
        allele.source = fields[0].strip() or None
    else:
        allele.functional_status = "unknown"
    if segment in ("V", "J"):
        result = resolve_allele_anchor(allele, segment, imgt_gapped=is_imgt and segment == "V")
        if segment == "J":
            allele.frame = None if result["position"] is None else result["position"] % 3
    else:
        allele.reference_evidence = {
            "schema_version": 1, "sequence_sha256": sequence_sha256(allele.ungapped_seq),
            "coding_phase": None, "anchor_status": "unknown",
            "boundary_5": "unknown", "boundary_3": "unknown",
            "source": None, "method": "fasta_header", "notes": [],
        }
    allele.reference_evidence["source"] = allele.source
    allele.reference_evidence["notes"].append("source_header=" + raw_header.lstrip(">"))
