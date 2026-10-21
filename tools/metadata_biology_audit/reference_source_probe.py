"""Reproduce a bounded, read-only study of pinned OGRDB IGL release 4 data.

This research probe supports the reference-curation design. It does not import,
repair, certify, or compile references, and it is not a general AIRR validator.
Downloaded source bytes remain in the explicitly selected local cache. Network
access occurs only with --fetch; subsequent runs can reproduce the observations
offline using the same cache. Only the Python standard library is required.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import urllib.request


MAX_SOURCE_BYTES = 16 * 1024 * 1024
ROOT_URL = "https://ogrdb.airr-community.org/api_v2/germline/set/9606.IGLambda_VJ/4"
SOURCE_HASH = "48aaf90d81377836901dda1013f4244cb00ecf1437092ab76d764534d647e9c5"
ANALYSIS_HASH = "9e69fdb59d28cdd133a4ca502a87135f73e3dbd90a8d46eb82af98d36fb138fa"
SOURCES = {
    "source": (ROOT_URL + "/airr_ex", "igl_v4_airr_ex.json", SOURCE_HASH),
    "analysis": (ROOT_URL + "/airr", "igl_v4_airr.json", ANALYSIS_HASH),
    "noformat": (ROOT_URL, "igl_v4_noformat.json", SOURCE_HASH),
}
DOCUMENTATION = {
    "ogrdb_formats": "https://wordpress.vdjbase.org/index.php/ogrdb/explanation-of-germline-set-formats/",
    "ogrdb_api": "https://github.com/airr-community/ogrdb/blob/master/schema/ogrdb_api_v2_guide.md",
    "airr_germline_schema": "https://docs.airr-community.org/en/stable/datarep/germline.html",
    "ogrdb_data_license": "https://ogrdb.airr-community.org/render_page/licensing_statement.html",
    "producer_j_conventions": (
        "https://github.com/williamdlees/receptor_utils/blob/"
        "e7983d47e6662449db100940c28ef83491c4f51e/"
        "src/receptor_utils/aux_formats.py#L163"
    ),
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_json(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False).encode("utf-8")


def verify_source(data: bytes, expected: str, label: str) -> None:
    if len(data) > MAX_SOURCE_BYTES:
        raise ValueError(f"{label}: source exceeds the 16 MiB limit")
    actual = digest(data)
    if actual != expected:
        raise ValueError(
            f"{label}: pinned SHA256 mismatch; expected {expected}, got {actual}. "
            "Do not update the pin without reviewing changed upstream content."
        )


def load_sources(cache: Path, fetch: bool) -> tuple[dict, dict]:
    """Acquire or read the exact pinned bytes; acquisition time is separate."""
    data_by_name = {}
    receipt_path = cache / "reference_source_probe_acquisition.json"
    receipt = {}
    if fetch:
        for name, (url, _, expected) in SOURCES.items():
            request = urllib.request.Request(
                url, headers={"User-Agent": "GenAIRR-reference-source-research-probe/1"}
            )
            with urllib.request.urlopen(request, timeout=30) as response:
                data = response.read(MAX_SOURCE_BYTES + 1)
                content_type = response.headers.get("Content-Type")
                final_url = response.geturl()
            verify_source(data, expected, name)
            data_by_name[name] = data
            receipt[name] = {
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                "url": url,
                "final_url": final_url,
                "content_type": content_type,
                "sha256": expected,
                "bytes": len(data),
            }
        # A failed download/hash check never replaces any existing cache files.
        cache.mkdir(parents=True, exist_ok=True)
        for name, (_, filename, _) in SOURCES.items():
            (cache / filename).write_bytes(data_by_name[name])
        receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    else:
        for name, (_, filename, expected) in SOURCES.items():
            with (cache / filename).open("rb") as stream:
                data = stream.read(MAX_SOURCE_BYTES + 1)
            verify_source(data, expected, name)
            data_by_name[name] = data
        if receipt_path.exists():
            with receipt_path.open("rb") as stream:
                raw_receipt = stream.read(65537)
            if len(raw_receipt) > 65536:
                raise ValueError("cache acquisition receipt exceeds 64 KiB")
            receipt = json.loads(raw_receipt)
            for name, (url, _, expected) in SOURCES.items():
                item = receipt.get(name, {})
                if (item.get("url") != url or item.get("sha256") != expected
                        or item.get("bytes") != len(data_by_name[name])):
                    raise ValueError(f"{name}: acquisition receipt does not match pinned source")
    return data_by_name, receipt


def v_boundary_observation(allele: dict) -> dict:
    """Compare coordinate claims inside this export, without choosing a boundary."""
    start, end = allele["gene_start"], allele["gene_end"]
    sequence, coding = allele["sequence"], allele["coding_sequence"]
    rss_start = allele.get("v_rs_start")
    rss_index = rss_start - 1 if isinstance(rss_start, int) else None
    coding_end = start - 1 + len(coding)
    return {
        "label": allele["label"],
        "gene_start": start, "gene_end": end, "v_rs_start": rss_start,
        "coordinate_convention_tested": "1-based inclusive in record sequence",
        "declared_gene_end_to_rss_gap_nt": rss_index - end if rss_index is not None else None,
        "coding_length_implied_end_to_rss_gap_nt": rss_index - coding_end if rss_index is not None else None,
        "declared_slice_equals_coding": sequence[start - 1:end] == coding,
        "coding_equals_record_slice_of_same_length": sequence[start - 1:coding_end] == coding,
        "coding_equals_record_slice_up_to_rss": sequence[start - 1:rss_index] == coding
            if rss_index is not None else None,
        "bases_between_declared_end_and_rss": sequence[end:rss_index]
            if rss_index is not None and end <= rss_index else None,
        "record_rss_first_7_nt": sequence[rss_index:rss_index + 7]
            if rss_index is not None else None,
    }


def sibling_tail_observations(allele: dict, alleles: list[dict], extension_length: int) -> dict:
    # This bounded probe uses the supplied OGRDB labels to group siblings. This
    # is not a general gene-name parser or evidence of distinct genomic origin.
    gene_label = allele["label"].split("*", 1)[0]
    siblings = []
    for other in sorted(alleles, key=lambda item: item["label"]):
        if other["label"] == allele["label"] or other["label"].split("*", 1)[0] != gene_label:
            continue
        terminal_length = min(extension_length, len(other["coding_sequence"]))
        sibling_boundary = v_boundary_observation(other)
        siblings.append({
            "label": other["label"],
            "coding_terminal_nt": other["coding_sequence"][-terminal_length:]
                if terminal_length else "",
            "terminal_extension_matches": other["coding_sequence"].endswith(
                allele["coding_sequence"][-extension_length:]) if extension_length else None,
            "coding_last_12": other["coding_sequence"][-12:],
            "last_12_matches": other["coding_sequence"][-12:] == allele["coding_sequence"][-12:],
            "coding_equals_own_record_slice_up_to_rss": sibling_boundary["coding_equals_record_slice_up_to_rss"],
            "embedded_unrearranged_support_count": len(other.get("unrearranged_support") or []),
        })
    return {
        "grouping_rule": "same label before first '*'; supplied OGRDB labels only",
        "gene_label": gene_label,
        "terminal_comparison_length_nt": extension_length,
        "siblings": siblings,
        "sibling_count": len(siblings),
        "all_siblings_terminal_extension_matches": all(x["terminal_extension_matches"] for x in siblings)
            if siblings else None,
        "all_siblings_last_12_matches": all(x["last_12_matches"] for x in siblings)
            if siblings else None,
        "evidence_scope": (
            "Within the same pinned export only. Embedded support counts report supplied "
            "entries, not independent verification, accession retrieval, or support quality. "
            "No independent genomic source was acquired for any sibling."
        ),
    }


def annotation_context(allele: dict) -> dict:
    """Retain historical affirmed-span notes alongside present coordinate claims."""
    context = {key: allele.get(key) for key in (
        "allele_description_id", "allele_description_ref", "release_version",
        "release_date", "inference_type", "status", "curational_tags", "curation",
    )}
    context["embedded_unrearranged_support"] = [
        {key: item.get(key) for key in (
            "sequence_id", "repository_name", "repository_ref", "patch_no",
            "gff_seqid", "gff_start", "gff_end", "strand", "curation",
        )}
        for item in allele.get("unrearranged_support") or []
    ]
    context["scope"] = (
        "Unmodified annotation text and support references in the pinned OGRDB record. "
        "Historical affirmation extent can differ from supplied genomic context or "
        "current coding_sequence extent; no external accession was retrieved here."
    )
    return context


def coding_interval_disagreements(alleles: list[dict]) -> list[dict]:
    result = []
    for allele in sorted(alleles, key=lambda item: item["label"]):
        start, end = allele["gene_start"], allele["gene_end"]
        coding = allele["coding_sequence"]
        declared_slice = allele["sequence"][start - 1:end]
        if declared_slice != coding:
            observation = {
                "label": allele["label"], "gene_start": start, "gene_end": end,
                "interval_convention_tested": "1-based inclusive in sequence",
                "declared_interval_length": end - start + 1,
                "extracted_slice_length": len(declared_slice),
                "coding_sequence_length": len(coding),
                "slice_last_12": declared_slice[-12:],
                "coding_last_12": coding[-12:],
                "curational_tags": allele.get("curational_tags"),
            }
            if allele["sequence_type"] == "V":
                boundary = v_boundary_observation(allele)
                extends_slice = coding.startswith(declared_slice) and len(coding) > len(declared_slice)
                extension = coding[len(declared_slice):] if extends_slice else None
                if extends_slice and boundary["coding_equals_record_slice_up_to_rss"]:
                    code = "coding_end_before_declared_v_rss"
                elif extends_slice and not boundary["coding_equals_record_slice_of_same_length"]:
                    code = "coding_extension_disagrees_with_record_context"
                else:
                    code = "coding_interval_sequence_disagreement"
                extension_length = len(extension) if extension is not None else 0
                observation.update({
                    "finding_code": code,
                    "coding_starts_with_declared_slice": extends_slice,
                    "coding_extension_nt": extension,
                    "record_sequence_after_declared_end_same_length":
                        allele["sequence"][end:end + extension_length],
                    "boundary_observation": boundary,
                    "same_gene_sibling_comparisons": sibling_tail_observations(
                        allele, alleles, extension_length),
                    "annotation_context": annotation_context(allele),
                    "evidence_scope": "Internal consistency of supplied sequence/coordinates/annotations; not independent boundary proof",
                    "automatic_repair_performed": False,
                })
            result.append(observation)
    return result


def v_boundary_summary(alleles: list[dict]) -> dict:
    records = [v_boundary_observation(allele)
               for allele in sorted(alleles, key=lambda item: item["label"])
               if allele["sequence_type"] == "V"]
    return {
        "record_count": len(records),
        "distance_definition": "(v_rs_start - 1) - end; positive means bases before RSS, negative means overlap",
        "declared_gene_end_to_rss_gap_distribution": [
            {"gap_nt": gap, "count": count}
            for gap, count in sorted(Counter(x["declared_gene_end_to_rss_gap_nt"] for x in records).items(),
                                     key=lambda item: (item[0] is None, item[0] or 0))
        ],
        "coding_length_implied_end_to_rss_gap_distribution": [
            {"gap_nt": gap, "count": count}
            for gap, count in sorted(Counter(x["coding_length_implied_end_to_rss_gap_nt"] for x in records).items(),
                                     key=lambda item: (item[0] is None, item[0] or 0))
        ],
        "records": records,
    }


def paralog_observations(indexed: dict[str, dict]) -> dict:
    result = {}
    for name in ("source", "analysis"):
        records = indexed[name]
        mentions = []
        identities = set(records)
        for label, allele in sorted(records.items()):
            for paralog in sorted(allele.get("paralogs") or []):
                identities.add(paralog)
                mentions.append({
                    "finding_code": "paralog_named_without_own_record"
                        if paralog not in records else None,
                    "record_label": label,
                    "paralog_label": paralog,
                    "paralog_has_own_record": paralog in records,
                    "present_as_record_in_other_variant": paralog in indexed[
                        "analysis" if name == "source" else "source"],
                })
        result[name] = {
            "record_count": len(records),
            "record_plus_paralog_identity_count": len(identities),
            "paralog_only_labels": sorted(identities - records.keys()),
            "mentions": mentions,
        }
    source_allele = indexed["source"]["IGLJ2*01"]
    analysis_allele = indexed["analysis"]["IGLJ3*01"]
    result["iglj2_iglj3_relationship"] = {
        "source_record_label": source_allele["label"],
        "source_paralog_label": "IGLJ3*01",
        "analysis_record_label": analysis_allele["label"],
        "source_allele_description_id": source_allele["allele_description_id"],
        "analysis_allele_description_id": analysis_allele["allele_description_id"],
        "source_record_annotation_context": annotation_context(source_allele),
        "analysis_record_annotation_context": annotation_context(analysis_allele),
        "full_sequence_bytes_equal": source_allele["sequence"] == analysis_allele["sequence"],
        "coding_sequence_bytes_equal": source_allele["coding_sequence"] == analysis_allele["coding_sequence"],
        "full_record_equal": source_allele == analysis_allele,
        "differing_record_fields": [key for key in sorted(source_allele.keys() | analysis_allele.keys())
                                    if source_allele.get(key) != analysis_allele.get(key)],
        "interpretation": (
            "Identical supplied sequence strings do not collapse distinct gene/allele identities. "
            "The Source variant names IGLJ3*01 only as a paralog; its own genomic annotations "
            "must not be copied from IGLJ2*01 without scoped evidence."
        ),
    }
    return result


def identical_sequence_groups(alleles: list[dict], field: str) -> list[dict]:
    groups = {}
    for allele in alleles:
        groups.setdefault(allele[field], []).append(allele["label"])
    return sorted(({
        "labels": sorted(labels), "sequence_length": len(sequence),
        "sequence_sha256": digest(sequence.encode("ascii")),
    } for sequence, labels in groups.items() if len(labels) > 1), key=lambda item: item["labels"])


def j_annotation_observations(alleles: list[dict]) -> list[dict]:
    result = []
    for allele in sorted(alleles, key=lambda item: item["label"]):
        if allele["sequence_type"] != "J":
            continue
        anchor = allele["j_cdr3_end"]
        frame = allele["j_codon_frame"]
        coding_index = anchor - 1
        full_triplet = allele["sequence"][coding_index:coding_index + 3]
        coding_triplet = allele["coding_sequence"][coding_index:coding_index + 3]
        # AIRR 2.0 prose defines the codon position of the first coding base.
        # The alternative describes the pinned producer utility's convention.
        airr_offset = (1 - frame) % 3
        producer_offset = frame - 1
        result.append({
            "label": allele["label"],
            "sequence_length": len(allele["sequence"]),
            "coding_sequence_length": len(allele["coding_sequence"]),
            "gene_start": allele["gene_start"], "gene_end": allele["gene_end"],
            "j_cdr3_end": anchor, "j_codon_frame": frame,
            "triplet_if_anchor_is_sequence_relative_1based": full_triplet,
            "triplet_if_anchor_is_coding_relative_1based": coding_triplet,
            "sequence_relative_triplet_is_phe_or_trp": full_triplet in ("TTT", "TTC", "TGG"),
            "coding_relative_triplet_is_phe_or_trp": coding_triplet in ("TTT", "TTC", "TGG"),
            "airr_prose_first_complete_codon_offset_0based": airr_offset,
            "producer_convention_first_complete_codon_offset_0based": producer_offset,
            "coding_relative_anchor_aligned_with_airr_prose": (coding_index - airr_offset) % 3 == 0,
            "coding_relative_anchor_aligned_with_producer_convention": (coding_index - producer_offset) % 3 == 0,
            "curational_tags": allele.get("curational_tags"),
        })
    return result


def provider_label_observations(alleles: list[dict]) -> dict:
    """Retain the producer's record labels without treating them as evidence."""
    records = [{
        "label": allele["label"],
        "allele_description_ref": allele.get("allele_description_ref"),
        "allele_description_ref_ends_with_ext": isinstance(
            allele.get("allele_description_ref"), str
        ) and allele["allele_description_ref"].endswith("_ext"),
    } for allele in sorted(alleles, key=lambda item: item["label"])]
    return {
        "record_count": len(records),
        "records_with_terminal_ext_label": sum(
            record["allele_description_ref_ends_with_ext"] for record in records
        ),
        "terminal_label_tested": "_ext",
        "records": records,
        "evidence_scope": (
            "Producer-supplied allele_description_ref strings only. The terminal _ext "
            "signal must be considered alongside endpoint documentation and sequence "
            "differences; it does not prove any particular extension or coding boundary."
        ),
    }


def v_alignment_coverage_observations(alleles: list[dict]) -> dict:
    """Report supplied gapped assets; their absence leaves span mapping unknown."""
    records = []
    for allele in sorted(alleles, key=lambda item: item["label"]):
        if allele["sequence_type"] != "V":
            continue
        delineations = []
        for delineation in allele.get("v_gene_delineations") or []:
            aligned = delineation.get("aligned_sequence")
            if "aligned_sequence" not in delineation:
                state = "missing_field"
            elif aligned is None:
                state = "null"
            elif isinstance(aligned, str):
                state = "nonempty_string" if aligned else "empty_string"
            else:
                state = "other_supplied_type"
            delineations.append({
                "sequence_delineation_id": delineation.get("sequence_delineation_id"),
                "delineation_scheme": delineation.get("delineation_scheme"),
                "aligned_sequence_state": state,
                "aligned_sequence_length": len(aligned) if isinstance(aligned, str) else None,
            })
        records.append({
            "label": allele["label"],
            "delineations": delineations,
            "has_nonempty_aligned_sequence": any(
                item["aligned_sequence_state"] == "nonempty_string"
                for item in delineations
            ),
        })
    supplied_count = sum(record["has_nonempty_aligned_sequence"] for record in records)
    return {
        "v_record_count": len(records),
        "v_records_with_nonempty_aligned_sequence": supplied_count,
        "v_records_without_nonempty_aligned_sequence": len(records) - supplied_count,
        "delineation_aligned_sequence_state_counts": dict(sorted(Counter(
            item["aligned_sequence_state"]
            for record in records for item in record["delineations"]
        ).items())),
        "records": records,
        "interpretation": (
            "Missing, null or empty aligned_sequence supplies no gapped sequence asset "
            "for mapping the historical IMGT-gapped affirmation spans. This mapping "
            "remains unknown here and requires a separately pinned compatible gapped "
            "asset; absence alone is not a finding that an allele is biologically invalid. "
            "Delineation coordinate fields alone do not reconstruct the omitted alignment."
        ),
    }


def observations(data_by_name: dict[str, bytes]) -> dict:
    sets = {name: json.loads(raw)["GermlineSet"][0]
            for name, raw in data_by_name.items()}
    alleles = {name: value["allele_descriptions"] for name, value in sets.items()}
    indexed = {name: {a["label"]: a for a in records}
               for name, records in alleles.items()}
    source, analysis = indexed["source"], indexed["analysis"]
    source_disagreements = coding_interval_disagreements(alleles["source"])
    disagreement_by_label = {item["label"]: item for item in source_disagreements}
    coding_changes = []
    for label in sorted(source.keys() & analysis.keys()):
        before, after = source[label], analysis[label]
        if before["coding_sequence"] == after["coding_sequence"]:
            continue
        coding_changes.append({
            "label": label,
            "source_coding_length": len(before["coding_sequence"]),
            "analysis_coding_length": len(after["coding_sequence"]),
            "source_coding_last_12": before["coding_sequence"][-12:],
            "analysis_coding_last_12": after["coding_sequence"][-12:],
            "sequence_identical": before["sequence"] == after["sequence"],
            "gene_coordinates_identical": (before["gene_start"], before["gene_end"])
                == (after["gene_start"], after["gene_end"]),
            "source_curational_tags": before.get("curational_tags"),
            "analysis_curational_tags": after.get("curational_tags"),
            "analysis_coding_equals_source_declared_slice": after["coding_sequence"]
                == before["sequence"][before["gene_start"] - 1:before["gene_end"]],
            "source_coding_starts_with_analysis_coding": before["coding_sequence"].startswith(after["coding_sequence"]),
            "removed_source_coding_suffix": before["coding_sequence"][len(after["coding_sequence"]):]
                if before["coding_sequence"].startswith(after["coding_sequence"]) else None,
            "source_finding_code": disagreement_by_label.get(label, {}).get("finding_code"),
            "analysis_boundary_observation": v_boundary_observation(after)
                if after["sequence_type"] == "V" else None,
        })
    analysis_disagreements = coding_interval_disagreements(alleles["analysis"])
    source_j = j_annotation_observations(alleles["source"])
    paralogs = paralog_observations(indexed)
    provider_labels = {name: provider_label_observations(alleles[name])
                       for name in ("source", "analysis")}
    v_alignments = {name: v_alignment_coverage_observations(alleles[name])
                    for name in ("source", "analysis")}
    result = {
        "source_equals_noformat_bytes": data_by_name["source"] == data_by_name["noformat"],
        "source_equals_analysis_bytes": data_by_name["source"] == data_by_name["analysis"],
        "source_label_interpretation": (
            "'source' and 'analysis' are endpoint/documentation labels, not guarantees "
            "of absence/presence of extensions or duplicate identities. The pinned "
            "airr_ex payload has both extension classes and a paralog-only identity; "
            "airr trims the six coding sequences to their declared intervals and "
            "contains a separate IGLJ3*01 record. Neither representation is certified here."
        ),
        "counts": {
            "source_alleles": len(alleles["source"]),
            "analysis_alleles": len(alleles["analysis"]),
            "source_coding_interval_disagreements": len(source_disagreements),
            "analysis_coding_interval_disagreements": len(analysis_disagreements),
            "coding_sequence_differences_between_variants": len(coding_changes),
            "source_j_alleles": len(source_j),
            "source_j_with_coding_relative_phe_trp": sum(x["coding_relative_triplet_is_phe_or_trp"] for x in source_j),
            "source_j_with_sequence_relative_phe_trp": sum(x["sequence_relative_triplet_is_phe_or_trp"] for x in source_j),
            "source_j_coding_anchor_aligned_with_airr_prose": sum(x["coding_relative_anchor_aligned_with_airr_prose"] for x in source_j),
            "source_j_coding_anchor_aligned_with_producer_convention": sum(x["coding_relative_anchor_aligned_with_producer_convention"] for x in source_j),
            "source_records_plus_paralog_identities": paralogs["source"]["record_plus_paralog_identity_count"],
            "analysis_records_plus_paralog_identities": paralogs["analysis"]["record_plus_paralog_identity_count"],
            "source_finding_codes": dict(sorted(Counter(
                item["finding_code"] for item in source_disagreements).items())),
            "source_allele_description_refs_with_terminal_ext":
                provider_labels["source"]["records_with_terminal_ext_label"],
            "analysis_allele_description_refs_with_terminal_ext":
                provider_labels["analysis"]["records_with_terminal_ext_label"],
            "source_v_records": v_alignments["source"]["v_record_count"],
            "analysis_v_records": v_alignments["analysis"]["v_record_count"],
            "source_v_records_with_nonempty_aligned_sequence":
                v_alignments["source"]["v_records_with_nonempty_aligned_sequence"],
            "analysis_v_records_with_nonempty_aligned_sequence":
                v_alignments["analysis"]["v_records_with_nonempty_aligned_sequence"],
        },
        "duplicate_labels": {
            name: {label: count for label, count in sorted(Counter(a["label"] for a in alleles[name]).items()) if count > 1}
            for name in ("source", "analysis")
        },
        "source_only_labels": sorted(source.keys() - analysis.keys()),
        "analysis_only_labels": sorted(analysis.keys() - source.keys()),
        "source_coding_interval_disagreements": source_disagreements,
        "analysis_coding_interval_disagreements": analysis_disagreements,
        "coding_sequence_variant_differences": coding_changes,
        "provider_record_label_observations": provider_labels,
        "v_alignment_coverage_observations": v_alignments,
        "source_j_annotations": source_j,
        "analysis_j_annotations": j_annotation_observations(alleles["analysis"]),
        "paralog_example": {"label": "IGLJ2*01", "source_paralogs": source["IGLJ2*01"].get("paralogs")},
        "paralog_identity_observations": paralogs,
        "v_boundary_observations": {name: v_boundary_summary(alleles[name])
                                    for name in ("source", "analysis")},
        "identical_sequence_groups": {
            name: {field: identical_sequence_groups(alleles[name], field)
                   for field in ("sequence", "coding_sequence")}
            for name in ("source", "analysis")
        },
        "api_metadata": {
            name: {"info": json.loads(data_by_name[name]).get("Info"),
                   "set_name": sets[name]["germline_set_name"],
                   "set_ref": sets[name]["germline_set_ref"],
                   "release_version": sets[name]["release_version"],
                   "release_date": sets[name]["release_date"]}
            for name in ("source", "analysis")
        },
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Fetch pinned bytes and write a report:\n"
            "  python3 tools/metadata_biology_audit/reference_source_probe.py "
            "--fetch --cache /tmp/genairr-reference-source-probe --output probe.json\n\n"
            "Reproduce offline (no network):\n"
            "  python3 tools/metadata_biology_audit/reference_source_probe.py "
            "--cache /tmp/genairr-reference-source-probe --output offline.json\n\n"
            "With the same acquisition receipt, both reports are byte-identical. "
            "findings_sha256 excludes retrieval times, which belong to acquisition only."
        ),
    )
    parser.add_argument("--fetch", action="store_true", help="explicitly fetch pinned upstream sources")
    parser.add_argument("--cache", required=True, type=Path, help="directory holding raw pinned JSON and acquisition receipt")
    parser.add_argument("--output", required=True, type=Path, help="write compact research observations as JSON")
    args = parser.parse_args()
    try:
        raw, acquisition = load_sources(args.cache, args.fetch)
        findings = observations(raw)
        report = {
            "schema": "genairr_reference_source_research_probe_v2",
            "schema_change": (
                "v2 retains v1 observations and adds per-record endpoint comparisons, "
                "two extension finding classes, raw affirmation/support annotations, "
                "sibling-tail comparisons, all-V RSS gaps including missing RSS, "
                "duplicate sequence groups, and paralog-only identity accounting. "
                "Additive v2 observations also retain provider _ext record labels and "
                "V aligned_sequence coverage, distinguishing absent gapped assets "
                "from biological invalidity."
            ),
            "scope": "Pinned OGRDB human IGL release 4 representations; design evidence only",
            "interpretation": (
                "These observations compare exported fields and coordinate conventions. "
                "They do not establish which biological sequence or boundary is correct, "
                "certify references, or authorize automatic repair. Motif consistency is a "
                "diagnostic check, not independent biological evidence. Supplied genomic "
                "slices, RSS coordinates and sibling sequences share this source and do not "
                "settle the scope of historically affirmed bases. Endpoint labels confer no "
                "per-record completeness or de-duplication guarantee."
            ),
            "documentation_urls": DOCUMENTATION,
            "sources": {name: {"url": url, "sha256": expected, "bytes": len(raw[name])}
                        for name, (url, _, expected) in SOURCES.items()},
            "acquisition": acquisition,
            "findings_sha256": digest(canonical_json(findings)),
            "findings": findings,
        }
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
                               + "\n", encoding="utf-8")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"reference source probe failed: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote {args.output}; findings SHA256 {report['findings_sha256']}")
    print(json.dumps(findings["counts"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
