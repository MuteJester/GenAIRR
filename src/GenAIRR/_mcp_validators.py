"""Per-record AIRR consistency checks used by validate_records.

Two tiers:

  Schema-only (run unconditionally):
    1. Required fields present: sequence, v_call, j_call, productive
    2. junction_start / junction_end within [0, sequence_length]
    3. junction_aa nucleotide-length math (when both junction + junction_aa present)
    4. productive=True => no '*' in junction_aa
    5. CIGAR strings (when present) only use M/I/D ops -- no S/H

  Refdata-driven (skipped when refdata is None):
    6. v_call / d_call / j_call alleles exist in the config's pools
    7. v_germline_end - v_germline_start <= claimed allele length
    8. anchor position from claimed v_call's allele is consistent with junction_start
    9. locus field matches the config's locus (skipped when config_name is None)

Returns a list of human-readable issue strings. Empty list = record is valid
under the active check set.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from GenAIRR._mcp_refdata import find_allele, locus_from_config_name


_REQUIRED_FIELDS = ("sequence", "v_call", "j_call", "productive")
_VALID_CIGAR_OPS = set("MID")
_CIGAR_FIELDS = ("v_cigar", "d_cigar", "j_cigar")

# Map the config-name locus token (e.g. 'TCRB' from 'human_tcrb') to the
# token the engine actually emits in record['locus'] (e.g. 'TRB'). The IG
# loci are emitted as-is (IGH, IGK, IGL); the TCR loci get rewritten from
# the IMGT-style 'TCRA/B/G/D' to the engine's 'TRA/B/G/D'.
_LOCUS_ALIAS = {
    "TCRA": "TRA",
    "TCRB": "TRB",
    "TCRG": "TRG",
    "TCRD": "TRD",
}


def _check_required_fields(rec: Dict[str, Any]) -> List[str]:
    return [f"missing required field: {f}" for f in _REQUIRED_FIELDS if f not in rec]


def _check_junction_bounds(rec: Dict[str, Any]) -> List[str]:
    seq_len = rec.get("sequence_length")
    if seq_len is None:
        return []
    issues = []
    for end_key in ("junction_start", "junction_end"):
        v = rec.get(end_key)
        if v is None:
            continue
        if not (0 <= v <= seq_len):
            issues.append(f"{end_key}={v} out of range [0, {seq_len}]")
    return issues


def _check_junction_aa_nt_length(rec: Dict[str, Any]) -> List[str]:
    junction = rec.get("junction")
    junction_aa = rec.get("junction_aa")
    if not junction or not junction_aa or not isinstance(junction, str):
        return []
    nt_len = len(junction)
    aa_len = len(junction_aa)
    if nt_len // 3 != aa_len:
        return [f"junction_aa length {aa_len} != junction nt length {nt_len} / 3"]
    return []


def _check_productive_implies_no_stop(rec: Dict[str, Any]) -> List[str]:
    if rec.get("productive") is not True:
        return []
    junction_aa = rec.get("junction_aa") or ""
    if "*" in junction_aa:
        return [f"productive=True but junction_aa contains stop codon: {junction_aa!r}"]
    return []


def _check_cigar_ops(rec: Dict[str, Any]) -> List[str]:
    issues = []
    for field in _CIGAR_FIELDS:
        cigar = rec.get(field) or ""
        if not cigar:
            continue
        # CIGAR is run-length-encoded as <count><op><count><op>... -- we only
        # need to verify the ops, so scan for any non-digit character.
        bad = sorted({c for c in cigar if not c.isdigit() and c not in _VALID_CIGAR_OPS})
        if bad:
            issues.append(
                f"{field} contains disallowed op(s) {bad} -- only M/I/D are valid"
            )
    return issues


def _check_alleles_exist(rec: Dict[str, Any], refdata: Any) -> List[str]:
    """Check 6: v_call / d_call / j_call alleles exist in their pools.

    Uses find_allele from _mcp_refdata; handles comma-separated tie sets.
    Missing pool (e.g. D on a VJ chain) is silently skipped -- the engine
    won't emit a d_call there in the first place.
    """
    issues: List[str] = []
    for segment, key in [("v", "v_call"), ("d", "d_call"), ("j", "j_call")]:
        call = rec.get(key)
        if not call:
            continue
        # Comma-separated tie sets -- every listed allele must exist.
        names = [c.strip() for c in str(call).split(",") if c.strip()]
        for name in names:
            if find_allele(refdata, segment, name) is None:
                issues.append(f"{key} allele not found in {segment} pool: {name}")
    return issues


def _check_v_germline_span_fits_allele(rec: Dict[str, Any], refdata: Any) -> List[str]:
    """Check 7: v_germline_end - v_germline_start <= claimed v_call allele length.

    Looks up the first v_call allele via find_allele; compares the engine's
    reported germline span against the allele's actual length. A span longer
    than the allele means the AIRR record is inconsistent with the germline
    pool it was supposedly drawn from.
    """
    v_call = rec.get("v_call")
    if not v_call:
        return []
    g_start = rec.get("v_germline_start")
    g_end = rec.get("v_germline_end")
    if g_start is None or g_end is None:
        return []
    first = str(v_call).split(",", 1)[0].strip()
    if not first:
        return []
    allele = find_allele(refdata, "v", first)
    if allele is None:
        # Already flagged by check 6; don't double-report.
        return []
    span = g_end - g_start
    raw_seq = allele.seq()
    allele_len = len(raw_seq)
    if span > allele_len:
        return [
            f"v_germline span {span} (= {g_end} - {g_start}) exceeds "
            f"claimed allele {first!r} length {allele_len}"
        ]
    return []


def _check_v_anchor_consistent_with_junction(
    rec: Dict[str, Any], refdata: Any
) -> List[str]:
    """Check 8: claimed v_call allele's anchor projected into pool space via
    v_sequence_start - v_germline_start + allele.anchor should land at
    junction_start (== the V Cys position).

    Skip when any of v_call, v_germline_start, v_sequence_start,
    allele.anchor, or junction_start are missing.
    """
    v_call = rec.get("v_call")
    if not v_call:
        return []
    v_seq_start = rec.get("v_sequence_start")
    v_germ_start = rec.get("v_germline_start")
    junction_start = rec.get("junction_start")
    if v_seq_start is None or v_germ_start is None or junction_start is None:
        return []
    first = str(v_call).split(",", 1)[0].strip()
    if not first:
        return []
    allele = find_allele(refdata, "v", first)
    if allele is None:
        return []
    anchor = getattr(allele, "anchor", None)
    if anchor is None:
        return []
    projected = v_seq_start - v_germ_start + anchor
    if projected != junction_start:
        return [
            f"v_call {first!r} anchor projected to {projected} "
            f"(= {v_seq_start} - {v_germ_start} + {anchor}) "
            f"does not match junction_start={junction_start}"
        ]
    return []


def _check_locus_matches(rec: Dict[str, Any], config_name: str) -> List[str]:
    """Check 9: rec['locus'] matches locus_from_config_name(config_name).

    Only fires when a record has a locus field. TCR loci get normalised
    from the alias-style 'TCRA/B/G/D' to the engine-emitted 'TRA/B/G/D'
    before comparison.
    """
    rec_locus = rec.get("locus")
    if not rec_locus:
        return []
    cfg_locus = locus_from_config_name(config_name)
    if not cfg_locus:
        return []
    expected = _LOCUS_ALIAS.get(cfg_locus, cfg_locus)
    if str(rec_locus).upper() != expected:
        return [
            f"locus={rec_locus!r} does not match config {config_name!r} "
            f"locus (expected {expected!r})"
        ]
    return []


def validate_one_record(
    rec: Dict[str, Any],
    refdata: Any = None,
    config_name: Optional[str] = None,
) -> List[str]:
    """Run every applicable check on a single record.

    Returns a list of human-readable issue strings. Empty list means the
    record passes all enabled checks. Refdata-driven checks (6-8) only run
    when refdata is supplied; check 9 (locus) only runs when config_name
    is supplied.
    """
    issues: List[str] = []
    issues.extend(_check_required_fields(rec))
    # If required fields are missing, downstream checks may KeyError --
    # short-circuit so we don't crash inside a single bad record.
    if any("missing required field" in i for i in issues):
        return issues

    issues.extend(_check_junction_bounds(rec))
    issues.extend(_check_junction_aa_nt_length(rec))
    issues.extend(_check_productive_implies_no_stop(rec))
    issues.extend(_check_cigar_ops(rec))

    if refdata is not None:
        issues.extend(_check_alleles_exist(rec, refdata))
        issues.extend(_check_v_germline_span_fits_allele(rec, refdata))
        issues.extend(_check_v_anchor_consistent_with_junction(rec, refdata))
    if config_name is not None:
        issues.extend(_check_locus_matches(rec, config_name))

    return issues
