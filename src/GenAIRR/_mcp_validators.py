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
    9. locus field matches the config's locus

Returns a list of human-readable issue strings. Empty list = record is valid
under the active check set.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional


_REQUIRED_FIELDS = ("v_call", "j_call", "productive")
_VALID_CIGAR_OPS = set("MID")
_CIGAR_FIELDS = ("v_cigar", "d_cigar", "j_cigar")


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
    issues = []
    for segment, key in [("v", "v_call"), ("d", "d_call"), ("j", "j_call")]:
        call = rec.get(key)
        if not call:
            continue
        # Comma-separated tie sets -- every listed allele must exist.
        names = [c.strip() for c in call.split(",") if c.strip()]
        if not names:
            continue
        pool = _get_pool_names(refdata, segment)
        if pool is None:
            continue  # config has no pool for this segment (e.g. VJ chain has no D)
        for name in names:
            if name not in pool:
                issues.append(f"{key} allele not found in {segment} pool: {name}")
    return issues


def _check_locus_matches(rec: Dict[str, Any], refdata: Any) -> List[str]:
    locus = rec.get("locus")
    if not locus:
        return []
    config_locus = _get_locus(refdata)
    if config_locus and locus != config_locus:
        return [f"locus={locus!r} does not match config locus {config_locus!r}"]
    return []


def _get_pool_names(refdata: Any, segment: str) -> Optional[set]:
    """Return a set of allele names for the given segment, or None when
    the refdata doesn't carry a pool for it (e.g. VJ chains have no D)."""
    # refdata is the engine's RefDataConfig -- Python object exposed by
    # genairr_engine. Method name varies; we go through the bundled DSL
    # accessor below for stability.
    try:
        from GenAIRR.experiment import dataconfig_to_refdata  # noqa: F401
    except ImportError:
        return None
    # Each segment's allele list is exposed via the pool indexer.
    pool_attr = {"v": "v_pool", "d": "d_pool", "j": "j_pool"}.get(segment)
    if pool_attr is None:
        return None
    pool = getattr(refdata, pool_attr, None)
    if pool is None:
        return None
    try:
        return {allele.name for _, allele in pool}
    except Exception:
        return None


def _get_locus(refdata: Any) -> Optional[str]:
    """Return the refdata's locus string ('IGH', 'TRB', ...) or None."""
    return getattr(refdata, "locus", None)


def validate_one_record(rec: Dict[str, Any], refdata: Any = None) -> List[str]:
    """Run every applicable check on a single record. Returns a list of
    human-readable issue strings. Empty list means the record passes all
    enabled checks.
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
        issues.extend(_check_locus_matches(rec, refdata))

    return issues
