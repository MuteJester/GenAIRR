"""
Helper functions for the GenAIRR MCP server.

Validation, scoring, analysis, and formatting logic used by the MCP
tool definitions in ``mcp_server.py``.  Every function is pure — takes
data in, returns JSON-serialisable data out.
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from ..alleles.allele import VAllele, DAllele, JAllele, CAllele


# ═══════════════════════════════════════════════════════════════════
# Allele lookup helpers
# ═══════════════════════════════════════════════════════════════════

def flatten_alleles(allele_dict) -> Dict[str, Any]:
    """Flatten ``{gene: [Allele, ...]}`` to ``{name: Allele}``."""
    out = {}
    if allele_dict is None:
        return out
    for alleles in allele_dict.values():
        for allele in alleles:
            out[allele.name] = allele
    return out


def find_allele(dc, allele_name: str):
    """Search all segments of a DataConfig for an allele by name."""
    for attr in ("v_alleles", "d_alleles", "j_alleles", "c_alleles"):
        ad = getattr(dc, attr, None)
        if ad is None:
            continue
        for alleles in ad.values():
            for a in alleles:
                if a.name == allele_name:
                    return a
    return None


def _segment_label(allele) -> str:
    if isinstance(allele, VAllele):
        return "V"
    if isinstance(allele, DAllele):
        return "D"
    if isinstance(allele, JAllele):
        return "J"
    if isinstance(allele, CAllele):
        return "C"
    return "?"


def format_allele_info(allele) -> Dict[str, Any]:
    """Return a compact JSON-safe dict for a single Allele."""
    seq = allele.ungapped_seq
    preview = seq[:80] + ("..." if len(seq) > 80 else "")
    info: Dict[str, Any] = {
        "name": allele.name,
        "segment": _segment_label(allele),
        "sequence_preview": preview,
        "ungapped_length": allele.ungapped_len,
        "family": allele.family,
        "gene": allele.gene,
    }
    anchor = getattr(allele, "anchor", None)
    if anchor is not None:
        info["anchor"] = anchor
    return info


# ═══════════════════════════════════════════════════════════════════
# Mutation parsing
# ═══════════════════════════════════════════════════════════════════

def parse_mutations(mutation_str: str) -> Dict[int, Tuple[str, str]]:
    """Parse ``"pos:X>Y,pos:X>Y,..."`` into ``{pos: (from_base, to_base)}``."""
    result: Dict[int, Tuple[str, str]] = {}
    if not mutation_str:
        return result
    for entry in mutation_str.split(","):
        entry = entry.strip()
        if not entry or ":" not in entry:
            continue
        pos_str, change = entry.split(":", 1)
        try:
            pos = int(pos_str)
        except ValueError:
            continue
        if ">" not in change:
            continue
        parts = change.split(">", 1)
        if len(parts) == 2:
            result[pos] = (parts[0], parts[1])
    return result


# ═══════════════════════════════════════════════════════════════════
# Allele scoring
# ═══════════════════════════════════════════════════════════════════

def count_matches(seq_region: str, allele_region: str) -> Tuple[int, int]:
    """Count N-aware base matches. Returns ``(n_matches, n_compared)``."""
    n_match = 0
    n_compared = 0
    for s, a in zip(seq_region, allele_region):
        su, au = s.upper(), a.upper()
        if su == "N" or au == "N":
            continue
        n_compared += 1
        if su == au:
            n_match += 1
    return n_match, n_compared


def score_allele_against_sequence(
    seq: str,
    allele_seq: str,
    seq_start: int,
    seq_end: int,
    germ_start: int,
    germ_end: int,
) -> Optional[Tuple[int, int]]:
    """Score how well an allele matches the sequence region."""
    if germ_start < 0 or germ_end <= germ_start:
        return None
    allele_len = len(allele_seq)
    region_len = min(seq_end - seq_start, germ_end - germ_start)
    if region_len <= 0:
        return None

    n_match = 0
    n_compared = 0
    for i in range(region_len):
        s = seq[seq_start + i].upper() if seq_start + i < len(seq) else "N"
        if germ_start + i < allele_len:
            a = allele_seq[germ_start + i].upper()
        else:
            n_compared += 1
            continue
        if s == "N" or a == "N":
            continue
        n_compared += 1
        if s == a:
            n_match += 1
    return n_match, n_compared


def score_all_alleles(dc, rec: Dict[str, Any], segment: str) -> List[Dict[str, Any]]:
    """Score every allele of *segment* against the record, return top results."""
    seg = segment.lower()
    allele_dict = getattr(dc, f"{seg}_alleles", None)
    if not allele_dict:
        return []

    flat = flatten_alleles(allele_dict)
    seq = rec.get("sequence", "")
    ss = rec.get(f"{seg}_sequence_start", 0)
    se = rec.get(f"{seg}_sequence_end", 0)
    gs = rec.get(f"{seg}_germline_start", 0)
    ge = rec.get(f"{seg}_germline_end", 0)

    reported_calls = set()
    call_str = rec.get(f"{seg}_call", "") or ""
    for c in call_str.split(","):
        c = c.strip()
        if c:
            reported_calls.add(c)

    results = []
    for name, allele in flat.items():
        score = score_allele_against_sequence(
            seq, allele.ungapped_seq, ss, se, gs, ge,
        )
        if score is None:
            continue
        matches, compared = score
        fraction = matches / compared if compared > 0 else 0.0
        results.append({
            "name": name,
            "matches": matches,
            "compared": compared,
            "fraction": round(fraction, 4),
            "is_reported_call": name in reported_calls,
        })

    results.sort(key=lambda r: (-r["fraction"], -r["matches"]))
    return results[:10]


# ═══════════════════════════════════════════════════════════════════
# Record validation
# ═══════════════════════════════════════════════════════════════════

def validate_record(rec: Dict[str, Any]) -> Dict[str, Any]:
    """Run comprehensive validation checks on a single AIRR record."""
    checks: List[Dict[str, Any]] = []
    seq = rec.get("sequence", "")
    seq_len = len(seq)
    germ = rec.get("germline_alignment", "")

    # 1. nucleotide_validity
    bad_chars = set(seq.upper()) - {"A", "T", "C", "G", "N", ""}
    checks.append({
        "name": "nucleotide_validity",
        "passed": len(bad_chars) == 0,
        "detail": f"invalid chars: {bad_chars}" if bad_chars else "ok",
    })

    # 2. coordinate_bounds
    coord_fields = [
        "v_sequence_start", "v_sequence_end",
        "d_sequence_start", "d_sequence_end",
        "j_sequence_start", "j_sequence_end",
    ]
    oob = []
    for f in coord_fields:
        val = rec.get(f, 0)
        if isinstance(val, int) and (val < 0 or val > seq_len):
            oob.append(f"{f}={val}")
    checks.append({
        "name": "coordinate_bounds",
        "passed": len(oob) == 0,
        "detail": f"out of bounds: {', '.join(oob)}" if oob else "ok",
    })

    # 3. segment_ordering
    vs = rec.get("v_sequence_start", 0)
    ve = rec.get("v_sequence_end", 0)
    ds = rec.get("d_sequence_start", 0)
    de = rec.get("d_sequence_end", 0)
    js = rec.get("j_sequence_start", 0)
    je = rec.get("j_sequence_end", 0)
    has_d = ds != de
    has_j = js != je  # J can be absent if 3' corruption removed it
    is_rc = rec.get("is_reverse_complement", False)
    order_ok = vs <= ve
    if not is_rc:
        if has_d and has_j:
            order_ok = order_ok and ve <= ds and ds <= de and de <= js and js <= je
        elif has_d:
            order_ok = order_ok and ve <= ds and ds <= de
        elif has_j:
            order_ok = order_ok and ve <= js and js <= je
    else:
        # In RC read orientation, segment order is mirrored.
        if has_d and has_j:
            order_ok = order_ok and js <= je and je <= ds and ds <= de and de <= vs
        elif has_d:
            order_ok = order_ok and ds <= de and de <= vs
        elif has_j:
            order_ok = order_ok and js <= je and je <= vs
    # If neither D nor J present, only V ordering matters
    checks.append({
        "name": "segment_ordering",
        "passed": order_ok,
        "detail": (
            f"V=[{vs},{ve}] D=[{ds},{de}] J=[{js},{je}] rc={is_rc}"
            if not order_ok else "ok"
        ),
    })

    # 4. junction_bounds
    jstart = rec.get("junction_start", 0)
    jend = rec.get("junction_end", 0)
    jb_ok = 0 <= jstart <= jend <= seq_len if isinstance(jstart, int) and isinstance(jend, int) else True
    checks.append({
        "name": "junction_bounds",
        "passed": jb_ok,
        "detail": f"junction [{jstart},{jend}] seq_len={seq_len}" if not jb_ok else "ok",
    })

    # 5. junction_length_match
    jlen = rec.get("junction_length", 0) or 0
    computed_jlen = jend - jstart if isinstance(jstart, int) and isinstance(jend, int) else 0
    jlm_ok = jlen == computed_jlen
    checks.append({
        "name": "junction_length_match",
        "passed": jlm_ok,
        "detail": f"field={jlen} computed={computed_jlen}" if not jlm_ok else "ok",
    })

    # 6. productive_consistency
    productive = rec.get("productive", False)
    stop_codon = rec.get("stop_codon", False)
    vj_in_frame = rec.get("vj_in_frame", False)
    prod_ok = True
    prod_detail = "ok"
    if productive:
        if stop_codon:
            prod_ok = False
            prod_detail = "productive=True but stop_codon=True"
        if not vj_in_frame:
            prod_ok = False
            prod_detail = "productive=True but vj_in_frame=False"
    checks.append({
        "name": "productive_consistency",
        "passed": prod_ok,
        "detail": prod_detail,
    })

    # 7. mutation_count
    mut_str = rec.get("mutations", "") or ""
    n_mut_field = rec.get("n_mutations", 0) or 0
    parsed = parse_mutations(mut_str)
    mc_ok = n_mut_field == len(parsed)
    checks.append({
        "name": "mutation_count",
        "passed": mc_ok,
        "detail": f"field={n_mut_field} parsed={len(parsed)}" if not mc_ok else "ok",
    })

    # 8. mutation_content
    mut_errors = []
    if germ and seq and parsed:
        for pos, (from_b, to_b) in parsed.items():
            if pos >= len(germ) or pos >= len(seq):
                mut_errors.append(f"pos {pos} out of bounds")
                continue
            g = germ[pos].upper()
            s = seq[pos].upper()
            if g != "N" and g != "-" and from_b.upper() != g:
                mut_errors.append(f"pos {pos}: germ={g} but from={from_b}")
            if s != "N" and to_b.upper() != s:
                mut_errors.append(f"pos {pos}: seq={s} but to={to_b}")
    checks.append({
        "name": "mutation_content",
        "passed": len(mut_errors) == 0,
        "detail": "; ".join(mut_errors[:5]) if mut_errors else "ok",
    })

    # 9. sequence_length_field
    sl_field = rec.get("sequence_length", None)
    sl_ok = sl_field is None or sl_field == seq_len
    checks.append({
        "name": "sequence_length_field",
        "passed": sl_ok,
        "detail": f"field={sl_field} actual={seq_len}" if not sl_ok else "ok",
    })

    n_passed = sum(1 for c in checks if c["passed"])
    n_failed = len(checks) - n_passed
    return {
        "valid": n_failed == 0,
        "n_passed": n_passed,
        "n_failed": n_failed,
        "checks": checks,
    }


# ═══════════════════════════════════════════════════════════════════
# Germline alignment analysis
# ═══════════════════════════════════════════════════════════════════

def align_to_germline(dc, rec: Dict[str, Any]) -> Dict[str, Any]:
    """Compare a record's sequence to its germline reference alleles."""
    seq = rec.get("sequence", "")
    germ = rec.get("germline_alignment", "")
    mut_str = rec.get("mutations", "") or ""
    reported = parse_mutations(mut_str)

    # Build position-level comparison
    details = []
    true_mut = 0
    unreported = 0
    phantom = 0

    compare_len = min(len(seq), len(germ))
    for pos in range(compare_len):
        s = seq[pos].upper()
        g = germ[pos].upper()
        if g in ("-", "N") or s == "N":
            continue
        is_reported = pos in reported
        if s != g:
            if is_reported:
                true_mut += 1
            else:
                unreported += 1
                if len(details) < 50:
                    details.append({
                        "pos": pos, "seq_base": s, "germ_base": g,
                        "reported": False, "type": "unreported",
                    })
        else:
            if is_reported:
                phantom += 1
                if len(details) < 50:
                    details.append({
                        "pos": pos, "seq_base": s, "germ_base": g,
                        "reported": True, "type": "phantom",
                    })

    return {
        "v_allele": rec.get("v_call", ""),
        "d_allele": rec.get("d_call", ""),
        "j_allele": rec.get("j_call", ""),
        "n_true_mutations": true_mut,
        "n_unreported": unreported,
        "n_phantom": phantom,
        "details": details,
    }


# ═══════════════════════════════════════════════════════════════════
# Region classification
# ═══════════════════════════════════════════════════════════════════

def classify_record_positions(dc, rec: Dict[str, Any]) -> Dict[str, Any]:
    """Classify every position in the sequence into IMGT regions."""
    from ..utilities.imgt_regions import compute_v_region_boundaries, classify_position

    seq = rec.get("sequence", "")
    v_call = (rec.get("v_call", "") or "").split(",")[0].strip()
    if not v_call:
        return {"error": "no v_call in record"}

    # Find V allele in config
    v_allele = find_allele(dc, v_call)
    if v_allele is None or not isinstance(v_allele, VAllele):
        return {"error": f"V allele {v_call} not found in config"}

    v_boundaries = compute_v_region_boundaries(v_allele)
    vs = rec.get("v_sequence_start", 0)
    ve = rec.get("v_sequence_end", 0)
    js = rec.get("j_sequence_start", 0)
    je = rec.get("j_sequence_end", 0)
    jstart = rec.get("junction_start", 0)
    jend = rec.get("junction_end", 0)

    # Classify each position
    region_counts: Dict[str, int] = defaultdict(int)
    region_ranges: Dict[str, List[int]] = defaultdict(list)
    for pos in range(len(seq)):
        region = classify_position(pos, vs, ve, v_boundaries, jstart, jend, je)
        region_counts[region] += 1
        region_ranges[region].append(pos)

    # Compute boundaries (min, max+1) for each region
    boundaries = {}
    for region, positions in region_ranges.items():
        if positions:
            boundaries[region] = [min(positions), max(positions) + 1]

    # Mutations by region
    mut_str = rec.get("mutations", "") or ""
    parsed = parse_mutations(mut_str)
    mut_by_region: Dict[str, int] = defaultdict(int)
    for pos in parsed:
        if pos < len(seq):
            region = classify_position(pos, vs, ve, v_boundaries, jstart, jend, je)
            mut_by_region[region] += 1

    # Sequence excerpts by region (truncated)
    seq_by_region = {}
    for region, positions in region_ranges.items():
        if positions:
            start, end = min(positions), max(positions) + 1
            excerpt = seq[start:end]
            seq_by_region[region] = excerpt[:80] + ("..." if len(excerpt) > 80 else "")

    return {
        "region_boundaries": boundaries,
        "mutations_by_region": dict(mut_by_region),
        "sequence_by_region": seq_by_region,
    }


# ═══════════════════════════════════════════════════════════════════
# Mutation pattern analysis
# ═══════════════════════════════════════════════════════════════════

_TRANSITIONS = {
    ("A", "G"), ("G", "A"),  # purine ↔ purine
    ("C", "T"), ("T", "C"),  # pyrimidine ↔ pyrimidine
}

_WRCY = re.compile(r"[AT][AG]C[CT]", re.IGNORECASE)
_RGYW = re.compile(r"[AG]G[CT][AT]", re.IGNORECASE)


def analyze_mutation_patterns(records: List[Dict], dc) -> Dict[str, Any]:
    """Analyze SHM patterns across a batch of records."""
    from ..utilities.imgt_regions import compute_v_region_boundaries, classify_position

    total_mutations = 0
    per_region: Dict[str, int] = defaultdict(int)
    n_transitions = 0
    n_transversions = 0
    wrcy_count = 0
    rgyw_count = 0

    for rec in records:
        seq = rec.get("sequence", "")
        germ = rec.get("germline_alignment", "")
        mut_str = rec.get("mutations", "") or ""
        parsed = parse_mutations(mut_str)
        if not parsed:
            continue

        # Get V allele for region classification
        v_call = (rec.get("v_call", "") or "").split(",")[0].strip()
        v_allele = find_allele(dc, v_call) if v_call else None
        v_boundaries = (
            compute_v_region_boundaries(v_allele)
            if v_allele and isinstance(v_allele, VAllele)
            else {}
        )

        vs = rec.get("v_sequence_start", 0)
        ve = rec.get("v_sequence_end", 0)
        je = rec.get("j_sequence_end", 0)
        jstart = rec.get("junction_start", 0)
        jend = rec.get("junction_end", 0)

        for pos, (from_b, to_b) in parsed.items():
            total_mutations += 1
            fb, tb = from_b.upper(), to_b.upper()

            # Transition / transversion
            if (fb, tb) in _TRANSITIONS:
                n_transitions += 1
            elif fb in "ACGT" and tb in "ACGT":
                n_transversions += 1

            # Region classification
            if v_boundaries and pos < len(seq):
                region = classify_position(pos, vs, ve, v_boundaries, jstart, jend, je)
                per_region[region] += 1

            # Hotspot check on germline context
            if germ and 1 <= pos < len(germ) - 2:
                context = germ[pos - 1:pos + 3]
                if len(context) == 4 and "-" not in context:
                    if _WRCY.fullmatch(context):
                        wrcy_count += 1
                    if _RGYW.fullmatch(context):
                        rgyw_count += 1

    # Compute fractions
    per_region_out = {}
    for region in ("FWR1", "CDR1", "FWR2", "CDR2", "FWR3", "CDR3", "FWR4", "NP"):
        count = per_region.get(region, 0)
        per_region_out[region] = {
            "count": count,
            "fraction": round(count / total_mutations, 4) if total_mutations else 0,
        }

    ti_tv = round(n_transitions / n_transversions, 3) if n_transversions > 0 else None
    motif_total = wrcy_count + rgyw_count

    return {
        "n_records": len(records),
        "total_mutations": total_mutations,
        "per_region": per_region_out,
        "transition_transversion_ratio": ti_tv,
        "hotspot_analysis": {
            "wrcy_count": wrcy_count,
            "rgyw_count": rgyw_count,
            "total_in_motif": motif_total,
            "fraction_in_motif": round(motif_total / total_mutations, 4) if total_mutations else 0,
        },
    }


# ═══════════════════════════════════════════════════════════════════
# Statistical summary
# ═══════════════════════════════════════════════════════════════════

def _stats(values: List[float]) -> Dict[str, Any]:
    """Compute mean/std/min/max for a list of floats."""
    if not values:
        return {"mean": 0, "std": 0, "min": 0, "max": 0}
    n = len(values)
    mean = sum(values) / n
    variance = sum((x - mean) ** 2 for x in values) / n if n > 1 else 0
    return {
        "mean": round(mean, 4),
        "std": round(math.sqrt(variance), 4),
        "min": round(min(values), 4),
        "max": round(max(values), 4),
    }


def _top_n(counter: Dict[str, int], n: int = 20) -> Dict[str, int]:
    """Return the top-N entries from a counter dict."""
    return dict(sorted(counter.items(), key=lambda kv: -kv[1])[:n])


def compute_dataset_summary(records: List[Dict]) -> Dict[str, Any]:
    """Compute aggregate statistics over a list of AIRR records."""
    n = len(records)
    if n == 0:
        return {"n_sequences": 0}

    productive_count = sum(1 for r in records if r.get("productive"))
    mutation_rates = [r["mutation_rate"] for r in records if isinstance(r.get("mutation_rate"), (int, float))]
    junction_lengths = [r["junction_length"] for r in records if isinstance(r.get("junction_length"), int) and r["junction_length"] > 0]
    np1_lengths = [r["np1_length"] for r in records if isinstance(r.get("np1_length"), int)]
    np2_lengths = [r["np2_length"] for r in records if isinstance(r.get("np2_length"), int)]

    v_usage: Dict[str, int] = defaultdict(int)
    d_usage: Dict[str, int] = defaultdict(int)
    j_usage: Dict[str, int] = defaultdict(int)
    for r in records:
        vc = (r.get("v_call", "") or "").split(",")[0].strip()
        dc_call = (r.get("d_call", "") or "").split(",")[0].strip()
        jc = (r.get("j_call", "") or "").split(",")[0].strip()
        if vc:
            v_usage[vc] += 1
        if dc_call:
            d_usage[dc_call] += 1
        if jc:
            j_usage[jc] += 1

    return {
        "n_sequences": n,
        "productive_rate": round(productive_count / n, 4),
        "mutation_rate": _stats(mutation_rates),
        "junction_length": _stats([float(x) for x in junction_lengths]),
        "v_gene_usage": _top_n(v_usage),
        "d_gene_usage": _top_n(d_usage),
        "j_gene_usage": _top_n(j_usage),
        "np1_length": _stats([float(x) for x in np1_lengths]),
        "np2_length": _stats([float(x) for x in np2_lengths]),
    }


# ═══════════════════════════════════════════════════════════════════
# Config section extraction
# ═══════════════════════════════════════════════════════════════════

def extract_config_section(dc, section: str) -> Dict[str, Any]:
    """Extract a specific section of config internals for inspection."""
    if section == "gene_use":
        raw = getattr(dc, "gene_use_dict", {}) or {}
        out = {}
        for seg, usage in raw.items():
            if isinstance(usage, dict):
                out[seg] = {k: round(v, 4) if isinstance(v, float) else v for k, v in usage.items()}
            else:
                out[seg] = usage
        return {"section": "gene_use", "data": out}

    if section == "trimming":
        raw = getattr(dc, "trim_dicts", {}) or {}
        out = {}
        for trim_key, families in raw.items():
            if isinstance(families, dict):
                summary = {}
                for family, genes in families.items():
                    if isinstance(genes, dict):
                        # Truncate to first 10 entries per gene
                        summary[family] = {
                            k: {str(kk): round(vv, 4) if isinstance(vv, float) else vv
                                for kk, vv in list(v.items())[:15]}
                            if isinstance(v, dict) else v
                            for k, v in list(genes.items())[:10]
                        }
                    else:
                        summary[family] = genes
                out[trim_key] = summary
            else:
                out[trim_key] = families
        return {"section": "trimming", "data": out}

    if section == "np_params":
        return {
            "section": "np_params",
            "data": {
                "NP_first_bases": getattr(dc, "NP_first_bases", {}),
                "NP_lengths": getattr(dc, "NP_lengths", {}),
                "NP_transitions_keys": list((getattr(dc, "NP_transitions", {}) or {}).keys()),
            },
        }

    if section == "p_nucleotides":
        return {
            "section": "p_nucleotides",
            "data": getattr(dc, "p_nucleotide_length_probs", {}),
        }

    return {"error": f"Unknown section: {section!r}. Use: gene_use, trimming, np_params, p_nucleotides"}


# ═══════════════════════════════════════════════════════════════════
# ASeq introspection helpers (for pipeline hook snapshots)
# ═══════════════════════════════════════════════════════════════════

def summarize_aseq(nodes: List[Dict]) -> Dict[str, Any]:
    """Compact summary of an ASeq node list from a snapshot."""
    if not nodes:
        return {"total_nodes": 0}

    seg_counts: Dict[str, int] = defaultdict(int)
    flag_counts: Dict[str, int] = defaultdict(int)
    anchor_positions: List[int] = []

    for n in nodes:
        seg_counts[n.get("seg", "?")] += 1
        for f in n.get("flags", []):
            flag_counts[f] += 1
        if "anchor" in n.get("flags", []):
            anchor_positions.append(n.get("pos", -1))

    return {
        "total_nodes": len(nodes),
        "segments": dict(seg_counts),
        "flag_counts": dict(flag_counts),
        "anchor_positions": anchor_positions,
    }


def detect_node_anomalies(nodes: List[Dict]) -> List[Dict[str, Any]]:
    """Flag suspicious node states in an ASeq snapshot."""
    anomalies = []

    for n in nodes:
        pos = n.get("pos", -1)
        cur = n.get("cur", "?")
        germ = n.get("germ", ".")
        seg = n.get("seg", "?")
        flags = n.get("flags", [])

        # Mutated but cur == germ
        if "mutated" in flags and cur == germ and germ != ".":
            anomalies.append({
                "pos": pos, "type": "phantom_mutation",
                "detail": f"flagged mutated but cur={cur}==germ={germ}",
            })

        # Indel insertion with mutated flag
        if "indel_ins" in flags and "mutated" in flags:
            anomalies.append({
                "pos": pos, "type": "mutated_indel",
                "detail": f"indel-inserted node also has mutated flag",
            })

        # Germline segment node with null germline (indel-inserted nodes are exempt)
        if seg in ("V", "D", "J", "C") and germ == "." and "indel_ins" not in flags:
            anomalies.append({
                "pos": pos, "type": "missing_germline",
                "detail": f"segment={seg} but germline is null",
            })

        # NP node with germline set (should be null)
        if seg in ("NP1", "NP2") and germ not in (".", ""):
            if "indel_ins" not in flags:
                anomalies.append({
                    "pos": pos, "type": "np_has_germline",
                    "detail": f"NP segment has germline={germ}",
                })

        # Non-ACGT current base in non-N node
        if cur.upper() not in ("A", "C", "G", "T", "N", "?"):
            anomalies.append({
                "pos": pos, "type": "invalid_base",
                "detail": f"current base is '{cur}'",
            })

    return anomalies


def validate_codon_rail_snapshot(codons: List[Dict]) -> Dict[str, Any]:
    """Validate codon rail from a snapshot."""
    if not codons:
        return {"valid": True, "n_codons": 0, "n_stop_codons": 0, "issues": []}

    issues = []
    n_stop = 0

    for c in codons:
        aa = c.get("aa", ".")
        is_stop = c.get("is_stop", False)
        bases = c.get("bases", "???")

        if is_stop:
            n_stop += 1

        # Verify amino acid matches bases (basic check)
        if len(bases) == 3 and "?" not in bases:
            from ..utilities.misc import translate
            translated = translate(bases)
            expected_aa = translated[0] if translated else None
            if expected_aa and aa != "." and aa != expected_aa:
                issues.append({
                    "codon_idx": c.get("idx", -1),
                    "bases": bases,
                    "expected_aa": expected_aa,
                    "actual_aa": aa,
                })

    return {
        "valid": len(issues) == 0,
        "n_codons": len(codons),
        "n_stop_codons": n_stop,
        "issues": issues[:20],
    }


def build_germline_diff(nodes: List[Dict]) -> List[Dict[str, Any]]:
    """Position-by-position diff where cur != germ."""
    diffs = []
    for n in nodes:
        cur = n.get("cur", "?")
        germ = n.get("germ", ".")
        if germ != "." and cur != germ:
            diffs.append({
                "pos": n.get("pos", -1),
                "cur": cur,
                "germ": germ,
                "seg": n.get("seg", "?"),
                "flags": n.get("flags", []),
            })
    return diffs


def diff_snapshots(
    snap_a_nodes: List[Dict],
    snap_b_nodes: List[Dict],
) -> Dict[str, Any]:
    """Compare two snapshots: what changed between them."""
    len_a = len(snap_a_nodes)
    len_b = len(snap_b_nodes)

    # Build position maps
    a_by_pos = {n["pos"]: n for n in snap_a_nodes}
    b_by_pos = {n["pos"]: n for n in snap_b_nodes}

    modified = []
    added_positions = sorted(set(b_by_pos) - set(a_by_pos))
    removed_positions = sorted(set(a_by_pos) - set(b_by_pos))

    # Find modifications at shared positions
    shared = sorted(set(a_by_pos) & set(b_by_pos))
    for pos in shared:
        a = a_by_pos[pos]
        b = b_by_pos[pos]
        changes = {}
        if a.get("cur") != b.get("cur"):
            changes["cur"] = {"from": a.get("cur"), "to": b.get("cur")}
        if a.get("germ") != b.get("germ"):
            changes["germ"] = {"from": a.get("germ"), "to": b.get("germ")}
        if set(a.get("flags", [])) != set(b.get("flags", [])):
            changes["flags"] = {
                "added": sorted(set(b.get("flags", [])) - set(a.get("flags", []))),
                "removed": sorted(set(a.get("flags", [])) - set(b.get("flags", []))),
            }
        if a.get("seg") != b.get("seg"):
            changes["seg"] = {"from": a.get("seg"), "to": b.get("seg")}
        if changes:
            modified.append({"pos": pos, "changes": changes})

    # Flag stats
    a_flags = defaultdict(int)
    b_flags = defaultdict(int)
    for n in snap_a_nodes:
        for f in n.get("flags", []):
            a_flags[f] += 1
    for n in snap_b_nodes:
        for f in n.get("flags", []):
            b_flags[f] += 1

    all_flag_names = sorted(set(a_flags) | set(b_flags))
    flag_changes = {}
    for f in all_flag_names:
        ac = a_flags.get(f, 0)
        bc = b_flags.get(f, 0)
        if ac != bc:
            flag_changes[f] = {"before": ac, "after": bc, "delta": bc - ac}

    return {
        "before_length": len_a,
        "after_length": len_b,
        "n_added": len(added_positions),
        "n_removed": len(removed_positions),
        "n_modified": len(modified),
        "added_positions": added_positions[:50],
        "removed_positions": removed_positions[:50],
        "modified": modified[:50],
        "flag_changes": flag_changes,
    }


def format_snapshot_timeline(snapshots: List[Dict]) -> str:
    """Human-readable timeline of how the sequence evolved through stages."""
    lines = []
    for snap in snapshots:
        hook = snap.get("hook", "?")
        nodes = snap.get("nodes", [])
        summary = summarize_aseq(nodes)

        seg_str = " ".join(
            f"{seg}:{count}"
            for seg, count in sorted(summary.get("segments", {}).items())
        )
        flags = summary.get("flag_counts", {})

        parts = [f"{hook:25s}  {summary['total_nodes']:4d} nodes"]
        if seg_str:
            parts.append(f"({seg_str})")
        if flags.get("mutated", 0):
            parts.append(f"+{flags['mutated']} mutations")
        if flags.get("indel_ins", 0):
            parts.append(f"+{flags['indel_ins']} indels")
        if flags.get("is_n", 0):
            parts.append(f"+{flags['is_n']} Ns")

        lines.append("  ".join(parts))

    return "\n".join(lines)


def aggregate_flag_stats(
    all_nodes_lists: List[List[Dict]],
) -> Dict[str, Any]:
    """Batch flag distribution across N sequences."""
    total_nodes = 0
    per_flag: Dict[str, int] = defaultdict(int)
    per_seg_flag: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    n_anomalies = 0

    for nodes in all_nodes_lists:
        total_nodes += len(nodes)
        anomalies = detect_node_anomalies(nodes)
        n_anomalies += len(anomalies)
        for n in nodes:
            seg = n.get("seg", "?")
            for f in n.get("flags", []):
                per_flag[f] += 1
                per_seg_flag[seg][f] += 1

    return {
        "n_sequences": len(all_nodes_lists),
        "total_nodes": total_nodes,
        "per_flag_counts": dict(per_flag),
        "per_segment_flag_counts": {
            seg: dict(flags) for seg, flags in per_seg_flag.items()
        },
        "total_anomalies": n_anomalies,
    }
