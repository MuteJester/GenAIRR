"""Summary statistics over a batch of AIRR record dicts.

Used by `simulate_repertoire` (always-on summary on every call) and
`summarize_dataset` (computes the same summary on user-supplied AIRR
records). Sharing the computation here means the two tools return
identical shapes -- the agent can compare simulated vs external data
apples-to-apples.

Histogram dict keys are JSON-stringified ints per the spec
(Section 4 rule 4). The agent does `int(k)` when computing on them.
"""
from __future__ import annotations

from collections import Counter
from statistics import median
from typing import Any, Dict, List


def _quartiles(values: List[int]) -> List[int]:
    """Return [min, q1, median, q3, max] as a 5-element list of ints.

    Empty input returns [0, 0, 0, 0, 0]. We use linear interpolation on
    sorted values so the result matches what a numpy-aware agent would
    expect from numpy.quantile(values, [0, 0.25, 0.5, 0.75, 1.0]).
    """
    if not values:
        return [0, 0, 0, 0, 0]
    sv = sorted(values)
    n = len(sv)

    def quantile(q: float) -> int:
        # Linear interpolation between two adjacent order statistics.
        pos = q * (n - 1)
        lo = int(pos)
        hi = min(lo + 1, n - 1)
        frac = pos - lo
        return int(round(sv[lo] + (sv[hi] - sv[lo]) * frac))

    return [sv[0], quantile(0.25), quantile(0.5), quantile(0.75), sv[-1]]


def _top_counter(values: List[str], top_n: int) -> Dict[str, int]:
    """Count occurrences and return the top-N as a descending-order dict."""
    counts = Counter(v for v in values if v)  # skip empty strings
    return dict(counts.most_common(top_n))


def compute_repertoire_summary(
    records: List[Dict[str, Any]],
    *,
    summary_top_n: int = 10,
) -> Dict[str, Any]:
    """Build the structured summary used by simulate_repertoire and
    summarize_dataset. Returns a dict containing the always-on fields
    plus optional fields (mutation / corruption quartiles) when the
    underlying counter is non-zero on at least one record.
    """
    n = len(records)

    # -- Productive --------------------------------------------------
    productive_count = sum(1 for r in records if r.get("productive") is True)
    productive_rate = productive_count / n if n > 0 else 0.0

    # -- Gene usage --------------------------------------------------
    # v_call / d_call / j_call may be comma-separated tie-sets; we
    # attribute usage to the first allele (the one projected_call_name
    # surfaced as the canonical pick). This matches what an aligner
    # reading the AIRR record would report as "the call".
    def first_call(rec: Dict[str, Any], key: str) -> str:
        value = rec.get(key) or ""
        return value.split(",", 1)[0].strip()

    v_values = [first_call(r, "v_call") for r in records]
    d_values = [first_call(r, "d_call") for r in records]
    j_values = [first_call(r, "j_call") for r in records]

    summary: Dict[str, Any] = {
        "productive_count": productive_count,
        "productive_rate": productive_rate,
        "v_usage_top": _top_counter(v_values, summary_top_n),
        "d_usage_top": _top_counter(d_values, summary_top_n),
        "j_usage_top": _top_counter(j_values, summary_top_n),
        "n_unique_v": len({v for v in v_values if v}),
        "n_unique_d": len({v for v in d_values if v}),
        "n_unique_j": len({v for v in j_values if v}),
    }

    # -- CDR3 length distribution ------------------------------------
    cdr3_lengths = [
        len(r["junction_aa"])
        for r in records
        if r.get("junction_aa") and isinstance(r["junction_aa"], str)
    ]
    cdr3_hist_counter = Counter(cdr3_lengths)
    summary["cdr3_length_histogram"] = {
        str(k): v for k, v in sorted(cdr3_hist_counter.items())
    }
    if cdr3_lengths:
        cdr3_sorted = sorted(cdr3_lengths)
        summary["cdr3_length_mean"] = round(sum(cdr3_lengths) / len(cdr3_lengths), 2)
        summary["cdr3_length_median"] = int(median(cdr3_sorted))
        summary["cdr3_length_p95"] = cdr3_sorted[
            min(int(0.95 * (len(cdr3_sorted) - 1)), len(cdr3_sorted) - 1)
        ]

    # -- Mutation distribution (when any record had mutations) -------
    mutation_counts = [int(r.get("n_mutations") or 0) for r in records]
    if any(c > 0 for c in mutation_counts):
        summary["mutation_count_quartiles"] = _quartiles(mutation_counts)
        # rate is per-record mutations / per-record sequence length, averaged.
        rates = []
        for r in records:
            n_mut = int(r.get("n_mutations") or 0)
            seq_len = int(r.get("sequence_length") or 0)
            if seq_len > 0:
                rates.append(n_mut / seq_len)
        if rates:
            summary["mutation_rate_per_record_mean"] = round(sum(rates) / len(rates), 4)

    # -- Corruption quartiles (each gated on non-zero presence) ------
    for field, key in [
        ("n_pcr_errors", "n_pcr_errors_quartiles"),
        ("n_indels", "n_indels_quartiles"),
        ("n_quality_errors", "n_quality_errors_quartiles"),
    ]:
        values = [int(r.get(field) or 0) for r in records]
        if any(v > 0 for v in values):
            summary[key] = _quartiles(values)

    return summary
