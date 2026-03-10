"""
visualize.py — Generate standalone HTML files showing an "exploding view"
dissection of a simulated AIRR sequence record.

Usage:
    from GenAIRR.utilities.visualize import visualize_sequence
    result = experiment.run(n=1)
    visualize_sequence(result[0], "sequence.html")
"""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any, Dict, Optional, Union

# ── Segment colors (matching the website DissectionBay) ──────

SEGMENT_COLORS = {
    "V": "#4A90D9",
    "NP1": "#5DC48C",
    "D": "#E8685A",
    "NP2": "#5DC48C",
    "J": "#E8A838",
}

RED = "#DC2626"
PURPLE = "#9B6FC4"


# ── Field extraction helpers ─────────────────────────────────

def _int(rec: dict, key: str, default: int = 0) -> int:
    v = rec.get(key)
    if v is None:
        return default
    try:
        return int(v)
    except (ValueError, TypeError):
        return default


def _str(rec: dict, key: str, default: str = "") -> str:
    v = rec.get(key)
    return str(v) if v is not None else default


def _bool(rec: dict, key: str) -> bool:
    v = rec.get(key)
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        return v.lower() in ("true", "1", "yes", "t")
    return bool(v) if v is not None else False


def _float(rec: dict, key: str, default: float = 0.0) -> float:
    v = rec.get(key)
    if v is None:
        return default
    try:
        return float(v)
    except (ValueError, TypeError):
        return default


def _parse_mutations(raw: Any) -> Dict[int, str]:
    """Parse the mutations field (string like '42:A>G;100:T>C' or dict)."""
    if not raw:
        return {}
    if isinstance(raw, dict):
        return {int(k): str(v) for k, v in raw.items()}
    if isinstance(raw, str):
        # Try JSON first
        try:
            parsed = json.loads(raw)
            if isinstance(parsed, dict):
                return {int(k): str(v) for k, v in parsed.items()}
        except (json.JSONDecodeError, ValueError):
            pass
        # Try comma or semicolon-separated "pos:from>to" format
        out = {}
        # Split on comma or semicolon
        sep = "," if "," in raw else ";"
        for part in raw.split(sep):
            part = part.strip()
            if not part:
                continue
            if ":" in part:
                pos_str, desc = part.split(":", 1)
                try:
                    out[int(pos_str)] = desc
                except ValueError:
                    pass
        return out
    return {}


# ── Build the HTML ────────────────────────────────────────────

def _esc(s: str) -> str:
    return html.escape(s, quote=True)


def _build_segment_bar_html(segments: list, seq_len: int, mutations: Dict[int, str]) -> str:
    """Build the color-coded assembled sequence bar."""
    parts = []
    for seg in segments:
        w = max(((seg["end"] - seg["start"]) / seq_len) * 100, 1.5)
        label = seg["label"] if w > 4 else ""
        start_label = seg["start"] + 1
        end_label = seg["end"]
        parts.append(
            f'<div class="seg-bar-part" style="width:{w:.2f}%;background:{seg["color"]}" '
            f'title="{_esc(seg["full_label"])}: {start_label}..{end_label}">'
            f'<span class="seg-bar-label">{_esc(label)}</span>'
            f'<span class="seg-bar-pos left">{start_label}</span>'
            + (f'<span class="seg-bar-pos right">{end_label}</span>' if w > 8 else "")
            + "</div>"
        )

    # Mutation dots
    dots = []
    for pos in mutations:
        left_pct = (pos / seq_len) * 100
        dots.append(f'<div class="mut-dot" style="left:{left_pct:.3f}%"></div>')

    return (
        '<div class="assembled-bar">'
        + "".join(parts)
        + "".join(dots)
        + "</div>"
    )


def _build_junction_bracket(junction_start: int, junction_end: int, junction_aa: str, seq_len: int) -> str:
    if seq_len == 0 or junction_end <= junction_start:
        return ""
    left_pct = (junction_start / seq_len) * 100
    width_pct = ((junction_end - junction_start) / seq_len) * 100
    return (
        '<div class="junction-bracket">'
        f'<div class="junction-inner" style="left:{left_pct:.2f}%;width:{width_pct:.2f}%">'
        '<svg style="width:100%;height:12px" preserveAspectRatio="none">'
        f'<line x1="0" y1="0" x2="0" y2="10" stroke="{PURPLE}" stroke-width="1.5"/>'
        f'<line x1="0" y1="10" x2="100%" y2="10" stroke="{PURPLE}" stroke-width="1.5"/>'
        f'<line x1="100%" y1="0" x2="100%" y2="10" stroke="{PURPLE}" stroke-width="1.5"/>'
        "</svg>"
        f'<span class="junction-aa">CDR3: {_esc(junction_aa)}</span>'
        "</div>"
        "</div>"
    )


def _sparkline_svg(start: int, end: int, mutations: Dict[int, str], height: int = 24) -> str:
    length = end - start
    if length <= 0:
        return ""
    bins = min(length, 60)
    bin_size = max(1, length // bins)
    counts = []
    for i in range(bins):
        b_start = start + i * bin_size
        b_end = b_start + bin_size
        c = sum(1 for p in mutations if b_start <= p < b_end)
        counts.append(c)
    mx = max(max(counts), 1)
    rects = []
    for i, c in enumerate(counts):
        fill = RED if c > 0 else "#555"
        opacity = 0.85 if c > 0 else 0.2
        h = c if c > 0 else 0.15
        rects.append(
            f'<rect x="{i}" y="{mx - h}" width="0.85" height="{h}" '
            f'fill="{fill}" opacity="{opacity}"/>'
        )
    return (
        f'<svg width="100%" height="{height}" viewBox="0 0 {bins} {mx}" '
        f'preserveAspectRatio="none" style="display:block;border-radius:2px;overflow:hidden">'
        + "".join(rects)
        + "</svg>"
    )


def _pnp_bar(prefix: str, n_region: str, suffix: str, color: str) -> str:
    total = len(prefix) + len(n_region) + len(suffix)
    if total == 0:
        return ""
    parts_data = [
        ("P", len(prefix), 0.6),
        ("N", len(n_region), 0.25),
        ("P", len(suffix), 0.6),
    ]
    parts = []
    for label, length, opacity in parts_data:
        if length > 0:
            w = max((length / total) * 100, 3)
            parts.append(
                f'<div class="pnp-part" style="width:{w:.1f}%;background:{color};opacity:{opacity}">'
                f'<span class="pnp-label">{label}</span></div>'
            )
    return '<div class="pnp-bar">' + "".join(parts) + "</div>"


def _nt_strip(sequence: str, start_pos: int, color: str, mutations: Dict[int, str]) -> str:
    """Nucleotide detail with mutation highlighting."""
    rows = []
    for row_idx in range(0, len(sequence), 80):
        chunk = sequence[row_idx : row_idx + 80]
        chars = []
        for ci, ch in enumerate(chunk):
            abs_pos = start_pos + row_idx + ci
            is_mut = abs_pos in mutations
            cls = "nt-char mut" if is_mut else "nt-char"
            dot = '<span class="mut-marker"></span>' if is_mut else ""
            chars.append(f'<span class="{cls}">{dot}{_esc(ch)}</span>')
        line_start = start_pos + row_idx + 1
        line_end = min(start_pos + row_idx + 80, start_pos + len(sequence))
        rows.append(
            f'<div class="nt-row">'
            f'<span class="nt-linenum">{line_start}</span>'
            f'<div class="nt-chars">{"".join(chars)}</div>'
            f'<span class="nt-linenum end">{line_end}</span>'
            f"</div>"
        )
    bg = f"color-mix(in srgb, {color} 5%, #1a1a2e)"
    return f'<div class="nt-strip" style="background:{bg}">{"".join(rows)}</div>'


def _metric(label: str, value: str, color: Optional[str] = None) -> str:
    style = f' style="color:{color}"' if color else ""
    return (
        f'<div class="metric">'
        f'<span class="metric-label">{_esc(label)}</span>'
        f'<span class="metric-value"{style}>{_esc(value)}</span>'
        f"</div>"
    )


def _trim_indicator(label: str, bases: int) -> str:
    if not bases:
        return ""
    return (
        f'<div class="trim-info">'
        f'<span class="trim-icon">&#9986;</span> '
        f'{_esc(label)}: <strong>{bases} bp</strong> trimmed'
        f"</div>"
    )


def _segment_panel(seg: dict, rec: dict, mutations: Dict[int, str], sequence: str) -> str:
    """Build a single exploded segment panel."""
    seg_start = seg["start"]
    seg_end = seg["end"]
    seg_len = seg_end - seg_start
    sub_seq = sequence[seg_start:seg_end]
    seg_muts = {p: v for p, v in mutations.items() if seg_start <= p < seg_end}
    mut_count = len(seg_muts)
    seg_id = seg["id"]
    color = seg["color"]

    header = (
        f'<div class="seg-panel-header" style="border-top:3px solid {color}">'
        f'<span class="seg-chip" style="background:{color}">{_esc(seg["label"])}</span>'
        f'<span class="seg-gene">{_esc(seg["full_label"])}</span>'
        f'<span class="seg-pos">{seg_start + 1}..{seg_end} ({seg_len} nt)</span>'
        f"</div>"
    )

    body_parts = []

    if seg_id == "V":
        body_parts.append(
            f'<div class="metric-grid g3">'
            + _metric("LENGTH", f"{seg_len} bp")
            + _metric("MUTATIONS", str(mut_count), RED if mut_count > 0 else None)
            + _metric("3' TRIM", f"{_int(rec, 'v_trim_3')} bp")
            + "</div>"
        )
        if seg_len > 0:
            body_parts.append(
                '<div class="spark-section">'
                '<span class="micro-label">Mutation Density</span>'
                + _sparkline_svg(seg_start, seg_end, mutations)
                + "</div>"
            )
        body_parts.append(_trim_indicator("V-gene 3'", _int(rec, "v_trim_3")))

    elif seg_id == "D":
        d_inv = _bool(rec, "d_inverted")
        body_parts.append(
            f'<div class="metric-grid g4">'
            + _metric("LENGTH", f"{seg_len} bp")
            + _metric("5' TRIM", f"{_int(rec, 'd_trim_5')} bp")
            + _metric("3' TRIM", f"{_int(rec, 'd_trim_3')} bp")
            + _metric("INVERTED", "YES" if d_inv else "NO", RED if d_inv else None)
            + "</div>"
        )
        if d_inv:
            body_parts.append(
                '<div class="inv-badge">'
                '<span class="inv-icon">&#8645;</span> '
                "D-gene inverted (reverse complement used)</div>"
            )
        body_parts.append(_trim_indicator("5'", _int(rec, "d_trim_5")))
        body_parts.append(_trim_indicator("3'", _int(rec, "d_trim_3")))

    elif seg_id == "J":
        body_parts.append(
            f'<div class="metric-grid g3">'
            + _metric("LENGTH", f"{seg_len} bp")
            + _metric("MUTATIONS", str(mut_count), RED if mut_count > 0 else None)
            + _metric("5' TRIM", f"{_int(rec, 'j_trim_5')} bp")
            + "</div>"
        )
        if seg_len > 0:
            body_parts.append(
                '<div class="spark-section">'
                '<span class="micro-label">Mutation Density</span>'
                + _sparkline_svg(seg_start, seg_end, mutations)
                + "</div>"
            )
        body_parts.append(_trim_indicator("J-gene 5'", _int(rec, "j_trim_5")))

    elif seg_id == "NP1":
        p_pre = _str(rec, "np1_p_prefix")
        n_reg = _str(rec, "np1_n_region")
        p_suf = _str(rec, "np1_p_suffix")
        body_parts.append(
            f'<div class="metric-grid g3">'
            + _metric("P-PREFIX", f"{len(p_pre)} bp")
            + _metric("N-ADDITION", f"{len(n_reg)} bp")
            + _metric("P-SUFFIX", f"{len(p_suf)} bp")
            + "</div>"
        )
        body_parts.append(
            '<span class="micro-label">P | N | P Composition</span>'
            + _pnp_bar(p_pre, n_reg, p_suf, color)
        )
        body_parts.append(
            '<div class="np-seq">'
            + f'<span style="color:{color};opacity:0.7">{_esc(p_pre)}</span>'
            + f'<span>{_esc(n_reg)}</span>'
            + f'<span style="color:{color};opacity:0.7">{_esc(p_suf)}</span>'
            + "</div>"
        )

    elif seg_id == "NP2":
        p_pre = _str(rec, "np2_p_prefix")
        n_reg = _str(rec, "np2_n_region")
        p_suf = _str(rec, "np2_p_suffix")
        body_parts.append(
            f'<div class="metric-grid g3">'
            + _metric("P-PREFIX", f"{len(p_pre)} bp")
            + _metric("N-ADDITION", f"{len(n_reg)} bp")
            + _metric("P-SUFFIX", f"{len(p_suf)} bp")
            + "</div>"
        )
        body_parts.append(
            '<span class="micro-label">P | N | P Composition</span>'
            + _pnp_bar(p_pre, n_reg, p_suf, color)
        )
        body_parts.append(
            '<div class="np-seq">'
            + f'<span style="color:{color};opacity:0.7">{_esc(p_pre)}</span>'
            + f'<span>{_esc(n_reg)}</span>'
            + f'<span style="color:{color};opacity:0.7">{_esc(p_suf)}</span>'
            + "</div>"
        )

    # Nucleotide detail (always shown in static view)
    if seg_len > 0 and seg_id in ("V", "D", "J"):
        mut_note = f' <span style="color:{RED}">({mut_count} mutations highlighted)</span>' if mut_count > 0 else ""
        body_parts.append(
            f'<div class="nt-section">'
            f'<span class="micro-label">Nucleotide Sequence{mut_note}</span>'
            + _nt_strip(sub_seq, seg_start, color, seg_muts)
            + "</div>"
        )

    return (
        f'<div class="seg-panel">'
        + header
        + '<div class="seg-panel-body">'
        + "".join(body_parts)
        + "</div></div>"
    )


# ── CSS ──────────────────────────────────────────────────────

CSS = """
:root {
  --bg: #0f1019;
  --bg2: #181825;
  --bg3: #1e1e32;
  --fg: #e0e0f0;
  --fg2: #8888aa;
  --border: #2a2a40;
  --accent: #6366f1;
  --red: #DC2626;
  --purple: #9B6FC4;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
  background: var(--bg);
  color: var(--fg);
  line-height: 1.5;
  padding: 2rem;
}
.container { max-width: 1200px; margin: 0 auto; }
h1 { font-size: 1.4rem; font-weight: 600; margin-bottom: 0.25rem; }
h2 { font-size: 1.1rem; font-weight: 600; margin-bottom: 0.75rem; color: var(--fg2); }
.subtitle { color: var(--fg2); font-size: 0.85rem; margin-bottom: 1.5rem; }

/* Top bar */
.top-bar {
  display: flex; align-items: center; gap: 0.75rem;
  padding: 1rem 1.25rem;
  background: var(--bg2); border: 1px solid var(--border);
  border-radius: 10px; margin-bottom: 1.25rem;
}
.top-bar .dna-icon { color: var(--accent); font-size: 1.2rem; }
.badge {
  padding: 0.2rem 0.6rem; border-radius: 999px;
  font-size: 0.7rem; font-weight: 600; text-transform: uppercase;
  letter-spacing: 0.05em;
}
.badge-pos { background: rgba(34,197,94,0.15); color: #22c55e; }
.badge-neg { background: rgba(220,38,38,0.15); color: #DC2626; }
.seq-id { color: var(--fg2); font-size: 0.8rem; margin-left: auto; font-family: monospace; }

/* Summary grid */
.summary-grid {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
  gap: 0.5rem; margin-bottom: 1.25rem;
}
.summary-cell {
  background: var(--bg2); border: 1px solid var(--border);
  border-radius: 8px; padding: 0.6rem 0.8rem;
}
.summary-cell .metric-label {
  display: block; font-size: 0.6rem; text-transform: uppercase;
  letter-spacing: 0.08em; color: var(--fg2); margin-bottom: 0.15rem;
}
.summary-cell .metric-value {
  font-size: 0.9rem; font-weight: 600;
}
.summary-cell .metric-value.mono { font-family: monospace; font-size: 0.75rem; word-break: break-all; }

/* Corruption */
.corruption-row {
  display: flex; gap: 0.75rem; margin-bottom: 1rem;
}
.corruption-badge {
  display: inline-flex; align-items: center; gap: 0.4rem;
  background: rgba(220,38,38,0.08); border: 1px solid rgba(220,38,38,0.25);
  border-radius: 6px; padding: 0.35rem 0.75rem;
  font-size: 0.75rem; color: var(--red);
}

/* Assembled bar */
.section-label { font-size: 0.7rem; text-transform: uppercase; letter-spacing: 0.08em; color: var(--fg2); margin-bottom: 0.4rem; }
.assembled-bar {
  display: flex; height: 38px; border-radius: 6px; overflow: hidden;
  position: relative; margin-bottom: 0;
  border: 1px solid var(--border);
}
.seg-bar-part {
  position: relative; display: flex; align-items: center; justify-content: center;
  cursor: default; transition: opacity 0.2s;
  min-width: 0;
}
.seg-bar-label {
  font-size: 0.7rem; font-weight: 700; color: #fff; text-shadow: 0 1px 2px rgba(0,0,0,0.4);
  pointer-events: none;
}
.seg-bar-pos {
  position: absolute; bottom: 2px; font-size: 0.55rem; color: rgba(255,255,255,0.7);
  pointer-events: none;
}
.seg-bar-pos.left { left: 3px; }
.seg-bar-pos.right { right: 3px; }
.mut-dot {
  position: absolute; top: -2px; width: 4px; height: 4px;
  background: var(--red); border-radius: 50%;
  transform: translateX(-50%);
  pointer-events: none;
}

/* Junction bracket */
.junction-bracket { position: relative; height: 28px; margin-bottom: 1.5rem; }
.junction-inner { position: absolute; text-align: center; }
.junction-aa {
  display: block; font-size: 0.65rem; color: var(--purple);
  font-family: monospace; margin-top: 1px; white-space: nowrap;
  overflow: hidden; text-overflow: ellipsis;
}

/* Segment panels */
.segments-grid {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
  gap: 0.75rem; margin-top: 1rem;
}
.seg-panel {
  background: var(--bg2); border: 1px solid var(--border);
  border-radius: 8px; overflow: hidden;
}
.seg-panel-header {
  display: flex; align-items: center; gap: 0.5rem;
  padding: 0.6rem 0.8rem; border-bottom: 1px solid var(--border);
}
.seg-chip {
  display: inline-block; padding: 0.15rem 0.5rem; border-radius: 4px;
  font-size: 0.65rem; font-weight: 700; color: #fff;
}
.seg-gene { font-size: 0.75rem; font-weight: 600; }
.seg-pos { font-size: 0.7rem; color: var(--fg2); margin-left: auto; }
.seg-panel-body { padding: 0.75rem 0.8rem; display: flex; flex-direction: column; gap: 0.6rem; }

/* Metrics inside panels */
.metric-grid { display: grid; gap: 0.4rem; }
.metric-grid.g3 { grid-template-columns: repeat(3, 1fr); }
.metric-grid.g4 { grid-template-columns: repeat(4, 1fr); }
.metric {
  background: var(--bg3); border-radius: 5px; padding: 0.35rem 0.5rem;
  text-align: center;
}
.metric-label {
  display: block; font-size: 0.55rem; text-transform: uppercase;
  letter-spacing: 0.06em; color: var(--fg2);
}
.metric-value { font-size: 0.8rem; font-weight: 600; }

/* Trim */
.trim-info {
  font-size: 0.7rem; color: var(--fg2);
  display: flex; align-items: center; gap: 0.3rem;
}
.trim-icon { font-size: 0.85rem; }

/* Inversion */
.inv-badge {
  display: flex; align-items: center; gap: 0.4rem;
  background: rgba(220,38,38,0.08); border: 1px solid rgba(220,38,38,0.2);
  border-radius: 5px; padding: 0.3rem 0.6rem;
  font-size: 0.7rem; color: var(--red);
}
.inv-icon { font-size: 1rem; }

/* PNP bar */
.pnp-bar {
  display: flex; height: 20px; border-radius: 4px; overflow: hidden;
  border: 1px solid var(--border);
}
.pnp-part {
  display: flex; align-items: center; justify-content: center;
  min-width: 0;
}
.pnp-label { font-size: 0.6rem; font-weight: 700; color: #fff; }
.np-seq {
  font-family: monospace; font-size: 0.7rem; word-break: break-all;
  padding: 0.3rem; background: var(--bg3); border-radius: 4px;
}

/* Nucleotide strip */
.nt-section { margin-top: 0.25rem; }
.nt-strip {
  border-radius: 6px; padding: 0.5rem; overflow-x: auto;
  font-family: 'Fira Code', 'JetBrains Mono', 'Cascadia Code', monospace;
  font-size: 0.65rem; line-height: 1.6;
}
.nt-row { display: flex; align-items: center; gap: 0.5rem; }
.nt-linenum { color: var(--fg2); min-width: 3ch; text-align: right; font-size: 0.6rem; user-select: none; }
.nt-linenum.end { text-align: left; }
.nt-chars { display: flex; flex-wrap: wrap; }
.nt-char {
  position: relative; display: inline-block; width: 0.8em; text-align: center;
}
.nt-char.mut { color: var(--red); font-weight: 700; }
.mut-marker {
  position: absolute; top: -3px; left: 50%; transform: translateX(-50%);
  width: 3px; height: 3px; border-radius: 50%; background: var(--red);
}
.micro-label {
  display: block; font-size: 0.6rem; text-transform: uppercase;
  letter-spacing: 0.06em; color: var(--fg2); margin-bottom: 0.3rem;
}
.spark-section { /* wrapper */ }

/* Germline alignment */
.germline-section { margin-top: 1.25rem; }
.germline-row {
  display: flex; font-family: monospace; font-size: 0.65rem; line-height: 1.6;
  gap: 0.5rem; align-items: center;
}
.germline-label { min-width: 5ch; color: var(--fg2); font-size: 0.6rem; text-align: right; }
.germline-chars { display: flex; flex-wrap: wrap; }
.germline-match { color: var(--fg2); opacity: 0.4; }
.germline-mismatch { color: var(--red); font-weight: 700; }
.germline-gap { color: var(--fg2); opacity: 0.2; }

/* Footer */
.footer {
  margin-top: 2rem; padding-top: 1rem; border-top: 1px solid var(--border);
  font-size: 0.65rem; color: var(--fg2); text-align: center;
}
"""


# ── Main function ────────────────────────────────────────────

def visualize_sequence(
    record: Dict[str, Any],
    path: Union[str, Path],
    title: Optional[str] = None,
) -> Path:
    """
    Generate a standalone HTML file with an "exploding view" dissection
    of a simulated AIRR sequence record.

    Parameters
    ----------
    record : dict
        A single AIRR record dict (one element from SimulationResult).
    path : str or Path
        Output file path for the HTML file.
    title : str, optional
        Custom title for the page. Defaults to "Sequence Dissection".

    Returns
    -------
    Path
        The path to the generated file.
    """
    path = Path(path)
    rec = record

    # Extract fields
    sequence = _str(rec, "sequence")
    seq_len = len(sequence)
    germline = _str(rec, "germline_alignment")
    v_call = _str(rec, "v_call")
    d_call = _str(rec, "d_call")
    j_call = _str(rec, "j_call")
    junction_aa = _str(rec, "junction_aa")
    productive = _bool(rec, "productive")
    mutation_rate = _float(rec, "mutation_rate")
    mutations = _parse_mutations(rec.get("mutations"))
    total_mutations = len(mutations)
    cdr3_length = len(junction_aa)

    v_seq_start = _int(rec, "v_sequence_start")
    v_seq_end = _int(rec, "v_sequence_end")
    d_seq_start = _int(rec, "d_sequence_start")
    d_seq_end = _int(rec, "d_sequence_end")
    j_seq_start = _int(rec, "j_sequence_start")
    j_seq_end = _int(rec, "j_sequence_end")
    junction_start = _int(rec, "junction_start", _int(rec, "junction_sequence_start"))
    junction_end = _int(rec, "junction_end", _int(rec, "junction_sequence_end"))

    corruption_5 = _int(rec, "corruption_5prime")
    corruption_3 = _int(rec, "corruption_3prime")

    seq_id = _str(rec, "sequence_id", "sequence")
    page_title = title or "Sequence Dissection"

    # Build segments
    segments = [
        {"id": "V", "label": "V", "full_label": v_call, "color": SEGMENT_COLORS["V"], "start": v_seq_start, "end": v_seq_end},
        {"id": "NP1", "label": "NP1", "full_label": "N/P Region 1", "color": SEGMENT_COLORS["NP1"], "start": v_seq_end, "end": d_seq_start},
        {"id": "D", "label": "D", "full_label": d_call, "color": SEGMENT_COLORS["D"], "start": d_seq_start, "end": d_seq_end},
        {"id": "NP2", "label": "NP2", "full_label": "N/P Region 2", "color": SEGMENT_COLORS["NP2"], "start": d_seq_end, "end": j_seq_start},
        {"id": "J", "label": "J", "full_label": j_call, "color": SEGMENT_COLORS["J"], "start": j_seq_start, "end": j_seq_end},
    ]

    # ── Build HTML sections ──────────────────────────────────

    # Top bar
    badge_cls = "badge badge-pos" if productive else "badge badge-neg"
    badge_text = "Productive" if productive else "Non-productive"
    top_bar = (
        f'<div class="top-bar">'
        f'<span class="dna-icon">&#x1F9EC;</span>'
        f'<div><h1>{_esc(page_title)}</h1></div>'
        f'<span class="{badge_cls}">{badge_text}</span>'
        f'<span class="seq-id">{_esc(seq_id)}</span>'
        f"</div>"
    )

    # Summary grid
    mut_color = RED if total_mutations > 0 else None
    prod_color = RED if not productive else None
    summary_items = [
        ("V-GENE", v_call, SEGMENT_COLORS["V"], False),
        ("D-GENE", d_call, SEGMENT_COLORS["D"], False),
        ("J-GENE", j_call, SEGMENT_COLORS["J"], False),
        ("JUNCTION AA", junction_aa, None, True),
        ("TOTAL LENGTH", f"{seq_len} nt", None, False),
        ("MUTATIONS", f"{total_mutations} ({mutation_rate * 100:.1f}%)", mut_color, False),
        ("CDR3 LENGTH", f"{cdr3_length} aa", None, False),
        ("PRODUCTIVE", "Yes" if productive else "No", prod_color, False),
    ]
    summary_cells = []
    for label, val, color, mono in summary_items:
        style = f' style="color:{color}"' if color else ""
        mono_cls = " mono" if mono else ""
        summary_cells.append(
            f'<div class="summary-cell">'
            f'<span class="metric-label">{_esc(label)}</span>'
            f'<span class="metric-value{mono_cls}"{style}>{_esc(str(val))}</span>'
            f"</div>"
        )
    summary_grid = '<div class="summary-grid">' + "".join(summary_cells) + "</div>"

    # Corruption warnings
    corruption_html = ""
    if corruption_5 or corruption_3:
        parts = []
        if corruption_5:
            parts.append(f'<div class="corruption-badge">&#9888; 5\' Corruption: +{corruption_5} nt added</div>')
        if corruption_3:
            parts.append(f'<div class="corruption-badge">&#9888; 3\' Corruption: &minus;{corruption_3} nt removed</div>')
        corruption_html = '<div class="corruption-row">' + "".join(parts) + "</div>"

    # Assembled bar + junction
    bar_border_left = f"border-left:3px dashed {RED};" if corruption_5 else ""
    bar_border_right = f"border-right:3px dashed {RED};" if corruption_3 else ""
    if bar_border_left or bar_border_right:
        bar_html = _build_segment_bar_html(segments, seq_len, mutations)
        bar_html = bar_html.replace(
            'class="assembled-bar"',
            f'class="assembled-bar" style="{bar_border_left}{bar_border_right}"',
        )
    else:
        bar_html = _build_segment_bar_html(segments, seq_len, mutations)

    junction_html = _build_junction_bracket(junction_start, junction_end, junction_aa, seq_len)

    # Exploded segment panels
    panels = []
    for seg in segments:
        panels.append(_segment_panel(seg, rec, mutations, sequence))

    # Germline alignment (if available)
    germline_html = ""
    if germline and len(germline) >= len(sequence):
        germline_html = _build_germline_alignment(sequence, germline)

    # Footer
    footer = '<div class="footer">Generated by GenAIRR &mdash; Synthetic AIRR Sequence Simulator</div>'

    # ── Assemble HTML ────────────────────────────────────────

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{_esc(page_title)} — {_esc(seq_id)}</title>
<style>{CSS}</style>
</head>
<body>
<div class="container">
{top_bar}
{summary_grid}
{corruption_html}

<div>
<span class="section-label">Assembled Sequence</span>
{bar_html}
{junction_html}
</div>

<div>
<span class="section-label">&#9889; Exploded Segments</span>
<div class="segments-grid">
{"".join(panels)}
</div>
</div>

{germline_html}
{footer}
</div>
</body>
</html>"""

    path.write_text(html_content, encoding="utf-8")
    return path


def _build_germline_alignment(sequence: str, germline: str) -> str:
    """Build a germline vs sequence alignment visualization."""
    # Show first 200 positions to keep it readable
    show_len = min(len(sequence), len(germline), 400)
    if show_len == 0:
        return ""

    rows = []
    for row_start in range(0, show_len, 80):
        chunk_seq = sequence[row_start : row_start + 80]
        chunk_germ = germline[row_start : row_start + 80]

        seq_chars = []
        germ_chars = []
        match_chars = []
        for s, g in zip(chunk_seq, chunk_germ):
            if g == "." or g == "N":
                germ_chars.append(f'<span class="germline-gap">{_esc(g)}</span>')
                seq_chars.append(f'<span class="germline-gap">{_esc(s)}</span>')
                match_chars.append(f'<span class="germline-gap"> </span>')
            elif s == g:
                germ_chars.append(f'<span class="germline-match">{_esc(g)}</span>')
                seq_chars.append(f'<span class="germline-match">{_esc(s)}</span>')
                match_chars.append(f'<span class="germline-match">|</span>')
            else:
                germ_chars.append(f'<span class="germline-mismatch">{_esc(g)}</span>')
                seq_chars.append(f'<span class="germline-mismatch">{_esc(s)}</span>')
                match_chars.append(f'<span class="germline-mismatch">*</span>')

        pos_label = str(row_start + 1)
        rows.append(
            f'<div class="germline-row">'
            f'<span class="germline-label">{pos_label}</span>'
            f'<div class="germline-chars">{"".join(germ_chars)}</div>'
            f'</div>'
            f'<div class="germline-row">'
            f'<span class="germline-label"></span>'
            f'<div class="germline-chars">{"".join(match_chars)}</div>'
            f'</div>'
            f'<div class="germline-row" style="margin-bottom:0.5rem">'
            f'<span class="germline-label"></span>'
            f'<div class="germline-chars">{"".join(seq_chars)}</div>'
            f'</div>'
        )

    return (
        '<div class="germline-section">'
        '<span class="section-label">Germline Alignment</span>'
        + "".join(rows)
        + "</div>"
    )
