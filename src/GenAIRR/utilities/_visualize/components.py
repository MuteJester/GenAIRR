"""HTML component builders for the exploded-view render."""

from __future__ import annotations

from typing import Dict, Optional

from .parse import _bool, _esc, _int, _str
from .styles import PURPLE, RED


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
