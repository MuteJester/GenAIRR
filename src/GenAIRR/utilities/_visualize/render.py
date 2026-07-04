"""Public renderer: assemble a full standalone HTML dissection page."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Union

from .alignment import _build_germline_alignment
from .components import (
    _build_junction_bracket,
    _build_segment_bar_html,
    _segment_panel,
)
from .parse import _bool, _esc, _float, _int, _parse_mutations, _str
from .styles import CSS, RED, SEGMENT_COLORS


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
