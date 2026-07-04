"""Germline vs. sequence alignment visualization."""

from __future__ import annotations

from .parse import _esc


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
