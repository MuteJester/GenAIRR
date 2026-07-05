"""Segment colors and the CSS stylesheet for the HTML render."""

from __future__ import annotations

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
