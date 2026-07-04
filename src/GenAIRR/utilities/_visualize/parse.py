"""Field extraction / parsing helpers and HTML escaping."""

from __future__ import annotations

import html
import json
from typing import Any, Dict


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


def _esc(s: str) -> str:
    return html.escape(s, quote=True)
