"""Canonical accessors for the v2.0.0 RefDataConfig surface.

The engine's RefDataConfig uses indexed accessors and methods (not
properties); this module hides those quirks behind a small set of
helpers used by the MCP tool implementations:

    resolve_refdata(name)              -> RefDataConfig (raises MCPError on unknown)
    iter_alleles(rd, segment)          -> Iterator[Allele]
    find_allele(rd, segment, name)     -> Optional[Allele]
    locus_from_config_name(name)       -> Optional[str]
"""
from __future__ import annotations

from typing import Any, Iterator, Optional

import GenAIRR as ga
from GenAIRR.mcp_errors import CONFIG_NOT_FOUND, MCPError


def resolve_refdata(config_name: str) -> Any:
    """Resolve a config alias (e.g. 'human_igh') to a RefDataConfig.

    Goes through ga.Experiment.on(name), which accepts the friendly
    aliases users see in the README and docs. Raises MCPError with
    CONFIG_NOT_FOUND when the alias isn't recognised.
    """
    try:
        return ga.Experiment.on(config_name).refdata
    except (ValueError, KeyError, FileNotFoundError) as exc:
        raise MCPError(
            CONFIG_NOT_FOUND,
            f"Config {config_name!r} not found.",
            hint="Call list_configs to see available configs.",
        ) from exc


def _pool_size(rd: Any, segment: str) -> int:
    """Return the per-segment pool size. Missing pool (e.g. D on VJ
    chain) returns 0 rather than raising."""
    method_name = f"{segment}_pool_size"
    method = getattr(rd, method_name, None)
    if method is None:
        return 0
    try:
        return int(method())
    except Exception:  # noqa: BLE001 -- defensive against future engine changes
        return 0


def _indexed_accessor(rd: Any, segment: str) -> Optional[Any]:
    """Return the per-segment indexed lookup method (rd.v_allele etc.)
    or None when the pool isn't present."""
    return getattr(rd, f"{segment}_allele", None)


def iter_alleles(rd: Any, segment: str) -> Iterator[Any]:
    """Yield every Allele in the segment's pool (V, D, J, or C).

    Returns an empty iterator when the pool isn't present in the
    config (e.g. asking for D on a VJ chain like human_igk).
    """
    size = _pool_size(rd, segment)
    if size == 0:
        return
    indexed = _indexed_accessor(rd, segment)
    if indexed is None:
        return
    for i in range(size):
        yield indexed(i)


def find_allele(rd: Any, segment: str, name: str) -> Optional[Any]:
    """Look up an allele by exact name in the segment's pool.

    Returns the Allele on match; returns None when not found (callers
    typically raise MCPError(ALLELE_NOT_FOUND) at the call site so the
    error message can include the config name + allele name).
    """
    for allele in iter_alleles(rd, segment):
        if allele.name == name:
            return allele
    return None


def locus_from_config_name(name: str) -> Optional[str]:
    """Derive the locus token from a config alias.

    Convention: GenAIRR aliases are 'species_locus' (e.g. 'human_igh',
    'mouse_tcrb', 'rabbit_igk'). Returns the locus part uppercased.
    Returns None when the name doesn't have the expected shape.

    Used for the validate_records check that record['locus'] matches
    the config's locus, since RefDataConfig itself doesn't expose a
    locus attribute.
    """
    if not name or "_" not in name:
        return None
    parts = name.split("_")
    if len(parts) < 2 or not parts[1]:
        return None
    return parts[1].upper()
