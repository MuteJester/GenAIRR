"""Hard-coded parameter bundles for the simulate_preset MCP tool.

Each preset is a dict of parameters consumed directly by simulate_repertoire
-- no other indirection. Adding a preset means appending an entry to PRESETS
and adding a happy-path test in tests/test_mcp_server.py.

Per the spec (section 3), presets give an LLM agent a way to call
"realistic scenario by name" without needing to know every parameter knob:

    naive_b_cell           -- human_igh, productive, no SHM, no corruption
    memory_b_cell_shm      -- human_igh, productive, s5f SHM count=(5, 15)
    tcr_beta_pool          -- human_tcrb, productive, no SHM (TCRs don't mutate)
    low_quality_sequencing -- human_igh, no constraint, heavy corruption load
    clonal_expansion       -- human_igh, productive, 20 x 25 clonal lineages
"""
from __future__ import annotations

from typing import Any, Dict

from GenAIRR.mcp_errors import INVALID_PRESET, MCPError


PRESETS: Dict[str, Dict[str, Any]] = {
    "naive_b_cell": {
        "config": "human_igh",
        "n": 500,
        "productive_only": True,
        # No SHM, no corruption -- naive B cells haven't seen antigen.
    },
    "memory_b_cell_shm": {
        "config": "human_igh",
        "n": 500,
        "productive_only": True,
        "mutation_model": "s5f",
        "mutation_count_min": 5,
        "mutation_count_max": 15,
    },
    "tcr_beta_pool": {
        "config": "human_tcrb",
        "n": 500,
        "productive_only": True,
        # TCRs don't undergo SHM in physiology, so no mutate clause.
    },
    "low_quality_sequencing": {
        "config": "human_igh",
        "n": 500,
        "productive_only": False,
        "mutation_model": "s5f",
        "mutation_count_min": 0,
        "mutation_count_max": 5,
        "five_prime_loss_max": 10,
        "three_prime_loss_max": 10,
        "pcr_error_count_max": 5,
        "indel_count_max": 3,
        "n_injection_count_max": 5,
        "quality_count_max": 8,
    },
    "clonal_expansion": {
        "config": "human_igh",
        "productive_only": True,
        "n_clones": 20,
        "clone_size": 25,
        # Per-descendant SHM after the clonal fork.
        "mutation_model": "s5f",
        "mutation_count_min": 3,
        "mutation_count_max": 10,
        "pcr_error_count_max": 2,
    },
}


def resolve_preset(name: str) -> Dict[str, Any]:
    """Look up a preset bundle by name. Raises MCPError(INVALID_PRESET)
    when the name isn't known."""
    if name not in PRESETS:
        raise MCPError(
            INVALID_PRESET,
            f"Unknown preset {name!r}.",
            hint=f"Available presets: {sorted(PRESETS.keys())}.",
        )
    # Return a shallow copy so the caller can override fields safely.
    return dict(PRESETS[name])
