"""
GenAIRR MCP Server — expose AIRR sequence simulation as tools for any agent.

Provides 21 tools via the Model Context Protocol:

Generation:
  1. list_configs      — discover available species/chain configurations
  2. simulate          — generate AIRR sequences with full parameter control
  3. simulate_preset   — generate using named experimental presets
  4. narrate           — human-readable execution trace of one simulation
  5. config_info       — metadata about a specific configuration

Introspection:
  6. inspect_allele    — look up a specific allele's properties
  7. list_alleles      — survey alleles in a config by segment
  8. query_config      — deep-dive into config probability distributions

Validation:
  9. validate_record   — run 9 consistency checks on one AIRR record
 10. validate_batch    — simulate + validate every record (audit tool)

Germline Analysis:
 11. align_to_germline — compare sequence to germline, find unreported/phantom mutations
 12. score_allele_calls— re-score allele calls against full reference

Sequence Analysis:
 13. analyze_mutations — SHM pattern analysis (regions, Ti/Tv, hotspots)
 14. classify_regions  — map every position to IMGT structural regions

Targeted Simulation:
 15. simulate_allele   — lock specific V/D/J alleles for focused testing

Statistics:
 16. summarize_dataset — aggregate stats without returning individual records

Deep Introspection (Pipeline Hooks):
 17. introspect_sequence — simulate with hooks to inspect ASeq at pipeline stages
 18. inspect_flags       — batch flag audit across N sequences
 19. inspect_codon_rail  — inspect the codon rail at a specific pipeline stage
 20. diff_snapshots      — compare ASeq state between two pipeline stages
 21. replay_with_trace   — full pipeline replay with trace + snapshot timeline

Run with::

    python -m GenAIRR.mcp_server          # STDIO transport
    genairr-mcp                           # console entry point

Connect from Claude Code (.claude/mcp.json)::

    {
        "mcpServers": [
            {
                "name": "genairr",
                "type": "stdio",
                "command": "genairr-mcp"
            }
        ]
    }
"""

from __future__ import annotations

import sys
from typing import Any, Dict, List, Optional

from fastmcp import FastMCP

# ── C backend availability ──────────────────────────────────────

_C_AVAILABLE = True
_C_ERROR = ""

try:
    from GenAIRR._native import CSimulator  # noqa: F401
except ImportError as _e:
    _C_AVAILABLE = False
    _C_ERROR = str(_e)


def _check_backend() -> Optional[str]:
    """Return an error message if the C backend is unavailable."""
    if not _C_AVAILABLE:
        return (
            f"GenAIRR C backend not available: {_C_ERROR}\n"
            "Build it with: cd src/GenAIRR/_native/csrc && "
            "mkdir -p build && cd build && cmake .. && make"
        )
    return None


# ── Preset definitions ──────────────────────────────────────────

PRESETS: Dict[str, Dict[str, Any]] = {
    "naive": {
        # Clean V(D)J rearrangements, no mutation, no artifacts
    },
    "memory": {
        "min_mutation_rate": 0.02,
        "max_mutation_rate": 0.08,
        "mutation_model": "s5f",
        "productive": True,
    },
    "illumina": {
        "min_mutation_rate": 0.02,
        "max_mutation_rate": 0.08,
        "mutation_model": "s5f",
        "productive": True,
        "corrupt_5prime": True,
        "corrupt_3prime": True,
        "paired_end_read_length": 300,
        "quality_base_error": 0.001,
        "quality_peak_error": 0.02,
    },
    "nanopore": {
        "min_mutation_rate": 0.02,
        "max_mutation_rate": 0.08,
        "mutation_model": "s5f",
        "productive": True,
        "long_read_error_rate": 0.03,
    },
    "noisy": {
        "min_mutation_rate": 0.02,
        "max_mutation_rate": 0.08,
        "mutation_model": "s5f",
        "productive": True,
        "corrupt_5prime": True,
        "corrupt_3prime": True,
        "paired_end_read_length": 300,
        "quality_base_error": 0.001,
        "quality_peak_error": 0.02,
        "contaminant_rate": 0.01,
        "indel_prob": 0.005,
        "n_prob": 0.005,
    },
    "tcr": {
        "productive": True,
    },
}


# ── Helpers ─────────────────────────────────────────────────────

def _build_experiment(
    config: str,
    *,
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    mutation_model: Optional[str] = None,
    with_isotype_rates: bool = False,
    selection_strength: Optional[float] = None,
    d_inversion_prob: Optional[float] = None,
    receptor_revision_prob: Optional[float] = None,
    corrupt_5prime: bool = False,
    corrupt_3prime: bool = False,
    paired_end_read_length: Optional[int] = None,
    long_read_error_rate: Optional[float] = None,
    reverse_complement_prob: Optional[float] = None,
    quality_base_error: Optional[float] = None,
    quality_peak_error: Optional[float] = None,
    umi_length: Optional[int] = None,
    pcr_error_rate: Optional[float] = None,
    pcr_cycles: Optional[int] = None,
    primer_mask_length: Optional[int] = None,
    contaminant_rate: Optional[float] = None,
    indel_prob: Optional[float] = None,
    n_prob: Optional[float] = None,
):
    """Translate flat MCP parameters into an Experiment object."""
    from GenAIRR.experiment import Experiment
    from GenAIRR.ops import (
        rate, model, with_antigen_selection, with_isotype_rates as isotype_rates,
        with_d_inversion, with_receptor_revision,
        with_5prime_loss, with_3prime_loss, paired_end, long_read,
        with_quality_profile, with_reverse_complement,
        with_umi, with_pcr, with_primer_mask,
        with_contaminants, with_indels, with_ns,
    )

    exp = Experiment.on(config)

    # ── Recombine phase ──
    recombine_clauses = []
    if d_inversion_prob is not None:
        recombine_clauses.append(with_d_inversion(d_inversion_prob))
    if receptor_revision_prob is not None:
        recombine_clauses.append(with_receptor_revision(receptor_revision_prob))
    if recombine_clauses:
        exp = exp.recombine(*recombine_clauses)

    # ── Mutate phase ──
    mutate_clauses = []
    if min_mutation_rate is not None or max_mutation_rate is not None:
        mutate_clauses.append(rate(
            min_mutation_rate if min_mutation_rate is not None else 0.01,
            max_mutation_rate if max_mutation_rate is not None else 0.05,
        ))
        mutate_clauses.append(model(mutation_model or "s5f"))
    if with_isotype_rates:
        mutate_clauses.append(isotype_rates())
    if selection_strength is not None:
        mutate_clauses.append(with_antigen_selection(selection_strength))
    if mutate_clauses:
        exp = exp.mutate(*mutate_clauses)

    # ── Prepare phase ──
    prepare_clauses = []
    if primer_mask_length is not None:
        prepare_clauses.append(with_primer_mask(primer_mask_length))
    if umi_length is not None:
        prepare_clauses.append(with_umi(umi_length))
    if pcr_error_rate is not None or pcr_cycles is not None:
        prepare_clauses.append(with_pcr(
            error_rate=pcr_error_rate if pcr_error_rate is not None else 1e-4,
            cycles=pcr_cycles if pcr_cycles is not None else 30,
        ))
    if prepare_clauses:
        exp = exp.prepare(*prepare_clauses)

    # ── Sequence phase ──
    sequence_clauses = []
    if corrupt_5prime:
        sequence_clauses.append(with_5prime_loss())
    if corrupt_3prime:
        sequence_clauses.append(with_3prime_loss())
    if paired_end_read_length is not None:
        sequence_clauses.append(paired_end(paired_end_read_length))
    if long_read_error_rate is not None:
        sequence_clauses.append(long_read(error_rate=long_read_error_rate))
    if quality_base_error is not None or quality_peak_error is not None:
        sequence_clauses.append(with_quality_profile(
            base=quality_base_error if quality_base_error is not None else 0.001,
            peak=quality_peak_error if quality_peak_error is not None else 0.02,
        ))
    if reverse_complement_prob is not None:
        sequence_clauses.append(with_reverse_complement(reverse_complement_prob))
    if sequence_clauses:
        exp = exp.sequence(*sequence_clauses)

    # ── Observe phase ──
    observe_clauses = []
    if contaminant_rate is not None:
        observe_clauses.append(with_contaminants(contaminant_rate))
    if indel_prob is not None:
        observe_clauses.append(with_indels(indel_prob))
    if n_prob is not None:
        observe_clauses.append(with_ns(n_prob))
    if observe_clauses:
        exp = exp.observe(*observe_clauses)

    return exp


def _filter_fields(
    records: List[Dict[str, Any]],
    fields: Optional[List[str]],
) -> List[Dict[str, Any]]:
    """Filter AIRR records to requested fields only."""
    if not fields:
        return records
    return [{k: r[k] for k in fields if k in r} for r in records]


def _run_simulation(
    config: str = "human_igh",
    n: int = 10,
    seed: Optional[int] = None,
    productive: bool = False,
    fields: Optional[List[str]] = None,
    **kwargs,
) -> Dict[str, Any]:
    """Shared simulation logic for simulate and simulate_preset."""
    err = _check_backend()
    if err:
        return {"error": err}

    n = max(1, min(n, 10000))

    try:
        exp = _build_experiment(config, **kwargs)
        result = exp.run(n=n, seed=seed, productive=productive)
        records = [dict(r) for r in result]
        records = _filter_fields(records, fields)
        field_names = list(records[0].keys()) if records else []
        return {
            "count": len(records),
            "fields": field_names,
            "records": records,
        }
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


# ── MCP Server ──────────────────────────────────────────────────

mcp = FastMCP(
    "GenAIRR",
    instructions=(
        "GenAIRR simulates synthetic immune receptor sequences in AIRR format. "
        "Use list_configs() to discover species/chains, then simulate() to "
        "generate sequences. Use simulate_preset() for common experimental "
        "scenarios without needing to know parameter details.\n\n"
        "For deeper analysis: inspect_allele/list_alleles to explore germline references, "
        "validate_record/validate_batch to audit record consistency, "
        "align_to_germline/score_allele_calls to verify metadata accuracy, "
        "analyze_mutations/classify_regions for SHM pattern analysis, "
        "simulate_allele to test specific alleles, "
        "summarize_dataset for aggregate statistics.\n\n"
        "For deep introspection: introspect_sequence to place hooks at pipeline stages "
        "and inspect the full ASeq linked list state, diff_snapshots to compare two "
        "pipeline stages, inspect_flags for batch flag audits, inspect_codon_rail for "
        "reading frame analysis, replay_with_trace for full pipeline replay."
    ),
)


@mcp.tool()
def list_configs(
    species_filter: Optional[str] = None,
    chain_filter: Optional[str] = None,
) -> List[str]:
    """List available species/chain configurations for AIRR sequence simulation.

    GenAIRR has 106+ configs across 23 species (human, mouse, rat, rabbit,
    rhesus, cow, dog, cat, pig, sheep, etc.) and chain types (IGH, IGK, IGL,
    TCRA, TCRB, TCRD, TCRG).

    Returns config names to pass as the 'config' parameter to simulate().

    Common shorthand aliases also work: "human_igh", "mouse_igk", "rabbit_tcrb".

    Args:
        species_filter: Case-insensitive substring to filter by species (e.g. "mouse").
        chain_filter: Case-insensitive substring to filter by chain type (e.g. "igh", "tcr").
    """
    from GenAIRR.data import list_configs as _list

    configs = _list()
    if species_filter:
        sf = species_filter.upper()
        configs = [c for c in configs if sf in c]
    if chain_filter:
        cf = chain_filter.upper()
        configs = [c for c in configs if cf in c]
    return configs


@mcp.tool()
def simulate(
    config: str = "human_igh",
    n: int = 10,
    seed: Optional[int] = None,
    productive: bool = False,
    # Mutation
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    mutation_model: Optional[str] = None,
    with_isotype_rates: bool = False,
    selection_strength: Optional[float] = None,
    # Recombination
    d_inversion_prob: Optional[float] = None,
    receptor_revision_prob: Optional[float] = None,
    # Sequencing
    corrupt_5prime: bool = False,
    corrupt_3prime: bool = False,
    paired_end_read_length: Optional[int] = None,
    reverse_complement_prob: Optional[float] = None,
    quality_base_error: Optional[float] = None,
    quality_peak_error: Optional[float] = None,
    # Library prep
    umi_length: Optional[int] = None,
    pcr_error_rate: Optional[float] = None,
    pcr_cycles: Optional[int] = None,
    primer_mask_length: Optional[int] = None,
    # Post-sequencing
    contaminant_rate: Optional[float] = None,
    indel_prob: Optional[float] = None,
    n_prob: Optional[float] = None,
    # Output
    fields: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Simulate synthetic AIRR immune receptor sequences.

    Generates realistic V(D)J rearranged sequences with optional somatic
    hypermutation, sequencing artifacts, and noise. Returns AIRR-format
    records with full annotations.

    Quick start: simulate(config="human_igh", n=10)

    For mutated sequences: simulate(config="human_igh", n=100,
        min_mutation_rate=0.02, max_mutation_rate=0.08)

    For realistic Illumina: simulate(config="human_igh", n=100,
        min_mutation_rate=0.02, max_mutation_rate=0.08,
        corrupt_5prime=True, paired_end_read_length=300)

    Args:
        config: Species/chain config name (use list_configs() to see options).
                Common: "human_igh", "human_igk", "mouse_igh", "human_tcrb".
        n: Number of sequences to generate (1-10000).
        seed: Random seed for reproducibility.
        productive: Only return productive rearrangements (in-frame, no stop codons).
        min_mutation_rate: Minimum SHM rate (e.g. 0.01). Setting this enables mutation.
        max_mutation_rate: Maximum SHM rate (e.g. 0.08).
        mutation_model: "s5f" (context-dependent, default) or "uniform".
        with_isotype_rates: Adjust mutation rates by isotype (class-switch recombination).
        selection_strength: Antigen selection pressure intensity (0-1).
        d_inversion_prob: Probability of D-gene inversion during recombination (0-1).
        receptor_revision_prob: Probability of V-gene replacement (0-1).
        corrupt_5prime: Simulate 5' end signal loss (truncation).
        corrupt_3prime: Simulate 3' end signal loss.
        paired_end_read_length: Enable paired-end reads (common: 150, 250, 300).
        reverse_complement_prob: Fraction of reads in antisense orientation (0-1).
        quality_base_error: Base sequencing error rate at 5' end (e.g. 0.001).
        quality_peak_error: Peak sequencing error rate at 3' end (e.g. 0.02).
        umi_length: Prepend random UMI barcode of this length (common: 8, 12, 16).
        pcr_error_rate: Per-base per-cycle PCR error rate (e.g. 1e-4).
        pcr_cycles: Number of PCR amplification cycles (e.g. 30).
        primer_mask_length: Overwrite N 5' bases with germline (0 = full FR1).
        contaminant_rate: Per-sequence contamination probability (0-1).
        indel_prob: Per-position indel probability (0-1).
        n_prob: Per-position ambiguous-N probability (0-1).
        fields: Subset of AIRR fields to return. None returns all ~47 fields.
                Useful subsets: ["sequence", "v_call", "d_call", "j_call",
                "junction_aa", "productive", "mutation_rate"]
    """
    return _run_simulation(
        config=config, n=n, seed=seed, productive=productive,
        fields=fields,
        min_mutation_rate=min_mutation_rate,
        max_mutation_rate=max_mutation_rate,
        mutation_model=mutation_model,
        with_isotype_rates=with_isotype_rates,
        selection_strength=selection_strength,
        d_inversion_prob=d_inversion_prob,
        receptor_revision_prob=receptor_revision_prob,
        corrupt_5prime=corrupt_5prime,
        corrupt_3prime=corrupt_3prime,
        paired_end_read_length=paired_end_read_length,
        reverse_complement_prob=reverse_complement_prob,
        quality_base_error=quality_base_error,
        quality_peak_error=quality_peak_error,
        umi_length=umi_length,
        pcr_error_rate=pcr_error_rate,
        pcr_cycles=pcr_cycles,
        primer_mask_length=primer_mask_length,
        contaminant_rate=contaminant_rate,
        indel_prob=indel_prob,
        n_prob=n_prob,
    )


@mcp.tool()
def simulate_preset(
    preset: str,
    config: str = "human_igh",
    n: int = 10,
    seed: Optional[int] = None,
    fields: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Generate sequences using a named experimental preset.

    Presets capture common real-world scenarios so you don't need to know
    individual parameter values. Available presets:

    - "naive":    Unmutated V(D)J rearrangements (germline B/T cells)
    - "memory":   Mutated sequences (2-8% SHM, S5F model, productive only)
    - "illumina": Realistic Illumina sequencing (memory + paired-end 300bp,
                  quality errors, 5'/3' loss)
    - "nanopore": Long-read sequencing (memory + homopolymer errors)
    - "noisy":    Kitchen-sink noise (illumina + contaminants + indels + Ns)
    - "tcr":      T-cell receptor preset (naive, productive only)

    Args:
        preset: One of: "naive", "memory", "illumina", "nanopore", "noisy", "tcr".
        config: Species/chain config name.
        n: Number of sequences (1-10000).
        seed: Random seed for reproducibility.
        fields: Subset of AIRR fields to return (None = all).
    """
    if preset not in PRESETS:
        return {
            "error": f"Unknown preset {preset!r}. "
                     f"Available: {sorted(PRESETS.keys())}"
        }

    params = dict(PRESETS[preset])
    productive = params.pop("productive", False)
    return _run_simulation(
        config=config, n=n, seed=seed, productive=productive,
        fields=fields, **params,
    )


@mcp.tool()
def narrate_simulation(
    config: str = "human_igh",
    seed: int = 42,
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    productive: bool = False,
) -> str:
    """Simulate one sequence and return a detailed human-readable narrative.

    Shows every internal step: allele selection, exonuclease trimming,
    N/P nucleotide addition, junction assembly, productivity assessment,
    somatic hypermutation events (each mutation position, context, base change),
    and any sequencing artifacts.

    Useful for understanding what GenAIRR does, debugging, or explaining
    a simulation result to a user.

    Args:
        config: Species/chain config name.
        seed: Random seed for reproducibility.
        min_mutation_rate: Enable SHM with this minimum rate.
        max_mutation_rate: Maximum SHM rate.
        productive: Only generate productive rearrangements.
    """
    err = _check_backend()
    if err:
        return err

    try:
        from GenAIRR.result import narrate as _narrate

        exp = _build_experiment(
            config,
            min_mutation_rate=min_mutation_rate,
            max_mutation_rate=max_mutation_rate,
        )
        if productive:
            # narrate() doesn't take productive directly; compile handles it
            sim = exp.compile(seed=seed, productive=True)
            csim = sim._sim
            csim.set_trace(True)
            csim.simulate_one()
            trace_text = csim.get_trace()
            csim.set_trace(False)
            # Return raw trace (narrate formats it, but we'd need the full
            # function; simpler to just call narrate with a non-productive exp)
            return _narrate(exp, seed=seed, color=False)

        return _narrate(exp, seed=seed, color=False)
    except Exception as e:
        return f"Error: {type(e).__name__}: {e}"


@mcp.tool()
def config_info(config: str) -> Dict[str, Any]:
    """Get detailed information about a specific species/chain configuration.

    Returns metadata including species, chain type, whether the chain has
    a D segment, the reference database source, and counts of V/D/J/C
    alleles in the germline pool.

    Args:
        config: Config name (e.g. "human_igh", "MOUSE_IGK_IMGT",
                "rabbit_tcrb"). Use list_configs() to see all options.
    """
    try:
        from GenAIRR.protocol import _resolve_config

        dc = _resolve_config(config)
        meta = dc.metadata

        def _count_alleles(allele_dict):
            if not allele_dict:
                return 0
            return sum(len(v) for v in allele_dict.values())

        return {
            "config": config,
            "species": meta.species.value if meta.species else "unknown",
            "chain_type": meta.chain_type.value if meta.chain_type else "unknown",
            "has_d": meta.has_d,
            "reference_set": meta.reference_set,
            "last_updated": str(meta.last_updated),
            "n_v_alleles": _count_alleles(dc.v_alleles),
            "n_d_alleles": _count_alleles(dc.d_alleles),
            "n_j_alleles": _count_alleles(dc.j_alleles),
            "n_c_alleles": _count_alleles(dc.c_alleles),
        }
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


# ── Category 1: Reference Introspection ────────────────────────


@mcp.tool()
def inspect_allele(config: str, allele_name: str) -> Dict[str, Any]:
    """Look up a specific allele's properties from the germline reference.

    Returns the allele's sequence (preview), ungapped length, anchor position,
    gene family, and segment type. Useful for investigating individual alleles
    (e.g. checking if TRBV17*01 has a valid anchor).

    Args:
        config: Species/chain config name.
        allele_name: Full allele name (e.g. "IGHV1-2*01", "TRBV17*01").
    """
    try:
        from GenAIRR.protocol import _resolve_config
        from GenAIRR.utilities.mcp_helpers import find_allele, format_allele_info

        dc = _resolve_config(config)
        allele = find_allele(dc, allele_name)
        if allele is None:
            return {"error": f"Allele {allele_name!r} not found in {config}"}
        return format_allele_info(allele)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


@mcp.tool()
def list_alleles(
    config: str,
    segment: str = "v",
    gene_filter: Optional[str] = None,
    limit: int = 50,
) -> Dict[str, Any]:
    """List alleles in a configuration by segment type.

    Survey the V, D, J, or C alleles with optional gene-name filtering.
    Returns a compact summary of each allele (name, length, anchor, family).

    Args:
        config: Species/chain config name.
        segment: Segment type: "v", "d", "j", or "c".
        gene_filter: Case-insensitive substring to filter allele names (e.g. "TRBV17").
        limit: Maximum number of alleles to return (default 50).
    """
    try:
        from GenAIRR.protocol import _resolve_config
        from GenAIRR.utilities.mcp_helpers import flatten_alleles, format_allele_info

        dc = _resolve_config(config)
        seg = segment.lower()
        allele_dict = getattr(dc, f"{seg}_alleles", None)
        if allele_dict is None:
            return {"error": f"No {seg.upper()} alleles in {config}"}

        flat = flatten_alleles(allele_dict)
        if gene_filter:
            gf = gene_filter.upper()
            flat = {k: v for k, v in flat.items() if gf in k.upper()}

        total = len(flat)
        alleles = [format_allele_info(a) for a in list(flat.values())[:limit]]
        return {"segment": seg.upper(), "total": total, "returned": len(alleles), "alleles": alleles}
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


@mcp.tool()
def query_config(config: str, section: str) -> Dict[str, Any]:
    """Deep-dive into a configuration's internal probability distributions.

    Inspect gene usage probabilities, trimming distributions, NP region
    Markov chain parameters, or P-nucleotide length probabilities.

    Args:
        config: Species/chain config name.
        section: One of: "gene_use", "trimming", "np_params", "p_nucleotides".
    """
    try:
        from GenAIRR.protocol import _resolve_config
        from GenAIRR.utilities.mcp_helpers import extract_config_section

        dc = _resolve_config(config)
        return extract_config_section(dc, section)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


# ── Category 2: Record Validation ─────────────────────────────


@mcp.tool()
def validate_record(record: Dict[str, Any]) -> Dict[str, Any]:
    """Run 9 consistency checks on a single AIRR record.

    Checks: nucleotide validity, coordinate bounds, segment ordering,
    junction bounds, junction length match, productive consistency,
    mutation count, mutation content (each pos:X>Y verified against
    germline_alignment and sequence), and sequence_length field.

    Pass a record dict from simulate() output.

    Args:
        record: A single AIRR record dictionary.
    """
    try:
        from GenAIRR.utilities.mcp_helpers import validate_record as _validate
        return _validate(record)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


@mcp.tool()
def validate_batch(
    config: str = "human_igh",
    n: int = 100,
    seed: Optional[int] = None,
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    productive: bool = False,
    corrupt_5prime: bool = False,
    indel_prob: Optional[float] = None,
) -> Dict[str, Any]:
    """Simulate sequences and validate every record — the audit tool.

    Generates n sequences and runs all 9 validation checks on each one.
    Returns a summary of pass/fail counts per check and details of any failures.

    Args:
        config: Species/chain config name.
        n: Number of sequences to generate and validate (1-10000).
        seed: Random seed for reproducibility.
        min_mutation_rate: Enable SHM with this minimum rate.
        max_mutation_rate: Maximum SHM rate.
        productive: Only generate productive rearrangements.
        corrupt_5prime: Enable 5' end loss.
        indel_prob: Per-position indel probability.
    """
    try:
        from GenAIRR.utilities.mcp_helpers import validate_record as _validate

        sim_result = _run_simulation(
            config=config, n=n, seed=seed, productive=productive,
            min_mutation_rate=min_mutation_rate,
            max_mutation_rate=max_mutation_rate,
            corrupt_5prime=corrupt_5prime,
            indel_prob=indel_prob,
        )
        if "error" in sim_result:
            return sim_result

        records = sim_result["records"]
        check_summary: Dict[str, Dict[str, int]] = {}
        failures = []

        for rec in records:
            result = _validate(rec)
            for check in result["checks"]:
                name = check["name"]
                if name not in check_summary:
                    check_summary[name] = {"passed": 0, "failed": 0}
                if check["passed"]:
                    check_summary[name]["passed"] += 1
                else:
                    check_summary[name]["failed"] += 1
            if not result["valid"] and len(failures) < 10:
                failures.append({
                    "v_call": rec.get("v_call", ""),
                    "failed_checks": [c for c in result["checks"] if not c["passed"]],
                })

        n_valid = sum(1 for r in records if _validate(r)["valid"])
        return {
            "n_simulated": len(records),
            "n_valid": n_valid,
            "n_failed": len(records) - n_valid,
            "check_summary": check_summary,
            "failures": failures,
        }
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


# ── Category 3: Germline Analysis ─────────────────────────────


@mcp.tool()
def align_to_germline(config: str, record: Dict[str, Any]) -> Dict[str, Any]:
    """Compare a simulated sequence to its germline reference alleles.

    Cross-references the sequence and germline_alignment fields with the
    reported mutations to find:
    - True mutations (seq != germline, reported)
    - Unreported mutations (seq != germline, NOT reported)
    - Phantom mutations (seq == germline, but reported)

    Args:
        config: Species/chain config name.
        record: A single AIRR record dictionary.
    """
    try:
        from GenAIRR.protocol import _resolve_config
        from GenAIRR.utilities.mcp_helpers import align_to_germline as _align

        dc = _resolve_config(config)
        return _align(dc, record)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


@mcp.tool()
def score_allele_calls(config: str, record: Dict[str, Any]) -> Dict[str, Any]:
    """Re-score allele calls against the full germline reference.

    For each segment (V, D, J), scores every allele in the reference
    against the sequence region and returns the top matches. Verifies
    the reported call is the best match (or ties for best).

    Uses N-aware scoring matching C bitmap semantics.

    Args:
        config: Species/chain config name.
        record: A single AIRR record dictionary.
    """
    try:
        from GenAIRR.protocol import _resolve_config
        from GenAIRR.utilities.mcp_helpers import score_all_alleles

        dc = _resolve_config(config)
        result = {}
        for seg in ("v", "d", "j"):
            scores = score_all_alleles(dc, record, seg)
            reported = (record.get(f"{seg}_call", "") or "").split(",")[0].strip()
            reported_is_best = False
            if scores:
                best_frac = scores[0]["fraction"]
                reported_is_best = any(
                    s["is_reported_call"] and s["fraction"] >= best_frac - 0.0001
                    for s in scores
                )
            result[seg] = {
                "reported": reported,
                "top_5": scores[:5],
                "reported_is_best": reported_is_best,
            }
        return result
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


# ── Category 4: Sequence Analysis ─────────────────────────────


@mcp.tool()
def analyze_mutations(
    config: str = "human_igh",
    n: int = 200,
    seed: Optional[int] = None,
    min_mutation_rate: float = 0.02,
    max_mutation_rate: float = 0.08,
) -> Dict[str, Any]:
    """Analyze somatic hypermutation patterns across simulated sequences.

    Generates mutated sequences and analyzes:
    - Per-IMGT-region mutation distribution (FWR1/CDR1/FWR2/CDR2/FWR3/CDR3/FWR4)
    - Transition/transversion ratio
    - SHM hotspot motif targeting (WRCY/RGYW)

    Args:
        config: Species/chain config name.
        n: Number of sequences to generate (1-10000).
        seed: Random seed for reproducibility.
        min_mutation_rate: Minimum SHM rate.
        max_mutation_rate: Maximum SHM rate.
    """
    err = _check_backend()
    if err:
        return {"error": err}

    try:
        from GenAIRR.protocol import _resolve_config
        from GenAIRR.utilities.mcp_helpers import analyze_mutation_patterns

        n = max(1, min(n, 10000))
        dc = _resolve_config(config)
        exp = _build_experiment(
            config,
            min_mutation_rate=min_mutation_rate,
            max_mutation_rate=max_mutation_rate,
        )
        result = exp.run(n=n, seed=seed)
        records = [dict(r) for r in result]
        return analyze_mutation_patterns(records, dc)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


@mcp.tool()
def classify_regions(config: str, record: Dict[str, Any]) -> Dict[str, Any]:
    """Classify every position in a sequence into IMGT structural regions.

    Maps each position to FWR1, CDR1, FWR2, CDR2, FWR3, CDR3, FWR4, or NP.
    Also breaks down mutations by region and provides sequence excerpts.

    Args:
        config: Species/chain config name.
        record: A single AIRR record dictionary.
    """
    try:
        from GenAIRR.protocol import _resolve_config
        from GenAIRR.utilities.mcp_helpers import classify_record_positions

        dc = _resolve_config(config)
        return classify_record_positions(dc, record)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


# ── Category 5: Targeted Simulation ───────────────────────────


@mcp.tool()
def simulate_allele(
    config: str = "human_igh",
    v_allele: Optional[str] = None,
    d_allele: Optional[str] = None,
    j_allele: Optional[str] = None,
    n: int = 10,
    seed: Optional[int] = None,
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    productive: bool = False,
    fields: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Simulate sequences with specific V/D/J alleles locked.

    Forces the simulator to use the specified allele(s) for every sequence.
    Useful for stress-testing specific alleles (e.g. "generate 100 sequences
    using only TRBV17*01 to verify they're all non-productive").

    Args:
        config: Species/chain config name.
        v_allele: Lock this V allele (e.g. "IGHV1-2*01"). None = random.
        d_allele: Lock this D allele. None = random.
        j_allele: Lock this J allele. None = random.
        n: Number of sequences (1-10000).
        seed: Random seed for reproducibility.
        min_mutation_rate: Enable SHM with this minimum rate.
        max_mutation_rate: Maximum SHM rate.
        productive: Only return productive rearrangements.
        fields: Subset of AIRR fields to return.
    """
    err = _check_backend()
    if err:
        return {"error": err}

    try:
        from GenAIRR._native import CSimulator
        from GenAIRR.protocol import _resolve_config
        from GenAIRR.dataconfig.gdc_io import to_gdc_bytes

        n = max(1, min(n, 10000))
        dc = _resolve_config(config)

        # Determine S5F chain category for GDC serialization
        s5f_cat = None
        if min_mutation_rate is not None:
            ct = dc.metadata.chain_type.value if dc.metadata and dc.metadata.chain_type else ""
            if "HEAVY" in ct or "BETA" in ct or "DELTA" in ct:
                s5f_cat = "heavy"
            else:
                s5f_cat = "light"

        gdc_bytes = to_gdc_bytes(dc, s5f_chain_category=s5f_cat)
        sim = CSimulator(gdc_bytes=gdc_bytes)

        # Lock alleles
        if v_allele:
            sim.lock_allele("v", v_allele)
        if d_allele:
            sim.lock_allele("d", d_allele)
        if j_allele:
            sim.lock_allele("j", j_allele)

        # Set features/params
        if seed is not None:
            sim.set_seed(seed)
        if productive:
            sim.set_feature("productive", True)
        if min_mutation_rate is not None:
            sim.set_feature("mutation", True)
            sim.set_param("mutation_rate_min", min_mutation_rate)
            sim.set_param("mutation_rate_max", max_mutation_rate or min_mutation_rate)

        records = sim.simulate(n)
        records = _filter_fields(records, fields)
        field_names = list(records[0].keys()) if records else []
        return {
            "count": len(records),
            "fields": field_names,
            "records": records,
        }
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


# ── Category 6: Statistical Summary ───────────────────────────


@mcp.tool()
def summarize_dataset(
    config: str = "human_igh",
    n: int = 500,
    seed: Optional[int] = None,
    preset: Optional[str] = None,
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    productive: bool = False,
) -> Dict[str, Any]:
    """Generate sequences and return aggregate statistics (no individual records).

    Computes: productive rate, mutation rate stats, junction length distribution,
    V/D/J gene usage (top 20), NP region length stats.

    Much more token-efficient than simulate() for dataset-level analysis.

    Args:
        config: Species/chain config name.
        n: Number of sequences to generate (1-10000).
        seed: Random seed for reproducibility.
        preset: Optional preset name ("naive", "memory", etc.) for defaults.
        min_mutation_rate: Minimum SHM rate (overrides preset).
        max_mutation_rate: Maximum SHM rate (overrides preset).
        productive: Only generate productive rearrangements.
    """
    try:
        from GenAIRR.utilities.mcp_helpers import compute_dataset_summary

        # Merge preset if provided
        kwargs: Dict[str, Any] = {}
        if preset and preset in PRESETS:
            kwargs = dict(PRESETS[preset])
            productive = kwargs.pop("productive", productive)
        if min_mutation_rate is not None:
            kwargs["min_mutation_rate"] = min_mutation_rate
        if max_mutation_rate is not None:
            kwargs["max_mutation_rate"] = max_mutation_rate

        sim_result = _run_simulation(
            config=config, n=n, seed=seed, productive=productive, **kwargs,
        )
        if "error" in sim_result:
            return sim_result

        return compute_dataset_summary(sim_result["records"])
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


# ── Category 7: Deep Introspection (Pipeline Hooks) ────────────


def _build_hooked_simulator(
    config: str,
    seed: Optional[int] = None,
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    productive: bool = False,
    v_allele: Optional[str] = None,
    indel_prob: Optional[float] = None,
    corrupt_5prime: bool = False,
    corrupt_3prime: bool = False,
    n_prob: Optional[float] = None,
):
    """Build a CSimulator configured for hooked simulation."""
    from GenAIRR._native import CSimulator
    from GenAIRR.protocol import _resolve_config
    from GenAIRR.dataconfig.gdc_io import to_gdc_bytes

    dc = _resolve_config(config)

    s5f_cat = None
    if min_mutation_rate is not None:
        ct = dc.metadata.chain_type.value if dc.metadata and dc.metadata.chain_type else ""
        if "HEAVY" in ct or "BETA" in ct or "DELTA" in ct:
            s5f_cat = "heavy"
        else:
            s5f_cat = "light"

    gdc_bytes = to_gdc_bytes(dc, s5f_chain_category=s5f_cat)
    sim = CSimulator(gdc_bytes=gdc_bytes)

    if seed is not None:
        sim.set_seed(seed)
    if productive:
        sim.set_feature("productive", True)
    if min_mutation_rate is not None:
        sim.set_feature("mutate", True)
        sim.set_param("min_mutation_rate", min_mutation_rate)
        sim.set_param("max_mutation_rate", max_mutation_rate or min_mutation_rate)
    if v_allele:
        sim.lock_allele("v", v_allele)
    if indel_prob is not None:
        sim.set_feature("indels", True)
        sim.set_param("indel_prob", indel_prob)
    if corrupt_5prime:
        sim.set_feature("corrupt_5_prime", True)
    if corrupt_3prime:
        sim.set_feature("corrupt_3_prime", True)
    if n_prob is not None:
        sim.set_feature("insert_ns", True)
        sim.set_param("n_prob", n_prob)

    return sim


@mcp.tool()
def introspect_sequence(
    config: str = "human_igh",
    seed: int = 42,
    hooks: Optional[List[str]] = None,
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    productive: bool = False,
    v_allele: Optional[str] = None,
    indel_prob: Optional[float] = None,
    corrupt_5prime: bool = False,
    corrupt_3prime: bool = False,
    compact: bool = True,
) -> Dict[str, Any]:
    """Simulate one sequence with pipeline hooks to inspect internal ASeq state.

    Places hooks at specified pipeline stages and captures a full snapshot
    of the ASeq linked list at each point. Each snapshot shows every
    nucleotide node with its current base, germline base, segment, flags,
    germline position, frame phase, amino acid, and productivity.

    Hook points:
      post_assembly, post_functionality, post_d_inversion, post_receptor_rev,
      post_mutation, post_selection, post_corrupt_5, post_corrupt_3,
      post_indels, post_ns, post_pcr, post_quality, final

    Args:
        config: Species/chain config name.
        seed: Random seed for reproducibility.
        hooks: List of hook point names (default: post_assembly, post_mutation, final).
        min_mutation_rate: Enable SHM with this minimum rate.
        max_mutation_rate: Maximum SHM rate.
        productive: Only generate productive rearrangements.
        v_allele: Lock a specific V allele.
        indel_prob: Enable indels with this probability.
        corrupt_5prime: Enable 5' corruption.
        corrupt_3prime: Enable 3' corruption.
        compact: If True, return summaries instead of full node lists.
    """
    err = _check_backend()
    if err:
        return {"error": err}

    if hooks is None:
        hooks = ["post_assembly", "post_mutation", "final"]

    try:
        from GenAIRR.utilities.mcp_helpers import (
            summarize_aseq, detect_node_anomalies, format_snapshot_timeline,
        )

        sim = _build_hooked_simulator(
            config, seed=seed,
            min_mutation_rate=min_mutation_rate,
            max_mutation_rate=max_mutation_rate,
            productive=productive, v_allele=v_allele,
            indel_prob=indel_prob,
            corrupt_5prime=corrupt_5prime,
            corrupt_3prime=corrupt_3prime,
        )

        rec, snapshots, trace = sim.simulate_one_hooked(hooks)

        # Process snapshots
        processed = []
        for snap in snapshots:
            nodes = snap["nodes"]
            entry: Dict[str, Any] = {
                "hook": snap["hook"],
                "summary": summarize_aseq(nodes),
                "anomalies": detect_node_anomalies(nodes),
            }
            if not compact:
                entry["nodes"] = nodes
            processed.append(entry)

        timeline = format_snapshot_timeline(snapshots)

        return {
            "airr_record": rec,
            "snapshots": processed,
            "timeline": timeline,
            "trace": trace,
        }
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


@mcp.tool()
def inspect_flags(
    config: str = "human_igh",
    n: int = 100,
    seed: Optional[int] = None,
    hooks: Optional[List[str]] = None,
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    productive: bool = False,
    indel_prob: Optional[float] = None,
) -> Dict[str, Any]:
    """Batch flag audit: inspect per-node flags across N sequences.

    Simulates N sequences with hooks and aggregates flag distributions.
    Shows per-flag counts, per-segment flag counts, and total anomalies.

    Args:
        config: Species/chain config name.
        n: Number of sequences (1-1000).
        seed: Random seed.
        hooks: Hook points to inspect (default: ["final"]).
        min_mutation_rate: Enable SHM.
        max_mutation_rate: Maximum SHM rate.
        productive: Only productive rearrangements.
        indel_prob: Enable indels.
    """
    err = _check_backend()
    if err:
        return {"error": err}

    if hooks is None:
        hooks = ["final"]

    n = max(1, min(n, 1000))

    try:
        from GenAIRR.utilities.mcp_helpers import aggregate_flag_stats

        sim = _build_hooked_simulator(
            config, seed=seed,
            min_mutation_rate=min_mutation_rate,
            max_mutation_rate=max_mutation_rate,
            productive=productive, indel_prob=indel_prob,
        )
        sim.set_hooks(hooks)

        all_nodes = []
        for _ in range(n):
            sim.simulate_one()
            snaps = sim.get_snapshots()
            for snap in snaps:
                all_nodes.append(snap["nodes"])

        sim.set_hooks([])
        return aggregate_flag_stats(all_nodes)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


@mcp.tool()
def inspect_codon_rail(
    config: str = "human_igh",
    seed: int = 42,
    hook: str = "final",
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
) -> Dict[str, Any]:
    """Inspect the codon rail at a specific pipeline stage.

    The codon rail is the reading-frame overlay on the ASeq: each codon
    is a triplet of nucleotides with a translated amino acid. Shows all
    codons with their bases, amino acid, stop codon status, and segment.

    Args:
        config: Species/chain config name.
        seed: Random seed.
        hook: Hook point name (default: "final").
        min_mutation_rate: Enable SHM.
        max_mutation_rate: Maximum SHM rate.
    """
    err = _check_backend()
    if err:
        return {"error": err}

    try:
        from GenAIRR.utilities.mcp_helpers import validate_codon_rail_snapshot

        sim = _build_hooked_simulator(
            config, seed=seed,
            min_mutation_rate=min_mutation_rate,
            max_mutation_rate=max_mutation_rate,
        )
        sim.set_hooks([hook])
        sim.simulate_one()
        snaps = sim.get_snapshots()
        sim.set_hooks([])

        if not snaps:
            return {"error": f"No snapshot captured at hook '{hook}'"}

        codons = sim.get_snapshot_codon_rail(0)
        validation = validate_codon_rail_snapshot(codons)

        return {
            "hook": hook,
            "codons": codons,
            "validation": validation,
        }
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


@mcp.tool()
def diff_snapshots(
    config: str = "human_igh",
    seed: int = 42,
    hook_before: str = "post_assembly",
    hook_after: str = "post_mutation",
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    productive: bool = False,
    indel_prob: Optional[float] = None,
    corrupt_5prime: bool = False,
    corrupt_3prime: bool = False,
) -> Dict[str, Any]:
    """Compare the ASeq state at two pipeline stages.

    Shows exactly what changed between two hook points: which nodes were
    added, removed, or modified, and how flags changed. Perfect for
    questions like "what exactly did S5F mutation do to this sequence?"

    Args:
        config: Species/chain config name.
        seed: Random seed.
        hook_before: First hook point (default: "post_assembly").
        hook_after: Second hook point (default: "post_mutation").
        min_mutation_rate: Enable SHM.
        max_mutation_rate: Maximum SHM rate.
        productive: Only productive rearrangements.
        indel_prob: Enable indels.
        corrupt_5prime: Enable 5' corruption.
        corrupt_3prime: Enable 3' corruption.
    """
    err = _check_backend()
    if err:
        return {"error": err}

    try:
        from GenAIRR.utilities.mcp_helpers import (
            summarize_aseq, diff_snapshots as _diff,
        )

        sim = _build_hooked_simulator(
            config, seed=seed,
            min_mutation_rate=min_mutation_rate,
            max_mutation_rate=max_mutation_rate,
            productive=productive, indel_prob=indel_prob,
            corrupt_5prime=corrupt_5prime,
            corrupt_3prime=corrupt_3prime,
        )

        rec, snapshots, trace = sim.simulate_one_hooked([hook_before, hook_after])

        # Find the two snapshots
        snap_a = next((s for s in snapshots if s["hook"] == hook_before), None)
        snap_b = next((s for s in snapshots if s["hook"] == hook_after), None)

        if not snap_a:
            return {"error": f"Hook '{hook_before}' not captured (step may not be active)"}
        if not snap_b:
            return {"error": f"Hook '{hook_after}' not captured (step may not be active)"}

        diff_result = _diff(snap_a["nodes"], snap_b["nodes"])
        diff_result["before_summary"] = summarize_aseq(snap_a["nodes"])
        diff_result["after_summary"] = summarize_aseq(snap_b["nodes"])
        diff_result["before_hook"] = hook_before
        diff_result["after_hook"] = hook_after

        return diff_result
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


@mcp.tool()
def replay_with_trace(
    config: str = "human_igh",
    seed: int = 42,
    hooks: Optional[List[str]] = None,
    min_mutation_rate: Optional[float] = None,
    max_mutation_rate: Optional[float] = None,
    productive: bool = False,
    indel_prob: Optional[float] = None,
    include_nodes: bool = False,
) -> Dict[str, Any]:
    """Full pipeline replay with trace log and snapshot timeline.

    Simulates one sequence with trace logging enabled and captures
    snapshots at multiple pipeline stages. Returns the AIRR record,
    execution trace, and a compact timeline showing how the sequence
    evolved through each stage.

    Args:
        config: Species/chain config name.
        seed: Random seed.
        hooks: Hook points (default: post_assembly, post_functionality,
               post_mutation, final).
        min_mutation_rate: Enable SHM.
        max_mutation_rate: Maximum SHM rate.
        productive: Only productive rearrangements.
        indel_prob: Enable indels.
        include_nodes: Include full node lists in snapshots (verbose).
    """
    err = _check_backend()
    if err:
        return {"error": err}

    if hooks is None:
        hooks = ["post_assembly", "post_functionality", "post_mutation", "final"]

    try:
        from GenAIRR.utilities.mcp_helpers import (
            summarize_aseq, detect_node_anomalies, format_snapshot_timeline,
        )

        sim = _build_hooked_simulator(
            config, seed=seed,
            min_mutation_rate=min_mutation_rate,
            max_mutation_rate=max_mutation_rate,
            productive=productive, indel_prob=indel_prob,
        )

        rec, snapshots, trace = sim.simulate_one_hooked(hooks)

        processed = []
        for snap in snapshots:
            nodes = snap["nodes"]
            entry: Dict[str, Any] = {
                "hook": snap["hook"],
                "summary": summarize_aseq(nodes),
                "anomalies": detect_node_anomalies(nodes),
            }
            if include_nodes:
                entry["nodes"] = nodes
            processed.append(entry)

        timeline = format_snapshot_timeline(snapshots)

        return {
            "airr_record": rec,
            "trace": trace,
            "timeline": timeline,
            "snapshots": processed,
        }
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


# ── Entry point ─────────────────────────────────────────────────

def main():
    """Run the GenAIRR MCP server (STDIO transport)."""
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
