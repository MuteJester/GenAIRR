"""GenAIRR MCP server.

Exposes 14 tools targeted at an LLM agent doing immunology research. Every
tool returns a uniform envelope:

    {"ok": True,  "tool": <name>, "elapsed_ms": int, "result": {...}}
    {"ok": False, "tool": <name>, "elapsed_ms": int, "error": {...}}

The envelope is applied via the @envelope decorator from GenAIRR.mcp_errors.
Tool implementations return only the inner `result` dict and raise MCPError
for structured failures.

Run: python -m fastmcp run src/GenAIRR/mcp_server.py
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional

from fastmcp import FastMCP

import GenAIRR as ga
from GenAIRR._mcp_refdata import (
    find_allele,
    iter_alleles,
    locus_from_config_name,
    resolve_refdata,
)
from GenAIRR.mcp_errors import (
    ALLELE_NOT_FOUND,
    INVALID_PARAMETER,
    MALFORMED_RECORD,
    SEED_REPLAY_MISMATCH,
    MCPError,
    envelope,
)


mcp = FastMCP("GenAIRR")


# -- Helpers (private) ----------------------------------------------


def _split_config_name(name: str) -> Optional[str]:
    """Extract the species token from a config name like 'human_igh' -> 'Human'.

    Returns None when the name doesn't match the species_chain shape.
    """
    if "_" not in name:
        return None
    species_token = name.split("_", 1)[0]
    return species_token.capitalize()


def _all_config_aliases() -> List[str]:
    """Return the user-facing config aliases (lowercased friendly names).

    These are the names accepted by ``Experiment.on(...)`` -- short forms
    like 'human_igh' as well as long forms like 'human_igh_imgt'. The
    engine's ``ga.list_configs()`` exposes only the uppercase technical
    names (e.g. 'HUMAN_IGH_OGRDB'), which are not what an MCP agent should
    be passing to the other tools. We surface the aliases instead so the
    output of ``list_configs`` round-trips directly into every other tool.
    """
    from GenAIRR.experiment import _CONFIG_ALIASES

    return sorted(_CONFIG_ALIASES.keys())


# -- Discovery tools ------------------------------------------------


@mcp.tool()
@envelope("list_configs")
def list_configs(species_filter: Optional[str] = None) -> Dict[str, Any]:
    """Enumerate available species/chain configurations.

    Args:
        species_filter: case-insensitive species name to narrow results
            (e.g. "Human", "Mouse"). When None, returns every config.

    Returns:
        {
          "configs": List[str],            # e.g. ["human_igh", "mouse_tcrb", ...]
          "n_total": int,
          "by_species": Dict[str, List[str]],  # e.g. {"Human": ["human_igh", ...]}
        }
    """
    all_configs = _all_config_aliases()
    if species_filter:
        norm = species_filter.lower()
        all_configs = [c for c in all_configs if c.startswith(f"{norm}_")]

    by_species: Dict[str, List[str]] = defaultdict(list)
    for cfg in all_configs:
        species = _split_config_name(cfg)
        if species is not None:
            by_species[species].append(cfg)

    return {
        "configs": all_configs,
        "n_total": len(all_configs),
        "by_species": dict(sorted(by_species.items())),
    }


# -- config_info tool ----------------------------------------------


@mcp.tool()
@envelope("config_info")
def config_info(config: str) -> Dict[str, Any]:
    """Per-config gene/allele counts, locus, chain type, anchor coverage.

    Args:
        config: config name (e.g. "human_igh", "mouse_tcrb").

    Returns:
        {
          "config": str,
          "chain_type": str,              # e.g. "vdj", "vj"
          "locus": str,                   # e.g. "IGH"
          "species": str,
          "v_pool_size": int,
          "d_pool_size": int,             # 0 for VJ chains
          "j_pool_size": int,
          "anchor_coverage": Dict[str, float],  # {"v": 1.0, "j": 0.98}
        }

    Raises:
        MCPError(CONFIG_NOT_FOUND): when the config name isn't recognised.
    """
    refdata = resolve_refdata(config)
    species = _split_config_name(config) or "Unknown"
    locus = locus_from_config_name(config) or ""

    def anchor_rate(segment: str) -> float:
        total = 0
        anchored = 0
        for allele in iter_alleles(refdata, segment):
            total += 1
            if getattr(allele, "anchor", None) is not None:
                anchored += 1
        return round(anchored / total, 4) if total > 0 else 0.0

    return {
        "config": config,
        "chain_type": str(getattr(refdata, "chain_type", "")),
        "locus": locus,
        "species": species,
        "v_pool_size": refdata.v_pool_size(),
        "d_pool_size": refdata.d_pool_size() if refdata.has_d() else 0,
        "j_pool_size": refdata.j_pool_size(),
        "anchor_coverage": {
            "v": anchor_rate("v"),
            "j": anchor_rate("j"),
        },
    }


# -- list_alleles tool ---------------------------------------------


@mcp.tool()
@envelope("list_alleles")
def list_alleles(
    config: str,
    segment: str = "v",
    limit: int = 100,
) -> Dict[str, Any]:
    """Allele names from a config's V/D/J/C pool.

    Args:
        config: config name (e.g. "human_igh").
        segment: one of "v", "d", "j", "c".
        limit: cap on the number of names returned. Set high to enumerate all.

    Returns:
        {
          "config": str,
          "segment": str,
          "n_total": int,
          "allele_names": List[str],
          "truncated": bool,
        }

    Raises:
        MCPError(CONFIG_NOT_FOUND), MCPError(INVALID_PARAMETER).
    """
    if segment not in {"v", "d", "j", "c"}:
        raise MCPError(
            INVALID_PARAMETER,
            f"segment={segment!r} must be one of 'v', 'd', 'j', 'c'.",
        )

    refdata = resolve_refdata(config)
    all_names = sorted({allele.name for allele in iter_alleles(refdata, segment)})
    truncated = len(all_names) > limit
    return {
        "config": config,
        "segment": segment,
        "n_total": len(all_names),
        "allele_names": all_names[:limit],
        "truncated": truncated,
    }


# -- inspect_allele tool -------------------------------------------


@mcp.tool()
@envelope("inspect_allele")
def inspect_allele(config: str, allele_name: str) -> Dict[str, Any]:
    """One allele's full detail -- sequence, anchor, length, gene family.

    Args:
        config: config name.
        allele_name: e.g. "IGHVF1-G1*01". Searches V, D, J pools.

    Returns:
        {
          "config": str,
          "allele_name": str,
          "segment": str,                 # "v" | "d" | "j"
          "length": int,
          "anchor": Optional[int],
          "gene_family": Optional[str],
          "sequence": str,                # raw germline bases
        }

    Raises:
        MCPError(CONFIG_NOT_FOUND), MCPError(ALLELE_NOT_FOUND).
    """
    refdata = resolve_refdata(config)

    for segment in ("v", "d", "j"):
        allele = find_allele(refdata, segment, allele_name)
        if allele is None:
            continue
        raw_seq = allele.seq()
        sequence = (
            raw_seq.decode("ascii")
            if isinstance(raw_seq, (bytes, bytearray))
            else str(raw_seq)
        )
        return {
            "config": config,
            "allele_name": allele_name,
            "segment": segment,
            "length": len(sequence),
            "anchor": int(allele.anchor) if allele.anchor is not None else None,
            "gene_family": getattr(allele, "gene", None),
            "sequence": sequence,
        }

    raise MCPError(
        ALLELE_NOT_FOUND,
        f"Allele {allele_name!r} not found in any pool of config {config!r}.",
        hint="Call list_alleles to enumerate available allele names.",
    )


# -- Simulation helpers (private) -----------------------------------


def _build_experiment_from_params(
    *,
    config: str,
    mutation_model: Optional[str],
    mutation_count_min: Optional[int],
    mutation_count_max: Optional[int],
    five_prime_loss_max: Optional[int],
    three_prime_loss_max: Optional[int],
    pcr_error_count_max: Optional[int],
    indel_count_max: Optional[int],
    n_injection_count_max: Optional[int],
    quality_count_max: Optional[int],
    contaminant_prob: Optional[float],
    rev_comp_prob: Optional[float],
    n_clones: Optional[int],
    clone_size: Optional[int],
    v_alleles: Optional[List[str]],
    d_alleles: Optional[List[str]],
    j_alleles: Optional[List[str]],
) -> Any:
    """Translate flat MCP parameters into a v2.0.0 Experiment object.

    The clonal-fork policy (per the README's realistic-pipeline section):
      passes BEFORE with_clonal_structure apply to the parent rearrangement
      (shared by every sister sequence in the clone); passes AFTER apply
      per-descendant (independent SHM + per-descendant sequencing artefacts).
    So we order: recombine -> using -> with_clonal_structure -> mutate -> corrupt*.

    Raises:
        MCPError(CONFIG_NOT_FOUND): unknown config alias.
        MCPError(INVALID_PARAMETER): half-specified clonal pair or mutation
            model without any count bound.
    """
    # resolve_refdata raises MCPError(CONFIG_NOT_FOUND) on unknown aliases.
    # Probe the config first so the error envelope is consistent with every
    # other tool, even though we don't need the refdata object itself here.
    resolve_refdata(config)

    # Validate clonal pair (half-specified is a user error)
    if (n_clones is None) != (clone_size is None):
        raise MCPError(
            INVALID_PARAMETER,
            "n_clones and clone_size must be provided together (or both omitted).",
        )

    exp = ga.Experiment.on(config).recombine()

    # -- Allele lock (parent-level) -----------------------------------
    using_kwargs: Dict[str, Any] = {}
    if v_alleles:
        using_kwargs["v"] = list(v_alleles)
    if d_alleles:
        using_kwargs["d"] = list(d_alleles)
    if j_alleles:
        using_kwargs["j"] = list(j_alleles)
    if using_kwargs:
        exp = exp.restrict_alleles(**using_kwargs)

    # -- Clonal fork (between parent and per-descendant noise) --------
    if n_clones is not None:
        # The (n_clones is None) != (clone_size is None) guard above
        # ensures both are set together. The assert pins the invariant
        # for the type checker.
        assert clone_size is not None
        exp = exp.with_clonal_structure(n_clones=int(n_clones), size=int(clone_size))

    # -- Per-descendant SHM -------------------------------------------
    if (
        mutation_model is not None
        or mutation_count_min is not None
        or mutation_count_max is not None
    ):
        # Resolve `count` shape: int when only one bound given, tuple when both.
        if mutation_count_min is not None and mutation_count_max is not None:
            mut_count: Any = (int(mutation_count_min), int(mutation_count_max))
        elif mutation_count_max is not None:
            mut_count = (0, int(mutation_count_max))
        elif mutation_count_min is not None:
            mut_count = int(mutation_count_min)
        else:
            # mutation_model set without any count -- engine requires count.
            raise MCPError(
                INVALID_PARAMETER,
                "mutation_model set without mutation_count_min or mutation_count_max.",
            )
        exp = exp.mutate(model=(mutation_model or "s5f"), count=mut_count)

    # -- Per-descendant sequencing artefacts --------------------------
    if five_prime_loss_max is not None:
        exp = exp.primer_trim_5prime(length=(0, int(five_prime_loss_max)))
    if three_prime_loss_max is not None:
        exp = exp.primer_trim_3prime(length=(0, int(three_prime_loss_max)))
    if indel_count_max is not None:
        exp = exp.library_indels(
            count=(0, int(indel_count_max)),
            insertion_prob=0.5,
        )
    if pcr_error_count_max is not None:
        exp = exp.pcr_amplify(count=(0, int(pcr_error_count_max)))
    if n_injection_count_max is not None:
        exp = exp.mask_low_quality(count=(0, int(n_injection_count_max)))
    if quality_count_max is not None:
        exp = exp.sequencing_errors(count=(0, int(quality_count_max)))
    if contaminant_prob is not None:
        exp = exp.contaminate(prob=float(contaminant_prob))
    if rev_comp_prob is not None:
        exp = exp.random_strand_orientation(prob=float(rev_comp_prob))

    # Note: productive_only is applied by the caller as .productive_only(),
    # not here, so this helper can be re-used by both productive and
    # non-productive simulation entry points.
    return exp


# -- simulate_repertoire tool --------------------------------------


def _simulate_repertoire_impl(
    *,
    config: str,
    n: int = 100,
    seed: int = 0,
    productive_only: bool = False,
    mutation_model: Optional[str] = None,
    mutation_count_min: Optional[int] = None,
    mutation_count_max: Optional[int] = None,
    five_prime_loss_max: Optional[int] = None,
    three_prime_loss_max: Optional[int] = None,
    pcr_error_count_max: Optional[int] = None,
    indel_count_max: Optional[int] = None,
    n_injection_count_max: Optional[int] = None,
    quality_count_max: Optional[int] = None,
    contaminant_prob: Optional[float] = None,
    rev_comp_prob: Optional[float] = None,
    n_clones: Optional[int] = None,
    clone_size: Optional[int] = None,
    v_alleles: Optional[List[str]] = None,
    d_alleles: Optional[List[str]] = None,
    j_alleles: Optional[List[str]] = None,
    return_records: bool = False,
    return_records_limit: int = 100,
    summary_top_n: int = 10,
) -> Dict[str, Any]:
    """Inner implementation of simulate_repertoire -- returns the raw result dict.

    Shared by simulate_repertoire, simulate_preset, and simulate_allele so that
    each tool's @envelope decorator fires exactly once per MCP call. fastmcp 3.x
    registers @mcp.tool() functions as plain Python functions, so we cannot
    bypass the envelope by introspecting `.fn.__wrapped__` -- a shared private
    helper is the clean way to delegate.

    See simulate_repertoire's docstring for parameter and return-shape semantics.
    """
    from GenAIRR._mcp_summary import compute_repertoire_summary

    exp = _build_experiment_from_params(
        config=config,
        mutation_model=mutation_model,
        mutation_count_min=mutation_count_min,
        mutation_count_max=mutation_count_max,
        five_prime_loss_max=five_prime_loss_max,
        three_prime_loss_max=three_prime_loss_max,
        pcr_error_count_max=pcr_error_count_max,
        indel_count_max=indel_count_max,
        n_injection_count_max=n_injection_count_max,
        quality_count_max=quality_count_max,
        contaminant_prob=contaminant_prob,
        rev_comp_prob=rev_comp_prob,
        n_clones=n_clones,
        clone_size=clone_size,
        v_alleles=v_alleles,
        d_alleles=d_alleles,
        j_alleles=j_alleles,
    )

    if productive_only:
        exp = exp.productive_only()

    # When clonal: n comes from n_clones x clone_size -- don't pass `n`.
    if n_clones is not None:
        result = exp.run_records(seed=seed)
    else:
        result = exp.run_records(n=n, seed=seed)

    records = list(result.records)
    summary = compute_repertoire_summary(records, summary_top_n=summary_top_n)

    out: Dict[str, Any] = {
        "config": config,
        "n_records": len(records),
        "seed": seed,
        "productive_only": productive_only,
        **summary,
    }

    if return_records:
        out["records"] = records[:return_records_limit]
        out["records_truncated"] = len(records) > return_records_limit

    return out


@mcp.tool()
@envelope("simulate_repertoire")
def simulate_repertoire(
    config: str,
    n: int = 100,
    seed: int = 0,
    productive_only: bool = False,
    mutation_model: Optional[str] = None,
    mutation_count_min: Optional[int] = None,
    mutation_count_max: Optional[int] = None,
    five_prime_loss_max: Optional[int] = None,
    three_prime_loss_max: Optional[int] = None,
    pcr_error_count_max: Optional[int] = None,
    indel_count_max: Optional[int] = None,
    n_injection_count_max: Optional[int] = None,
    quality_count_max: Optional[int] = None,
    contaminant_prob: Optional[float] = None,
    rev_comp_prob: Optional[float] = None,
    n_clones: Optional[int] = None,
    clone_size: Optional[int] = None,
    v_alleles: Optional[List[str]] = None,
    d_alleles: Optional[List[str]] = None,
    j_alleles: Optional[List[str]] = None,
    return_records: bool = False,
    return_records_limit: int = 100,
    summary_top_n: int = 10,
) -> Dict[str, Any]:
    """Headline tool. Generate a repertoire and return a structured summary.

    Translates 22 flat parameters into a v2.0.0 Experiment via the private
    _build_experiment_from_params helper, runs the engine, and returns the
    summary computed by compute_repertoire_summary. AIRR records themselves
    are opt-in via return_records=True and capped at return_records_limit.

    See spec section 2 for parameter and return-shape semantics.

    Args:
        config: species/chain config alias (e.g. "human_igh").
        n: number of records to generate. Ignored when n_clones is set
            (total = n_clones x clone_size).
        seed: deterministic RNG seed.
        productive_only: when True, attach the productive-sequence
            contract bundle via :meth:`Experiment.productive_only` so
            every returned record satisfies it by construction.
        mutation_model: SHM model (e.g. "s5f"); requires a count bound.
        mutation_count_min / mutation_count_max: SHM count bounds.
        five_prime_loss_max / three_prime_loss_max: max bases trimmed
            from each end.
        pcr_error_count_max / indel_count_max / n_injection_count_max /
            quality_count_max: per-record sequencing artefact caps.
        contaminant_prob / rev_comp_prob: probability knobs.
        n_clones / clone_size: clonal lineage fork (both or neither).
        v_alleles / d_alleles / j_alleles: optional allele locks.
        return_records: include the actual AIRR record dicts in the
            response (default False -- agents see only the summary).
        return_records_limit: cap on records returned when opted in.
        summary_top_n: cap on each top-N usage table in the summary.

    Returns:
        {
          "config": str,
          "n_records": int,
          "seed": int,
          "productive_only": bool,
          ...summary fields...,
          "records": Optional[List[Dict]],
          "records_truncated": Optional[bool],
        }

    Raises:
        MCPError(CONFIG_NOT_FOUND): when the config alias is unknown.
        MCPError(INVALID_PARAMETER): when the n_clones / clone_size pair is
            half-specified, or mutation_model is set without any count bound.
    """
    return _simulate_repertoire_impl(
        config=config,
        n=n,
        seed=seed,
        productive_only=productive_only,
        mutation_model=mutation_model,
        mutation_count_min=mutation_count_min,
        mutation_count_max=mutation_count_max,
        five_prime_loss_max=five_prime_loss_max,
        three_prime_loss_max=three_prime_loss_max,
        pcr_error_count_max=pcr_error_count_max,
        indel_count_max=indel_count_max,
        n_injection_count_max=n_injection_count_max,
        quality_count_max=quality_count_max,
        contaminant_prob=contaminant_prob,
        rev_comp_prob=rev_comp_prob,
        n_clones=n_clones,
        clone_size=clone_size,
        v_alleles=v_alleles,
        d_alleles=d_alleles,
        j_alleles=j_alleles,
        return_records=return_records,
        return_records_limit=return_records_limit,
        summary_top_n=summary_top_n,
    )


# -- simulate_preset tool ------------------------------------------


@mcp.tool()
@envelope("simulate_preset")
def simulate_preset(
    preset: str,
    n: Optional[int] = None,
    seed: int = 0,
    return_records: bool = False,
    return_records_limit: int = 100,
) -> Dict[str, Any]:
    """Pre-baked simulation scenarios by name.

    See spec section 3 for the 5 preset definitions:
    naive_b_cell, memory_b_cell_shm, tcr_beta_pool,
    low_quality_sequencing, clonal_expansion.

    Args:
        preset: preset name.
        n: override the preset's default n. Ignored for clonal presets
            (where total = n_clones x clone_size).
        seed: deterministic seed.
        return_records / return_records_limit: passed through to the
            underlying simulate_repertoire implementation.

    Raises:
        MCPError(INVALID_PRESET): when the name is unknown.
    """
    from GenAIRR._mcp_presets import resolve_preset

    params = resolve_preset(preset)  # raises INVALID_PRESET
    # Caller's n overrides the preset's n unless the preset is clonal
    # (where total = n_clones x clone_size).
    if n is not None and "n_clones" not in params:
        params["n"] = int(n)
    params["seed"] = seed
    params["return_records"] = return_records
    params["return_records_limit"] = return_records_limit

    # Delegate to the shared inner implementation so the @envelope decorator
    # fires only once (on simulate_preset itself). fastmcp 3.x registers
    # tools as plain functions, so we cannot use simulate_repertoire.fn.
    return _simulate_repertoire_impl(**params)


# -- simulate_allele tool ------------------------------------------


@mcp.tool()
@envelope("simulate_allele")
def simulate_allele(
    config: str,
    v_allele: str,
    j_allele: str,
    d_allele: Optional[str] = None,
    n: int = 50,
    seed: int = 0,
    mutation_count_max: int = 10,
    productive_only: bool = True,
) -> Dict[str, Any]:
    """Focused simulation under a fixed V/D/J allele lock.

    Internally delegates to the shared simulate_repertoire implementation
    with v_alleles=[v_allele] / j_alleles=[j_allele] (and d_alleles=[d_allele]
    when supplied) plus a moderate SHM range. Returns the same summary shape
    -- useful for "what does this exact allele combination produce?" questions.

    Args:
        config: species/chain config alias (e.g. "human_igh").
        v_allele / j_allele: allele names to lock (e.g. "IGHVF1-G1*01").
        d_allele: optional D allele lock for VDJ chains.
        n: number of records to generate.
        seed: deterministic RNG seed.
        mutation_count_max: per-record SHM upper bound. Set to 0 to disable
            SHM entirely (so v_call / j_call stay exactly == the locked truth).
        productive_only: default True -- focused queries usually want clean
            productive records.
    """
    return _simulate_repertoire_impl(
        config=config,
        n=n,
        seed=seed,
        productive_only=productive_only,
        mutation_model="s5f" if mutation_count_max > 0 else None,
        mutation_count_min=0 if mutation_count_max > 0 else None,
        mutation_count_max=mutation_count_max if mutation_count_max > 0 else None,
        v_alleles=[v_allele],
        d_alleles=[d_allele] if d_allele else None,
        j_alleles=[j_allele],
    )


# -- validate_records tool -----------------------------------------


@mcp.tool()
@envelope("validate_records")
def validate_records(
    records: List[Dict[str, Any]],
    config: Optional[str] = None,
) -> Dict[str, Any]:
    """Per-record AIRR consistency checks.

    See spec section 3 + GenAIRR/_mcp_validators.py for the 9 checks
    (5 schema-only, 4 refdata-driven). Returns per-record issue lists
    plus aggregate counts.

    Args:
        records: list of AIRR record dicts.
        config: when provided, enables refdata-driven checks (checks 6-9).

    Raises:
        MCPError(CONFIG_NOT_FOUND): when config is supplied but unknown.
    """
    from GenAIRR._mcp_validators import validate_one_record

    refdata = resolve_refdata(config) if config is not None else None

    per_record: List[Dict[str, Any]] = []
    issue_counts: Dict[str, int] = defaultdict(int)
    n_valid = 0
    for i, rec in enumerate(records):
        if not isinstance(rec, dict):
            issues = [f"record at index {i} is not a dict"]
        else:
            issues = validate_one_record(rec, refdata=refdata, config_name=config)

        sequence_id = rec.get("sequence_id") if isinstance(rec, dict) else None
        per_record.append({
            "sequence_id": sequence_id if sequence_id else f"<index {i}>",
            "valid": len(issues) == 0,
            "issues": issues,
        })
        if not issues:
            n_valid += 1
        for issue in issues:
            issue_counts[issue] += 1

    return {
        "n_records": len(records),
        "n_valid": n_valid,
        "n_invalid": len(records) - n_valid,
        "per_record": per_record,
        "issue_counts": dict(sorted(issue_counts.items(), key=lambda kv: -kv[1])),
    }


# -- align_to_germline tool ----------------------------------------


def _align_one_segment(
    sequence: str,
    seq_start: Optional[int],
    seq_end: Optional[int],
    germline: str,
    ref_start: Optional[int],
    ref_end: Optional[int],
) -> Dict[str, Any]:
    """Compare two equal-length string slices base by base.

    The engine has already emitted matching start/end spans (CIGAR M+I-aware
    on the sequence side, M+D-aware on the germline side); here we only
    compute match/mismatch on the M positions, ignoring I and D ops by
    simple slicing.
    """
    if seq_start is None or seq_end is None or ref_start is None or ref_end is None:
        return {
            "identity": None,
            "match_count": 0,
            "mismatch_count": 0,
            "mismatch_positions": [],
        }

    obs = sequence[seq_start:seq_end].upper()
    ref = germline[ref_start:ref_end].upper()
    n = min(len(obs), len(ref))

    matches = 0
    mismatches: List[int] = []
    for i in range(n):
        if obs[i] == ref[i]:
            matches += 1
        else:
            mismatches.append(seq_start + i)
    return {
        "match_count": matches,
        "mismatch_count": len(mismatches),
        "identity": round(matches / n, 4) if n > 0 else None,
        "mismatch_positions": mismatches,
    }


def _germline_seq_for_call(
    refdata: Any, segment: str, call: str, config: str
) -> Optional[str]:
    """Resolve the first allele name in a comma-separated tie set to its
    germline string. Returns None when the call is empty or the segment's
    pool is missing (e.g. d_call on a VJ chain). Raises ALLELE_NOT_FOUND
    when the call is non-empty but doesn't match any allele in the pool.
    """
    if not call:
        return None
    first = str(call).split(",", 1)[0].strip()
    if not first:
        return None
    allele = find_allele(refdata, segment, first)
    if allele is None:
        # Distinguish "missing pool" from "missing allele": iter_alleles
        # is empty for missing pools, so a None result there means the
        # segment isn't populated for this config -- treat as no-germline
        # rather than as an error.
        if not any(True for _ in iter_alleles(refdata, segment)):
            return None
        raise MCPError(
            ALLELE_NOT_FOUND,
            f"Allele {first!r} (from {segment}_call) not found in config {config!r}.",
            hint="Call list_alleles to enumerate available allele names.",
        )
    raw_seq = allele.seq()
    return (
        raw_seq.decode("ascii")
        if isinstance(raw_seq, (bytes, bytearray))
        else str(raw_seq)
    )


@mcp.tool()
@envelope("align_to_germline")
def align_to_germline(config: str, record: Dict[str, Any]) -> Dict[str, Any]:
    """Per-segment alignment of a sequence against its claimed germline.

    Compares record["sequence"][v_sequence_start:v_sequence_end] against
    the V allele named in record["v_call"] at [v_germline_start:v_germline_end].
    Returns per-segment identity + mismatch positions for V/D/J.

    Args:
        config: config name (used to resolve allele names to germline).
        record: an AIRR record dict.

    Raises:
        MCPError(CONFIG_NOT_FOUND): unknown config alias.
        MCPError(MALFORMED_RECORD): record missing the sequence field.
        MCPError(ALLELE_NOT_FOUND): a claimed allele isn't in the pool.
    """
    if "sequence" not in record:
        raise MCPError(MALFORMED_RECORD, "record missing required field: sequence")

    refdata = resolve_refdata(config)
    sequence = record.get("sequence") or ""

    out: Dict[str, Any] = {}
    for segment in ("v", "d", "j"):
        call = record.get(f"{segment}_call") or ""
        germline = _germline_seq_for_call(refdata, segment, call, config)
        claimed_call = str(call).split(",", 1)[0].strip() if call else ""
        if germline is None:
            out[f"{segment}_alignment"] = {
                "claimed_call": claimed_call,
                "identity": None,
                "match_count": 0,
                "mismatch_count": 0,
                "mismatch_positions": [],
            }
            continue
        align = _align_one_segment(
            sequence=sequence,
            seq_start=record.get(f"{segment}_sequence_start"),
            seq_end=record.get(f"{segment}_sequence_end"),
            germline=germline,
            ref_start=record.get(f"{segment}_germline_start"),
            ref_end=record.get(f"{segment}_germline_end"),
        )
        align["claimed_call"] = claimed_call
        out[f"{segment}_alignment"] = align

    return out


# -- score_allele_calls tool ---------------------------------------


@mcp.tool()
@envelope("score_allele_calls")
def score_allele_calls(
    config: str,
    record: Dict[str, Any],
    top_k: int = 5,
) -> Dict[str, Any]:
    """Re-score the claimed v_call / d_call / j_call against the full pool.

    For each segment, compute match count between
    record['sequence'][seg_sequence_start:seg_sequence_end] and EVERY allele
    in the pool at the same projected offsets; return the top-K alleles
    with their match counts + identity, plus an is_claimed flag for the
    one(s) named in record[f"{segment}_call"]. Ties broken by allele name
    (alphabetical) for determinism.

    Args:
        config: config name.
        record: AIRR record dict.
        top_k: how many top-ranked alleles to return per segment.

    Raises:
        MCPError(CONFIG_NOT_FOUND): unknown config alias.
        MCPError(MALFORMED_RECORD): record missing the sequence field.
    """
    if "sequence" not in record:
        raise MCPError(MALFORMED_RECORD, "record missing required field: sequence")

    refdata = resolve_refdata(config)
    sequence = record.get("sequence") or ""

    def _score_segment(segment: str) -> List[Dict[str, Any]]:
        seq_start = record.get(f"{segment}_sequence_start")
        seq_end = record.get(f"{segment}_sequence_end")
        ref_start = record.get(f"{segment}_germline_start")
        ref_end = record.get(f"{segment}_germline_end")
        if any(v is None for v in (seq_start, seq_end, ref_start, ref_end)):
            return []
        obs = sequence[seq_start:seq_end].upper()

        claimed_names = {
            n.strip()
            for n in str(record.get(f"{segment}_call") or "").split(",")
            if n.strip()
        }

        scores: List[Dict[str, Any]] = []
        for allele in iter_alleles(refdata, segment):
            raw = allele.seq()
            ref_str = (
                raw.decode("ascii")
                if isinstance(raw, (bytes, bytearray))
                else str(raw)
            ).upper()
            ref_slice = ref_str[ref_start:ref_end]
            n = min(len(obs), len(ref_slice))
            if n == 0:
                continue
            matches = sum(1 for i in range(n) if obs[i] == ref_slice[i])
            scores.append({
                "allele": allele.name,
                "score": matches,
                "identity": round(matches / n, 4),
                "is_claimed": allele.name in claimed_names,
            })

        scores.sort(key=lambda e: (-e["score"], e["allele"]))
        return scores[:top_k]

    rankings = {f"{seg}_score_ranking": _score_segment(seg) for seg in ("v", "d", "j")}
    claimed_is_best = {
        seg: bool(rankings[f"{seg}_score_ranking"]
                  and rankings[f"{seg}_score_ranking"][0]["is_claimed"])
        for seg in ("v", "d", "j")
    }

    return {**rankings, "claimed_call_is_best": claimed_is_best}


# -- SHM analysis ----------------------------------------------------

# WRC and RGYW are SHM hotspot motifs from Yaari et al. We check both the
# forward (WRC) and reverse-complement (GYW) plus the canonical 5-mer
# motif (RGYW) seen by an aligner. Position refers to the C (forward) or
# G (reverse) base -- that's the one AID targets.
_HOTSPOT_NUCS = {"A", "T"}  # W = A or T

_PURINES = {"A", "G"}
_PYRIMIDINES = {"C", "T"}


def _is_transition(from_base: str, to_base: str) -> bool:
    f, t = from_base.upper(), to_base.upper()
    return (f in _PURINES and t in _PURINES) or (f in _PYRIMIDINES and t in _PYRIMIDINES)


def _is_hotspot_context(sequence: str, pos: int) -> bool:
    """Return True when the base at `pos` sits in a WRC (A/T-A/G-C) or
    GYW (G-C/T-A/T) trinucleotide context. Hotspot definition follows
    AID's preferred motifs; this is a coarse proxy used to flag mutations
    that landed in canonical SHM-favoured positions."""
    if pos < 1 or pos >= len(sequence) - 1:
        return False
    a = sequence[pos - 1].upper()
    b = sequence[pos].upper()
    c = sequence[pos + 1].upper()
    # WRC: a in W, b in R (A or G), c == C
    if a in _HOTSPOT_NUCS and b in {"A", "G"} and c == "C":
        return True
    # GYW: a == G, b in Y (C or T), c in W
    if a == "G" and b in {"C", "T"} and c in _HOTSPOT_NUCS:
        return True
    return False


@mcp.tool()
@envelope("analyze_mutations")
def analyze_mutations(config: str, record: Dict[str, Any]) -> Dict[str, Any]:
    """SHM pattern analysis on a single record.

    Compares record["sequence"] against the V germline (claimed v_call)
    over [v_sequence_start, v_sequence_end] and produces:

      - n_mutations
      - transition_count, transversion_count, transition_transversion_ratio
      - hotspot_targeted_count, hotspot_targeting_rate
      - per_mutation: List[{position, from, to, is_hotspot}]

    Args:
        config: config name.
        record: AIRR record dict.

    Raises:
        MCPError(CONFIG_NOT_FOUND), MCPError(MALFORMED_RECORD),
        MCPError(ALLELE_NOT_FOUND).
    """
    if "sequence" not in record:
        raise MCPError(MALFORMED_RECORD, "record missing required field: sequence")

    refdata = resolve_refdata(config)

    sequence = (record.get("sequence") or "").upper()
    v_call = record.get("v_call") or ""
    first_v = v_call.split(",", 1)[0].strip() if v_call else ""

    seq_start = record.get("v_sequence_start")
    seq_end = record.get("v_sequence_end")
    ref_start = record.get("v_germline_start")
    ref_end = record.get("v_germline_end")

    germline = None
    if first_v:
        allele = find_allele(refdata, "v", first_v)
        if allele is not None:
            raw = allele.seq()
            germline = (
                raw.decode("ascii")
                if isinstance(raw, (bytes, bytearray))
                else str(raw)
            ).upper()

    per_mutation: List[Dict[str, Any]] = []
    transitions = 0
    transversions = 0
    hotspots = 0

    if (
        germline
        and seq_start is not None
        and seq_end is not None
        and ref_start is not None
        and ref_end is not None
    ):
        obs = sequence[seq_start:seq_end]
        ref = germline[ref_start:ref_end]
        n = min(len(obs), len(ref))
        for i in range(n):
            if obs[i] != ref[i] and obs[i] != "N" and ref[i] != "N":
                from_b, to_b = ref[i], obs[i]
                is_hotspot = _is_hotspot_context(sequence, seq_start + i)
                per_mutation.append({
                    "position": seq_start + i,
                    "from": from_b,
                    "to": to_b,
                    "is_hotspot": is_hotspot,
                })
                if _is_transition(from_b, to_b):
                    transitions += 1
                else:
                    transversions += 1
                if is_hotspot:
                    hotspots += 1

    n_mutations = len(per_mutation)
    return {
        "n_mutations": n_mutations,
        "transition_count": transitions,
        "transversion_count": transversions,
        "transition_transversion_ratio": (
            round(transitions / transversions, 2) if transversions > 0 else None
        ),
        "hotspot_targeted_count": hotspots,
        "hotspot_targeting_rate": (
            round(hotspots / n_mutations, 3) if n_mutations > 0 else 0.0
        ),
        "per_mutation": per_mutation,
    }


# -- Region classification ------------------------------------------


@mcp.tool()
@envelope("classify_regions")
def classify_regions(config: str, record: Dict[str, Any]) -> Dict[str, Any]:
    """IMGT region breakdown of a record.

    In v2.0.0 the engine's Allele type carries no `imgt_regions` attribute
    (only anchor/gene/name/segment/seq), so per-allele FWR1/CDR1/FWR2/CDR2/FWR3
    boundaries are unavailable to project into the record's coordinate space.
    We degrade gracefully and emit only the regions derivable from the AIRR
    coordinates the engine does emit:

      - CDR3: junction_start + 3 .. junction_end - 3 (IMGT convention:
        junction minus the flanking anchor codons).
      - FWR4: junction_end .. j_sequence_end (J anchor sits at junction_end-3;
        FWR4 extends from there to j_sequence_end).

    If/when the Allele type grows an imgt_regions attribute downstream,
    this tool should be extended to project those onto FWR1/CDR1/FWR2/CDR2/FWR3.

    Args:
        config: config name.
        record: AIRR record dict.

    Returns:
        {
          "sequence_length": int,
          "regions": Dict[str, [int, int]],   # 0-based half-open
          "region_lengths": Dict[str, int],
          "junction_in_cdr3": bool,
        }

    Raises:
        MCPError(MALFORMED_RECORD), MCPError(CONFIG_NOT_FOUND).
    """
    if "sequence" not in record:
        raise MCPError(MALFORMED_RECORD, "record missing required field: sequence")

    # Resolve refdata for validation that the config exists, even if we don't
    # currently project per-allele boundaries (see docstring).
    refdata = resolve_refdata(config)

    sequence = record.get("sequence") or ""
    seq_len = len(sequence)
    v_call = (record.get("v_call") or "").split(",", 1)[0].strip()
    v_seq_start = record.get("v_sequence_start")
    j_seq_end = record.get("j_sequence_end")

    regions: Dict[str, List[int]] = {}

    # Look up V allele -- preserved as a hook for future imgt_regions support.
    # Today the v2.0.0 Allele type has no `imgt_regions` attribute, so we
    # cannot project FWR1/CDR1/FWR2/CDR2/FWR3 from allele coords onto sequence
    # coords. We still resolve the allele to surface ALLELE_NOT_FOUND-style
    # mismatches early if v_call references something not in the pool, but
    # we don't currently raise -- find_allele returning None is tolerated.
    if v_call and v_seq_start is not None:
        v_allele = find_allele(refdata, "v", v_call)
        if v_allele is not None:
            imgt = getattr(v_allele, "imgt_regions", None)
            if isinstance(imgt, dict):
                v_germline_start = record.get("v_germline_start") or 0
                for name, span in imgt.items():
                    try:
                        a, b = int(span[0]), int(span[1])
                    except (TypeError, IndexError, ValueError):
                        continue
                    proj_start = v_seq_start + max(0, a - v_germline_start)
                    proj_end = v_seq_start + max(0, b - v_germline_start)
                    proj_start = max(0, min(seq_len, proj_start))
                    proj_end = max(0, min(seq_len, proj_end))
                    regions[name] = [proj_start, proj_end]

    # CDR3: derived from junction_start/junction_end (excluding the
    # flanking anchor codons -- AIRR junction includes them, IMGT CDR3
    # is junction minus first 3 and last 3).
    junction_start = record.get("junction_start")
    junction_end = record.get("junction_end")
    cdr3_in_junction = False
    if junction_start is not None and junction_end is not None:
        cdr3_start = junction_start + 3
        cdr3_end = junction_end - 3
        if cdr3_end >= cdr3_start:
            regions["CDR3"] = [cdr3_start, cdr3_end]
            cdr3_in_junction = True

    # FWR4: from the J anchor + 3 to sequence end (the J anchor sits at
    # junction_end - 3; FWR4 extends from there to j_sequence_end).
    if j_seq_end is not None and junction_end is not None and j_seq_end > junction_end:
        regions["FWR4"] = [junction_end, j_seq_end]

    region_lengths = {name: span[1] - span[0] for name, span in regions.items()}

    return {
        "sequence_length": seq_len,
        "regions": regions,
        "region_lengths": region_lengths,
        "junction_in_cdr3": cdr3_in_junction,
    }


# -- Dataset summarisation -------------------------------------------


@mcp.tool()
@envelope("summarize_dataset")
def summarize_dataset(
    records: List[Dict[str, Any]],
    summary_top_n: int = 10,
) -> Dict[str, Any]:
    """Compute the same summary shape as simulate_repertoire on user-supplied records.

    Lets the agent compare external AIRR data against simulated data
    apples-to-apples: simulate_repertoire returns {v_usage_top, ...} and
    so does this tool.

    Args:
        records: list of AIRR record dicts.
        summary_top_n: how many top V/D/J calls to surface.
    """
    from GenAIRR._mcp_summary import compute_repertoire_summary

    return {
        "n_records": len(records),
        **compute_repertoire_summary(records, summary_top_n=summary_top_n),
    }


# -- Reproducibility: replay_seed -----------------------------------


_TRACE_HEADLINE_ADDRESSES = (
    "sample_allele.v",
    "sample_allele.d",
    "sample_allele.j",
    "trim.v_5", "trim.v_3",
    "trim.d_5", "trim.d_3",
    "trim.j_5", "trim.j_3",
    "np.np1.length",
    "np.np2.length",
)


def _summarize_trace(outcome: Any) -> Dict[str, Any]:
    """Pull the headline draws out of a full Outcome.trace() -- allele picks,
    trim lengths, NP lengths -- and stringify the values for JSON safety.
    The agent gets the meaningful per-record randomness without the full
    byte-level dump."""
    trace = outcome.trace()
    key_addresses: Dict[str, str] = {}
    for addr in _TRACE_HEADLINE_ADDRESSES:
        rec = trace.find(addr)
        if rec is None:
            continue
        # ChoiceRecord.value is a discriminated value (Int/Base/etc.);
        # str-cast is the safest JSON projection.
        key_addresses[addr] = str(rec.value)
    return {
        "n_passes": len(outcome.pass_names()),
        "pass_names": list(outcome.pass_names()),
        "key_addresses": key_addresses,
    }


@mcp.tool()
@envelope("replay_seed")
def replay_seed(
    config: str,
    seed: int,
    params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Reproduce a specific simulation and return the record + trace summary.

    Args:
        config: config name.
        seed: deterministic seed.
        params: optional subset of simulate_repertoire params (productive_only,
            mutation_model, mutation_count_min/max, corruption knobs, etc.).
            Defaults to no constraints / no SHM / no corruption.

    Returns:
        {
          "config": str, "seed": int, "params_used": Dict[str, Any],
          "record": Dict[str, Any],   # full AIRR fields
          "trace_summary": {
            "n_passes": int,
            "pass_names": List[str],
            "key_addresses": Dict[str, str],   # allele picks + trim + np lengths
          },
        }

    Raises:
        MCPError(CONFIG_NOT_FOUND).
    """
    params = dict(params or {})
    # Reuse simulate_repertoire's parameter translation but ALSO need the
    # Outcome for the trace summary. simulate_repertoire returns records
    # not outcomes. So we rebuild the Experiment + call .run() (Outcome path).
    productive_only = bool(params.pop("productive_only", False))

    # Translate to an Experiment using the same helper as simulate_repertoire.
    # _build_experiment_from_params no longer accepts productive_only -- that
    # is applied by attaching .productive_only() to the returned Experiment.
    build_kwargs = {
        "config": config,
        "mutation_model": params.get("mutation_model"),
        "mutation_count_min": params.get("mutation_count_min"),
        "mutation_count_max": params.get("mutation_count_max"),
        "five_prime_loss_max": params.get("five_prime_loss_max"),
        "three_prime_loss_max": params.get("three_prime_loss_max"),
        "pcr_error_count_max": params.get("pcr_error_count_max"),
        "indel_count_max": params.get("indel_count_max"),
        "n_injection_count_max": params.get("n_injection_count_max"),
        "quality_count_max": params.get("quality_count_max"),
        "contaminant_prob": params.get("contaminant_prob"),
        "rev_comp_prob": params.get("rev_comp_prob"),
        "n_clones": params.get("n_clones"),
        "clone_size": params.get("clone_size"),
        "v_alleles": params.get("v_alleles"),
        "d_alleles": params.get("d_alleles"),
        "j_alleles": params.get("j_alleles"),
    }
    exp = _build_experiment_from_params(**build_kwargs)
    if productive_only:
        exp = exp.productive_only()

    if build_kwargs["n_clones"] is not None:
        outcomes = exp.run(seed=seed)
    else:
        outcomes = exp.run(n=1, seed=seed)
    if not outcomes:
        raise MCPError(
            SEED_REPLAY_MISMATCH,
            f"Replay with seed={seed} produced no outcomes.",
        )

    # Build the AIRR record from the first outcome (consistent with what
    # simulate_repertoire would return at index 0). Use resolve_refdata
    # rather than exp.refdata for consistency with every other tool that
    # touches the refdata surface.
    from GenAIRR._airr_record import outcome_to_airr_record
    refdata = resolve_refdata(config)
    record = outcome_to_airr_record(outcomes[0], refdata, sequence_id="seq0")

    return {
        "config": config,
        "seed": seed,
        "params_used": {**build_kwargs, "productive_only": productive_only},
        "record": record,
        "trace_summary": _summarize_trace(outcomes[0]),
    }


# -- Entry point ---------------------------------------------------


def main() -> None:
    """Run the GenAIRR MCP server over the STDIO transport.

    Called by `python -m GenAIRR.mcp_server`, which is how the Claude
    Code / Claude Desktop / Cursor MCP client launches this process
    per the .mcp.json config in the repo root.
    """
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
