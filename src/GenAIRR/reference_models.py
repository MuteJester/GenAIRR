"""Reference empirical models — the cartridge's defaults plane.

A reference cartridge ships with empirical distributions the engine
uses as defaults when the user doesn't override them — NP lengths,
exonuclease trim amounts, etc. Historically these lived inside
``DataConfig`` as nested ``{family: {gene: {value: prob}}}`` dicts,
extracted into the engine's flat ``[(value, weight), ...]`` shape on
demand by ``GenAIRR._dataconfig_extract``.

This module adds an explicit, typed cartridge-side authoring surface.
A user can attach a :class:`ReferenceEmpiricalModels` to a DataConfig
and the loader will prefer it over the legacy extraction path. The
shape mirrors what the engine ultimately consumes — flat
``(value, weight)`` lists — so there's no marginalisation magic in
the middle.

Scope of v1
~~~~~~~~~~~

- NP1 / NP2 length distributions.
- V_3 / D_5 / D_3 / J_5 trim distributions.

Validation is shape-only here: weights finite and positive, values
non-negative ints, trim/NP keys recognised, D trims rejected on VJ
chains. The Rust ``RefDataConfig`` does not yet carry empirical
models; this plane is pure Python for now (the engine consumes
already-lowered pass distributions). Later slices may surface a
normalised summary into the Rust cartridge for trace identity, but
that's deliberately out of scope here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# Canonical key sets for each model dict. Keep small + explicit so
# users can't author typo'd keys ("v3" / "v_3" / "V3") and have them
# silently disappear into a "no defaults" fallback.
NP_KEYS: Tuple[str, ...] = ("NP1", "NP2")
TRIM_KEYS: Tuple[str, ...] = ("V_3", "D_5", "D_3", "J_5")
TRIM_KEYS_VJ: Tuple[str, ...] = ("V_3", "J_5")
# P-nucleotide per-end length-distribution keys. Mirror the trim
# keys (same four ends), with VJ rejecting D_5 / D_3 the same way
# the trim validation does.
P_NUCLEOTIDE_END_KEYS: Tuple[str, ...] = ("V_3", "D_5", "D_3", "J_5")
P_NUCLEOTIDE_END_KEYS_VJ: Tuple[str, ...] = ("V_3", "J_5")


@dataclass
class EmpiricalDistributionSpec:
    """One empirical distribution as a flat ``(value, weight)`` list.

    ``values`` is a list of ``(int, float)`` pairs. Values are
    non-negative ints (lengths, trim amounts). Weights are positive
    finite floats; the engine normalises them at sample time, so
    callers don't need to normalise themselves — they just need to
    be positive. Duplicate values aren't enforced at this layer (the
    engine accumulates) but discourage them stylistically.
    """

    values: List[Tuple[int, float]]

    def validate(self, name: str) -> None:
        """Shape-check this distribution. Raises ``ValueError`` on
        the first violation with a ``name``-tagged message so the
        caller knows which distribution failed (e.g. ``"NP1"`` or
        ``"V_3"``).
        """
        if not isinstance(self.values, list):
            raise ValueError(
                f"{name}: values must be a list of (int, float) pairs, "
                f"got {type(self.values).__name__}"
            )
        if not self.values:
            raise ValueError(f"{name}: empirical distribution must be non-empty")
        for i, pair in enumerate(self.values):
            if not isinstance(pair, tuple) or len(pair) != 2:
                raise ValueError(
                    f"{name}: entry {i} must be a (value, weight) tuple, got {pair!r}"
                )
            value, weight = pair
            if not isinstance(value, int) or isinstance(value, bool):
                raise ValueError(
                    f"{name}: entry {i} value must be a non-negative int, got {value!r}"
                )
            if value < 0:
                raise ValueError(
                    f"{name}: entry {i} value must be non-negative, got {value}"
                )
            if not isinstance(weight, (int, float)) or isinstance(weight, bool):
                raise ValueError(
                    f"{name}: entry {i} weight must be a finite positive float, "
                    f"got {weight!r}"
                )
            wf = float(weight)
            if not math.isfinite(wf):
                raise ValueError(
                    f"{name}: entry {i} weight must be finite, got {weight!r}"
                )
            if wf <= 0.0:
                raise ValueError(
                    f"{name}: entry {i} weight must be > 0, got {weight}"
                )


# Canonical DNA bases the NP base model supports. Bases outside
# this set are rejected at construction time — the engine only
# knows how to sample A/C/G/T at NP positions.
_NP_BASE_ALPHABET: Tuple[str, ...] = ("A", "C", "G", "T")

# Recognised `NpBaseModelSpec.kind` values. All three lower into
# concrete engine generators (`UniformNpGenerator` /
# `CategoricalNpGenerator` / `MarkovBaseGenerator`) and fold into
# the plan signature via `GenerateNPPass.parameter_signature`. See
# `docs/junction_n_addition_audit.md`.
_NP_BASE_MODEL_KINDS: Tuple[str, ...] = (
    "uniform",
    "empirical_first_base",
    "markov",
)


@dataclass
class NpBaseModelSpec:
    """One V(D)J N-addition base sampling model.

    Three kinds are recognised by the Python validator:

    - ``"uniform"`` (default-equivalent) — every NP position
      samples uniformly from A/C/G/T. Byte-identical to the
      pre-slice engine when no model is configured. ``first_base``
      and ``transitions`` must be ``None``.
    - ``"empirical_first_base"`` — every NP position samples
      independently from the supplied ``first_base`` categorical
      distribution. ``first_base`` is required; ``transitions``
      must be ``None``. This is the v1 cartridge-owned biology
      surface for the typed NP base model.
    - ``"markov"`` — true Markov chain where each base is sampled
      conditional on the previously emitted base via the
      ``transitions`` matrix. Both ``first_base`` and
      ``transitions`` are required.

      **Wired end-to-end as of the NP-Markov slice.** The
      engine ships an ``NpBaseGenerator`` trait whose
      ``MarkovBaseGenerator`` implementation threads the
      previous base through ``sample_base(rng,
      previous_base, …)``; the typed Python spec lowers
      directly into it via ``push_generate_np(...,
      markov_transitions=...)``. Plan-signature replay
      folds the full Markov payload (first-base row + 4
      transition rows in canonical A/C/G/T order). See
      ``docs/junction_n_addition_audit.md`` (Markov
      shipped) and the validation matrix's "NP base
      models / Markov N-addition" row.

    Validation rules:

    - ``kind`` must be one of ``{"uniform",
      "empirical_first_base", "markov"}``.
    - ``first_base`` (when present) is ``dict[str, float]``
      with keys from ``{"A", "C", "G", "T"}``. Values are
      finite, non-negative, with at least one strictly
      positive. Unknown bases are rejected.
    - ``transitions`` (when present) is ``dict[str, dict[str,
      float]]`` — outer key is the previously-emitted base,
      inner is the next-base categorical. Same per-base
      validation rules; each row must be non-empty and have at
      least one positive weight. Outer-key coverage of A/C/G/T
      is enforced — if any from-base is missing, the validator
      raises (a partial Markov matrix is almost certainly an
      authoring bug).
    """

    kind: str
    first_base: Optional[Dict[str, float]] = None
    transitions: Optional[Dict[str, Dict[str, float]]] = None

    def __post_init__(self) -> None:
        # Validate eagerly so a malformed cartridge fails at
        # construction time, not at simulation start.
        self.validate(name="<NpBaseModelSpec>")

    def validate(self, name: str = "<NpBaseModelSpec>") -> None:
        """Shape-check the spec. Raises ``ValueError`` on the
        first violation with a ``name``-tagged message so the
        caller knows which model failed (e.g. ``"NP1"`` /
        ``"NP2"``).
        """
        if not isinstance(self.kind, str):
            raise ValueError(
                f"{name}: kind must be a string, got {type(self.kind).__name__}"
            )
        if self.kind not in _NP_BASE_MODEL_KINDS:
            raise ValueError(
                f"{name}: unknown kind {self.kind!r}; expected one of "
                f"{list(_NP_BASE_MODEL_KINDS)}"
            )

        if self.kind == "uniform":
            if self.first_base is not None:
                raise ValueError(
                    f"{name}: kind='uniform' must not carry first_base; got {self.first_base!r}"
                )
            if self.transitions is not None:
                raise ValueError(
                    f"{name}: kind='uniform' must not carry transitions; got {self.transitions!r}"
                )
            return

        if self.kind == "empirical_first_base":
            if self.first_base is None:
                raise ValueError(
                    f"{name}: kind='empirical_first_base' requires first_base"
                )
            if self.transitions is not None:
                raise ValueError(
                    f"{name}: kind='empirical_first_base' must not carry transitions; "
                    f"use kind='markov' for previous-base conditioning"
                )
            _validate_base_categorical(self.first_base, name=f"{name}.first_base")
            return

        if self.kind == "markov":
            if self.first_base is None:
                raise ValueError(
                    f"{name}: kind='markov' requires first_base (the categorical "
                    "for the first NP position before any previous base exists)"
                )
            if self.transitions is None:
                raise ValueError(
                    f"{name}: kind='markov' requires transitions (the from-base → "
                    "to-base matrix for positions 1..n)"
                )
            _validate_base_categorical(self.first_base, name=f"{name}.first_base")
            if not isinstance(self.transitions, dict):
                raise ValueError(
                    f"{name}.transitions: must be dict[str, dict[str, float]], "
                    f"got {type(self.transitions).__name__}"
                )
            if not self.transitions:
                raise ValueError(
                    f"{name}.transitions: must be non-empty"
                )
            seen_from = set()
            for from_base, row in self.transitions.items():
                if not isinstance(from_base, str) or from_base not in _NP_BASE_ALPHABET:
                    raise ValueError(
                        f"{name}.transitions: from-base {from_base!r} is not one of "
                        f"{list(_NP_BASE_ALPHABET)}"
                    )
                seen_from.add(from_base)
                if not isinstance(row, dict) or not row:
                    raise ValueError(
                        f"{name}.transitions[{from_base!r}]: must be a non-empty "
                        f"dict[str, float] row"
                    )
                _validate_base_categorical(
                    row, name=f"{name}.transitions[{from_base!r}]"
                )
            missing = set(_NP_BASE_ALPHABET) - seen_from
            if missing:
                raise ValueError(
                    f"{name}.transitions: missing rows for from-base(s) {sorted(missing)!r}; "
                    f"a Markov NP base model must define a transition row for every "
                    f"canonical base — partial matrices are almost certainly an "
                    f"authoring bug"
                )


def _validate_base_categorical(weights: Dict[str, float], *, name: str) -> None:
    """Per-base categorical validator shared between
    ``first_base`` and per-row ``transitions`` checks. The dict
    must have ``str`` keys from the canonical A/C/G/T alphabet
    with finite non-negative float values, at least one strictly
    positive."""
    if not isinstance(weights, dict):
        raise ValueError(
            f"{name}: must be dict[str, float], got {type(weights).__name__}"
        )
    if not weights:
        raise ValueError(f"{name}: must be non-empty")
    total = 0.0
    for base, weight in weights.items():
        if not isinstance(base, str) or base not in _NP_BASE_ALPHABET:
            raise ValueError(
                f"{name}: base {base!r} is not one of {list(_NP_BASE_ALPHABET)}"
            )
        if isinstance(weight, bool) or not isinstance(weight, (int, float)):
            raise ValueError(
                f"{name}[{base!r}]: weight must be a finite non-negative number, "
                f"got {type(weight).__name__}"
            )
        wf = float(weight)
        if not math.isfinite(wf):
            raise ValueError(
                f"{name}[{base!r}]: weight must be finite, got {weight!r}"
            )
        if wf < 0.0:
            raise ValueError(
                f"{name}[{base!r}]: weight must be non-negative, got {wf}"
            )
        total += wf
    if total <= 0.0:
        raise ValueError(
            f"{name}: at least one weight must be strictly positive; sum was {total}"
        )


@dataclass
class AlleleUsageSpec:
    """Per-segment allele-usage weights authored on a cartridge.

    Each segment field is a ``{allele_name: weight}`` dict — names
    are AIRR-convention allele names (``"IGHV1-2*02"`` etc.); weights
    are positive finite floats. Weights need NOT sum to 1.0 at
    authoring time; the bridge resolver normalises per segment when
    it lowers into the dense pool-aligned ``Tuple[float, ...]`` the
    engine consumes. An empty segment dict means "fall back to
    uniform for that segment" — useful when the author only has
    data for one or two segments and wants the others uniform.

    Validation:

    - Each name is a non-empty string.
    - Each weight is a finite, strictly positive float.
    - When ``chain_type="vj"`` is supplied to
      :meth:`validate`, non-empty D-segment entries are rejected
      (D on a VJ cartridge is an authoring bug, mirrored from the
      trims / NP-base / P-nucleotide planes).
    - Allele-name existence in the cartridge's V/D/J pools is NOT
      checked here — the pool lookup lives in the bridge resolver
      (`_dataconfig_extract._allele_usage_from_models`) so the spec
      layer stays decoupled from a specific catalogue.

    Slice — Allele Usage Estimation v1. Produced either by hand or
    by :meth:`GenAIRR.ReferenceCartridgeBuilder.estimate_allele_usage`.
    """

    v: Dict[str, float] = field(default_factory=dict)
    d: Dict[str, float] = field(default_factory=dict)
    j: Dict[str, float] = field(default_factory=dict)

    def validate(self, chain_type: Optional[str] = None, *, name: str = "allele_usage") -> None:
        """Shape-check the spec. Raises :class:`ValueError` on the
        first violation tagged with the ``name`` argument so the
        caller can route the error to a clear field path
        (``"reference_models.allele_usage"``)."""
        ct = chain_type.lower() if isinstance(chain_type, str) else None
        for segment_label, mapping in (("V", self.v), ("D", self.d), ("J", self.j)):
            if not isinstance(mapping, dict):
                raise ValueError(
                    f"{name}.{segment_label.lower()}: must be a dict[str, float], "
                    f"got {type(mapping).__name__}"
                )
            if ct == "vj" and segment_label == "D" and mapping:
                raise ValueError(
                    f"{name}.d: non-empty on a VJ chain — D-segment "
                    "allele usage is meaningless without a D pool. Drop "
                    "the D entries or use a VDJ cartridge."
                )
            for allele_name, weight in mapping.items():
                if not isinstance(allele_name, str) or not allele_name:
                    raise ValueError(
                        f"{name}.{segment_label.lower()}: allele names must "
                        f"be non-empty strings, got {allele_name!r}"
                    )
                if isinstance(weight, bool) or not isinstance(
                    weight, (int, float)
                ):
                    raise ValueError(
                        f"{name}.{segment_label.lower()}[{allele_name!r}]: "
                        f"weight must be a finite positive number, got "
                        f"{type(weight).__name__}"
                    )
                wf = float(weight)
                if not math.isfinite(wf):
                    raise ValueError(
                        f"{name}.{segment_label.lower()}[{allele_name!r}]: "
                        f"weight must be finite, got {weight!r}"
                    )
                if wf <= 0.0:
                    raise ValueError(
                        f"{name}.{segment_label.lower()}[{allele_name!r}]: "
                        f"weight must be strictly positive, got {wf}"
                    )

    def nonempty_segments(self) -> Tuple[str, ...]:
        """Return the segment labels whose dict is non-empty.

        Used by the manifest block to surface which segments the
        cartridge actually overrides — empty segments fall through
        to uniform at the bridge."""
        out = []
        for label, mapping in (("V", self.v), ("D", self.d), ("J", self.j)):
            if mapping:
                out.append(label)
        return tuple(out)


@dataclass
class ReferenceEmpiricalModels:
    """The **empirical models** plane of a reference cartridge —
    default distributions the engine samples from when the user
    doesn't override them at recombine time.

    A cartridge has four planes:

    - identity (species, locus, reference set, name, source)
    - catalogue (V/D/J alleles)
    - rules (see :class:`ReferenceRulesSpec`) — how to interpret
      alleles
    - **empirical models** — this object — defaults for NP lengths
      + exonuclease trim distributions

    Two dicts, both keyed by stable string identifiers:

    - ``np_lengths``: ``{"NP1", "NP2"}`` → length distribution.
    - ``trims``: ``{"V_3", "D_5", "D_3", "J_5"}`` → trim distribution.

    Each value is an :class:`EmpiricalDistributionSpec` — a flat list
    of ``(value, weight)`` pairs the engine consumes directly (no
    marginalisation magic in between, unlike the legacy nested-dict
    extraction path).

    Empty dicts are valid (model absent → fall back to the legacy
    ``DataConfig.NP_lengths``/``trim_dicts`` extraction, then to the
    uniform placeholder). Partially-populated dicts are valid (NP1
    set, NP2 absent → NP1 used; NP2 falls through to legacy).

    ``validate(chain_type=...)`` enforces shape: known keys only,
    non-empty distributions, positive finite weights. The
    ``chain_type`` argument (``"vj"`` or ``"vdj"``) rejects D-trim
    models on VJ cartridges — D trims are nonsense without a D
    segment, and silently ignoring them would mask the user's
    authoring mistake.

    See ``docs/reference_cartridge.md`` for end-to-end examples.
    """

    np_lengths: Dict[str, EmpiricalDistributionSpec] = field(default_factory=dict)
    trims: Dict[str, EmpiricalDistributionSpec] = field(default_factory=dict)
    # V(D)J N-addition base sampling models, keyed by NP region
    # (``"NP1"`` / ``"NP2"``). Slice — Typed NP base model. Empty
    # dict (default) means every NP position samples uniformly
    # A/C/G/T — byte-identical to the pre-slice engine.
    # Authored shape: ``{"NP1": NpBaseModelSpec(kind="empirical_first_base",
    # first_base={"A": 0.2, "C": 0.3, "G": 0.3, "T": 0.2})}``. See
    # ``docs/junction_n_addition_audit.md`` for the cascade
    # discipline (typed → uniform; legacy `NP_first_bases` /
    # `NP_transitions` auto-lift is deliberately deferred).
    np_bases: Dict[str, "NpBaseModelSpec"] = field(default_factory=dict)
    # Per-end P-nucleotide length distributions, keyed by junction
    # side label (``"V_3"`` / ``"D_5"`` / ``"D_3"`` / ``"J_5"``).
    # Slice — P-nucleotide v1. Empty dict (default) means the
    # pipeline omits every P-pass — byte-identical to the
    # pre-slice baseline. The cartridge author opts in by
    # populating per-end distributions; bases are templated
    # palindromic complements of the source allele's post-trim,
    # post-orientation coding flank (deterministic, not sampled).
    # See ``docs/p_nucleotide_design.md`` for the full design.
    p_nucleotide_lengths: Dict[str, EmpiricalDistributionSpec] = field(
        default_factory=dict
    )
    # Per-segment allele-usage weights. Slice — Allele Usage
    # Estimation v1. ``None`` means uniform sampling across each
    # segment's pool (byte-identical to the pre-slice baseline).
    # When present, the bridge resolver narrows / reweights the
    # `SampleAllelePass` distribution per segment; an explicit
    # ``Experiment.recombine(v_allele_weights=...)`` kwarg still
    # takes precedence per the documented precedence order. Legacy
    # ``DataConfig.gene_use_dict`` is NOT auto-lifted — that orphan
    # boundary is preserved per
    # ``docs/allele_usage_estimation_design.md`` §2.3.
    allele_usage: Optional["AlleleUsageSpec"] = None

    def validate(self, chain_type: Optional[str] = None) -> None:
        """Shape-check the whole bundle. Raises ``ValueError`` on
        the first violation. ``chain_type`` is ``"vj"`` or ``"vdj"``;
        when supplied, D-trim entries on a VJ cartridge are rejected.
        """
        if not isinstance(self.np_lengths, dict):
            raise ValueError(
                f"np_lengths must be a dict, got {type(self.np_lengths).__name__}"
            )
        if not isinstance(self.trims, dict):
            raise ValueError(
                f"trims must be a dict, got {type(self.trims).__name__}"
            )

        for key, spec in self.np_lengths.items():
            if key not in NP_KEYS:
                raise ValueError(
                    f"np_lengths key {key!r} is not recognised; expected one of {NP_KEYS}"
                )
            if not isinstance(spec, EmpiricalDistributionSpec):
                raise ValueError(
                    f"np_lengths[{key!r}] must be an EmpiricalDistributionSpec, "
                    f"got {type(spec).__name__}"
                )
            spec.validate(name=f"np_lengths[{key}]")

        ct = chain_type.lower() if isinstance(chain_type, str) else None
        for key, spec in self.trims.items():
            if key not in TRIM_KEYS:
                raise ValueError(
                    f"trims key {key!r} is not recognised; expected one of {TRIM_KEYS}"
                )
            if ct == "vj" and key not in TRIM_KEYS_VJ:
                raise ValueError(
                    f"trims[{key!r}] is set but the cartridge is a VJ chain "
                    f"with no D segment — D trims must not be authored on VJ "
                    f"reference models. Drop {key!r} or use a VDJ cartridge."
                )
            if not isinstance(spec, EmpiricalDistributionSpec):
                raise ValueError(
                    f"trims[{key!r}] must be an EmpiricalDistributionSpec, "
                    f"got {type(spec).__name__}"
                )
            spec.validate(name=f"trims[{key}]")

        # NP base sampling models. Same key vocabulary as
        # `np_lengths` — `"NP1"` / `"NP2"`. Each value is an
        # `NpBaseModelSpec` (already self-validates in
        # `__post_init__`; we re-validate here under the typed
        # name to surface a clear `np_bases[…]` error path).
        if not isinstance(self.np_bases, dict):
            raise ValueError(
                f"np_bases must be a dict, got {type(self.np_bases).__name__}"
            )
        for key, spec in self.np_bases.items():
            if key not in NP_KEYS:
                raise ValueError(
                    f"np_bases key {key!r} is not recognised; expected one of {NP_KEYS}"
                )
            if not isinstance(spec, NpBaseModelSpec):
                raise ValueError(
                    f"np_bases[{key!r}] must be an NpBaseModelSpec, "
                    f"got {type(spec).__name__}"
                )
            spec.validate(name=f"np_bases[{key}]")

        # Per-end P-nucleotide length distributions. Same key
        # vocabulary as `trims`: V_3 / D_5 / D_3 / J_5. Empty
        # dict is the byte-identical default. D-end keys are
        # rejected on VJ chains the same way trim D-end keys
        # are — a VJ cartridge with `p_nucleotide_lengths["D_5"]`
        # is an authoring bug, surfaced loudly rather than
        # silently dropped.
        if not isinstance(self.p_nucleotide_lengths, dict):
            raise ValueError(
                "p_nucleotide_lengths must be a dict, got "
                f"{type(self.p_nucleotide_lengths).__name__}"
            )
        for key, spec in self.p_nucleotide_lengths.items():
            if key not in P_NUCLEOTIDE_END_KEYS:
                raise ValueError(
                    f"p_nucleotide_lengths key {key!r} is not recognised; "
                    f"expected one of {P_NUCLEOTIDE_END_KEYS}"
                )
            if ct == "vj" and key not in P_NUCLEOTIDE_END_KEYS_VJ:
                raise ValueError(
                    f"p_nucleotide_lengths[{key!r}] is set but the cartridge "
                    f"is a VJ chain with no D segment — D-end P-nucleotides "
                    f"must not be authored on VJ reference models. Drop "
                    f"{key!r} or use a VDJ cartridge."
                )
            if not isinstance(spec, EmpiricalDistributionSpec):
                raise ValueError(
                    f"p_nucleotide_lengths[{key!r}] must be an "
                    f"EmpiricalDistributionSpec, got {type(spec).__name__}"
                )
            spec.validate(name=f"p_nucleotide_lengths[{key}]")

        # Per-segment allele-usage spec (Slice — Allele Usage
        # Estimation v1). ``None`` is the byte-identical default;
        # when present, the spec is delegated to
        # :meth:`AlleleUsageSpec.validate` which enforces the
        # per-name + per-weight shape AND the VJ-rejects-D
        # boundary.
        if self.allele_usage is not None:
            if not isinstance(self.allele_usage, AlleleUsageSpec):
                raise ValueError(
                    f"allele_usage must be an AlleleUsageSpec, got "
                    f"{type(self.allele_usage).__name__}"
                )
            self.allele_usage.validate(
                chain_type=ct, name="allele_usage"
            )
