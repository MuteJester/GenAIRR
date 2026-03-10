"""
Clause functions for the phase-call DSL.

Each function returns a typed, frozen dataclass that is passed as an
argument to one of the five Experiment phase methods::

    from GenAIRR import Experiment
    from GenAIRR.ops import rate, with_antigen_selection

    result = (
        Experiment.on("human_igh")
        .mutate(rate(0.02, 0.08), with_antigen_selection(0.5))
        .run(n=1000, seed=42)
    )

Clause types are checked at phase-call time — passing a mutation clause
to ``.recombine()`` raises ``TypeError`` immediately.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, Union


# =====================================================================
# Base types — one per phase
# =====================================================================

class Clause:
    """Base for all DSL clauses. Never instantiated directly."""

    def __label__(self) -> str:
        return self.__class__.__name__


class RecombineClause(Clause):
    """Clause valid inside ``.recombine(...)``."""
    pass


class MutateClause(Clause):
    """Clause valid inside ``.mutate(...)``."""
    pass


class PrepareClause(Clause):
    """Clause valid inside ``.prepare(...)``."""
    pass


class SequenceClause(Clause):
    """Clause valid inside ``.sequence(...)``."""
    pass


class ObserveClause(Clause):
    """Clause valid inside ``.observe(...)``."""
    pass


# =====================================================================
# Recombine clauses
# =====================================================================

@dataclass(frozen=True)
class UsingClause(RecombineClause):
    v: Optional[Tuple[str, ...]] = None
    d: Optional[Tuple[str, ...]] = None
    j: Optional[Tuple[str, ...]] = None
    c: Optional[Tuple[str, ...]] = None

    def __label__(self) -> str:
        parts = []
        for seg, names in [("V", self.v), ("D", self.d),
                           ("J", self.j), ("C", self.c)]:
            if names:
                parts.append(f"{seg}={','.join(names)}")
        return f"Using {', '.join(parts)}"


@dataclass(frozen=True)
class DInversionClause(RecombineClause):
    prob: float = 0.15

    def __label__(self) -> str:
        return f"D-gene inversion (prob={self.prob})"


@dataclass(frozen=True)
class ReceptorRevisionClause(RecombineClause):
    prob: float = 0.05
    footprint_min: int = 5
    footprint_max: int = 20

    def __label__(self) -> str:
        return (f"Receptor revision (prob={self.prob}, "
                f"footprint={self.footprint_min}-{self.footprint_max})")


# =====================================================================
# Mutate clauses
# =====================================================================

@dataclass(frozen=True)
class ModelClause(MutateClause):
    name: str = "s5f"

    def __label__(self) -> str:
        return f"Model: {self.name.upper()}"


@dataclass(frozen=True)
class RateClause(MutateClause):
    min_rate: float = 0.01
    max_rate: float = 0.05

    def __label__(self) -> str:
        return f"Rate: {self.min_rate} – {self.max_rate}"


@dataclass(frozen=True)
class IsotypeRatesClause(MutateClause):

    def __label__(self) -> str:
        return "Isotype-adjusted rates (CSR)"


@dataclass(frozen=True)
class AntigenSelectionClause(MutateClause):
    strength: float = 0.5
    cdr_r_acceptance: float = 0.85
    fwr_r_acceptance: float = 0.40

    def __label__(self) -> str:
        return f"Antigen selection (strength={self.strength})"


# =====================================================================
# Prepare clauses
# =====================================================================

@dataclass(frozen=True)
class UMIClause(PrepareClause):
    length: int = 12

    def __label__(self) -> str:
        return f"UMI ({self.length} nt)"


@dataclass(frozen=True)
class PCRClause(PrepareClause):
    error_rate: float = 1e-4
    cycles: int = 30

    def __label__(self) -> str:
        return f"PCR (rate={self.error_rate}, cycles={self.cycles})"


@dataclass(frozen=True)
class PrimerMaskClause(PrepareClause):
    length: int = 0  # 0 = full FR1

    def __label__(self) -> str:
        if self.length == 0:
            return "Primer mask (full FR1)"
        return f"Primer mask ({self.length} bp)"


# =====================================================================
# Sequence clauses
# =====================================================================

@dataclass(frozen=True)
class PairedEndClause(SequenceClause):
    read_length: int = 300

    def __label__(self) -> str:
        return f"Paired-end ({self.read_length} bp)"


@dataclass(frozen=True)
class LongReadClause(SequenceClause):
    error_rate: float = 0.03
    min_run_length: int = 3
    insertion_bias: float = 0.6

    def __label__(self) -> str:
        return f"Long-read errors (rate={self.error_rate})"


@dataclass(frozen=True)
class FivePrimeLossClause(SequenceClause):
    min_remove: int = 1
    max_remove: int = 20
    min_add: int = 1
    max_add: int = 10

    def __label__(self) -> str:
        return f"5' signal loss (remove={self.min_remove}-{self.max_remove})"


@dataclass(frozen=True)
class ThreePrimeLossClause(SequenceClause):
    min_remove: int = 1
    max_remove: int = 20
    min_add: int = 1
    max_add: int = 10

    def __label__(self) -> str:
        return f"3' signal loss (remove={self.min_remove}-{self.max_remove})"


@dataclass(frozen=True)
class QualityProfileClause(SequenceClause):
    base: float = 0.001
    peak: float = 0.02

    def __label__(self) -> str:
        return f"Quality profile (base={self.base}, peak={self.peak})"


@dataclass(frozen=True)
class ReverseComplementClause(SequenceClause):
    prob: float = 0.5

    def __label__(self) -> str:
        return f"Reverse complement (prob={self.prob})"


# =====================================================================
# Observe clauses
# =====================================================================

@dataclass(frozen=True)
class ContaminantsClause(ObserveClause):
    rate: float = 0.01
    source: str = "random"

    def __label__(self) -> str:
        return f"Contaminants ({self.source}, rate={self.rate})"


@dataclass(frozen=True)
class IndelsClause(ObserveClause):
    prob: float = 0.01

    def __label__(self) -> str:
        return f"Indels (prob={self.prob})"


@dataclass(frozen=True)
class NsClause(ObserveClause):
    prob: float = 0.01

    def __label__(self) -> str:
        return f"N bases (prob={self.prob})"


# =====================================================================
# Constructor functions — the public API
# =====================================================================

# ── Helpers ───────────────────────────────────────────────────────

def _check_probability(value: float, name: str) -> None:
    """Validate a probability parameter is in [0, 1]."""
    if value < 0 or value > 1:
        raise ValueError(
            f"{name} must be between 0 and 1, got {value}"
        )


def _check_positive_int(value: int, name: str) -> None:
    """Validate a positive integer parameter."""
    if value < 1:
        raise ValueError(
            f"{name} must be at least 1, got {value}"
        )


def _check_non_negative_float(value: float, name: str) -> None:
    """Validate a non-negative float parameter."""
    if value < 0:
        raise ValueError(
            f"{name} must be non-negative, got {value}"
        )


# ── Recombine ─────────────────────────────────────────────────────

def using(*, v=None, d=None, j=None, c=None) -> UsingClause:
    """Restrict rearrangement to specific allele(s).

    Args:
        v: V allele name(s) — string or list of strings.
        d: D allele name(s).
        j: J allele name(s).
        c: C allele name(s).

    Example::

        using(v="IGHV1-2*01")
        using(v=["IGHV1-2*01", "IGHV1-2*02"], j="IGHJ4*02")
    """
    def _norm(x):
        if x is None:
            return None
        return (x,) if isinstance(x, str) else tuple(x)
    return UsingClause(v=_norm(v), d=_norm(d), j=_norm(j), c=_norm(c))


def with_d_inversion(prob: float = 0.15) -> DInversionClause:
    """D-gene segment inversion during recombination.

    Args:
        prob: Probability of inversion per sequence (0–1).
    """
    _check_probability(prob, "prob")
    return DInversionClause(prob=prob)


def with_receptor_revision(
    prob: float = 0.05,
    footprint: Tuple[int, int] = (5, 20),
) -> ReceptorRevisionClause:
    """V-gene replacement with old-V footprint.

    Args:
        prob: Probability of revision per sequence (0–1).
        footprint: (min, max) footprint length in nucleotides.
    """
    _check_probability(prob, "prob")
    fmin, fmax = footprint
    if fmin < 0:
        raise ValueError(f"footprint min must be non-negative, got {fmin}")
    if fmax < fmin:
        raise ValueError(
            f"footprint max ({fmax}) must be >= footprint min ({fmin})"
        )
    return ReceptorRevisionClause(
        prob=prob, footprint_min=fmin, footprint_max=fmax,
    )


# ── Mutate ────────────────────────────────────────────────────────

def model(name: str = "s5f") -> ModelClause:
    """Select the somatic hypermutation model.

    Args:
        name: ``"s5f"`` (context-dependent, default) or
              ``"uniform"`` (position-independent, future).

    Note:
        ``model("uniform")`` is accepted but not yet implemented
        in the C backend. A ``RuntimeWarning`` is emitted at compile
        time. Sequences will use S5F mutation until the uniform
        model is available.
    """
    valid = ("s5f", "uniform")
    if name not in valid:
        raise ValueError(
            f"Unknown mutation model {name!r}. Choose from: {valid}"
        )
    return ModelClause(name=name)


def rate(min_rate: float = 0.01, max_rate: float = 0.05) -> RateClause:
    """Set the per-sequence mutation rate range.

    Each simulated sequence draws its mutation rate uniformly
    from ``[min_rate, max_rate]``.

    Args:
        min_rate: Minimum mutation rate (fraction of sequence mutated).
        max_rate: Maximum mutation rate.
    """
    _check_non_negative_float(min_rate, "min_rate")
    _check_non_negative_float(max_rate, "max_rate")
    if min_rate > max_rate:
        raise ValueError(
            f"min_rate ({min_rate}) must be <= max_rate ({max_rate})"
        )
    return RateClause(min_rate=min_rate, max_rate=max_rate)


def with_antigen_selection(
    strength: float = 0.5,
    cdr_r_acceptance: float = 0.85,
    fwr_r_acceptance: float = 0.40,
) -> AntigenSelectionClause:
    """Simulate antigen-driven selection pressure.

    Models germinal center selection by accepting or rejecting
    replacement mutations based on IMGT region (CDR vs FWR).

    Args:
        strength: Overall selection intensity (0–1).
        cdr_r_acceptance: Acceptance rate for CDR replacement mutations.
        fwr_r_acceptance: Acceptance rate for FWR replacement mutations.
    """
    _check_probability(strength, "strength")
    _check_probability(cdr_r_acceptance, "cdr_r_acceptance")
    _check_probability(fwr_r_acceptance, "fwr_r_acceptance")
    return AntigenSelectionClause(
        strength=strength,
        cdr_r_acceptance=cdr_r_acceptance,
        fwr_r_acceptance=fwr_r_acceptance,
    )


def with_isotype_rates() -> IsotypeRatesClause:
    """Apply isotype-specific mutation rate adjustments (CSR).

    Simulates class switch recombination by adjusting per-sequence
    mutation rates based on the switched isotype.
    """
    return IsotypeRatesClause()


# ── Prepare ───────────────────────────────────────────────────────

def with_umi(length: int = 12) -> UMIClause:
    """Prepend a random UMI barcode.

    Args:
        length: UMI length in nucleotides (common: 8, 10, 12, 16).
    """
    _check_positive_int(length, "length")
    return UMIClause(length=length)


def with_pcr(error_rate: float = 1e-4, cycles: int = 30) -> PCRClause:
    """Simulate PCR amplification errors.

    Args:
        error_rate: Per-base per-cycle error rate.
        cycles: Number of PCR cycles.
    """
    _check_non_negative_float(error_rate, "error_rate")
    _check_positive_int(cycles, "cycles")
    return PCRClause(error_rate=error_rate, cycles=cycles)


def with_primer_mask(length: int = 0) -> PrimerMaskClause:
    """Overwrite FR1 with germline primer sequence.

    Args:
        length: Number of bases to mask. 0 means full FR1 region.
    """
    if length < 0:
        raise ValueError(f"length must be non-negative, got {length}")
    return PrimerMaskClause(length=length)


# ── Sequence ──────────────────────────────────────────────────────

def paired_end(read_length: int = 300) -> PairedEndClause:
    """Paired-end sequencing mode.

    Args:
        read_length: Cycles per read (common: 150, 250, 300).
    """
    _check_positive_int(read_length, "read_length")
    return PairedEndClause(read_length=read_length)


def long_read(
    error_rate: float = 0.03,
    min_run_length: int = 3,
    insertion_bias: float = 0.6,
) -> LongReadClause:
    """Nanopore/PacBio long-read error model.

    Args:
        error_rate: Base error rate (scales with homopolymer run length).
        min_run_length: Minimum homopolymer length to target.
        insertion_bias: Fraction of errors that are insertions (0–1).
    """
    _check_non_negative_float(error_rate, "error_rate")
    _check_positive_int(min_run_length, "min_run_length")
    _check_probability(insertion_bias, "insertion_bias")
    return LongReadClause(
        error_rate=error_rate,
        min_run_length=min_run_length,
        insertion_bias=insertion_bias,
    )


def with_5prime_loss(
    min_remove: int = 1,
    max_remove: int = 20,
    min_add: int = 1,
    max_add: int = 10,
) -> FivePrimeLossClause:
    """Simulate 5' end signal loss (truncation or adapter read-through).

    Each sequence randomly gets one of three events: removal only (40%),
    addition only (30%), or removal followed by addition (30%). The
    amount removed/added is drawn uniformly from the given ranges.

    Args:
        min_remove: Minimum bases to remove from the 5' end.
        max_remove: Maximum bases to remove from the 5' end.
        min_add: Minimum random bases to add at the 5' end.
        max_add: Maximum random bases to add at the 5' end.
    """
    if min_remove < 0 or max_remove < min_remove:
        raise ValueError(
            f"Need 0 <= min_remove <= max_remove, "
            f"got min_remove={min_remove}, max_remove={max_remove}"
        )
    if min_add < 0 or max_add < min_add:
        raise ValueError(
            f"Need 0 <= min_add <= max_add, "
            f"got min_add={min_add}, max_add={max_add}"
        )
    return FivePrimeLossClause(
        min_remove=min_remove, max_remove=max_remove,
        min_add=min_add, max_add=max_add,
    )


def with_3prime_loss(
    min_remove: int = 1,
    max_remove: int = 20,
    min_add: int = 1,
    max_add: int = 10,
) -> ThreePrimeLossClause:
    """Simulate 3' end signal loss (truncation or adapter read-through).

    Each sequence randomly gets one of three events: removal only (40%),
    addition only (30%), or removal followed by addition (30%). The
    amount removed/added is drawn uniformly from the given ranges.

    Args:
        min_remove: Minimum bases to remove from the 3' end.
        max_remove: Maximum bases to remove from the 3' end.
        min_add: Minimum random bases to add at the 3' end.
        max_add: Maximum random bases to add at the 3' end.
    """
    if min_remove < 0 or max_remove < min_remove:
        raise ValueError(
            f"Need 0 <= min_remove <= max_remove, "
            f"got min_remove={min_remove}, max_remove={max_remove}"
        )
    if min_add < 0 or max_add < min_add:
        raise ValueError(
            f"Need 0 <= min_add <= max_add, "
            f"got min_add={min_add}, max_add={max_add}"
        )
    return ThreePrimeLossClause(
        min_remove=min_remove, max_remove=max_remove,
        min_add=min_add, max_add=max_add,
    )


def with_quality_profile(
    base: float = 0.001,
    peak: float = 0.02,
) -> QualityProfileClause:
    """Position-dependent Illumina sequencing error profile.

    Error rate increases linearly from 5' (base) to 3' (peak).

    Args:
        base: Error rate at 5' end (Q30 ~ 0.001).
        peak: Error rate at 3' end (Q17 ~ 0.02).
    """
    _check_non_negative_float(base, "base")
    _check_non_negative_float(peak, "peak")
    return QualityProfileClause(base=base, peak=peak)


def with_reverse_complement(prob: float = 0.5) -> ReverseComplementClause:
    """Reverse-complement a fraction of reads.

    Args:
        prob: Probability of reverse-complementing each read (0–1).
    """
    _check_probability(prob, "prob")
    return ReverseComplementClause(prob=prob)


# ── Observe ───────────────────────────────────────────────────────

def with_contaminants(
    rate: float = 0.01,
    source: str = "random",
) -> ContaminantsClause:
    """Spike in contaminant sequences.

    Args:
        rate: Per-sequence contamination probability (0–1).
        source: ``"random"`` or ``"phix"``.
    """
    _check_probability(rate, "rate")
    valid_sources = ("random", "phix")
    if source not in valid_sources:
        raise ValueError(
            f"Unknown contaminant source {source!r}. "
            f"Choose from: {valid_sources}"
        )
    return ContaminantsClause(rate=rate, source=source)


def with_indels(prob: float = 0.01) -> IndelsClause:
    """Random insertion/deletion events in the observed data.

    Args:
        prob: Per-position indel probability (0–1).
    """
    _check_probability(prob, "prob")
    return IndelsClause(prob=prob)


def with_ns(prob: float = 0.01) -> NsClause:
    """Ambiguous N bases in the observed data.

    Args:
        prob: Per-position N-insertion probability (0–1).
    """
    _check_probability(prob, "prob")
    return NsClause(prob=prob)


# =====================================================================
# Public exports
# =====================================================================

__all__ = [
    # Recombine
    "using", "with_d_inversion", "with_receptor_revision",
    # Mutate
    "model", "rate", "with_antigen_selection", "with_isotype_rates",
    # Prepare
    "with_umi", "with_pcr", "with_primer_mask",
    # Sequence
    "paired_end", "long_read",
    "with_5prime_loss", "with_3prime_loss",
    "with_quality_profile", "with_reverse_complement",
    # Observe
    "with_contaminants", "with_indels", "with_ns",
]
