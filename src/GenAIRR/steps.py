"""
Simulation step descriptors — internal IR for the clause→step compiler.

Each step maps directly to a C pipeline feature flag and its
associated parameters. **Steps are not the user API**: pass clauses
from :mod:`GenAIRR.ops` to the Experiment phase methods (``.mutate``,
``.sequence``, ``.observe``, etc.) and the compiler in
``protocol._clauses_to_steps`` produces the matching Step descriptors
internally.

Example::

    from GenAIRR import Experiment, Productivity
    from GenAIRR.ops import rate, model

    result = (
        Experiment.on("human_igh")
        .mutate(rate(0.01, 0.05), model("s5f"))
        .run(n=1000, seed=42, productivity=Productivity.PRODUCTIVE_ONLY)
    )
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional


# ── Base ─────────────────────────────────────────────────────────

@dataclass
class Step:
    """Base class for all simulation steps."""

    # Subclasses set these class-level attributes
    _feature: str = ""        # C feature flag name (empty = always on)
    _category: str = "core"   # For display: core, mutation, simulation, corruption

    def configure(self, sim) -> None:
        """Apply this step's feature flags and parameters to a CSimulator."""
        if self._feature:
            sim.set_feature(self._feature, True)

    def __repr__(self):
        params = {k: v for k, v in self.__dict__.items()
                  if not k.startswith('_') and v is not None}
        if params:
            args = ", ".join(f"{k}={v!r}" for k, v in params.items())
            return f"{self.__class__.__name__}({args})"
        return f"{self.__class__.__name__}()"


# ── Core rearrangement (always present) ──────────────────────────

@dataclass
class Rearrange(Step):
    """V(D)J rearrangement with trimming, NP insertion, and assembly."""
    _feature: str = field(default="", init=False, repr=False)
    _category: str = field(default="core", init=False, repr=False)

    def configure(self, sim) -> None:
        pass  # Always included — no feature flag needed


# ── Allele locking ───────────────────────────────────────────────

@dataclass
class LockAlleles(Step):
    """Lock rearrangement to use specific allele(s).

    Each key is a segment ("v", "d", "j", "c") mapped to a single
    allele name or a list of allele names. The simulator picks
    uniformly among the locked alleles for each segment.
    """
    locks: dict = field(default_factory=dict)

    _feature: str = field(default="", init=False, repr=False)
    _category: str = field(default="core", init=False, repr=False)

    def configure(self, sim) -> None:
        for segment, names in self.locks.items():
            if isinstance(names, str):
                names = [names]
            for name in names:
                sim.lock_allele(segment, name)


# ── Mutation ─────────────────────────────────────────────────────

@dataclass
class Mutate(Step):
    """Somatic hypermutation. ``model`` selects the kernel —
    ``"s5f"`` (default, context-dependent) or ``"uniform"`` (T2-4,
    position-independent)."""
    min_rate: float = 0.01
    max_rate: float = 0.05
    model: str = "s5f"

    _feature: str = field(default="mutate", init=False, repr=False)
    _category: str = field(default="mutation", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("mutate", True)
        sim.set_mutation_model(self.model)
        sim.set_param("min_mutation_rate", self.min_rate)
        sim.set_param("max_mutation_rate", self.max_rate)


@dataclass
class SimulateCSR(Step):
    """Class-switch recombination — isotype-aware mutation rates.
    Honors ``model`` (s5f or uniform) the same way ``Mutate`` does;
    CSR rate adjustment applies in either kernel."""
    min_rate: float = 0.01
    max_rate: float = 0.05
    model: str = "s5f"

    _feature: str = field(default="csr", init=False, repr=False)
    _category: str = field(default="mutation", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("mutate", True)
        sim.set_feature("csr", True)
        sim.set_mutation_model(self.model)
        sim.set_param("min_mutation_rate", self.min_rate)
        sim.set_param("max_mutation_rate", self.max_rate)


# ── Simulation ops ───────────────────────────────────────────────

@dataclass
class SelectionPressure(Step):
    """Antigen-driven selection pressure after SHM.

    Walks every R-mutation in the coding region (V/NP1/D/NP2/J) and
    accepts/reverts each based on its region:
      - Anchor codon (V Cys, J W/F): kept with anchor_r_acceptance
        (default 0 — anchor disruption breaks the V/J fold).
      - CDR (CDR1/CDR2/CDR3 incl. NP1/D/NP2): kept with cdr_r_acceptance.
      - FWR (FR1/FR2/FR3/FR4): kept with fwr_r_acceptance.
    Anchorless V or J alleles are skipped on that segment.
    """
    strength: float = 0.5
    cdr_r_acceptance: float = 0.85
    fwr_r_acceptance: float = 0.40
    anchor_r_acceptance: float = 0.0

    _feature: str = field(default="selection_pressure", init=False, repr=False)
    _category: str = field(default="simulation", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("selection_pressure", True)
        sim.set_param("selection_strength", self.strength)
        sim.set_param("cdr_r_acceptance", self.cdr_r_acceptance)
        sim.set_param("fwr_r_acceptance", self.fwr_r_acceptance)
        sim.set_param("anchor_r_acceptance", self.anchor_r_acceptance)


@dataclass
class SimulateDGeneInversion(Step):
    """D gene segment inversion (reverse-complement)."""
    probability: float = 0.15

    _feature: str = field(default="d_inversion", init=False, repr=False)
    _category: str = field(default="simulation", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("d_inversion", True)
        sim.set_param("d_inversion_prob", self.probability)


@dataclass
class SimulateReceptorRevision(Step):
    """V-gene replacement with old-V footprint."""
    probability: float = 0.05
    footprint_min: int = 5
    footprint_max: int = 20

    _feature: str = field(default="receptor_revision", init=False, repr=False)
    _category: str = field(default="simulation", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("receptor_revision", True)
        sim.set_param("revision_prob", self.probability)
        sim.set_param("footprint_min", self.footprint_min)
        sim.set_param("footprint_max", self.footprint_max)


# ── Corruption ops ───────────────────────────────────────────────

@dataclass
class Corrupt5Prime(Step):
    """5' end truncation or adapter read-through."""
    min_remove: int = 1
    max_remove: int = 20
    min_add: int = 1
    max_add: int = 10
    _feature: str = field(default="corrupt_5_prime", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("corrupt_5_prime", True)
        sim.set_param("corrupt_5_remove_min", self.min_remove)
        sim.set_param("corrupt_5_remove_max", self.max_remove)
        sim.set_param("corrupt_5_add_min", self.min_add)
        sim.set_param("corrupt_5_add_max", self.max_add)


@dataclass
class Corrupt3Prime(Step):
    """3' end truncation or adapter read-through."""
    min_remove: int = 1
    max_remove: int = 20
    min_add: int = 1
    max_add: int = 10
    _feature: str = field(default="corrupt_3_prime", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("corrupt_3_prime", True)
        sim.set_param("corrupt_3_remove_min", self.min_remove)
        sim.set_param("corrupt_3_remove_max", self.max_remove)
        sim.set_param("corrupt_3_add_min", self.min_add)
        sim.set_param("corrupt_3_add_max", self.max_add)


@dataclass
class CorruptQuality(Step):
    """Position-dependent Illumina sequencing errors."""
    base_error_rate: float = 0.001
    peak_error_rate: float = 0.02

    _feature: str = field(default="quality_errors", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("quality_errors", True)
        sim.set_param("base_error_rate", self.base_error_rate)
        sim.set_param("peak_error_rate", self.peak_error_rate)


@dataclass
class SimulatePairedEnd(Step):
    """Paired-end merge artifacts with bathtub error profile."""
    read_length: int = 300

    _feature: str = field(default="paired_end", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("paired_end", True)
        sim.set_param("pe_read_length", self.read_length)


@dataclass
class PCRAmplification(Step):
    """Uniform PCR polymerase errors."""
    error_rate: float = 1e-4
    n_cycles: int = 30

    _feature: str = field(default="pcr", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("pcr", True)
        sim.set_param("pcr_error_rate", self.error_rate)
        sim.set_param("pcr_cycles", self.n_cycles)


@dataclass
class SimulateUMI(Step):
    """Random UMI barcode prepend."""
    umi_length: int = 12

    _feature: str = field(default="umi", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("umi", True)
        sim.set_param("umi_length", self.umi_length)


@dataclass
class PrimerMask(Step):
    """FR1 overwritten with germline primer sequence."""
    mask_length: int = 0  # 0 = full FR1

    _feature: str = field(default="primer_mask", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("primer_mask", True)
        sim.set_param("primer_mask_length", self.mask_length)


@dataclass
class ReverseComplement(Step):
    """Antisense read orientation (~50% of reads)."""
    probability: float = 0.5

    _feature: str = field(default="reverse_complement", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("reverse_complement", True)
        sim.set_param("rc_prob", self.probability)


@dataclass
class SpikeContaminants(Step):
    """Cross-sample or phiX contamination."""
    probability: float = 0.01
    contaminant_type: str = "random"  # "random" or "phix"

    _feature: str = field(default="contaminants", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("contaminants", True)
        sim.set_param("contamination_prob", self.probability)
        ct = 1 if self.contaminant_type == "phix" else 0
        sim.set_param("contaminant_type", ct)


@dataclass
class SkewBaseComposition(Step):
    """Nanopore/PacBio homopolymer-targeted indels."""
    error_rate: float = 0.03
    min_run_length: int = 3
    insertion_bias: float = 0.6

    _feature: str = field(default="long_read_errors", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("long_read_errors", True)
        sim.set_param("long_read_error_rate", self.error_rate)
        sim.set_param("min_run_length", self.min_run_length)
        sim.set_param("insertion_bias", self.insertion_bias)


@dataclass
class InsertIndels(Step):
    """Random insertion/deletion simulation."""
    probability: float = 0.01

    _feature: str = field(default="indels", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("indels", True)
        sim.set_param("indel_prob", self.probability)


@dataclass
class InsertNs(Step):
    """Ambiguous N base insertion."""
    probability: float = 0.01

    _feature: str = field(default="insert_ns", init=False, repr=False)
    _category: str = field(default="corruption", init=False, repr=False)

    def configure(self, sim) -> None:
        sim.set_feature("insert_ns", True)
        sim.set_param("n_prob", self.probability)
