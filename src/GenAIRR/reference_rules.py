"""Reference-rules cartridge spec (Python side).

``DataConfig`` can carry an optional :class:`ReferenceRulesSpec`
describing the rules the engine should consult when interpreting the
catalogue — V/J anchor expectations, allowed nucleotide alphabet,
severity policy. The loader in
:func:`GenAIRR._refdata_resolver.dataconfig_to_refdata` transfers the
spec into the Rust ``RefDataConfig.rules`` slice; the Rust side
remains the authority on whether the final cartridge is usable.

This module only enforces *shape* (single-character bases / amino
acids, valid severity strings, a minimum sensible alphabet) so the
bridge fails fast with a clear Python ``ValueError`` instead of
crossing the PyO3 boundary with garbage. Biological policy lives in
the Rust validator.

Default values mirror the Rust ``ReferenceRules::default()`` —
construct a ``ReferenceRulesSpec()`` with no arguments to get the
same lenient defaults the loader applies for unknown loci.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple


# Allowed severity strings — match the Rust ``RefDataIssueSeverity``
# variants, lowercased to keep the Python surface readable.
_VALID_SEVERITIES: Tuple[str, ...] = ("fatal", "curatable")

# Minimum required alphabet — any A/C/G/T/N reference catalogue needs
# the four canonical bases. ``N`` is recommended but not required (a
# strict spec may forbid it). Custom alphabets can ADD letters
# (e.g. IUPAC ambiguity codes) but not remove the four canonical bases.
_MIN_REQUIRED_BASES: Tuple[str, ...] = ("A", "C", "G", "T")

# Single-letter amino acids the spec recognises. Stop codon (`*`) is
# included so callers can model unusual anchor rules; productivity
# semantics live elsewhere.
_VALID_AA = set("ACDEFGHIKLMNPQRSTVWY*")


@dataclass
class AnchorRuleSpec:
    """Anchor rule for one segment (V or J) — one of the **rules**
    in a reference cartridge.

    Rules define **how the engine interprets** the catalogue's
    alleles. The V/J anchor rule says which amino acids the anchor
    codon may translate to, whether an anchor is required at all, and
    how anchor-related issues are classified (Fatal vs Curatable) for
    the strict-vs-AllowCuratable validation gate.

    Fields mirror the Rust ``AnchorRule`` struct one-to-one. Severity
    strings are ``"fatal"`` or ``"curatable"`` (lower-case) — the same
    vocabulary the validator's issue dicts use.

    The default constructor is intentionally restrictive (anchor
    required, both severities curatable). Pass an empty
    ``expected_aa`` list with ``required=False`` to model an
    anchorless rule (legacy/exploratory catalogues).

    See ``docs/reference_cartridge.md`` for the cartridge model and
    the full validation/curation/compile flow.
    """

    expected_aa: List[str]
    required: bool = True
    missing_severity: str = "curatable"
    mismatch_severity: str = "curatable"

    def validate(self, segment: str) -> None:
        """Shape-check this rule. Raises ``ValueError`` on the first
        violation. ``segment`` is ``"V"`` or ``"J"`` for the diagnostic
        message — the validator itself is segment-agnostic.
        """
        if not isinstance(self.expected_aa, (list, tuple)):
            raise ValueError(
                f"{segment} anchor expected_aa must be a list of single-character "
                f"strings, got {type(self.expected_aa).__name__}"
            )
        if self.required and not self.expected_aa:
            raise ValueError(
                f"{segment} anchor expected_aa must be non-empty when required=True"
            )
        for aa in self.expected_aa:
            if not isinstance(aa, str) or len(aa) != 1:
                raise ValueError(
                    f"{segment} anchor expected_aa entries must be one-character "
                    f"strings, got {aa!r}"
                )
            if aa.upper() not in _VALID_AA:
                raise ValueError(
                    f"{segment} anchor expected_aa entry {aa!r} is not a recognised "
                    f"amino-acid letter (one of {sorted(_VALID_AA)})"
                )
        if not isinstance(self.required, bool):
            raise ValueError(
                f"{segment} anchor required must be a bool, "
                f"got {type(self.required).__name__}"
            )
        for field_name in ("missing_severity", "mismatch_severity"):
            value = getattr(self, field_name)
            if value not in _VALID_SEVERITIES:
                raise ValueError(
                    f"{segment} anchor {field_name}={value!r} must be one of "
                    f"{_VALID_SEVERITIES}"
                )


@dataclass
class ReferenceRulesSpec:
    """The **rules** plane of a reference cartridge — how the engine
    should interpret the catalogue's alleles.

    A cartridge has four planes:

    - **identity** (species, locus, reference set, name, source)
    - **catalogue** (V/D/J alleles)
    - **rules** — this object
    - **empirical models** (NP-length + trim distributions, see
      :class:`ReferenceEmpiricalModels`)

    Plus an orthogonal concept:

    - **curation** — which subset of the catalogue participates
      in simulation (``Experiment.curate_refdata(...)``).

    This spec mirrors the Rust ``ReferenceRules`` struct one-to-one
    and is shape-validated here (single-char bases / amino acids,
    severity strings, A/C/G/T present in the alphabet). The Rust
    ``RefDataConfig.validate()`` remains the authority on whether the
    catalogue is **usable** under these rules.

    Default constructor reproduces the Rust ``ReferenceRules::default()``:
    A/C/G/T/N alphabet, V expects ``C``, J accepts ``W`` or ``F``,
    both severities curatable. When the spec is attached to a
    :class:`DataConfig`, ``dataconfig_to_refdata`` ships it verbatim
    into the cartridge — overriding the loader's bundled-locus
    inference (which narrows J to ``["W"]`` for IGH, ``["F"]`` for
    IGK/IGL/TR*).

    See ``docs/reference_cartridge.md`` for the cartridge model and
    end-to-end examples (custom J anchor ``Y``, extended alphabet,
    non-standard species).
    """

    allowed_bases: List[str] = field(
        default_factory=lambda: ["A", "C", "G", "T", "N"]
    )
    v_anchor: AnchorRuleSpec = field(default_factory=lambda: AnchorRuleSpec(["C"]))
    j_anchor: AnchorRuleSpec = field(default_factory=lambda: AnchorRuleSpec(["W", "F"]))

    def validate(self) -> None:
        """Shape-check the spec. Raises ``ValueError`` on the first
        violation. Designed to fail fast before the spec crosses the
        PyO3 boundary."""
        if not isinstance(self.allowed_bases, (list, tuple)):
            raise ValueError(
                f"allowed_bases must be a list, got {type(self.allowed_bases).__name__}"
            )
        if not self.allowed_bases:
            raise ValueError("allowed_bases must be a non-empty list")
        for b in self.allowed_bases:
            if not isinstance(b, str) or len(b) != 1:
                raise ValueError(
                    f"allowed_bases entries must be single-character strings, got {b!r}"
                )
            if not b.isascii() or not b.isalpha():
                raise ValueError(
                    f"allowed_bases entries must be ASCII letters, got {b!r}"
                )
        canonical = {b.upper() for b in self.allowed_bases}
        missing = [base for base in _MIN_REQUIRED_BASES if base not in canonical]
        if missing:
            raise ValueError(
                f"allowed_bases must include the canonical DNA bases; "
                f"missing {missing}. Custom alphabets may ADD letters "
                f"(IUPAC ambiguity codes, etc.) but not remove A/C/G/T."
            )
        if not isinstance(self.v_anchor, AnchorRuleSpec):
            raise ValueError(
                f"v_anchor must be an AnchorRuleSpec, "
                f"got {type(self.v_anchor).__name__}"
            )
        if not isinstance(self.j_anchor, AnchorRuleSpec):
            raise ValueError(
                f"j_anchor must be an AnchorRuleSpec, "
                f"got {type(self.j_anchor).__name__}"
            )
        self.v_anchor.validate("V")
        self.j_anchor.validate("J")
