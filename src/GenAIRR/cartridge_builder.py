"""Reference cartridge authoring API.

Public surface for constructing a :class:`~GenAIRR.DataConfig` from
FASTA inputs with an auditable build trail. v1 ships the facade
and report infrastructure; statistical estimators
(:meth:`estimate_allele_usage`, :meth:`estimate_trim_distributions`,
:meth:`estimate_np_length_distributions`,
:meth:`estimate_np_base_model`,
:meth:`estimate_p_nucleotide_lengths`, :meth:`estimate_shm_rates`)
are deferred to a follow-up slice and are explicitly absent from
this module — users hit `AttributeError` if they try to call them,
which the audit
[`docs/reference_cartridge_authoring_audit.md`](../../docs/reference_cartridge_authoring_audit.md)
§11 documents as the intentional v1 boundary.

The v1 contract:

- :class:`ReferenceCartridgeBuilder` is the only public class.
- :meth:`ReferenceCartridgeBuilder.from_fasta` is the only
  constructor.
- :meth:`infer_identity` / :meth:`infer_v_subregions` /
  :meth:`with_rules` / :meth:`with_models` are the only stages.
- :meth:`build` returns a plain :class:`~GenAIRR.DataConfig` —
  drop-in replaceable with the 106 bundled cartridges.
- :meth:`report` returns a :class:`CartridgeBuildReport` whose
  :meth:`~CartridgeBuildReport.to_dict` is JSON-clean.

See
[`docs/reference_cartridge_authoring_audit.md`](../../docs/reference_cartridge_authoring_audit.md)
for the full design.
"""
from __future__ import annotations

import io
from dataclasses import asdict, dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, Dict, IO, List, Optional, Union

from .alleles.allele import DAllele, JAllele, VAllele


class _SafeAnchorMixin:
    """`_find_anchor` shim that degrades gracefully when the
    ``GenAIRR._native._anchor`` C resolver isn't available.

    The bundled cartridges' alleles are pickled with anchors
    already populated, so the C resolver is only needed when
    constructing fresh alleles from FASTA — which is exactly
    what the builder does. Some development environments (and
    minimal CI installs) don't ship the native module; in
    those, we want :class:`ReferenceCartridgeBuilder` to still
    succeed with `anchor = None` and a warning in the report.

    Mixin precedence: the mixin's ``_find_anchor`` is called
    first; it delegates to ``super()._find_anchor()`` and
    catches :class:`ImportError`. Any other exception bubbles.
    """

    def _find_anchor(self):  # type: ignore[override]
        try:
            return super()._find_anchor()  # type: ignore[misc]
        except ImportError:
            self.anchor = None
            self.anchor_meta = None
            return None


class _BuilderVAllele(_SafeAnchorMixin, VAllele):
    """V allele used by :class:`ReferenceCartridgeBuilder`.

    Identical to :class:`~GenAIRR.alleles.allele.VAllele` in
    every respect except that ``_find_anchor`` degrades
    gracefully when the native C resolver isn't built (see
    :class:`_SafeAnchorMixin`).
    """


class _BuilderJAllele(_SafeAnchorMixin, JAllele):
    """J allele used by :class:`ReferenceCartridgeBuilder` —
    same shape as the V variant."""


class _BuilderDAllele(DAllele):
    """D allele used by :class:`ReferenceCartridgeBuilder`.

    D's ``_find_anchor`` is a no-op in the base class so no
    mixin is needed.
    """
from .dataconfig.config_info import ConfigInfo
from .dataconfig.data_config import DataConfig
from .dataconfig.enums import ChainType, Species
from .genotype_priors import PopulationGenotypeModel
from .reference_models import (
    AlleleUsageSpec,
    EmpiricalDistributionSpec,
    NP_KEYS,
    NpBaseModelSpec,
    P_NUCLEOTIDE_END_KEYS,
    P_NUCLEOTIDE_END_KEYS_VJ,
    ReferenceEmpiricalModels,
    TRIM_KEYS,
    TRIM_KEYS_VJ,
)

_NP_CANONICAL_BASES = ("A", "C", "G", "T")
from .reference_rules import ReferenceRulesSpec
from .utilities.imgt_regions import compute_v_region_boundaries
from .utilities.misc import parse_fasta


# Accepted input shapes for FASTA inputs: a filesystem path, a path-
# like, or an already-open text file. Strings that contain a newline
# are treated as raw FASTA content (useful for tests).
FastaInput = Union[str, Path, IO[str]]


def _resolve_chain_type(value: Union[ChainType, str]) -> ChainType:
    """Normalise the ``chain_type`` argument into a :class:`ChainType`.

    Accepts either the enum directly or the canonical string label
    (``"BCR_HEAVY"`` etc.). Raises :class:`ValueError` on an
    unrecognised value with a clear catalogue of accepted strings.
    """
    if isinstance(value, ChainType):
        return value
    if isinstance(value, str):
        normalised = value.strip().upper()
        for member in ChainType:
            if member.name == normalised or member.value.upper() == normalised:
                return member
    raise ValueError(
        f"chain_type must be a ChainType enum or one of "
        f"{[m.name for m in ChainType]}, got {value!r}"
    )


def _resolve_species(value: Union[Species, str, None]) -> Optional[Species]:
    """Normalise the ``species`` argument into a :class:`Species`.

    ``None`` passes through. Accepts the enum or the canonical name
    label (case-insensitive); raises :class:`ValueError` on an
    unrecognised string.
    """
    if value is None or isinstance(value, Species):
        return value
    if isinstance(value, str):
        normalised = value.strip()
        for member in Species:
            if member.name == normalised.upper() or member.value == normalised:
                return member
    raise ValueError(
        f"species must be a Species enum or recognised name string, "
        f"got {value!r}"
    )


def _locus_label_from_chain_type(chain_type: ChainType) -> str:
    """Conventional locus label derived from a chain type.

    Mirrors the AIRR-C / IMGT locus naming so the cartridge
    manifest's ``identity.locus`` matches what downstream tools
    expect.
    """
    mapping = {
        ChainType.BCR_HEAVY: "IGH",
        ChainType.BCR_LIGHT_KAPPA: "IGK",
        ChainType.BCR_LIGHT_LAMBDA: "IGL",
        ChainType.TCR_ALPHA: "TRA",
        ChainType.TCR_BETA: "TRB",
        ChainType.TCR_GAMMA: "TRG",
        ChainType.TCR_DELTA: "TRD",
    }
    return mapping[chain_type]


def _open_fasta(source: FastaInput):
    """Open a FASTA input. Returns a `(file, close_on_exit)` tuple."""
    if hasattr(source, "read"):
        return source, False
    if isinstance(source, (str, Path)):
        text = str(source)
        # Raw FASTA content shortcut for tests: a string with a `>`
        # somewhere and at least one newline. Path-like strings that
        # happen to contain a `>` are vanishingly rare (no POSIX
        # filename starts with `>`), so this heuristic is safe.
        if "\n" in text and text.lstrip().startswith(">"):
            return io.StringIO(text), True
        path = Path(text)
        return open(path, "r"), True
    raise TypeError(
        f"FASTA input must be a path, path-like, or open text file, "
        f"got {type(source).__name__}"
    )


def _split_tie_set(raw: str) -> List[str]:
    """Parse a comma-separated AIRR tie-set call into a clean list.

    Empty / whitespace-only entries drop out. The returned list
    preserves insertion order — important for the
    ``ambiguous="truth_first"`` policy which keeps only the first
    entry."""
    return [token.strip() for token in raw.split(",") if token.strip()]


def _allele_names_for_segment(
    buckets: Dict[str, List[Any]],
) -> "set[str]":
    """Flatten a builder's per-gene allele bucket dict into a set
    of canonical allele names. Used by
    :meth:`ReferenceCartridgeBuilder.estimate_allele_usage` to
    decide which AIRR-record allele names are "known" to the
    cartridge under construction."""
    out: "set[str]" = set()
    for allele_list in buckets.values():
        for allele in allele_list:
            name = getattr(allele, "name", None)
            if isinstance(name, str) and name:
                out.add(name)
    return out


def _load_rearrangements(
    source: Any,
) -> "tuple[List[Dict[str, Any]], str]":
    """Normalise a ``rearrangements`` argument into a
    ``(records, source_label)`` tuple.

    Accepted inputs:

    - ``list[dict]`` — used verbatim; source label
      ``"records:N"`` carrying the row count.
    - path-like (``str`` / ``Path``) — parsed via
      :class:`csv.DictReader` with ``delimiter='\\t'`` (AIRR-C TSV
      convention); source label is the path.
    - open text file handle — parsed via :class:`csv.DictReader`;
      source label is ``"file:N"``.

    Raises :class:`TypeError` for unrecognised shapes."""
    import csv

    if isinstance(source, list):
        return list(source), f"records:{len(source)}"
    if hasattr(source, "read"):
        rows = list(csv.DictReader(source, delimiter="\t"))
        return rows, f"file:{len(rows)}"
    if isinstance(source, (str, Path)):
        path = Path(source)
        with open(path, "r", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return rows, str(path)
    raise TypeError(
        f"rearrangements must be a list[dict], path, or open text "
        f"file, got {type(source).__name__}"
    )


def _parse_allele_name(raw_header: str) -> str:
    """Pull the allele name from a FASTA header.

    Accepts:
    - bare names (``>IGHV1-2*02``).
    - IMGT-style pipe-delimited headers
      (``>X62106|IGHV1-2*02|Homo sapiens|F|...``) — the first field
      that looks like ``GeneFamily*allele`` wins, falling back to the
      second pipe-token.
    - Whitespace-separated headers — first token after the ``>``.
    """
    header = raw_header.lstrip(">").strip()
    if "|" in header:
        parts = [p.strip() for p in header.split("|") if p.strip()]
        for part in parts:
            if "*" in part:
                return part
        if len(parts) >= 2:
            return parts[1]
        return parts[0]
    return header.split()[0]


@dataclass
class CartridgeBuildReport:
    """Audit trail for a :class:`ReferenceCartridgeBuilder` run.

    The report accumulates one entry per builder stage in
    :attr:`stages`. Each entry is a JSON-clean dict carrying the
    stage's ``inputs``, ``inferred`` payload, and any per-stage
    warnings. Cross-stage rejections (alleles dropped by the
    parser, etc.) accumulate in :attr:`rejected`; cross-stage
    warnings (e.g. the cartridge has no V alleles with
    ``gapped_seq``) accumulate in :attr:`warnings`.

    :attr:`manifest_snapshot` and :attr:`checksum_at_build_time`
    are populated by :meth:`ReferenceCartridgeBuilder.build` only
    — they capture the cartridge state at the moment ``build()``
    finalised.
    """

    stages: List[Dict[str, Any]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    rejected: List[Dict[str, Any]] = field(default_factory=list)
    manifest_snapshot: Optional[Dict[str, Any]] = None
    checksum_at_build_time: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-clean dict of the report.

        Equivalent to :func:`dataclasses.asdict` but exposed as a
        method so a future revision can shim non-JSON-clean
        sub-objects without breaking the call sites.
        """
        return asdict(self)


class ReferenceCartridgeBuilder:
    """Stage-based builder for a :class:`~GenAIRR.DataConfig`.

    See :mod:`GenAIRR.cartridge_builder` module docs and
    `docs/reference_cartridge_authoring_audit.md` for the design.

    Usage::

        builder = ReferenceCartridgeBuilder.from_fasta(
            v_fasta="igh_v.fasta",
            j_fasta="igh_j.fasta",
            d_fasta="igh_d.fasta",
            chain_type="BCR_HEAVY",
        )
        builder.infer_identity(species="HUMAN", reference_set="MY_REF")
        builder.infer_v_subregions()
        builder.with_rules(my_rules_spec)
        builder.with_models(my_empirical_models)
        cfg = builder.build()
        report = builder.report()

    Each stage method mutates the builder in place and returns
    ``self`` so calls can be chained.
    """

    def __init__(self, chain_type: ChainType) -> None:
        self._chain_type: ChainType = chain_type
        self._v_alleles: Dict[str, List[VAllele]] = {}
        self._d_alleles: Dict[str, List[DAllele]] = {}
        self._j_alleles: Dict[str, List[JAllele]] = {}
        self._name: Optional[str] = None
        self._metadata: Optional[ConfigInfo] = None
        self._reference_rules: Optional[ReferenceRulesSpec] = None
        self._reference_models: Optional[ReferenceEmpiricalModels] = None
        self._genotype_priors: Optional[PopulationGenotypeModel] = None
        self._report = CartridgeBuildReport()

    # ──────────────────────────────────────────────────────────
    # Constructors
    # ──────────────────────────────────────────────────────────

    @classmethod
    def from_fasta(
        cls,
        *,
        v_fasta: FastaInput,
        j_fasta: FastaInput,
        d_fasta: Optional[FastaInput] = None,
        c_fasta: Optional[FastaInput] = None,
        chain_type: Union[ChainType, str],
    ) -> "ReferenceCartridgeBuilder":
        """Build a builder from FASTA inputs.

        ``v_fasta`` and ``j_fasta`` are required; ``d_fasta`` is
        required only when the chain type has a D segment
        (``BCR_HEAVY`` / ``TCR_BETA`` / ``TCR_DELTA``). ``c_fasta``
        is accepted but the parsed C alleles are not yet wired into
        the resulting cartridge — v1 limits to V/D/J authoring.
        Records dropped during parsing land in
        :attr:`CartridgeBuildReport.rejected`.

        Each ``*_fasta`` argument can be a path-like, an open text
        file, or a raw FASTA string (a string containing a newline
        and starting with ``>`` is parsed directly).
        """
        chain = _resolve_chain_type(chain_type)
        if not chain.has_d and d_fasta is not None:
            raise ValueError(
                f"d_fasta supplied for chain_type {chain.name}, which has "
                f"no D segment. Drop d_fasta or use a chain type with a "
                f"D segment ({[m.name for m in ChainType if m.has_d]})."
            )
        if chain.has_d and d_fasta is None:
            raise ValueError(
                f"chain_type {chain.name} has a D segment — d_fasta is "
                f"required."
            )

        builder = cls(chain_type=chain)

        # Per-segment parse with structured rejection tracking. Each
        # rejected allele lands in a dict explaining which input it
        # came from and why it was dropped — auditable downstream.
        inputs_summary: Dict[str, Any] = {
            "chain_type": chain.name,
        }
        v_parsed, v_rejected = _parse_segment_fasta(
            v_fasta, segment="V", builder=builder, allele_cls=_BuilderVAllele
        )
        inputs_summary["v_alleles_parsed"] = sum(len(g) for g in builder._v_alleles.values())
        inputs_summary["v_alleles_rejected"] = len(v_rejected)

        j_parsed, j_rejected = _parse_segment_fasta(
            j_fasta, segment="J", builder=builder, allele_cls=_BuilderJAllele
        )
        inputs_summary["j_alleles_parsed"] = sum(len(g) for g in builder._j_alleles.values())
        inputs_summary["j_alleles_rejected"] = len(j_rejected)

        if d_fasta is not None:
            d_parsed, d_rejected = _parse_segment_fasta(
                d_fasta, segment="D", builder=builder, allele_cls=_BuilderDAllele
            )
            inputs_summary["d_alleles_parsed"] = sum(
                len(g) for g in builder._d_alleles.values()
            )
            inputs_summary["d_alleles_rejected"] = len(d_rejected)

        if c_fasta is not None:
            # v1 boundary: parsing C alleles is accepted but the
            # built cartridge does not yet populate `c_alleles`. The
            # report records the count so a future slice can wire
            # them in without breaking the API.
            cf_handle, cf_close = _open_fasta(c_fasta)
            try:
                c_count = sum(1 for _ in parse_fasta(cf_handle))
            finally:
                if cf_close:
                    cf_handle.close()
            inputs_summary["c_alleles_parsed_but_unused_in_v1"] = c_count
            builder._report.warnings.append(
                "c_fasta parsed but ignored: v1 does not populate "
                "DataConfig.c_alleles; supply via with_models / manual "
                "DataConfig editing if you need them"
            )

        # First stage entry: from_fasta.
        builder._report.stages.append(
            {
                "stage": "from_fasta",
                "inputs": inputs_summary,
                "inferred": {
                    "v_genes": len(builder._v_alleles),
                    "j_genes": len(builder._j_alleles),
                    "d_genes": len(builder._d_alleles),
                },
                "warnings": [],
            }
        )
        # Per-segment rejected alleles flow to the top-level
        # rejected list with their segment + reason annotated.
        builder._report.rejected.extend(v_rejected + j_rejected)
        if d_fasta is not None:
            builder._report.rejected.extend(d_rejected)
        return builder

    # ──────────────────────────────────────────────────────────
    # Identity / annotations / programmable surfaces
    # ──────────────────────────────────────────────────────────

    def infer_identity(
        self,
        *,
        species: Optional[Union[Species, str]] = None,
        locus: Optional[str] = None,
        reference_set: Optional[str] = None,
        name: Optional[str] = None,
        source: str = "ReferenceCartridgeBuilder",
    ) -> "ReferenceCartridgeBuilder":
        """Populate the cartridge's identity / metadata plane.

        ``species`` / ``locus`` / ``reference_set`` are user-supplied
        — none of these can be reliably inferred from FASTA alone, so
        the audit (§7) requires explicit input. When ``locus`` is
        omitted the conventional label is derived from the chain
        type (``BCR_HEAVY → "IGH"``, etc.).

        ``source`` is recorded on every parsed allele's ``.source``
        field so downstream consumers can identify cartridges
        produced by this builder vs the bundled OGRDB / IMGT
        loaders.
        """
        species_enum = _resolve_species(species)
        chain = self._chain_type
        locus_value = locus or _locus_label_from_chain_type(chain)
        cartridge_name = name or f"USER_{locus_value}"

        inputs: Dict[str, Any] = {
            "species": species_enum.name if species_enum else None,
            "locus": locus_value,
            "reference_set": reference_set,
            "name": cartridge_name,
            "source": source,
        }
        warnings: List[str] = []
        if species_enum is None:
            warnings.append(
                "species not provided — manifest identity.species will "
                "be null. Cartridges without species cannot use the "
                "bundled S5F kernel selection heuristics."
            )

        self._name = cartridge_name
        self._metadata = ConfigInfo(
            species=species_enum or Species.HUMAN,
            chain_type=chain,
            reference_set=reference_set or "user-supplied",
            last_updated=date.today(),
            has_d=chain.has_d,
        )

        # Tag every parsed allele with the source string so the
        # manifest's per-allele provenance reflects the build origin.
        for buckets in (self._v_alleles, self._j_alleles, self._d_alleles):
            for allele_list in buckets.values():
                for allele in allele_list:
                    allele.source = source
                    if species_enum is not None:
                        allele.species = species_enum.value
                    allele.locus = locus_value

        self._report.stages.append(
            {
                "stage": "infer_identity",
                "inputs": inputs,
                "inferred": {
                    "name": cartridge_name,
                    "metadata_set": True,
                },
                "warnings": warnings,
            }
        )
        return self

    def infer_v_subregions(self) -> "ReferenceCartridgeBuilder":
        """Derive per-V-allele IMGT subregion intervals.

        Walks every V allele with a non-empty ``gapped_seq`` and
        runs :func:`compute_v_region_boundaries`, attaching the
        resulting ``{label: (start, end)}`` dict to
        ``allele.subregions``. Alleles without ``gapped_seq`` (or
        with a derivation that raises) are left unset and counted
        in the per-stage report.

        Idempotent — calling twice re-derives and overwrites the
        previous output, with a ``replaced=True`` flag on the
        report entry.
        """
        replaced = any(
            entry.get("stage") == "infer_v_subregions"
            for entry in self._report.stages
        )

        annotated = 0
        skipped_no_gapped: List[str] = []
        skipped_derivation_failed: List[Dict[str, Any]] = []
        for gene_alleles in self._v_alleles.values():
            for allele in gene_alleles:
                gapped = getattr(allele, "gapped_seq", None)
                # `compute_v_region_boundaries` walks IMGT-gapped
                # positions assuming `.` is the gap marker. A
                # sequence with no dots isn't IMGT-numbered and
                # would produce degenerate boundaries; treat
                # that case as "skip + warn" the same way an
                # absent `gapped_seq` is treated.
                if not gapped or "." not in gapped:
                    skipped_no_gapped.append(allele.name)
                    continue
                try:
                    bounds = compute_v_region_boundaries(allele)
                except Exception as exc:
                    skipped_derivation_failed.append(
                        {"allele": allele.name, "reason": str(exc)}
                    )
                    continue
                # Store as the dict shape `_resolve_v_subregions`
                # accepts: ``{label: (start, end)}``.
                allele.subregions = {
                    label: (int(s), int(e)) for label, (s, e) in bounds.items()
                }
                annotated += 1

        warnings: List[str] = []
        if skipped_no_gapped:
            warnings.append(
                f"{len(skipped_no_gapped)} V allele(s) lack gapped_seq — "
                f"subregion attribution will route to "
                f"n_v_unannotated_mutations"
            )
        if skipped_derivation_failed:
            warnings.append(
                f"{len(skipped_derivation_failed)} V allele(s) had "
                f"compute_v_region_boundaries failures — see "
                f"report.rejected for details"
            )

        self._report.stages.append(
            {
                "stage": "infer_v_subregions",
                "inputs": {
                    "derivation": "imgt_gapped",
                    "replaced": replaced,
                },
                "inferred": {
                    "alleles_annotated": annotated,
                    "alleles_skipped_no_gapped": len(skipped_no_gapped),
                    "alleles_skipped_derivation_failed": len(
                        skipped_derivation_failed
                    ),
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(
            [
                {
                    "stage": "infer_v_subregions",
                    "allele": entry["allele"],
                    "reason": entry["reason"],
                }
                for entry in skipped_derivation_failed
            ]
        )
        return self

    def with_rules(
        self, reference_rules: ReferenceRulesSpec
    ) -> "ReferenceCartridgeBuilder":
        """Attach an explicit :class:`ReferenceRulesSpec`.

        The user-supplied spec is validated immediately so an
        author's bug surfaces at build-stage time, not at compile
        time. The report records the rule fields the spec defines.
        """
        if not isinstance(reference_rules, ReferenceRulesSpec):
            raise TypeError(
                f"with_rules expects a ReferenceRulesSpec instance, got "
                f"{type(reference_rules).__name__}"
            )
        reference_rules.validate()
        self._reference_rules = reference_rules
        self._report.stages.append(
            {
                "stage": "with_rules",
                "inputs": {
                    "allowed_bases": getattr(
                        reference_rules, "allowed_bases", None
                    ),
                    "v_anchor_required": getattr(
                        getattr(reference_rules, "v_anchor", None),
                        "required",
                        None,
                    ),
                    "j_anchor_required": getattr(
                        getattr(reference_rules, "j_anchor", None),
                        "required",
                        None,
                    ),
                },
                "inferred": {"rules_attached": True},
                "warnings": [],
            }
        )
        return self

    def with_models(
        self, reference_models: ReferenceEmpiricalModels
    ) -> "ReferenceCartridgeBuilder":
        """Attach an explicit :class:`ReferenceEmpiricalModels`.

        The supplied bundle is validated against the chain type
        (D-end keys rejected on VJ chains, etc.). The report records
        which typed planes are populated.
        """
        if not isinstance(reference_models, ReferenceEmpiricalModels):
            raise TypeError(
                f"with_models expects a ReferenceEmpiricalModels "
                f"instance, got {type(reference_models).__name__}"
            )
        chain_label = "vdj" if self._chain_type.has_d else "vj"
        reference_models.validate(chain_type=chain_label)
        self._reference_models = reference_models
        self._report.stages.append(
            {
                "stage": "with_models",
                "inputs": {
                    "chain_type_label": chain_label,
                },
                "inferred": {
                    "np_length_keys": sorted(
                        (reference_models.np_lengths or {}).keys()
                    ),
                    "trim_keys": sorted(
                        (reference_models.trims or {}).keys()
                    ),
                    "np_base_model_keys": sorted(
                        (reference_models.np_bases or {}).keys()
                    ),
                    "p_nucleotide_length_keys": sorted(
                        (reference_models.p_nucleotide_lengths or {}).keys()
                    ),
                    "allele_usage_segments": list(
                        reference_models.allele_usage.nonempty_segments()
                    ) if reference_models.allele_usage is not None else [],
                },
                "warnings": [],
            }
        )
        return self

    # ──────────────────────────────────────────────────────────
    # Estimators
    # ──────────────────────────────────────────────────────────

    def _build_working_cfg(self):
        """A lightweight DataConfig carrying only the parsed catalogues — used by
        genotype-prior validation/estimation to resolve gene/allele names and
        construct throwaway Genotypes for novel functional validation."""
        return DataConfig(
            name=self._name or "WORKING",
            metadata=self._metadata,
            v_alleles=self._v_alleles,
            d_alleles=self._d_alleles or None,
            j_alleles=self._j_alleles,
            c_alleles=None,
        )

    def set_genotype_priors(
        self, model: PopulationGenotypeModel
    ) -> "ReferenceCartridgeBuilder":
        """Attach a hand-authored population genotype prior, validated against
        this cartridge's chain type and catalogue. Chainable."""
        from GenAIRR.genotype import Genotype

        if not isinstance(model, PopulationGenotypeModel):
            raise TypeError(
                f"set_genotype_priors expects a PopulationGenotypeModel, "
                f"got {type(model).__name__}")
        chain_label = "vdj" if self._chain_type.has_d else "vj"
        model.validate(chain_type=chain_label)  # catalogue-free

        # Catalogue-aware checks against the builder's pools.
        cfg = self._build_working_cfg()
        by_seg = {"V": cfg.v_alleles, "D": cfg.d_alleles or {}, "J": cfg.j_alleles}
        for table_name, table in (("allele_frequencies", model.allele_frequencies),
                                  ("haplotype_deletion_prob", model.haplotype_deletion_prob)):
            for seg, genes in (table or {}).items():
                catalogue = by_seg[seg]
                for gene, payload in genes.items():
                    if gene not in catalogue:
                        raise ValueError(
                            f"genotype_priors.{table_name}: {seg} gene {gene!r} is not "
                            f"in the cartridge")
                    if table_name == "allele_frequencies":
                        names = {a.name for a in catalogue[gene]}
                        for allele in payload:
                            if allele not in names:
                                raise ValueError(
                                    f"genotype_priors.allele_frequencies: {allele!r} is "
                                    f"not a known allele of {gene!r}")
        # Novels: reuse Genotype.add_novel_allele for full functional validation.
        helper = Genotype.from_dataconfig(cfg)
        for nv in model.novel_alleles:
            helper.add_novel_allele(
                nv.name, base=nv.base_allele, sequence=nv.sequence.upper(),
                segment=nv.segment, allow_nonfunctional=nv.allow_nonfunctional)

        self._genotype_priors = model
        self._report.stages.append({
            "stage": "set_genotype_priors",
            "inputs": {"model_id": model.model_id, "source": model.source,
                       "chain_type_label": chain_label},
            "inferred": {"model_checksum": model.content_checksum(),
                         "novel_allele_count": len(model.novel_alleles)},
            "warnings": [],
        })
        return self

    def estimate_genotype_priors(
        self, genotypes, **kwargs
    ) -> "ReferenceCartridgeBuilder":
        """Estimate a population genotype prior from observed ``Genotype`` objects
        and attach it (chainable). Thin wrapper over
        :meth:`PopulationGenotypeModel.from_genotypes` followed by
        :meth:`set_genotype_priors`."""
        model = PopulationGenotypeModel.from_genotypes(
            genotypes, cfg=self._build_working_cfg(), **kwargs)
        return self.set_genotype_priors(model)

    def estimate_allele_usage(
        self,
        rearrangements: Any,
        *,
        min_count: float = 1.0,
        ambiguous: str = "fractional",
        replace: bool = True,
    ) -> "ReferenceCartridgeBuilder":
        """Estimate per-segment allele-usage weights from observed
        AIRR rearrangement records.

        ``rearrangements`` accepts:

        - a list of dicts (each row is one AIRR record),
        - a path-like (filesystem path to an AIRR TSV — parsed via
          ``csv.DictReader`` with tab delimiter), or
        - an open text file handle pointing at AIRR TSV.

        ``ambiguous`` selects the tie-set policy (column values
        like ``"IGHV1*01,IGHV2*01"``):

        - ``"fractional"`` (default): split one row's credit
          ``1.0`` evenly across all known alleles in the tie set.
          Unknown allele names in the tie set are excluded; if
          NO names in the tie set are known to the cartridge, the
          row is recorded as ``unknown_allele`` and skipped.
        - ``"truth_first"``: credit only the first comma-separated
          allele in the tie set. Matches the existing
          :func:`GenAIRR._mcp_summary` convention.
        - ``"reject"``: drop ambiguous (multi-call) rows entirely
          and record them in ``report.rejected``.

        ``min_count`` drops alleles whose final per-segment count
        is strictly below the threshold. Pre-normalisation; the
        per-segment weights remaining after the drop are then
        renormalised to sum to ``1.0``.

        ``replace`` (default ``True``) controls idempotency. When
        ``True``, calling :meth:`estimate_allele_usage` twice
        overwrites the previous spec and writes ``replaced=True``
        on the new stage entry. When ``False``, the second call
        raises :class:`ValueError`.

        Updates ``self._reference_models.allele_usage`` so a
        downstream :meth:`build` carries the estimated spec into
        the cartridge.
        """
        if ambiguous not in ("fractional", "truth_first", "reject"):
            raise ValueError(
                f"ambiguous must be one of 'fractional' / 'truth_first' / "
                f"'reject', got {ambiguous!r}"
            )
        previously_estimated = any(
            entry.get("stage") == "estimate_allele_usage"
            for entry in self._report.stages
        )
        if previously_estimated and not replace:
            raise ValueError(
                "estimate_allele_usage already ran; pass replace=True to "
                "overwrite the previous spec"
            )

        records, source_label = _load_rearrangements(rearrangements)
        chain_has_d = self._chain_type.has_d
        v_pool = _allele_names_for_segment(self._v_alleles)
        d_pool = _allele_names_for_segment(self._d_alleles)
        j_pool = _allele_names_for_segment(self._j_alleles)

        v_counts: Dict[str, float] = {}
        d_counts: Dict[str, float] = {}
        j_counts: Dict[str, float] = {}
        skipped = {
            "missing_required_column": 0,
            "unknown_allele": {"V": 0, "D": 0, "J": 0},
            "missing_d_call_on_vdj": 0,
            "ambiguous_rejected": 0,
        }
        rejected_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []
        d_on_vj_warned = False

        for row_idx, row in enumerate(records):
            v_raw = (row.get("v_call") or "").strip()
            j_raw = (row.get("j_call") or "").strip()
            d_raw = (row.get("d_call") or "").strip()

            # Required column check.
            if not v_raw or not j_raw:
                skipped["missing_required_column"] += 1
                rejected_entries.append(
                    {
                        "stage": "estimate_allele_usage",
                        "row_index": row_idx,
                        "reason": "missing_required_column",
                    }
                )
                continue

            # VDJ requires d_call; VJ ignores any d_call.
            if chain_has_d and not d_raw:
                skipped["missing_d_call_on_vdj"] += 1
                rejected_entries.append(
                    {
                        "stage": "estimate_allele_usage",
                        "row_index": row_idx,
                        "reason": "missing_d_call_on_vdj",
                    }
                )
                continue
            if not chain_has_d and d_raw:
                if not d_on_vj_warned:
                    warnings.append(
                        "d_call column present on a VJ cartridge — D "
                        "contribution ignored for every row"
                    )
                    d_on_vj_warned = True
                d_raw = ""  # ignore D contribution silently after warning

            v_tie = _split_tie_set(v_raw)
            j_tie = _split_tie_set(j_raw)
            d_tie = _split_tie_set(d_raw) if d_raw else []

            if ambiguous == "reject":
                if len(v_tie) > 1 or len(j_tie) > 1 or (chain_has_d and len(d_tie) > 1):
                    skipped["ambiguous_rejected"] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_allele_usage",
                            "row_index": row_idx,
                            "reason": "ambiguous_rejected",
                        }
                    )
                    continue

            # truth_first collapses the tie set to its first entry.
            if ambiguous == "truth_first":
                v_tie = v_tie[:1]
                j_tie = j_tie[:1]
                d_tie = d_tie[:1] if d_tie else []

            # Resolve known alleles in each tie set.
            v_known = [n for n in v_tie if n in v_pool]
            j_known = [n for n in j_tie if n in j_pool]
            d_known = (
                [n for n in d_tie if n in d_pool] if d_tie else []
            )

            # Unknown-allele bookkeeping. We record one rejection
            # entry per unknown allele per segment per row so the
            # report names the actual unknown names.
            for n in v_tie:
                if n not in v_pool:
                    skipped["unknown_allele"]["V"] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_allele_usage",
                            "row_index": row_idx,
                            "segment": "V",
                            "allele_name": n,
                            "reason": "unknown_allele",
                        }
                    )
            for n in j_tie:
                if n not in j_pool:
                    skipped["unknown_allele"]["J"] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_allele_usage",
                            "row_index": row_idx,
                            "segment": "J",
                            "allele_name": n,
                            "reason": "unknown_allele",
                        }
                    )
            if d_tie:
                for n in d_tie:
                    if n not in d_pool:
                        skipped["unknown_allele"]["D"] += 1
                        rejected_entries.append(
                            {
                                "stage": "estimate_allele_usage",
                                "row_index": row_idx,
                                "segment": "D",
                                "allele_name": n,
                                "reason": "unknown_allele",
                            }
                        )

            # Fractional credit (the default and the truth_first
            # collapsed path; both use the same accumulation
            # logic now that truth_first reduced the tie set).
            if v_known:
                share = 1.0 / len(v_known)
                for n in v_known:
                    v_counts[n] = v_counts.get(n, 0.0) + share
            if j_known:
                share = 1.0 / len(j_known)
                for n in j_known:
                    j_counts[n] = j_counts.get(n, 0.0) + share
            if d_known:
                share = 1.0 / len(d_known)
                for n in d_known:
                    d_counts[n] = d_counts.get(n, 0.0) + share

        # min_count filter + normalisation per segment.
        below_min: Dict[str, int] = {"V": 0, "D": 0, "J": 0}

        def _filter_and_normalise(
            counts: Dict[str, float], segment_label: str
        ) -> Dict[str, float]:
            kept = {n: c for n, c in counts.items() if c >= min_count}
            dropped = len(counts) - len(kept)
            below_min[segment_label] = dropped
            total = sum(kept.values())
            if total <= 0.0:
                return {}
            return {n: w / total for n, w in kept.items()}

        v_weights = _filter_and_normalise(v_counts, "V")
        d_weights = _filter_and_normalise(d_counts, "D")
        j_weights = _filter_and_normalise(j_counts, "J")

        if below_min["V"] or below_min["D"] or below_min["J"]:
            warnings.append(
                f"alleles below min_count={min_count} dropped — "
                f"V={below_min['V']}, D={below_min['D']}, J={below_min['J']}"
            )

        spec = AlleleUsageSpec(v=v_weights, d=d_weights, j=j_weights)
        chain_label = "vdj" if chain_has_d else "vj"
        spec.validate(chain_type=chain_label, name="allele_usage")

        # Attach to the existing reference_models (creating it if
        # absent). The cartridge built downstream carries the
        # spec automatically.
        if self._reference_models is None:
            self._reference_models = ReferenceEmpiricalModels()
        self._reference_models = ReferenceEmpiricalModels(
            np_lengths=self._reference_models.np_lengths,
            trims=self._reference_models.trims,
            np_bases=self._reference_models.np_bases,
            p_nucleotide_lengths=self._reference_models.p_nucleotide_lengths,
            allele_usage=spec,
        )

        self._report.stages.append(
            {
                "stage": "estimate_allele_usage",
                "inputs": {
                    "record_count": len(records),
                    "ambiguous": ambiguous,
                    "min_count": float(min_count),
                    "source": source_label,
                    "replaced": previously_estimated,
                },
                "inferred": {
                    "V": v_weights,
                    "D": d_weights,
                    "J": j_weights,
                    "skipped": skipped,
                    "below_min_count": below_min,
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(rejected_entries)
        return self

    def estimate_trim_distributions(
        self,
        rearrangements: Any,
        *,
        min_count: int = 1,
        pseudocount: float = 0.0,
        replace: bool = True,
    ) -> "ReferenceCartridgeBuilder":
        """Estimate per-key trim distributions from observed
        AIRR rearrangement records.

        Writes ``EmpiricalDistributionSpec`` instances into
        ``self._reference_models.trims`` keyed by ``"V_3"`` /
        ``"D_5"`` / ``"D_3"`` / ``"J_5"`` (see
        :data:`GenAIRR.reference_models.TRIM_KEYS`). VJ cartridges
        only produce the ``V_3`` and ``J_5`` keys (see
        :data:`GenAIRR.reference_models.TRIM_KEYS_VJ`); the D-end
        keys are skipped silently because the cartridge has no D
        segment to trim.

        ``rearrangements`` accepts the same shapes as
        :meth:`estimate_allele_usage`:

        - a list of dicts (each row is one AIRR record),
        - a path-like to an AIRR TSV (parsed via ``csv.DictReader``
          with tab delimiter), or
        - an open text file handle pointing at AIRR TSV.

        AIRR columns consumed (per audit §1.3):

        - **VJ:** ``v_trim_3``, ``j_trim_5``.
        - **VDJ:** ``v_trim_3``, ``d_trim_5``, ``d_trim_3``, ``j_trim_5``.

        Columns deliberately ignored: ``v_trim_5`` / ``j_trim_3``
        (no engine pass — hard-zero in projection), and the
        observation-stage ``end_loss_5_length`` / ``end_loss_3_length``
        (separate corruption surface — see
        ``docs/primer_trim_end_loss_audit.md``).

        Per-row validation is **field-local**: a malformed or
        missing value in one trim column drops that field's
        contribution only — the row still feeds its other
        well-formed fields. Per-row structured entries land in
        ``report.rejected`` with reasons
        ``"missing_required_column"`` / ``"malformed_trim_value"``
        / ``"negative_trim_value"``, each carrying the AIRR
        column name.

        ``min_count`` (int, default 1) drops trim values whose
        observed integer count is strictly below the threshold
        before normalisation.

        ``pseudocount`` (float, default 0.0) adds a uniform prior
        to every **observed** trim value before normalisation
        (no support expansion; values never observed stay
        unobserved). Applied AFTER the ``min_count`` filter so
        the filter looks at raw observations.

        ``replace`` (default ``True``) controls idempotency. When
        ``True``, calling this method twice overwrites the
        previous specs and writes ``replaced=True`` on the new
        stage entry. When ``False`` AND a prior typed-plane
        ``trims`` is already attached to ``self._reference_models``,
        the call raises :class:`ValueError` before consuming any
        records.
        """
        if isinstance(min_count, bool) or not isinstance(min_count, int):
            raise ValueError(
                f"min_count must be an int, got {min_count!r}"
            )
        if min_count < 0:
            raise ValueError(
                f"min_count must be non-negative, got {min_count!r}"
            )
        if isinstance(pseudocount, bool) or not isinstance(
            pseudocount, (int, float)
        ):
            raise ValueError(
                f"pseudocount must be a non-negative number, got {pseudocount!r}"
            )
        if pseudocount < 0.0:
            raise ValueError(
                f"pseudocount must be non-negative, got {pseudocount!r}"
            )

        if not replace:
            existing_trims = (
                self._reference_models.trims
                if self._reference_models is not None
                else None
            )
            if existing_trims:
                raise ValueError(
                    "trim distributions already attached to this cartridge; "
                    "pass replace=True to overwrite"
                )
        previously_estimated = any(
            entry.get("stage") == "estimate_trim_distributions"
            for entry in self._report.stages
        )

        records, source_label = _load_rearrangements(rearrangements)
        chain_has_d = self._chain_type.has_d
        active_keys = TRIM_KEYS if chain_has_d else TRIM_KEYS_VJ

        # Map plane key → AIRR column.
        key_to_column = {
            "V_3": "v_trim_3",
            "D_5": "d_trim_5",
            "D_3": "d_trim_3",
            "J_5": "j_trim_5",
        }
        # AIRR columns the estimator deliberately does NOT consume.
        ignored_columns = ("v_trim_5", "j_trim_3")

        counters: Dict[str, Dict[int, int]] = {k: {} for k in active_keys}
        skipped = {
            "missing_required_column": {col: 0 for col in
                                        (key_to_column[k] for k in active_keys)},
            "malformed_trim_value": {col: 0 for col in
                                     (key_to_column[k] for k in active_keys)},
            "negative_trim_value": {col: 0 for col in
                                    (key_to_column[k] for k in active_keys)},
        }
        dropped_columns: Dict[str, int] = {col: 0 for col in ignored_columns}
        rejected_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []
        ignored_warned = {col: False for col in ignored_columns}
        vj_d_warned = False

        for row_idx, row in enumerate(records):
            # Surface non-zero `v_trim_5` / `j_trim_3` columns: track
            # the count but never consume them. One warning per
            # column across the dataset.
            for col in ignored_columns:
                raw = row.get(col)
                if raw not in (None, "", "0"):
                    try:
                        if int(str(raw).strip()) != 0:
                            dropped_columns[col] += 1
                            if not ignored_warned[col]:
                                warnings.append(
                                    f"{col} column non-zero in input — "
                                    f"contribution dropped (no V_5 / J_3 "
                                    f"trim pass in the engine)"
                                )
                                ignored_warned[col] = True
                    except (TypeError, ValueError):
                        # Malformed in an unused column is uninteresting.
                        pass

            # Warn once per dataset if a VJ cartridge sees populated
            # D-trim columns in the input. Same boundary as
            # `estimate_allele_usage`'s D-call ignore.
            if not chain_has_d:
                for col in ("d_trim_5", "d_trim_3"):
                    raw = row.get(col)
                    if raw not in (None, "", "0") and not vj_d_warned:
                        try:
                            if int(str(raw).strip()) != 0:
                                warnings.append(
                                    "VJ cartridge: d_trim_5 / d_trim_3 "
                                    "columns present in records; "
                                    "contribution ignored"
                                )
                                vj_d_warned = True
                                break
                        except (TypeError, ValueError):
                            pass

            # Field-local validation: each key/column independently.
            for key in active_keys:
                col = key_to_column[key]
                raw = row.get(col)
                if raw is None or str(raw).strip() == "":
                    skipped["missing_required_column"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_trim_distributions",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "reason": "missing_required_column",
                        }
                    )
                    continue
                try:
                    value = int(str(raw).strip())
                except (TypeError, ValueError):
                    skipped["malformed_trim_value"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_trim_distributions",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": raw,
                            "reason": "malformed_trim_value",
                        }
                    )
                    continue
                if value < 0:
                    skipped["negative_trim_value"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_trim_distributions",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": value,
                            "reason": "negative_trim_value",
                        }
                    )
                    continue
                counters[key][value] = counters[key].get(value, 0) + 1

        # min_count + pseudocount + per-key normalisation.
        below_min: Dict[str, int] = {k: 0 for k in active_keys}
        inferred_pairs: Dict[str, List[tuple]] = {k: [] for k in active_keys}

        for key in active_keys:
            raw_counts = counters[key]
            # Drop values strictly below min_count.
            kept = {v: c for v, c in raw_counts.items() if c >= min_count}
            below_min[key] = len(raw_counts) - len(kept)
            # Apply pseudocount to observed values only.
            if pseudocount > 0.0:
                kept = {v: (c + pseudocount) for v, c in kept.items()}
            # Normalise to sum 1.0 if anything survived.
            total = float(sum(kept.values()))
            if total <= 0.0 or not kept:
                inferred_pairs[key] = []
                continue
            normalised = [(v, kept[v] / total) for v in sorted(kept.keys())]
            inferred_pairs[key] = normalised

        if any(below_min[k] for k in active_keys):
            warnings.append(
                f"trim values below min_count={min_count} dropped — "
                + ", ".join(f"{k}={below_min[k]}" for k in active_keys)
            )

        # Build the per-key EmpiricalDistributionSpec instances.
        new_trims: Dict[str, EmpiricalDistributionSpec] = {}
        for key, pairs in inferred_pairs.items():
            if not pairs:
                continue
            spec = EmpiricalDistributionSpec(pairs)
            spec.validate(name=f"trims[{key}]")
            new_trims[key] = spec

        # Attach to the existing reference_models (creating it if
        # absent). Other typed planes are preserved.
        if self._reference_models is None:
            self._reference_models = ReferenceEmpiricalModels()
        self._reference_models = ReferenceEmpiricalModels(
            np_lengths=self._reference_models.np_lengths,
            trims=new_trims,
            np_bases=self._reference_models.np_bases,
            p_nucleotide_lengths=self._reference_models.p_nucleotide_lengths,
            allele_usage=self._reference_models.allele_usage,
        )
        # Validate the full container under the cartridge's chain
        # type so D-on-VJ etc. raise at attach time.
        chain_label = "vdj" if chain_has_d else "vj"
        self._reference_models.validate(chain_type=chain_label)

        self._report.stages.append(
            {
                "stage": "estimate_trim_distributions",
                "inputs": {
                    "record_count": len(records),
                    "min_count": int(min_count),
                    "pseudocount": float(pseudocount),
                    "source": source_label,
                    "replaced": previously_estimated,
                },
                "inferred": {
                    "V_3": inferred_pairs.get("V_3", []),
                    "D_5": inferred_pairs.get("D_5", []),
                    "D_3": inferred_pairs.get("D_3", []),
                    "J_5": inferred_pairs.get("J_5", []),
                    "skipped": skipped,
                    "below_min_count": {k: below_min.get(k, 0) for k in TRIM_KEYS},
                    "dropped_columns": dropped_columns,
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(rejected_entries)
        return self

    def estimate_np_length_distributions(
        self,
        rearrangements: Any,
        *,
        min_count: int = 1,
        pseudocount: float = 0.0,
        replace: bool = True,
    ) -> "ReferenceCartridgeBuilder":
        """Estimate per-key NP length distributions from
        observed AIRR rearrangement records.

        Writes ``EmpiricalDistributionSpec`` instances into
        ``self._reference_models.np_lengths`` keyed by
        ``"NP1"`` and (on VDJ cartridges) ``"NP2"`` (see
        :data:`GenAIRR.reference_models.NP_KEYS`). VJ
        cartridges produce only the ``NP1`` key; the
        cartridge has no NP2 region. The boundary is
        enforced at estimation time because the typed-plane
        validator does NOT chain-type-reject ``NP2`` on
        VJ at attach time (see audit §2.2 of
        ``docs/np_length_estimation_design.md``).

        AIRR columns consumed (per audit §1.2):

        - ``np1_length`` (always).
        - ``np2_length`` (VDJ only — VJ rows with a
          non-zero ``np2_length`` raise a one-time warning
          and the contribution is skipped).

        Columns deliberately ignored: ``np1`` / ``np2``
        (sequence-derived length is sensitive to
        post-claim reabsorption — see audit §5.3),
        ``p_v_3_length`` / ``p_d_5_length`` /
        ``p_d_3_length`` / ``p_j_5_length`` (separate
        P-nucleotide biology — see
        ``docs/p_nucleotide_design.md``),
        and ``junction_length`` (aggregate arithmetic too
        fragile across simulators — audit §5.2).

        Per-row validation is **field-local**: a malformed
        or missing value in one NP column drops only that
        key's contribution; the row's other column still
        feeds. Per-row structured entries land in
        ``report.rejected`` with reasons
        ``"missing_required_column"`` /
        ``"malformed_length_value"`` /
        ``"negative_length_value"``, each carrying the
        AIRR column name.

        ``min_count`` (int, default 1) drops length values
        whose observed integer count is strictly below
        the threshold before normalisation.

        ``pseudocount`` (float, default 0.0) adds a uniform
        prior to every **observed** length value before
        normalisation (no support expansion). Applied
        AFTER the ``min_count`` filter.

        ``replace`` (default ``True``) controls idempotency.
        When ``False`` AND a prior typed-plane
        ``np_lengths`` is attached to
        ``self._reference_models``, the call raises
        :class:`ValueError` before consuming any records.
        """
        if isinstance(min_count, bool) or not isinstance(min_count, int):
            raise ValueError(
                f"min_count must be an int, got {min_count!r}"
            )
        if min_count < 0:
            raise ValueError(
                f"min_count must be non-negative, got {min_count!r}"
            )
        if isinstance(pseudocount, bool) or not isinstance(
            pseudocount, (int, float)
        ):
            raise ValueError(
                f"pseudocount must be a non-negative number, got {pseudocount!r}"
            )
        if pseudocount < 0.0:
            raise ValueError(
                f"pseudocount must be non-negative, got {pseudocount!r}"
            )

        if not replace:
            existing = (
                self._reference_models.np_lengths
                if self._reference_models is not None
                else None
            )
            if existing:
                raise ValueError(
                    "np_length distributions already attached to this "
                    "cartridge; pass replace=True to overwrite"
                )
        previously_estimated = any(
            entry.get("stage") == "estimate_np_length_distributions"
            for entry in self._report.stages
        )

        records, source_label = _load_rearrangements(rearrangements)
        chain_has_d = self._chain_type.has_d
        active_keys = NP_KEYS if chain_has_d else ("NP1",)

        key_to_column = {"NP1": "np1_length", "NP2": "np2_length"}

        counters: Dict[str, Dict[int, int]] = {k: {} for k in active_keys}
        skipped = {
            "missing_required_column": {
                key_to_column[k]: 0 for k in active_keys
            },
            "malformed_length_value": {
                key_to_column[k]: 0 for k in active_keys
            },
            "negative_length_value": {
                key_to_column[k]: 0 for k in active_keys
            },
        }
        # VJ chains track `np2_length` contributions as a dropped
        # column with a single one-time warning.
        dropped_columns: Dict[str, int] = {}
        if not chain_has_d:
            dropped_columns["np2_length"] = 0
        rejected_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []
        vj_np2_warned = False

        for row_idx, row in enumerate(records):
            # On VJ, surface non-zero `np2_length` as a dropped
            # column with one warning across the dataset.
            if not chain_has_d:
                raw_np2 = row.get("np2_length")
                if raw_np2 not in (None, "", "0"):
                    try:
                        if int(str(raw_np2).strip()) != 0:
                            dropped_columns["np2_length"] += 1
                            if not vj_np2_warned:
                                warnings.append(
                                    "VJ cartridge: np2_length column "
                                    "present in records; contribution "
                                    "ignored (no NP2 region on a VJ "
                                    "chain)"
                                )
                                vj_np2_warned = True
                    except (TypeError, ValueError):
                        pass

            # Field-local validation per active key.
            for key in active_keys:
                col = key_to_column[key]
                raw = row.get(col)
                if raw is None or str(raw).strip() == "":
                    skipped["missing_required_column"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_np_length_distributions",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "reason": "missing_required_column",
                        }
                    )
                    continue
                try:
                    value = int(str(raw).strip())
                except (TypeError, ValueError):
                    skipped["malformed_length_value"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_np_length_distributions",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": raw,
                            "reason": "malformed_length_value",
                        }
                    )
                    continue
                if value < 0:
                    skipped["negative_length_value"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_np_length_distributions",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": value,
                            "reason": "negative_length_value",
                        }
                    )
                    continue
                counters[key][value] = counters[key].get(value, 0) + 1

        # min_count + pseudocount + per-key normalisation.
        below_min: Dict[str, int] = {k: 0 for k in active_keys}
        inferred_pairs: Dict[str, List[tuple]] = {k: [] for k in active_keys}

        for key in active_keys:
            raw_counts = counters[key]
            kept = {v: c for v, c in raw_counts.items() if c >= min_count}
            below_min[key] = len(raw_counts) - len(kept)
            if pseudocount > 0.0:
                kept = {v: (c + pseudocount) for v, c in kept.items()}
            total = float(sum(kept.values()))
            if total <= 0.0 or not kept:
                inferred_pairs[key] = []
                continue
            normalised = [(v, kept[v] / total) for v in sorted(kept.keys())]
            inferred_pairs[key] = normalised

        if any(below_min[k] for k in active_keys):
            warnings.append(
                f"NP-length values below min_count={min_count} dropped — "
                + ", ".join(f"{k}={below_min[k]}" for k in active_keys)
            )

        # Build per-key EmpiricalDistributionSpec instances.
        new_np_lengths: Dict[str, EmpiricalDistributionSpec] = {}
        for key, pairs in inferred_pairs.items():
            if not pairs:
                continue
            spec = EmpiricalDistributionSpec(pairs)
            spec.validate(name=f"np_lengths[{key}]")
            new_np_lengths[key] = spec

        # Attach to existing reference_models (creating if absent).
        # Other typed planes are preserved.
        if self._reference_models is None:
            self._reference_models = ReferenceEmpiricalModels()
        self._reference_models = ReferenceEmpiricalModels(
            np_lengths=new_np_lengths,
            trims=self._reference_models.trims,
            np_bases=self._reference_models.np_bases,
            p_nucleotide_lengths=self._reference_models.p_nucleotide_lengths,
            allele_usage=self._reference_models.allele_usage,
        )
        chain_label = "vdj" if chain_has_d else "vj"
        self._reference_models.validate(chain_type=chain_label)

        self._report.stages.append(
            {
                "stage": "estimate_np_length_distributions",
                "inputs": {
                    "record_count": len(records),
                    "min_count": int(min_count),
                    "pseudocount": float(pseudocount),
                    "source": source_label,
                    "replaced": previously_estimated,
                },
                "inferred": {
                    "NP1": inferred_pairs.get("NP1", []),
                    "NP2": inferred_pairs.get("NP2", []),
                    "skipped": skipped,
                    "below_min_count": {k: below_min.get(k, 0) for k in NP_KEYS},
                    "dropped_columns": dropped_columns,
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(rejected_entries)
        return self

    def estimate_np_base_model(
        self,
        rearrangements: Any,
        *,
        kind: str = "markov",
        min_count: int = 1,
        pseudocount: float = 0.0,
        replace: bool = True,
    ) -> "ReferenceCartridgeBuilder":
        """Estimate per-key NP base sampling model from
        observed AIRR rearrangement records.

        Writes ``NpBaseModelSpec`` instances into
        ``self._reference_models.np_bases`` keyed by
        ``"NP1"`` (always) and (on VDJ cartridges)
        ``"NP2"``. VJ cartridges produce only the ``NP1``
        key.

        ``kind`` selects the model family:

        - ``"empirical_first_base"``: estimates a single
          categorical over A/C/G/T from the **full base
          composition** of every observed NP string (every
          position, not just position 0). The engine
          samples every NP position independently from
          this distribution — the name is preserved for
          API stability with existing cartridges; the
          biologically correct estimate is the full-base
          composition.
        - ``"markov"`` (default): estimates a first-base
          row from position 0 of each NP string plus a
          4×4 transition matrix from every observed
          (prev, next) pair.

        AIRR columns consumed (per audit §1.2):

        - ``np1`` (always).
        - ``np2`` (VDJ only — VJ rows with a non-empty
          ``np2`` raise a one-time warning and the
          contribution is skipped).

        Columns deliberately ignored: ``junction`` (audit
        §1.2), ``p_v_3_length`` / ``p_d_5_length`` /
        ``p_d_3_length`` / ``p_j_5_length`` (separate
        P-nucleotide biology), ``np1_length`` /
        ``np2_length`` (length-only — owned by
        :meth:`estimate_np_length_distributions`).

        Per-row validation is **field-local**: a malformed
        or missing value in one NP column drops only that
        key's contribution. Per-row structured entries
        land in ``report.rejected`` with reasons
        ``"missing_required_column"`` (empty / missing
        string) or ``"noncanonical_base"`` (any character
        outside ``{A,C,G,T}`` after uppercasing).

        ``min_count`` (int, default 1) drops first-base
        categories whose observed count is strictly below
        the threshold before normalisation. For ``markov``,
        ``min_count`` applies to the **first-base row
        only**: dropping transition cells could leave a
        from-base row with no positive weights, which the
        :class:`NpBaseModelSpec` validator rejects. v1
        keeps transition rows intact and surfaces the
        first-base drops in ``below_min_count.first_base``.

        ``pseudocount`` (float, default 0.0) adds a uniform
        prior:

        - ``empirical_first_base``: added to every A/C/G/T
          base category before normalisation.
        - ``markov``: added to every A/C/G/T first-base
          category AND to every cell of the 4×4 transition
          matrix.

        ``replace`` (default ``True``) controls idempotency.
        When ``False`` AND a prior typed-plane
        ``np_bases`` is already attached, the call raises
        :class:`ValueError` before consuming any records.
        """
        if kind not in ("empirical_first_base", "markov"):
            raise ValueError(
                f"kind must be one of 'empirical_first_base' / "
                f"'markov', got {kind!r}"
            )
        if isinstance(min_count, bool) or not isinstance(min_count, int):
            raise ValueError(
                f"min_count must be an int, got {min_count!r}"
            )
        if min_count < 0:
            raise ValueError(
                f"min_count must be non-negative, got {min_count!r}"
            )
        if isinstance(pseudocount, bool) or not isinstance(
            pseudocount, (int, float)
        ):
            raise ValueError(
                f"pseudocount must be a non-negative number, got {pseudocount!r}"
            )
        if pseudocount < 0.0:
            raise ValueError(
                f"pseudocount must be non-negative, got {pseudocount!r}"
            )

        if not replace:
            existing = (
                self._reference_models.np_bases
                if self._reference_models is not None
                else None
            )
            if existing:
                raise ValueError(
                    "np_bases model already attached to this cartridge; "
                    "pass replace=True to overwrite"
                )
        previously_estimated = any(
            entry.get("stage") == "estimate_np_base_model"
            for entry in self._report.stages
        )

        records, source_label = _load_rearrangements(rearrangements)
        chain_has_d = self._chain_type.has_d
        active_keys = ("NP1", "NP2") if chain_has_d else ("NP1",)

        key_to_column = {"NP1": "np1", "NP2": "np2"}
        bases = _NP_CANONICAL_BASES  # ("A","C","G","T")

        # Per-key tallies. first_base[key] is a dict over A/C/G/T;
        # transitions[key] is a 4×4 dict-of-dict over A/C/G/T from→to.
        first_base_counts: Dict[str, Dict[str, int]] = {
            key: {b: 0 for b in bases} for key in active_keys
        }
        transition_counts: Dict[str, Dict[str, Dict[str, int]]] = {
            key: {b: {t: 0 for t in bases} for b in bases}
            for key in active_keys
        }
        skipped = {
            "missing_required_column": {
                key_to_column[k]: 0 for k in active_keys
            },
            "noncanonical_base": {
                key_to_column[k]: 0 for k in active_keys
            },
        }
        # VJ chains track non-empty `np2` strings as a dropped
        # column with a single one-time warning.
        dropped_columns: Dict[str, int] = {}
        if not chain_has_d:
            dropped_columns["np2"] = 0
        rejected_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []
        vj_np2_warned = False

        for row_idx, row in enumerate(records):
            # Surface non-empty `np2` on VJ as a dropped column
            # with one warning across the dataset.
            if not chain_has_d:
                raw_np2 = row.get("np2")
                if isinstance(raw_np2, str) and raw_np2.strip():
                    dropped_columns["np2"] += 1
                    if not vj_np2_warned:
                        warnings.append(
                            "VJ cartridge: np2 column non-empty in "
                            "records; contribution ignored (no NP2 "
                            "region on a VJ chain)"
                        )
                        vj_np2_warned = True

            # Field-local validation per active key.
            for key in active_keys:
                col = key_to_column[key]
                raw = row.get(col)
                if raw is None or not isinstance(raw, str) or not raw.strip():
                    skipped["missing_required_column"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_np_base_model",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "reason": "missing_required_column",
                        }
                    )
                    continue
                upper = raw.strip().upper()
                non_canonical = sorted(set(upper) - set(bases))
                if non_canonical:
                    skipped["noncanonical_base"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_np_base_model",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "unknown_chars": non_canonical,
                            "reason": "noncanonical_base",
                        }
                    )
                    continue
                # Tally bases.
                #   - For `empirical_first_base` the engine samples
                #     every NP position independently from the same
                #     distribution, so the biologically correct
                #     estimate is the FULL base composition (every
                #     position).
                #   - For `markov` the first-base row models position
                #     0 only (the seed of the chain); position 1+
                #     are modelled by transitions.
                if kind == "empirical_first_base":
                    for b in upper:
                        first_base_counts[key][b] += 1
                else:
                    first_base_counts[key][upper[0]] += 1
                for i in range(len(upper) - 1):
                    prev_b = upper[i]
                    next_b = upper[i + 1]
                    transition_counts[key][prev_b][next_b] += 1

        # min_count + pseudocount + per-kind spec construction.
        below_min: Dict[str, Dict[str, int]] = {
            key: {"first_base": 0} for key in active_keys
        }
        inferred: Dict[str, Dict[str, Any]] = {
            key: {} for key in active_keys
        }

        new_np_bases: Dict[str, NpBaseModelSpec] = {}
        from_base_errors: List[str] = []

        for key in active_keys:
            fb_counts = first_base_counts[key]
            # Apply min_count filter to first-base row (drop bases
            # whose count is strictly below threshold).
            fb_kept = {b: c for b, c in fb_counts.items() if c >= min_count}
            below_min[key]["first_base"] = len(fb_counts) - len(fb_kept)
            # Apply pseudocount AFTER filter — observed-only support.
            if pseudocount > 0.0:
                for b in bases:
                    if b in fb_kept:
                        fb_kept[b] = fb_kept[b] + pseudocount
                    else:
                        fb_kept[b] = pseudocount
            fb_total = float(sum(fb_kept.values()))
            if fb_total <= 0.0:
                # Nothing observed AND no pseudocount → skip this key.
                inferred[key] = {"first_base": None, "transitions": None}
                continue
            first_base_weights = {
                b: fb_kept[b] / fb_total
                for b in sorted(fb_kept.keys())
                if fb_kept[b] > 0
            }

            spec_payload: Dict[str, Any] = {
                "first_base": first_base_weights,
            }
            if kind == "markov":
                # Apply pseudocount to every cell of every transition row.
                tr = transition_counts[key]
                tr_weights: Dict[str, Dict[str, float]] = {}
                for from_b in bases:
                    row_counts = {t: float(tr[from_b][t]) for t in bases}
                    if pseudocount > 0.0:
                        for t in bases:
                            row_counts[t] += pseudocount
                    row_total = sum(row_counts.values())
                    if row_total <= 0.0:
                        # Pseudocount=0 and from-base never observed: surface a
                        # tagged error so the caller knows to pass pseudocount.
                        from_base_errors.append(
                            f"{key}: transition row for from-base "
                            f"{from_b!r} has no observed transitions "
                            f"and pseudocount=0; pass pseudocount > 0 "
                            f"or supply more data"
                        )
                        continue
                    tr_weights[from_b] = {
                        t: row_counts[t] / row_total
                        for t in bases
                        if row_counts[t] > 0
                    }
                if from_base_errors:
                    # Bail before constructing a partial spec.
                    raise ValueError("; ".join(from_base_errors))
                spec_payload["transitions"] = tr_weights
            else:
                spec_payload["transitions"] = None

            # Construct + validate the spec.
            spec = NpBaseModelSpec(
                kind=kind,
                first_base=spec_payload["first_base"],
                transitions=spec_payload["transitions"],
            )
            spec.validate(name=f"np_bases[{key}]")
            new_np_bases[key] = spec
            inferred[key] = {
                "first_base": dict(spec_payload["first_base"]),
                "transitions": (
                    {k: dict(v) for k, v in spec_payload["transitions"].items()}
                    if spec_payload["transitions"] is not None
                    else None
                ),
            }

        if any(below_min[k]["first_base"] for k in active_keys):
            warnings.append(
                f"first-base values below min_count={min_count} dropped — "
                + ", ".join(
                    f"{k}={below_min[k]['first_base']}"
                    for k in active_keys
                )
            )

        # Attach to existing reference_models (creating if absent);
        # other typed planes preserved.
        if self._reference_models is None:
            self._reference_models = ReferenceEmpiricalModels()
        self._reference_models = ReferenceEmpiricalModels(
            np_lengths=self._reference_models.np_lengths,
            trims=self._reference_models.trims,
            np_bases=new_np_bases,
            p_nucleotide_lengths=self._reference_models.p_nucleotide_lengths,
            allele_usage=self._reference_models.allele_usage,
        )
        chain_label = "vdj" if chain_has_d else "vj"
        self._reference_models.validate(chain_type=chain_label)

        # NP2 entry surfaces in inferred even on VJ (as empty dict)
        # for shape stability — same discipline as the NP-length
        # estimator's empty `D_5` / `D_3` entries.
        np1_inferred = inferred.get("NP1", {"first_base": None, "transitions": None})
        np2_inferred = inferred.get("NP2", {"first_base": None, "transitions": None})

        self._report.stages.append(
            {
                "stage": "estimate_np_base_model",
                "inputs": {
                    "record_count": len(records),
                    "kind": kind,
                    "min_count": int(min_count),
                    "pseudocount": float(pseudocount),
                    "source": source_label,
                    "replaced": previously_estimated,
                },
                "inferred": {
                    "NP1": np1_inferred,
                    "NP2": np2_inferred,
                    "skipped": skipped,
                    "below_min_count": {
                        k: below_min.get(k, {"first_base": 0})
                        for k in NP_KEYS
                    },
                    "dropped_columns": dropped_columns,
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(rejected_entries)
        return self

    def estimate_p_nucleotide_lengths(
        self,
        rearrangements: Any,
        *,
        min_count: int = 1,
        pseudocount: float = 0.0,
        replace: bool = True,
    ) -> "ReferenceCartridgeBuilder":
        """Estimate per-end P-nucleotide length distributions
        from observed AIRR rearrangement records.

        Writes ``EmpiricalDistributionSpec`` instances into
        ``self._reference_models.p_nucleotide_lengths``
        keyed by ``"V_3"`` and ``"J_5"`` (always) plus
        ``"D_5"`` and ``"D_3"`` (VDJ only). VJ cartridges
        produce only the V_3 and J_5 keys; the typed-plane
        validator chain-type-rejects D-end keys on VJ at
        attach time.

        **Provenance warning.** This estimator requires
        AIRR-like records that ALREADY carry GenAIRR's
        P-length fields (``p_v_3_length`` /
        ``p_d_5_length`` / ``p_d_3_length`` /
        ``p_j_5_length``). It does NOT infer P-lengths
        from generic AIRR junction sequences, NP strings,
        or trim arithmetic — that inference problem is
        out of scope for v1. External AIRR tools (IgBLAST,
        MiXCR, …) do not model P-nucleotide additions, so
        their output will either omit these columns
        entirely (rejection storm) or populate them as
        zero (degenerate ``[(0, 1.0)]`` distribution).
        The estimator emits a stage-level warning per key
        when ≥ 95% of contributing rows reported zero —
        diagnostic of P-naïve input. See
        ``docs/p_nucleotide_length_estimation_design.md``
        §5 for the realistic-input enumeration.

        AIRR columns consumed:

        - ``p_v_3_length`` (always).
        - ``p_j_5_length`` (always).
        - ``p_d_5_length`` (VDJ only — VJ rows with a
          non-zero value raise a one-time warning per
          column and the contribution is skipped).
        - ``p_d_3_length`` (same).

        Columns deliberately ignored: every other AIRR
        column. v1 does NOT derive P-lengths from
        junction-arithmetic, NP strings, trim fields, or
        end-loss fields.

        Per-row validation is **field-local**: a malformed
        or missing value in one P-length column drops only
        that key's contribution. Per-row structured
        entries land in ``report.rejected`` with reasons
        ``"missing_required_column"`` /
        ``"malformed_length_value"`` /
        ``"negative_length_value"``.

        ``min_count`` (int, default 1) drops length values
        whose observed count is strictly below the
        threshold before normalisation.

        ``pseudocount`` (float, default 0.0) adds a uniform
        prior to every **observed** length value before
        normalisation (no support expansion). Applied
        AFTER the ``min_count`` filter.

        ``replace`` (default ``True``) controls idempotency.
        When ``False`` AND a prior typed-plane
        ``p_nucleotide_lengths`` is already attached,
        the call raises :class:`ValueError` before
        consuming any records.
        """
        if isinstance(min_count, bool) or not isinstance(min_count, int):
            raise ValueError(
                f"min_count must be an int, got {min_count!r}"
            )
        if min_count < 0:
            raise ValueError(
                f"min_count must be non-negative, got {min_count!r}"
            )
        if isinstance(pseudocount, bool) or not isinstance(
            pseudocount, (int, float)
        ):
            raise ValueError(
                f"pseudocount must be a non-negative number, got {pseudocount!r}"
            )
        if pseudocount < 0.0:
            raise ValueError(
                f"pseudocount must be non-negative, got {pseudocount!r}"
            )

        if not replace:
            existing = (
                self._reference_models.p_nucleotide_lengths
                if self._reference_models is not None
                else None
            )
            if existing:
                raise ValueError(
                    "p_nucleotide_lengths model already attached to "
                    "this cartridge; pass replace=True to overwrite"
                )
        previously_estimated = any(
            entry.get("stage") == "estimate_p_nucleotide_lengths"
            for entry in self._report.stages
        )

        records, source_label = _load_rearrangements(rearrangements)
        chain_has_d = self._chain_type.has_d
        active_keys = (
            P_NUCLEOTIDE_END_KEYS if chain_has_d
            else P_NUCLEOTIDE_END_KEYS_VJ
        )

        key_to_column = {
            "V_3": "p_v_3_length",
            "D_5": "p_d_5_length",
            "D_3": "p_d_3_length",
            "J_5": "p_j_5_length",
        }
        vj_ignored_columns = ("p_d_5_length", "p_d_3_length")

        counters: Dict[str, Dict[int, int]] = {k: {} for k in active_keys}
        contributing_counts: Dict[str, int] = {k: 0 for k in active_keys}
        zero_counts: Dict[str, int] = {k: 0 for k in active_keys}
        skipped = {
            "missing_required_column": {
                key_to_column[k]: 0 for k in active_keys
            },
            "malformed_length_value": {
                key_to_column[k]: 0 for k in active_keys
            },
            "negative_length_value": {
                key_to_column[k]: 0 for k in active_keys
            },
        }
        # VJ chains track nonzero `p_d_*_length` columns as
        # dropped columns with one warning per column.
        dropped_columns: Dict[str, int] = {}
        if not chain_has_d:
            for col in vj_ignored_columns:
                dropped_columns[col] = 0
        rejected_entries: List[Dict[str, Any]] = []
        warnings: List[str] = []
        vj_ignored_warned = {col: False for col in vj_ignored_columns}

        for row_idx, row in enumerate(records):
            # VJ: surface nonzero D-end columns as dropped + warn once per column.
            if not chain_has_d:
                for col in vj_ignored_columns:
                    raw = row.get(col)
                    if raw in (None, "", "0"):
                        continue
                    try:
                        value = int(str(raw).strip())
                    except (TypeError, ValueError):
                        continue
                    if value != 0:
                        dropped_columns[col] += 1
                        if not vj_ignored_warned[col]:
                            warnings.append(
                                f"VJ cartridge: {col} non-zero in "
                                f"input — contribution ignored (no "
                                f"D segment on a VJ chain)"
                            )
                            vj_ignored_warned[col] = True

            # Field-local validation per active key.
            for key in active_keys:
                col = key_to_column[key]
                raw = row.get(col)
                if raw is None or str(raw).strip() == "":
                    skipped["missing_required_column"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_p_nucleotide_lengths",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "reason": "missing_required_column",
                        }
                    )
                    continue
                try:
                    value = int(str(raw).strip())
                except (TypeError, ValueError):
                    skipped["malformed_length_value"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_p_nucleotide_lengths",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": raw,
                            "reason": "malformed_length_value",
                        }
                    )
                    continue
                if value < 0:
                    skipped["negative_length_value"][col] += 1
                    rejected_entries.append(
                        {
                            "stage": "estimate_p_nucleotide_lengths",
                            "row_index": row_idx,
                            "column": col,
                            "key": key,
                            "value": value,
                            "reason": "negative_length_value",
                        }
                    )
                    continue
                counters[key][value] = counters[key].get(value, 0) + 1
                contributing_counts[key] += 1
                if value == 0:
                    zero_counts[key] += 1

        # min_count + pseudocount + per-key normalisation.
        below_min: Dict[str, int] = {k: 0 for k in active_keys}
        inferred_pairs: Dict[str, List[tuple]] = {k: [] for k in active_keys}
        zero_fraction: Dict[str, float] = {}

        for key in active_keys:
            raw_counts = counters[key]
            kept = {v: c for v, c in raw_counts.items() if c >= min_count}
            below_min[key] = len(raw_counts) - len(kept)
            if pseudocount > 0.0:
                kept = {v: (c + pseudocount) for v, c in kept.items()}
            total = float(sum(kept.values()))
            if total <= 0.0 or not kept:
                inferred_pairs[key] = []
            else:
                inferred_pairs[key] = [
                    (v, kept[v] / total) for v in sorted(kept.keys())
                ]
            # zero_fraction tracks observed-row provenance before
            # normalisation — useful diagnostic for P-naïve input.
            n_contrib = contributing_counts[key]
            if n_contrib > 0:
                zero_fraction[key] = zero_counts[key] / n_contrib
            else:
                zero_fraction[key] = 0.0

        if any(below_min[k] for k in active_keys):
            warnings.append(
                f"P-nucleotide length values below min_count={min_count} "
                f"dropped — "
                + ", ".join(f"{k}={below_min[k]}" for k in active_keys)
            )

        # Per-key provenance auto-warning (audit §5.2).
        for key in active_keys:
            if (
                contributing_counts[key] > 0
                and zero_fraction[key] >= 0.95
            ):
                warnings.append(
                    f"p_nucleotide_lengths[{key}] is >=95% zero; "
                    f"input may be P-naive or lack P annotations"
                )

        # Build per-key EmpiricalDistributionSpec instances.
        new_p_lengths: Dict[str, EmpiricalDistributionSpec] = {}
        for key, pairs in inferred_pairs.items():
            if not pairs:
                continue
            spec = EmpiricalDistributionSpec(pairs)
            spec.validate(name=f"p_nucleotide_lengths[{key}]")
            new_p_lengths[key] = spec

        # Attach to existing reference_models (creating if absent);
        # other typed planes preserved.
        if self._reference_models is None:
            self._reference_models = ReferenceEmpiricalModels()
        self._reference_models = ReferenceEmpiricalModels(
            np_lengths=self._reference_models.np_lengths,
            trims=self._reference_models.trims,
            np_bases=self._reference_models.np_bases,
            p_nucleotide_lengths=new_p_lengths,
            allele_usage=self._reference_models.allele_usage,
        )
        chain_label = "vdj" if chain_has_d else "vj"
        self._reference_models.validate(chain_type=chain_label)

        self._report.stages.append(
            {
                "stage": "estimate_p_nucleotide_lengths",
                "inputs": {
                    "record_count": len(records),
                    "min_count": int(min_count),
                    "pseudocount": float(pseudocount),
                    "source": source_label,
                    "replaced": previously_estimated,
                },
                "inferred": {
                    "V_3": inferred_pairs.get("V_3", []),
                    "D_5": inferred_pairs.get("D_5", []),
                    "D_3": inferred_pairs.get("D_3", []),
                    "J_5": inferred_pairs.get("J_5", []),
                    "skipped": skipped,
                    "below_min_count": {
                        k: below_min.get(k, 0) for k in P_NUCLEOTIDE_END_KEYS
                    },
                    "zero_fraction": {
                        k: zero_fraction.get(k, 0.0)
                        for k in P_NUCLEOTIDE_END_KEYS
                    },
                    "dropped_columns": dropped_columns,
                },
                "warnings": warnings,
            }
        )
        self._report.rejected.extend(rejected_entries)
        return self

    # ──────────────────────────────────────────────────────────
    # Finalisation
    # ──────────────────────────────────────────────────────────

    def build(self) -> DataConfig:
        """Assemble and return a validated :class:`DataConfig`.

        The cartridge is finalised with:

        - allele dicts populated from FASTA.
        - ``metadata`` from :meth:`infer_identity` (or a built-in
          stub when :meth:`infer_identity` was never called).
        - ``reference_rules`` / ``reference_models`` when those
          stages ran.
        - ``build_report`` populated with the audit trail + a
          manifest snapshot + the post-build checksum.

        Final step calls :meth:`DataConfig.verify_integrity` so a
        malformed cartridge surfaces at build time, not at
        :func:`Experiment.on` time. The build report is attached
        before integrity check so :attr:`build_report` is
        consistent regardless of whether the check passes.
        """
        if not self._v_alleles:
            raise ValueError(
                "build(): no V alleles parsed — call from_fasta(...) first"
            )
        if not self._j_alleles:
            raise ValueError(
                "build(): no J alleles parsed — call from_fasta(...) first"
            )
        if self._chain_type.has_d and not self._d_alleles:
            raise ValueError(
                f"build(): chain_type {self._chain_type.name} requires "
                f"D alleles, but none were parsed — verify d_fasta input"
            )

        # Fill in a default ConfigInfo if infer_identity was never
        # called. The default uses HUMAN / "user-supplied" so the
        # cartridge survives downstream consumers that read
        # `metadata.species` without an `is None` guard. The
        # build-report stage entry surfaces the omission.
        if self._metadata is None:
            self._metadata = ConfigInfo(
                species=Species.HUMAN,
                chain_type=self._chain_type,
                reference_set="user-supplied",
                last_updated=date.today(),
                has_d=self._chain_type.has_d,
            )
            self._report.warnings.append(
                "infer_identity() was not called — cartridge metadata "
                "defaults to HUMAN species + 'user-supplied' reference_set"
            )

        cfg = DataConfig(
            name=self._name or f"USER_{self._chain_type.name}",
            metadata=self._metadata,
            v_alleles=self._v_alleles,
            d_alleles=self._d_alleles or None,
            j_alleles=self._j_alleles,
            c_alleles=None,  # v1 boundary
            reference_rules=self._reference_rules,
            reference_models=self._reference_models,
            genotype_priors=self._genotype_priors,
        )
        # Pull a manifest snapshot + checksum. The manifest call
        # runs before verify_integrity so the report carries the
        # final state even when integrity blows up.
        try:
            manifest = cfg.cartridge_manifest()
        except Exception as exc:  # pragma: no cover — defensive
            manifest = {"error": f"manifest_failed: {exc!s}"}
        self._report.manifest_snapshot = manifest

        # Stamp the canonical checksum onto the cartridge so
        # `verify_integrity` succeeds. `compute_checksum` is
        # transient-surgery-safe and tolerates the unfilled
        # `schema_sha256` field (it zeroes the field for the
        # hash computation regardless of its prior value).
        try:
            cfg.schema_sha256 = cfg.compute_checksum()
            self._report.checksum_at_build_time = cfg.schema_sha256
        except Exception as exc:  # pragma: no cover — defensive
            self._report.checksum_at_build_time = None
            self._report.warnings.append(
                f"compute_checksum failed at build-time: {exc!s}"
            )

        self._report.stages.append(
            {
                "stage": "build",
                "inputs": {},
                "inferred": {
                    "name": cfg.name,
                    "checksum": self._report.checksum_at_build_time,
                    "manifest_snapshot_attached": True,
                },
                "warnings": [],
            }
        )
        cfg.build_report = self._report

        # Final gate — surfaces malformed cartridges before any
        # downstream consumer (Experiment.on, pickle persist, etc.).
        cfg.verify_integrity()
        return cfg

    def report(self) -> CartridgeBuildReport:
        """Return the accumulated :class:`CartridgeBuildReport`."""
        return self._report


# ──────────────────────────────────────────────────────────────────
# Private helpers
# ──────────────────────────────────────────────────────────────────


def _parse_segment_fasta(
    source: FastaInput,
    *,
    segment: str,
    builder: ReferenceCartridgeBuilder,
    allele_cls: Any,
) -> "tuple[int, list[dict]]":
    """Parse one FASTA file into the builder's per-segment bucket.

    Returns ``(parsed_count, rejected_list)``. Rejected entries are
    dicts ``{stage, segment, allele_name, reason}`` ready for the
    report.

    Each parsed sequence is normalised to upper-case; ``.`` gaps
    are preserved so V alleles can later route through
    :meth:`infer_v_subregions`. Allele names duplicate-protect:
    the second occurrence of an exact name is rejected with the
    reason ``"duplicate_name"``.
    """
    file_handle, close_on_exit = _open_fasta(source)
    rejected: List[Dict[str, Any]] = []
    seen_names: set = set()
    parsed = 0
    try:
        for raw_header, raw_seq in parse_fasta(file_handle):
            name = _parse_allele_name(raw_header)
            if not name:
                rejected.append(
                    {
                        "stage": "from_fasta",
                        "segment": segment,
                        "allele_name": None,
                        "reason": "empty_or_unparseable_header",
                    }
                )
                continue
            if name in seen_names:
                rejected.append(
                    {
                        "stage": "from_fasta",
                        "segment": segment,
                        "allele_name": name,
                        "reason": "duplicate_name",
                    }
                )
                continue
            seen_names.add(name)
            gapped = raw_seq.strip().upper()
            if not gapped:
                rejected.append(
                    {
                        "stage": "from_fasta",
                        "segment": segment,
                        "allele_name": name,
                        "reason": "empty_sequence",
                    }
                )
                continue
            try:
                # `length` is the gapped length per the Allele
                # contract; ungapped length is derived inside
                # `Allele.__init__` from the gapped sequence.
                allele = allele_cls(
                    name=name, gapped_sequence=gapped, length=len(gapped)
                )
            except Exception as exc:
                rejected.append(
                    {
                        "stage": "from_fasta",
                        "segment": segment,
                        "allele_name": name,
                        "reason": f"allele_constructor_failed: {exc!s}",
                    }
                )
                continue
            bucket = {
                "V": builder._v_alleles,
                "D": builder._d_alleles,
                "J": builder._j_alleles,
            }[segment]
            bucket.setdefault(allele.gene, []).append(allele)
            parsed += 1
    finally:
        if close_on_exit:
            file_handle.close()
    return parsed, rejected
