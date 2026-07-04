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
    ReferenceEmpiricalModels,
)
from .reference_rules import ReferenceRulesSpec
from .utilities.imgt_regions import compute_v_region_boundaries
from .utilities.misc import parse_fasta
from ._cartridge_estimators import _CartridgeEstimators


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


class ReferenceCartridgeBuilder(_CartridgeEstimators):
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
