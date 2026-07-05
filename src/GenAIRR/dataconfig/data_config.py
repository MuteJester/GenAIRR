from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional, Dict, List, Any

from GenAIRR.dataconfig.cartridge_views import (
    CartridgeCatalogueView,
    CartridgeIdentityView,
    CartridgeModelsView,
    CartridgeRulesView,
)
from GenAIRR.dataconfig.config_info import ConfigInfo
from GenAIRR.dataconfig.enums import ChainType
from GenAIRR.dataconfig._manifest import build_manifest
# Backward-compat re-export: these two documented-gap tuples moved to
# ``_manifest`` but are still imported from ``data_config`` by callers
# (e.g. tests/test_cartridge_manifest.py). Re-exported to preserve the
# import surface without duplicating the definitions.
from GenAIRR.dataconfig._manifest import (
    _DOCUMENTED_DROPPED_ALLELE_FIELDS,
    _DOCUMENTED_ORPHAN_DATACONFIG_FIELDS,
)
import copy
import hashlib
import pickle
from GenAIRR.alleles.allele import Allele
if TYPE_CHECKING:
    # Authoring-layer types used only as field annotations. Deferred to
    # TYPE_CHECKING so loading a pickled DataConfig (every bundled
    # cartridge is one) does not eagerly force-import the entire
    # cartridge-authoring subtree.
    from GenAIRR.reference_models import ReferenceEmpiricalModels
    from GenAIRR.reference_rules import ReferenceRulesSpec
    from GenAIRR.genotype_priors import PopulationGenotypeModel


DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS = {0: 0.50, 1: 0.25, 2: 0.15, 3: 0.07, 4: 0.03}

# Bump on any breaking change to DataConfig fields. Old pickles without
# this version (or with an older one) refuse to load and prompt the user
# to re-migrate. See _check_schema() / verify_integrity() below.
SCHEMA_VERSION = 1


class DataConfigError(Exception):
    """Raised when a DataConfig fails validation."""


@dataclass
class DataConfig:
    """
    Configuration class for storing data related to sequence generation, allele usage, trimming, and mutation rates.

    This class encapsulates various dictionaries and settings used in the simulation and analysis of immunoglobulin sequences.
    """
    # --- Attributes ---
    name: Optional[str] = None
    metadata: Optional[ConfigInfo] = None

    # Gene usage frequencies
    gene_use_dict: Dict[str, Any] = field(default_factory=dict)

    # Allele dictionaries (grouped by gene name → list of allele objects)
    v_alleles: Optional[Dict[str, List[Allele]]] = None
    d_alleles: Optional[Dict[str, List[Allele]]] = None
    j_alleles: Optional[Dict[str, List[Allele]]] = None
    c_alleles: Optional[Dict[str, List[Allele]]] = None

    # Trimming and NP region parameters
    trim_dicts: Dict[str, Any] = field(default_factory=dict)
    NP_transitions: Dict[str, Any] = field(default_factory=dict)
    NP_first_bases: Dict[str, Any] = field(default_factory=dict)
    NP_lengths: Dict[str, Any] = field(default_factory=dict)

    # Correction maps and ASC tables
    correction_maps: Dict[str, Any] = field(default_factory=dict)
    asc_tables: Dict[str, Any] = field(default_factory=dict)

    # P-nucleotide length distribution (0-4 bp, geometric decay)
    p_nucleotide_length_probs: Dict[int, float] = field(
        default_factory=lambda: dict(DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS)
    )

    # D-J pairing map: {D_gene_name: [compatible_J_gene_name, ...]}
    # When None, no D-J family pairing constraint is applied.
    dj_pairing_map: Optional[Dict[str, list]] = None

    # Schema version + integrity checksum (T2-1). Populated by the migration
    # script for all builtins; pickles without these fields fail load-time
    # verification with a clear migration hint. schema_sha256 is sha256 over
    # pickle.dumps(self, protocol=4) with schema_sha256 set to "" — see
    # compute_checksum() below.
    schema_version: int = SCHEMA_VERSION
    schema_sha256: str = ""

    # Build-time provenance. Populated by
    # :class:`GenAIRR.cartridge_builder.ReferenceCartridgeBuilder.build`
    # with a :class:`~GenAIRR.cartridge_builder.CartridgeBuildReport`
    # carrying per-stage inputs / inferred / warnings / rejected
    # entries, a manifest snapshot, and the post-build checksum.
    # ``None`` for the 106 bundled cartridges (predate the builder)
    # and for manual `DataConfig(...)` constructions; downstream
    # consumers can treat absence as "legacy / unaudited cartridge".
    build_report: Optional[Any] = None

    # Optional reference rules spec — programmable interpretation layer
    # (allowed alphabet, V/J anchor expectations + severities). When
    # set, transferred verbatim into ``RefDataConfig.rules`` by
    # ``dataconfig_to_refdata``; when ``None``, the loader falls back
    # to the bundled locus-derived defaults.
    #
    # Soft-transition checksum policy: ``None`` is excluded from the
    # checksum so legacy pickles (which lack this key in ``__dict__``)
    # continue to verify. Any non-``None`` value is folded into the
    # checksum because it materially changes simulation semantics.
    reference_rules: Optional[ReferenceRulesSpec] = None

    # Optional empirical-models bundle — typed NP-length and trim
    # distributions that override the legacy nested-dict extraction
    # path. When set, ``recombine()`` consumes these directly; when
    # ``None``, the loader falls back to the legacy
    # ``NP_lengths`` / ``trim_dicts`` extraction and finally the
    # uniform placeholder. Same soft-transition checksum policy as
    # ``reference_rules``.
    reference_models: Optional[ReferenceEmpiricalModels] = None

    # Donor-population germline prior plane (Slice — Cartridge genotype plane).
    # ``None`` means no population prior; ``Genotype.sample(cfg)`` then falls
    # back to a uniform synthetic prior. A non-``None`` plane is cartridge
    # identity (folds into compute_checksum). See
    # ``site_docs/guides/genotype.md`` ("Population genotype models on a
    # cartridge"). Same soft-transition checksum policy as ``reference_rules`` /
    # ``reference_models``.
    genotype_priors: Optional[PopulationGenotypeModel] = None

    def __getattr__(self, name):
        # Backward-compat shim for pickled DataConfigs missing post-v1
        # fields. Note: schema_version / schema_sha256 fall through to
        # the dataclass class-level defaults (SCHEMA_VERSION and "")
        # when the instance __dict__ doesn't have them, so they don't
        # need entries here — verify_integrity() then sees version=1
        # but checksum="", which trips the missing-checksum branch and
        # produces a clean migration error.
        if name == 'p_nucleotide_length_probs':
            return dict(DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS)
        if name == 'dj_pairing_map':
            return None
        if name == 'build_report':
            return None
        if name == 'reference_rules':
            return None
        if name == 'reference_models':
            return None
        if name == 'genotype_priors':
            return None
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def compute_checksum(self) -> str:
        """Compute the canonical sha256 of this DataConfig's contents.

        The hash is taken over ``pickle.dumps(self, protocol=4)`` with
        the following transient surgery on ``__dict__``:

        - ``schema_sha256`` is temporarily zeroed so the field cannot
          include itself in its own hash.
        - ``build_report`` is removed entirely — it's diagnostic
          provenance with no semantic effect, and legacy pickles
          predating the field never had it in ``__dict__`` at all.
          Removing it keeps the wire shape identical for both lineages.
        - ``reference_rules`` is removed entirely **only when its value
          is None**. Legacy pickles don't have the key; new
          ``DataConfig()`` instances do (default ``None`` populates
          ``__dict__``). Popping the None case lets legacy and default-
          new pickles share a checksum. Any non-``None`` value is
          retained because it materially changes simulation semantics
          and should produce a distinct hash.
        """
        saved_sha = self.schema_sha256
        had_report = 'build_report' in self.__dict__
        saved_report = self.__dict__.get('build_report')
        # Soft-transition shim for post-v1 optional fields. For each
        # field we ALSO pop only when its value is ``None`` so legacy
        # pickles (which never had the key in __dict__) and default-
        # new instances (which do have the key set to None) produce
        # the same checksum. Non-``None`` values are retained because
        # they materially change simulation semantics and should be
        # part of cartridge identity.
        rr_value = self.__dict__.get('reference_rules')
        pop_rules = 'reference_rules' in self.__dict__ and rr_value is None
        rm_value = self.__dict__.get('reference_models')
        pop_models = 'reference_models' in self.__dict__ and rm_value is None
        gp_value = self.__dict__.get('genotype_priors')
        pop_priors = 'genotype_priors' in self.__dict__ and gp_value is None

        self.schema_sha256 = ""
        if had_report:
            del self.__dict__['build_report']
        if pop_rules:
            del self.__dict__['reference_rules']
        if pop_models:
            del self.__dict__['reference_models']
        if pop_priors:
            del self.__dict__['genotype_priors']
        try:
            blob = pickle.dumps(self, protocol=4)
            return hashlib.sha256(blob).hexdigest()
        finally:
            self.schema_sha256 = saved_sha
            if had_report:
                self.__dict__['build_report'] = saved_report
            if pop_rules:
                self.__dict__['reference_rules'] = rr_value
            if pop_models:
                self.__dict__['reference_models'] = rm_value
            if pop_priors:
                self.__dict__['genotype_priors'] = gp_value

    def verify_integrity(self) -> None:
        """Validate schema_version and schema_sha256.

        Raises ``DataConfigError`` if the pickle predates the current
        schema version or if its stored checksum does not match the
        recomputed checksum (corrupted/tampered file).
        """
        version = getattr(self, "schema_version", 0)
        if version != SCHEMA_VERSION:
            raise DataConfigError(
                f"DataConfig '{self.name or 'Unnamed'}' has schema_version="
                f"{version}, but the current code expects "
                f"{SCHEMA_VERSION}. Re-migrate or rebuild this pickle "
                f"from the IMGT source via dataconfig builders."
            )
        stored = getattr(self, "schema_sha256", "") or ""
        if not stored:
            raise DataConfigError(
                f"DataConfig '{self.name or 'Unnamed'}' is missing "
                f"schema_sha256. Re-migrate this pickle to populate the "
                f"checksum field."
            )
        actual = self.compute_checksum()
        if actual != stored:
            raise DataConfigError(
                f"DataConfig '{self.name or 'Unnamed'}' checksum mismatch: "
                f"stored={stored[:16]}…, actual={actual[:16]}…. "
                f"The pickle has been modified or corrupted since it was "
                f"saved. Reinstall GenAIRR or rebuild this config from "
                f"its IMGT source."
            )

    def validate(self):
        """
        Validate that this DataConfig has the minimum required fields for simulation.

        Raises:
            DataConfigError: If any required field is missing or malformed.
        """
        errors = []

        # Metadata
        if self.metadata is not None and not isinstance(self.metadata, ConfigInfo):
            errors.append(f"metadata must be a ConfigInfo instance, got {type(self.metadata).__name__}")

        # Alleles
        if not self.v_alleles or len(self.v_alleles) == 0:
            errors.append("v_alleles is required and must be non-empty")
        if not self.j_alleles or len(self.j_alleles) == 0:
            errors.append("j_alleles is required and must be non-empty")

        # D alleles: required when metadata says has_d
        if self.metadata and self.metadata.has_d:
            if not self.d_alleles or len(self.d_alleles) == 0:
                errors.append("d_alleles is required for chains with D segments (metadata.has_d=True)")
        # Consistency: has_d=False but d_alleles provided
        if self.metadata and not self.metadata.has_d:
            if self.d_alleles and len(self.d_alleles) > 0:
                errors.append(
                    "d_alleles is non-empty but metadata.has_d=False. "
                    "Either remove d_alleles or set has_d=True."
                )

        # Gene usage
        if not self.gene_use_dict:
            errors.append("gene_use_dict is required and must be non-empty")
        else:
            for key in ("V", "J"):
                if key not in self.gene_use_dict:
                    errors.append(f"gene_use_dict must contain '{key}' key")

        # Trim dicts
        if not self.trim_dicts:
            errors.append("trim_dicts is required and must be non-empty")
        else:
            for key in ("V_3", "J_5"):
                if key not in self.trim_dicts:
                    errors.append(f"trim_dicts must contain '{key}' key")
            if self.metadata and self.metadata.has_d:
                for key in ("D_5", "D_3"):
                    if key not in self.trim_dicts:
                        errors.append(f"trim_dicts must contain '{key}' key for chains with D segments")

        if errors:
            raise DataConfigError(
                f"DataConfig '{self.name or 'Unnamed'}' failed validation:\n  - " + "\n  - ".join(errors)
            )

    def _unfold_alleles(self, gene_segment: str) -> List[str]:
        """Unfolds the alleles for a given gene segment (v, d, j, c) into a flat list."""
        alleles_dict = getattr(self, f"{gene_segment}_alleles")
        if alleles_dict is None:
            return []
        # This comprehension is clearer: iterate through the lists of alleles and flatten them.
        return [allele for allele_list in alleles_dict.values() for allele in allele_list]

    def _count_alleles(self, gene_segment: str) -> int:
        """Counts the number of alleles for a given gene segment."""
        return len(self._unfold_alleles(gene_segment))

    # ──────────────────────────────────────────────────────────────
    # Cartridge plane views
    # ──────────────────────────────────────────────────────────────
    #
    # The reference cartridge model groups DataConfig fields into
    # four planes — identity, catalogue, rules, empirical models.
    # These properties expose each plane as a frozen view dataclass
    # without moving the underlying fields. Same pickle, same
    # checksum; the views are a documentation/discovery surface.
    # See ``docs/reference_cartridge.md``.

    @property
    def cartridge_identity(self) -> "CartridgeIdentityView":
        """Read-only view of the **identity** plane (name +
        metadata). See :class:`CartridgeIdentityView`."""
        return CartridgeIdentityView(name=self.name, metadata=self.metadata)

    @property
    def cartridge_catalogue(self) -> "CartridgeCatalogueView":
        """Read-only view of the **catalogue** plane (V/D/J/C
        allele dicts). See :class:`CartridgeCatalogueView`."""
        return CartridgeCatalogueView(
            v_alleles=self.v_alleles,
            d_alleles=self.d_alleles,
            j_alleles=self.j_alleles,
            c_alleles=self.c_alleles,
        )

    @property
    def cartridge_rules(self) -> "CartridgeRulesView":
        """Read-only view of the **rules** plane
        (``reference_rules``). See :class:`CartridgeRulesView`.

        ``getattr(self, "reference_rules", None)`` is used so legacy
        pickles missing the field continue to surface ``None``
        rather than raising ``AttributeError``.
        """
        return CartridgeRulesView(
            reference_rules=getattr(self, "reference_rules", None),
        )

    @property
    def cartridge_models(self) -> "CartridgeModelsView":
        """Read-only view of the **empirical models** plane
        (``reference_models`` plus legacy ``NP_lengths`` /
        ``trim_dicts``). See :class:`CartridgeModelsView`."""
        return CartridgeModelsView(
            reference_models=getattr(self, "reference_models", None),
            legacy_np_lengths=self.NP_lengths,
            legacy_trim_dicts=self.trim_dicts,
        )

    # ──────────────────────────────────────────────────────────────
    # Cartridge manifest — single inspectable export surface
    # ──────────────────────────────────────────────────────────────

    def functional_status_counts(self) -> Dict[str, Dict[str, int]]:
        """Per-segment histogram of allele ``functional_status``
        values.

        Returns ``{segment: {status: count}}`` where ``segment`` is
        ``"v"`` / ``"d"`` / ``"j"`` / ``"c"`` and ``status`` is one
        of ``"functional"``, ``"orf"``, ``"pseudogene"``,
        ``"unknown"`` (the four canonical IMGT classifications) plus
        ``"unannotated"`` for alleles whose status is ``None``.

        Status normalisation matches the Python→Rust bridge: case-
        insensitive matching against
        ``{"functional", "orf", "pseudogene", "unknown"}`` plus the
        ``"F"`` / ``"P"`` aliases. Unknown / unrecognised strings
        collapse to ``"unannotated"`` (mirrors the bridge dropping
        them to ``None``).

        Empty / missing segment dicts contribute an all-zeros entry
        so the output shape is stable across cartridges.
        """
        from GenAIRR._refdata_resolver import _normalise_functional_status

        buckets = ("functional", "orf", "pseudogene", "unknown", "unannotated")
        out: Dict[str, Dict[str, int]] = {}
        segment_dicts = (
            ("v", self.v_alleles),
            ("d", self.d_alleles),
            ("j", self.j_alleles),
            ("c", self.c_alleles),
        )
        for seg, alleles_by_gene in segment_dicts:
            counts = {b: 0 for b in buckets}
            if alleles_by_gene:
                for _gene, allele_list in alleles_by_gene.items():
                    for allele in allele_list:
                        raw = getattr(allele, "functional_status", None)
                        normalised = _normalise_functional_status(raw)
                        if normalised is None:
                            counts["unannotated"] += 1
                        else:
                            counts[normalised] = counts.get(normalised, 0) + 1
            out[seg] = counts
        return out

    def cartridge_manifest(
        self, refdata: Optional[Any] = None
    ) -> Dict[str, Any]:
        """Return a stable, JSON-serialisable summary of this
        cartridge — the **single inspection surface** for the four
        planes plus curation, identity, hashes, and the documented
        completeness gaps.

        Output shape (audit §11 Slice 1):

        ::

            {
              "schema_version": int,
              "identity": {name, species, locus, reference_set, source},
              "catalogue": {v_count, d_count, j_count, c_count,
                            functional_status_counts: {v: {...}, d: {...}, ...}},
              "rules": {has_explicit_rules, allowed_bases, v_anchor, j_anchor},
              "models": {has_reference_models, np_length_keys, trim_keys,
                         legacy_np_lengths_present, legacy_trim_dicts_present},
              "curation": {source_tag, policies: [str, ...]},
              "hashes": {data_config_checksum, refdata_content_hash},
              "dropped_allele_fields": [str, ...],
              "orphan_dataconfig_fields": [str, ...],
              "errors": [str, ...],  # only when refdata bridge fails
            }

        Read-only — does NOT mutate this ``DataConfig``. Calling
        ``cartridge_manifest`` twice yields equal dicts and
        ``compute_checksum`` stays stable across calls.

        **Bridge cost.** Building the ``rules`` and ``hashes`` planes
        requires running ``dataconfig_to_refdata(self)`` (unless an
        ``refdata`` is passed in) — same cost as a real engine
        compile setup. For tooling that already has a refdata in
        hand (e.g. after ``refdata.curated(...)``), pass it via
        ``refdata`` to skip the rebuild and to surface curation
        tagging that lives on the refdata identity.

        **Error handling.** Bridge failures (invalid cartridge,
        missing fields, etc.) DO NOT raise. The failing planes
        fall back to safe defaults and the error message is
        appended to the ``"errors"`` list. Callers can inspect
        ``manifest["errors"]`` to decide whether the cartridge is
        actually usable.
        """
        return build_manifest(self, refdata, schema_version=SCHEMA_VERSION)

    # --- Public Properties ---
    @property
    def number_of_v_alleles(self) -> int:
        return self._count_alleles('v')

    @property
    def number_of_d_alleles(self) -> int:
        return self._count_alleles('d')

    @property
    def number_of_j_alleles(self) -> int:
        return self._count_alleles('j')

    @property
    def number_of_c_alleles(self) -> int:
        return self._count_alleles('c')

    # --- Public Methods ---
    def allele_list(self, gene_segment: str) -> List[str]:
        """Returns a flattened list of all alleles for a given gene segment."""
        return self._unfold_alleles(gene_segment)

    def copy(self):
        """
        Creates a deep, independent copy of this DataConfig object.

        Returns:
            DataConfig: A new DataConfig object with all attributes and nested
                        data structures duplicated.
        """
        return copy.deepcopy(self)

    def __repr__(self) -> str:
        parts = [f"<{self.name or 'Unnamed'} - Data Config>"]
        # Check each allele type before adding it to the representation
        if self.v_alleles is not None:
            parts.append(f"<{self.number_of_v_alleles} V Alleles>")
        if self.d_alleles is not None:
            parts.append(f"<{self.number_of_d_alleles} D Alleles>")
        if self.j_alleles is not None:
            parts.append(f"<{self.number_of_j_alleles} J Alleles>")
        if self.c_alleles is not None:
            parts.append(f"<{self.number_of_c_alleles} C Alleles>")
        return "-".join(parts)