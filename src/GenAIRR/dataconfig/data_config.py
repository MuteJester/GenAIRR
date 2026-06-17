from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any

from GenAIRR.dataconfig.cartridge_views import (
    CartridgeCatalogueView,
    CartridgeIdentityView,
    CartridgeModelsView,
    CartridgeRulesView,
)
from GenAIRR.dataconfig.config_info import ConfigInfo
from GenAIRR.dataconfig.enums import ChainType
import copy
import hashlib
import pickle
from GenAIRR.alleles.allele import Allele
from GenAIRR.reference_models import ReferenceEmpiricalModels
from GenAIRR.reference_rules import ReferenceRulesSpec
from GenAIRR.genotype_priors import PopulationGenotypeModel


DEFAULT_P_NUCLEOTIDE_LENGTH_PROBS = {0: 0.50, 1: 0.25, 2: 0.15, 3: 0.07, 4: 0.03}

# Bump on any breaking change to DataConfig fields. Old pickles without
# this version (or with an older one) refuse to load and prompt the user
# to re-migrate. See _check_schema() / verify_integrity() below.
SCHEMA_VERSION = 1


# Documented completeness gaps surfaced verbatim by
# ``DataConfig.cartridge_manifest``. See
# ``docs/reference_cartridge_completeness_audit.md`` §1 (orphan
# fields) and §2 (dropped allele fields) for the rationale on each
# entry. These are tuples so manifest output stays stable across
# calls and across Python sessions.
_DOCUMENTED_DROPPED_ALLELE_FIELDS = (
    "aliases",
    "anchor_meta",
    "gapped_seq",
    "family",
    "locus",
    "species",
    "source",
)
_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS = (
    "gene_use_dict",
    "NP_transitions",
    "NP_first_bases",
    "correction_maps",
    "asc_tables",
    "p_nucleotide_length_probs",
    "dj_pairing_map",
)


def _jsonable(value):
    """Coerce a value to a JSON-clean form. Enums become their
    ``.name`` (uppercase canonical) for stability; tuples become
    lists; everything else passes through. Used by
    :meth:`DataConfig.cartridge_manifest` to make the identity
    plane safe for ``json.dumps``."""
    if value is None:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if hasattr(value, "name") and not isinstance(value, type):
        # Enum-shaped values (Species, ChainType, etc.).
        return value.name
    if isinstance(value, (tuple, list)):
        return [_jsonable(v) for v in value]
    return str(value)


def _normalise_anchor_rule(rule):
    """Coerce an anchor rule dict from the Rust ``rules()`` accessor
    into a JSON-clean shape. Returns ``None`` for missing rules.

    The Rust side may surface the rule as a dict with ``required``,
    ``expected_aa``, ``missing_severity``, ``mismatch_severity``
    keys (string severities). We pass it through but coerce any
    iterable fields to plain ``list`` to keep ``json.dumps`` happy.
    """
    if rule is None:
        return None
    if isinstance(rule, dict):
        out = dict(rule)
        if "expected_aa" in out and out["expected_aa"] is not None:
            out["expected_aa"] = list(out["expected_aa"])
        return out
    # Unknown shape — drop silently rather than raise; the manifest
    # error list isn't the right channel for "this is unexpected
    # but harmless".
    return None


def _allele_usage_manifest_block(cfg):
    """Build the ``models.allele_usage`` manifest block per the
    Allele Usage Estimation v1 audit (§8 of
    ``docs/allele_usage_estimation_design.md``).

    Returns a JSON-clean dict carrying:

    - ``available``: whether ``reference_models.allele_usage`` is
      a non-``None`` spec.
    - ``segments``: the canonical V/D/J label list.
    - ``nonempty_segments``: subset of ``segments`` whose
      author-supplied dict is non-empty.
    - ``legacy_gene_use_dict_present``: documents the orphan
      legacy dict's bundled-cartridge presence so a downstream
      tool can decide whether to author a typed plane from it.
    - ``legacy_fallback``: always ``False`` in v1 (no auto-lift).
    - ``in_plan_signature``: ``False`` — inherits soft gap 1
      from the plan-signature completeness audit. Tightening is
      a separate slice covering BOTH the kwarg AND the
      cartridge-driven path in lockstep.
    - ``source``: which Python field the manifest block reads
      from. Always ``"ReferenceEmpiricalModels.allele_usage"``.
    """
    rm = getattr(cfg, "reference_models", None)
    spec = getattr(rm, "allele_usage", None) if rm is not None else None
    available = spec is not None
    nonempty: list = []
    if available:
        try:
            nonempty = list(spec.nonempty_segments())
        except Exception:
            nonempty = []
    return {
        "available": available,
        "segments": ["V", "D", "J"],
        "nonempty_segments": nonempty,
        "legacy_gene_use_dict_present": bool(
            getattr(cfg, "gene_use_dict", None)
        ),
        "legacy_fallback": False,
        "in_plan_signature": False,  # documented soft gap 1
        "source": "ReferenceEmpiricalModels.allele_usage",
    }


def _genotype_priors_manifest_block(cfg):
    """Build the ``models.genotype_priors`` manifest block (Slice — Cartridge
    genotype plane). Audit-sized: counts and identity, never the full tables.
    Reads the top-level ``DataConfig.genotype_priors`` plane (independent of
    ``reference_models``)."""
    def _empty(available, valid, error=None):
        b = {
            "available": available,
            "valid": valid,
            "model_id": None,
            "source": None,
            "version": None,
            "model_checksum": None,
            "segments_with_frequencies": [],
            "freq_gene_counts": {"V": 0, "D": 0, "J": 0},
            "deletion_gene_counts": {"V": 0, "D": 0, "J": 0},
            "novel_allele_count": 0,
            "chromosome_weights": None,
            "source_field": "DataConfig.genotype_priors",
        }
        if error is not None:
            b["validation_error"] = error
        return b

    model = getattr(cfg, "genotype_priors", None)
    if model is None:
        return _empty(False, None)
    if not isinstance(model, PopulationGenotypeModel):
        # Field is typed Optional[PopulationGenotypeModel] but Python won't enforce
        # it; a garbage value must not crash the manifest (called from build()).
        return _empty(True, False,
                      f"genotype_priors is not a PopulationGenotypeModel "
                      f"(got {type(model).__name__})")

    # A plane may have been attached directly (bypassing builder validation).
    # Report validity rather than leaking non-JSON-clean numerics (e.g. NaN
    # chromosome weights) into the manifest.
    try:
        model.validate(
            chain_type=getattr(getattr(cfg, "metadata", None), "chain_type", None))
        valid, validation_error = True, None
    except ValueError as exc:
        valid, validation_error = False, str(exc)

    # content_checksum / float() can themselves raise on a malformed model; never
    # let that crash the manifest — report None and the validity flag instead.
    try:
        checksum = model.content_checksum()
    except Exception:
        checksum = None
    cw = None
    if valid:
        try:
            cw = [float(model.chromosome_weights[0]), float(model.chromosome_weights[1])]
        except Exception:
            cw = None
    freq = model.allele_frequencies if isinstance(model.allele_frequencies, dict) else {}
    dele = model.haplotype_deletion_prob if isinstance(model.haplotype_deletion_prob, dict) else {}
    novels = model.novel_alleles if isinstance(model.novel_alleles, (list, tuple)) else []
    block = {
        "available": True,
        "valid": valid,
        "model_id": (model.model_id or None) if isinstance(model.model_id, str) else None,
        "source": (model.source or None) if isinstance(model.source, str) else None,
        "version": (model.version or None) if isinstance(model.version, str) else None,
        "model_checksum": checksum,
        "segments_with_frequencies": [s for s in ("V", "D", "J") if freq.get(s)],
        "freq_gene_counts": {s: len(freq.get(s, {})) for s in ("V", "D", "J")},
        "deletion_gene_counts": {s: len(dele.get(s, {})) for s in ("V", "D", "J")},
        "novel_allele_count": len(novels),
        "chromosome_weights": cw,
        "source_field": "DataConfig.genotype_priors",
    }
    if validation_error is not None:
        block["validation_error"] = validation_error
    return block


def _np_length_models_manifest_block(cfg):
    """Build the ``models.np_length_models`` manifest block
    per the NP Length Distribution Estimation v1 audit
    (§8 of ``docs/np_length_estimation_design.md``).

    Returns a JSON-clean dict carrying:

    - ``keys``: sorted list of NP-key strings present on
      the typed ``ReferenceEmpiricalModels.np_lengths``
      plane (subset of ``("NP1", "NP2")``).
    - ``source``: which Python field the manifest block
      reads from. Always ``"ReferenceEmpiricalModels.np_lengths"``.
    - ``in_plan_signature``: ``True`` — NP-length
      distributions fold into the plan signature via
      ``GenerateNPPass.parameter_signature`` /
      ``fmt_int_dist`` (verified at audit time). **No
      soft gap inherited.**
    - ``legacy_np_lengths_present``: documents the bundled
      cartridges' legacy nested-dict presence so a
      downstream tool can decide whether to author a
      typed plane from it.
    - ``legacy_fallback``: always ``False`` in v1 (no
      auto-lift).
    """
    rm = getattr(cfg, "reference_models", None)
    np_lengths = getattr(rm, "np_lengths", None) if rm is not None else None
    keys = sorted(np_lengths.keys()) if np_lengths else []
    return {
        "keys": keys,
        "source": "ReferenceEmpiricalModels.np_lengths",
        "in_plan_signature": True,
        "legacy_np_lengths_present": bool(getattr(cfg, "NP_lengths", None)),
        "legacy_fallback": False,
    }


def _trim_models_manifest_block(cfg):
    """Build the ``models.trim_models`` manifest block per the
    Trim Distribution Estimation v1 audit (§8 of
    ``docs/trim_distribution_estimation_design.md``).

    Returns a JSON-clean dict carrying:

    - ``keys``: sorted list of trim-key strings present on the
      typed ``ReferenceEmpiricalModels.trims`` plane (subset of
      ``("V_3", "D_5", "D_3", "J_5")``).
    - ``source``: which Python field the manifest block reads
      from. Always ``"ReferenceEmpiricalModels.trims"``.
    - ``in_plan_signature``: ``True`` — trim distributions fold
      into the plan signature via ``fmt_int_dist`` (verified at
      audit time). **No soft gap inherited** — unlike the
      allele-usage slice.
    - ``legacy_trim_dicts_present``: documents the bundled
      cartridges' legacy nested-dict presence so a downstream
      tool can decide whether to author a typed plane from it.
    - ``legacy_fallback``: always ``False`` in v1 (no auto-lift).
    """
    rm = getattr(cfg, "reference_models", None)
    trims = getattr(rm, "trims", None) if rm is not None else None
    keys = sorted(trims.keys()) if trims else []
    return {
        "keys": keys,
        "source": "ReferenceEmpiricalModels.trims",
        "in_plan_signature": True,
        "legacy_trim_dicts_present": bool(getattr(cfg, "trim_dicts", None)),
        "legacy_fallback": False,
    }


def _p_nucleotide_length_keys(cfg):
    """Return the sorted list of P-end keys that have a typed
    `EmpiricalDistributionSpec` on
    `cfg.reference_models.p_nucleotide_lengths`. Returns an empty
    list when no typed P-plane is authored — the bundled
    cartridges' default.
    """
    rm = getattr(cfg, "reference_models", None)
    if rm is None:
        return []
    plane = getattr(rm, "p_nucleotide_lengths", None) or {}
    return sorted(plane.keys())


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
                f"{SCHEMA_VERSION}. Re-migrate this pickle (see "
                f".private/scripts/migrate_pkl_to_v1.py) or rebuild it "
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
        errors: List[str] = []

        # Identity — assembled from BOTH the Python-side ConfigInfo
        # (species, reference_set, chain-type-derived data) and the
        # bridged refdata.identity() (which exposes locus, source,
        # plus the curation source-tag suffix when applicable). The
        # bridged identity is the authoritative answer for fields
        # the engine consumes; the ConfigInfo fields fill gaps the
        # bridge doesn't surface.
        meta = self.metadata
        identity_dict: Dict[str, Any] = {
            "name": _jsonable(self.name),
            "species": _jsonable(getattr(meta, "species", None) if meta else None),
            "locus": None,  # populated from bridged identity below
            "reference_set": _jsonable(
                getattr(meta, "reference_set", None) if meta else None
            ),
            "source": None,  # populated from bridged identity below
        }

        # Catalogue counts + functional-status histogram. These read
        # only the Python-side allele lists; no bridge needed.
        catalogue_dict: Dict[str, Any] = {
            "v_count": self.number_of_v_alleles,
            "d_count": self.number_of_d_alleles,
            "j_count": self.number_of_j_alleles,
            "c_count": self.number_of_c_alleles,
            "functional_status_counts": self.functional_status_counts(),
        }

        # Models plane — pure Python; no bridge.
        ref_models = getattr(self, "reference_models", None)
        has_reference_models = ref_models is not None
        np_length_keys: List[str] = []
        trim_keys: List[str] = []
        # Typed NP base model inventory (Slice — Typed NP base
        # model). ``np_base_models`` is keyed by NP region
        # (``"NP1"`` / ``"NP2"``) and surfaces the per-region
        # ``kind``. Cartridges without typed NP base models report
        # ``kind="uniform"`` implicitly via the empty / absent
        # entry; the manifest summary block below makes the
        # opt-in / fallback explicit so a consumer can detect
        # which cartridges author non-uniform NP biology.
        np_base_models: Dict[str, Any] = {}
        if has_reference_models:
            np_length_keys = sorted((ref_models.np_lengths or {}).keys())
            trim_keys = sorted((ref_models.trims or {}).keys())
            np_bases_dict = getattr(ref_models, "np_bases", None) or {}
            for key in sorted(np_bases_dict.keys()):
                spec = np_bases_dict[key]
                np_base_models[key] = {"kind": spec.kind}
        # SHM model inventory — documents which mutation models the
        # engine ships and which S5F kernels are bundled. The kernel
        # choice itself is a per-experiment parameter (passed to
        # ``Experiment.mutate(...)``), NOT a per-cartridge property,
        # so the manifest reports what's *available* rather than
        # what was used. Explicitly carries ``in_content_hash=False``
        # so a user reading the manifest sees that two simulation
        # runs using different S5F kernels would be indistinguishable
        # by ``content_hash`` alone (the v1 boundary pinned by the
        # SHM model audit §3).
        from GenAIRR._s5f_loader import (
            DEFAULT_S5F_KERNEL,
            available_s5f_kernels,
            builtin_s5f_kernel_digest,
        )

        shm_dict: Dict[str, Any] = {
            "available_models": ["uniform", "s5f"],
            "s5f_kernels_available": available_s5f_kernels(),
            "default_s5f_kernel": DEFAULT_S5F_KERNEL,
            # Digest of the default kernel's bundled pickle bytes.
            # ``None`` if the bytes can't be read (corrupted
            # install) — the manifest doesn't raise.
            "s5f_kernel_digest": builtin_s5f_kernel_digest(DEFAULT_S5F_KERNEL),
            # Hard-coded ``False`` because the bridge doesn't fold
            # the kernel choice into Rust ``content_hash``; the v1
            # boundary pinned in the audit. A future slice that
            # closes the boundary updates BOTH this field and the
            # companion audit pin.
            "in_content_hash": False,
            # Per-segment SHM rate scalars — the segment-rates
            # slice. The manifest advertises the capability + the
            # flat default; rate vectors are per-experiment (passed
            # to ``Experiment.mutate(segment_rates=...)``) and
            # therefore not part of cartridge identity. The
            # ``in_content_hash`` field documents that the rate
            # vector itself doesn't enter ``content_hash`` — same
            # v1 boundary as the kernel choice.
            "segment_rate_support": {
                "available": True,
                "buckets": ["V", "D", "J", "NP"],
                "default": {"V": 1.0, "D": 1.0, "J": 1.0, "NP": 1.0},
                "in_content_hash": False,
            },
            # V-region substructure annotation surface — the Slice 1
            # output of the V-subregion cartridge work. ``available``
            # is True as soon as the bridge derives subregions (so it
            # toggles on for cartridges with IMGT ``gapped_seq``);
            # ``annotated_v_count`` / ``total_v_count`` quantify
            # coverage. The ``derivation`` tag identifies which
            # source produced the annotations
            # (``"bridge_imgt_gapped_seq"`` for the default helper,
            # ``"explicit"`` if the user set ``allele.subregions``
            # directly). ``in_content_hash`` is True — the
            # ``refdata_content_hash`` does fold subregion intervals
            # in (see ``engine_rs/src/trace_file.rs``). Sampling
            # rates (``v_subregion_rates``) and CDR/FR mutation
            # counters do NOT live here; those are deliberately
            # absent and pinned by the contract suite.
            "v_subregion_support": {
                "available": False,  # populated below from bridged refdata
                "labels": ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"],
                "annotated_v_count": 0,
                "total_v_count": self.number_of_v_alleles,
                "derivation": "bridge_imgt_gapped_seq",
                "in_content_hash": True,
            },
            # V-subregion SHM rate-vector capability — the
            # ``v_subregion_rates`` kwarg on ``Experiment.mutate``
            # (Slice B). The manifest advertises the kwarg vocabulary
            # (five canonical labels + two aliases ``FWR`` / ``CDR``)
            # plus the flat default. ``in_plan_signature=True``
            # because the rate vector flows into
            # ``pass_plan_signature`` via Slice A's mutate-pass
            # ``parameter_signature`` (mismatched rates fail replay
            # at the signature gate). ``in_content_hash=False``
            # because rates are a per-experiment parameter, not a
            # cartridge property — same boundary the
            # ``segment_rate_support`` block documents.
            "v_subregion_rate_support": {
                "available": True,
                "labels": ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"],
                "aliases": {
                    "FWR": ["FWR1", "FWR2", "FWR3"],
                    "CDR": ["CDR1", "CDR2"],
                },
                "default": {
                    "FWR1": 1.0,
                    "CDR1": 1.0,
                    "FWR2": 1.0,
                    "CDR2": 1.0,
                    "FWR3": 1.0,
                },
                "in_plan_signature": True,
                "in_content_hash": False,
            },
            # V-subregion realised SHM mutation counters (Slice —
            # ``docs/v_subregion_mutation_counters_audit.md``). Six
            # AIRR fields that partition ``n_v_mutations`` by the
            # assigned V allele's IMGT subregion intervals: five
            # canonical labels plus ``n_v_unannotated_mutations``
            # for V events that can't be attributed. Source of
            # truth is the engine's event ledger filtered to the
            # mutate.{uniform,s5f} passes — the validator surfaces
            # six per-field mismatch issues + a cross-field
            # ``VSubregionMutationCountSumMismatch`` invariant.
            # ``requires_annotations=True`` documents the
            # dependency on the Slice-1 cartridge annotation
            # surface — without subregion intervals, every V SHM
            # event routes to the unannotated bucket.
            # ``in_content_hash=False`` because counters are a
            # per-record observable, not a cartridge property.
            "v_subregion_counter_support": {
                "available": True,
                "fields": [
                    "n_fwr1_mutations",
                    "n_cdr1_mutations",
                    "n_fwr2_mutations",
                    "n_cdr2_mutations",
                    "n_fwr3_mutations",
                    "n_v_unannotated_mutations",
                ],
                "partition_of": "n_v_mutations",
                "requires_annotations": True,
                "unannotated_bucket": "n_v_unannotated_mutations",
                "in_content_hash": False,
            },
        }

        models_dict: Dict[str, Any] = {
            "has_reference_models": has_reference_models,
            "np_length_keys": np_length_keys,
            "trim_keys": trim_keys,
            "legacy_np_lengths_present": bool(self.NP_lengths),
            "legacy_trim_dicts_present": bool(self.trim_dicts),
            "shm": shm_dict,
            # V(D)J N-addition base sampling models (Slice — Typed
            # NP base model). The block advertises which NP regions
            # carry typed cartridge-owned base distributions; an
            # empty ``np_base_models`` dict means every NP position
            # samples uniformly A/C/G/T (the pre-slice default).
            # ``legacy_fallback`` documents whether the resolver
            # auto-lifted the legacy ``NP_transitions`` /
            # ``NP_first_bases`` orphan fields — v1 keeps this
            # ``False`` because auto-lift would silently change
            # output bytes vs the pre-slice baseline; a follow-up
            # slice may enable it under an explicit opt-in. The
            # ``legacy_np_*_present`` flags surface the bundled
            # cartridges' orphan data so a downstream tool can
            # decide whether to author a typed spec from it.
            # ``in_plan_signature=True`` because the base
            # distribution's ``support()`` enters
            # ``GenerateNPPass::parameter_signature`` via
            # ``fmt_byte_dist`` (Slice A discipline).
            # ``in_content_hash=False`` because the model is a
            # per-experiment defaults plane, not cartridge
            # identity — same boundary the per-segment SHM rates
            # block documents.
            "np_base_models": {
                "models": np_base_models,
                "legacy_fallback": False,
                "legacy_np_transitions_present": bool(
                    getattr(self, "NP_transitions", None)
                ),
                "legacy_np_first_bases_present": bool(
                    getattr(self, "NP_first_bases", None)
                ),
                "supported_kinds": [
                    "uniform",
                    "empirical_first_base",
                    "markov",
                ],
                "deferred_kinds": [],
                "in_plan_signature": True,
                "in_content_hash": False,
            },
            # P-nucleotide per-end length distribution inventory
            # (Slice — P-nucleotide v1). `length_keys` advertises
            # which V(D)J coding-end junction sides the cartridge
            # authored a typed `EmpiricalDistributionSpec` for.
            # Empty list (the bundled-cartridge default) means
            # the pipeline omits every P-pass — byte-identical
            # to the pre-slice baseline.
            # `legacy_p_nucleotide_length_probs_present` surfaces
            # the bundled cartridges' orphan dict so a downstream
            # tool can decide whether to author a typed plane
            # from it; `legacy_fallback=False` documents that
            # auto-lift is NOT enabled (same boundary the Markov
            # slice held for `NP_transitions`).
            # `in_plan_signature=True` because each `PAdditionPass`
            # folds its per-end length distribution via
            # `fmt_int_dist` into `parameter_signature`, so
            # replay against a different P-length distribution
            # fails the signature gate.
            "p_nucleotide_models": {
                "length_keys": _p_nucleotide_length_keys(self),
                "legacy_p_nucleotide_length_probs_present": bool(
                    getattr(self, "p_nucleotide_length_probs", None)
                ),
                "legacy_fallback": False,
                "supported_ends": ["V_3", "D_5", "D_3", "J_5"],
                "in_plan_signature": True,
                "in_content_hash": False,
            },
            # Per-segment allele-usage weights (Slice — Allele
            # Usage Estimation v1). ``available`` is ``True`` only
            # when the cartridge author / builder attached a
            # non-``None`` ``allele_usage`` spec on
            # ``reference_models``. ``nonempty_segments`` lists
            # the segments whose dicts are non-empty so a
            # downstream consumer can decide whether the
            # cartridge actually overrides V / D / J or only some
            # of them. ``legacy_gene_use_dict_present`` surfaces
            # the bundled cartridges' orphan dict so a downstream
            # tool can decide whether to author a typed plane
            # from it; ``legacy_fallback=False`` documents that
            # auto-lift is NOT enabled (same boundary the Markov /
            # P-nucleotide slices respected). ``in_plan_signature``
            # is the documented soft gap 1 from the plan-signature
            # completeness audit — the per-experiment kwarg
            # surface AND this cartridge-driven path both fold
            # through `SampleAllelePass.parameter_signature` which
            # returns the empty string; tightening is a separate
            # slice.
            "allele_usage": _allele_usage_manifest_block(self),
            # Trim Distribution Estimation v1 — typed-plane
            # `trims` block. ``in_plan_signature=True`` because
            # trim distributions already fold into the plan
            # signature via ``fmt_int_dist`` (unlike the
            # allele-usage soft gap). The minimal top-level
            # ``trim_keys`` / ``legacy_trim_dicts_present``
            # entries above remain for backwards compatibility;
            # this block is the structured surface.
            "trim_models": _trim_models_manifest_block(self),
            # NP Length Distribution Estimation v1 —
            # typed-plane `np_lengths` block. Same shape /
            # discipline as `trim_models`: ``in_plan_signature
            # =True`` (folded via ``GenerateNPPass.parameter_signature``),
            # ``legacy_fallback=False`` (no auto-lift of the
            # legacy ``DataConfig.NP_lengths`` orphan dict),
            # additive on top of the minimal top-level
            # ``np_length_keys`` / ``legacy_np_lengths_present``
            # entries above.
            "np_length_models": _np_length_models_manifest_block(self),
            # Donor-population germline prior plane (Slice — Cartridge genotype
            # plane). Read from the top-level ``DataConfig.genotype_priors``
            # field (NOT ``reference_models``); ``source_field`` records that.
            "genotype_priors": _genotype_priors_manifest_block(self),
        }

        # Bridge once (or accept the provided refdata) to read the
        # rules plane + content_hash + identity.source curation tag.
        bridged_refdata = refdata
        if bridged_refdata is None:
            try:
                from GenAIRR._refdata_resolver import dataconfig_to_refdata
                bridged_refdata = dataconfig_to_refdata(self)
            except Exception as exc:
                errors.append(
                    f"dataconfig_to_refdata failed: "
                    f"{type(exc).__name__}: {exc}"
                )
                bridged_refdata = None

        rules_dict: Dict[str, Any] = {
            "has_explicit_rules": getattr(self, "reference_rules", None)
            is not None,
            "allowed_bases": None,
            "v_anchor": None,
            "j_anchor": None,
        }
        refdata_content_hash: Optional[str] = None
        curation_source_tag: Optional[str] = None
        curation_policies: List[str] = []

        if bridged_refdata is not None:
            # Rules plane — three Rust accessors:
            # ``allowed_bases()`` + ``v_anchor_rule()`` +
            # ``j_anchor_rule()``. Each returns a JSON-friendly dict
            # / list; ``_normalise_anchor_rule`` coerces tuples to
            # lists for ``json.dumps`` safety.
            try:
                rules_dict["allowed_bases"] = list(
                    bridged_refdata.allowed_bases()
                )
            except Exception as exc:
                errors.append(
                    f"refdata.allowed_bases() failed: "
                    f"{type(exc).__name__}: {exc}"
                )
            try:
                rules_dict["v_anchor"] = _normalise_anchor_rule(
                    bridged_refdata.v_anchor_rule()
                )
            except Exception as exc:
                errors.append(
                    f"refdata.v_anchor_rule() failed: "
                    f"{type(exc).__name__}: {exc}"
                )
            try:
                rules_dict["j_anchor"] = _normalise_anchor_rule(
                    bridged_refdata.j_anchor_rule()
                )
            except Exception as exc:
                errors.append(
                    f"refdata.j_anchor_rule() failed: "
                    f"{type(exc).__name__}: {exc}"
                )
            try:
                refdata_content_hash = bridged_refdata.content_hash()
            except Exception as exc:
                errors.append(
                    f"refdata.content_hash() failed: "
                    f"{type(exc).__name__}: {exc}"
                )
            # V-subregion coverage — walk the bridged V pool and
            # count alleles with non-empty subregion lists. Failures
            # here go through the manifest's standard ``errors``
            # channel rather than raising — the rest of the manifest
            # is still useful even if coverage stats can't be read.
            try:
                annotated_v = 0
                v_total = bridged_refdata.v_pool_size()
                for v_id in range(v_total):
                    if bridged_refdata.v_allele(v_id).subregions:
                        annotated_v += 1
                shm_dict["v_subregion_support"]["available"] = (
                    annotated_v > 0
                )
                shm_dict["v_subregion_support"]["annotated_v_count"] = (
                    annotated_v
                )
                shm_dict["v_subregion_support"]["total_v_count"] = v_total
            except Exception as exc:
                errors.append(
                    f"v_subregion_support coverage failed: "
                    f"{type(exc).__name__}: {exc}"
                )
            # Identity (locus + source) from the bridged side, plus
            # curation parsing. The Rust curation path appends
            # ``|curated:<policy_tag>`` per applied policy; the
            # policy tag itself may contain ``|`` (e.g.
            # ``functional_status:functional|keep_unannotated=true``).
            # So we split on the ``|curated:`` separator after the
            # first base-source piece. The base source becomes
            # ``curation.source_tag``; the trailing policy tags
            # become ``curation.policies``.
            try:
                ident = bridged_refdata.identity() or {}
                bridged_locus = ident.get("locus")
                if bridged_locus is not None:
                    identity_dict["locus"] = _jsonable(bridged_locus)
                source_tag = ident.get("source")
                curation_source_tag = source_tag
                if source_tag is not None:
                    identity_dict["source"] = _jsonable(source_tag)
                if source_tag and "|curated:" in source_tag:
                    pieces = source_tag.split("|curated:")
                    curation_policies = pieces[1:]
            except Exception as exc:
                errors.append(
                    f"refdata.identity() failed: {type(exc).__name__}: {exc}"
                )

        curation_dict: Dict[str, Any] = {
            "source_tag": curation_source_tag,
            "policies": curation_policies,
        }

        hashes_dict: Dict[str, Any] = {
            "data_config_checksum": self.compute_checksum(),
            "refdata_content_hash": refdata_content_hash,
        }

        # Documented completeness gaps from the audit — listed
        # verbatim so a user inspecting a manifest sees which
        # surfaces aren't covered by the views / the bridge.
        dropped_allele_fields = list(_DOCUMENTED_DROPPED_ALLELE_FIELDS)
        orphan_dataconfig_fields = list(_DOCUMENTED_ORPHAN_DATACONFIG_FIELDS)

        manifest: Dict[str, Any] = {
            "schema_version": getattr(self, "schema_version", SCHEMA_VERSION),
            "identity": identity_dict,
            "catalogue": catalogue_dict,
            "rules": rules_dict,
            "models": models_dict,
            "curation": curation_dict,
            "hashes": hashes_dict,
            "dropped_allele_fields": dropped_allele_fields,
            "orphan_dataconfig_fields": orphan_dataconfig_fields,
        }
        if errors:
            manifest["errors"] = errors
        return manifest

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