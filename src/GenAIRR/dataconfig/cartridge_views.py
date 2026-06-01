"""Read-only views over the four cartridge planes of a
:class:`GenAIRR.dataconfig.DataConfig`.

A :class:`DataConfig` is, structurally, a flat bag of historical
fields (``v_alleles``, ``NP_lengths``, ``trim_dicts``,
``reference_rules``, ``reference_models``, ``metadata``, ...). The
cartridge architecture organises these into four planes:

1. **identity** — what cartridge is this? (name, metadata)
2. **catalogue** — which alleles exist? (V/D/J/C pools)
3. **rules** — how should the engine interpret alleles?
   (``reference_rules``)
4. **empirical models** — what are the default sampling
   distributions? (``reference_models`` plus legacy
   ``NP_lengths`` / ``trim_dicts``)

These view dataclasses give docs, tests, and tooling a place to talk
about the planes without moving any stored fields. The underlying
``DataConfig`` is the source of truth; the views are cheap snapshots
of what the cartridge currently exposes. No behaviour change — same
pickle, same checksum, same load path.

See ``docs/reference_cartridge.md`` for the full cartridge model
and the relationship between these planes, curation, and compile.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, List, Optional

if TYPE_CHECKING:  # pragma: no cover — typing-only
    from GenAIRR.alleles.allele import Allele
    from GenAIRR.dataconfig.config_info import ConfigInfo
    from GenAIRR.reference_models import ReferenceEmpiricalModels
    from GenAIRR.reference_rules import ReferenceRulesSpec


@dataclass(frozen=True)
class CartridgeIdentityView:
    """The **identity** plane — *what cartridge is this?*

    Mirrors the Rust ``ReferenceIdentity`` fields conceptually, but
    sourced directly from the Python-side ``DataConfig.name`` and
    ``DataConfig.metadata`` (a :class:`ConfigInfo`). Read-only;
    construct via :meth:`DataConfig.cartridge_identity`.
    """

    name: Optional[str]
    metadata: Optional["ConfigInfo"]


@dataclass(frozen=True)
class CartridgeCatalogueView:
    """The **catalogue** plane — *which alleles exist?*

    Each field is the underlying ``DataConfig`` allele dict
    (``{gene: [Allele, ...]}``) or ``None`` when absent. Reads
    are cheap (references, no copy). Read-only.
    """

    v_alleles: Optional[Dict[str, List["Allele"]]]
    d_alleles: Optional[Dict[str, List["Allele"]]]
    j_alleles: Optional[Dict[str, List["Allele"]]]
    c_alleles: Optional[Dict[str, List["Allele"]]]

    def _count(self, alleles: Optional[Dict[str, List["Allele"]]]) -> int:
        if alleles is None:
            return 0
        return sum(len(lst) for lst in alleles.values())

    @property
    def number_of_v_alleles(self) -> int:
        return self._count(self.v_alleles)

    @property
    def number_of_d_alleles(self) -> int:
        return self._count(self.d_alleles)

    @property
    def number_of_j_alleles(self) -> int:
        return self._count(self.j_alleles)

    @property
    def number_of_c_alleles(self) -> int:
        return self._count(self.c_alleles)


@dataclass(frozen=True)
class CartridgeRulesView:
    """The **rules** plane — *how should the engine interpret
    alleles?*

    Carries the cartridge's :class:`ReferenceRulesSpec`, or
    ``None`` when the cartridge relies on the loader's
    bundled-locus defaults. Read-only.
    """

    reference_rules: Optional["ReferenceRulesSpec"]


@dataclass(frozen=True)
class CartridgeModelsView:
    """The **empirical models** plane — *what are the default
    sampling distributions?*

    Exposes both surfaces:

    - ``reference_models`` is the typed, validated
      :class:`ReferenceEmpiricalModels` spec (or ``None``).
    - ``legacy_np_lengths`` / ``legacy_trim_dicts`` are the
      historical nested-dict shapes (``DataConfig.NP_lengths`` /
      ``DataConfig.trim_dicts``). The extractor in
      ``_dataconfig_extract`` consults the typed spec first and
      falls back to the legacy dicts; this view shows both so users
      can see what their cartridge actually carries.

    Read-only; both surfaces hold direct references to the
    underlying ``DataConfig`` storage (no copy).
    """

    reference_models: Optional["ReferenceEmpiricalModels"]
    legacy_np_lengths: Dict[str, Any]
    legacy_trim_dicts: Dict[str, Any]
