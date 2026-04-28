import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional

from GenAIRR.dataconfig.make.base_dataconfig_builder import BaseDataConfigGenerator
from GenAIRR.dataconfig.make.auxiliary_builders.trimming_proportion_builder import TrimmingProbabilityGenerator
from GenAIRR.dataconfig.make.auxiliary_builders.np_markov_chain_builder import NPMarkovParameterBuilder
from GenAIRR.utilities.data_utilities import create_allele_dict
from GenAIRR.utilities.asc_utilities import create_asc_germline_set

logger = logging.getLogger(__name__)


# ── BuildReport (T2-8 Phase 3b — pure-Python aggregator) ────────────


@dataclass
class BuildReportRejection:
    """One rejected allele's diagnostic record."""
    allele_name: str
    reason: str
    position: int = -1
    codon: str = "???"


@dataclass
class BuildReport:
    """Aggregate outcome of a `make_from_reference` invocation.

    Counts and per-allele rejection records collected as we route
    each LoadedAlleleRecord through the AnchorResolver. Use this to
    diagnose "why didn't allele X end up in my config?" without
    having to re-derive the answer from logs.

    Three rejection categories tracked separately:
      - Anchor rejections (`rejections`): the AnchorResolver couldn't
        find a valid anchor (typically pseudogenes or non-IMGT
        sequences without proper alignment).
      - Filter exclusions (`filtered_out`): the record was excluded
        by the F/ORF/P/partial filter policy (T2-9).
      - Notes: free-form annotations the user or builder added.

    Attached to ``DataConfig.build_report`` by ``make_from_reference``
    so diagnostics survive pickle round-trip.
    """
    total_seen: int = 0
    accepted: int = 0
    by_segment: dict = field(default_factory=lambda: defaultdict(int))
    by_status: dict = field(default_factory=lambda: defaultdict(int))
    by_confidence: dict = field(default_factory=lambda: defaultdict(int))
    rejections: list = field(default_factory=list)
    # T2-9: separate bucket for filter-policy exclusions. Kept apart
    # from anchor rejections because the two have different semantics
    # (filter = "user said don't include this category"; rejection =
    # "we tried and failed to resolve an anchor").
    filtered_out: list = field(default_factory=list)
    by_filter_status: dict = field(default_factory=lambda: defaultdict(int))
    notes: list = field(default_factory=list)

    def record_accepted(self, rec, anchor_result) -> None:
        self.total_seen += 1
        self.accepted += 1
        self.by_segment[rec.segment.name] += 1
        self.by_status[rec.functional_status.name] += 1
        self.by_confidence[anchor_result.confidence.name] += 1

    def record_rejected(self, rec, anchor_result) -> None:
        self.total_seen += 1
        self.by_segment[rec.segment.name] += 1
        self.by_status[rec.functional_status.name] += 1
        self.by_confidence[anchor_result.confidence.name] += 1
        self.rejections.append(BuildReportRejection(
            allele_name=rec.name,
            reason=anchor_result.reason or "(no reason)",
            position=anchor_result.position,
            codon=anchor_result.codon,
        ))

    def record_filtered(self, rec) -> None:
        """T2-9: log a record that was excluded by the filter policy
        (not an anchor failure — the policy chose to skip this
        functional-status category)."""
        self.total_seen += 1
        self.by_filter_status[rec.functional_status.name] += 1
        # Keep the list bounded — we don't want a 50K-allele filter
        # exclusion to balloon a pickle. Names only, no per-record
        # detail beyond name + status.
        self.filtered_out.append((rec.name, rec.functional_status.name))

    def summary(self) -> str:
        """One-line human summary, useful for logs."""
        parts = [
            f"BuildReport: accepted={self.accepted}/{self.total_seen}",
            f"rejected={len(self.rejections)}",
        ]
        if self.filtered_out:
            parts.append(f"filtered={len(self.filtered_out)} "
                         f"({dict(self.by_filter_status)})")
        parts.append(f"by_confidence={dict(self.by_confidence)}")
        return ", ".join(parts)


class RandomDataConfigBuilder(BaseDataConfigGenerator):
    """
    Generator that creates a DataConfig object with random (uniform/decaying) distributions
    for gene usage, trimming, and NP region parameters.
    """

    def __init__(self, convert_to_asc=False, *, species=None, chain_type=None, reference_set=None):
        super().__init__(
            convert_to_asc,
            species=species, chain_type=chain_type, reference_set=reference_set,
        )

    def _load_random_gene_usage(self):
        """
        Assign uniform probabilities to gene usage across all known alleles.
        """
        pointer_to_reference = self._get_reference_pointers()

        gene_use_dict = {allele: {} for allele in self.alleles}
        if self.has_d:
            gene_use_dict['D'] = {}

        for allele_type in gene_use_dict:
            current_alleles = list(pointer_to_reference[allele_type])
            n = len(current_alleles)
            gene_use_dict[allele_type] = {name: 1 / n for name in current_alleles}

        self.dataconfig.gene_use_dict = gene_use_dict

    def make(self, v_reference_path, j_reference_path, c_reference_path=None,
             d_reference_path=None, *, v_anchor_finder=None, j_anchor_finder=None,
             keep_anchorless=False):
        """
        Construct and return a DataConfig with all randomly simulated components.

        Args:
            v_reference_path: V gene FASTA file path or preloaded dict.
            j_reference_path: J gene FASTA file path or preloaded dict.
            c_reference_path: C gene FASTA file path or preloaded dict (optional).
            d_reference_path: D gene FASTA file path or preloaded dict (optional).
            v_anchor_finder: Custom V anchor callable (see ``create_allele_dict``).
            j_anchor_finder: Custom J anchor callable (see ``create_allele_dict``).
            keep_anchorless: If True, keep V/J alleles without anchors
                (they can be sampled but never produce productive sequences).

        Returns:
            DataConfig with random distributions and validated metadata.
        """
        # Read and load reference sequences
        v, j, c, d = self._read_reference_files(
            v_reference_path, j_reference_path, c_reference_path, d_reference_path,
            v_anchor_finder=v_anchor_finder, j_anchor_finder=j_anchor_finder,
            keep_anchorless=keep_anchorless,
        )
        self._load_alleles(v_alleles=v, j_alleles=j, c_alleles=c, d_alleles=d)
        logger.info('Alleles mounted to DataConfig')

        # Random gene usage and trimming proportions
        self._load_random_gene_usage()
        logger.info('Random gene usage mounted to DataConfig')

        self._load_trimming_probs()
        logger.info('Random trimming proportions mounted to DataConfig')

        np_generator = NPMarkovParameterBuilder(self.dataconfig, has_d=self.has_d)
        np_generator.generate_all_random(max_length=50)
        logger.info('Random NP parameters mounted to DataConfig')

        # Correction maps are computed lazily during graph compilation
        self._finalize()
        logger.info('DataConfig build complete')

        return self.dataconfig

    # ── Phase 2/3 public API: make_from_reference ─────────────────

    def make_from_reference(self,
                            v_reference: str,
                            j_reference: str,
                            d_reference: Optional[str] = None,
                            c_reference: Optional[str] = None,
                            *,
                            v_format: Optional[str] = None,
                            j_format: Optional[str] = None,
                            v_locus_hint=None,
                            j_locus_hint=None,
                            v_sidecar: Optional[str] = None,
                            j_sidecar: Optional[str] = None,
                            keep_anchorless: bool = False,
                            strict: bool = False,
                            include_orf: bool = True,
                            include_pseudo: bool = False,
                            include_partial: bool = False):
        """Build a DataConfig from any of the supported reference formats.

        T2-8 Phase 2/3 entry point. Each segment can be a different
        format — common patterns:

        - IMGT V-QUEST: pass FASTA paths, no other args needed.
        - AIRR-C: pass ``.json`` paths.
        - OGRDB: pass FASTA paths plus ``v_sidecar=`` / ``j_sidecar=``
          (the corresponding ``.json`` files).
        - IgBLAST: pass FASTAs plus ``.aux`` / ``.ndm.imgt`` sidecars.
        - Plain FASTA: pass FASTA paths and locus hints.

        ``v_format`` / ``j_format`` are optional explicit overrides
        (``"imgt"``, ``"airrc"``, ``"ogrdb"``, ``"igblast"``,
        ``"plain"``); if omitted, format auto-detection runs.

        Returns:
            (DataConfig, BuildReport) — the report aggregates per-
            allele resolution outcomes for diagnostics.

        Raises:
            ValueError: if both V and J resolution end up empty
                (almost certainly a wrong path or wrong format).
        """
        from .._native_helpers import (   # local import to avoid cycle
            load_segment_alleles,
        )

        report = BuildReport()
        # T2-9: filter flags apply to every segment so users get
        # consistent F/ORF/P/partial policy across V/D/J/C.
        filter_kwargs = dict(
            include_orf=include_orf,
            include_pseudo=include_pseudo,
            include_partial=include_partial,
        )
        v_alleles = load_segment_alleles(
            v_reference, segment="V", format=v_format,
            locus_hint=v_locus_hint, sidecar_path=v_sidecar,
            keep_anchorless=keep_anchorless, strict=strict,
            report=report, **filter_kwargs)
        j_alleles = load_segment_alleles(
            j_reference, segment="J", format=j_format,
            locus_hint=j_locus_hint, sidecar_path=j_sidecar,
            keep_anchorless=keep_anchorless, strict=strict,
            report=report, **filter_kwargs)
        d_alleles = (load_segment_alleles(
            d_reference, segment="D", report=report, **filter_kwargs)
            if d_reference else None)
        c_alleles = (load_segment_alleles(
            c_reference, segment="C", report=report, **filter_kwargs)
            if c_reference else None)

        if not v_alleles or not j_alleles:
            raise ValueError(
                f"make_from_reference: empty V or J after loading "
                f"(report: {report.summary()}). Check the file paths "
                f"and format hints, or pass keep_anchorless=True.")

        # has_d toggles allele list construction + downstream validation.
        # Setting it from the input keeps the builder honest about
        # which chain class we're producing.
        self.has_d = d_alleles is not None and len(d_alleles) > 0

        self._load_alleles(v_alleles=v_alleles, j_alleles=j_alleles,
                           c_alleles=c_alleles, d_alleles=d_alleles)
        logger.info("Alleles mounted to DataConfig (%s)", report.summary())

        self._load_random_gene_usage()
        self._load_trimming_probs()

        np_generator = NPMarkovParameterBuilder(self.dataconfig,
                                                has_d=self.has_d)
        np_generator.generate_all_random(max_length=50)

        # T2-8 closure: attach the build report to the DataConfig
        # BEFORE _finalize() so the report's diagnostic data is
        # carried in any downstream pickle. Excluded from the
        # checksum (see compute_checksum), so it won't invalidate
        # integrity if mutated post-load.
        self.dataconfig.build_report = report

        self._finalize()
        # T2-8 closure: stamp the integrity checksum after validation
        # so a custom-built DataConfig round-trips cleanly through
        # the same verify_integrity pathway as shipped builtins.
        # `compute_checksum` is hash-stable across legacy/new pickle
        # shapes (see DataConfig.compute_checksum docstring).
        self.dataconfig.schema_sha256 = self.dataconfig.compute_checksum()
        logger.info("DataConfig build complete; checksum stamped: %s...",
                    self.dataconfig.schema_sha256[:16])
        return self.dataconfig, report
