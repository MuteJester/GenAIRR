"""
Helpers that bridge the C-side ReferenceLoader / AnchorResolver
pipeline to the Python-side DataConfig builder (T2-8 Phase 2/3).

`load_segment_alleles` reads a reference file via the appropriate
loader, runs the AnchorResolver per record, populates a BuildReport,
and returns the standard ``Dict[gene_name, List[Allele]]`` structure
the existing builder code consumes.

This is deliberately a thin module — the heavy lifting is on the C
side; here we just orchestrate the per-record decisions.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Optional

from ..alleles.allele import VAllele, DAllele, JAllele, CAllele
from .._native._anchor import (
    AnchorConfidence,
    AnchorMethod,
    FunctionalStatus,
    Locus,
    Segment,
    ReferenceLoader,
    resolve_anchor,
)


# Map "V"/"D"/"J"/"C" string segments to (Allele class, Segment enum).
_SEGMENT_MAP = {
    "V": (VAllele, Segment.V),
    "D": (DAllele, Segment.D),
    "J": (JAllele, Segment.J),
    "C": (CAllele, Segment.C),
}


def _status_in_policy(status: FunctionalStatus,
                      *,
                      include_orf: bool,
                      include_pseudo: bool,
                      include_partial: bool) -> bool:
    """T2-9: returns True if a record with the given status passes
    the filter policy. F is always kept; F/ORF/P/partial controlled
    by flags; UNKNOWN (no status info available, e.g. plain FASTA)
    is always kept — caller can set status manually if desired."""
    if status == FunctionalStatus.F:
        return True
    if status == FunctionalStatus.ORF:
        return include_orf
    if status == FunctionalStatus.PSEUDO:
        return include_pseudo
    if status == FunctionalStatus.PARTIAL:
        return include_partial
    # FUNC_UNKNOWN — loader didn't surface a status; always keep.
    return True


def _open_loader(path: str,
                 *,
                 format: Optional[str],
                 locus_hint: Optional[Locus],
                 segment_hint: Segment,
                 sidecar_path: Optional[str]) -> ReferenceLoader:
    """Dispatch based on explicit format or auto-detection.

    Format names: "imgt" / "airrc" / "ogrdb" / "igblast" / "plain" / None
    (None means auto-detect)."""
    locus = locus_hint or Locus.UNKNOWN

    if format == "imgt":
        return ReferenceLoader.open_imgt(path, segment_hint=segment_hint)
    if format == "airrc":
        return ReferenceLoader.open_airrc(path)
    if format == "ogrdb":
        return ReferenceLoader.open_ogrdb(
            path, sidecar_path,
            locus_hint=locus, segment_hint=segment_hint)
    if format == "igblast":
        return ReferenceLoader.open_igblast(
            path, sidecar_path,
            locus_hint=locus, segment_hint=segment_hint)
    if format == "plain":
        return ReferenceLoader.open_plain(
            path, locus_hint=locus, segment_hint=segment_hint)
    if format is None:
        return ReferenceLoader.open_auto(
            path, locus_hint=locus, segment_hint=segment_hint)
    raise ValueError(
        f"Unknown reference format {format!r}. Choose from: "
        f"imgt, airrc, ogrdb, igblast, plain, or None for auto-detect.")


def load_segment_alleles(path: str,
                         *,
                         segment: str,                      # "V" / "D" / "J" / "C"
                         format: Optional[str] = None,
                         locus_hint: Optional[Locus] = None,
                         sidecar_path: Optional[str] = None,
                         keep_anchorless: bool = False,
                         strict: bool = False,
                         include_orf: bool = True,
                         include_pseudo: bool = False,
                         include_partial: bool = False,
                         report=None) -> dict:
    """Load all alleles of a single segment from a reference file.

    Returns the standard ``{gene_name: [Allele, ...]}`` dict the
    existing builder code expects, with anchors fully resolved by
    the C-side AnchorResolver.

    Filter policy (T2-9):
      - F alleles always kept.
      - ORF alleles kept iff ``include_orf=True`` (default).
      - Pseudogenes (P) kept iff ``include_pseudo=True`` (default False).
      - Partial sequences kept iff ``include_partial=True`` (default False).
      - UNKNOWN-status records (typically plain FASTA / IgBLAST FASTA
        without sidecar) always kept — caller's responsibility to
        annotate if filtering is needed.

    Rejection policy (anchor resolver):
      - V or J with `keep_anchorless=False`: the allele is DROPPED
        and recorded in `report.rejections`.
      - V or J with `keep_anchorless=True`: the allele is kept with
        anchor=-1; downstream simulator code skips them in the
        productive retry loop.
      - D or C: always kept (no anchor concept for these segments).
    """
    if segment not in _SEGMENT_MAP:
        raise ValueError(f"segment must be V/D/J/C, got {segment!r}")
    allele_cls, segment_enum = _SEGMENT_MAP[segment]

    out = defaultdict(list)
    loader = _open_loader(
        path, format=format, locus_hint=locus_hint,
        segment_hint=segment_enum, sidecar_path=sidecar_path)

    with loader:
        for rec in loader:
            gene = rec.name.split("*")[0]

            # T2-9: apply F/ORF/P/partial filter BEFORE anchor
            # resolution. Filtering first keeps the BuildReport's
            # accept/reject counts focused on anchor outcomes for
            # records the user actually wanted to load.
            if not _status_in_policy(rec.functional_status,
                                     include_orf=include_orf,
                                     include_pseudo=include_pseudo,
                                     include_partial=include_partial):
                if report is not None:
                    report.record_filtered(rec)
                continue

            # D and C segments don't have anchors — pass through directly.
            if segment_enum in (Segment.D, Segment.C):
                # Use the gapped sequence if available (preserves
                # IMGT alignment); otherwise the ungapped one.
                seq = rec.gapped_sequence or rec.sequence
                a = allele_cls(rec.name, seq, len(rec.sequence))
                _attach_metadata(a, rec)
                out[gene].append(a)
                if report is not None:
                    # No anchor result for D/C — synthesize a
                    # CONF_CANONICAL placeholder so the report's
                    # confidence histogram stays meaningful.
                    from .._native._anchor import AnchorResult, AnchorResidue
                    fake = AnchorResult(
                        position=-1, codon="???",
                        residue=AnchorResidue.UNKNOWN,
                        confidence=AnchorConfidence.CANONICAL,
                        method=AnchorMethod.NONE,
                        reason=None)
                    report.record_accepted(rec, fake)
                continue

            # V or J: run the resolver.
            result = resolve_anchor(rec,
                                    segment=segment_enum,
                                    locus=rec.locus,
                                    strict=strict)

            if result.confidence == AnchorConfidence.REJECTED:
                if report is not None:
                    report.record_rejected(rec, result)
                if not keep_anchorless:
                    continue
                anchor = -1
            else:
                if report is not None:
                    report.record_accepted(rec, result)
                anchor = result.position

            seq = rec.gapped_sequence or rec.sequence
            a = allele_cls(rec.name, seq, len(rec.sequence),
                           anchor_override=anchor)
            # Surface the AnchorResult on the Allele for diagnostics.
            try:
                a.anchor_meta = result
            except AttributeError:
                pass  # if class declares __slots__ without anchor_meta
            _attach_metadata(a, rec)
            out[gene].append(a)

    return dict(out)


def _attach_metadata(allele, rec) -> None:
    """Copy provenance fields from a LoadedAlleleRecord onto an
    Allele instance (T2-8 closure step 1 — preserve metadata that
    was previously dropped on the way into VAllele/JAllele/DAllele/
    CAllele)."""
    try:
        allele.functional_status = rec.functional_status
        allele.locus = rec.locus
        allele.aliases = rec.aliases   # already a tuple of str
        allele.species = rec.species
        allele.source = rec.source
    except AttributeError:
        # Class declares __slots__ without these names — silently
        # tolerate. Today's Allele subclasses don't use __slots__.
        pass
