"""Verification that already exists in the Rust suites.

``engine_rs/tests/metadata_truth/`` holds the same charter as this layer and
has done so for longer: 86 tests across eight category files, whose README
states the guarantee as "metadata that exactly reflects what a perfect aligner
/ biologist / oracle would extract from the sequence given only its bases".

Each record below is derived from the category file's own module docstring, not
from reading its individual tests, and claims only what that docstring states
it covers. The mechanism and field claims are validated against the catalog, so
a record cannot credit a mechanism with a field it does not own.
"""
from __future__ import annotations

from dataclasses import dataclass

ROOT = "engine_rs/tests/metadata_truth"


@dataclass(frozen=True)
class RustSource:
    """One Rust category file, and the mechanisms its docstring claims."""

    path: str
    scope: str
    mechanism_ids: tuple[str, ...] = ()
    # No field list: which fields a file covers is derived from the file in
    # `rust_covered_fields`. The hand-written table it replaced had 17 claims
    # naming a field absent from the crediting file.


RUST_SOURCES = (
    RustSource(
        path=f"{ROOT}/trim_ambiguity.rs",
        scope="Trimming can leave too few bases to distinguish the truth allele, "
              "so the call must surface every indistinguishable allele.",
        mechanism_ids=("germline.allele_choice", "recombination.trimming"),
    ),
    RustSource(
        path=f"{ROOT}/np_extension.rs",
        scope="NP evidence may extend an inferred placement; the standard calls "
              "keep the generating assignments and every generated NP byte stays "
              "inside its actual junction interval.",
        mechanism_ids=("junction.n_addition", "germline.allele_choice"),
    ),
    RustSource(
        path=f"{ROOT}/cross_segment.rs",
        scope="One segment's suffix matching the next segment's prefix lets a "
              "live call extend past the structural boundary.",
        mechanism_ids=("recombination.assembly", "germline.allele_choice"),
    ),
    RustSource(
        path=f"{ROOT}/mutation_effects.rs",
        scope="SHM and PCR / quality / contaminant base edits flip bases in the "
              "pool; the call set must reflect post-mutation evidence.",
        mechanism_ids=("maturation.hypermutation", "assay.pcr_error",
                       "assay.quality_error", "assay.ambiguous_calls",
                       "germline.allele_choice"),
    ),
    RustSource(
        path=f"{ROOT}/indel_effects.rs",
        scope="Indels shift the whole pool layout below the affected position; "
              "coordinates must follow and the truth allele must be retained.",
        mechanism_ids=("assay.polymerase_indels", "recombination.assembly"),
    ),
    RustSource(
        path=f"{ROOT}/junction_anchor.rs",
        scope="Junction bounds, translation, in-frame and stop-codon flags and "
              "the productive flag stay self-consistent under every corruption.",
        mechanism_ids=("junction.assembly", "junction.coding_end_processing",
                       "expression.productivity"),
    ),
    RustSource(
        path=f"{ROOT}/corruption_stack.rs",
        scope="Every corruption source applied to one record at once; bugs hide "
              "in the interactions rather than in any single source.",
        mechanism_ids=("assay.end_loss", "assay.orientation", "assay.contamination"),
    ),
)

#: Cross-cutting: the alignment-truth surface is the projection of whatever
#: every mechanism did, so no single mechanism owns it. Recorded separately
#: rather than attributed, so nothing claims ownership it does not have.
PROJECTION_SOURCES = (
    RustSource(
        path=f"{ROOT}/airr_projection.rs",
        scope="The live-call to AIRR-record projection path, across 69 named "
              "fields: fallbacks, empty cases and format rules.",
    ),
)


def all_rust_sources() -> tuple[RustSource, ...]:
    """Every recorded Rust source that claims a mechanism."""
    return RUST_SOURCES


def rust_coverage() -> dict[str, tuple[str, ...]]:
    """mechanism id -> the Rust source paths covering it."""
    out: dict[str, list[str]] = {}
    for source in RUST_SOURCES:
        for mechanism_id in source.mechanism_ids:
            out.setdefault(mechanism_id, []).append(source.path)
    return {k: tuple(v) for k, v in out.items()}


def _file_text(path: str) -> str:
    from tests.support import paths
    candidate = paths.repo_root() / path
    return candidate.read_text(encoding="utf-8") if candidate.exists() else ""


def _mentions(text: str, field: str) -> bool:
    """Whether ``field`` appears in ``text`` as a whole word.

    Word-bounded on purpose: a substring test credits ``np_extension.rs`` with
    ``j_call`` on the strength of ``j_call_names``, which is a different
    identifier.
    """
    import re
    return re.search(rf"\b{re.escape(field)}\b", text) is not None


def rust_covered_fields() -> dict[str, tuple[str, ...]]:
    """mechanism id -> the truth fields the Rust suite covers.

    Derived from what each file actually mentions, not from a hand-written
    table. The table was transcribed from module docstrings and 17 of its 51
    claims named a field that appears nowhere in the file crediting it -- six
    in cross_segment.rs, five in mutation_effects.rs, three in
    indel_effects.rs, two in corruption_stack.rs, one in np_extension.rs. None
    of the three guards over this module read the files, so crediting
    trim_ambiguity.rs with the end-loss fields passed all of them.

    A mention is weaker evidence than a reading: a file that names a field is
    not thereby proved to verify it. That is the price of a claim that cannot
    drift, and it is why this column still cannot move a status.
    """
    from ._stages import MECHANISMS
    owned = {m.id: m.truth_fields for m in MECHANISMS}
    out: dict[str, set[str]] = {}
    for source in RUST_SOURCES:
        text = _file_text(source.path)
        if not text:
            continue
        for mechanism_id in source.mechanism_ids:
            for field in owned.get(mechanism_id, ()):
                if _mentions(text, field):
                    out.setdefault(mechanism_id, set()).add(field)
    return {k: tuple(sorted(v)) for k, v in out.items()}
