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

from dataclasses import dataclass, field as _field

ROOT = "engine_rs/tests/metadata_truth"


@dataclass(frozen=True)
class RustSource:
    """One Rust category file, and the mechanisms its docstring claims."""

    path: str
    scope: str
    mechanism_ids: tuple[str, ...] = ()
    fields: dict = _field(default_factory=dict)


RUST_SOURCES = (
    RustSource(
        path=f"{ROOT}/trim_ambiguity.rs",
        scope="Trimming can leave too few bases to distinguish the truth allele, "
              "so the call must surface every indistinguishable allele.",
        mechanism_ids=("germline.allele_choice", "recombination.trimming"),
        fields={"germline.allele_choice": ("v_call", "d_call", "j_call"),
                "recombination.trimming": ("v_trim_3", "d_trim_5", "d_trim_3", "j_trim_5")},
    ),
    RustSource(
        path=f"{ROOT}/np_extension.rs",
        scope="NP evidence may extend an inferred placement; the standard calls "
              "keep the generating assignments and every generated NP byte stays "
              "inside its actual junction interval.",
        mechanism_ids=("junction.n_addition", "germline.allele_choice"),
        fields={"junction.n_addition": ("np1", "np2", "np1_length", "np2_length"),
                "germline.allele_choice": ("v_call", "d_call", "j_call")},
    ),
    RustSource(
        path=f"{ROOT}/cross_segment.rs",
        scope="One segment's suffix matching the next segment's prefix lets a "
              "live call extend past the structural boundary.",
        mechanism_ids=("recombination.assembly", "germline.allele_choice"),
        fields={"recombination.assembly": ("v_germline_start", "v_germline_end",
                                           "d_germline_start", "d_germline_end",
                                           "j_germline_start", "j_germline_end"),
                "germline.allele_choice": ("v_call", "d_call", "j_call")},
    ),
    RustSource(
        path=f"{ROOT}/mutation_effects.rs",
        scope="SHM and PCR / quality / contaminant base edits flip bases in the "
              "pool; the call set must reflect post-mutation evidence.",
        mechanism_ids=("maturation.hypermutation", "assay.pcr_error",
                       "assay.quality_error", "assay.ambiguous_calls",
                       "germline.allele_choice"),
        fields={"maturation.hypermutation": ("n_mutations", "n_v_mutations",
                                             "n_d_mutations", "n_j_mutations"),
                "assay.pcr_error": ("n_pcr_errors",),
                "assay.quality_error": ("n_quality_errors",),
                "assay.ambiguous_calls": ("sequence",),
                "germline.allele_choice": ("v_call", "d_call", "j_call")},
    ),
    RustSource(
        path=f"{ROOT}/indel_effects.rs",
        scope="Indels shift the whole pool layout below the affected position; "
              "coordinates must follow and the truth allele must be retained.",
        mechanism_ids=("assay.polymerase_indels", "recombination.assembly"),
        fields={"assay.polymerase_indels": ("n_indels", "n_v_indels", "n_d_indels",
                                            "n_j_indels"),
                "recombination.assembly": ("sequence", "sequence_length")},
    ),
    RustSource(
        path=f"{ROOT}/junction_anchor.rs",
        scope="Junction bounds, translation, in-frame and stop-codon flags and "
              "the productive flag stay self-consistent under every corruption.",
        mechanism_ids=("junction.assembly", "junction.coding_end_processing",
                       "expression.productivity"),
        fields={"junction.assembly": ("junction", "junction_aa", "junction_length"),
                "junction.coding_end_processing": ("junction_start", "junction_end"),
                "expression.productivity": ("productive", "stop_codon", "vj_in_frame")},
    ),
    RustSource(
        path=f"{ROOT}/corruption_stack.rs",
        scope="Every corruption source applied to one record at once; bugs hide "
              "in the interactions rather than in any single source.",
        mechanism_ids=("assay.end_loss", "assay.orientation", "assay.contamination"),
        fields={"assay.end_loss": ("end_loss_5_length", "end_loss_3_length"),
                "assay.orientation": ("rev_comp",),
                "assay.contamination": ("is_contaminant",)},
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


def rust_covered_fields() -> dict[str, tuple[str, ...]]:
    """mechanism id -> the truth fields the Rust suite covers."""
    out: dict[str, set[str]] = {}
    for source in RUST_SOURCES:
        for mechanism_id, fields in source.fields.items():
            out.setdefault(mechanism_id, set()).update(fields)
    return {k: tuple(sorted(v)) for k, v in out.items()}
