"""Interpretation profiles: what a provider's coordinate fields mean, and on what evidence.

A profile is chosen explicitly (it is named in the source lock's adapter reference) and is
never selected or switched from the data.  Exactly one profile is implemented:
``ogrdb_api_v2_human_ig_observed_v1``, the conventions observed in the pinned OGRDB API v2
human IG exports (IGH_VDJ release 10, IGKappa_VJ release 5, IGLambda_VJ release 4; the exact
validated bytes are listed in :data:`~GenAIRR._reference_sources.common.PROFILE_SCOPES`).
It is GenAIRR's reading of what that producer emits, not a claim that every API v2 export
follows the AIRR 2.0 prose; where the two differ the producer's code decides.

The profile records, as data:

* the coordinate rules and the producer source lines and documents each rests on (pinned
  by commit, git blob and SHA-256; nothing is read at run time);
* the landmark residues used only to *report* noncanonical codons (never to place one);
* the conventions that remain unresolved, which stay unresolved in every inventory:
  J phase independence, V phase independence beyond anchor derivation, the meaning of a
  Boolean ``functional``, physical coding-end completeness, historical affirmed spans
  and RSS coordinates;
* the historical curation notes known to limit an affirmed span, each bound to its exact
  text and to the coding sequences it was read against, with no reviewed interpretation.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping, Optional

from GenAIRR._reference_sources._fields import bad, require_type
from GenAIRR._reference_sources.common import (
    COORDINATE_PROFILE_ID, COORDINATE_PROFILE_VERSION, PROFILE_SCOPES, CoordinateProfileRef,
)


@dataclass(frozen=True)
class PinnedSource:
    """A producer source file at an exact commit, with the line ranges relied on."""

    repository: str
    commit: str
    path: str
    git_blob: str
    sha256: str
    lines: tuple            # ((line range, what it establishes), ...)


@dataclass(frozen=True)
class PinnedDocument:
    """A documentation page by the digest of the exact bytes retrieved."""

    id: str
    url: str
    sha256: str
    bytes: int


@dataclass(frozen=True)
class ProfileEvidence:
    repository: str
    commit: str
    files: tuple            # PinnedSource of the producer, at ``commit``
    related_tool: PinnedSource
    documents: tuple        # PinnedDocument


@dataclass(frozen=True)
class NoteInterpretation:
    """A reviewed reading of a historical note: the note affirms gapped columns 1 to
    ``last_affirmed_column`` (1-based, inclusive) of the record's own gapped alignment."""

    last_affirmed_column: int


@dataclass(frozen=True)
class HistoricalNote:
    """A curation note that limits the historically affirmed span of one allele.

    Bound to the exact text of the record field (``text_sha256``) and to the coding
    sequence of each record export it was read against (``sequences``: export ->
    sequence digest).  ``numbering`` is the numbering its positions use.  Without a
    reviewed ``interpretation`` no affirmed interval is ever derived from it.
    """

    set_id: str
    release: int
    label: str
    field: str
    text_sha256: str
    quotes: tuple
    numbering: str
    sequences: tuple        # ((export, sequence sha256), ...)
    interpretation: Optional[NoteInterpretation]


@dataclass(frozen=True)
class CoordinateProfile:
    """One interpretation profile (see the module documentation)."""

    id: str
    version: int
    rules: Mapping
    basis: Mapping
    landmark_residues: Mapping
    v_delineation_scheme: str
    missing_delineation_value: int
    landmark_columns: tuple           # 1-based inclusive gapped columns of IMGT codon 104
    j_phase_independence: str
    v_phase_derivation: str
    functionality_classes: Mapping
    unresolved: tuple
    historical_notes: tuple
    evidence: ProfileEvidence

    @property
    def ref(self) -> CoordinateProfileRef:
        return CoordinateProfileRef(self.id, self.version)

    @property
    def scope(self) -> tuple:
        """The exact sources the profile has been validated against (the frozen registry entry)."""
        return PROFILE_SCOPES[(self.id, self.version)]


# ── ogrdb_api_v2_human_ig_observed_v1 ───────────────────────────────

_OGRDB = "https://github.com/airr-community/ogrdb"
_OGRDB_COMMIT = "c572507728789eb6436d2747a872ce434077ae5f"
_TO_AIRR = "ogrdb/germline_set/to_airr.py"
_API = "api_v2/api.py"
_MODELS = "api_v2/models.py"
_ROUTES = "ogrdb/sequence/sequence_routes.py"
_FASTA = "ogrdb/germline_set/descs_to_fasta.py"
_MAINT = "ogrdb/maint/maint_routes.py"
_DB = "db/gene_description_db.py"
_TOOL = "src/receptor_utils/aux_formats.py"


def _ogrdb(path: str, git_blob: str, sha256: str, lines: tuple) -> PinnedSource:
    return PinnedSource(_OGRDB, _OGRDB_COMMIT, path, git_blob, "sha256:" + sha256, lines)


_EVIDENCE = ProfileEvidence(
    repository=_OGRDB,
    commit=_OGRDB_COMMIT,
    files=(
        _ogrdb(_TO_AIRR, "ff193fbebee416d07fbdd2074c6f2f7545b66d5e",
               "3f849c7f3bb0c72dc9dff91064334c17ee496b0779474b7ff81aab7eb5cf0fb2", (
                   ("108-110", "coding_sequence is the gapped coding sequence without gaps, plus the 3' extension "
                               "in '_ex' exports"),
                   ("119", "functional is whether the stored functionality text contains 'F'"),
                   ("143-155", "the IMGT delineation's unaligned sequence is the coding sequence"),
                   ("166", "cdr3_start is converted to a 1-based coordinate in the coding sequence"),
                   ("179-180", "j_codon_frame and j_cdr3_end are passed through unchanged"))),
        _ogrdb(_API, "c2c897ba143d8a9a9396c07db7daebfcc647653a",
               "7389f2ef709f9f285ed5bc10637f9f6a6b9b4af18dc00e27bd3d3f9f1f5cda19", (
                   ("339-343", "'_ex' formats are built with extensions"),
                   ("349", "JSON exports are built by germline_set_to_airr"),
                   ("360", "FASTA exports are built by descs_to_fasta"),
                   ("542-578", "served delineations carry the numeric fields only"),
                   ("560-570", "a delineation coordinate stored as None is served as 0"))),
        _ogrdb(_MODELS, "725de65e18fdd4e161e8d532fbe83b678b760332",
               "95fafa1c24b5bbe39a6b6a51123a7e56891a9f5d0df8ba49b4318996c6600448", (
                   ("246-247", "cdr3_start is a coordinate in the unaligned sequence"),)),
        _ogrdb(_ROUTES, "dc925622a3663152b66cb0241003422fe19bd333",
               "b2b996b96597a69cae895376b99389b5c3e490c6c9d883c4e1d814beef4a4802", (
                   ("1821-1825", "default IMGT ranges over the gapped sequence; FR3 ends at column 312"),
                   ("1828-1868", "cdr3_start = 1 + non-gap characters before gapped column 313"),
                   ("1203", "bulk upload stores cdr3_start as a full-sequence coordinate"),
                   ("1107-1108", "uploaders can set the gapped feature ranges"),
                   ("2085-2090", "cdr3_start can be edited by hand"))),
        _ogrdb(_FASTA, "8fbb09a2d3e825fff13b161d9f363280ed273571",
               "8e8cc5991b2ba9b59f1b32aea76d0a4419a5107fd61e5202a0ae851fb15b4bb0", (
                   ("87-91", "FASTA relabelling of names without '*' (no fixture label affected)"),
                   ("98-105", "gapped exports write the gapped coding sequence, plus the extension for '_ex'"))),
        _ogrdb(_MAINT, "9a92570e6c0b2e791415685fe36f3264340da156",
               "2bc46f504435f55d7f340e5458e573200b925ff3907b1df3c660db79f7ba2d8c", (
                   ("503-509", "a maintenance check expects cdr3_start at aligned position 313"),
                   ("552-562", "aligned positions 310-312 are checked as the second cysteine"))),
        _ogrdb(_DB, "63205b5b736bc774386820966e31b4db4ed2212c",
               "48cc37fafc8002f19f0fb7162b41a87458e5cc2db9eb15653a5afe8dd60081f2", (
                   ("394", "stored cdr3_start is a full-sequence coordinate"),)),
    ),
    related_tool=PinnedSource(
        "https://github.com/williamdlees/receptor_utils", "e7983d47e6662449db100940c28ef83491c4f51e", _TOOL,
        "cabb6203799cf3b4e11f3b84cf6d7ea625534ee3",
        "sha256:d2da33d0e5c216b4d82d580dd397287db5b55fa6406241dece8646c21f2eb545", (
            ("163", "J records need (j_cdr3_end - j_codon_frame) % 3 == 0"),
            ("176", "coding frame = j_codon_frame - 1"),
            ("178", "zero-based last CDR3 nucleotide = j_cdr3_end - 2"))),
    documents=(
        PinnedDocument("ogrdb_api_v2_guide",
                       "https://raw.githubusercontent.com/airr-community/ogrdb/"
                       "c572507728789eb6436d2747a872ce434077ae5f/schema/ogrdb_api_v2_guide.md",
                       "sha256:ea7d9c78b34916cd8103a0d130a1c23453b0b8909794fb70f1cb4b072e35d8a5", 6990),
        PinnedDocument("ogrdb_germline_set_formats",
                       "https://wordpress.vdjbase.org/index.php/ogrdb/explanation-of-germline-set-formats/",
                       "sha256:9201dceff16a465338e3ac667a49648299e9fef1de8536fcbd5120ff26fb3a5b", 57094),
        PinnedDocument("airr_germline_schema", "https://docs.airr-community.org/en/stable/datarep/germline.html",
                       "sha256:8b83d03e62cbba14f2f715b7bc861207c75dd8802cd41b519433da803f34a54b", 82899),
        PinnedDocument("ogrdb_openapi", "https://ogrdb.airr-community.org/schema/ogrdb-api-openapi3.yaml",
                       "sha256:c06ff8c2e0790825956fdaa3a7b43e74ef7bad2777f05846a5d5012035586cb3", 10258),
        PinnedDocument("ogrdb_licensing_statement",
                       "https://ogrdb.airr-community.org/render_page/licensing_statement.html",
                       "sha256:1cf20853b6302f5ea3d55f322592821c513842f79ba8b61604931a2cdec7dbe7", 27047),
    ),
)

# The rules, in the notation of the producer's fields; coordinates.py implements them.
_RULES = {
    "j_anchor_coding0": "j_cdr3_end - 1",
    "j_reported_phase": "j_codon_frame - 1",
    "v_anchor_coding0": "cdr3_start - 4",
    "v_cdr3_start_from_gapped": "1 + count of non-gap characters in gapped columns 1..312",
    "coding_interval": "sequence[gene_start - 1 : gene_end]",
}
_BASIS = {
    "j_anchor_coding0": ((_TO_AIRR, "179-180"), (_TOOL, "178")),
    "j_reported_phase": ((_TO_AIRR, "179-180"), (_TOOL, "163"), (_TOOL, "176")),
    "v_anchor_coding0": ((_TO_AIRR, "143-155"), (_TO_AIRR, "166"), (_MODELS, "246-247"), (_ROUTES, "1828-1868")),
    "v_cdr3_start_from_gapped": ((_ROUTES, "1821-1825"), (_ROUTES, "1828-1868"), (_ROUTES, "1107-1108"),
                                 (_ROUTES, "2085-2090"), (_MAINT, "552-562")),
    "coding_interval": ((_TO_AIRR, "108-110"), (_API, "339-343")),
    "gapped_alignment": ((_FASTA, "87-91"), (_FASTA, "98-105"), (_API, "360")),
    "missing_delineation_value": ((_API, "560-570"),),
}
_UNRESOLVED = (
    ("v_coding_phase_independence", "No source field states a V coding phase; a phase from the anchor is "
                                    "anchor-derived, never independent."),
    ("j_phase_derivation", "j_codon_frame converts to a reported phase value; whether the producer derived it "
                           "independently of j_cdr3_end is not documented."),
    ("functional_false_meaning", "functional is whether the stored functionality text contains 'F'; it does not "
                                 "classify an allele as functional, ORF or pseudogene."),
    ("coding_end_completeness", "Curational tags, functional, RSS adjacency and exact context slices do not "
                                "establish complete physical coding ends; '_ex' extensions are provider "
                                "inferences."),
    ("historical_affirmed_spans", "Affirmation limits in curation notes use gapped numbering and have no "
                                  "reviewed interpretation."),
    ("v_rss_coordinates", "RSS coordinates are kept as observations only; some contradict the record's coding "
                          "interval."),
    ("server_revision", "The deployed producer revision is not observable; each record is checked against its "
                        "own exports instead."),
)


def _note(label: str, text_sha256: str, quotes: tuple, airr: str, airr_ex: str) -> HistoricalNote:
    return HistoricalNote(set_id="9606.IGLambda_VJ", release=4, label=label, field="curation",
                          text_sha256="sha256:" + text_sha256, quotes=quotes, numbering="imgt_gapped",
                          sequences=(("airr", "sha256:" + airr), ("airr_ex", "sha256:" + airr_ex)),
                          interpretation=None)


_UNCERTAIN_TWO = ("Plots of the final 3\u2019 nucleotides were considered, but the variability seen made it "
                  "impossible to consider the final two nucleotides of the sequence.")
_NOTES = (
    _note("IGLV10-54*04", "12408e2acae180953310dd6010927ba2abd467b8d556556d97f9c2e1504a5f70", (
        "Plots of the final 3\u2019 nucleotides showed a high level of variability, making it impossible to "
        "determine the final three nucleotides with certainty.",
        "The sequence, up to and including nucleotide 339, was affirmed as a Level 1 sequence.",
        "Uncertainty regarding nucleotides 340-341 will be indicated in IARC and OGRDB publications by two dots at "
        "the end of the affirmed sequence."),
        "054d5b997ca16c9313b56e814b1604ce669dd2d5d42b8d93922a7686d1b8eec1",
        "03800d1a35cc28cd4586817b201cb54f09dbefe54b17dd15e73ca059675d7c68"),
    _note("IGLV2-14*03", "3c06429abcbabf8c891cb8c95db7e2970060b8f4e94a60afe38d61220e20dafe", (
        "There was support for the sequence being affirmed up to and including nucleotide 337, in line with IARC "
        "policy not to infer bases substantially affected by trimming during the rearrangement process, as a "
        "Level 1 sequence, with two additional 3\u2019 nucleotide being likely present in the sequence.",),
        "5ca70268a23d9f1ba3e8381a3d9d1cd95743f6f6e3a7a7136057305e8d0db9d9",
        "a09488a27413bd19d5366eedd5c07dfe124e7ec952af381551427d1e1690d71e"),
    _note("IGLV3-21*04", "d72b8ac22784b8260deb5bcccb4b49de0f34e7a4424b727aac6e69b0c3b0044e", (
        _UNCERTAIN_TWO,
        "The sequence, up to and including nucleotide 339, was affirmed as the Level 1 sequence, IGLV3-21*i01.",
        "Uncertainty regarding nucleotides 340-341 will be indicated in IARC and OGRDB publications by two dots at "
        "the end of the affirmed sequence."),
        "d0b983b1223e55be23ae8339d1e5060814900433fa90fd9a8d1e7cb4f3c8d285",
        "c1207f8ac1954e1610e6b6dbf0b7bc38301594c032dfc32ce293dc18f60c2b92"),
    _note("IGLV3-25*03", "73959a3aa995cbb0cb6f807584e709cf441f3745f0a06d8f2252196e5b80ad44", (
        _UNCERTAIN_TWO,
        "The sequence, up to and including nucleotide 339, was affirmed as the Level 1 sequence, IGLV3-25*i01."),
        "47c8e1d37a1beed3957c0d059baf17113bd146f473d40ba4cd46b270dba43213",
        "59b875956dd72a662d02f278efb5a0f9537252de8ee84f09f6cb57d76e671d44"),
    _note("IGLV6-57*04", "ba5d8ad8612f577e490f97bf07f874f3c4414d93454ec2bc602beb08259e57a3", (
        "Plots of the final 3\u2019 nucleotides were considered, but the variability seen made it impossible to be "
        "certain of the final two nucleotides of the germline sequence.",
        "The sequence, up to and including nucleotide 333, was affirmed as the Level 1 sequence, IGLV6-57*i01.",
        "Uncertainty regarding nucleotides 334-335 will be indicated in IARC and OGRDB publications by two dots at "
        "the end of the affirmed sequence."),
        "b7961ecc7ec1e4cd2707b66aaeb34e6bc24cb62e25f6f14245eb0ad3d83fe7eb",
        "9e9dcda340a220c8a7e4a15d71b9640b520ca1c6d2350fb35e5f989f6c6b886d"),
)

PROFILE = CoordinateProfile(
    id=COORDINATE_PROFILE_ID,
    version=COORDINATE_PROFILE_VERSION,
    rules=MappingProxyType(_RULES),
    basis=MappingProxyType(_BASIS),
    landmark_residues=MappingProxyType({"V": {"IGH": ("C",), "IGK": ("C",), "IGL": ("C",)},
                                        "J": {"IGH": ("W",), "IGK": ("F",), "IGL": ("F",)}}),
    v_delineation_scheme="IMGT",
    missing_delineation_value=0,
    landmark_columns=(310, 312),
    j_phase_independence="unresolved",
    v_phase_derivation="anchor_derived",
    functionality_classes=MappingProxyType({}),
    unresolved=_UNRESOLVED,
    historical_notes=_NOTES,
    evidence=_EVIDENCE,
)

# (profile id, version) -> profile.  Only profiles with reviewed evidence are listed.
PROFILES = MappingProxyType({(PROFILE.id, PROFILE.version): PROFILE})


def profile_for(ref: CoordinateProfileRef) -> CoordinateProfile:
    """The implemented profile named by ``ref`` (``unsupported_profile`` if none is)."""
    require_type(ref, CoordinateProfileRef, "ref")
    profile = PROFILES.get((ref.id, ref.version))
    if profile is None:
        raise bad("/coordinate_profile", f"interpretation profile {ref.id} version {ref.version} is not implemented "
                                         f"(implemented: {sorted(PROFILES)})", "unsupported_profile")
    return profile
