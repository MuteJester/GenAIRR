"""Coordinate rules of the OGRDB interpretation profile, as pure functions of field values.

Coordinate spaces are always named.  A *record sequence* position is 1-based in the
provider's fields (``gene_start``, ``gene_end``, ``v_rs_start``); a *coding sequence*
anchor is a zero-based index of the first nucleotide of the whole landmark codon, with all
three bases in bounds; *gapped columns* are IMGT alignment columns, 1-based in the
producer's numbering and zero-based half-open in :func:`gapped_to_coding`.  Intervals
returned here are zero-based and half-open.

Every function validates meaning first and reports a problem instead of guessing: a
missing, Boolean, non-integer or out-of-range value gives no position, and nothing is
ever moved to make a codon look canonical.  The AIRR 2.0 prose reading of the J fields
is deliberately not implemented: the producer's convention is the only one applied, and
no other reading is tried as a fallback.

* J: ``anchor = j_cdr3_end - 1`` in the coding sequence; the reported phase value is
  ``j_codon_frame - 1`` (a source-stated value; its independence from the landmark is
  not established here).
* V: ``anchor = cdr3_start - 4`` in the coding sequence, bound to the record's gapped
  alignment: the alignment ungaps to the coding sequence, the landmark codon at gapped
  columns 310-312 holds no gap, and ``1 +`` the non-gap count in columns 1-312 equals
  ``cdr3_start``.
* Coding interval: ``sequence[gene_start - 1 : gene_end]``, applied only where it
  reproduces the supplied coding sequence.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from GenAIRR._airr_record import _translate_codon
from GenAIRR._reference.model import sequence_digest
from GenAIRR._reference_sources.inventory import IMGT_GAP


class _Absent:
    """A field the record does not have (distinct from a present JSON null)."""

    def __repr__(self) -> str:
        return "ABSENT"


ABSENT = _Absent()

# IMGT codon 104 occupies gapped columns 310-312 (1-based); the CDR3 starts at column 313.
LANDMARK_FIRST_COLUMN = 310
LANDMARK_LAST_COLUMN = 312
# Zero-based half-open gapped columns of the landmark codon.
_LANDMARK = (LANDMARK_FIRST_COLUMN - 1, LANDMARK_LAST_COLUMN)

_JSON_TYPES = {dict: "object", list: "array", str: "string", int: "number", float: "number", type(None): "null"}


def residue(codon: str) -> str:
    """The amino acid of ``codon`` in the standard code (``X`` if it is not three of A/C/G/T)."""
    return _translate_codon(codon)


def _non_gap(text: str) -> int:
    return sum(1 for char in text if char != IMGT_GAP)


# ── stated positions ────────────────────────────────────────────────


@dataclass(frozen=True)
class Position:
    """An exact integer field value, or why there is none."""

    value: Optional[int]
    problem: Optional[str]


def stated_position(value: Any, *, minimum: Optional[int] = 1, zero_is_missing: bool = False) -> Position:
    """Read an integer coordinate field exactly: ``true`` and ``8.0`` are not positions.

    ``zero_is_missing``: the provider serves a missing value as 0 (delineation coordinates).
    """
    if value is ABSENT:
        return Position(None, "absent")
    if value is None:
        return Position(None, "null")
    if type(value) is bool:
        return Position(None, f"a Boolean ({str(value).lower()}), not an integer")
    if type(value) is not int:
        return Position(None, f"a JSON {_JSON_TYPES.get(type(value), type(value).__name__)}, not an integer")
    if zero_is_missing and value == 0:
        return Position(None, "0, which the provider serves for a missing coordinate")
    if minimum is not None and value < minimum:
        return Position(None, f"{value} is below {minimum}")
    return Position(value, None)


def _named(name: str, position: Position) -> Optional[str]:
    return None if position.problem is None else f"{name}: {position.problem}"


# ── J ───────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class JLandmark:
    """A J record's anchor (coding sequence, zero-based) and reported phase value."""

    anchor: Optional[int]
    anchor_problem: Optional[str]
    phase: Optional[int]
    phase_problem: Optional[str]
    codon: Optional[str]

    @property
    def consistent(self) -> Optional[bool]:
        """Whether the reported phase is the anchor modulo 3 (``None`` unless both are known)."""
        if self.anchor is None or self.phase is None:
            return None
        return self.phase == self.anchor % 3


def j_landmark(coding: str, j_cdr3_end: Any, j_codon_frame: Any) -> JLandmark:
    """Apply ``anchor = j_cdr3_end - 1`` and ``phase = j_codon_frame - 1`` to a J coding sequence."""
    anchor = codon = phase = None
    end = stated_position(j_cdr3_end)
    anchor_problem = _named("j_cdr3_end", end)
    if end.value is not None:
        if end.value - 1 + 3 > len(coding):
            anchor_problem = (f"j_cdr3_end: {end.value} puts the landmark codon at {end.value - 1}-{end.value + 1}, "
                              f"beyond the {len(coding)} bases of the coding sequence")
        else:
            anchor = end.value - 1
            codon = coding[anchor:anchor + 3]
    frame = stated_position(j_codon_frame, minimum=None)
    phase_problem = _named("j_codon_frame", frame)
    if frame.value is not None:
        if frame.value in (1, 2, 3):
            phase = frame.value - 1
        else:
            phase_problem = f"j_codon_frame: {frame.value} is not 1, 2 or 3"
    return JLandmark(anchor, anchor_problem, phase, phase_problem, codon)


# ── V ───────────────────────────────────────────────────────────────


def gapped_cdr3_start(gapped: str) -> Optional[int]:
    """``1 +`` the non-gap characters in gapped columns 1-312, or ``None`` if there are fewer columns."""
    if len(gapped) < LANDMARK_LAST_COLUMN:
        return None
    return 1 + _non_gap(gapped[:LANDMARK_LAST_COLUMN])


def gapped_to_coding(gapped: str, start: int, end: int) -> Optional[tuple]:
    """Map gapped columns ``[start, end)`` (zero-based) to ``[u(start), u(end))`` in the ungapped
    sequence, where ``u(i)`` counts the non-gap characters before column ``i``.  ``None``
    unless every column of the interval exists and holds a base (a gap or a column beyond
    the alignment is never guessed)."""
    if not 0 <= start < end <= len(gapped) or IMGT_GAP in gapped[start:end]:
        return None
    before = _non_gap(gapped[:start])
    return before, before + (end - start)


@dataclass(frozen=True)
class VLandmark:
    """A V record's anchor (coding sequence, zero-based), or why it is unresolved.

    ``kind`` names the failing step: ``field`` (``cdr3_start`` itself), ``alignment`` (no
    gapped record) or ``binding`` (the gapped record does not bind the number).
    """

    anchor: Optional[int]
    problem: Optional[str]
    kind: Optional[str]
    codon: Optional[str]
    gapped_cdr3_start: Optional[int]


def v_landmark(coding: str, cdr3_start: Any, gapped: Optional[str]) -> VLandmark:
    """Apply ``anchor = cdr3_start - 4``, bound to the record's own gapped alignment."""
    stated = stated_position(cdr3_start, zero_is_missing=True)
    if stated.value is None:
        return VLandmark(None, _named("cdr3_start", stated), "field", None, None)
    if stated.value < 4:
        return VLandmark(None, f"cdr3_start: {stated.value} is below 4, so the landmark codon would start before "
                               "the coding sequence", "field", None, None)
    anchor = stated.value - 4
    if anchor + 3 > len(coding):
        return VLandmark(None, f"cdr3_start: {stated.value} puts the landmark codon at {anchor}-{anchor + 2}, beyond "
                               f"the {len(coding)} bases of the coding sequence", "field", None, None)
    if gapped is None:
        return VLandmark(None, "no gapped record is bound to this record", "alignment", None, None)
    if gapped.replace(IMGT_GAP, "") != coding:
        return VLandmark(None, "the gapped record does not ungap to the coding sequence", "binding", None, None)
    found = gapped_cdr3_start(gapped)
    if found is None:
        return VLandmark(None, f"the gapped record has {len(gapped)} columns, fewer than "
                               f"{LANDMARK_LAST_COLUMN}", "binding", None, None)
    if gapped_to_coding(gapped, *_LANDMARK) is None:
        return VLandmark(None, f"the landmark codon at gapped columns {LANDMARK_FIRST_COLUMN}-{LANDMARK_LAST_COLUMN} "
                               "holds a gap", "binding", None, found)
    if found != stated.value:
        return VLandmark(None, f"cdr3_start {stated.value} differs from {found}, 1 + the non-gap characters in gapped "
                               f"columns 1-{LANDMARK_LAST_COLUMN}", "binding", None, found)
    return VLandmark(anchor, None, None, coding[anchor:anchor + 3], found)


# ── coding interval ─────────────────────────────────────────────────


@dataclass(frozen=True)
class CodingInterval:
    """The declared coding interval ``[start, end)`` of the record sequence and whether
    ``sequence[start:end]`` is the coding sequence; if not, how they relate."""

    start: Optional[int]
    end: Optional[int]
    problem: Optional[str]
    reproduces: Optional[bool]
    relation: Optional[str]
    extra: Optional[str]


def _relation(declared: str, coding: str) -> tuple:
    for longer, shorter, who in ((coding, declared, "candidate_extends_interval"),
                                 (declared, coding, "interval_extends_candidate")):
        if len(longer) > len(shorter):
            if longer.startswith(shorter):
                return f"{who}_3prime", longer[len(shorter):]
            if longer.endswith(shorter):
                return f"{who}_5prime", longer[:len(longer) - len(shorter)]
    return "different", None


def coding_interval(sequence: Any, gene_start: Any, gene_end: Any, coding: str) -> CodingInterval:
    """Apply ``sequence[gene_start - 1 : gene_end]`` and compare it with the coding sequence."""
    def unusable(problem):
        return CodingInterval(None, None, problem, None, None, None)
    if sequence is ABSENT or type(sequence) is not str:
        found = "absent" if sequence is ABSENT else f"a JSON {_JSON_TYPES.get(type(sequence), 'value')}, not text"
        return unusable(f"sequence: {found}")
    start, end = stated_position(gene_start), stated_position(gene_end)
    for name, position in (("gene_start", start), ("gene_end", end)):
        if position.problem is not None:
            return unusable(_named(name, position))
    if start.value > end.value:
        return unusable(f"gene_start {start.value} is after gene_end {end.value}")
    if end.value > len(sequence):
        return unusable(f"gene_end {end.value} is beyond the {len(sequence)} bases of the record sequence")
    declared = sequence[start.value - 1:end.value]
    if declared == coding:
        return CodingInterval(start.value - 1, end.value, None, True, None, None)
    relation, extra = _relation(declared, coding)
    return CodingInterval(start.value - 1, end.value, None, False, relation, extra)


# ── extensions and RSS ──────────────────────────────────────────────


@dataclass(frozen=True)
class ExtensionContext:
    """How a coding extension compares with the record bases at the same positions.

    ``status``: ``matches`` (every extension base is present and equal),
    ``disagrees`` (a present base differs) or ``without_record_context`` (the record ends
    before the extension does and every present base agrees).  ``bases`` are the record
    bases found at the extension's positions.
    """

    status: str
    bases: str


def extension_context(sequence: str, interval: tuple, side: str, extra: str) -> ExtensionContext:
    """Compare ``extra`` with the record bases beyond the declared ``interval`` on ``side``
    (``"3"``: after the interval; ``"5"``: before it)."""
    start, end = interval
    if side == "3":
        bases = sequence[end:end + len(extra)]
        expected = extra[:len(bases)]
    else:
        bases = sequence[max(0, start - len(extra)):start]
        expected = extra[len(extra) - len(bases):]
    if bases != expected:
        return ExtensionContext("disagrees", bases)
    return ExtensionContext("matches" if len(bases) == len(extra) else "without_record_context", bases)


@dataclass(frozen=True)
class RssPosition:
    """Where a declared V RSS starts relative to the declared coding interval.

    ``gap`` = ``v_rs_start - gene_end - 1`` (bases between them, negative when the RSS
    starts inside or before the interval); ``position``: ``after_coding_end``,
    ``inside_coding_interval`` or ``before_coding_interval``.
    """

    gap: Optional[int]
    position: Optional[str]
    problem: Optional[str]


def v_rss_position(gene_start: int, gene_end: int, v_rs_start: Any) -> RssPosition:
    """Compare a declared ``v_rs_start`` with a usable declared coding interval (1-based fields)."""
    stated = stated_position(v_rs_start)
    if stated.value is None:
        return RssPosition(None, None, _named("v_rs_start", stated))
    gap = stated.value - gene_end - 1
    if gap >= 0:
        position = "after_coding_end"
    elif stated.value >= gene_start:
        position = "inside_coding_interval"
    else:
        position = "before_coding_interval"
    return RssPosition(gap, position, None)


# ── historically affirmed spans ─────────────────────────────────────


@dataclass(frozen=True)
class AffirmedSpan:
    """The coding interval a historical note affirms, or why none is derived."""

    interval: Optional[tuple]
    reason: Optional[str]


def affirmed_span(note, export: str, coding: str, gapped: Optional[str]) -> AffirmedSpan:
    """Map a note's reviewed affirmed columns onto the coding sequence of one record.

    Only a reviewed interpretation, the exact sequence the note was read against and an
    exact gapped alignment of that sequence give an interval; anything else stays
    unresolved with its reason.
    """
    def unresolved(reason):
        return AffirmedSpan(None, reason)
    if note.interpretation is None:
        return unresolved("no reviewed interpretation of this note exists")
    if dict(note.sequences).get(export) != sequence_digest(coding):
        return unresolved(f"the {export} coding sequence is not the sequence this note was reviewed against")
    if gapped is None:
        return unresolved("no gapped record is bound to this record")
    if gapped.replace(IMGT_GAP, "") != coding:
        return unresolved("the gapped record does not ungap to the coding sequence")
    column = note.interpretation.last_affirmed_column
    if not 1 <= column <= len(gapped):
        return unresolved(f"affirmed column {column} is beyond the {len(gapped)} columns of the gapped record")
    return AffirmedSpan((0, _non_gap(gapped[:column])), None)
