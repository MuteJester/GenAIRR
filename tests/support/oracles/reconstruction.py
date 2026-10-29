"""Check the recorded sequence against the reference tables it claims to come from.

The provenance oracle next door derives the junction's untemplated content from
the engine's own record. That is useful but it cannot contradict the engine: if
provenance were corrupt, the oracle would agree with it. This module closes
that gap by reading the *reference* -- the allele sequences the simulation says
it drew from -- and asserting the recorded germline matches.

Taking the reference tables is not circular. Calling the production junction
builder or region projector would be; reading the allele a record names and
looking up the base at the coordinate it claims is the independent step.

What that step establishes is **allele identity**, and only that. It is not an
independent check on the coordinates, and must not be read as one: production
lays each germline byte down beside its coordinate from the same two numbers
(``Nucleotide::germline(slice[i], slice_start + i)``), so a trim that slices
one base too many shifts the byte and the coordinate together and every
comparison here still agrees. Measured: shifting both by one over 366
positions yields zero mismatches. The coordinates are checked by the
``coordinates`` oracle beside this one, which never compares a base and
instead requires the trim amounts, the end-loss lengths, the exported
coordinates and the reference allele's own length to add up.

Two further things this module knows that the provenance oracle does not:

- **P nucleotides are templated.** A P base is the reverse complement of the
  segment's own post-trim tip, so "no germline position" does not mean
  "untemplated" -- it means "not at a germline coordinate". P and N are
  distinguished in the output by letter case, lowercase for P. That is an
  undocumented truth channel which any consumer calling ``.upper()`` destroys.
- **A recorded substitution should change the base.** Comparing the final
  bases to the recorded germline bases gives the observable substitutions,
  which is what a mutation caller can actually recover.
"""
from __future__ import annotations

from dataclasses import dataclass

from . import provenance

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


@dataclass(frozen=True)
class Mismatch:
    """A position whose recorded germline base is not what the reference says."""

    index: int
    germline_position: int
    segment: str
    recorded: str
    reference: str


@dataclass(frozen=True)
class UntemplatedRun:
    """A contiguous stretch carrying no germline coordinate."""

    start: int
    bases: str
    kind: str                      # "P" | "N" | "mixed"

    @property
    def end(self) -> int:
        return self.start + len(self.bases)


def segments_by_position(simulation) -> dict[int, str]:
    """Position -> the segment whose region contains it."""
    out: dict[int, str] = {}
    for region in simulation.regions():
        segment = region.segment() if callable(region.segment) else region.segment
        start = region.start() if callable(region.start) else region.start
        end = region.end() if callable(region.end) else region.end
        for i in range(start, end):
            out[i] = segment
    return out


def _allele_sequence(refdata, simulation, segment: str) -> bytes | None:
    getter, allele_id = {
        "V": (refdata.v_allele, simulation.v_allele_id),
        "D": (refdata.d_allele, simulation.d_allele_id),
        "J": (refdata.j_allele, simulation.j_allele_id),
    }.get(segment, (None, None))
    if getter is None:
        return None
    identifier = allele_id() if callable(allele_id) else allele_id
    if identifier is None:
        return None
    allele = getter(identifier)
    seq = allele.seq() if callable(allele.seq) else allele.seq
    return seq


@dataclass(frozen=True)
class TemplatedCheck:
    """The result of comparing recorded germline bases to the reference.

    ``mismatches`` alone is ambiguous: an empty tuple means either "every
    position agreed" or "no position was examined". The second happens whenever
    a simulation carries no germline coordinates at all, or whose regions are
    empty, or whose only templated segment is an inverted D -- and a test that
    reads only the emptiness would report such a case as verified.

    So the count travels with the verdict, and callers assert on both.
    """

    mismatches: tuple[Mismatch, ...]
    checked: int
    skipped_inverted: int
    #: Positions whose segment names an allele the reference does not have.
    unresolved_allele: int = 0
    #: Positions whose germline coordinate falls outside the allele it names.
    #: Either is a defect in its own right, and both used to be a bare
    #: ``continue`` -- not compared, not counted, not reported. A reference
    #: truncated by thirty bases produced 0 mismatches, a lower ``checked``,
    #: and ``clean``; nothing asserted on ``checked``'s absolute value, so the
    #: count added here did not close this on its own.
    out_of_range: int = 0

    @property
    def clean(self) -> bool:
        """Something was compared, it all agreed, and nothing was unexaminable."""
        return (not self.mismatches and self.checked > 0
                and self.unresolved_allele == 0 and self.out_of_range == 0)


def templated_check(simulation, refdata, *, skip_inverted: bool = True) -> TemplatedCheck:
    """Compare every recorded germline base to the reference's base.

    This catches a wrong allele assembled with internally consistent
    provenance -- a failure the engine's own record cannot reveal, because its
    provenance would agree with itself.

    It does not catch a mis-trimmed segment. A coordinate is only visible here
    when it moves independently of the byte beside it, which a record-level
    mutant can do and a trimming defect cannot. Use the ``coordinates`` oracle
    for that.

    Inverted D is skipped by default: its coordinates are forward-reference
    while the sequence carries the reverse complement, so a direct lookup does
    not apply and comparing it needs its own claim.
    """
    germline = simulation.germline_bases()
    segments = segments_by_position(simulation)
    cache: dict[str, bytes | None] = {}
    out: list[Mismatch] = []
    checked = skipped = unresolved = out_of_range = 0
    for index in range(len(germline)):
        position = simulation.germline_position(index)
        if position is None:
            continue
        segment = segments.get(index, "?")
        if segment not in ("V", "D", "J"):
            continue
        if skip_inverted and segment == "D" and _is_inverted(simulation):
            skipped += 1
            continue
        if segment not in cache:
            cache[segment] = _allele_sequence(refdata, simulation, segment)
        reference = cache[segment]
        if reference is None:
            unresolved += 1
            continue
        if not 0 <= position < len(reference):
            out_of_range += 1
            continue
        checked += 1
        recorded = chr(germline[index]).upper()
        expected = chr(reference[position]).upper()
        if recorded != expected:
            out.append(Mismatch(index=index, germline_position=position, segment=segment,
                                recorded=recorded, reference=expected))
    return TemplatedCheck(mismatches=tuple(out), checked=checked, skipped_inverted=skipped,
                          unresolved_allele=unresolved, out_of_range=out_of_range)


def _is_inverted(simulation) -> bool:
    for region in simulation.regions():
        segment = region.segment() if callable(region.segment) else region.segment
        if segment == "D":
            inverted = getattr(simulation, "d_inverted", None)
            if inverted is not None:
                return bool(inverted() if callable(inverted) else inverted)
    return False


def untemplated_runs(simulation) -> tuple[UntemplatedRun, ...]:
    """Contiguous stretches with no germline coordinate, classified P or N.

    Case is the channel: lowercase is a P nucleotide, uppercase an N addition.
    A run reported as ``mixed`` spans both, which happens where a P abuts an N.

    The splitting itself lives in ``provenance``, which owns per-base origin.
    This adds only the P/N classification, so there is one definition of what
    a run is rather than two that could drift apart.
    """
    runs = []
    for run in provenance.untemplated_runs(simulation):
        text = run.bases.decode()
        lower = sum(1 for character in text if character.islower())
        kind = "P" if lower == len(text) else "N" if lower == 0 else "mixed"
        runs.append(UntemplatedRun(start=run.start, bases=text, kind=kind))
    return tuple(runs)


def split_by_case(run: UntemplatedRun) -> tuple[str, str]:
    """``(p_part, n_part)`` of a run, by case, preserving order."""
    return ("".join(c for c in run.bases if c.islower()),
            "".join(c for c in run.bases if c.isupper()))


def observable_substitutions(simulation) -> tuple[int, ...]:
    """Templated positions whose base differs from the germline base recorded there.

    Compared case-insensitively on purpose: case carries the P/N channel, not
    base identity, and a case-sensitive comparison silently reports every P
    nucleotide as a substitution.
    """
    bases, germline = simulation.bases(), simulation.germline_bases()
    return tuple(i for i in range(len(bases))
                 if simulation.germline_position(i) is not None
                 and chr(bases[i]).upper() != chr(germline[i]).upper())


def reverse_complement(sequence: str) -> str:
    return "".join(COMPLEMENT.get(base.upper(), "N") for base in reversed(sequence))


def np_fields_disagreement(record, runs) -> str | None:
    """How the recorded np fields differ from the untemplated runs, or ``None``.

    A zero-length NP region contributes no run, so the comparison is against
    the non-empty ones: records with ``np1_length == 0`` do occur, which is why
    this is not a positional comparison.

    Shared between the truth test and its mutant deliberately. A sensitivity
    test that reimplements the comparison it is testing proves only that two
    pieces of test code agree.
    """
    pairs = [(record["np1"], record["np1_length"]), (record["np2"], record["np2_length"])]
    expected = [(bases, length) for bases, length in pairs if length]
    found = [(run.bases, len(run.bases)) for run in runs]
    if found == expected:
        return None
    return (f"untemplated runs {[bases for bases, _ in found]} do not match "
            f"np1={record['np1']!r}/{record['np1_length']} "
            f"np2={record['np2']!r}/{record['np2_length']}")


def p_tip_disagreement(sequence: str, boundary: int, p_length: int) -> tuple[str, str] | None:
    """``(found, expected)`` if the P run at ``boundary`` is not a tip palindrome.

    A P addition of length ``n`` sitting at ``boundary`` must be the reverse
    complement of the ``n`` bases immediately before it -- the segment's own
    post-trim tip, copied back across the break. ``None`` means it is.

    Written as a function over the sequence rather than over a ``Simulation``
    so the same comparison serves the truth test and its mutant: the mutant
    hands it a sequence with the P run overwritten and must see it fire.

    Raises on input it cannot compare rather than returning ``None``. Returning
    ``None`` for both "agreed" and "there was nothing to compare" is the
    ambiguity ``templated_check`` exists to remove; re-introducing it in a
    helper written during that same fix would be the defect, not an exception
    to it. The caller decides whether a record has a P run to check.
    """
    if p_length <= 0:
        raise ValueError(
            f"p_length is {p_length}; there is no P run here to compare. Gate "
            "on the recorded p_*_length before calling.")
    if boundary - p_length < 0 or boundary + p_length > len(sequence):
        raise ValueError(
            f"a {p_length}-base P run at {boundary} needs {p_length} bases "
            f"either side, and the sequence is {len(sequence)} long")
    found = sequence[boundary:boundary + p_length].upper()
    expected = reverse_complement(sequence[boundary - p_length:boundary])
    return None if found == expected else (found, expected)
