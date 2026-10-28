"""Check the recorded sequence against the reference tables it claims to come from.

The provenance oracle next door derives the junction's untemplated content from
the engine's own record. That is useful but it cannot contradict the engine: if
provenance were corrupt, the oracle would agree with it. This module closes
that gap by reading the *reference* -- the allele sequences the simulation says
it drew from -- and asserting the recorded germline matches.

Taking the reference tables is not circular. Calling the production junction
builder or region projector would be; reading the allele a record names and
looking up the base at the coordinate it claims is the independent step.

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


def templated_mismatches(simulation, refdata, *, skip_inverted: bool = True
                         ) -> tuple[Mismatch, ...]:
    """Positions where the recorded germline base is not the reference's base.

    This is the check that catches a wrong allele assembled with internally
    consistent provenance, or a germline coordinate shifted by one -- failures
    the engine's own record cannot reveal.

    Inverted D is skipped by default: its coordinates are forward-reference
    while the sequence carries the reverse complement, so a direct lookup does
    not apply and comparing it needs its own claim.
    """
    germline = simulation.germline_bases()
    segments = segments_by_position(simulation)
    cache: dict[str, bytes | None] = {}
    out: list[Mismatch] = []
    for index in range(len(germline)):
        position = simulation.germline_position(index)
        if position is None:
            continue
        segment = segments.get(index, "?")
        if segment not in ("V", "D", "J"):
            continue
        if skip_inverted and segment == "D" and _is_inverted(simulation):
            continue
        if segment not in cache:
            cache[segment] = _allele_sequence(refdata, simulation, segment)
        reference = cache[segment]
        if reference is None or not 0 <= position < len(reference):
            continue
        recorded = chr(germline[index]).upper()
        expected = chr(reference[position]).upper()
        if recorded != expected:
            out.append(Mismatch(index=index, germline_position=position, segment=segment,
                                recorded=recorded, reference=expected))
    return tuple(out)


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
    """
    bases = simulation.bases()
    runs: list[UntemplatedRun] = []
    current: list[tuple[int, str]] = []

    def flush() -> None:
        if not current:
            return
        text = "".join(ch for _, ch in current)
        lower = sum(1 for ch in text if ch.islower())
        kind = "P" if lower == len(text) else "N" if lower == 0 else "mixed"
        runs.append(UntemplatedRun(start=current[0][0], bases=text, kind=kind))

    for index in range(len(bases)):
        if simulation.germline_position(index) is None:
            char = chr(bases[index])
            if current and index != current[-1][0] + 1:
                flush()
                current = []
            current.append((index, char))
        elif current:
            flush()
            current = []
    flush()
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
