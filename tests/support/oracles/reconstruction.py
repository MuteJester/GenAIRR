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


#: Every identity this module can report, in the order the comparisons appear
#: below. Declared so ``CIRCULARITY`` can be checked to cover them. Unlike the
#: sibling oracles, the result types here carry no identity field -- a
#: ``Mismatch`` is a position, and two of the comparisons return a bare string
#: or tuple -- so there is no construction site to validate against this
#: tuple. :data:`REPORTS` is the substitute: it names, per public function,
#: which of these identities that function can report.
IDENTITIES = (
    "recorded germline base == reference base at the claimed coordinate",
    "recorded germline base == base of the allele the CALL names",
    "an inverted segment's byte == complement of the reference base",
    "observable_substitutions == the engine's n_observable_mutations",
    "a visible substitution is bucketed by the region or segment carrying it",
    "the np fields are the untemplated runs' bases and lengths",
    "a P run is the reverse complement of the segment's post-trim tip",
)

#: Every public function here, against the identities whose disagreements it
#: can report -- ``()`` where it reports none and only supplies material to
#: another comparison. Exhaustive by test: the meta test requires the keys to
#: be exactly this module's public functions, so a new comparison added here
#: fails until its author says which identity it reports, and an identity has
#: to be in :data:`IDENTITIES` and therefore probed.
REPORTS = {
    "segments_by_position": (),
    "templated_check": (
        "recorded germline base == reference base at the claimed coordinate",
        "recorded germline base == base of the allele the CALL names",
        "an inverted segment's byte == complement of the reference base",
    ),
    "untemplated_runs": (),
    "split_by_case": (),
    "observable_substitutions": (
        "observable_substitutions == the engine's n_observable_mutations",),
    "reverse_complement": (),
    "np_fields_disagreement": (
        "the np fields are the untemplated runs' bases and lengths",),
    "p_tip_disagreement": (
        "a P run is the reverse complement of the segment's post-trim tip",),
    "observable_substitutions_by_v_subregion": (
        "a visible substitution is bucketed by the region or segment carrying it",),
    "observable_substitutions_by_segment": (
        "a visible substitution is bucketed by the region or segment carrying it",),
}


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
    #: Positions in a segment assembled in reverse orientation. Counted rather
    #: than skipped: the IR stores the complement of the reference base at the
    #: original allele coordinate, so the comparison applies with that
    #: expectation. Measured 417 of 417 against the complement and 0 against
    #: the plain reference, on the assembly truth test's own fixture
    #: (invert_d(prob=1.0), n=25, seed=0).
    inverted_compared: int
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


def templated_check(simulation, refdata, *, inverted=(),
                    resolve=None) -> TemplatedCheck:
    """Compare every recorded germline base to the reference's base.

    This catches a wrong allele assembled with internally consistent
    provenance -- a failure the engine's own record cannot reveal, because its
    provenance would agree with itself.

    It does not catch a mis-trimmed segment. A coordinate is only visible here
    when it moves independently of the byte beside it, which a record-level
    mutant can do and a trimming defect cannot. Use the ``coordinates`` oracle
    for that.

    ``resolve`` chooses WHICH allele each segment's bases are compared
    against. It defaults to the allele the simulation's own assignment names,
    which is what establishes allele identity. Passing a resolver that looks
    the allele up by the record's exported ``*_call`` instead establishes
    something this walk could not otherwise see: that the name a benchmark
    reads is the name of the allele whose bases were actually laid down. The
    two differ -- the engine's re-projection comparison recomputes the call
    from the same assignment and so agrees with itself either way.

    ``inverted`` names the segments the simulation assembled in reverse
    orientation, which the caller reads from the record's ``d_inverted``. It
    has to be passed in: a ``Simulation`` exposes no orientation at all, which
    is why the previous ``skip_inverted`` could not work. It read
    ``getattr(simulation, "d_inverted", None)``, always found ``None``, and so
    skipped nothing while its docstring said inverted D was skipped by default.
    Nothing caught it because no registered test ran an inverting stack.

    An inverted segment is now compared rather than passed over. The IR stores
    the complement of the reference base at the original allele coordinate --
    ``Nucleotide::germline(complement(seq[pos]), pos)`` -- so the comparison
    applies with the expectation complemented. Measured on the assembly truth test's own fixture (n=25, seed=0):
    417 inverted positions compared, all agreeing with the complement, and
    0 agreeing with the plain reference.
    """
    germline = simulation.germline_bases()
    segments = segments_by_position(simulation)
    cache: dict[str, bytes | None] = {}
    out: list[Mismatch] = []
    checked = inverted_seen = unresolved = out_of_range = 0
    inverted = {s.upper() for s in inverted}
    for index in range(len(germline)):
        position = simulation.germline_position(index)
        if position is None:
            continue
        segment = segments.get(index, "?")
        if segment not in ("V", "D", "J"):
            continue
        if segment not in cache:
            cache[segment] = (resolve(segment) if resolve is not None
                              else _allele_sequence(refdata, simulation, segment))
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
        if segment in inverted:
            inverted_seen += 1
            expected = COMPLEMENT.get(expected, "N")
        if recorded != expected:
            out.append(Mismatch(index=index, germline_position=position, segment=segment,
                                recorded=recorded, reference=expected))
    return TemplatedCheck(mismatches=tuple(out), checked=checked,
                          inverted_compared=inverted_seen,
                          unresolved_allele=unresolved, out_of_range=out_of_range)


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

    **This counts mismatches, not hypermutation.** It cannot see which pass
    wrote a base, so a PCR misincorporation, a sequencing miscall or an N mask
    lands in the result exactly as a mutation does. As the expected value for
    ``n_observable_mutations`` it is therefore valid only where hypermutation
    is the only pass that writes a base; see this module's ``CIRCULARITY``
    entry for the measurement.
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

    **Valid only where the NP regions are the only untemplated bases.** A
    polymerase insertion carries no germline coordinate either, so it forms a
    run of its own and this comparison reports a disagreement on a record
    whose np fields are correct -- 194 of 200 on a pipeline with
    ``polymerase_indels(insertion_prob=1.0)``. Those are false alarms, not
    findings. The np fields are the interval between the segment bounds, which
    ``test_np_fields_are_the_interval_between_the_segments`` checks directly
    and which holds under insertions.

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


#: What a both-sides perturbation produces, per identity. See the package
#: docstring for why this is required and what it does not guarantee.
#:
#: Keys are checked against :data:`IDENTITIES` by the meta test, so an
#: identity this module reports cannot stay unprobed. The two that used to be
#: outside the mapping -- ``np_fields_disagreement`` and
#: ``p_tip_disagreement`` -- were in it only as a sentence saying they were
#: not, which is how an unprobed comparison survives a review. Both are
#: probed below.
#:
#: Still outside: the ``unresolved_allele`` and ``out_of_range`` counters,
#: which can each sink a verdict without a mismatch being reported. They are
#: checked by their own mutants in tests/biology/test_oracle_sensitivity.py,
#: which is weaker than a probe -- a mutant shows the check can fire, a probe
#: shows what it cannot see.
CIRCULARITY = {
    "recorded germline base == reference base at the claimed coordinate":
        "CIRCULAR. Shift the germline coordinate by one and re-read each byte "
        "from its own segment's reference: 0 mismatches of 366 compared. The "
        "engine lays byte and coordinate down from the same two numbers "
        "(Nucleotide::germline(slice[i], slice_start + i)), so a mis-trimmed "
        "segment moves both together. What this identity establishes is "
        "ALLELE IDENTITY -- a rotated reference breaks it -- and nothing about "
        "the coordinates. They are checked by the coordinates oracle against "
        "the allele's length, which this walk does not produce.",
    "observable_substitutions == the engine's n_observable_mutations":
        "CIRCULAR IN POSITION, BLIND TO ATTRIBUTION, and not a tautology "
        "about a single byte. Three probes, all measured over 30 records. "
        "(1) Read both the base and the germline byte one position over, "
        "together: 535 positions reported either way -- identical, because "
        "both bytes come from the same nucleotide record, so moving the read "
        "moves both. (2) Shift only the germline read: 8203 positions, so the "
        "comparison does fire once the two sides are decoupled. (3) Run a "
        "pipeline with no mutation pass at all and five sequencing errors: "
        "the engine reports 0 and this oracle reports 111. The oracle cannot "
        "see WHICH pass wrote a base, and the field counts only the "
        "hypermutation passes, so the identity holds exactly where "
        "hypermutation is the only pass that writes a base -- which is the "
        "condition the registered test runs under, and the reason the "
        "corruption-bearing cases are asserted against the error counters "
        "instead. Deriving the field from this comparison rather than from "
        "the mutation passes' own ledger is the defect that was shipped and "
        "then fixed: it reported library noise as hypermutation.",
    "recorded germline base == base of the allele the CALL names":
        "NOT CIRCULAR, and it sees what the default resolution cannot. "
        "Resolving each segment's reference by the record's exported "
        "`*_call` rather than by the simulation's own assignment: swapping "
        "v_call to a same-length allele is reported on 40 of 40 records, "
        "where the assignment-resolved walk reports 0 of 40 -- it is still "
        "comparing against the allele that was actually used. Swapping "
        "d_call fires 40 of 40. Swapping v_call to another allele of the "
        "same gene fires on 37 of 37 records the probe could perturb; the "
        "other 3 have a V gene with no second allele, so nothing was "
        "changed in them. 195 of 198 IGH V alleles share an exact length "
        "with another, so the length relation that was the only external "
        "check on v_call cannot separate them -- this can. "
        "Read `clean` rather than `mismatches`: a swap to a shorter allele "
        "puts positions past its end, which the walk counts as out-of-range "
        "(17 of 40 on the d_call probe) and not as a mismatch. Probing "
        "`mismatches` alone made a shorter-allele swap look invisible.",
    "the np fields are the untemplated runs' bases and lengths":
        "CIRCULAR IN CONTENT. Rewriting every untemplated base to A in the "
        "sequence AND in the np fields together leaves the comparison clean "
        "on 60 of 60 records: both sides are views of the same bytes, so the "
        "np bases the NP pass actually drew are not what this establishes. "
        "Rewriting the sequence alone fires on 59 of 60 -- the 60th record's "
        "only untemplated run already read 'AA', so the rewrite changed "
        "nothing there. What the identity does establish is that the exported "
        "np1/np2 strings and lengths are the untemplated stretches of the "
        "read, in order, which is what a consumer slicing the junction gets. "
        "Where the untemplated bases came from is the np composition "
        "question, and no test here answers it.",
    "a P run is the reverse complement of the segment's post-trim tip":
        "CIRCULAR IN POSITION, independent in palindromy. Overwriting the tip "
        "and the P run together with a matched pair (A*n and its reverse "
        "complement) leaves the comparison clean on 60 of 60 records, so this "
        "says nothing about the tip being the real post-trim germline tip -- "
        "that is the coordinates oracle's end-loss and trim arithmetic. "
        "Overwriting the P run alone is caught, measured by "
        "test_a_p_nucleotide_that_is_not_a_palindrome_is_caught. The "
        "palindrome relation itself is independent: it is nowhere in the "
        "record, so the engine cannot satisfy it except by having written the "
        "bases the P pass is defined to write.",
    "a visible substitution is bucketed by the region or segment carrying it":
        "NOT CIRCULAR in the bucketing, which is the step the base comparison "
        "above does not cover. Unperturbed, the buckets equal the record's "
        "own counters on 60 of 60 at one mutation per record and at twenty.\n\n"
        "        Subregions, at TWENTY mutations per record: shifting every "
        "declared interval by +15 stops matching on 58 of 60, emptying the "
        "table on 60 of 60, and rotating the LABELS across the same intervals "
        "on 48 of 60. At ONE mutation the same shift stops matching on only "
        "11 of the 49 records whose mutation lands in V, and the emptied "
        "table on 48 of 49 -- a single substitution usually sits far enough "
        "inside its region that a 15-base displacement does not move it out. "
        "The earlier version of this note carried the twenty-mutation numbers "
        "under a sentence saying the probe was run at one, which executing it "
        "as test_a_displaced_subregion_table_no_longer_matches_the_recorded_"
        "counters corrected. Reversing the order of the interval list does "
        "nothing at all (0 of 60) and is not a perturbation: the intervals "
        "are disjoint, so the same one still matches whichever end the search "
        "starts from. That measurement is recorded because it looked like a "
        "blind spot and is not one.\n\n"
        "        Segments: swapping the recorded V and D counts makes it "
        "disagree on 55 of 60, so a mutation routed to the wrong segment has "
        "somewhere to show up -- which the partition identity does not "
        "provide, since moving one mutation from D to J leaves the sum "
        "unchanged. The 5 of 60 that survive the swap are records whose V and "
        "D counts happen to be equal.\n\n"
        "        At twenty mutations per record the comparison is an "
        "inequality rather than an equality in both families, because an "
        "event whose base was overwritten or lost is still an event. That is "
        "the semantics, not a weakness, and it is why the observable "
        "counters exist beside the event ones.",
    "an inverted segment's byte == complement of the reference base":
        "NOT CIRCULAR in orientation. Declaring the segment forward when it is "
        "inverted gives mismatches on every inverted position (16 of 16 on the "
        "measured record), and declaring it inverted agrees on 417 of 417 "
        "across the assembly test's fixture. The "
        "coordinate circularity above still applies.",
}


V_SUBREGION_BUCKETS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "v_unannotated")


def observable_substitutions_by_v_subregion(simulation, v_subregions) -> dict[str, int]:
    """Recoverable V substitutions, bucketed by the allele's own IMGT table.

    ``v_subregions`` is the assigned V allele's ``subregions`` list, read by
    the caller from the reference and handed in -- the same arrangement the
    junction oracle uses, because this module may not import GenAIRR and so
    cannot resolve an allele itself. Each entry is ``(label, start, end)`` in
    reference coordinates, half-open.

    Restricted to V-segment positions. An earlier version of this bucketed
    every templated position against the V table, which put D and J
    coordinates into V subregions and reported 3 of 60 agreement where the
    truth was 60 of 60 -- a false failure that looked exactly like a defect
    in the counters.

    Same limit as :func:`observable_substitutions`: it sees a difference, not
    which pass wrote it, so it is valid as an expected value only where
    hypermutation is the only pass that writes a base.
    """
    segments = segments_by_position(simulation)
    bases, germline = simulation.bases(), simulation.germline_bases()
    counts = {bucket: 0 for bucket in V_SUBREGION_BUCKETS}
    for index in range(len(bases)):
        if segments.get(index) != "V":
            continue
        position = simulation.germline_position(index)
        if position is None:
            continue
        observed, reference = chr(bases[index]).upper(), chr(germline[index]).upper()
        if observed == reference:
            continue
        if observed not in "ACGT" or reference not in "ACGT":
            continue
        label = next((name.lower() for name, start, end in v_subregions
                      if start <= position < end), None)
        counts[label if label in counts else "v_unannotated"] += 1
    return counts


SEGMENT_BUCKETS = ("v", "d", "j", "np")


def observable_substitutions_by_segment(simulation) -> dict[str, int]:
    """Visible substitutions bucketed by the segment carrying them.

    The NP bucket collects both untemplated regions and anything the
    provenance does not label V, D or J, which is what the record's
    ``n_np_mutations`` does.

    Same limit as the other two in this family: it sees a difference, not
    which pass wrote it, so it is valid as an expected value only where
    hypermutation is the only pass that writes a base.
    """
    segments = segments_by_position(simulation)
    bases, germline = simulation.bases(), simulation.germline_bases()
    counts = {bucket: 0 for bucket in SEGMENT_BUCKETS}
    for index in range(len(bases)):
        observed, reference = chr(bases[index]).upper(), chr(germline[index]).upper()
        if observed == reference:
            continue
        if observed not in "ACGT" or reference not in "ACGT":
            continue
        segment = (segments.get(index) or "").upper()
        counts["v" if segment == "V" else
               "d" if segment == "D" else
               "j" if segment == "J" else "np"] += 1
    return counts
