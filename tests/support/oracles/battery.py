"""Every truth check this suite knows how to run, over one record.

Each biology test picks its oracle, builds that oracle's reference material by
hand, and chooses a fixture the oracle is valid on. That is right for a test
about one mechanism, and it leaves two things undone: no record is ever checked
by *all* of them, and the conditions under which a check is invalid live in
docstrings and in each test's choice of fixture rather than anywhere a program
can read.

This module closes both. It runs every applicable oracle over a record and
returns what they found, tagged with which oracle and which identity. It is
what makes a generated pipeline worth running: the sweep draws a stack nobody
wrote by hand, and the battery is the only thing that knows what should be
true of the result.

**Applicability is the load-bearing part.** An oracle applied outside its
stated conditions reports disagreements that are not defects -- the NP-field
comparison treats a polymerase insertion as an untemplated run and cries on
194 of 200 records; the observable-substitution identity counts a sequencing
miscall as hypermutation; the identity fraction becomes an inequality once an
indel adds alignment columns. A sweep that reported those would bury a real
finding in noise and get switched off, which is how a discovery campaign dies.

So each check declares, as a predicate over the pipeline that produced the
record, when it may be believed. Two consequences worth stating:

- ``Result.skipped`` reports what was not run and why, because a check that
  silently sits out is indistinguishable from a check that passed.
- Where a condition only weakens a comparison rather than voiding it, there
  are two checks instead of one: ``identity_exact`` is skipped under indels
  and ``identity_bound`` holds the inequality that survives them. Skipping is
  the last resort, not the first.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from . import coordinates, expression, junction, projection, reconstruction

#: Oracle modules deliberately NOT in the battery, with the reason. Checked by
#: the battery's own test, so a new oracle is either wired in or declared.
NOT_IN_BATTERY = {
    "duplex": "only applies to a record carrying a duplex_repair plane, which "
              "the authored duplex surface produces and no drawn pipeline "
              "does; tests/biology/assay/test_duplex_truth.py runs it",
}


@dataclass(frozen=True)
class Finding:
    """One disagreement, with enough context to act on it."""

    oracle: str
    identity: str
    sequence_id: str
    recorded: object
    expected: object
    detail: str = ""

    def __str__(self) -> str:
        where = f" [{self.detail}]" if self.detail else ""
        return (f"{self.oracle}: {self.identity}{where}\n"
                f"    recorded {self.recorded!r}\n"
                f"    expected {self.expected!r}")


@dataclass(frozen=True)
class Result:
    """What the battery found, and what it did not look at."""

    findings: tuple[Finding, ...]
    ran: tuple[str, ...]
    #: check name -> why it was not run on this record.
    skipped: dict[str, str] = field(default_factory=dict)
    #: check name -> how much it compared, so a vacuous pass is visible.
    compared: dict[str, int] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return bool(self.findings)

    @property
    def vacuous(self) -> tuple[str, ...]:
        """Checks that ran and compared nothing."""
        return tuple(name for name in self.ran if not self.compared.get(name))

    def report(self) -> str:
        lines = [f"{len(self.findings)} finding(s) from {len(self.ran)} check(s)"]
        lines.extend(str(finding) for finding in self.findings)
        if self.skipped:
            lines.append("  not run: " + "; ".join(
                f"{name} ({why})" for name, why in sorted(self.skipped.items())))
        if self.vacuous:
            lines.append(f"  compared nothing: {list(self.vacuous)}")
        return "\n".join(lines)


# ── applicability ────────────────────────────────────────────────────────────
#
# Each predicate answers one question: may this check be believed on a record
# the given pipeline produced? Written against the pipeline rather than the
# record wherever possible, because a pipeline is what a finding is reported
# against and what a reader reproduces.

def _writes_bases_other_than_hypermutation(pipeline) -> bool:
    """Any pass that changes a base without being somatic hypermutation."""
    return any(pipeline.has(name) for name in
               ("pcr_amplify", "sequencing_errors", "ambiguous_base_calls",
                "contaminate"))


def _may_insert(pipeline) -> bool:
    """Whether polymerase indels could add an untemplated base.

    Conservative on purpose: a drawn ``insertion_prob`` of 0 still counts,
    because the parameter may be a distribution and the draw is per record.
    """
    return pipeline.has("polymerase_indels")


def _changes_alignment_columns(pipeline) -> bool:
    return pipeline.has("polymerase_indels")


def _may_overwrite_the_molecule(pipeline) -> bool:
    return pipeline.has("contaminate")


@dataclass(frozen=True)
class Check:
    """One oracle call, with the conditions under which it may be believed."""

    name: str
    oracle: str
    #: (record, simulation, refdata, pipeline) -> (findings, how much compared)
    run: Callable
    #: (reason, predicate over the pipeline) -- all must hold to run.
    applies_when: tuple[tuple[str, Callable], ...] = ()

    def skip_reason(self, pipeline) -> str | None:
        for reason, predicate in self.applies_when:
            if not predicate(pipeline):
                return reason
        return None


# ── reference material ───────────────────────────────────────────────────────

def _allele_lengths(refdata, record) -> dict[str, int]:
    """Segment letter -> the length of the allele this record's call names."""
    from tests.support import alleles

    out = {}
    for letter in ("v", "d", "j"):
        name = record.get(f"{letter}_call")
        if not name:
            continue
        out[letter] = len(alleles.sequence(
            alleles.by_name(refdata, letter.upper(), name)))
    return out


def _reference_sequences(refdata, simulation) -> dict[str, str]:
    """Segment letter -> the assigned allele's bases, for the projection."""
    getters = {"v": refdata.v_allele, "d": refdata.d_allele, "j": refdata.j_allele}
    identifiers = {"v": simulation.v_allele_id(), "d": simulation.d_allele_id(),
                   "j": simulation.j_allele_id()}
    out = {}
    for segment, identifier in identifiers.items():
        if identifier is None:
            continue
        allele = getters[segment](identifier)
        raw = allele.seq() if callable(allele.seq) else allele.seq
        out[segment] = raw.decode().upper()
    return out


def _finding(oracle, disagreement) -> Finding:
    """One oracle's ``Disagreement`` as a battery finding."""
    return Finding(
        oracle, disagreement.identity,
        getattr(disagreement, "sequence_id", "?"),
        disagreement.recorded, disagreement.expected,
        detail=getattr(disagreement, "detail", "")
        or getattr(disagreement, "segment", ""))


# ── the checks ───────────────────────────────────────────────────────────────

def _coordinate_check(record, simulation, refdata, pipeline):
    lengths = _allele_lengths(refdata, record)
    found = coordinates.disagreements(record, lengths)
    return tuple(_finding("coordinates", one) for one in found), len(lengths)


def _projection_check(record, simulation, refdata, pipeline):
    inverted = {"d"} if record.get("d_inverted") else set()
    sequences = _reference_sequences(refdata, simulation)
    found = projection.disagreements(record, sequences, inverted=inverted)
    columns = projection.columns_compared(record, sequences, inverted=inverted)
    return (tuple(_finding("projection", one) for one in found),
            sum(columns.values()))


def _junction_check(record, simulation, refdata, pipeline):
    from tests.support import alleles

    if not record.get("v_call") or not record.get("j_call"):
        return (), 0
    v_allele = alleles.by_name(refdata, "V", record["v_call"])
    j_allele = alleles.by_name(refdata, "J", record["j_call"])
    found = junction.disagreements(
        record, simulation, reconstruction.segments_by_position(simulation),
        v_anchor=v_allele.anchor, j_anchor=j_allele.anchor,
        v_sequence=alleles.sequence(v_allele),
        j_sequence=alleles.sequence(j_allele))
    return tuple(_finding("junction", one) for one in found), 1


def _expression_check(record, simulation, refdata, pipeline):
    if record.get("junction_length") is None:
        # No junction, no resolvable frame. The junction oracle owns that
        # contract; inventing an expectation here would be worse.
        return (), 0
    found = expression.disagreements(record)
    return tuple(_finding("expression", one) for one in found), 3


def _assembly_check(record, simulation, refdata, pipeline):
    """The base-level walk: every recorded germline base against the reference."""
    inverted = {"D"} if record.get("d_inverted") else set()
    walk = reconstruction.templated_check(simulation, refdata, inverted=inverted)
    identifier = record.get("sequence_id", "?")
    findings = []
    if walk.unresolved_allele or walk.out_of_range:
        findings.append(Finding(
            "reconstruction",
            "every segment resolves and every coordinate is inside its allele",
            identifier,
            f"unresolved={walk.unresolved_allele} out_of_range={walk.out_of_range}",
            "0 of each"))
    for mismatch in walk.mismatches:
        findings.append(Finding(
            "reconstruction",
            "recorded germline base == reference base at the claimed coordinate",
            identifier, mismatch.recorded, mismatch.reference,
            detail=f"{mismatch.segment} at germline {mismatch.germline_position}"))
    return tuple(findings), walk.checked


def _np_fields_check(record, simulation, refdata, pipeline):
    """The exported np fields against the untemplated stretches of the read."""
    runs = reconstruction.untemplated_runs(simulation)
    message = reconstruction.np_fields_disagreement(record, runs)
    if message is None:
        return (), len(runs)
    return (Finding("reconstruction",
                    "the np fields are the untemplated runs' bases and lengths",
                    record.get("sequence_id", "?"),
                    f"np1={record['np1']!r} np2={record['np2']!r}",
                    message),), len(runs)


def _observable_mutation_check(record, simulation, refdata, pipeline):
    """``n_observable_mutations`` against the differences a caller can see.

    The identity that caught the shipped defect where the field was derived
    from a germline comparison and so reported library noise as hypermutation.
    """
    recorded = record.get("n_observable_mutations")
    if recorded is None:
        return (), 0
    visible = reconstruction.observable_substitutions(simulation)
    if recorded != len(visible):
        return (Finding("reconstruction",
                        "observable_substitutions == the engine's "
                        "n_observable_mutations",
                        record.get("sequence_id", "?"), recorded, len(visible),
                        detail=f"{len(visible)} visible differences"),
                ), max(len(visible), 1)
    return (), max(len(visible), 1)


def _identity_exact_check(record, simulation, refdata, pipeline):
    """Each segment's identity as matches over its templated positions."""
    counts = reconstruction.germline_match_counts(simulation)
    identifier = record.get("sequence_id", "?")
    findings = []
    checked = 0
    for segment, name in (("V", "v_identity"), ("D", "d_identity"),
                          ("J", "j_identity")):
        recorded = record.get(name)
        if segment not in counts:
            if recorded is not None:
                findings.append(Finding(
                    "reconstruction",
                    "a segment with no templated base reports a null identity",
                    identifier, recorded, None, detail=name))
            continue
        matching, total = counts[segment]
        checked += 1
        if recorded is None or abs(recorded - matching / total) > 1e-9:
            findings.append(Finding(
                "reconstruction",
                "each segment's identity is its matching fraction",
                identifier, recorded, matching / total,
                detail=f"{name}: {matching} of {total} match"))
    return tuple(findings), checked


def _identity_bound_check(record, simulation, refdata, pipeline):
    """The inequality that survives an indel.

    An inserted or deleted column has no matching germline base and the engine
    counts it in its own denominator, so the recorded identity sits at or below
    the substitution-only fraction. Always applicable, which is what keeps the
    indel case covered rather than skipped.
    """
    counts = reconstruction.germline_match_counts(simulation)
    identifier = record.get("sequence_id", "?")
    findings = []
    checked = 0
    for segment, name in (("V", "v_identity"), ("D", "d_identity"),
                          ("J", "j_identity")):
        recorded = record.get(name)
        if segment not in counts or recorded is None:
            continue
        matching, total = counts[segment]
        checked += 1
        if recorded > matching / total + 1e-9:
            findings.append(Finding(
                "reconstruction",
                "an indel lowers the identity it touches and never raises it",
                identifier, recorded, matching / total,
                detail=f"{name}: {matching} of {total} templated positions match"))
    return tuple(findings), checked


def _p_tip_check(record, simulation, refdata, pipeline):
    """Each recorded V 3' P run against the segment tip it mirrors.

    Only the V 3' run: the helper compares the bases after a boundary against
    the reverse complement of those before it, which is the arrangement at the
    V end. The other three P ends sit the other way round and are checked by
    the registered P-nucleotide tests, which know which side each one mirrors.
    """
    length = record.get("p_v_3_length") or 0
    boundary = record.get("v_sequence_end")
    sequence = record["sequence"]
    if not length or boundary is None:
        return (), 0
    if boundary - length < 0 or boundary + length > len(sequence):
        return (), 0
    result = reconstruction.p_tip_disagreement(sequence, boundary, length)
    if result is None:
        return (), 1
    found, expected = result
    return (Finding("reconstruction",
                    "a P run is the reverse complement of the segment's "
                    "post-trim tip", record.get("sequence_id", "?"),
                    found, expected, detail="p_v_3_length"),), 1


CHECKS = (
    Check("coordinates", "coordinates", _coordinate_check),
    Check("projection", "projection", _projection_check),
    Check("junction", "junction", _junction_check,
          applies_when=(("a contaminated read's bases are overwritten, so its "
                         "anchor codons are not the reference's and the bound "
                         "identities do not apply",
                         lambda pipeline: not _may_overwrite_the_molecule(pipeline)),)),
    Check("expression", "expression", _expression_check,
          applies_when=(("a contaminated read's coding span is uniform noise, "
                         "so a stop codon in it says nothing about the "
                         "productivity assessment",
                         lambda pipeline: not _may_overwrite_the_molecule(pipeline)),)),
    Check("assembly", "reconstruction", _assembly_check,
          applies_when=(("a contaminated read's bases are overwritten, so they "
                         "are deliberately not the germline's",
                         lambda pipeline: not _may_overwrite_the_molecule(pipeline)),)),
    Check("np_fields", "reconstruction", _np_fields_check,
          applies_when=(("a polymerase insertion carries no germline "
                         "coordinate either, so it forms an untemplated run of "
                         "its own and this comparison reports a record whose "
                         "np fields are correct -- 194 of 200 measured",
                         lambda pipeline: not _may_insert(pipeline)),
                        ("a contaminated read's untemplated stretches are "
                         "whatever the overwrite left",
                         lambda pipeline: not _may_overwrite_the_molecule(pipeline)),)),
    Check("observable_mutations", "reconstruction", _observable_mutation_check,
          applies_when=(("the comparison sees a difference, not which pass "
                         "wrote it, so a PCR misincorporation, a sequencing "
                         "miscall or an N mask lands in it exactly as a "
                         "mutation does; valid only where hypermutation is the "
                         "only pass that writes a base",
                         lambda pipeline:
                         not _writes_bases_other_than_hypermutation(pipeline)),
                        ("an indel shifts the positions the comparison walks",
                         lambda pipeline: not _changes_alignment_columns(pipeline)),)),
    Check("identity_exact", "reconstruction", _identity_exact_check,
          applies_when=(("an indel adds alignment columns the engine counts in "
                         "its denominator and the provenance has no position "
                         "for, so the equality becomes the inequality "
                         "identity_bound checks instead",
                         lambda pipeline: not _changes_alignment_columns(pipeline)),)),
    Check("identity_bound", "reconstruction", _identity_bound_check),
    Check("p_tip", "reconstruction", _p_tip_check,
          applies_when=(("an indel or an overwrite inside a P run breaks the "
                         "palindrome without the record being wrong",
                         lambda pipeline:
                         not _changes_alignment_columns(pipeline)
                         and not _may_overwrite_the_molecule(pipeline)),)),
)


def run(outcome, record, refdata, pipeline) -> Result:
    """Every applicable check over one record.

    ``pipeline`` decides which checks may be believed. A pipeline with no steps
    runs everything, which is what a hand-built fixture wants.
    """
    simulation = outcome.final_simulation()
    findings: list[Finding] = []
    ran: list[str] = []
    skipped: dict[str, str] = {}
    compared: dict[str, int] = {}
    for check in CHECKS:
        reason = check.skip_reason(pipeline)
        if reason is not None:
            skipped[check.name] = reason
            continue
        found, counted = check.run(record, simulation, refdata, pipeline)
        ran.append(check.name)
        compared[check.name] = counted
        findings.extend(found)
    return Result(tuple(findings), tuple(ran), skipped, compared)
