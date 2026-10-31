"""Check the duplex repair plane against the molecule it claims to describe.

``duplex_repair`` is the only record of what a duplex repair actually did. The
exported ``sequence`` cannot carry it: a base resected from one strand and then
restored by template-directed fill leaves a molecule identical to the one a
blunt ligation produces, so a read is blind to the difference and this plane is
where it lives.

What this module reads is that plane and the sequence beside it. Two things
follow from that, and both are recorded in ``CIRCULARITY``: the comparison
between the product's top strand and the exported sequence is between two
projections of one state, while the comparisons across planes -- the
untemplated run against ``np1``, the physical counts against the columns that
produced them -- are between things written by different code paths.

The authored script is the independent side, and it is not available here: a
script is a GenAIRR object and an oracle may not import GenAIRR. The tests that
register this mechanism author the expected molecule and the expected counts by
hand and compare those; this module checks the plane's internal arithmetic.
"""
from __future__ import annotations

from dataclasses import dataclass

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

#: Every identity this module can report. A ``Disagreement`` may not name one
#: that is not here -- see the package docstring.
IDENTITIES = (
    "the product's top strand is the exported sequence",
    "each column's bottom base is the complement of its top base",
    "a column's kind accounts for its origin",
    "the physical counts are what the product shows",
    "the untemplated additions on the top strand are the np1 field",
    "the selected history is one the model offers",
)

#: What each ``kind`` in a product column means about that base's origin.
#: ``original_strand`` is a base that was already there and names where it came
#: from; ``n_addition`` is untemplated and names nothing; ``templated_synthesis``
#: was copied back from the opposite strand and names the template it came off.
ORIGIN_REQUIRED = {"original_strand": "original", "templated_synthesis": None,
                   "n_addition": None}


@dataclass(frozen=True)
class Disagreement:
    """One identity that does not hold, and where."""

    sequence_id: str
    identity: str
    recorded: object
    expected: object
    detail: str = ""

    def __post_init__(self) -> None:
        if self.identity not in IDENTITIES:
            raise ValueError(
                f"{self.identity!r} is not a declared identity. Add it to "
                "IDENTITIES and probe it for circularity, or reuse the key it "
                "belongs to and put the specifics in `detail`.")


def columns(record) -> list:
    """The product's columns, or ``[]`` when the record carries no duplex plane."""
    plane = record.get("duplex_repair")
    if not isinstance(plane, dict):
        return []
    return list((plane.get("product") or {}).get("columns") or ())


def top_strand(record) -> str:
    """The top strand of the repaired product, 5' to 3'."""
    return "".join((column.get("top") or {}).get("base", "-")
                   for column in columns(record))


def kind_counts(record) -> dict[str, int]:
    """How many bases of each kind the product carries, over both strands."""
    out: dict[str, int] = {}
    for column in columns(record):
        for side in ("top", "bottom"):
            base = column.get(side)
            if base:
                out[base["kind"]] = out.get(base["kind"], 0) + 1
    return out


def untemplated_run(record) -> str:
    """The top strand's untemplated bases, in order.

    One run in the fixtures measured, because a single join adds at one place;
    returned as a string rather than a list of runs because the np fields it is
    compared against are strings too.
    """
    return "".join(column["top"]["base"] for column in columns(record)
                   if (column.get("top") or {}).get("kind") == "n_addition")


def disagreements(record, *, offered_histories=()) -> tuple[Disagreement, ...]:
    """Every identity above that does not hold for this record.

    ``offered_histories`` is the set of history names the authored model
    declares, read by the caller from the model it built. Left empty, the
    selection identity is not checked -- there is nothing to check it against,
    and inventing an expectation would be worse than declaring the gap.
    """
    out: list[Disagreement] = []
    sequence_id = record.get("sequence_id", "?")
    cells = columns(record)
    if not cells:
        return ()
    plane = record["duplex_repair"]

    top = top_strand(record)
    if top.upper() != (record["sequence"] or "").upper():
        out.append(Disagreement(
            sequence_id, "the product's top strand is the exported sequence",
            top, record["sequence"]))

    for position, column in enumerate(cells):
        top_base, bottom_base = column.get("top"), column.get("bottom")
        if top_base and bottom_base:
            expected = COMPLEMENT.get(top_base["base"].upper())
            if bottom_base["base"].upper() != expected:
                out.append(Disagreement(
                    sequence_id,
                    "each column's bottom base is the complement of its top base",
                    bottom_base["base"], expected,
                    detail=f"column {position}, top {top_base['base']}"))
        for side in ("top", "bottom"):
            base = column.get(side)
            if not base:
                continue
            if base["kind"] not in ORIGIN_REQUIRED:
                out.append(Disagreement(
                    sequence_id, "a column's kind accounts for its origin",
                    base["kind"], sorted(ORIGIN_REQUIRED),
                    detail=f"column {position} {side}"))
                continue
            required = ORIGIN_REQUIRED[base["kind"]]
            if required and not base.get(required):
                out.append(Disagreement(
                    sequence_id, "a column's kind accounts for its origin",
                    base.get(required), f"{base['kind']} names its {required}",
                    detail=f"column {position} {side}"))
            if base["kind"] == "n_addition" and base.get("original"):
                out.append(Disagreement(
                    sequence_id, "a column's kind accounts for its origin",
                    base["original"], "an untemplated base has no origin",
                    detail=f"column {position} {side}"))

    counts = kind_counts(record)
    physical = (plane.get("counts") or {}).get("physical") or {}
    expected_physical = {
        # Every base the product carries, both strands.
        "retained_nucleotides": sum(counts.values()),
        "untemplated_additions": counts.get("n_addition", 0),
        "templated_synthesis": counts.get("templated_synthesis", 0),
        # A resected base is gone from the product, so the original count is
        # what survived plus what was cut.
        "original_nucleotides": (counts.get("original_strand", 0)
                                 + physical.get("resected_nucleotides", 0)),
    }
    for name, expected in expected_physical.items():
        if physical.get(name) != expected:
            out.append(Disagreement(
                sequence_id, "the physical counts are what the product shows",
                physical.get(name), expected, detail=name))

    run = untemplated_run(record)
    if (record.get("np1") or "").upper() != run.upper():
        out.append(Disagreement(
            sequence_id,
            "the untemplated additions on the top strand are the np1 field",
            record.get("np1"), run))

    if offered_histories:
        chosen = (plane.get("selection") or {}).get("history_id")
        if chosen not in set(offered_histories):
            out.append(Disagreement(
                sequence_id, "the selected history is one the model offers",
                chosen, sorted(offered_histories)))
    return tuple(out)


#: What was perturbed, and what this oracle said. Measured on the authored
#: three-allele fixture in tests/support/cases/duplex.py, which produces zero
#: disagreements unperturbed on each of the blunt, N-addition and
#: resect-and-fill scripts.
CIRCULARITY = {
    "the product's top strand is the exported sequence":
        "CIRCULAR IN CONTENT. Rewriting the first base of the exported "
        "sequence AND of the product's first top column together leaves this "
        "silent on 3 of 3 scripts; rewriting the sequence alone fires on 3 of "
        "3. Both sides are projections of one repaired duplex, so what this "
        "establishes is that the read a consumer gets is the strand the plane "
        "describes -- which is the thing a benchmark needs and is not nothing "
        "-- and not that either is the molecule the script built. The authored "
        "comparison in tests/biology/assay/test_duplex_truth.py is what "
        "establishes that, from a script written by hand outside the engine.",
    "each column's bottom base is the complement of its top base":
        "CIRCULAR IN BASE IDENTITY, independent in complementarity. Replacing "
        "a column's top and bottom with a different complementary pair is "
        "silent on 3 of 3 scripts; replacing the top alone fires on 3 of 3. "
        "Complementarity is not recorded anywhere in the plane -- nothing "
        "stores 'these two are a pair' -- so the engine cannot satisfy it "
        "except by having written a real duplex, which is why a single-strand "
        "perturbation is caught.",
    "a column's kind accounts for its origin":
        "NOT CIRCULAR, and it reads one field against another the engine "
        "writes separately. Relabelling an original_strand base as n_addition "
        "while leaving its origin in place fires on 3 of 3 scripts, and "
        "clearing the origin of an original_strand base fires on 3 of 3. The "
        "kind is what a consumer filters on to tell templated from "
        "untemplated content; the origin is what it would use to place that "
        "base in a germline. A record where they disagree would mislead "
        "either way round.",
    "the physical counts are what the product shows":
        "NOT CIRCULAR in the direction that matters. The counts are a summary "
        "the engine writes beside the columns, not a recount of them, so "
        "incrementing untemplated_additions alone fires on 3 of 3 scripts "
        "while the product stays as it was. The relation that makes this more "
        "than arithmetic is original_nucleotides == surviving original bases "
        "+ resected: on the resect-and-fill script 29 original bases survive "
        "of 30, and the restored base is counted as templated_synthesis "
        "rather than as original, which is the distinction the sequence "
        "cannot show.",
    "the untemplated additions on the top strand are the np1 field":
        "NOT CIRCULAR -- two planes written by different projections. The "
        "legacy np1 field comes from the AIRR junction projection and the run "
        "from the duplex product's own columns; rewriting np1 alone fires on "
        "the N-addition script, and both agree on 'GG' and 'ACGT' for the two "
        "authored additions measured. A duplex repair that added bases "
        "without reporting them in np1 would leave a consumer reading the "
        "junction fields blind to untemplated content that is in the read.",
    "the selected history is one the model offers":
        "NOT CIRCULAR, and only checked when the caller supplies the offered "
        "set -- which it reads from the model it authored, outside the engine. "
        "Declaring a one-row model and asserting the selection names that row "
        "holds on 3 of 3 scripts; naming a history the model does not offer "
        "fires. Left unchecked when the set is empty, because an expectation "
        "invented here would be worse than the declared gap.",
}
