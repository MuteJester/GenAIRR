"""How much of the junction TdT actually added.

`np1` and `np2` are the stretch of molecule between the segment alignment
bounds. That identity holds exactly -- it is what the fields promise and
what `test_np_fields_are_the_interval_between_the_segments` pins -- but the
stretch is not the TdT contribution. It also holds:

- the **P nucleotides**, templated palindromic complements of the segment
  tips, which the documentation does say; and
- any **polymerase insertion** that landed between the two segments, which
  nothing said. Measured by running one seed with and without
  `polymerase_indels`, 24 of 200 records had an np length move by the
  insertion -- np1 27 -> 28 and np2 2 -> 3 among them.

So a consumer reading `np1_length` as the TdT amount was wrong by the
number of insertions that fell in the junction, with nothing in the record
to say by how much. These six fields partition each np interval by what
wrote each base, which is the only thing that answers the mechanism's own
claim.

No new biology: `P_NUC`, `N_NUC` and `INDEL_INSERTED` are already recorded
per nucleotide. This is a counting projection over provenance the engine
had all along.

The oracle here reads those flags through the masks exposed on the
simulation, so it derives the partition independently of the fields it
checks -- see `test_base_provenance_exposure.py` for why letter case could
not serve that purpose.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder

pytestmark = pytest.mark.behavior

RECORDS = 80
SEED = 9


def _interval(record, which):
    """The `[start, end)` pool span of an np region, or `None`.

    Taken from the same alignment bounds the builder uses, so the oracle
    reads the interval the fields are about rather than re-deriving a
    different one.
    """
    if which == 1:
        start, end = record["v_sequence_end"], record["d_sequence_start"]
        if end is None:
            end = record["j_sequence_start"]
    else:
        start, end = record["d_sequence_end"], record["j_sequence_start"]
    if start is None or end is None or end <= start:
        return None
    return int(start), int(end)


def _origin_counts(simulation, span):
    """`(tdt, p, inserted)` over a span, counted from the engine's flags."""
    start, end = span
    p_mask = simulation.p_nucleotide_mask()
    n_mask = simulation.tdt_addition_mask()
    inserted_mask = simulation.indel_inserted_mask()
    tdt = sum(1 for i in range(start, end) if n_mask[i])
    p = sum(1 for i in range(start, end) if p_mask[i])
    inserted = sum(1 for i in range(start, end) if inserted_mask[i])
    return tdt, p, inserted


@verifies("junction.n_addition",
          claim="each np interval is partitioned by what wrote each base: "
                "np1_n_tdt counts the TdT additions, np1_n_p the palindromic "
                "additions and np1_n_inserted the polymerase insertions that "
                "landed between the segments, and the three sum to "
                "np1_length. Counted from the engine's per-base provenance "
                "flags, which is why the partition still holds once an indel "
                "pass has moved the interval",
          fields=("np1_n_tdt", "np1_n_p", "np1_n_inserted",
                  "np2_n_tdt", "np2_n_p", "np2_n_inserted"),
          basis="provenance")
@pytest.mark.parametrize("with_indels", [False, True], ids=["plain", "indels"])
def test_the_np_partition_sums_to_the_np_length(with_indels):
    """The partition identity, on both sides of an indel pass.

    Both sides matter because the indel pass is what made `np1_length` stop
    meaning "TdT added this many". If the identity held only without it, the
    fields would be reporting a quantity that is correct exactly where it
    was never in doubt.
    """
    exp = experiment_builder.with_p_nucleotides_at_every_end()
    if with_indels:
        exp = exp.polymerase_indels(count=2)

    checked = 0
    for outcome in exp.run(n=RECORDS, seed=SEED):
        record = outcome.airr_record(exp.refdata)
        # Each field named literally rather than built with an f-string: the
        # catalogue checks that a registration reads the fields it claims,
        # and a field reached through an interpolated key is invisible to it.
        for region, length, tdt, p, inserted in (
            ("np1", record["np1_length"], record["np1_n_tdt"],
             record["np1_n_p"], record["np1_n_inserted"]),
            ("np2", record["np2_length"], record["np2_n_tdt"],
             record["np2_n_p"], record["np2_n_inserted"]),
        ):
            assert tdt + p + inserted == length, (
                f"{region}: {tdt} TdT + {p} P + {inserted} inserted "
                f"!= {length} bases of {region}_length")
            assert min(tdt, p, inserted) >= 0, (
                f"{region} carries a negative count: "
                f"tdt={tdt} p={p} inserted={inserted}")
            checked += 1
    assert checked, "no record was examined"


@pytest.mark.parametrize("with_indels", [False, True], ids=["plain", "indels"])
def test_the_np_partition_matches_the_recorded_provenance(with_indels):
    """The counts against the flags, position by position.

    The sum identity above would be satisfied by any three numbers adding
    to the length, including a projection that put every base in the wrong
    bucket. This is the check that each count is the bucket it claims.
    """
    exp = experiment_builder.with_p_nucleotides_at_every_end()
    if with_indels:
        exp = exp.polymerase_indels(count=2)

    compared = 0
    for outcome in exp.run(n=RECORDS, seed=SEED):
        record = outcome.airr_record(exp.refdata)
        simulation = outcome.final_simulation()
        for which in (1, 2):
            span = _interval(record, which)
            if span is None:
                assert record[f"np{which}_n_tdt"] == 0
                assert record[f"np{which}_n_p"] == 0
                assert record[f"np{which}_n_inserted"] == 0
                continue
            tdt, p, inserted = _origin_counts(simulation, span)
            assert (record[f"np{which}_n_tdt"],
                    record[f"np{which}_n_p"],
                    record[f"np{which}_n_inserted"]) == (tdt, p, inserted), (
                f"np{which} partition disagrees with the per-base flags over "
                f"{span}: record says "
                f"{record[f'np{which}_n_tdt']}/{record[f'np{which}_n_p']}/"
                f"{record[f'np{which}_n_inserted']}, flags say "
                f"{tdt}/{p}/{inserted}")
            compared += tdt + p + inserted
    assert compared, "no np base was compared against its flags"


def test_the_tdt_count_is_not_merely_the_np_length_again():
    """The fields have to *differ* from `np1_length` to be worth exporting.

    A projection that set `n_tdt = np_length` and the other two to zero
    would satisfy the sum identity and agree with the flags on any record
    whose junction happens to be pure TdT. This pins that the distinction
    is real: with P addition forced at every coding end and an indel pass
    running, some record's TdT count must come out below its np length.
    """
    exp = (experiment_builder.with_p_nucleotides_at_every_end()
           .polymerase_indels(count=2))

    below = 0
    with_p = 0
    with_inserted = 0
    for outcome in exp.run(n=RECORDS, seed=SEED):
        record = outcome.airr_record(exp.refdata)
        for which in (1, 2):
            if record[f"np{which}_length"] == 0:
                continue
            if record[f"np{which}_n_tdt"] < record[f"np{which}_length"]:
                below += 1
            if record[f"np{which}_n_p"]:
                with_p += 1
            if record[f"np{which}_n_inserted"]:
                with_inserted += 1
    assert below, (
        "every np interval was pure TdT, so this configuration cannot tell "
        "the partition from a copy of np_length")
    assert with_p, "no P nucleotide landed in an np interval"
    assert with_inserted, (
        "no polymerase insertion landed in an np interval, which is the case "
        "the mechanism's truth gap was specifically about")
