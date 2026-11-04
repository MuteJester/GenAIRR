"""`count` means the same thing in every substitution model.

Three of the library's five substitution knobs draw the replacement base from
all four bases **including the base already at the site**, so about a quarter
of draws write what was already there. Those writes are correctly not counted
-- a sequence identical to its germline at that position is not a mutation --
which makes `count` an *attempt* count for those three and an *effect* count
for the other two. Measured on human_igh, 300 records, seed 11, before the fix:

| knob | realized/requested, count=1 | count=8 |
|---|---|---|
| `pcr_amplify` | 0.757 | 0.758 |
| `sequencing_errors` | 0.757 | 0.758 |
| `mutate(uniform)` | 0.780 | 0.767 |
| `ambiguous_base_calls` | 1.000 | 0.991 |
| `mutate(s5f)` | 1.000 | 1.000 |

So `pcr_amplify(count=8)` delivered about six PCR errors and nothing said so.
Three of the five corruption budgets in the library sat ~24% below their
nominal value, which is a quiet bias in every benchmark calibrated against
them.

The two that already agreed do so for different reasons, and neither is the
substitution draw being careful: `ambiguous_base_calls` *pins* its write to
`N` rather than sampling, and the S5F kernel's substitution row carries zero
mass on the base at the site, so its candidate list never contains it.

It is also wrong biologically. AID deamination and its repair yield
transitions and transversions; C->C is not a mutation that happens and fails
to be counted, it is not a mutation. A model whose "mutation" draw can return
the base already present is modelling nothing at that draw.

**How it surfaced**, which is the part worth keeping: a Rust test asked for
one PCR error and guarded its own premise by asserting the counter was
non-zero, so that tampering with it would not merely test the zero case. On
seed 17 the single draw was an identity write and `n_pcr_errors` came out 0.
A guard on a test's premise found a library-wide finding.

The expected value in every test here is the `count=` the experiment was
configured with -- a declared parameter, not another projection of the
engine's own output -- which is why the basis is `reference`.
"""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder

pytestmark = pytest.mark.behavior

RECORDS = 60
SEED = 11

#: Both a small and a large budget. The small one is where an attempt count
#: and an effect count differ by a whole event often enough to notice; the
#: large one is where the ~24% shortfall is unmistakable.
COUNTS = (1, 8)


def _realized(counter, count, configure):
    """`(realized, requested)` totals for one knob at one count.

    Returns totals rather than asserting, so the failure message can report
    the ratio -- the measurement is the finding, and a bare "3 != 4" would
    hide it.
    """
    exp = configure(experiment_builder.recombined(), count)
    realized = 0
    requested = 0
    short = []
    for outcome in exp.run(n=RECORDS, seed=SEED):
        record = outcome.airr_record(exp.refdata)
        got = record[counter]
        realized += got
        requested += count
        if got != count:
            short.append(got)
    return realized, requested, short


def _assert_every_record_realized(counter, count, configure):
    realized, requested, short = _realized(counter, count, configure)
    assert requested == RECORDS * count, (
        f"the fixture requested {requested} events over {RECORDS} records at "
        f"count={count}; the premise of this test is that it asked for "
        f"{RECORDS * count}")
    assert not short, (
        f"{counter} came out below the configured count={count} on "
        f"{len(short)} of {RECORDS} records (values {sorted(set(short))}); "
        f"realized/requested = {realized / requested:.3f}. A count the "
        f"library accepts and then under-delivers is a silent bias in every "
        f"budget calibrated against it")
    assert realized == requested


@verifies("assay.pcr_error",
          claim="n_pcr_errors equals the count the experiment asked for: "
                "every amplification draw changes the base at its site, so "
                "the configured budget is delivered rather than attempted",
          fields=("n_pcr_errors",),
          basis="reference")
@pytest.mark.parametrize("count", COUNTS)
def test_pcr_amplify_delivers_the_count_it_was_given(count):
    _assert_every_record_realized(
        "n_pcr_errors", count,
        lambda exp, n: exp.pcr_amplify(count=n))


@verifies("assay.quality_error",
          claim="n_quality_errors equals the count the experiment asked for. "
                "The comparison is case-insensitive on purpose: this pass "
                "writes its replacement lowercase, so a draw landing on the "
                "letter already there changes the byte without changing the "
                "base and must not count",
          fields=("n_quality_errors",),
          basis="reference")
@pytest.mark.parametrize("count", COUNTS)
def test_sequencing_errors_delivers_the_count_it_was_given(count):
    _assert_every_record_realized(
        "n_quality_errors", count,
        lambda exp, n: exp.sequencing_errors(count=n))


@verifies("maturation.hypermutation",
          claim="n_mutations equals the count the experiment asked for under "
                "the uniform kernel, as it already did under S5F: the "
                "replacement is drawn from the bases other than the one at "
                "the site, so count is a number of mutations rather than a "
                "number of attempts",
          fields=("n_mutations",),
          basis="reference")
@pytest.mark.parametrize("count", COUNTS)
def test_uniform_hypermutation_delivers_the_count_it_was_given(count):
    _assert_every_record_realized(
        "n_mutations", count,
        lambda exp, n: exp.mutate(model="uniform", count=n))


@pytest.mark.parametrize("count", COUNTS)
def test_the_s5f_kernel_already_delivered_its_count(count):
    """The anchor the other three are measured against.

    S5F was already at 1.000 before the fix, so it is what "count means
    effects" looked like when only one model did it. If this ever fails the
    comparison the module is built on has gone, and the uniform result above
    stops meaning the two models agree.
    """
    _assert_every_record_realized(
        "n_mutations", count,
        lambda exp, n: exp.mutate(model="s5f", count=n))
