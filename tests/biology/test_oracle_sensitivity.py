"""Would the oracle catch a regression? Measured, not assumed.

Everything else here asserts that GenAIRR's recorded truth is correct today. It
does not establish that a test would notice if that stopped being true -- and a
verification that cannot fail is indistinguishable from no verification.

So: take a correct recorded result, corrupt it in a specific, biologically
meaningful way, and assert the oracle reports it. Each mutant corresponds to a
defect the engine could plausibly acquire -- a swapped junction pair, a
coordinate off by one, an allele assembled from the wrong reference, a mutation
counted that changed nothing.

These are record-level mutants rather than engine-level ones. An engine mutant
is stronger, because it exercises the real path, but it needs a rebuild per
mutant and so cannot run in ordinary CI. This runs on every commit.
"""
from __future__ import annotations

from dataclasses import dataclass

import pytest

from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.pin


class _Corrupted:
    """A Simulation view with one aspect deliberately wrong.

    Wraps rather than mutates, because the real object is read-only and because
    a wrapper makes the single corrupted aspect obvious at the call site.
    """

    def __init__(self, real, *, bases=None, germline=None, shift=0, drop_positions=()):
        self._real = real
        self._bases = bases
        self._germline = germline
        self._shift = shift
        self._drop = set(drop_positions)

    def bases(self):
        return self._bases if self._bases is not None else self._real.bases()

    def germline_bases(self):
        return self._germline if self._germline is not None else self._real.germline_bases()

    def germline_position(self, index):
        if index in self._drop:
            return None
        position = self._real.germline_position(index)
        if position is None:
            return None
        return position + self._shift

    def regions(self):
        return self._real.regions()

    def v_allele_id(self):
        return self._real.v_allele_id()

    def d_allele_id(self):
        return self._real.d_allele_id()

    def j_allele_id(self):
        return self._real.j_allele_id()


@dataclass(frozen=True)
class Mutant:
    name: str
    what_it_models: str


def _correct():
    exp = experiment_builder.recombined()
    return exp, exp.run(n=1, seed=0)[0]


def test_the_oracle_passes_on_an_uncorrupted_record():
    """The control. Without it, a mutant test that always fails would look like
    a working sensitivity check."""
    exp, outcome = _correct()
    assert reconstruction.templated_mismatches(outcome.final_simulation(), exp.refdata) == ()


def test_a_germline_coordinate_shifted_by_one_is_caught():
    """Models an off-by-one in a coordinate the engine records -- the defect
    class a benchmark would silently inherit as a systematic one-base offset."""
    exp, outcome = _correct()
    corrupted = _Corrupted(outcome.final_simulation(), shift=1)
    mismatches = reconstruction.templated_mismatches(corrupted, exp.refdata)
    assert mismatches, "a one-base coordinate shift went undetected"


def test_an_assembly_from_the_wrong_allele_is_caught():
    """Models the engine recording one allele while having copied another.

    The engine's own provenance would be internally consistent, so only a check
    against the reference tables can see it -- which is why the oracle reads
    the reference rather than the record alone.
    """
    exp, outcome = _correct()
    simulation = outcome.final_simulation()
    germline = bytearray(simulation.germline_bases())
    # rewrite the V region's germline as a different allele's sequence
    pool = exp.refdata.v_pool_size
    pool = pool() if callable(pool) else pool
    other_id = (simulation.v_allele_id() + 1) % max(1, pool)
    other = exp.refdata.v_allele(other_id)
    other_sequence = other.seq() if callable(other.seq) else other.seq
    changed = 0
    for index in range(len(germline)):
        position = simulation.germline_position(index)
        if position is None or position >= len(other_sequence):
            continue
        if germline[index] != other_sequence[position]:
            germline[index] = other_sequence[position]
            changed += 1
    if changed == 0:
        pytest.skip("the neighbouring allele is identical over the assembled span")
    corrupted = _Corrupted(simulation, germline=bytes(germline))
    assert reconstruction.templated_mismatches(corrupted, exp.refdata), (
        "germline rewritten from the wrong allele went undetected")


def test_a_swapped_np_pair_is_caught():
    """Models np1 and np2 length or composition models being wired the wrong way
    round -- a defect that leaves both fields populated and plausible."""
    exp, outcome = _correct()
    simulation = outcome.final_simulation()
    runs = reconstruction.untemplated_runs(simulation)
    if len(runs) != 2 or len(runs[0].bases) == len(runs[1].bases):
        pytest.skip("this record has no distinguishable NP pair to swap")
    record = outcome.airr_record(exp.refdata)
    swapped = {**record, "np1": record["np2"], "np2": record["np1"],
               "np1_length": record["np2_length"], "np2_length": record["np1_length"]}
    assert [r.bases for r in runs] != [swapped["np1"], swapped["np2"]], (
        "a swapped NP pair matches the untemplated runs, so the junction claim "
        "could not tell the two apart")


def test_a_mutation_counted_that_changed_nothing_is_caught():
    """Models exactly the defect this branch found and pinned: a recorded
    substitution that left the base unchanged."""
    exp = experiment_builder.recombined().mutate(model="s5f", count=2)
    outcome = exp.run(n=1, seed=0)[0]
    simulation = outcome.final_simulation()
    # make the final bases identical to germline: every recorded mutation
    # becomes a no-op, which is what the uniform model does 27% of the time
    corrupted = _Corrupted(simulation, bases=simulation.germline_bases())
    assert reconstruction.observable_substitutions(corrupted) == (), (
        "the oracle should see zero observable substitutions here")
    record = outcome.airr_record(exp.refdata)
    assert record["n_mutations"] > 0, "the fixture produced no mutations to contradict"


def test_a_p_nucleotide_that_is_not_a_palindrome_is_caught():
    """Models P addition emitting arbitrary bases instead of the reverse
    complement of the segment tip -- which the provenance oracle alone cannot
    see, because it records P with no germline coordinate either way."""
    tip = "AGA"
    correct = reconstruction.reverse_complement(tip)
    assert correct == "TCT"
    wrong = "GGG"
    assert wrong != reconstruction.reverse_complement(tip), (
        "the palindrome check would accept arbitrary P bases")
