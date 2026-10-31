"""Every simulation pass the discovery sweep may draw, and what it may draw.

The suite's truth claims rest on near-single-pass pipelines. In tests/biology
``.mutate`` appears twenty-five times, ``end_loss`` seventeen, indels fourteen
-- almost always alone. Users stack five passes, and both of the hard defects
this project has found were interactions: the junction end wrong by the indel
length when an indel lands inside the J anchor codon, and ``end_loss_5_length``
counting only the first truncation. Nobody had tested those pairs.

This module is the input to a generator that builds stacks nobody wrote by
hand. Each pass declares four things, and each is load-bearing:

- ``stage``: where it sits in a legal pipeline. A drawn stack is sorted by
  stage, so the generator produces programs the compiler accepts instead of
  mostly ValueError. The order here is the DSL's own, not a guess --
  ``test_the_declared_stage_order_is_the_order_the_compiler_accepts`` holds it
  to that.
- ``parameters``: a Hypothesis strategy per keyword, drawing values inside the
  range the DSL documents. Deliberately narrow where width buys nothing: a
  drawn ``count`` of 400 mutations on a 300-base molecule tests the clamp, not
  the biology.
- ``requires``: preconditions, as predicates over the locus. ``mutate`` is
  B-cell only because the engine refuses to simulate hypermutation it does not
  claim; ``invert_d`` needs a D segment or it is a no-op that would pass by
  having nothing to do.
- ``mechanisms``: the catalog ids this pass exercises. It ties a campaign back
  to the mechanism table: a run can say which mechanisms it actually reached
  rather than how many examples it drew.

**Out of scope for v1**, and the reason is the same for all three: they change
the shape of the output rather than the content of a record, so the single
record battery does not apply to them unchanged. ``clonal_lineage`` /
``clonal_repertoire`` / ``expand_clones`` fork into families; ``recombine``
with ``duplex_repair`` replaces junction authoring wholesale and has its own
oracle; ``with_genotype`` changes which alleles exist, which is a different
question from whether a record's metadata is true.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from hypothesis import strategies as st

#: Loci the sweep draws from, with whether each has a D segment and whether it
#: is a B-cell locus. TCRB is in because a TCR locus is where a pass that
#: assumes SHM shows up, and the two IG loci differ in having a D at all.
LOCI = {
    "human_igh": {"has_d": True, "b_cell": True},
    "human_igk": {"has_d": False, "b_cell": True},
    "human_tcrb": {"has_d": True, "b_cell": False},
}

#: Pipeline stages in the order a legal pipeline applies them. `recombine` is
#: implicit -- every drawn pipeline starts there -- and `paired_end` is last
#: because the compiler's pre-pass pushes it to the end of the plan whatever
#: the call order, so generating it anywhere else would be generating a program
#: that does not mean what it reads like.
#:
#: `post_recombine` is a stage of its own for one pass: the DSL requires
#: `trim()` *immediately* after `recombine()`, with nothing between, and
#: refuses it even after an `invert_d()`. Found by the first run of
#: test_every_drawn_pipeline_builds_and_runs, which is what that test is for.
#: `process_coding_ends()` does NOT share the rule -- it only has to precede
#: mutation and corruption -- so it sits in `recombination` and is kept apart
#: from `trim` by `excludes` instead.
STAGES = ("post_recombine", "recombination", "mutation", "corruption",
          "readout", "constraint")


@dataclass(frozen=True)
class Pass:
    """One DSL method the sweep may put in a pipeline."""

    name: str
    stage: str
    #: keyword -> strategy for its value.
    parameters: dict[str, st.SearchStrategy]
    #: Catalog mechanism ids this pass exercises.
    mechanisms: tuple[str, ...] = ()
    #: Preconditions on the locus, as (reason, predicate).
    requires: tuple[tuple[str, Callable[[dict], bool]], ...] = ()
    #: Passes that may not appear beside this one, with the reason.
    excludes: dict[str, str] = field(default_factory=dict)
    notes: str = ""

    def __post_init__(self) -> None:
        if self.stage not in STAGES:
            raise ValueError(
                f"{self.name}: stage {self.stage!r} is not one of {STAGES}")
        if not self.parameters and not self.notes:
            raise ValueError(
                f"{self.name}: a pass with no parameter strategy can only ever "
                "be drawn as its default call, so say in `notes` why that is "
                "the whole of it")

    def applies_to(self, locus: str) -> bool:
        """Whether this pass can run on ``locus`` at all."""
        facts = LOCI[locus]
        return all(predicate(facts) for _reason, predicate in self.requires)

    def why_not(self, locus: str) -> tuple[str, ...]:
        """The reasons it cannot, for a message that explains itself."""
        facts = LOCI[locus]
        return tuple(reason for reason, predicate in self.requires
                     if not predicate(facts))


def _counts(low: int, high: int) -> st.SearchStrategy:
    """A count the DSL accepts: a fixed int, a uniform range, or weights.

    All three shapes because they are three code paths in the normalizer, and
    a sweep that only ever drew integers would never reach the other two.
    """
    fixed = st.integers(min_value=low, max_value=high)
    ranges = st.tuples(st.integers(min_value=low, max_value=high),
                       st.integers(min_value=low, max_value=high)).map(
        lambda pair: (min(pair), max(pair)))
    weighted = st.lists(
        st.tuples(st.integers(min_value=low, max_value=high),
                  st.floats(min_value=0.1, max_value=10, allow_nan=False,
                            allow_infinity=False)),
        min_size=1, max_size=3)
    return st.one_of(fixed, ranges, weighted)


def _probabilities() -> st.SearchStrategy:
    """A probability, including both endpoints.

    Zero is in deliberately: a pass at prob=0 must be a no-op that still
    records its coin flip, and that is a claim worth drawing.
    """
    return st.floats(min_value=0.0, max_value=1.0, allow_nan=False,
                     allow_infinity=False)


#: Trim and coding-end distributions are weight lists over a small range.
#: Negative values are the DSL's convention for "extend rather than trim" at a
#: coding end, which is how P nucleotides are reached.
def _trim_distribution() -> st.SearchStrategy:
    return st.lists(
        st.tuples(st.integers(min_value=0, max_value=8),
                  st.floats(min_value=0.1, max_value=5, allow_nan=False,
                            allow_infinity=False)),
        min_size=1, max_size=4)


def _coding_end_distribution() -> st.SearchStrategy:
    return st.lists(
        st.tuples(st.integers(min_value=-4, max_value=4),
                  st.floats(min_value=0.1, max_value=5, allow_nan=False,
                            allow_infinity=False)),
        min_size=1, max_size=4)


#: Mechanisms every drawn pipeline exercises through ``recombine()`` itself,
#: which is implicit and so has no ``Pass`` entry: a record cannot exist
#: without an assembly, a junction, coordinates and a projection. Declared so
#: the sweep's reach is counted honestly -- without this, the blind-spot list
#: would name mechanisms every single example covers.
ALWAYS_EXERCISED = (
    "recombination.assembly",
    "junction.assembly",
    "junction.n_addition",
    "projection.segment_coordinates",
    "projection.cigar",
    "projection.germline_alignment",
    "projection.translation",
    "projection.alignment_scores",
    "germline.allele_choice",
)


PASSES = (
    # ── recombination-stage modifiers ───────────────────────────────────────
    Pass(
        name="trim",
        stage="post_recombine",
        parameters={"v_3": _trim_distribution(), "j_5": _trim_distribution()},
        mechanisms=("recombination.trimming",),
        excludes={"process_coding_ends":
                  "the DSL refuses both: a signed coding-end outcome already "
                  "specifies reference loss, so trimming on top would be two "
                  "answers to one question"},
        notes="exonuclease chew-back at the coding ends, which is what moves "
              "the germline coordinates and can eat an anchor codon",
    ),
    Pass(
        name="process_coding_ends",
        stage="recombination",
        parameters={"v_3": _coding_end_distribution(),
                    "j_5": _coding_end_distribution(),
                    "reference_ends_are_coding_ends": st.just(True)},
        mechanisms=("junction.p_nucleotide", "junction.coding_end_processing"),
        excludes={"trim": "see trim's own entry: mutually exclusive"},
        notes="a negative draw extends rather than trims, which is how a P "
              "nucleotide is reached at all",
    ),
    Pass(
        name="invert_d",
        stage="recombination",
        parameters={"prob": _probabilities()},
        mechanisms=("recombination.d_inversion",),
        requires=(("the locus has no D segment, so inverting one is a no-op "
                   "that would pass by having nothing to do",
                   lambda facts: facts["has_d"]),),
    ),
    Pass(
        name="receptor_revision",
        stage="recombination",
        parameters={"prob": _probabilities(),
                    "same_haplotype": st.booleans()},
        mechanisms=("expression.receptor_revision", "germline.allele_choice"),
        requires=(("receptor revision is only valid for VDJ chains; the DSL "
                   "refuses it on a VJ locus, which the first run of "
                   "test_each_pass_runs_on_every_locus_it_claims found",
                   lambda facts: facts["has_d"]),),
    ),
    # ── mutation ────────────────────────────────────────────────────────────
    Pass(
        name="mutate",
        stage="mutation",
        parameters={"model": st.sampled_from(["s5f", "uniform"]),
                    "count": _counts(0, 25)},
        mechanisms=("maturation.hypermutation",),
        requires=(("somatic hypermutation does not occur in TCR sequences and "
                   "the engine refuses to simulate it",
                   lambda facts: facts["b_cell"]),),
    ),
    # ── corruption / library prep ───────────────────────────────────────────
    Pass(
        name="pcr_amplify",
        stage="corruption",
        parameters={"count": _counts(0, 10)},
        mechanisms=("assay.pcr_error",),
    ),
    Pass(
        name="sequencing_errors",
        stage="corruption",
        parameters={"count": _counts(0, 10)},
        mechanisms=("assay.quality_error",),
    ),
    Pass(
        name="ambiguous_base_calls",
        stage="corruption",
        parameters={"count": _counts(0, 8)},
        mechanisms=("assay.ambiguous_calls",),
    ),
    Pass(
        name="polymerase_indels",
        stage="corruption",
        parameters={"count": _counts(0, 4),
                    "insertion_prob": _probabilities(),
                    # Bounded deliberately. Under a productive constraint an
                    # indel count can be unsatisfiable -- four deletions
                    # cannot leave the coding frame intact -- and the
                    # constrained sampler is documented to search until its
                    # budget is spent and then raise StrictSamplingError. At
                    # the 20,000,000-node default that search takes six
                    # minutes, which a campaign cannot spend per example.
                    #
                    # The limitation this buys: a pipeline the campaign
                    # records as refused might succeed for a user on the
                    # default budget. The budget bounds how long the sampler
                    # looks, not what it finds, so a pipeline that produces
                    # records here produces the same records there.
                    "exact_node_budget": st.just(200_000)},
        mechanisms=("assay.polymerase_indels",),
        notes="the pass that breaks the one-column-one-base assumption every "
              "coordinate identity rests on",
    ),
    Pass(
        name="end_loss_5prime",
        stage="corruption",
        parameters={"length": _counts(0, 20)},
        mechanisms=("assay.end_loss",),
    ),
    Pass(
        name="end_loss_3prime",
        stage="corruption",
        parameters={"length": _counts(0, 20)},
        mechanisms=("assay.end_loss",),
    ),
    Pass(
        name="random_strand_orientation",
        stage="corruption",
        parameters={"prob": _probabilities()},
        mechanisms=("assay.orientation",),
    ),
    Pass(
        name="contaminate",
        stage="corruption",
        parameters={"prob": _probabilities()},
        mechanisms=("assay.contamination",),
        notes="overwrites the molecule with uniform bases, so a contaminated "
              "record's germline claims are deliberately false and the "
              "battery must not read them as findings",
    ),
    # ── readout ─────────────────────────────────────────────────────────────
    Pass(
        name="paired_end",
        stage="readout",
        parameters={"r1_length": st.integers(min_value=40, max_value=150),
                    "r2_length": st.integers(min_value=40, max_value=150),
                    "insert_size": st.integers(min_value=150, max_value=260)},
        mechanisms=("assay.paired_end",),
        notes="at most once per experiment, and the compiler pushes it to the "
              "end of the plan whatever the call order",
    ),
    # ── constraints ─────────────────────────────────────────────────────────
    Pass(
        name="productive_only",
        stage="constraint",
        parameters={},
        mechanisms=("expression.productivity",),
        notes="takes no arguments: it either constrains the repertoire or it "
              "is absent, so the default call is the whole of it. Worth "
              "drawing because it constrains every slot before proposing "
              "rather than rejection-sampling, so it changes what the other "
              "passes see",
    ),
)


def by_name() -> dict[str, Pass]:
    return {one.name: one for one in PASSES}


def for_locus(locus: str) -> tuple[Pass, ...]:
    """Every pass that can legally run on ``locus``."""
    return tuple(one for one in PASSES if one.applies_to(locus))
