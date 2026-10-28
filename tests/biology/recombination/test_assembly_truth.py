"""The assembled sequence must come from the alleles it says it came from."""
from __future__ import annotations

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder
from tests.support.oracles import reconstruction

pytestmark = pytest.mark.behavior

LOCI = ("human_igh", "human_igk", "human_tcrb")


@verifies("recombination.assembly",
          claim="every base carrying a germline coordinate matches the reference "
                "allele's base at that coordinate, so the recorded germline is the "
                "germline the record names",
          fields=("v_germline_start", "v_germline_end", "d_germline_start",
                  "d_germline_end", "j_germline_start", "j_germline_end"))
@pytest.mark.parametrize("locus", LOCI)
def test_every_templated_base_matches_the_reference_at_its_claimed_coordinate(locus):
    """Checked against the reference tables, not against the engine's own record.

    This is what catches a wrong allele assembled with internally consistent
    provenance, or a germline coordinate shifted by one -- neither of which the
    engine's own bookkeeping can reveal.
    """
    exp = experiment_builder.recombined(locus=locus)
    for outcome in exp.run(n=25, seed=0):
        simulation = outcome.final_simulation()
        mismatches = reconstruction.templated_mismatches(simulation, exp.refdata)
        assert mismatches == (), (
            f"{locus}: recorded germline disagrees with the reference at "
            f"{len(mismatches)} position(s); first: {mismatches[0]}")
