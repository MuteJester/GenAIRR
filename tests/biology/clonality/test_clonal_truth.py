"""What a clonal record says about the cell it came from.

A lineage-reconstruction benchmark is scored against these fields. All four
clonality mechanisms had no registered test, and until this module the
legacy clonal path had no validation either: every one of the eight fields
could be edited and `validate_records` still returned ok, because the check
that resolves a record back to its tree node ran only in duplex mode.

The identities here are deliberately the ones that hold. Measurements that
turned up export-shape problems -- `lineage_parent_id` resolving to no
exported record, `lineage_node_id` unique only within a clone,
`lineage_affinity` inert at the default selection strength -- are recorded
in the mechanisms' truth gaps rather than pinned as if they were intended.
"""
from __future__ import annotations

import collections

import pytest

from tests.biology.catalog import verifies
from tests.support.builders import experiment as experiment_builder

pytestmark = pytest.mark.behavior

LINEAGE = dict(n_clones=4, max_generations=3, n_max=24, n_sample=5, rate=0.05)


def _lineage(**overrides):
    config = dict(LINEAGE)
    config.update(overrides)
    return experiment_builder.recombined().clonal_lineage(**config)


@verifies("clonality.founder",
          claim="every lineage field is resolved back to the tree node the "
                "record came from, so an edited one is caught. The addressing "
                "is the (tree, clone, node) origin rather than the record's "
                "own contents, because two cells can carry identical DNA and "
                "history and an edited namespace tag is invisible to every "
                "other check. sequence_id is validated the same way and is "
                "deliberately not claimed here: the catalog records it as no "
                "mechanism's truth, which is right for a plain record where it "
                "is an arbitrary tag",
          fields=("clone_id",),
          basis="provenance")
def test_an_edited_clonal_identity_is_caught():
    """The hole this closes: on the legacy path nothing validated these.

    Each field is edited in turn and the record set re-validated, so the test
    fails if any one of them stops being checked -- which is how the gap
    existed in the first place, the whole check being gated on duplex mode.
    """
    exp = _lineage()
    result = exp.run_records(seed=21)
    assert result.validate_records(exp.refdata).ok, (
        "the unedited result must validate, or the rest of this proves nothing")
    for field, replacement in (("clone_id", 9999),
                               ("sequence_id", "clone9999_node9999")):
        original = result[0][field]
        result[0][field] = replacement
        report = result.validate_records(exp.refdata)
        result[0][field] = original
        assert not report.ok, f"an edited {field} was not caught"
        kinds = {issue["kind"] for failure in report.failures
                 for issue in failure["issues"]}
        assert "LineageMetadataMismatch" in kinds, (
            f"{field}: caught, but not as a lineage mismatch: {sorted(kinds)}")


@verifies("clonality.sampling",
          claim="duplicate_count and lineage_abundance are the node's "
                "sampling draw count, and the per-clone sum recovers n_sample "
                "exactly. A single record's count is not a clone size",
          fields=("duplicate_count", "lineage_abundance"),
          basis="provenance")
@pytest.mark.parametrize("n_sample", [1, 5, 20])
def test_abundance_sums_to_the_sample_size_per_clone(n_sample):
    """The identity a consumer needs in order to use these fields at all.

    Both fields carry the same number -- the node's abundance, stamped twice --
    so neither is a read count or a tree-cell count. What is recoverable is
    the configured sample size, as the per-clone sum.
    """
    exp = _lineage(n_sample=n_sample)
    result = exp.run_records(seed=21)
    assert len(result), "no records were produced"
    by_clone = collections.defaultdict(int)
    for record in result:
        assert record["duplicate_count"] == record["lineage_abundance"], (
            "the two fields disagree, so they are not one number stamped twice")
        by_clone[record["clone_id"]] += record["lineage_abundance"]
    for clone_id, total in sorted(by_clone.items()):
        assert total == n_sample, (
            f"clone {clone_id} sums to {total} cells, not the {n_sample} drawn")


@verifies("clonality.founder",
          claim="clone membership is exactly recoverable from the per-base "
                "germline provenance on a substitution-only pipeline: every "
                "record of one clone carries the founder's germline bases and "
                "coordinates, and distinct clones differ. It is NOT "
                "recoverable once a length-changing pass is configured, which "
                "the mechanism's truth gap records",
          fields=("clone_id",),
          basis="provenance")
def test_clone_membership_is_recoverable_from_the_germline_provenance():
    """Why this matters: it is the only clone key that does not depend on
    trusting the namespace tag, so it is what an independent check of
    `clone_id` can be built on.

    Measured as exact on substitution-only pipelines and as fragmenting under
    polymerase indels, which is why the claim says substitution-only rather
    than claiming it generally.
    """
    exp = _lineage()
    result = exp.run_records(seed=21)
    keys = collections.defaultdict(set)
    # Paired by index rather than rebuilding the record from the outcome: the
    # lineage metadata is stamped by the Python layer, so a record rebuilt
    # from the outcome alone carries no clone_id at all.
    assert len(result) == len(result.outcomes), (
        "records and outcomes are not index-aligned, so they cannot be paired")
    for record, outcome in zip(result, result.outcomes):
        simulation = outcome.final_simulation()
        provenance = (bytes(simulation.germline_bases()),
                      tuple(simulation.germline_position(index)
                            for index in range(len(simulation.bases()))))
        keys[record["clone_id"]].add(provenance)
    assert len(keys) > 1, "only one clone was produced; distinctness is untested"
    for clone_id, distinct in sorted(keys.items()):
        assert len(distinct) == 1, (
            f"clone {clone_id} spans {len(distinct)} distinct germline "
            f"provenances, so membership is not recoverable from it")
    collapsed = {next(iter(distinct)) for distinct in keys.values()}
    assert len(collapsed) == len(keys), (
        "two clones share a germline provenance, so the key does not separate them")
