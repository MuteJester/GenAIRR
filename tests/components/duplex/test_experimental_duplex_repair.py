import json

import pytest

from GenAIRR.experimental import DuplexRepair

pytestmark = pytest.mark.behavior


def opened(left="AAG", right="CAA", offset=0, **kwargs):
    return DuplexRepair(left, right, **kwargs).open("left", offset=offset).open("right", offset=offset)


def test_hairpin_products_preserve_two_strand_ancestries_and_exact_script_replay():
    script = opened(offset=1).synapse(right_start=3).ligate()
    result = script.run()
    assert result["sequence"] == "AAGCAA"
    assert result["complementary_sequence"] == "TTGCTT"
    assert result["columns"][2]["top"]["original"]["side"] == "left"
    assert result["columns"][2]["bottom"]["original"]["side"] == "right"
    assert result["columns"][3]["top"]["original"]["strand"] == "bottom"
    assert len({b["id"] for col in result["columns"] for b in (col["top"], col["bottom"])}) == 12
    assert DuplexRepair.from_dict(json.loads(json.dumps(script.to_dict()))).run() == result


def test_resection_is_recorded_even_after_template_restores_same_sequence():
    script = (opened("ACG", "TTA").resect("left", "top", count=1)
              .synapse(right_start=3).fill("top").ligate())
    result = script.run()
    assert result["sequence"] == "ACGTTA"
    assert len(result["removed"]) == 1
    assert result["columns"][2]["top"]["kind"] == "templated_synthesis"
    assert result["columns"][2]["top"]["template_id"] == result["columns"][2]["bottom"]["id"]


def test_two_n_additions_share_one_product_column_without_double_counting():
    result = (opened("AAA", "CCC").add_n("left", "G").add_n("right", "C")
              .synapse(right_start=4).ligate().run())
    assert result["sequence"] == "AAAGCCC"
    assert result["columns"][3]["top"]["kind"] == "n_addition"
    assert result["columns"][3]["bottom"]["kind"] == "n_addition"


@pytest.mark.parametrize("script,reason", [
    (opened(), "ligation"),
    (opened("AAA", "CCC", offset=1).synapse(right_start=3).ligate(), "incompatible"),
    (opened("AAA", "CCC").synapse(right_start=4).fill("top").ligate(), "template"),
    (opened().synapse(right_start=4).ligate(), "gaps"),
    (DuplexRepair("AAA", "CCC").synapse(right_start=3).ligate(), "opened"),
    (opened().synapse(right_start=2).ligate(), "flap"),
    (opened().synapse(right_start=3).ligate().fill("top"), "follow"),
])
def test_unsupported_or_unresolved_products_fail_explicitly(script, reason):
    with pytest.raises(ValueError, match=reason):
        script.run()


@pytest.mark.parametrize("value", [True, -1, 1.5, 2**64])
def test_python_coordinates_and_operation_counts_do_not_coerce(value):
    with pytest.raises(ValueError):
        DuplexRepair("AAA", "CCC").open("left", offset=value)
    if value != -1:
        with pytest.raises(ValueError):
            DuplexRepair("AAA", "CCC").synapse(right_start=value)
    else:
        assert DuplexRepair("AAA", "CCC").synapse(right_start=value).operations[-1]["right_start"] == -1


def test_script_is_immutable_and_work_bounds_are_enforced():
    base = DuplexRepair("AAA", "CCC")
    first = base.open("left", offset=0)
    assert base.operations == []
    operations = first.operations
    operations[0]["offset"] = 99
    assert first.operations[0]["offset"] == 0
    with pytest.raises(ValueError, match="budget"):
        DuplexRepair("AAA", "CCC", max_operations=0).open("left", offset=0)
    with pytest.raises(ValueError, match="budget"):
        opened(max_bases=5).synapse(right_start=3).ligate().run()


def test_opposite_strand_restoration_and_orientation_symmetry():
    left = (opened("ACG", "TTA").resect("left", "top", count=1)
            .synapse(right_start=3).fill("top").ligate().run())
    right = (opened("TAA", "CGT").resect("right", "bottom", count=1)
             .synapse(right_start=3).fill("bottom").ligate().run())
    complement = str.maketrans("ACGT", "TGCA")
    assert right["sequence"] == left["sequence"].translate(complement)[::-1]
    assert right["columns"][3]["bottom"]["kind"] == "templated_synthesis"


def test_unknown_operations_and_unit_operation_fields_are_rejected():
    script = opened().synapse(right_start=3).ligate()
    data = script.to_dict()
    data["operations"][-1]["invented"] = True
    with pytest.raises(ValueError, match="invalid duplex repair history"):
        DuplexRepair.from_dict(data).run()


def test_unpaired_three_prime_terminus_cannot_bootstrap_gap_synthesis():
    script = (DuplexRepair("AAA", "CCC").open("left", offset=0)
              .open("right", offset=1).resect("left", "bottom", count=2)
              .synapse(right_start=4).fill("top").fill("bottom").ligate())
    with pytest.raises(ValueError, match="paired 3-prime primer"):
        script.run()


def test_exhaustive_dinucleotide_opening_preserves_physical_bases_and_symmetry():
    from itertools import product

    complement = str.maketrans("ACGT", "TGCA")
    substrates = ["".join(bases) for bases in product("ACGT", repeat=2)]
    accepted = rejected = 0
    for left, right in product(substrates, repeat=2):
        # A blunt pair always joins. Symmetric one-base opening creates two
        # original-strand overhangs, compatible exactly when the facing bases
        # complement. This expectation is specified independently of the kernel.
        for offset in (0, 1):
            script = opened(left, right, offset).synapse(right_start=2).ligate()
            if offset and left[-1].translate(complement) != right[0]:
                with pytest.raises(ValueError, match="incompatible"):
                    script.run()
                rejected += 1
                continue
            result = script.run()
            assert result["sequence"] == left + right
            assert result["complementary_sequence"] == (left + right).translate(complement)[::-1]
            bases = [base for column in result["columns"] for base in (column["top"], column["bottom"])]
            assert {base["id"] for base in bases} == set(range(8))
            assert all(base["kind"] == "original_strand" and base["template_id"] is None for base in bases)
            assert result["removed"] == []
            mirrored = (opened(right.translate(complement)[::-1], left.translate(complement)[::-1], offset)
                        .synapse(right_start=2).ligate().run())
            assert mirrored["sequence"] == result["complementary_sequence"]
            accepted += 1
    assert (accepted, rejected) == (320, 192)
