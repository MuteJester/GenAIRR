"""Independent physical-base oracles for two joins sharing one finite D duplex."""

from itertools import product
import json

import pytest

from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin

pytestmark = pytest.mark.behavior


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def opened(name, left=0, right=0):
    return DuplexJoin(name).open("left", offset=left).open("right", offset=right)


def blunt(name, right_start):
    return opened(name).synapse(right_start=right_start).ligate()


def assembled(v, d, j, *, inverted=False, order=("vd", "dj"), scripts=None, **limits):
    if scripts is None:
        scripts = {"vd": blunt("vd", len(v)), "dj": blunt("dj", len(d))}
    script = CoupledDuplexRepair(v, d, j, d_inverted=inverted, **limits)
    for name in order:
        script = script.join(scripts[name])
    return script


def original_id(references, segment, strand, position):
    offset = sum(len(references[s]) for s in ("v", "d", "j")[:("v", "d", "j").index(segment)])
    return 2 * (offset + position) + (strand == "bottom")


def assert_conservation(result, references):
    """Every physical nucleotide is either retained once or explicitly removed."""
    columns = result["columns"]
    current = [base for column in columns for base in (column["top"], column["bottom"])]
    removed = result["removed"]
    ledger = result["ledger"]
    by_id = {base["id"]: base for base in ledger}
    assert len(by_id) == len(ledger)
    assert set(by_id) == set(range(len(ledger)))
    current_ids = [base["id"] for base in current]
    removed_ids = [base["id"] for base in removed]
    assert len(set(current_ids)) == len(current_ids)
    assert len(set(removed_ids)) == len(removed_ids)
    assert not set(current_ids).intersection(removed_ids)
    assert set(current_ids).union(removed_ids) == set(by_id)
    assert all(base == by_id[base["id"]] for base in current + removed)

    initial_count = 2 * sum(map(len, references.values()))
    for identifier in range(initial_count):
        base = by_id[identifier]
        origin = base["original"]
        assert base["kind"] == "original_strand"
        assert base["template_id"] is None
        expected = references[origin["segment"]][origin["reference_position"]]
        if origin["strand"] == "bottom":
            expected = expected.translate(COMPLEMENT)
        assert base["base"] == expected
        assert identifier == original_id(
            references, origin["segment"], origin["strand"], origin["reference_position"]
        )
    for identifier in range(initial_count, len(ledger)):
        base = by_id[identifier]
        assert base["original"] is None
        if base["kind"] == "templated_synthesis":
            template = by_id[base["template_id"]]
            assert template["id"] < identifier  # Acyclic historical template graph.
            assert base["base"] == template["base"].translate(COMPLEMENT)
        else:
            assert base["kind"] == "n_addition"
            assert base["template_id"] is None

    assert result["sequence"] == "".join(column["top"]["base"] for column in columns)
    assert result["complementary_sequence"] == rc(result["sequence"])
    assert [column["position"] for column in columns] == list(
        range(columns[0]["position"], columns[0]["position"] + len(columns))
    )
    assert all(
        column["top"]["base"].translate(COMPLEMENT) == column["bottom"]["base"]
        for column in columns
    )


@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
@pytest.mark.parametrize("inverted", [False, True])
def test_blunt_coupled_joins_share_exactly_one_d_and_preserve_original_orientation(order, inverted):
    references = {"v": "ACG", "d": "AAT", "j": "TCC"}
    script = assembled(*references.values(), inverted=inverted, order=order)
    result = script.run()
    assert result["sequence"] == ("ACGATTTCC" if inverted else "ACGAATTCC")
    assert len(result["ledger"]) == 18
    assert result["removed"] == []
    assert result["origin_positions"] == {"v": 0, "d": 3, "j": 6}
    seams = {seam["join"]: (seam["top"], seam["bottom"]) for seam in result["ligation_seams"]}
    assert seams == {"vd": (3, 3), "dj": (6, 6)}
    assert_conservation(result, references)
    for position in range(3):
        column = result["columns"][3 + position]
        source_position = 2 - position if inverted else position
        for strand, original_strand in (
            ("top", "bottom" if inverted else "top"),
            ("bottom", "top" if inverted else "bottom"),
        ):
            assert column[strand]["original"] == {
                "segment": "d", "strand": original_strand, "reference_position": source_position
            }


def test_exhaustive_short_blunt_products_and_both_join_orders_match_literal_oracle():
    # All single-base V/J substrates and every one- or two-base D, with both
    # orientations/orders. Expected identities are direct canonical coordinates;
    # this oracle never invokes a local repair product to construct the answer.
    d_substrates = list("ACGT") + ["".join(bases) for bases in product("ACGT", repeat=2)]
    checked = 0
    for v, d, j, inverted in product("ACGT", d_substrates, "ACGT", (False, True)):
        references = {"v": v, "d": d, "j": j}
        answers = []
        for order in (("vd", "dj"), ("dj", "vd")):
            result = assembled(v, d, j, inverted=inverted, order=order).run()
            assert result["sequence"] == v + (rc(d) if inverted else d) + j
            assert len(result["ledger"]) == 2 * (2 + len(d))
            assert result["removed"] == []
            assert_conservation(result, references)
            expected_top = [0]
            expected_bottom = [1]
            positions = reversed(range(len(d))) if inverted else range(len(d))
            for p in positions:
                expected_top.append(original_id(references, "d", "bottom" if inverted else "top", p))
                expected_bottom.append(original_id(references, "d", "top" if inverted else "bottom", p))
            expected_top.append(original_id(references, "j", "top", 0))
            expected_bottom.append(original_id(references, "j", "bottom", 0))
            assert [column["top"]["id"] for column in result["columns"]] == expected_top
            assert [column["bottom"]["id"] for column in result["columns"]] == expected_bottom
            answers.append(result["columns"])
            checked += 1
        assert answers[0] == answers[1]
    assert checked == 1280


@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
@pytest.mark.parametrize("inverted", [False, True])
def test_first_join_can_supply_the_stem_for_second_opening_without_an_original_d_core(order, inverted):
    # h5+h3 == D length. A blanket original-core arithmetic test would reject
    # this valid successive history; the partner strands supply pairing.
    references = {"v": "AAG", "d": "TG" if inverted else "CA", "j": "TAA"}
    scripts = {
        "vd": opened("vd", 1, 1).synapse(right_start=3).ligate(),
        "dj": opened("dj", 1, 1).synapse(right_start=2).ligate(),
    }
    result = assembled(*references.values(), inverted=inverted, order=order, scripts=scripts).run()
    assert result["sequence"] == "AAGCATAA"
    assert result["complementary_sequence"] == "TTATGCTT"
    assert len(result["ledger"]) == 16  # Opening itself synthesizes nothing.
    assert result["removed"] == []
    assert_conservation(result, references)
    expected_top = [0, 2, 4, 5, 7, 6, 12, 14] if inverted else [0, 2, 4, 5, 8, 9, 12, 14]
    expected_bottom = [1, 3, 9, 8, 10, 11, 13, 15] if inverted else [1, 3, 6, 7, 10, 11, 13, 15]
    assert [column["top"]["id"] for column in result["columns"]] == expected_top
    assert [column["bottom"]["id"] for column in result["columns"]] == expected_bottom
    seams = {seam["join"]: (seam["top"], seam["bottom"]) for seam in result["ligation_seams"]}
    assert seams == {"vd": (4, 2), "dj": (6, 4)}


@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
def test_inverted_d_one_sided_opening_uses_the_opposite_canonical_end(order):
    references = {"v": "AAG", "d": "TG", "j": "TAA"}
    scripts = {
        "vd": opened("vd", 1, 1).synapse(right_start=3).ligate(),
        "dj": blunt("dj", 2),
    }
    result = assembled(*references.values(), inverted=True, order=order, scripts=scripts).run()
    assert result["sequence"] == "AAGCATAA"
    # Physical D5 is canonical D3 after inversion. Its transferred physical
    # base keeps the canonical bottom-strand identity at the last D coordinate.
    assert result["columns"][2]["bottom"]["original"] == {
        "segment": "d", "strand": "bottom", "reference_position": 1
    }
    assert result["columns"][5]["top"]["original"] == {
        "segment": "j", "strand": "top", "reference_position": 0
    }
    assert_conservation(result, references)


@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
def test_two_opposing_n_pairs_add_two_columns_not_four(order):
    references = {"v": "AAA", "d": "CCC", "j": "GGG"}
    scripts = {
        "vd": opened("vd").add_n("left", "G").add_n("right", "C").synapse(right_start=4).ligate(),
        "dj": opened("dj").add_n("left", "A").add_n("right", "T").synapse(right_start=4).ligate(),
    }
    result = assembled(*references.values(), order=order, scripts=scripts).run()
    assert result["sequence"] == "AAAGCCCAGGG"
    assert result["origin_positions"] == {"v": 0, "d": 4, "j": 8}
    assert len(result["ledger"]) == 22
    assert sum(base["kind"] == "n_addition" for base in result["ledger"]) == 4
    for p in (3, 7):
        assert result["columns"][p]["top"]["kind"] == "n_addition"
        assert result["columns"][p]["bottom"]["kind"] == "n_addition"
    assert_conservation(result, references)


@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
def test_restored_d_base_keeps_synthesis_ancestry_through_the_other_join(order):
    references = {"v": "ACG", "d": "CAT", "j": "TAA"}
    scripts = {
        "vd": (opened("vd").resect("right", "top", count=1)
               .synapse(right_start=3).fill("top").ligate()),
        "dj": blunt("dj", 3),
    }
    result = assembled(*references.values(), order=order, scripts=scripts).run()
    assert result["sequence"] == "ACGCATTAA"
    restored = result["columns"][3]["top"]
    assert restored["kind"] == "templated_synthesis"
    assert restored["id"] == 18
    assert restored["template_id"] == original_id(references, "d", "bottom", 0)
    assert [base["id"] for base in result["removed"]] == [original_id(references, "d", "top", 0)]
    assert_conservation(result, references)


def test_second_hairpin_opening_can_transfer_a_first_join_synthesis_without_relabeling_it():
    references = {"v": "ATG", "d": "CAA", "j": "TAA"}
    scripts = {
        "dj": (opened("dj").resect("left", "top", count=2)
               .synapse(right_start=3).fill("top").ligate()),
        "vd": opened("vd", 2, 2).synapse(right_start=3).ligate(),
    }
    result = assembled(*references.values(), order=("dj", "vd"), scripts=scripts).run()
    assert result["sequence"] == "ATGCAATAA"
    transferred = result["columns"][1]["bottom"]
    assert transferred["kind"] == "templated_synthesis"
    assert transferred["original"] is None
    assert transferred["template_id"] == original_id(references, "d", "bottom", 1)
    assert len(result["ledger"]) == 20
    assert {base["id"] for base in result["removed"]} == {
        original_id(references, "d", "top", 1), original_id(references, "d", "top", 2)
    }
    assert_conservation(result, references)


def test_template_links_survive_template_removal_and_second_generation_synthesis():
    references = {"v": "ACG", "d": "CAA", "j": "TAA"}
    scripts = {
        "dj": (opened("dj").resect("left", "top", count=2)
               .synapse(right_start=3).fill("top").ligate()),
        "vd": (opened("vd").resect("right", "bottom", count=2)
               .synapse(right_start=3).fill("bottom").ligate()),
    }
    result = assembled(*references.values(), order=("dj", "vd"), scripts=scripts).run()
    assert result["sequence"] == "ACGCAATAA"
    top = result["columns"][4]["top"]
    bottom = result["columns"][4]["bottom"]
    assert top["kind"] == bottom["kind"] == "templated_synthesis"
    assert top["template_id"] == original_id(references, "d", "bottom", 1)
    assert top["template_id"] in {base["id"] for base in result["removed"]}
    assert bottom["template_id"] == top["id"]
    assert len(result["ledger"]) == 22
    assert_conservation(result, references)


def test_first_join_synthesis_can_support_opening_a_single_base_d():
    # The first join adds a C/G pair immediately before D. Its synthesized
    # bottom G supplies the retained pairing while D's sole original bottom
    # base moves around its remaining hairpin. Requiring offset < D length at
    # every join would wrongly reject this current-component substrate.
    references = {"v": "AAG", "d": "C", "j": "GAA"}
    scripts = {
        "vd": (opened("vd").add_n("left", "C")
               .synapse(right_start=4).fill("bottom").ligate()),
        "dj": opened("dj", 1, 1).synapse(right_start=1).ligate(),
    }
    result = assembled(*references.values(), scripts=scripts).run()
    assert result["sequence"] == "AAGCCGAA"
    assert [column["top"]["id"] for column in result["columns"]] == [0, 2, 4, 14, 6, 7, 10, 12]
    assert result["columns"][3]["bottom"]["template_id"] == 14
    assert len(result["ledger"]) == 16
    assert_conservation(result, references)


@pytest.mark.parametrize("order", [("vd", "dj"), ("dj", "vd")])
def test_coupled_history_round_trips_exactly_and_serialized_copies_do_not_mutate_scripts(order):
    script = assembled("ACG", "AAT", "TCC", inverted=True, order=order)
    encoded = json.loads(json.dumps(script.to_dict()))
    result = script.run()
    assert [join["join"] for join in result["joins"]] == list(order)
    assert CoupledDuplexRepair.from_dict(encoded).run() == result
    encoded["joins"][0]["operations"][0]["offset"] = 99
    assert script.to_dict()["joins"][0]["operations"][0]["offset"] == 0
    with pytest.raises(ValueError):
        CoupledDuplexRepair.from_dict(encoded).run()
    invalid = script.to_dict()
    invalid["joins"][1]["operations"][-1]["invented"] = True
    with pytest.raises(ValueError):
        CoupledDuplexRepair.from_dict(invalid).run()


def test_incomplete_duplicate_and_unresolved_joins_do_not_export_a_product():
    with pytest.raises(ValueError):
        CoupledDuplexRepair("AAA", "CCC", "GGG").join(blunt("vd", 3)).run()
    with pytest.raises(ValueError):
        CoupledDuplexRepair("AAA", "CCC", "GGG").join(blunt("vd", 3)).join(blunt("vd", 3)).run()
    with pytest.raises(ValueError):
        assembled("AAA", "CCC", "GGG", scripts={
            "vd": blunt("vd", 3),
            "dj": opened("dj").synapse(right_start=4).ligate(),
        }).run()


def test_overprocessing_is_rejected_instead_of_regrowing_or_clamping_the_finite_d():
    for first in ("vd", "dj"):
        d_side = "right" if first == "vd" else "left"
        impossible = DuplexJoin(first).open(d_side, offset=1)
        with pytest.raises(ValueError):
            CoupledDuplexRepair("AAA", "C", "GGG").join(impossible).join(
                blunt("dj" if first == "vd" else "vd", 1 if first == "vd" else 3)
            ).run()
    with pytest.raises(ValueError):
        assembled("AAA", "CC", "GGG", scripts={
            "vd": opened("vd").resect("right", "top", count=2).synapse(right_start=3).ligate(),
            "dj": blunt("dj", 2),
        }).run()


@pytest.mark.parametrize("first", ["vd", "dj"])
def test_second_join_resection_cannot_cut_through_a_previous_ligation(first):
    if first == "vd":
        scripts = {
            "vd": blunt("vd", 3),
            "dj": (opened("dj").resect("left", "bottom", count=2)
                   .synapse(right_start=2).fill("bottom").ligate()),
        }
        order = ("vd", "dj")
    else:
        scripts = {
            "dj": blunt("dj", 2),
            "vd": (opened("vd").resect("right", "top", count=2)
                   .synapse(right_start=3).fill("top").ligate()),
        }
        order = ("dj", "vd")
    with pytest.raises(ValueError, match="previous ligation seam"):
        assembled("AAA", "CC", "GGG", order=order, scripts=scripts).run()


def test_operation_and_synthesis_budgets_cover_both_joins_and_removed_bases():
    with pytest.raises(ValueError):
        assembled("AAA", "CCC", "GGG", max_operations=7).run()
    scripts = {
        "vd": opened("vd").add_n("left", "G").add_n("right", "C").synapse(right_start=4).ligate(),
        "dj": opened("dj").add_n("left", "A").add_n("right", "T").synapse(right_start=4).ligate(),
    }
    with pytest.raises(ValueError):
        assembled("AAA", "CCC", "GGG", scripts=scripts, max_bases=21).run()
    assert assembled("AAA", "CCC", "GGG", scripts=scripts, max_bases=22).run()["sequence"] == "AAAGCCCAGGG"
    restoration = {
        "vd": opened("vd").resect("right", "top", count=1).synapse(right_start=3).fill("top").ligate(),
        "dj": blunt("dj", 3),
    }
    # There are still18 bases in the product, but19 physical bases have existed.
    with pytest.raises(ValueError):
        assembled("ACG", "CAT", "TAA", scripts=restoration, max_bases=18).run()


@pytest.mark.parametrize("value", [0, 1, "false", None])
def test_inversion_requires_an_explicit_boolean(value):
    with pytest.raises(ValueError):
        CoupledDuplexRepair("AAA", "CCC", "GGG", d_inverted=value)
