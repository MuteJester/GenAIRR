"""Literal edit/ancestry oracles for native duplex founder family continuations.

These controls use equal authored history weights and existing neutral family
parameters. They do not calibrate biochemical frequencies from measured data.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import warnings

from GenAIRR import Experiment, RefDataConfig, _engine
from GenAIRR.experimental import CoupledDuplexRepair, DuplexHistory, DuplexJoin, DuplexRepairModel


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
FOUNDER_FIELDS = ("product", "evidence", "bindings", "selection", "founder_column_body_claims")
CALLS = {"v_call": "V*01", "d_call": "D*01", "j_call": "J*01"}


def source(inverted, order):
    """A restored V landmark, including nonpalindromic inverted D controls."""
    ref = RefDataConfig.vdj()
    ref.add_v_allele("V*01", "V", b"AAATGT", anchor=3)
    ref.add_d_allele("D*01", "D", b"ACT")
    ref.add_j_allele("J*01", "J", b"TTTAAA", anchor=0)
    repair = CoupledDuplexRepair("AAATGT", "ACT", "TTTAAA", d_inverted=inverted)
    for name in order:
        join = DuplexJoin(name).open("left", offset=0).open("right", offset=0)
        if name == "vd":
            join = join.resect("left", "top", count=1).synapse(right_start=6).fill("top")
        else:
            join = join.synapse(right_start=3)
        repair = repair.join(join.ligate())
    # Identical physical scripts with distinct IDs remain distinct equal-mass
    # authored choices; no default is inferred from the resulting DNA.
    rows = [DuplexHistory(name, repair, "V*01", "D*01", "J*01")
            for name in ("restored-a", "restored-b")]
    model = DuplexRepairModel.bind(ref, rows, reference_ends_are_coding_ends=True)
    expected = "AAATGT" + ("AGT" if inverted else "ACT") + "TTTAAA"
    return ref, model, expected


def claim_for_origin(origin, inverted):
    if origin < 6:
        return "v", origin
    if origin < 9:
        offset = origin - 6
        return "d", 2 - offset if inverted else offset
    return "j", origin - 9


def check_record(record, outcome, expected_founder, inverted, founder=None):
    """Reconstruct bytes, links and counters without using native projectors."""
    truth = record["duplex_repair"]
    observation = truth["observation"]
    assert truth["product"]["sequence"] == expected_founder
    assert truth["counts"]["founder_query"]["columns"] == len(expected_founder)
    assert truth["counts"]["physical"]["resected_nucleotides"] == 1
    assert truth["counts"]["physical"]["templated_synthesis"] == 1
    restored = truth["product"]["columns"][5]["top"]
    assert restored["kind"] == "templated_synthesis" and restored["original"] is None
    assert truth["evidence"]["placements"][str(restored["id"])]["role"] == "restored_body"
    assert all(record[key] == value for key, value in CALLS.items())
    assert record["d_inverted"] is inverted
    assert truth["selection"]["natural_weight"] == 1.0
    assert truth["selection"]["history_id"] in {"restored-a", "restored-b"}
    if founder is not None:
        for field in FOUNDER_FIELDS:
            assert truth[field] == founder[field]
        assert truth["counts"]["physical"] == founder["counts"]["physical"]
        assert truth["observation"]["founder_id"] == founder["observation"]["founder_id"]

    original_claims = truth["founder_column_body_claims"]
    assert len(original_claims) == len(expected_founder)
    for origin, column in enumerate(original_claims):
        assert column["query_position"] == origin
        assert [(claim["segment"], claim["reference_position"]) for claim in column["claims"]] == [claim_for_origin(origin, inverted)]

    query = list(expected_founder)
    links = [{"kind": "founder_column", "origin": i} for i in range(len(query))]
    removed = []
    counters = Counter()
    rev_comp = False
    for ordinal, edit in enumerate(observation["edits"]):
        assert edit["event_ordinal"] == ordinal
        action = edit["action"]
        kind = action["kind"]
        source_pass = edit["source_pass"]
        if kind == "substitution":
            at = action["at"]
            assert ord(query[at]) == action["old_base"]
            assert links[at] == action["link"]
            query[at] = chr(action["new_base"])
            if links[at]["kind"] == "founder_column":
                segment, germline = claim_for_origin(links[at]["origin"], inverted)
                assert action["segment"] == segment.upper()
                assert action["germline_pos"] == germline
            if source_pass in {"mutate.uniform", "mutate.s5f"}:
                counters["n_mutations"] += 1
                counters[f"n_{action['segment'].lower()}_mutations"] += 1
            elif source_pass == "corrupt.pcr":
                counters["applied_pcr"] += 1
            elif source_pass == "corrupt.quality":
                counters["applied_quality"] += 1
            else:
                raise AssertionError(f"unexpected substitution mechanism {source_pass}")
        elif kind == "insertion":
            at = action["at"]
            assert source_pass == "corrupt.indel"
            expected_link = {"kind": "insertion", "origin": {
                "event_ordinal": ordinal, "offset": 0, "original_base": action["base"]}}
            assert action["link"] == expected_link
            query.insert(at, chr(action["base"]))
            links.insert(at, expected_link)
            counters["n_indels"] += 1
            counters[f"n_{action['segment'].lower()}_indels"] += 1
        elif kind == "deletion":
            at = action["at"]
            assert ord(query.pop(at)) == action["base"]
            link = links.pop(at)
            assert link == action["link"]
            removed.append({"link": link, "query_position_at_removal": at,
                            "cause": action["cause"], "event_ordinal": ordinal})
            if action["cause"] == "observation_loss":
                assert source_pass in {"corrupt.end_loss.5", "corrupt.end_loss.3"}
                counters[f"end_loss_{source_pass[-1]}_length"] += 1
            else:
                assert action["cause"] == "structural_deletion" and source_pass == "corrupt.indel"
                counters["n_indels"] += 1
                counters[f"n_{action['segment'].lower()}_indels"] += 1
        else:
            assert kind == "reverse_complement" and source_pass == "corrupt.rev_comp"
            rev_comp = action["applied"]
    coding = "".join(query)
    assert outcome.final_simulation().bases().decode() == coding
    assert record["sequence"] == (coding.translate(COMPLEMENT)[::-1] if rev_comp else coding)
    assert observation["links"] == links and observation["removed"] == removed
    assert observation["rev_comp"] is rev_comp and record["rev_comp"] is rev_comp
    if record["sequence_alignment"] is not None:
        assert record["sequence_alignment"].replace("-", "") == coding
    surviving = [link["origin"] for link in links if link["kind"] == "founder_column"]
    lost = [entry["link"]["origin"] for entry in removed if entry["link"]["kind"] == "founder_column"]
    assert sorted(surviving + lost) == list(range(len(expected_founder)))
    assert len(truth["column_body_claims"]) == len(links)
    for q, column in enumerate(truth["column_body_claims"]):
        assert column["query_position"] == q
        if links[q]["kind"] == "founder_column":
            expected_claims = original_claims[links[q]["origin"]]["claims"]
            assert column["claims"] == expected_claims
        else:
            assert not column["claims"]

    for field in ("n_mutations", "n_v_mutations", "n_d_mutations", "n_j_mutations", "n_np_mutations",
                  "n_indels", "n_v_indels", "n_d_indels", "n_j_indels", "end_loss_5_length", "end_loss_3_length"):
        assert record[field] == counters[field], (field, record[field], counters[field])
    # In these controls every SHM stage precedes observation loss and uses the
    # same 15-base molecule, including repeated changes and reversions.
    assert math.isclose(record["mutation_rate"], counters["n_mutations"] / len(expected_founder), abs_tol=1e-12)
    for stage, field in (("pcr", "n_pcr_errors"), ("quality", "n_quality_errors")):
        attempts = sum(entry["sampled_attempts"] for entry in observation["attempt_stages"]
                       if entry["stage"] == stage)
        assert record[field] == attempts
        assert truth["applied_error_substitutions"][stage] == counters[f"applied_{stage}"]
    return coding, counters


def experiment_for(workflow, ref, model, inverted, productive, edited):
    experiment = Experiment.on(ref).recombine(duplex_repair=model)
    if inverted:
        experiment.invert_d(prob=1)
    if workflow == "fixed":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            experiment.expand_clones(n_clones=2, per_clone=4)
        experiment.mutate(model="uniform", count=2)
    elif workflow == "repertoire":
        experiment.clonal_repertoire(n_clones=2, size_distribution="lognormal", mu=2,
                                     sigma=0, max_size=7)
        experiment.mutate(model="uniform", count=2)
    else:
        experiment.clonal_lineage(n_clones=2, max_generations=3, n_max=18,
                                  n_sample=6, rate=0.08, lambda_base=2.5,
                                  selection_strength=0)
    if edited:
        (experiment.pcr_amplify(count=1).sequencing_errors(count=1)
         .end_loss_5prime(length=1).end_loss_3prime(length=1)
         .random_strand_orientation(prob=1))
    if productive:
        experiment.productive_only()
    return experiment


def run():
    counts = Counter()
    plans = []
    selections = Counter()
    dimensions = product(("fixed", "repertoire", "tree"), (("vd", "dj"), ("dj", "vd")),
                         (False, True), (False, True), (False, True), (False, True))
    for workflow, order, inverted, strict, productive, edited in dimensions:
        ref, model, expected = source(inverted, order)
        compiled = experiment_for(workflow, ref, model, inverted, productive, edited).compile()
        result = compiled.run_records(seed=23, strict=strict, validate_records=True)
        assert len(result.parents) == 2
        assert result.validate_families_with_parents(ref).ok
        founders = []
        founder_snapshots = []
        for parent in result.parents:
            parent_record = parent.airr_record(ref)
            check_record(parent_record, parent, expected, inverted)
            founders.append(parent_record["duplex_repair"])
            founder_snapshots.append(deepcopy(parent_record))
            selections[parent_record["duplex_repair"]["selection"]["history_id"]] += 1
            counts["founders"] += 1

        branch_maps = {}
        represented = Counter()
        for record, outcome in zip(result, result.outcomes, strict=True):
            parent_id = record["parent_id"]
            assert parent_id == record["clone_id"] and parent_id in (0, 1)
            check_record(record, outcome, expected, inverted, founders[parent_id])
            if workflow != "tree":
                assert record["n_mutations"] == 2
            if productive:
                assert record["productive"] is True
            if edited:
                assert record["n_pcr_errors"] == record["n_quality_errors"] == 1
                assert record["end_loss_5_length"] == record["end_loss_3_length"] == 1
                assert record["rev_comp"] is True
            branch_maps.setdefault(parent_id, set()).add(json.dumps(record["duplex_repair"]["observation"], sort_keys=True))
            represented[parent_id] += record.get("duplicate_count", 1)
            counts["records"] += 1
        if workflow == "fixed":
            assert len(result) == 8 and represented == {0: 4, 1: 4}
        elif workflow == "repertoire":
            assert represented == {0: 7, 1: 7}
        else:
            assert len(result.families) == len(result.lineage_trees) == 2
            assert represented == {0: 6, 1: 6}
            for family_id, (family, tree) in enumerate(zip(result.families, result.lineage_trees, strict=True)):
                nodes = tree.nodes()
                assert nodes[0].sequence == expected
                replayed_records = {}
                for node in nodes:
                    replayed = family.replay_node(node.id, ref)
                    replay_record = replayed.airr_record(ref)
                    coding, counters = check_record(replay_record, replayed, expected, inverted, founders[family_id])
                    assert coding == node.sequence
                    assert replayed.validate_record(ref) == []
                    if node.parent_id is None:
                        assert node.id == 0 and node.generation == 0 and node.mutations_from_parent == 0
                    else:
                        parent = replayed_records[node.parent_id]
                        old_edits = parent["duplex_repair"]["observation"]["edits"]
                        edits = replay_record["duplex_repair"]["observation"]["edits"]
                        assert edits[:len(old_edits)] == old_edits
                        edge = edits[len(old_edits):]
                        assert all(edit["source_pass"] == "mutate.s5f" and edit["action"]["kind"] == "substitution"
                                   for edit in edge)
                        assert len(edge) == node.mutations_from_parent
                        assert counters["n_mutations"] - parent["n_mutations"] == node.mutations_from_parent
                        assert node.generation == nodes[node.parent_id].generation + 1
                    replayed_records[node.id] = replay_record
                    counts["node_replays"] += 1
                for record in (record for record in result if record["parent_id"] == family_id):
                    node_id = record["lineage_node_id"]
                    replay_record = replayed_records[node_id]
                    assert nodes[node_id].observed
                    assert record["lineage_abundance"] == nodes[node_id].abundance
                    assert record["lineage_mutations_from_parent"] == nodes[node_id].mutations_from_parent
                    assert record["n_mutations"] == replay_record["n_mutations"]
                    if not edited:
                        assert record["sequence"] == replay_record["sequence"]
                        assert record["duplex_repair"] == replay_record["duplex_repair"]
                    else:
                        edits = record["duplex_repair"]["observation"]["edits"]
                        node_edits = replay_record["duplex_repair"]["observation"]["edits"]
                        assert edits[:len(node_edits)] == node_edits
                    counts["observed_node_comparisons"] += 1
        if workflow != "tree":
            assert all(len(maps) > 1 for maps in branch_maps.values())
        # No descendant may mutate the founder's query or serialized history.
        for parent, snapshot in zip(result.parents, founder_snapshots, strict=True):
            assert parent.airr_record(ref) == snapshot
        counts["plans"] += 1
        plans.append({"workflow": workflow, "join_order": list(order), "d_inverted": inverted,
                      "strict": strict, "productive_only": productive, "observation_suffix": edited,
                      "records": len(result), "represented_cells": dict(represented),
                      "distinct_observation_maps": {key: len(value) for key, value in branch_maps.items()}})
    assert set(selections) == {"restored-a", "restored-b"}
    return {"schema": "duplex_family_sweep_v1", "seed": 23,
            "extension_sha256": hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
            "counts": dict(counts), "errors": 0, "selected_history_ids": dict(selections),
            "selection_count_scope": "control plans reuse founder seeds; these correlated counts do not estimate or test population frequencies",
            "checks": ["literal coding-query and output-orientation reconstruction", "exact current/removed founder-link replay",
                       "canonical inverted-D placement oracle", "cumulative applied mutation and observation counter oracle",
                       "mutation-stage denominator before read loss", "one immutable repair founder per clone",
                       "actual parent-edge edit prefixes and counts", "all lineage nodes replay including unsampled nodes",
                       "represented-cell abundance conservation", "independent observed branch maps"],
            "scope": "authored 15-base restored-landmark controls; omitted unit history weights and existing neutral family controls; no empirical biochemical fitting",
            "plans": plans}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
