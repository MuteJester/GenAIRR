"""Literal edit replay and ancestry conservation across native duplex continuations."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path

from GenAIRR import RefDataConfig, _engine
from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin, project_duplex_simulation


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def check(seed, sim, projected):
    rec = projected["record"]
    truth = rec["duplex_repair"]
    obs = truth["observation"]
    query = list(seed)
    for ordinal, edit in enumerate(obs["edits"]):
        assert edit["event_ordinal"] == ordinal
        action = edit["action"]
        kind = action["kind"]
        if kind == "substitution":
            assert ord(query[action["at"]]) == action["old_base"]
            query[action["at"]] = chr(action["new_base"])
        elif kind == "insertion":
            query.insert(action["at"], chr(action["base"]))
        elif kind == "deletion":
            assert ord(query.pop(action["at"])) == action["base"]
        else:
            assert kind == "reverse_complement"
    current = "".join(query)
    assert sim.bases().decode() == current
    assert rec["sequence"] == (current.translate(COMPLEMENT)[::-1] if obs["rev_comp"] else current)
    if rec["sequence_alignment"] is not None:
        assert rec["sequence_alignment"].replace("-", "") == current
    assert len(obs["links"]) == len(current)
    founder = [x["origin"] for x in obs["links"] if x["kind"] == "founder_column"]
    founder += [x["link"]["origin"] for x in obs["removed"] if x["link"]["kind"] == "founder_column"]
    assert sorted(founder) == list(range(len(seed)))
    for column in truth["column_body_claims"]:
        q = column["query_position"]
        for claim in column["claims"]:
            link = obs["links"][q]
            assert link["kind"] == "founder_column"
            original_claims = truth["founder_column_body_claims"][link["origin"]]["claims"]
            normalized = dict(claim, support_ids=sorted(claim["support_ids"]))
            assert normalized in original_claims
    assert truth["counts"]["founder_query"]["columns"] == len(seed)
    assert truth["product"]["sequence"] == seed
    assert all(rec[f"{s}_call"] == f"{s.upper()}*01" for s in ("v", "d", "j"))


def run():
    # Distinct literal bases cover orientation, GC/AT ends and repeated motifs.
    substrates = ("AA", "AC", "CA", "CG", "GT", "TG", "TC", "TT")
    counts = Counter()
    methods = [
        ("push_mutate_uniform", ([(2, 1.0)],), {}),
        ("push_corrupt_pcr", ([(2, 1.0)],), {}),
        ("push_corrupt_quality", ([(2, 1.0)],), {}),
        ("push_corrupt_ns", ([(2, 1.0)],), {}),
        ("push_corrupt_indel", ([(2, 1.0)],), {"insertion_prob": 1.0}),
        ("push_corrupt_indel", ([(2, 1.0)],), {"insertion_prob": 0.0}),
        ("push_corrupt_5prime_loss", ([(3, 1.0)],), {}),
        ("push_corrupt_3prime_loss", ([(3, 1.0)],), {}),
        ("push_corrupt_rev_comp", (1.0,), {}),
    ]
    for v, d, j in product(substrates, repeat=3):
        ref = RefDataConfig.vdj()
        ref.add_v_allele("V*01", "V", v.encode())
        ref.add_d_allele("D*01", "D", d.encode())
        ref.add_j_allele("J*01", "J", j.encode())
        plans = []
        for method, args, kwargs in methods:
            plan = _engine.PassPlan()
            getattr(plan, method)(*args, **kwargs)
            plans.append(plan.compile(refdata=ref, strict=True, allow_curatable_refdata=True))
        for inverted, order in product((False, True), (("vd", "dj"), ("dj", "vd"))):
            script = CoupledDuplexRepair(v, d, j, d_inverted=inverted)
            for name in order:
                script = script.join(DuplexJoin(name).open("left", offset=0)
                    .open("right", offset=0).synapse(right_start=2).ligate())
            seed = script.seed(ref, v_call="V*01", d_call="D*01", j_call="J*01", reference_ends_are_coding_ends=True)
            expected = v + (d.translate(COMPLEMENT)[::-1] if inverted else d) + j
            initial = project_duplex_simulation(seed, ref)
            check(expected, seed, initial)
            counts["founders"] += 1
            for i, compiled in enumerate(plans):
                child = compiled.run_from(seed, counts["founders"] * 13 + i).final_simulation()
                check(expected, child, project_duplex_simulation(child, ref))
                counts["continuations"] += 1
            assert project_duplex_simulation(seed, ref) == initial
    return {"schema":"duplex_observation_sweep_v1", "counts":dict(counts),"errors":0,
        "extension_sha256":hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
        "checks":["independent literal edit replay", "founder-column conservation", "current certified claim remapping",
                  "query/alignment conservation", "immutable repair product and counters", "both join orders and D orientations"],
        "scope":"short authored controls; no fitted biochemical probabilities"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
