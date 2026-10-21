"""Finite authored law and metadata oracle using literal sequence arithmetic."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path

from GenAIRR import Experiment, RefDataConfig, _engine
from GenAIRR.experimental import CoupledDuplexRepair, DuplexHistory, DuplexJoin, DuplexRepairModel

MASK = (1 << 64) - 1
GAMMA = 0x9E3779B97F4A7C15
COMP = str.maketrans("ACGT", "TGCA")


def uniforms(seed):
    """Independent published SplitMix64 arithmetic, isolated from native RNG."""
    state = seed or GAMMA
    while True:
        state = (state + GAMMA) & MASK
        word = state
        word = ((word ^ (word >> 30)) * 0xBF58476D1CE4E5B9) & MASK
        word = ((word ^ (word >> 27)) * 0x94D049BB133111EB) & MASK
        yield ((word ^ (word >> 31)) >> 11) / (1 << 53)


def choose(choices, draw):
    target = draw * sum(weight for _, weight in choices)
    cumulative = 0
    for value, weight in choices:
        cumulative += weight
        if target < cumulative:
            return value
    raise AssertionError("literal positive-support selection was incomplete")


def script(dna, inverted, order, kind):
    s = CoupledDuplexRepair("AAATGT", dna, "TTTAAA", d_inverted=inverted)
    for name in order:
        join = DuplexJoin(name).open("left", offset=0).open("right", offset=0)
        if name == "vd":
            if kind == "restore":
                join = join.resect("left", "top", count=1)
            extra = {"triplet":"GCT", "frame_shift":"A"}.get(kind, "")
            if extra:
                join = join.add_n("left", extra).add_n("right", extra.translate(COMP)[::-1])
            join = join.synapse(right_start=6 + len(extra))
            if kind == "restore":
                join = join.fill("top")
        else:
            join = join.synapse(right_start=3)
        s = s.join(join.ligate())
    return s


def run(records_per_plan=512):
    r = RefDataConfig.vdj()
    for i in (0, 1):
        r.add_v_allele(f"V{i}*01", f"V{i}", b"AAATGT", anchor=3)
    for i, dna in enumerate(("ACT", "CCC")):
        r.add_d_allele(f"D{i}*01", f"D{i}", dna.encode())
    r.add_j_allele("J*01", "J", b"TTTAAA", anchor=0)
    count = Counter()
    plan_reports = []
    for order, productive, inversion_probability in product((("vd", "dj"), ("dj", "vd")), (False, True), (0.0, 0.37, 1.0)):
        rows, facts = [], []
        for v, d, inverted in product((0, 1), (0, 1), (False, True)):
            for kind, weight in (("blunt", 1), ("restore", 3), ("triplet", 2), ("frame_shift", 4)):
                dna = ("ACT", "CCC")[d]
                repair = script(dna, inverted, order, kind)
                row = DuplexHistory(f"{v}-{d}-{inverted}-{kind}", repair, f"V{v}*01", f"D{d}*01", "J*01", weight)
                rows.append(row)
                extra = {"triplet":"GCT", "frame_shift":"A"}.get(kind, "")
                expected = "AAATGT" + extra + (dna.translate(COMP)[::-1] if inverted else dna) + "TTTAAA"
                facts.append((v, d, inverted, kind, weight, expected))
        m = DuplexRepairModel.bind(r, rows, reference_ends_are_coding_ends=True)
        exp = Experiment.on(r).recombine(duplex_repair=m, v_allele_weights={"V1*01":3}, d_allele_weights={"D1*01":2}).invert_d(prob=inversion_probability)
        if productive:
            exp.productive_only()
        compiled = exp.compile()
        selected = Counter()
        for seed, outcome in enumerate(compiled.stream(n=records_per_plan, seed=1, strict=True), 1):
            draws = uniforms(seed)
            v = choose(((0, 1), (1, 3)), next(draws))
            d = choose(((0, 1), (1, 2)), next(draws))
            next(draws)  # J's one-element natural distribution.
            inverted = next(draws) < inversion_probability
            eligible = [(i, fact[4]) for i, fact in enumerate(facts)
                        if fact[:3] == (v, d, inverted) and (not productive or fact[3] != "frame_shift")]
            index = choose(eligible, next(draws))
            trace = outcome.trace()
            assert trace.find("sample_allele.v").value == v
            assert trace.find("sample_allele.d").value == d
            assert trace.find("sample_allele.d.inverted").value == inverted
            assert trace.find("duplex_repair.v1.history").value == index
            rec = outcome.airr_record(r)
            fact = facts[index]
            assert rec["sequence"] == fact[5]
            assert rec["sequence_alignment"] == rec["sequence"]
            assert rec["productive"] is (fact[3] != "frame_shift")
            assert rec["np1"] == {"triplet":"GCT", "frame_shift":"A"}.get(fact[3], "")
            assert rec["np1_length"] == len(rec["np1"]) and rec["np2_length"] == 0
            assert rec["n_mutations"] == rec["n_indels"] == 0
            assert rec["p_v_3_length"] is None
            truth = rec["duplex_repair"]
            assert truth["selection"]["history_id"] == rows[index].id
            assert truth["selection"]["natural_weight"] == fact[4]
            assert truth["selection"]["conditioning_law"] == "finite_authored_history_local_conditioning_v1"
            assert truth["counts"]["physical"]["resected_nucleotides"] == (1 if fact[3] == "restore" else 0)
            assert truth["current_paths"]["d"]["germline_alignment"] == (("ACT", "CCC")[d].translate(COMP)[::-1] if inverted else ("ACT", "CCC")[d])
            selected[index] += 1
            count["records"] += 1
        count["plans"] += 1
        plan_reports.append({"join_order":order, "productive_only":productive,
                             "inversion_probability":inversion_probability,
                             "selected_rows":{rows[i].id:n for i,n in sorted(selected.items())}})
    return {"schema":"duplex_experiment_sweep_v1", "extension_sha256":hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
            "counts":dict(count),"errors":0,"records_per_plan":records_per_plan,
            "checks":["independent exact PRNG/weighted-choice arithmetic", "local conditioning hierarchy", "generating calls and canonical inverted D", "literal query, body alignment and NP intervals", "same-DNA histories remain distinct", "repair resection separate from downstream counters"],
            "scope":"authored 15/16/18-base controls; no empirical biochemical rates",
            "plans":plan_reports}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--records-per-plan", type=int, default=512)
    args = parser.parse_args()
    if args.records_per_plan < 1:
        parser.error("--records-per-plan must be positive")
    args.output.write_text(json.dumps(run(args.records_per_plan), indent=2, allow_nan=False) + "\n")
