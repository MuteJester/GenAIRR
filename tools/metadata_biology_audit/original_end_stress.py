"""Reproducible original-end metadata sweep using only GenAIRR and stdlib.

Run from the repository: .venv/bin/python -m tools.metadata_biology_audit.original_end_stress
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import GenAIRR as ga
from GenAIRR import _engine

CONFIGS = [f"{species}_{locus}" for species in ("human", "mouse")
           for locus in ("igh", "igk", "igl", "tcra", "tcrb", "tcrd", "tcrg")]
COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def check(record, reference):
    meta = record["coding_end_processing"]
    assert meta["model"] == "original_coding_end_v1"
    assert meta["recombination_complete"]
    assert meta["coordinate_system"] == "zero_based_coding_oriented"
    coding = record["sequence"]
    if record["rev_comp"]:
        coding = coding.translate(COMPLEMENT)[::-1]
    historical = dict.fromkeys(("V_3", "D_5", "D_3", "J_5"), 0)
    observed = historical.copy()
    source_pools = {
        "V": [reference.v_allele(i) for i in range(reference.v_pool_size())],
        "D": [reference.d_allele(i) for i in range(reference.d_pool_size())],
        "J": [reference.j_allele(i) for i in range(reference.j_pool_size())],
    }
    for origin in meta["p_origins"]:
        source = next(a for a in source_pools[origin["source_segment"]]
                      if a.name == origin["original_allele"])
        sequence = source.seq().decode()
        base = sequence[origin["template_position"]]
        if origin["source_orientation"] == "forward":
            base = base.translate(COMPLEMENT)
        assert base.upper() == origin["original_base"].upper()
        historical[origin["assembled_end"]] += 1
        position = origin["query_position"]
        if position is None:
            assert origin["current_base"] is None
        else:
            observed[origin["assembled_end"]] += 1
            assert coding[position].upper() == origin["current_base"].upper()
    assert historical == meta["retained_p_at_join"]
    assert observed == meta["observed_p"]
    for end, count in historical.items():
        assert record[f"p_{end.lower()}_length"] == count
    for choice in meta["choices"]:
        segment = choice["segment"]
        assert record[f"{segment.lower()}_trim_5"] == max(choice["five_prime_outcome"], 0)
        assert record[f"{segment.lower()}_trim_3"] == max(choice["three_prime_outcome"], 0)
        for canonical, outcome in (("5", choice["five_prime_outcome"]), ("3", choice["three_prime_outcome"])):
            origins = [o for o in meta["p_origins"]
                       if o["choice_id"] == choice["choice_id"] and o["canonical_end"] == canonical]
            assert len(origins) == max(-outcome, 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20261020)
    parser.add_argument("--productive-indels", action="store_true", help="include strict two-event indel batches in productive profiles")
    parser.add_argument("--output", type=Path, default=Path("audit-docs/original_coding_end_stress_2026_10_20.json"))
    args = parser.parse_args()
    results = []
    rejections = []
    for name in CONFIGS:
        for productive in (False, True):
            exp = ga.Experiment.on(name).allow_curatable_refdata().recombine(
                np1_lengths=[(n, 1) for n in range(7)],
                **({"np2_lengths": [(n, 1) for n in range(7)]}
                   if ga.Experiment.on(name).chain_type == "vdj" else {}),
            )
            ends = dict(v_3=[(-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1)],
                        j_5=[(-1, 1), (0, 1), (1, 1)])
            if exp.chain_type == "vdj":
                ends.update(d_5=[(-1, 1), (0, 1), (1, 1)], d_3=[(-2, 1), (0, 1), (2, 1)])
            exp.process_coding_ends(reference_ends_are_coding_ends=True, **ends)
            if exp.chain_type == "vdj":
                exp.invert_d(prob=0.5)
            if productive:
                exp.productive_only()
                if "_ig" in name:
                    exp.mutate(model="uniform", rate=0.12)
                if args.productive_indels:
                    exp.polymerase_indels(count=2)
            else:
                if "_ig" in name:
                    exp.mutate(model="uniform", rate=0.15)
                else:
                    exp.sequencing_errors(rate=0.15)
                (exp.polymerase_indels(count=4)
                 .ambiguous_base_calls(count=2).end_loss_5prime(length=10)
                 .end_loss_3prime(length=5).random_strand_orientation(prob=0.5))
            try:
                compiled = exp.compile()
            except ValueError as error:
                # Record existing cartridge integrity failures explicitly.
                # Never bypass fatal validation to obtain a sweep result.
                reference_error = "[fatal]" in str(error)
                support_error = productive and "no productive completion" in str(error)
                if not reference_error and not support_error:
                    raise
                rejection = {"config": name, "productive": productive,
                             "status": "reference_rejected" if reference_error else "productive_support_rejected", "reason": str(error)}
                rejections.append(rejection)
                print(json.dumps(rejection), flush=True)
                continue
            result = compiled.run_records(n=args.records, seed=args.seed, strict=True)
            report = result.validate_records(exp.refdata)
            assert report.ok, (name, productive, report.failures[:1])
            for record in result:
                check(record, exp.refdata)
                if productive:
                    assert record["productive"] and record["vj_in_frame"] and not record["stop_codon"]
                    assert record["junction_aa"].startswith("C")
                    assert record["junction_aa"].endswith(("F", "W"))
                    assert "*" not in record["sequence_aa"]
            fresh = compiled.simulator.run(seed=args.seed, strict=True)
            trace = compiled.simulator.trace_file_from(fresh, seed=args.seed)
            replay = compiled.simulator.replay_from_trace_file(trace, strict=True)
            assert fresh.airr_record(compiled.refdata) == replay.airr_record(compiled.refdata)
            item = {"config": name, "productive": productive, "records": len(result), "errors": 0, "replay_equal": True}
            results.append(item)
            print(json.dumps(item), flush=True)
    artifact = {"seed": args.seed, "records": sum(r["records"] for r in results),
                "errors": 0, "profiles": results, "extension_path": _engine.__file__,
                "guarded_rejections": rejections,
                "productive_indel_count": 2 if args.productive_indels else 0,
                "reference_policy": "allow curatable reference annotations; productive support still requires valid retained anchors and coding sequence",
                "mutation_profiles": "IGH/IGK/IGL use SHM; unconstrained TCR profiles use sequencing errors; constrained TCR profiles have no SHM",
                "extension_sha256": hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
                "scope": "authored outcome distributions; metadata and replay correctness, not empirical frequency validation"}
    args.output.write_text(json.dumps(artifact, indent=2) + "\n")


if __name__ == "__main__":
    main()
