"""Recover authored retained outcomes from independent native simulations.

All support is feasible for the single literal allele per segment, so there is
no allele-selection or productive-conditioning correction hidden in this check.
This checks outcome estimation; it does not calibrate biochemical enzyme rates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import GenAIRR as ga
from GenAIRR import _engine


def run(training_size=2000, held_out_size=500):
    reference = ga.RefDataConfig.vdj()
    reference.add_v_allele("v*01", "v", b"ATGAAATGTAAA", anchor=6)
    reference.add_d_allele("d*01", "d", b"AACCGGTT")
    reference.add_j_allele("j*01", "j", b"AAATGGAAA", anchor=3)
    identity = reference.content_hash()
    authored = ga.CodingEndOutcomeModel(
        v=((0, -1, 1), (0, 0, 3), (0, 1, 1)),
        d=((-1, 1, 3), (1, -1, 1)),
        j=((0, 0, 3), (1, 0, 2)),
        reference_hashes=(identity,),
    )
    compiled = (ga.Experiment.on(reference)
                .recombine(np1_lengths=[(0, 1)], np2_lengths=[(0, 1)])
                .process_coding_ends(model=authored, reference_ends_are_coding_ends=True)
                .invert_d(prob=.5).compile())
    seed_start = 2026102000

    def rows(start, count):
        for seed in range(start, start + count):
            yield compiled.simulator.run(seed=seed, strict=True).airr_record(compiled.refdata)

    fitted = ga.CodingEndOutcomeModel.fit_generating_records(
        rows(seed_start, training_size), chain_type="vdj", reference_hash=identity,
        source="independent native synthetic founders; one allele per segment; all outcomes feasible",
        training_sampling_law="unconditioned_generating_choices",
        independent_generating_events=True, support=authored,
    )
    comparisons = {}
    maximum_error = 0.0
    for segment in ("v", "d", "j"):
        source = getattr(authored, segment)
        total = math.fsum(w for _, _, w in source)
        recovered = {(a, b): w for a, b, w in getattr(fitted, segment)}
        observed = {(r["five"], r["three"]): r for r in fitted.provenance["counts"][segment.upper()]}
        comparisons[segment.upper()] = []
        for five, three, weight in source:
            target, estimate = weight / total, recovered.get((five, three), 0.)
            error = abs(target - estimate)
            maximum_error = max(maximum_error, error)
            comparisons[segment.upper()].append({
                "five": five, "three": three, "target": target, "estimate": estimate,
                "absolute_error": error, **observed[(five, three)],
            })
    held_out = list(rows(seed_start + training_size, held_out_size))
    fitted_score = fitted.score_generating_records(held_out)
    authored_score = authored.score_generating_records(held_out)
    # A broad predetermined recovery tolerance, not a claim of simultaneous
    # confidence coverage. Marginal intervals are reported without cherry-picking.
    assert maximum_error < .04, (maximum_error, comparisons)
    assert fitted_score["unseen_outcomes"] == 0
    assert abs(fitted_score["log_likelihood"] - authored_score["log_likelihood"]) / held_out_size < .03
    return {
        "schema": "coding_end_synthetic_recovery_v1",
        "reference_hash": identity,
        "extension_sha256": hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
        "training_seed_range": [seed_start, seed_start + training_size - 1],
        "held_out_seed_range": [seed_start + training_size, seed_start + training_size + held_out_size - 1],
        "training_records": training_size, "held_out_records": held_out_size,
        "sampling_law": "unconditioned_generating_choices",
        "maximum_absolute_error": maximum_error, "recovery_tolerance": .04,
        "outcomes": comparisons,
        "fitted_held_out_score": fitted_score, "authored_held_out_score": authored_score,
        "fit_provenance": fitted.provenance,
        "limitations": ["synthetic retained outcomes, not measured biochemical rates",
                        "single feasible allele per segment removes feasibility filtering",
                        "marginal intervals are not simultaneous confidence intervals",
                        "measured-data donor/clone validation remains pending"],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(report)
    else:
        print(report, end="")
