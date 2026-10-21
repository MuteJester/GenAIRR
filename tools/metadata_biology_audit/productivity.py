"""Independent biological checks for productivity, stops and coding frame.

These probes report observed defects without changing the simulator. Coordinates
in observations are zero-based, half-open. They use the public Python API and an
independent standard genetic-code table, not GenAIRR's translation helper.
"""
from itertools import product

import GenAIRR as ga

_CODE = dict(zip(
    ("".join(codon) for codon in product("TCAG", repeat=3)),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
))


def _translate(sequence, frame=0):
    sequence = sequence.upper()[frame:]
    return "".join(
        _CODE.get(sequence[i:i + 3], "X")
        for i in range(0, len(sequence) - 2, 3)
    )


def _stops(sequence, frame=0):
    sequence = sequence.upper()
    return [
        {"start": i, "end": i + 3, "codon": sequence[i:i + 3]}
        for i in range(frame, len(sequence) - 2, 3)
        if _CODE.get(sequence[i:i + 3]) == "*"
    ]


def _run(exp, n, row):
    result = exp.run_records(n=n, seed=42)
    return result[row], result.validate_records(exp.refdata).summary()


def _fields(record):
    return {
        key: record[key] for key in (
            "v_call", "j_call", "productive", "stop_codon", "vj_in_frame",
            "v_sequence_start", "v_germline_start", "v_cigar", "junction_start",
            "junction_end", "junction", "junction_aa", "sequence_aa",
        )
    }


def run_checks():
    """Return machine-readable findings; reproduced=True means the bug persists."""
    findings = []

    exp = (ga.Experiment.on("human_igh").recombine().productive_only()
           .mutate(rate=0.03))
    record, issues = _run(exp, 1, 0)
    all_stops = _stops(record["sequence"])
    junction_stops = _stops(record["junction"])
    findings.append({
        "id": "productivity_outside_junction_shm",
        "title": "SHM creates an upstream stop but metadata still says productive",
        "expected": {
            "productive": False, "stop_codon": True,
            "reason": "A translated stop in the aligned V coding region is "
                      "outside the junction but still interrupts the receptor.",
        },
        "observed": {
            **_fields(record), "coding_stops": all_stops,
            "junction_stops": junction_stops,
            "sequence_translation_matches_independent_table":
                record["sequence_aa"] == _translate(record["sequence"]),
        },
        "reproduced": bool(
            record["productive"] is True and record["stop_codon"] is False
            and all_stops and not junction_stops
        ),
        "validator_issues": issues,
        "reproduction": "Experiment.on('human_igh').recombine().productive_only()"
                        ".mutate(rate=0.03).run_records(n=1, seed=42)[0]",
        "source": "engine_rs/src/airr_record/builder.rs:464",
    })

    exp = ga.Experiment.on("human_igh").recombine().productive_only()
    record, issues = _run(exp, 120, 119)
    allele = next(
        exp.refdata.v_allele(i)
        for i in range(exp.refdata.v_pool_size())
        if exp.refdata.v_allele(i).name == "IGHVF9-G32*02"
    )
    refseq = allele.seq().decode()
    findings.append({
        "id": "productivity_outside_junction_germline",
        "title": "A germline stop also escapes clean productive-only generation",
        "expected": {
            "productive": False, "stop_codon": True,
            "reason": "The retained reference TAG at [282,285) is an in-frame "
                      "V-region stop immediately before the Cys anchor.",
        },
        "observed": {
            **_fields(record), "coding_stops": _stops(record["sequence"]),
            "reference_allele": allele.name,
            "reference_functional_status": allele.functional_status,
            "reference_anchor": allele.anchor,
            "reference_codon_282_285": refseq[282:285],
            "record_codon_282_285": record["sequence"][282:285],
        },
        "reproduced": bool(
            record["productive"] is True and record["stop_codon"] is False
            and record["sequence"][282:285].upper()
            == refseq[282:285].upper() == "TAG"
        ),
        "validator_issues": issues,
        "reproduction": "Experiment.on('human_igh').recombine().productive_only()"
                        ".run_records(n=120, seed=42)[119]",
        "source": "engine_rs/src/airr_record/builder.rs:464",
    })

    exp = ga.Experiment.on("human_igh").recombine()
    record, issues = _run(exp, 2, 1)
    findings.append({
        "id": "stop_codon_out_of_frame",
        "title": "An out-of-frame junction forces stop_codon=False despite stops",
        "expected": {
            "stop_codon": True,
            "reason": "Loss of V/J frame agreement does not remove the stop "
                      "codons encountered when translating from the V frame.",
        },
        "observed": {
            **_fields(record), "coding_stops": _stops(record["sequence"]),
            "independent_junction_translation": _translate(record["junction"]),
        },
        "reproduced": bool(
            record["vj_in_frame"] is False and record["stop_codon"] is False
            and "*" in _translate(record["junction"])
        ),
        "validator_issues": issues,
        "reproduction": "Experiment.on('human_igh').recombine()"
                        ".run_records(n=2, seed=42)[1]",
        "source": "engine_rs/src/airr_record/builder.rs:484",
    })

    exp = (ga.Experiment.on("human_igh").recombine()
           .polymerase_indels(count=1, insertion_prob=1))
    record, issues = _run(exp, 20, 19)
    true_translation = _translate(record["sequence"])
    cfg = ga.RefDataConfig.vj()
    cfg.add_v_allele("v_test*01", "v_test", b"ATGAAATGTAAA", anchor=6)
    cfg.add_j_allele("j_test*01", "j_test", b"AAATGGGGG", anchor=3)
    minimal_exp = (ga.Experiment.on(cfg).recombine(np1_lengths=[(0, 1)])
                   .trim(enabled=False)
                   .polymerase_indels(count=1, insertion_prob=1))
    minimal, minimal_issues = _run(minimal_exp, 2, 1)
    findings.append({
        "id": "sequence_aa_upstream_frameshift",
        "title": "An upstream insertion resets translation to the junction frame",
        "expected": {
            "productive": False,
            "sequence_aa": true_translation,
            "reason": "The complete retained V starts in frame zero. An internal "
                      "one-base insertion shifts downstream translation; moving "
                      "the Cys anchor does not move the original coding origin.",
        },
        "observed": {
            **_fields(record), "sequence": record["sequence"],
            "independent_v_origin_translation": true_translation,
            "independent_v_origin_stops": _stops(record["sequence"]),
            "matches_wrong_frame_one":
                record["sequence_aa"] == _translate(record["sequence"], 1),
        },
        "reproduced": bool(
            record["v_sequence_start"] == 0
            and record["v_germline_start"] == 0
            and record["v_cigar"] == "136M1I163M"
            and record["sequence_aa"] != true_translation
            and record["sequence_aa"] == _translate(record["sequence"], 1)
            and record["productive"] is True
        ),
        "validator_issues": issues,
        "reproduction": "Experiment.on('human_igh').recombine()"
                        ".polymerase_indels(count=1, insertion_prob=1)"
                        ".run_records(n=20, seed=42)[19]",
        "minimal_example": {
            "v_reference": "ATGAAATGTAAA", "v_anchor": 6,
            "j_reference": "AAATGGGGG", "j_anchor": 3,
            "np_length": 0, "trim": False, "seed": 42, "n": 2, "row": 1,
            "sequence": minimal["sequence"],
            "reported_sequence_aa": minimal["sequence_aa"],
            "expected_sequence_aa": _translate(minimal["sequence"]),
            "productive": minimal["productive"],
            "v_cigar": minimal["v_cigar"],
            "validator_issues": minimal_issues,
        },
        "source": "engine_rs/src/airr_record/builder.rs:501",
    })

    for name, exp, finding_index in (
        ("clean", ga.Experiment.on("human_igh").recombine().productive_only(), 1),
        ("shm_3_percent", ga.Experiment.on("human_igh").recombine()
         .productive_only().mutate(rate=0.03), 0),
    ):
        result = exp.run_records(n=1000, seed=42)
        failures = [
            i for i, r in enumerate(result)
            if r["productive"] is True and "*" in _translate(r["sequence"])
        ]
        findings[finding_index]["sample"] = {
            "name": name, "n": 1000, "seed": 42,
            "productive_records_with_coding_stops": len(failures),
            "first_five_row_indexes": failures[:5],
            "note": "Seed-specific demonstration, not a population estimate.",
        }
    return findings
