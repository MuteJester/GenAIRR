"""Exhaustive short-duplex annotation oracle, without a sequence aligner."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path

from GenAIRR import RefDataConfig, _engine
from GenAIRR.experimental import CoupledDuplexRepair
from .coupled_duplex_sweep import COMPLEMENT, joins, origins, reverse_complement


def run():
    substrates = ["".join(bases) for bases in product("ACGT", repeat=2)]
    scripts = {(offset, order): joins(offset, order) for offset in (0, 1)
               for order in (("vd", "dj"), ("dj", "vd"))}
    accepted, rejected = Counter(), Counter()
    for v, d, j in product(substrates, repeat=3):
        r = RefDataConfig.vdj()
        r.add_v_allele("V*01", "V", v.encode())
        r.add_d_allele("D*01", "D", d.encode())
        r.add_j_allele("J*01", "J", j.encode())
        for inverted in (False, True):
            placed_d = reverse_complement(d) if inverted else d
            expected = v + placed_d + j
            for offset in (0, 1):
                compatible = not offset or (v[-1].translate(COMPLEMENT) == placed_d[0]
                                             and placed_d[-1].translate(COMPLEMENT) == j[0])
                for order in (("vd", "dj"), ("dj", "vd")):
                    script = CoupledDuplexRepair(v, d, j, d_inverted=inverted,
                                                joins_json=scripts[(offset, order)])
                    try:
                        p = script.project(r, v_call="V*01", d_call="D*01", j_call="J*01",
                                           reference_ends_are_coding_ends=True)
                    except ValueError as error:
                        assert not compatible and "incompatible" in str(error), (v, d, j, inverted, offset, order, str(error))
                        rejected[f"offset_{offset}"] += 1
                        continue
                    assert compatible
                    rec = p["record"]
                    truth = rec["duplex_repair"]
                    assert rec["sequence"] == rec["sequence_alignment"] == rec["germline_alignment"] == expected
                    assert rec["d_inverted"] == inverted
                    assert rec["np1"] == rec["np2"] == ""
                    assert rec["np1_length"] == rec["np2_length"] == 0
                    assert p["capabilities"]["simulation_seed"] is True
                    origins(truth["product"], {"v": v, "d": d, "j": j})
                    for index, name in enumerate(("v", "d", "j")):
                        path = p["paths"][name]
                        assert (path["query_start"], path["query_end"], path["reference_start"], path["reference_end"], path["cigar"]) == (2*index, 2*index+2, 0, 2, "2M")
                        assert rec[f"{name}_trim_5"] == rec[f"{name}_trim_3"] == 0
                        for point in path["body"]:
                            g = point["reference_position"]
                            expected_query = 2*index + (1-g if name == "d" and inverted else g)
                            assert point["query_position"] == expected_query
                            assert len(point["support_ids"]) in (1, 2)
                            if not offset:
                                assert len(point["support_ids"]) == 2
                    assert truth["overlaps"] == []
                    assert sum(len(p["paths"][s]["body"]) for s in ("v", "d", "j")) == 6
                    assert truth["counts"]["physical"]["hairpin_transfer_events"] == 4*offset
                    assert truth["counts"]["physical"]["untemplated_additions"] == 0
                    assert rec["n_mutations"] == rec["n_indels"] == 0
                    accepted[f"offset_{offset}"] += 1
    assert dict(accepted) == {"offset_0": 16384, "offset_1": 1024}
    assert dict(rejected) == {"offset_1": 15360}
    return {"schema": "duplex_projection_sweep_v1", "projection_model": "duplex_airr_projection_v1",
            "extension_sha256": hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
            "cases": 32768, "accepted": dict(accepted), "incompatible_rejected": dict(rejected),
            "errors": 0,
            "checks": ["literal sequence and full body alignment", "canonical inverted D paths",
                       "unique body claims supported by either strand", "hairpin transfers separate from body placement",
                       "empty NP intervals and zero net trims", "paired physical origin conservation",
                       "no repair operations counted as SHM or downstream indels"],
            "scope": "deterministic dinucleotide controls, not empirical frequencies"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
