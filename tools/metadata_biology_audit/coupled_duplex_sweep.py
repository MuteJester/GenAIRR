"""Exhaustive dinucleotide V/D/J controls for coupled duplex repair.

Expected products follow literal strand complementation and conservation, not
the simulation's own projection or classification helpers. No biochemical
frequencies are fitted or inferred by this check.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path

from GenAIRR import _engine
from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def joins(offset, order):
    scripts = {name: (DuplexJoin(name).open("left", offset=offset)
                     .open("right", offset=offset).synapse(right_start=2).ligate())
               for name in ("vd", "dj")}
    return json.dumps([scripts[name].to_dict() for name in order])


def origins(result, references):
    bases = [base for column in result["columns"] for base in (column["top"], column["bottom"])]
    assert result["removed"] == []
    assert len(bases) == len(result["ledger"]) == 2*sum(map(len, references.values()))
    identities = {}
    for base in bases:
        origin = base["original"]
        key = (origin["segment"], origin["strand"], origin["reference_position"])
        assert key not in identities
        identities[key] = base["id"]
        reference = references[origin["segment"]][origin["reference_position"]]
        expected = reference if origin["strand"] == "top" else reference.translate(COMPLEMENT)
        assert base["base"] == expected
        assert base["kind"] == "original_strand" and base["template_id"] is None
    assert len(set(identities.values())) == len(bases)
    expected_origins = {(segment, strand, position) for segment, sequence in references.items()
                        for strand in ("top", "bottom") for position in range(len(sequence))}
    assert set(identities) == expected_origins
    assert {base["id"] for base in result["ledger"]} == set(identities.values())
    return identities


def run():
    substrates = ["".join(bases) for bases in product("ACGT", repeat=2)]
    scripts = {(offset, order): joins(offset, order) for offset in (0, 1)
               for order in (("vd", "dj"), ("dj", "vd"))}
    accepted, rejected = Counter(), Counter()
    for v, d, j in product(substrates, repeat=3):
        identity_map = None
        references = {"v": v, "d": d, "j": j}
        for inverted in (False, True):
            placed_d = reverse_complement(d) if inverted else d
            expected = v + placed_d + j
            for offset in (0, 1):
                compatible = not offset or (v[-1].translate(COMPLEMENT) == placed_d[0]
                                             and placed_d[-1].translate(COMPLEMENT) == j[0])
                previous = None
                for order in (("vd", "dj"), ("dj", "vd")):
                    script = CoupledDuplexRepair(v, d, j, d_inverted=inverted,
                                                joins_json=scripts[(offset, order)])
                    try:
                        result = script.run()
                    except ValueError as error:
                        assert not compatible, (v, d, j, inverted, offset, order, str(error))
                        assert "incompatible" in str(error), str(error)
                        rejected[f"offset_{offset}"] += 1
                        continue
                    assert compatible, (v, d, j, inverted, offset, order)
                    assert result["sequence"] == expected
                    assert result["complementary_sequence"] == reverse_complement(expected)
                    identities = origins(result, references)
                    if identity_map is None:
                        identity_map = identities
                    assert identities == identity_map, "canonical original identities changed with order/orientation"
                    projection = (result["columns"], result["removed"])
                    if previous is not None:
                        assert projection == previous, "independent compatible joins did not commute"
                    previous = projection
                    accepted[f"offset_{offset}"] += 1
    assert dict(accepted) == {"offset_0": 16384, "offset_1": 1024}, accepted
    assert dict(rejected) == {"offset_1": 15360}, rejected
    return {"schema": "coupled_duplex_sweep_v1", "model": "experimental_coupled_duplex_repair_v1",
            "extension_sha256": hashlib.sha256(Path(_engine.__file__).read_bytes()).hexdigest(),
            "substrate_length": 2, "substrate_triples": 4096, "d_orientations": 2,
            "join_orders": [["vd", "dj"], ["dj", "vd"]], "opening_offsets": [0, 1],
            "accepted": dict(accepted), "incompatible_rejected": dict(rejected),
            "errors": 0,
            "checks": ["literal product and complementary strand", "canonical source nucleotide identities",
                       "single allocation of each D original", "no creation during opening",
                       "stable canonical IDs under inversion", "commutation for independent compatible fixtures"],
            "scope": "deterministic short-substrate controls, not empirical repair frequencies"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
