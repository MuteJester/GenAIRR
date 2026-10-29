"""Independent sequence checks and explicit-anchor P-addition fixtures."""

from itertools import product
import re

import GenAIRR as ga
from GenAIRR.alleles import DAllele, JAllele, VAllele
from GenAIRR.reference_models import (
    EmpiricalDistributionSpec,
    NpBaseModelSpec,
    ReferenceEmpiricalModels,
)


_CODONS = dict(zip(
    ("".join(codon) for codon in product("TCAG", repeat=3)),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
))
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def translate(sequence):
    sequence = sequence.upper()
    return "".join(_CODONS.get(sequence[i:i + 3], "X") for i in range(0, len(sequence) - 2, 3))


def reverse_complement(sequence):
    return sequence.translate(_COMPLEMENT)[::-1]


def synthetic_cartridge(p_lengths, *, v="AAACCCTGT", v_anchor=6, d=None, np_base="C"):
    """Use public cartridge/allele authoring surfaces with known coding frames."""
    has_d = d is not None
    prefix = "IGH" if has_d else "IGK"
    v_allele = VAllele(f"{prefix}V-PTEST*01", v, len(v), anchor_override=v_anchor)
    j = "TGGAAA" if has_d else "TTTGGG"
    j_allele = JAllele(f"{prefix}J-PTEST*01", j, len(j), anchor_override=0)
    d_alleles = None
    if has_d:
        d_allele = DAllele("IGHD-PTEST*01", d, len(d))
        d_alleles = {d_allele.gene: [d_allele]}
    return ga.DataConfig(
        name="productive-p-fixture",
        metadata=(ga.HUMAN_IGH_OGRDB if has_d else ga.HUMAN_IGK_OGRDB).metadata,
        v_alleles={v_allele.gene: [v_allele]},
        d_alleles=d_alleles,
        j_alleles={j_allele.gene: [j_allele]},
        reference_models=ReferenceEmpiricalModels(
            p_nucleotide_lengths={
                end: EmpiricalDistributionSpec(pairs) for end, pairs in p_lengths.items()
            },
            np_bases={
                region: NpBaseModelSpec(kind="empirical_first_base", first_base={np_base: 1.0})
                for region in (("NP1", "NP2") if has_d else ("NP1",))
            } if np_base is not None else {},
        ),
    )


def assert_productive_record(record, cfg):
    """Check biology and alignment arithmetic without the engine validator."""
    sequence = record["sequence"]
    if record["rev_comp"]:
        sequence = reverse_complement(sequence)
    sequence = sequence.upper()
    v_name = record.get("truth_v_call", record["v_call"]).split(",")[0]
    v = next(allele for alleles in cfg.v_alleles.values() for allele in alleles if allele.name == v_name)
    # These fixtures have no 5-prime read loss. Structural insertions must
    # not move the retained V coding origin to the displaced Cys anchor.
    assert record["end_loss_5_length"] == 0
    phase = (v.anchor - record["v_trim_5"]) % 3
    assert record["sequence_aa"] == translate(sequence[phase:])
    assert "*" not in translate(sequence[phase:record["j_sequence_end"]])
    assert record["junction_aa"] == translate(record["junction"])
    assert record["junction_aa"].startswith("C")
    assert record["junction_aa"][-1] in "FW"
    assert record["productive"] is True
    assert record["stop_codon"] is False
    assert record["vj_in_frame"] is True

    alignment = record["sequence_alignment"].upper()
    assert alignment.replace("-", "") == sequence
    assert len(alignment) == len(record["germline_alignment"])
    has_d = bool(record["d_call"])
    intervals = [
        ("np1", record["v_sequence_end"], record["d_sequence_start"] if has_d else record["j_sequence_start"]),
    ]
    if has_d:
        intervals.append(("np2", record["d_sequence_end"], record["j_sequence_start"]))
    for name, start, end in intervals:
        assert record[name].upper() == sequence[start:end]
        assert record[f"{name}_length"] == end - start
    for segment in ("v", "d", "j") if has_d else ("v", "j"):
        runs = [(int(count), op) for count, op in re.findall(r"(\d+)([MID=X])", record[f"{segment}_cigar"])]
        assert sum(n for n, op in runs if op in "MI=X") == (
            record[f"{segment}_sequence_end"] - record[f"{segment}_sequence_start"]
        )
        assert sum(n for n, op in runs if op in "MD=X") == (
            record[f"{segment}_germline_end"] - record[f"{segment}_germline_start"]
        )
