"""T2-14: receptor revision must keep codon state and productivity
annotations aligned with the revised V template."""
from __future__ import annotations

from GenAIRR import Experiment, Productivity
from GenAIRR.ops import with_receptor_revision
from GenAIRR.utilities.mcp_helpers import validate_codon_rail_snapshot


_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _translate(seq: str) -> str:
    aa = []
    for i in range(0, len(seq) - 2, 3):
        aa.append(_CODON.get(seq[i:i + 3].upper(), "X"))
    return "".join(aa)


def _actual_status(rec):
    seq = rec["sequence"].upper()
    junction = rec["junction_nt"].upper()
    stop_codon = "*" in _translate(seq)
    vj_in_frame = (
        rec["junction_start"] % 3 == 0
        and rec["junction_end"] % 3 == 0
        and len(junction) % 3 == 0
    )
    junction_aa = _translate(junction) if vj_in_frame else ""
    productive = (
        bool(junction_aa)
        and not stop_codon
        and vj_in_frame
        and "*" not in junction_aa
        and "X" not in junction_aa
        and junction_aa[0] == "C"
        and junction_aa[-1] in {"F", "W"}
    )
    return productive, stop_codon, vj_in_frame


class TestCodonRailAfterRevision:
    """Receptor revision rewrites live V bases, so the codon rail must
    be rebuilt immediately afterward."""

    def test_post_receptor_revision_codon_rail_matches_rewritten_bases(self):
        sim = (Experiment.on("human_igh")
               .recombine(with_receptor_revision(prob=1.0, footprint=(5, 20)))
               .compile(seed=42))

        _, snapshots, _ = sim._sim.simulate_one_hooked(
            ["post_functionality", "post_receptor_rev"]
        )
        rails = {
            snap["hook"]: sim._sim.get_snapshot_codon_rail(i)
            for i, snap in enumerate(snapshots)
        }

        before_v = [c["bases"] for c in rails["post_functionality"]
                    if c["seg"] == "V"]
        after_v = [c["bases"] for c in rails["post_receptor_rev"]
                   if c["seg"] == "V"]

        assert before_v, "Expected at least one V-overlapping codon before revision"
        assert after_v, "Expected at least one V-overlapping codon after revision"
        assert before_v != after_v, (
            "Receptor revision should change at least one V-overlapping "
            "codon; otherwise this test is not exercising the rail "
            "rebuild path"
        )

        validation = validate_codon_rail_snapshot(rails["post_receptor_rev"])
        assert validation["valid"], (
            f"post_receptor_rev codon rail is stale or inconsistent: "
            f"{validation['issues'][:3]}"
        )


class TestFinalProductivityAnnotations:
    """Final AIRR productivity fields must agree with the revised
    emitted sequence after receptor revision."""

    def test_final_productivity_matches_revised_sequence(self):
        seeds = [45, 51, 62]

        for seed in seeds:
            rec = (Experiment.on("human_igh")
                   .recombine(with_receptor_revision(prob=1.0, footprint=(5, 20)))
                   .run(n=1, seed=seed))[0]

            assert rec["receptor_revised"] is True, (
                f"Seed {seed} should exercise the revision path"
            )

            productive, stop_codon, vj_in_frame = _actual_status(rec)
            assert rec["productive"] == productive, (
                f"Seed {seed} final productive annotation disagrees with "
                f"the revised sequence"
            )
            assert rec["stop_codon"] == stop_codon, (
                f"Seed {seed} final stop_codon annotation disagrees with "
                f"the revised sequence"
            )
            assert rec["vj_in_frame"] == vj_in_frame, (
                f"Seed {seed} final vj_in_frame annotation disagrees with "
                f"the revised sequence"
            )


class TestProductiveGuardedRevision:
    """When the productive contract is active, an unsafe receptor
    revision candidate should be skipped locally."""

    def test_productive_only_skips_unsafe_revision_seed1(self):
        sim = (Experiment.on("human_igh")
               .recombine(with_receptor_revision(prob=1.0, footprint=(1, 1)))
               .compile(seed=1, productivity=Productivity.PRODUCTIVE_ONLY))
        sim._sim.set_param("max_productive_attempts", 200)

        rec, snapshots, _ = sim._sim.simulate_one_hooked(
            ["post_functionality", "post_receptor_rev"]
        )
        by_hook = {snap["hook"]: snap["nodes"] for snap in snapshots}

        before_v = "".join(n["cur"] for n in by_hook["post_functionality"]
                            if n["seg"] == "V")
        after_v = "".join(n["cur"] for n in by_hook["post_receptor_rev"]
                           if n["seg"] == "V")

        assert rec["productive"] is True
        assert rec["receptor_revised"] is False, (
            "Seed 1 should exercise the productive guard skip path for "
            "an unsafe receptor revision candidate"
        )
        assert before_v
        assert before_v == after_v, (
            "V bases changed even though productive-only receptor "
            "revision should have been skipped"
        )
