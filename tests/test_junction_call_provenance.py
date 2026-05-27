"""Junction-call provenance golden tests.

Companion to [docs/junction_call_audit.md](../docs/junction_call_audit.md).
Pins the derivation rules for the AIRR junction-related fields
(`junction`, `junction_aa`, `junction_length`, `junction_start`,
`junction_end`, `vj_in_frame`, `stop_codon`, `productive`) under
each evidence-changing mechanism.

Fixture conventions:

- Baseline VJ refdata: V=``AAACCCTGT`` (anchor at pos 6, anchor
  codon TGT = Cys), J=``TGGAAA`` (anchor at pos 0, anchor codon
  TGG = Trp). With NP1=3 (in-frame) the assembled sequence is
  18 bases and the junction is ``TGTACATGG`` (= ``CTW`` after
  translation, length 9 = 3 codons).
- VDJ refdata: V/J as above, D=``GGGCCC`` (6 bases). With
  NP1=NP2=3 the junction spans V_anchor → NP1(3) → D(6) →
  NP2(3) → J_anchor+3, total 18 bases / 6 codons.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge
from GenAIRR import CompiledExperiment
from GenAIRR.utilities.misc import reverse_complement


# ──────────────────────────────────────────────────────────────────
# Refdata factories
# ──────────────────────────────────────────────────────────────────


def _vj_basic() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vj()
    # V anchor codon = V[6:9] = TGT (Cys)
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    # J anchor codon = J[0:3] = TGG (Trp)
    cfg.add_j_allele("j1*01", "j1", b"TGGAAA", anchor=0)
    return cfg


def _vdj_basic() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCTGT", anchor=6)
    cfg.add_d_allele("d1*01", "d1", b"GGGCCC")
    cfg.add_j_allele("j1*01", "j1", b"TGGAAA", anchor=0)
    return cfg


def _baseline_vj() -> "ga.Experiment":
    return (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )


def _baseline_vdj() -> "ga.Experiment":
    return (
        ga.Experiment.on(_vdj_basic())
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )


def _record(exp, seed=0, **kwargs):
    return exp.run_records(n=1, seed=seed, **kwargs).records[0]


# ──────────────────────────────────────────────────────────────────
# 1. Field exposure (§1)
# ──────────────────────────────────────────────────────────────────


def test_junction_fields_present_on_airr_record() -> None:
    """Audit §1: enumerate the junction fields the AIRR projection
    exposes, plus the explicit non-exposure of `cdr3` / `cdr3_aa`."""
    rec = _record(_baseline_vj())
    for field in (
        "junction", "junction_aa", "junction_length",
        "junction_start", "junction_end",
        "vj_in_frame", "stop_codon", "productive",
    ):
        assert field in rec, f"AIRR record missing junction field {field!r}"
    # Non-exposed: junction IS the CDR3 window in IMGT convention.
    for field in ("cdr3", "cdr3_aa", "v_anchor_codon_pos", "j_anchor_codon_pos"):
        assert field not in rec, (
            f"unexpected field {field!r} appeared — audit §1 pinned it absent"
        )


# ──────────────────────────────────────────────────────────────────
# 2. Junction window boundaries (§2)
# ──────────────────────────────────────────────────────────────────


def test_vj_baseline_junction_window_spans_anchor_to_anchor_plus_3() -> None:
    """Audit §2: junction = [V_anchor_pool, J_anchor_pool + 3).

    Baseline fixture:
      sequence = AAACCCTGT (V) + NP1(3) + TGGAAA (J) = 18 bases.
      V_anchor_pool = 6 (V_region.start=0 + V_anchor=6 - V_trim_5=0).
      J_anchor_pool = 12 (J_region.start=12 + J_anchor=0 - J_trim_5=0).
      junction = sequence[6:15] = TGT + NP1 + TGG.
      junction_length = 9, three codons (CTW).
    """
    rec = _record(_baseline_vj())
    assert rec["junction_start"] == 6
    assert rec["junction_end"] == 15
    assert rec["junction_length"] == 9
    # First codon is the V anchor codon (TGT = Cys).
    assert rec["junction"][:3] == "TGT"
    # Last codon is the J anchor codon (TGG = Trp).
    assert rec["junction"][-3:] == "TGG"
    # junction equals the slice of sequence at those coordinates.
    assert rec["junction"] == rec["sequence"][rec["junction_start"]:rec["junction_end"]]


def test_vj_zero_np_junction_collapses_to_anchor_codons() -> None:
    """Audit §5.1: with NP1=0, the V and J anchor codons sit
    adjacent. Junction = V_anchor + J_anchor = TGT + TGG = TGTTGG.
    Length 6, two codons (CW). Productive."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(0, 1.0)])
        .trim(enabled=False)
    )
    rec = _record(exp)
    assert rec["junction"] == "TGTTGG"
    assert rec["junction_length"] == 6
    assert rec["junction_aa"] == "CW"
    assert rec["productive"] is True


def test_vdj_junction_spans_v_np1_d_np2_j() -> None:
    """Audit §5.1: VDJ junction crosses V → NP1 → D → NP2 → J.
    With NP1=NP2=3 and D=6: junction length = 3 (V_anchor_to_end)
    + 3 (NP1) + 6 (D) + 3 (NP2) + 3 (J_anchor) = 18 → 6 codons."""
    rec = _record(_baseline_vdj())
    assert rec["junction_length"] == 18
    assert rec["junction"][:3] == "TGT"   # V anchor codon
    assert rec["junction"][-3:] == "TGG"  # J anchor codon
    # D bases (GGGCCC) should appear somewhere in the middle.
    assert "GGGCCC" in rec["junction"]
    assert rec["vj_in_frame"] is True
    assert len(rec["junction_aa"]) == 6


# ──────────────────────────────────────────────────────────────────
# 3. Translation (§3)
# ──────────────────────────────────────────────────────────────────


def test_junction_translation_invariant_codon_to_aa() -> None:
    """Audit §3.2: in-frame junction translation is codon-aligned;
    len(junction)/3 == len(junction_aa). Each non-overlapping
    triplet maps to exactly one amino acid character."""
    rec = _record(_baseline_vj())
    junction = rec["junction"]
    junction_aa = rec["junction_aa"]
    assert rec["vj_in_frame"] is True
    assert len(junction) // 3 == len(junction_aa)
    # Codon-by-codon spot check against the known genetic code:
    # TGT=C, ACA=T, TGG=W on the baseline.
    expected_codons = {"TGT": "C", "ACA": "T", "TGG": "W"}
    for i in range(0, len(junction), 3):
        codon = junction[i:i+3].upper()
        if codon in expected_codons:
            assert junction_aa[i // 3] == expected_codons[codon], (
                f"codon {codon!r} should translate to "
                f"{expected_codons[codon]!r}, got {junction_aa[i // 3]!r}"
            )


def test_out_of_frame_junction_emits_empty_aa_and_false_stop() -> None:
    """Audit §3.1 + §6.1: when junction_length is not divisible by
    3, junction_aa is the empty string `""` (NOT None). stop_codon
    is also `False` (vacuously, not checked when OOF). productive
    is `False`.

    NP1=1 → junction_length = 3 + 1 + 3 = 7, OOF."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(1, 1.0)])
        .trim(enabled=False)
    )
    rec = _record(exp)
    assert rec["junction_length"] == 7
    assert rec["vj_in_frame"] is False
    assert rec["junction_aa"] == ""  # empty string, not None
    assert rec["stop_codon"] is False  # vacuously, OOF
    assert rec["productive"] is False


# ──────────────────────────────────────────────────────────────────
# 4. Productive triad (§4)
# ──────────────────────────────────────────────────────────────────


def test_productive_true_for_in_frame_no_stop_anchors_preserved() -> None:
    """Audit §4: productive = in_frame ∧ !stop ∧ anchors_preserved.
    Baseline fixture has all three: triad fires."""
    rec = _record(_baseline_vj())
    assert rec["vj_in_frame"] is True
    assert rec["stop_codon"] is False
    assert rec["productive"] is True


def test_in_frame_junction_with_stop_codon_is_not_productive() -> None:
    """Audit §4 row 2: in_frame=True, stop_codon=True → productive=False.
    Seed=0 with mutate(rate=0.5) lands a TGA inside the junction
    codon frame."""
    exp = (
        _baseline_vj().mutate(rate=0.5)
    )
    rec = _record(exp, seed=0)
    assert rec["vj_in_frame"] is True
    assert rec["stop_codon"] is True
    assert "*" in rec["junction_aa"]  # stop codon in AA translation
    assert rec["productive"] is False


def test_v_anchor_mutation_to_non_cys_breaks_productive_only() -> None:
    """Audit §4.2: in-frame + no-stop + V anchor amino acid
    changed → productive=False. Seed=2 with mutate(rate=0.5)
    produces junction_aa='HCW' (V anchor C → H = His)."""
    exp = _baseline_vj().mutate(rate=0.5)
    rec = _record(exp, seed=2)
    assert rec["vj_in_frame"] is True
    assert rec["stop_codon"] is False
    assert rec["junction_aa"][0] != "C"  # V anchor no longer Cys
    assert rec["productive"] is False


def test_j_anchor_mutation_to_non_trp_breaks_productive_only() -> None:
    """Audit §4.2 mirror: J anchor amino acid changed →
    productive=False. Seed=3 with mutate(rate=0.5) produces
    junction_aa='YEG' (J anchor W → G = Gly)."""
    exp = _baseline_vj().mutate(rate=0.5)
    rec = _record(exp, seed=3)
    assert rec["vj_in_frame"] is True
    assert rec["stop_codon"] is False
    assert rec["junction_aa"][-1] != "W"  # J anchor no longer Trp
    assert rec["productive"] is False


def test_stop_codon_is_false_when_out_of_frame() -> None:
    """Audit §4.1 / §6.3: stop_codon is `False` (not `None`) when
    vj_in_frame is False. The check is only meaningful in-frame;
    OOF returns the vacuous False."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(2, 1.0)])  # OOF
        .trim(enabled=False)
    )
    rec = _record(exp)
    assert rec["vj_in_frame"] is False
    assert rec["stop_codon"] is False
    assert rec["productive"] is False


# ──────────────────────────────────────────────────────────────────
# 5. Indels (§5.3)
# ──────────────────────────────────────────────────────────────────


def test_single_indel_in_middle_zone_breaks_in_frame() -> None:
    """Audit §5.3 row 2: an indel in [V_region.start, J_region.start)
    shifts J.region.start but not V.region.start → junction length
    ± 1 → out-of-frame.

    Seed=3 with polymerase_indels(count=1, insertion_prob=1.0)
    lands the insertion inside NP1 (between V and J), bumping
    junction_length from 9 to 10."""
    exp = _baseline_vj().polymerase_indels(count=1, insertion_prob=1.0)
    rec = _record(exp, seed=3)
    assert rec["junction_length"] == 10
    assert rec["vj_in_frame"] is False
    assert rec["productive"] is False


def test_indel_before_v_anchor_shifts_window_but_preserves_length() -> None:
    """Audit §5.3 row 1: an indel at site < V_region.start shifts
    both V and J region starts → junction window slides as a whole,
    length unchanged.

    Note: in our fixture V_region.start=0 so "site < V_region.start"
    is empty. But an insertion *inside* V before the anchor codon
    behaves similarly for the J side — V_anchor_pool shifts by +1
    AND J_anchor_pool shifts by +1, so junction length stays 9.

    Seed=1 with polymerase_indels(count=1, insertion_prob=1.0)
    yields v_cigar='5M1I4M' (insertion at V position 5, before the
    anchor codon at V[6:9]). junction_start shifts from 6 → 7,
    junction_end from 15 → 16, length stays 9."""
    exp = _baseline_vj().polymerase_indels(count=1, insertion_prob=1.0)
    rec = _record(exp, seed=1)
    assert rec["junction_length"] == 9
    assert rec["junction_start"] == 7  # shifted by +1
    assert rec["junction_end"] == 16  # shifted by +1
    assert rec["vj_in_frame"] is True  # length still in-frame


# ──────────────────────────────────────────────────────────────────
# 6. End-loss (§5.4)
# ──────────────────────────────────────────────────────────────────


def test_end_loss_5_prime_shifts_junction_coordinates_not_content() -> None:
    """Audit §5.4: 5' primer trim of N bases shifts junction
    coordinates left by N but doesn't change junction content
    (the anchor codons survive end-loss because they're not at
    the pool boundary).

    Baseline junction = TGTACATGG at [6, 15). After
    primer_trim_5prime(length=2): sequence is 16 bases, junction
    content unchanged but coordinates [4, 13)."""
    baseline_rec = _record(_baseline_vj())
    assert baseline_rec["junction_start"] == 6
    assert baseline_rec["junction_end"] == 15

    exp = _baseline_vj().primer_trim_5prime(length=2)
    rec = _record(exp)
    assert rec["junction"] == baseline_rec["junction"]  # content unchanged
    assert rec["junction_length"] == baseline_rec["junction_length"]
    assert rec["junction_start"] == 4  # shifted by -2
    assert rec["junction_end"] == 13   # shifted by -2
    assert rec["end_loss_5_length"] == 2
    assert rec["productive"] is True


# ──────────────────────────────────────────────────────────────────
# 7. Reverse-complement (§5.5)
# ──────────────────────────────────────────────────────────────────


def test_rev_comp_flips_junction_string_and_coordinates() -> None:
    """Audit §5.5: under rev-comp, the junction nucleotide string
    is reverse-complemented in place, and the (start, end)
    coordinates flip via new_start = seq_len - old_end,
    new_end = seq_len - old_start."""
    baseline = _record(_baseline_vj())
    rc = _record(_baseline_vj().random_strand_orientation(prob=1.0))

    assert baseline["rev_comp"] is False
    assert rc["rev_comp"] is True

    seq_len = len(baseline["sequence"])
    assert rc["junction"] == reverse_complement(baseline["junction"])
    assert rc["junction_start"] == seq_len - baseline["junction_end"]
    assert rc["junction_end"] == seq_len - baseline["junction_start"]
    # Length preserved.
    assert rc["junction_length"] == baseline["junction_length"]


def test_rev_comp_re_translates_junction_aa_from_rc_string() -> None:
    """Audit §5.5: junction_aa under rev-comp is the translation
    of the reverse-complemented junction string, NOT the reverse
    of the original junction_aa.

    Original junction TGTACATGG → CTW.
    rev-comp(TGTACATGG) = CCATGTACA → CCA(P) + TGT(C) + ACA(T) = PCT.
    So rc_record['junction_aa'] = 'PCT', not 'WTC'."""
    baseline = _record(_baseline_vj())
    rc = _record(_baseline_vj().random_strand_orientation(prob=1.0))
    assert baseline["junction_aa"] == "CTW"
    assert rc["junction_aa"] == "PCT"
    # And NOT the reverse string of the original:
    assert rc["junction_aa"] != baseline["junction_aa"][::-1]


# ──────────────────────────────────────────────────────────────────
# 8. Replay round-trip
# ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "configure,seed",
    [
        pytest.param(lambda e: e, 0, id="vj_baseline"),
        pytest.param(
            lambda e: e.mutate(rate=0.3), 0, id="shm"
        ),
        pytest.param(
            lambda e: e.polymerase_indels(count=1, insertion_prob=1.0),
            3,
            id="indel_in_middle",
        ),
        pytest.param(
            lambda e: e.primer_trim_5prime(length=2),
            0,
            id="end_loss",
        ),
        pytest.param(
            lambda e: e.random_strand_orientation(prob=1.0),
            0,
            id="rev_comp",
        ),
    ],
)
def test_replay_reproduces_all_junction_fields(configure, seed) -> None:
    """Trace replay reproduces every junction field exactly across
    every evidence-changing mechanism. If a future change leaks RNG
    state that the trace doesn't capture, this test diverges first."""
    compiled = configure(_baseline_vj()).compile()
    assert isinstance(compiled, CompiledExperiment)

    out_a = compiled.simulator.run(seed=seed)
    tf = compiled.simulator.trace_file_from(out_a, seed=seed)
    out_b = compiled.simulator.replay_from_trace_file(tf)

    from GenAIRR._airr_record import outcome_to_airr_record

    rec_a = outcome_to_airr_record(out_a, compiled.refdata, sequence_id="a")
    rec_b = outcome_to_airr_record(out_b, compiled.refdata, sequence_id="b")
    for field in [
        "junction", "junction_aa", "junction_length",
        "junction_start", "junction_end",
        "vj_in_frame", "stop_codon", "productive",
        "sequence", "sequence_aa", "rev_comp",
    ]:
        assert rec_a.get(field) == rec_b.get(field), (
            f"replay diverged on {field!r}: "
            f"{rec_a.get(field)!r} vs {rec_b.get(field)!r}"
        )


# ──────────────────────────────────────────────────────────────────
# 9. Pin-the-drift (§6)
# ──────────────────────────────────────────────────────────────────


def test_drift_pinned_junction_aa_is_empty_string_not_none_when_oof() -> None:
    """Audit §6.1: out-of-frame junction's junction_aa is the
    empty string `""`, not `None`. Pin the convention — a future
    fix to return `None` would flip this test, signalling the
    schema change."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(2, 1.0)])  # OOF
        .trim(enabled=False)
    )
    rec = _record(exp)
    assert rec["vj_in_frame"] is False
    assert rec["junction_aa"] == ""  # not None
    assert rec["junction_aa"] is not None


def test_drift_pinned_no_cdr3_aliases_on_airr_record() -> None:
    """Audit §6.2: no `cdr3` / `cdr3_aa` denormalised fields
    today — `junction` is the CDR3 window (IMGT convention).
    Pin absence; a future fix to add the aliases flips this."""
    rec = _record(_baseline_vj())
    for field in ("cdr3", "cdr3_aa"):
        assert field not in rec, (
            f"AIRR record now exposes {field!r} — flip this test "
            f"to assert cdr3 == junction[3:-3] when in-frame"
        )


def test_drift_pinned_stop_codon_is_false_not_none_when_oof() -> None:
    """Audit §6.3: out-of-frame junction's stop_codon is False
    (vacuously), not None. Pin the convention."""
    exp = (
        ga.Experiment.on(_vj_basic())
        .recombine(np1_lengths=[(2, 1.0)])
        .trim(enabled=False)
    )
    rec = _record(exp)
    assert rec["vj_in_frame"] is False
    assert rec["stop_codon"] is False
    assert rec["stop_codon"] is not None


def test_drift_pinned_no_productive_fail_reason_field() -> None:
    """Audit §6.4: no per-record reason field for why
    productive=False today. Users reconstruct from the triad
    (vj_in_frame, stop_codon, anchors_preserved). Pin absence."""
    exp = _baseline_vj().mutate(rate=0.5)
    rec = _record(exp, seed=0)  # stop-codon case
    assert rec["productive"] is False
    assert "productive_fail_reason" not in rec
    assert "productive_reason" not in rec
