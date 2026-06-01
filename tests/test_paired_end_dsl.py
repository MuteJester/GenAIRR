"""End-to-end tests for `Experiment.paired_end(...)` — Slice D.

Pins the user-facing DSL surface that exposes paired-end / read-
layout sampling as a fluent step:

- Both VDJ and VJ chains accepted.
- Fixed-value `r1_length` / `r2_length` / `insert_size` emit three
  trace records in the documented order
  (`r1_length`, `r2_length`, `insert_size`) at the END of the
  trace, after every biology / corruption record.
- `r2_length=None` mirrors `r1_length`.
- R1 / R2 sequences match the Slice B projection rules.
- Composes with mutation / end-loss / rev-comp without changing
  their outputs.
- Replay round-trip preserves the eight AIRR paired-end fields.
- Duplicate call rejects; bad geometry rejects at the DSL
  boundary for fixed values.
- DataFrame / CSV exports carry the eight populated fields.

No two-row FASTQ, no per-base quality scores, no read IDs — Slice E
is the release-consolidation slice.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Deterministic fixtures (matching the contract / schema test
# fixtures so future Slice E can reuse the same harness).
# ──────────────────────────────────────────────────────────────────


def _vdj_refdata() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vdj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_d_allele("d1*01", "d1", b"ACGTTA")
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def _vj_refdata() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"TGTAAACCC", anchor=0)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def _baseline(*, paired: bool, r1: int = 5, r2=None, insert: int = 12) -> "ga.Experiment":
    exp = (
        ga.Experiment.on(_vdj_refdata())
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )
    if paired:
        kwargs = {"r1_length": r1, "insert_size": insert}
        if r2 is not None:
            kwargs["r2_length"] = r2
        exp = exp.paired_end(**kwargs)
    return exp


# ──────────────────────────────────────────────────────────────────
# 1. Fixed values record the three trace records in order
# ──────────────────────────────────────────────────────────────────


def test_paired_end_fixed_values_emit_three_trace_records_in_order() -> None:
    """The three `paired_end.*` records appear in the documented
    `r1_length` → `r2_length` → `insert_size` order, with the
    sampled values matching the fixed inputs."""
    exp = _baseline(paired=True, r1=5, r2=5, insert=12)
    addrs_and_values = [
        (r.address, r.value)
        for r in exp.run(n=1, seed=0)[0].trace().choices()
        if r.address.startswith("paired_end.")
    ]
    assert [a for a, _ in addrs_and_values] == [
        "paired_end.r1_length",
        "paired_end.r2_length",
        "paired_end.insert_size",
    ], (
        f"paired_end.* trace records out of order: {addrs_and_values}"
    )
    assert [v for _, v in addrs_and_values] == [5, 5, 12]


# ──────────────────────────────────────────────────────────────────
# 2. r2_length=None mirrors r1_length
# ──────────────────────────────────────────────────────────────────


def test_r2_length_default_mirrors_r1_length() -> None:
    """When the user omits ``r2_length``, the recorded value
    equals ``r1_length``."""
    exp = _baseline(paired=True, r1=5, r2=None, insert=12)
    trace = exp.run(n=1, seed=0)[0].trace()
    r1_rec = next(r for r in trace.choices() if r.address == "paired_end.r1_length")
    r2_rec = next(r for r in trace.choices() if r.address == "paired_end.r2_length")
    assert r2_rec.value == r1_rec.value == 5


# ──────────────────────────────────────────────────────────────────
# 3. R1/R2 sequences match the Slice B projection rules
# ──────────────────────────────────────────────────────────────────


def test_r1_and_r2_sequences_match_projection_kernel_rules() -> None:
    """``r1_sequence == sequence[r1_start:r1_end]`` (forward) and
    ``r2_sequence == reverse_complement(sequence[r2_start:r2_end])``.
    Pinned end-to-end through the DSL so a Slice C refactor that
    silently swapped the R2 transform would surface here too."""
    rec = _baseline(paired=True, r1=5, r2=5, insert=12).run_records(n=1, seed=0).records[0]
    seq = rec["sequence"]
    assert rec["read_layout"] == "paired_end"
    assert rec["r1_sequence"] == seq[rec["r1_start"]:rec["r1_end"]]
    # R2 is the reverse-complement of the 3'-end window.
    def _py_revcomp(s: str) -> str:
        table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return s.translate(table)[::-1]
    assert rec["r2_sequence"] == _py_revcomp(seq[rec["r2_start"]:rec["r2_end"]])
    assert rec["insert_size"] == rec["r2_end"]


# ──────────────────────────────────────────────────────────────────
# 4. Works on VJ and VDJ
# ──────────────────────────────────────────────────────────────────


def test_paired_end_works_on_vdj_chain() -> None:
    """VDJ chains accept paired_end and emit the eight populated
    fields end-to-end."""
    rec = _baseline(paired=True, r1=5, r2=5, insert=12).run_records(n=1, seed=0).records[0]
    assert rec["read_layout"] == "paired_end"
    assert len(rec["r1_sequence"]) == 5
    assert len(rec["r2_sequence"]) == 5


def test_paired_end_works_on_vj_chain() -> None:
    """VJ chains are also supported — paired-end is sequencing-
    stage, not biology-stage."""
    exp = (
        ga.Experiment.on(_vj_refdata())
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .paired_end(r1_length=4, insert_size=10)
    )
    rec = exp.run_records(n=1, seed=0).records[0]
    assert rec["read_layout"] == "paired_end"
    assert len(rec["r1_sequence"]) == 4
    assert rec["insert_size"] == 10


# ──────────────────────────────────────────────────────────────────
# 5. Composes after end_loss + random_strand_orientation
# ──────────────────────────────────────────────────────────────────


def test_paired_end_composes_after_end_loss_and_strand_orientation() -> None:
    """End-loss shortens the molecule before paired-end carves
    windows; strand orientation flips the molecule. R1/R2 must
    still index into the projected (possibly rev-comped, possibly
    end-loss-shortened) `rec.sequence`."""
    exp = (
        ga.Experiment.on(_vdj_refdata())
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)], np2_lengths=[(3, 1.0)])
        .trim(enabled=False)
        .end_loss_5prime(length=[(2, 1.0)])
        .random_strand_orientation(prob=1.0)
        .paired_end(r1_length=4, insert_size=10)
    )
    rec = exp.run_records(n=1, seed=0).records[0]
    assert rec["read_layout"] == "paired_end"
    assert rec["end_loss_5_length"] == 2
    assert rec["rev_comp"] is True
    seq = rec["sequence"]
    # R1 = forward window into the post-rev-comp, post-end-loss seq.
    assert rec["r1_sequence"] == seq[rec["r1_start"]:rec["r1_end"]]


def test_paired_end_trace_records_land_after_biology_and_corruption() -> None:
    """Trace-order pin (per design doc §6 / the
    `_extract_paired_end_step` rationale): the three
    `paired_end.*` records appear AFTER every biology /
    corruption choice. We use a mutation step's count record as
    the proxy for "the last biology/corruption choice."""
    exp = (
        _baseline(paired=True, r1=5, r2=5, insert=12)
        .mutate(rate=0.01)
        .end_loss_5prime(length=[(1, 1.0)])
        .random_strand_orientation(prob=1.0)
    )
    addrs = [r.address for r in exp.run(n=1, seed=0)[0].trace().choices()]
    paired_end_indices = [
        i for i, a in enumerate(addrs) if a.startswith("paired_end.")
    ]
    assert paired_end_indices == [
        addrs.index("paired_end.r1_length"),
        addrs.index("paired_end.r2_length"),
        addrs.index("paired_end.insert_size"),
    ]
    # Every biology/corruption record sits before the first
    # paired_end record.
    mut_idx = addrs.index("mutate.s5f.count")
    end_loss_idx = addrs.index("corrupt.end_loss.5")
    rev_comp_idx = addrs.index("corrupt.rev_comp.applied")
    first_paired_end = paired_end_indices[0]
    assert mut_idx < first_paired_end
    assert end_loss_idx < first_paired_end
    assert rev_comp_idx < first_paired_end
    # And the paired-end records are contiguous + at the end.
    last_paired_end = paired_end_indices[-1]
    assert last_paired_end == len(addrs) - 1, (
        f"paired_end records should be the last in the trace; "
        f"got addrs[-1]={addrs[-1]!r}"
    )


def test_paired_end_composes_with_mutation_and_pcr_quality_errors() -> None:
    """Smoke pin that the DSL composes with mutation + PCR + quality
    errors without surfacing a paired-end-specific issue."""
    refdata = ga.Experiment.on("human_igh").allow_curatable_refdata().refdata
    exp = (
        ga.Experiment.on("human_igh")
        .allow_curatable_refdata()
        .recombine()
        .mutate(rate=0.01)
        .pcr_amplify(rate=1e-4)
        .sequencing_errors(rate=1e-4)
        .paired_end(r1_length=80, insert_size=200)
    )
    result = exp.run_records(n=3, seed=0)
    report = result.validate_records(refdata)
    # No paired-end issues. Other issues (D-tie under inversion,
    # etc.) are out of scope for this Slice D smoke test.
    paired_end_kinds = {
        "PairedEndFieldWithoutLayout",
        "ReadWindowOutOfBounds",
        "ReadSequenceMismatch",
        "ReadInsertSizeMismatch",
        "ReadLayoutMismatch",
    }
    bad_kinds = {
        issue["kind"]
        for failure in report.failures
        for issue in failure["issues"]
    }
    offenders = bad_kinds & paired_end_kinds
    assert not offenders, (
        f"paired-end issues surfaced under realistic stack: {offenders}"
    )


# ──────────────────────────────────────────────────────────────────
# 6. Replay round-trip preserves the eight AIRR fields
# ──────────────────────────────────────────────────────────────────


def test_paired_end_trace_replays_bit_for_bit() -> None:
    """A trace recorded by ``paired_end(r1=5, r2=5, insert=12)``
    replays through ``rerun_from_trace_file`` and reproduces all
    eight paired-end AIRR fields."""
    from GenAIRR._airr_record import outcome_to_airr_record

    exp = _baseline(paired=True, r1=5, r2=5, insert=12)
    compiled = exp.compile()
    fresh = compiled.simulator.run(seed=0)
    tf = compiled.simulator.trace_file_from(fresh, seed=0)
    replayed = compiled.simulator.rerun_from_trace_file(tf)

    fresh_rec = outcome_to_airr_record(fresh, exp._refdata, sequence_id="fresh")
    replayed_rec = outcome_to_airr_record(replayed, exp._refdata, sequence_id="replay")

    for field in (
        "read_layout",
        "r1_sequence",
        "r2_sequence",
        "r1_start",
        "r1_end",
        "r2_start",
        "r2_end",
        "insert_size",
    ):
        assert fresh_rec[field] == replayed_rec[field], (
            f"paired-end field {field!r} desynced under replay: "
            f"{fresh_rec[field]!r} → {replayed_rec[field]!r}"
        )


# ──────────────────────────────────────────────────────────────────
# 7. Duplicate / validation guards
# ──────────────────────────────────────────────────────────────────


def test_paired_end_called_twice_raises_value_error() -> None:
    """v1 picks a single paired-end layout per pipeline; the DSL
    rejects a second call so an over-eager builder fails loudly."""
    exp = _baseline(paired=False)
    exp.paired_end(r1_length=5, insert_size=12)
    with pytest.raises(ValueError, match="already configured"):
        exp.paired_end(r1_length=4, insert_size=10)


def test_paired_end_r1_length_must_be_positive() -> None:
    exp = _baseline(paired=False)
    with pytest.raises(ValueError, match="r1_length must be positive"):
        exp.paired_end(r1_length=0, insert_size=10)


def test_paired_end_r2_length_must_be_positive() -> None:
    exp = _baseline(paired=False)
    with pytest.raises(ValueError, match="r2_length must be positive"):
        exp.paired_end(r1_length=5, r2_length=0, insert_size=10)


def test_paired_end_insert_size_must_be_non_negative() -> None:
    """Insert-size sentinel is 0 (`no fragment geometry`), so the
    DSL accepts 0 here; the negative case is what
    `_normalize_count` rejects. Pin the actual rejection
    message comes from the normalize layer."""
    exp = _baseline(paired=False)
    with pytest.raises(ValueError, match="non-negative"):
        exp.paired_end(r1_length=5, insert_size=-1)


def test_paired_end_fixed_r1_larger_than_insert_rejects_at_dsl() -> None:
    """When all three are fixed values, the DSL pre-checks
    `r1_length > insert_size` and rejects rather than waiting
    for the engine."""
    exp = _baseline(paired=False)
    with pytest.raises(ValueError, match="r1_length .* > insert_size"):
        exp.paired_end(r1_length=15, insert_size=10)


def test_paired_end_fixed_r2_larger_than_insert_rejects_at_dsl() -> None:
    exp = _baseline(paired=False)
    with pytest.raises(ValueError, match="r2_length .* > insert_size"):
        exp.paired_end(r1_length=5, r2_length=15, insert_size=10)


def test_paired_end_rejects_unknown_kwargs() -> None:
    """Aliases like ``read_length=`` or ``r1=`` are not exposed
    yet; Python's TypeError surfaces the unknown kwarg."""
    exp = _baseline(paired=False)
    with pytest.raises(TypeError):
        exp.paired_end(read_length=150, insert_size=300)
    with pytest.raises(TypeError):
        exp.paired_end(r1=150, insert_size=300)


def test_paired_end_accepts_range_tuple_shape() -> None:
    """``(low, high)`` tuple syntax is normalised to a uniform
    integer distribution over the closed interval."""
    exp = _baseline(paired=False).paired_end(
        r1_length=(4, 6), r2_length=(4, 6), insert_size=(10, 14)
    )
    rec = exp.run_records(n=1, seed=0).records[0]
    assert rec["read_layout"] == "paired_end"
    assert 4 <= rec["r1_end"] - rec["r1_start"] <= 6
    assert 4 <= rec["r2_end"] - rec["r2_start"] <= 6
    assert 10 <= rec["insert_size"] <= 14


def test_paired_end_accepts_empirical_pairs_shape() -> None:
    """List-of-pairs syntax is normalised verbatim."""
    exp = _baseline(paired=False).paired_end(
        r1_length=[(5, 1.0)],
        r2_length=[(5, 1.0)],
        insert_size=[(12, 1.0)],
    )
    rec = exp.run_records(n=1, seed=0).records[0]
    assert rec["read_layout"] == "paired_end"
    assert rec["insert_size"] == 12


# ──────────────────────────────────────────────────────────────────
# 8. Baseline without `.paired_end()` keeps defaults
# ──────────────────────────────────────────────────────────────────


def test_baseline_without_paired_end_emits_no_paired_end_records() -> None:
    """The Slice A baseline-clean invariant survives Slice D. No
    `paired_end.*` records, all eight AIRR fields at defaults."""
    rec = _baseline(paired=False).run_records(n=1, seed=0).records[0]
    assert rec["read_layout"] == ""
    assert rec["r1_sequence"] == ""
    assert rec["r2_sequence"] == ""
    assert rec["r1_start"] is None
    assert rec["r1_end"] is None
    assert rec["r2_start"] is None
    assert rec["r2_end"] is None
    assert rec["insert_size"] == 0
    addrs = {r.address for r in _baseline(paired=False).run(n=1, seed=0)[0].trace().choices()}
    assert not any(a.startswith("paired_end.") for a in addrs)


# ──────────────────────────────────────────────────────────────────
# 9. DataFrame / CSV exports include the eight populated fields
# ──────────────────────────────────────────────────────────────────


def test_dataframe_export_carries_populated_paired_end_columns() -> None:
    """The eight columns are populated on every record when the
    DSL step is configured."""
    df = (
        _baseline(paired=True, r1=5, r2=5, insert=12)
        .run_records(n=3, seed=0)
        .to_dataframe()
    )
    for field in (
        "read_layout",
        "r1_sequence",
        "r2_sequence",
        "r1_start",
        "r1_end",
        "r2_start",
        "r2_end",
        "insert_size",
    ):
        assert field in df.columns
    # Per-row values reflect the fixed inputs.
    assert (df["read_layout"] == "paired_end").all()
    assert (df["insert_size"] == 12).all()
    assert (df["r1_sequence"].str.len() == 5).all()
    assert (df["r2_sequence"].str.len() == 5).all()
