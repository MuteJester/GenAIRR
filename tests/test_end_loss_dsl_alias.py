"""DSL naming: `end_loss_*prime` is canonical, `primer_trim_*prime` is alias.

Pins the *naming-cleanup* slice that promotes
:meth:`Experiment.end_loss_5prime` / :meth:`end_loss_3prime` as the
biologically-named public API for the engine's observation-stage
end-loss pass. The legacy :meth:`primer_trim_5prime` /
:meth:`primer_trim_3prime` survive as backwards-compatible aliases.

Both spellings must:

- route to the same trace address (``corrupt.end_loss.{5,3}``);
- populate the same AIRR field (``end_loss_{5,3}_length``);
- produce byte-identical AIRR records for the same seed;
- not break when chained or interleaved.

Plus a small ``describe()`` check that the new wording surfaces in
the narrative for users browsing experiments interactively.
"""
from __future__ import annotations

import pytest

import GenAIRR as ga
from GenAIRR import _engine as ge


# ──────────────────────────────────────────────────────────────────
# Shared deterministic VJ fixture (matches the provenance suite).
# ──────────────────────────────────────────────────────────────────


def _vj_refdata() -> "ge.RefDataConfig":
    cfg = ge.RefDataConfig.vj()
    cfg.add_v_allele("v1*01", "v1", b"AAACCCGGG", anchor=6)
    cfg.add_j_allele("j1*01", "j1", b"TTTAAA", anchor=0)
    return cfg


def _baseline() -> "ga.Experiment":
    return (
        ga.Experiment.on(_vj_refdata())
        .allow_curatable_refdata()
        .recombine(np1_lengths=[(3, 1.0)])
        .trim(enabled=False)
    )


def _record(exp: "ga.Experiment", seed: int = 0) -> dict:
    result = exp.run_records(n=1, seed=seed)
    assert len(result.records) == 1
    return result.records[0]


# ──────────────────────────────────────────────────────────────────
# 1. Trace address — same regardless of which DSL method was used
# ──────────────────────────────────────────────────────────────────


def test_end_loss_5prime_emits_corrupt_end_loss_5_trace_record() -> None:
    """The canonical DSL method routes to the same trace address the
    engine uses (``corrupt.end_loss.5``), so trace consumers don't
    need to know about the DSL surface naming at all."""
    outcomes = _baseline().end_loss_5prime(length=2).run(n=1, seed=0)
    records = outcomes[0].trace().choices()
    addrs = {r.address for r in records}
    assert "corrupt.end_loss.5" in addrs
    # And the canonical method does NOT emit anything under a
    # legacy ``corrupt.primer_trim.*`` address — there is no such
    # address in the engine.
    for r in records:
        assert not r.address.startswith("corrupt.primer_trim"), (
            f"unexpected legacy trace address: {r.address}"
        )


def test_end_loss_3prime_emits_corrupt_end_loss_3_trace_record() -> None:
    outcomes = _baseline().end_loss_3prime(length=2).run(n=1, seed=0)
    addrs = {r.address for r in outcomes[0].trace().choices()}
    assert "corrupt.end_loss.3" in addrs


# ──────────────────────────────────────────────────────────────────
# 2. AIRR field plumbing — `end_loss_{5,3}_length` is populated
# ──────────────────────────────────────────────────────────────────


def test_end_loss_5prime_populates_airr_end_loss_5_length() -> None:
    rec = _record(_baseline().end_loss_5prime(length=2))
    assert rec["end_loss_5_length"] == 2
    assert rec["end_loss_3_length"] == 0


def test_end_loss_3prime_populates_airr_end_loss_3_length() -> None:
    rec = _record(_baseline().end_loss_3prime(length=2))
    assert rec["end_loss_5_length"] == 0
    assert rec["end_loss_3_length"] == 2


def test_end_loss_fields_default_to_zero_with_no_end_loss_step() -> None:
    rec = _record(_baseline())
    assert rec["end_loss_5_length"] == 0
    assert rec["end_loss_3_length"] == 0


# ──────────────────────────────────────────────────────────────────
# 3. Alias contract — new name ≡ old name under the same seed
# ──────────────────────────────────────────────────────────────────


_COMPARABLE_FIELDS = (
    "sequence",
    "v_call",
    "j_call",
    "v_trim_5",
    "v_trim_3",
    "j_trim_5",
    "j_trim_3",
    "junction",
    "junction_length",
    "vj_in_frame",
    "stop_codon",
    "productive",
    "end_loss_5_length",
    "end_loss_3_length",
    "n_indels",
)


def _assert_records_equal(left: dict, right: dict) -> None:
    """Two records must agree on every load-bearing AIRR column."""
    for field in _COMPARABLE_FIELDS:
        assert left[field] == right[field], (
            f"records diverge on {field!r}: "
            f"end_loss spelling={left[field]!r} vs primer_trim spelling={right[field]!r}"
        )


def test_end_loss_5prime_byte_identical_to_primer_trim_5prime_alias() -> None:
    """The two spellings must produce the same AIRR record for the
    same seed. Anything else means the alias is leaking parameters
    or routing differently."""
    new = _record(_baseline().end_loss_5prime(length=2))
    old = _record(_baseline().primer_trim_5prime(length=2))
    _assert_records_equal(new, old)


def test_end_loss_3prime_byte_identical_to_primer_trim_3prime_alias() -> None:
    new = _record(_baseline().end_loss_3prime(length=2))
    old = _record(_baseline().primer_trim_3prime(length=2))
    _assert_records_equal(new, old)


def test_chained_end_loss_byte_identical_to_chained_primer_trim() -> None:
    """Chaining both ends must produce the same record under either
    spelling — including the case where the user mixes the two
    spellings in a single chain (legacy migration scenario)."""
    new = _record(
        _baseline().end_loss_5prime(length=2).end_loss_3prime(length=1)
    )
    old = _record(
        _baseline().primer_trim_5prime(length=2).primer_trim_3prime(length=1)
    )
    mixed = _record(
        _baseline().end_loss_5prime(length=2).primer_trim_3prime(length=1)
    )
    _assert_records_equal(new, old)
    _assert_records_equal(new, mixed)


@pytest.mark.parametrize(
    "length",
    [0, 1, 3, (0, 4), [(1, 1.0), (3, 1.0)]],
    ids=["zero", "one", "three", "range", "empirical"],
)
def test_alias_accepts_every_length_shape_the_canonical_accepts(length) -> None:
    """The alias is a thin pass-through, so every ``length=`` shape
    accepted by the canonical method must work on the alias too."""
    # Just exercising the construction path is enough — if the alias
    # mis-handled any shape, this would raise at recording time.
    new = _baseline().end_loss_5prime(length=length)
    old = _baseline().primer_trim_5prime(length=length)
    new.run(n=1, seed=0)
    old.run(n=1, seed=0)


# ──────────────────────────────────────────────────────────────────
# 4. describe() narrative — surfaces "end-loss" wording
# ──────────────────────────────────────────────────────────────────


def test_describe_uses_end_loss_wording_for_canonical_method() -> None:
    text = _baseline().end_loss_5prime(length=2).describe()
    assert "5'-end loss" in text
    # And — load-bearing for users browsing the narrative — the
    # describe text should NOT carry the legacy "(primer/adapter
    # trim)" parenthetical that used to confuse the surface
    # vocabulary. The engine pass is end-loss; the narrative says
    # end-loss.
    assert "primer/adapter trim" not in text


def test_describe_alias_renders_same_narrative_as_canonical() -> None:
    """The alias compiles to the same pipeline step, so its
    describe() output must match the canonical method's output."""
    new = _baseline().end_loss_5prime(length=2).describe()
    old = _baseline().primer_trim_5prime(length=2).describe()
    assert new == old


# ──────────────────────────────────────────────────────────────────
# 5. Lockstep with the audit doc
# ──────────────────────────────────────────────────────────────────


def test_audit_doc_lists_end_loss_as_canonical_dsl_surface() -> None:
    """Section 6.2 of ``docs/primer_trim_end_loss_audit.md`` must
    describe the rename as **resolved**, with ``end_loss_*prime``
    listed as canonical. A future slice that reverts the rename
    would also flip this test, surfacing the doc/test drift."""
    from pathlib import Path

    doc = (
        Path(__file__).resolve().parent.parent
        / "docs"
        / "primer_trim_end_loss_audit.md"
    )
    assert doc.is_file()
    text = doc.read_text(encoding="utf-8")
    assert "end_loss_5prime" in text
    assert "end_loss_3prime" in text
    assert "Resolved" in text or "resolved" in text
