"""T2-11: cross-platform reproducibility regression test.

The CI matrix runs the test suite on Linux, macOS, and Windows in
parallel — but pre-T2-11 the simulation outputs were never compared,
so a platform-specific divergence (different RNG bytes, different
float formatting in the libc, different allele sampling order) would
produce subtly-wrong sequences without any test catching it.

This test asserts that running ``Experiment.on(chain).run_to_file(
n=10, seed=42)`` produces output structurally identical to the
checked-in golden files at ``tests/golden/seed42_n10_*.tsv``. Any
field divergence — value, type, missing key — fails the test, and
the failure message points at the exact (row, column) for triage.

Three loci are covered to exercise both code paths:
  * ``human_igh``  — heavy BCR (V/D/J/C, full junction)
  * ``human_igk``  — light BCR (V/J only — no D segment)
  * ``human_tcrb`` — TCR (V/D/J — productive rate ~20%)

When the AIRR output format intentionally changes (e.g. a new column
or a deliberate algorithm change), regenerate the goldens with::

    python .private/scripts/regenerate_golden_files.py
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from GenAIRR import Experiment


GOLDEN_DIR = Path(__file__).resolve().parent / "golden"

# Fields the simulator emits as floats. Compared with `pytest.approx`
# (tight tolerance) instead of strict equality so the test tolerates
# benign libc-level float-repr differences while still catching real
# compute divergence.
_FLOAT_FIELDS = {"mutation_rate"}

# Tolerance for float comparisons: tight enough that any real numeric
# divergence (different RNG sample, different rate calc) trips, loose
# enough that "0.024999999" vs "0.025" (libc rounding) doesn't.
_FLOAT_REL = 1e-9
_FLOAT_ABS = 1e-12


_GOLDEN_CONFIGS = [
    ("human_igh",  "BCR-heavy with D"),
    ("human_igk",  "BCR-light without D"),
    ("human_tcrb", "TCR with D"),
]


def _read_tsv(path: Path) -> list[dict]:
    """Read a TSV file as a list of dicts (one per row)."""
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


@pytest.mark.parametrize("chain,label", _GOLDEN_CONFIGS)
def test_golden_reproducibility(chain: str, label: str, tmp_path):
    """Run the simulation against the same chain + seed + n that
    produced the checked-in golden, and assert structural identity.

    A failure here means this platform produced different sequences
    than the platform that generated the golden. Common causes:
      - C compiler used non-IEEE float math (extremely rare, would
        affect mutation_rate)
      - libc `printf`-formatting of doubles differs (handled by
        the float-tolerance branch below)
      - real RNG / sampling / allele-pool ordering bug
      - format change without regenerating goldens (run
        ``python .private/scripts/regenerate_golden_files.py``)
    """
    golden = GOLDEN_DIR / f"seed42_n10_{chain}.tsv"
    assert golden.exists(), (
        f"missing golden file {golden}; regenerate with "
        f"`python .private/scripts/regenerate_golden_files.py`")
    expected = _read_tsv(golden)

    # Reproduce the simulation under test.
    out = tmp_path / f"reproduced_{chain}.tsv"
    Experiment.on(chain).run_to_file(
        n=10, output_path=str(out), seed=42)
    actual = _read_tsv(out)

    # Sanity: same number of rows.
    assert len(actual) == len(expected), (
        f"[{label}] row count differs — golden has {len(expected)}, "
        f"this run produced {len(actual)}")
    # Sanity: same column set.
    assert actual[0].keys() == expected[0].keys(), (
        f"[{label}] column set differs:\n"
        f"  golden-only: {set(expected[0]) - set(actual[0])}\n"
        f"  actual-only: {set(actual[0]) - set(expected[0])}\n"
        f"  (regenerate goldens if a column was added intentionally)")

    # Field-by-field structural compare. First mismatch surfaces
    # exact (row, column) — much easier to triage than a diff dump.
    for row_idx, (exp, act) in enumerate(zip(expected, actual)):
        for col in exp.keys():
            exp_val = exp[col]
            act_val = act[col]
            if col in _FLOAT_FIELDS:
                _assert_float_equal(exp_val, act_val, row_idx, col, label)
            else:
                assert exp_val == act_val, (
                    f"[{label}] row {row_idx} col {col!r}: "
                    f"golden={exp_val!r}, actual={act_val!r}\n"
                    f"  (cross-platform divergence — see "
                    f"tests/test_reproducibility_golden.py docstring)")


def _assert_float_equal(expected_str: str, actual_str: str,
                        row_idx: int, col: str, label: str) -> None:
    """Compare two TSV-stringified floats with `_FLOAT_REL`/`_FLOAT_ABS`
    tolerance. Empty / missing values must match strictly (an empty
    field on one side and a populated float on the other is a real
    bug, not a rounding artifact)."""
    if expected_str == "" or actual_str == "":
        assert expected_str == actual_str, (
            f"[{label}] row {row_idx} col {col!r}: "
            f"golden={expected_str!r}, actual={actual_str!r} "
            f"(one side is empty)")
        return
    exp_f = float(expected_str)
    act_f = float(actual_str)
    # NaN handling — NaN != NaN by IEEE 754; treat NaN-NaN as match,
    # NaN-anything-else as mismatch.
    if math.isnan(exp_f) and math.isnan(act_f):
        return
    assert act_f == pytest.approx(exp_f, rel=_FLOAT_REL, abs=_FLOAT_ABS), (
        f"[{label}] row {row_idx} col {col!r} (float): "
        f"golden={exp_f}, actual={act_f}, "
        f"diff={act_f - exp_f:.3e} (tol rel={_FLOAT_REL}, abs={_FLOAT_ABS})\n"
        f"  (real numeric divergence — likely RNG or compute path bug)")
