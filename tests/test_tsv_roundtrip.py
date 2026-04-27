"""T1-11 Python-side round-trip: every TSV row that the C engine
emits must parse back through `csv.DictReader` with exactly
N_COLUMNS keys, no leakage of stray fields.

The C engine sanitizes embedded tabs/newlines in string fields
inline (see `TsvWriter` in airr.c). This test exercises the full
write_to_file → read_with_csv path on real simulator output, plus
a regression for the column-count invariant on a heavy-feature
configuration (where `note` would most likely accumulate text).
"""
from __future__ import annotations

import csv
import os
import tempfile

from GenAIRR import Experiment, Productivity
from GenAIRR.ops import (
    rate, model, with_indels, with_5prime_loss, with_3prime_loss,
)


def _read_tsv_rows(path: str) -> tuple[list[str], list[dict]]:
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = list(reader.fieldnames or [])
        rows = list(reader)
    return header, rows


def test_normal_run_produces_well_formed_tsv():
    """Baseline: a vanilla simulation produces a TSV every row of which
    has exactly the same number of fields as the header."""
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as fp:
        path = fp.name
    try:
        n_written = (Experiment.on("human_igh")
                     .run_to_file(n=200, output_path=path, seed=42))
        assert n_written == 200

        header, rows = _read_tsv_rows(path)
        n_cols = len(header)
        assert n_cols >= 50, f"Header has only {n_cols} columns?"
        assert len(rows) == 200
        for i, row in enumerate(rows):
            # csv.DictReader puts excess fields under None when a row
            # has too many columns; missing columns leave keys at None.
            assert None not in row, \
                f"row {i}: extra columns leaked → {row.get(None)!r}"
            assert all(v is not None for v in row.values()), \
                f"row {i}: missing fields"
            assert len(row) == n_cols, \
                f"row {i}: {len(row)} fields vs {n_cols} header columns"
    finally:
        os.unlink(path)


def test_heavy_feature_run_column_count_invariant():
    """Heavy config: SHM + indels + 5'/3' loss all emit notes and
    annotation strings. Verify the column-count invariant still
    holds and the truncation note from T1-2 doesn't leak tabs."""
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as fp:
        path = fp.name
    try:
        n = 300
        (Experiment.on("human_igh")
         .mutate(rate(0.03, 0.08), model("s5f"))
         .sequence(with_5prime_loss(min_remove=10, max_remove=50),
                    with_3prime_loss(min_remove=10, max_remove=50))
         .observe(with_indels(prob=0.05))
         .run_to_file(n=n, output_path=path, seed=123,
                      productivity=Productivity.PRODUCTIVE_MIXED))

        header, rows = _read_tsv_rows(path)
        n_cols = len(header)
        assert len(rows) == n
        seen_partial_note = False
        for i, row in enumerate(rows):
            assert None not in row, \
                f"row {i}: extra columns leaked → {row.get(None)!r}"
            assert len(row) == n_cols
            # Sanity: note column should never contain a literal tab
            # or newline (sanitized to space).
            note = row.get("note", "")
            assert "\t" not in note, f"row {i} note has tab: {note!r}"
            assert "\n" not in note, f"row {i} note has newline: {note!r}"
            if "junction partial" in note or "junction fully removed" in note:
                seen_partial_note = True
        # Heavy 3' trim → at least some records hit the T1-2 truncation
        # path. This is the high-risk path for tab/newline leakage,
        # so make sure it actually fired.
        assert seen_partial_note, \
            "expected at least one junction-truncation note under heavy 3' loss"
    finally:
        os.unlink(path)


def test_truncated_run_produces_well_formed_tsv():
    """Productive-only mode emits records with fewer notes but the
    column invariant must still hold."""
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as fp:
        path = fp.name
    try:
        (Experiment.on("human_igh")
         .mutate(rate(0.02, 0.05), model("s5f"))
         .run_to_file(n=100, output_path=path, seed=7,
                      productivity=Productivity.PRODUCTIVE_ONLY))
        header, rows = _read_tsv_rows(path)
        n_cols = len(header)
        for i, row in enumerate(rows):
            assert None not in row
            assert len(row) == n_cols
    finally:
        os.unlink(path)
