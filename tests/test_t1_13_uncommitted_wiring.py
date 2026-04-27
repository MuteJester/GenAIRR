"""T1-13: regression tests for AIRR fields and entry points that
shipped without Python-level coverage.

Specifically:
  - ``n_insertions`` / ``n_deletions`` AIRR fields under the indels
    feature flag.  The C side has ``test_indel_pool_exhaustion_count``
    (T0-6 regression) but no test exercised these from the user-facing
    Experiment.run path.
  - ``run_to_file`` row count and schema parity with ``run()``.
  - ``run_to_file`` reproducibility (byte-identical output for same seed).
"""
from __future__ import annotations

import csv
import filecmp
import os
import tempfile

from GenAIRR import Experiment
from GenAIRR.ops import with_indels


# ── n_insertions / n_deletions correctness ────────────────────────


class TestIndelCountFields:
    """Per-record `n_insertions` / `n_deletions` AIRR fields must be
    consistent with whether the indels feature is enabled."""

    def test_no_indels_when_disabled(self):
        """Default config does not enable indels — every record must
        report zero insertions and zero deletions."""
        result = Experiment.on("human_igh").run(n=20, seed=42)
        for i, r in enumerate(result):
            assert r["n_insertions"] == 0, \
                f"row {i}: n_insertions={r['n_insertions']} without indels enabled"
            assert r["n_deletions"] == 0, \
                f"row {i}: n_deletions={r['n_deletions']} without indels enabled"

    def test_indels_emitted_when_enabled(self):
        """With indels(prob=0.05), aggregate counts across many records
        must be > 0 for both insertions and deletions. Catches the
        T0-6 over-count regression and any regression that disables
        the indel step entirely."""
        result = (Experiment.on("human_igh")
                  .observe(with_indels(prob=0.05))
                  .run(n=200, seed=42))
        records = list(result)
        total_ins = sum(r["n_insertions"] for r in records)
        total_del = sum(r["n_deletions"] for r in records)
        assert total_ins > 0, \
            f"with_indels(0.05) produced 0 total insertions across 200 records"
        assert total_del > 0, \
            f"with_indels(0.05) produced 0 total deletions across 200 records"

    def test_indel_counts_non_negative_and_bounded(self):
        """Per-record sanity: counts must be ≥ 0 and not exceed
        ``sequence_length`` (a hard upper bound — every indel either
        inserts a base or removes one)."""
        result = (Experiment.on("human_igh")
                  .observe(with_indels(prob=0.10))
                  .run(n=100, seed=7))
        for i, r in enumerate(result):
            assert r["n_insertions"] >= 0, \
                f"row {i}: negative n_insertions={r['n_insertions']}"
            assert r["n_deletions"] >= 0, \
                f"row {i}: negative n_deletions={r['n_deletions']}"
            # n_insertions can't exceed sequence_length — every inserted
            # node lives in the sequence.
            assert r["n_insertions"] <= r["sequence_length"], \
                f"row {i}: n_insertions={r['n_insertions']} > seq_len={r['sequence_length']}"

    def test_indel_count_matches_germline_alignment_n_count(self):
        """An inserted indel node has segment=V/D/J/C but germline='\\0',
        so build_germline_alignment writes 'N' at that position. NP
        regions and adapter nodes also produce 'N'. With no corruption
        and known NP lengths, the invariant is::

            count('N' in germline_alignment)
                == np1_length + np2_length + n_insertions

        This catches any regression where the C engine miscounts
        inserted nodes (the T0-6 bug was exactly this — pool-exhausted
        inserts incremented n_insertions but didn't actually create
        the node)."""
        result = (Experiment.on("human_igh")
                  .observe(with_indels(prob=0.05))
                  .run(n=50, seed=42))
        for i, r in enumerate(result):
            germ = r["germline_alignment"]
            seq  = r["sequence"]
            assert len(seq) == len(germ), \
                f"row {i}: seq/germ length mismatch ({len(seq)} vs {len(germ)})"
            # Case-insensitive 'N' count.
            n_count = sum(1 for c in germ if c.upper() == "N")
            expected = (r["np1_length"] + r["np2_length"]
                        + r["n_insertions"])
            assert n_count == expected, (
                f"row {i}: germline-N count {n_count} != "
                f"np1({r['np1_length']}) + np2({r['np2_length']}) + "
                f"ins({r['n_insertions']}) = {expected}"
            )


# ── run_to_file: row count, schema parity, reproducibility ────────


class TestRunToFile:
    """`run_to_file` is the streaming entrypoint for n > 1M. It bypasses
    Python materialization and writes AIRR TSV directly. Must be in
    lockstep with `run()` on schema and reproducibility."""

    def test_row_count_matches_n(self):
        """run_to_file returns the number of rows written; the file
        must contain exactly that many data rows (plus header)."""
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as fp:
            path = fp.name
        try:
            for n in (1, 50, 250):
                written = (Experiment.on("human_igh")
                           .run_to_file(n=n, output_path=path, seed=42))
                assert written == n
                with open(path, encoding="utf-8") as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    rows = list(reader)
                assert len(rows) == n, \
                    f"run_to_file(n={n}): file has {len(rows)} rows"
        finally:
            os.unlink(path)

    def test_schema_parity_with_run(self):
        """run() returns dicts whose keys are the AIRR schema; run_to_file
        writes a TSV header. Both must enumerate the same column set
        in the same order."""
        # In-memory path
        recs = list(Experiment.on("human_igh").run(n=2, seed=42))
        run_keys = list(recs[0].keys())

        # Streaming path
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as fp:
            path = fp.name
        try:
            Experiment.on("human_igh").run_to_file(
                n=2, output_path=path, seed=42)
            with open(path, encoding="utf-8") as f:
                header = f.readline().rstrip("\n").split("\t")
            assert header == run_keys, (
                f"Schema drift between run() and run_to_file().\n"
                f"  run() keys:  {run_keys}\n"
                f"  TSV header:  {header}"
            )
        finally:
            os.unlink(path)

    def test_deterministic_with_same_seed(self):
        """Two run_to_file calls with the same seed and same Experiment
        must produce byte-identical TSV files."""
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as fp1, \
             tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as fp2:
            path1, path2 = fp1.name, fp2.name
        try:
            Experiment.on("human_igh").run_to_file(
                n=20, output_path=path1, seed=4242)
            Experiment.on("human_igh").run_to_file(
                n=20, output_path=path2, seed=4242)
            assert filecmp.cmp(path1, path2, shallow=False), \
                "run_to_file with same seed produced different bytes"
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_run_and_run_to_file_produce_equivalent_records(self):
        """Same Experiment + same seed via run() and run_to_file()
        must produce semantically equal records on shared columns."""
        n = 5
        recs = list(Experiment.on("human_igh").run(n=n, seed=7))

        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as fp:
            path = fp.name
        try:
            Experiment.on("human_igh").run_to_file(
                n=n, output_path=path, seed=7)
            with open(path, encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                tsv_rows = list(reader)
        finally:
            os.unlink(path)

        assert len(recs) == len(tsv_rows) == n
        for i, (rec, tsv_row) in enumerate(zip(recs, tsv_rows)):
            # Sequence is the most fundamental cross-check.
            assert rec["sequence"] == tsv_row["sequence"], \
                f"row {i}: run() sequence != run_to_file sequence"
            # productive flag — TSV stores "T"/"F", run() returns bool.
            tsv_prod = (tsv_row["productive"] == "T")
            assert rec["productive"] == tsv_prod, \
                f"row {i}: productive flag drift"
