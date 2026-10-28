"""End-to-end implementation tests for **paired-end FASTQ
export** — `SimulationResult.to_paired_fastq(r1_path, r2_path,
…)`.

The slice's headline contract:

- Two files, one R1 record + one R2 record per AIRR record.
- Headers ``@{sequence_id}/1`` and ``@{sequence_id}/2`` —
  Illumina-portable suffix convention, no metadata pipe.
- R1 / R2 bodies are the AIRR ``r1_sequence`` / ``r2_sequence``
  fields uppercased; R2 is **already** reverse-complemented at
  projection time and is written verbatim.
- Quality strings are produced by the existing pluggable models
  (``constant`` / ``illumina``) consumed by the single-end
  exporter; each read is scored independently.
- Refuses to clobber existing files unless ``overwrite=True``.
- Refuses to write non-paired records, empty windows, or
  length-mismatched quality arrays.

See [`docs/fastq_export_design.md`](../docs/fastq_export_design.md)
for the audit and per-spec rationale.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

import GenAIRR as ga


# ──────────────────────────────────────────────────────────────────
# Shared fixtures
# ──────────────────────────────────────────────────────────────────


def _paired_records(n: int = 3, seed: int = 4242, *, r1_length: int = 140,
                    insert_size: int = 280, **mutate_kwargs):
    """Run a small human-IGH paired-end batch and return the
    :class:`SimulationResult`. Optional kwargs are forwarded to
    ``.mutate(...)`` so individual tests can opt into corruption
    / SHM."""
    exp = ga.Experiment.on("human_igh").recombine()
    if mutate_kwargs:
        exp = exp.mutate(**mutate_kwargs)
    exp = exp.paired_end(r1_length=r1_length, insert_size=insert_size)
    return exp.run_records(n=n, seed=seed)


# ──────────────────────────────────────────────────────────────────
# Spec 1 — fixed paired-end result writes two FASTQ files
# ──────────────────────────────────────────────────────────────────


def test_writes_two_fastq_files_with_matching_record_counts(tmp_path):
    """Each AIRR record produces exactly four lines per file
    (`@header` / `sequence` / `+` / `quality`). Two records ⇒ 8
    lines per file across both."""
    n = 3
    result = _paired_records(n=n)
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    result.to_paired_fastq(str(r1_path), str(r2_path), quality="constant", q=30)

    r1_lines = r1_path.read_text(encoding="utf-8").splitlines()
    r2_lines = r2_path.read_text(encoding="utf-8").splitlines()
    assert len(r1_lines) == n * 4
    assert len(r2_lines) == n * 4
    # Every fourth line starts with `@` (the FASTQ record header).
    for i in range(n):
        assert r1_lines[i * 4].startswith("@")
        assert r2_lines[i * 4].startswith("@")
        assert r1_lines[i * 4 + 2] == "+"
        assert r2_lines[i * 4 + 2] == "+"


# ──────────────────────────────────────────────────────────────────
# Spec 2 — headers use sequence_id/1 and sequence_id/2
# ──────────────────────────────────────────────────────────────────


def test_headers_use_sequence_id_with_illumina_suffix(tmp_path):
    """Header format: ``@{sequence_id}/1`` for R1 and
    ``@{sequence_id}/2`` for R2 — the audit's §5 canonical
    shape. No `|v_call=…|j_call=…` metadata pipe (tool
    portability — STAR < 2.7 splits on `|`)."""
    result = _paired_records(n=2)
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    result.to_paired_fastq(str(r1_path), str(r2_path), quality="constant", q=30)
    r1_lines = r1_path.read_text(encoding="utf-8").splitlines()
    r2_lines = r2_path.read_text(encoding="utf-8").splitlines()
    for i, rec in enumerate(result.records):
        sequence_id = rec["sequence_id"]
        assert r1_lines[i * 4] == f"@{sequence_id}/1", (
            f"R1 header mismatch at record {i}: got {r1_lines[i * 4]!r}, "
            f"expected @{sequence_id}/1"
        )
        assert r2_lines[i * 4] == f"@{sequence_id}/2", (
            f"R2 header mismatch at record {i}: got {r2_lines[i * 4]!r}, "
            f"expected @{sequence_id}/2"
        )
        # The header must NOT carry a `|v_call=...` metadata pipe.
        assert "|" not in r1_lines[i * 4]
        assert "|" not in r2_lines[i * 4]


# ──────────────────────────────────────────────────────────────────
# Spec 3 — R1/R2 sequence lines match AIRR fields exactly
# ──────────────────────────────────────────────────────────────────


def test_r1_and_r2_sequences_match_airr_fields_uppercased(tmp_path):
    """The body line at position ``4i + 1`` in each file equals
    the corresponding AIRR field uppercased. R2 is already
    reverse-complemented at projection time — the writer outputs
    it verbatim (no second RC)."""
    result = _paired_records(n=4)
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    result.to_paired_fastq(str(r1_path), str(r2_path), quality="constant", q=30)
    r1_lines = r1_path.read_text(encoding="utf-8").splitlines()
    r2_lines = r2_path.read_text(encoding="utf-8").splitlines()
    for i, rec in enumerate(result.records):
        assert r1_lines[i * 4 + 1] == rec["r1_sequence"].upper()
        assert r2_lines[i * 4 + 1] == rec["r2_sequence"].upper()


# ──────────────────────────────────────────────────────────────────
# Spec 4 — quality string length matches read length
# ──────────────────────────────────────────────────────────────────


def test_quality_string_length_matches_read_length(tmp_path):
    """Every quality string (`4i + 3`) is exactly the same length
    as its corresponding sequence body (`4i + 1`). Pinned per
    read, per file."""
    result = _paired_records(n=5)
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    result.to_paired_fastq(str(r1_path), str(r2_path), quality="illumina")
    r1_lines = r1_path.read_text(encoding="utf-8").splitlines()
    r2_lines = r2_path.read_text(encoding="utf-8").splitlines()
    for i in range(5):
        assert len(r1_lines[i * 4 + 1]) == len(r1_lines[i * 4 + 3]), (
            f"R1 record {i}: seq len {len(r1_lines[i * 4 + 1])} != "
            f"quality len {len(r1_lines[i * 4 + 3])}"
        )
        assert len(r2_lines[i * 4 + 1]) == len(r2_lines[i * 4 + 3]), (
            f"R2 record {i}: seq len {len(r2_lines[i * 4 + 1])} != "
            f"quality len {len(r2_lines[i * 4 + 3])}"
        )


# ──────────────────────────────────────────────────────────────────
# Spec 5 — Constant-quality model works
# ──────────────────────────────────────────────────────────────────


def test_constant_quality_writes_q30_phred(tmp_path):
    """``quality="constant", q=30`` produces a quality string
    whose every non-N base is Q30 (Phred+33 ASCII ``'?'`` →
    chr(33 + 30) == '?')."""
    result = _paired_records(n=1)
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    result.to_paired_fastq(str(r1_path), str(r2_path), quality="constant", q=30)
    r1_lines = r1_path.read_text(encoding="utf-8").splitlines()
    seq = r1_lines[1]
    q = r1_lines[3]
    expected = "".join("#" if c in ("N", "n") else "?" for c in seq)
    assert q == expected, (
        f"R1 quality string didn't match expected constant-Q30 "
        f"pattern: q={q[:40]!r} expected={expected[:40]!r}"
    )


# ──────────────────────────────────────────────────────────────────
# Spec 6 — Illumina model smoke (independent shape per read)
# ──────────────────────────────────────────────────────────────────


def test_illumina_quality_runs_and_quality_chars_are_printable(tmp_path):
    """The Illumina model produces a printable Phred+33 quality
    string (ASCII 33–126) of the right length. The model is
    consulted **independently** for R1 and R2 — both reads get
    their own ramp-up + tail-down shape starting from position 0."""
    result = _paired_records(n=1)
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    result.to_paired_fastq(str(r1_path), str(r2_path), quality="illumina")
    for path in (r1_path, r2_path):
        lines = path.read_text(encoding="utf-8").splitlines()
        seq = lines[1]
        q = lines[3]
        # Every Phred+33 char is in the printable range.
        for ch in q:
            assert 33 <= ord(ch) <= 126, (
                f"non-printable Phred+33 byte {ord(ch)} in {path.name}"
            )
        assert len(q) == len(seq)


# ──────────────────────────────────────────────────────────────────
# Spec 7 — Non-paired result raises
# ──────────────────────────────────────────────────────────────────


def test_non_paired_result_raises_value_error(tmp_path):
    """A `SimulationResult` produced without ``.paired_end(...)``
    has ``read_layout=""`` on every record. The writer raises
    ``ValueError`` on the first record with a hint pointing at
    ``Experiment.paired_end(...)``."""
    result = (
        ga.Experiment.on("human_igh")
        .recombine()
        .run_records(n=2, seed=4242)
    )
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    with pytest.raises(ValueError, match="paired-end FASTQ export requires"):
        result.to_paired_fastq(
            str(r1_path), str(r2_path), quality="constant", q=30
        )


# ──────────────────────────────────────────────────────────────────
# Spec 8 — Empty R1/R2 raises
# ──────────────────────────────────────────────────────────────────


def test_empty_r1_sequence_raises(tmp_path):
    """A consumer who hand-edits the record dict to blank out
    ``r1_sequence`` triggers a ``ValueError`` from the writer —
    we never silently emit a zero-length FASTQ read."""
    result = _paired_records(n=1)
    # Tamper the record dict's r1_sequence directly.
    result.records[0]["r1_sequence"] = ""
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    with pytest.raises(ValueError, match="empty r1_sequence"):
        result.to_paired_fastq(
            str(r1_path), str(r2_path), quality="constant", q=30
        )


def test_empty_r2_sequence_raises(tmp_path):
    result = _paired_records(n=1)
    result.records[0]["r2_sequence"] = ""
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    with pytest.raises(ValueError, match="empty r2_sequence"):
        result.to_paired_fastq(
            str(r1_path), str(r2_path), quality="constant", q=30
        )


# ──────────────────────────────────────────────────────────────────
# Spec 9 — Existing path overwrite discipline
# ──────────────────────────────────────────────────────────────────


def test_existing_path_without_overwrite_raises(tmp_path):
    """The default ``overwrite=False`` refuses to clobber existing
    files. Both R1 and R2 paths get the same guard; the writer
    raises before any quality model is constructed."""
    result = _paired_records(n=1)
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    # First write succeeds.
    result.to_paired_fastq(str(r1_path), str(r2_path), quality="constant", q=30)
    # Second write with default overwrite=False refuses.
    with pytest.raises(FileExistsError, match="r1_path"):
        result.to_paired_fastq(
            str(r1_path), str(r2_path), quality="constant", q=30
        )


def test_overwrite_true_allows_replacement(tmp_path):
    """``overwrite=True`` replaces existing files. The new
    content reflects the second write."""
    result_a = _paired_records(n=1, seed=1111)
    result_b = _paired_records(n=2, seed=2222)
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    # First write — 1 record.
    result_a.to_paired_fastq(str(r1_path), str(r2_path), quality="constant", q=30)
    assert len(r1_path.read_text().splitlines()) == 4
    # Second write — 2 records — replaces.
    result_b.to_paired_fastq(
        str(r1_path), str(r2_path), quality="constant", q=30, overwrite=True
    )
    assert len(r1_path.read_text().splitlines()) == 8
    assert len(r2_path.read_text().splitlines()) == 8


def test_overwrite_guard_fires_when_only_r2_path_exists(tmp_path):
    """If only one of the two output paths already exists, the
    writer still refuses — the guard runs over both paths
    independently."""
    result = _paired_records(n=1)
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    # Pre-create only r2_path.
    r2_path.write_text("preexisting r2 data\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="r2_path"):
        result.to_paired_fastq(
            str(r1_path), str(r2_path), quality="constant", q=30
        )
    # And r1_path was not partially written — the writer raised
    # before opening either file handle.
    assert not r1_path.exists()


# ──────────────────────────────────────────────────────────────────
# Replay round-trip (parallels the single-end test discipline)
# ──────────────────────────────────────────────────────────────────


def test_paired_fastq_decodes_as_two_synchronized_files(tmp_path):
    """For each AIRR record at index ``i``, the R1 file's
    `(4i, 4i+1, 4i+2, 4i+3)` block and the R2 file's same block
    refer to the SAME ``sequence_id`` — `/1` for R1 and `/2` for
    R2. Parallels the single-end
    ``test_to_fastq_pairs_with_constant_q_is_decodable``
    discipline."""
    result = _paired_records(n=5)
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    result.to_paired_fastq(str(r1_path), str(r2_path), quality="constant", q=30)
    r1_lines = r1_path.read_text(encoding="utf-8").splitlines()
    r2_lines = r2_path.read_text(encoding="utf-8").splitlines()
    assert len(r1_lines) == len(r2_lines) == 5 * 4
    for i in range(5):
        r1_id = r1_lines[i * 4][1:].rsplit("/", 1)[0]
        r2_id = r2_lines[i * 4][1:].rsplit("/", 1)[0]
        assert r1_id == r2_id, (
            f"R1/R2 desynchronised at index {i}: "
            f"R1 id={r1_id!r} R2 id={r2_id!r}"
        )
        # And the suffix is correct on each side.
        assert r1_lines[i * 4].endswith("/1")
        assert r2_lines[i * 4].endswith("/2")
