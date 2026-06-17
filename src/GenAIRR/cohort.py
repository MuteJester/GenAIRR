"""Cohort orchestration results — N subjects each with their own diploid
genotype, produced by :meth:`GenAIRR.Experiment.run_cohort`.

``CohortResult`` is built from per-subject ``CohortSubjectResult`` entries and
exposes per-subject access (``result_for`` / ``refdata_for``) plus combined
export. It explicitly stores each subject's refdata because ``SimulationResult``
does not preserve it (needed for ``validate_records`` and novel-allele subjects).
See ``.private/specs/2026-06-17-genotype-cohorts-design.md``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List


def _resolve_subject_ids(raw_ids: List[Any]) -> List[str]:
    """Resolve per-subject IDs: all-None -> subject_0..N-1; all-present ->
    require unique after str() normalization; mixed -> raise."""
    none_count = sum(1 for i in raw_ids if i is None)
    n = len(raw_ids)
    if none_count == n:
        return [f"subject_{i}" for i in range(n)]
    if none_count != 0:
        raise ValueError(
            "run_cohort: some genotypes have a subject_id and others don't; set "
            "subject_id on all genotypes or none")
    ids = [str(i) for i in raw_ids]
    if len(ids) != len(set(ids)):
        raise ValueError("run_cohort: duplicate subject_id after normalization")
    return ids


def _check_count(c, where: str) -> int:
    if isinstance(c, bool) or not isinstance(c, int) or c < 0:
        raise ValueError(f"{where} must be an int >= 0, got {c!r}")
    return c


def _resolve_counts(n_genotypes: int, n_per_subject, counts) -> List[int]:
    """Resolve per-subject record counts. ``counts`` (a parallel sequence)
    overrides ``n_per_subject`` when supplied; entries are validated as
    non-bool ints >= 0."""
    _check_count(n_per_subject, "n_per_subject")
    if counts is None:
        return [int(n_per_subject)] * n_genotypes
    counts = list(counts)
    if len(counts) != n_genotypes:
        raise ValueError(
            f"run_cohort: counts length {len(counts)} != number of genotypes "
            f"{n_genotypes}")
    return [_check_count(c, f"counts[{i}]") for i, c in enumerate(counts)]


@dataclass(frozen=True)
class CohortSubjectResult:
    """One subject's slice of a cohort run."""

    subject_id: str
    genotype: Any            # GenAIRR.genotype.Genotype (snapshot)
    result: Any              # GenAIRR.result.SimulationResult
    refdata: Any             # the refdata this subject ran against (base or effective)
    seed: int                # derived per-subject sub-seed
    count: int               # records requested for this subject


class CohortResult:
    """Combined result of a :meth:`Experiment.run_cohort` run.

    ``subjects`` is the single source of truth; ``subject_ids`` / ``genotypes`` /
    ``results`` are derived so parallel lists can't drift."""

    def __init__(self, subjects: List[CohortSubjectResult]):
        self._subjects = list(subjects)
        self._by_id = {s.subject_id: s for s in self._subjects}

    @property
    def subjects(self) -> List[CohortSubjectResult]:
        return list(self._subjects)

    @property
    def subject_ids(self) -> List[str]:
        return [s.subject_id for s in self._subjects]

    @property
    def genotypes(self) -> List[Any]:
        return [s.genotype for s in self._subjects]

    @property
    def results(self) -> List[Any]:
        return [s.result for s in self._subjects]

    def result_for(self, subject_id: str):
        return self._by_id[subject_id].result

    def refdata_for(self, subject_id: str):
        return self._by_id[subject_id].refdata

    @property
    def records(self) -> List[Dict[str, Any]]:
        """A fresh concatenated list of every subject's record dicts (each
        already ``subject_id``-tagged and ``sequence_id``-namespaced). Mutating
        the returned list does not affect the cohort."""
        out: List[Dict[str, Any]] = []
        for s in self._subjects:
            out.extend(s.result.records)
        return out

    def __len__(self) -> int:
        return sum(len(s.result) for s in self._subjects)

    def __repr__(self) -> str:
        return f"<CohortResult subjects={len(self._subjects)} records={len(self)}>"

    # ── combined export ─────────────────────────────────────────────
    def to_dataframe(self, *, airr_strict: bool = False):
        """Combined DataFrame over every subject's records. Column set is the
        stable union across subjects (pandas fills missing keys with NaN)."""
        import pandas as pd

        records = self.records
        if airr_strict:
            from .result import _to_airr_strict
            records = [_to_airr_strict(r) for r in records]
        return pd.DataFrame(records)

    def to_csv(self, path: str, *, airr_strict: bool = False) -> None:
        """Write the combined records as CSV (union columns guaranteed by the
        DataFrame)."""
        self.to_dataframe(airr_strict=airr_strict).to_csv(path, index=False)

    def to_fasta(self, path: str) -> None:
        """Write one FASTA record per concatenated record. Header is exactly
        ``>{sequence_id}`` (the namespaced id); body is the record's
        ``sequence``. Unlike ``SimulationResult.to_fasta`` (whose headers come
        from the enumerate index), this keeps cohort headers globally unique."""
        with open(path, "w", encoding="utf-8") as fh:
            for rec in self.records:
                sid = rec.get("sequence_id", "")
                seq = rec.get("sequence", "")
                fh.write(f">{sid}\n{seq}\n")
