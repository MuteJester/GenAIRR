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
