"""Quality-score models for FASTQ output.

This module provides quality models that synthesise per-base Phred
scores from the assembled sequence. The models are deliberately
*position-driven*, not case-driven:

- GenAIRR uses a mixed-case convention internally (lowercase for
  the IMGT germline source plus `corrupt_quality` substitutes;
  uppercase for SHM mutations, PCR errors, and indel insertions).
  Those two pools overlap by case alone, so case is not a reliable
  "low quality" signal.
- For FASTQ output, the standard convention is uppercase bases +
  Phred quality strings. The Q-string carries confidence; we don't
  encode it in case.

Two signals the models *do* honour:

1. **Position within the read** — Illumina reads start at moderate
   quality, peak in the middle, and drop at the 3' end. Constant
   models ignore position; the Illumina model produces a smoothed
   trapezoid.
2. **Ambiguous base** — `N` / `n` → low Q (Phred convention is Q2,
   the lowest legal Illumina quality).

Phred+33 ASCII encoding (the Sanger / Illumina-1.8+ format) is
used downstream by the FASTQ writer:
    ASCII = Q + 33   →   Q=30 → ASCII 63 → "?"

Models in this file are **pure Python** and stateless: they take
a sequence string and return the same-length list of Q-values.
No simulator state, no RNG (tests can rely on deterministic
outputs).
"""

from __future__ import annotations

from typing import List


_AMBIGUOUS_Q = 2  # Phred convention for `N` bases (Illumina)


class ConstantQualityModel:
    """Every position gets the same `q`, except ``N`` bases which
    drop to `n_q`. Useful for sanity-test pipelines and when the
    user just wants Q30 across the board.

    Parameters:
        q: Phred quality assigned to non-N bases.
        n_q: Phred quality assigned to `N` / `n` bases. Defaults
            to Q2 per Illumina spec.
    """

    def __init__(self, q: int = 30, n_q: int = _AMBIGUOUS_Q) -> None:
        self._validate("q", q)
        self._validate("n_q", n_q)
        self.q = int(q)
        self.n_q = int(n_q)

    @staticmethod
    def _validate(name: str, value: int) -> None:
        if not isinstance(value, int) or isinstance(value, bool):
            raise TypeError(f"{name} must be int, got {type(value).__name__}")
        if not (0 <= value <= 93):
            raise ValueError(f"{name} must be in [0, 93] (Phred+33 range), got {value}")

    def quality_array(self, sequence: str) -> List[int]:
        out: List[int] = []
        for c in sequence:
            if c in ("N", "n"):
                out.append(self.n_q)
            else:
                out.append(self.q)
        return out


class IlluminaQualityModel:
    """Smoothed trapezoid that approximates an Illumina MiSeq /
    NextSeq quality profile: ramps up over the first `ramp_len`
    bases from `start_q` to `peak_q`, holds `peak_q` through the
    middle, then drops linearly to `end_q` over the last
    `tail_len` bases.

    `N` bases override the position-based score with `n_q`.

    Parameters:
        peak_q: peak quality in the middle of the read.
        start_q: quality at position 0 (read start, 5' end).
        end_q: quality at the final position (read end, 3' end).
        ramp_len: number of leading bases over which Q ramps from
            `start_q` to `peak_q`. Default 10.
        tail_len: number of trailing bases over which Q drops
            from `peak_q` to `end_q`. Default 30.
        n_q: Phred quality for `N` bases.
    """

    def __init__(
        self,
        peak_q: int = 35,
        start_q: int = 25,
        end_q: int = 18,
        ramp_len: int = 10,
        tail_len: int = 30,
        n_q: int = _AMBIGUOUS_Q,
    ) -> None:
        for name, value in [
            ("peak_q", peak_q),
            ("start_q", start_q),
            ("end_q", end_q),
            ("n_q", n_q),
        ]:
            ConstantQualityModel._validate(name, value)
        if ramp_len < 0:
            raise ValueError(f"ramp_len must be >= 0, got {ramp_len}")
        if tail_len < 0:
            raise ValueError(f"tail_len must be >= 0, got {tail_len}")
        self.peak_q = int(peak_q)
        self.start_q = int(start_q)
        self.end_q = int(end_q)
        self.ramp_len = int(ramp_len)
        self.tail_len = int(tail_len)
        self.n_q = int(n_q)

    def quality_array(self, sequence: str) -> List[int]:
        n = len(sequence)
        out: List[int] = []
        for i, c in enumerate(sequence):
            base = self._position_q(i, n)
            if c in ("N", "n"):
                out.append(self.n_q)
            else:
                out.append(base)
        return out

    def _position_q(self, i: int, n: int) -> int:
        if n == 0:
            return self.peak_q
        # Ramp: positions [0, ramp_len) interpolate start → peak.
        if i < self.ramp_len:
            if self.ramp_len == 0:
                return self.peak_q
            frac = (i + 1) / max(self.ramp_len, 1)
            return int(round(self.start_q + frac * (self.peak_q - self.start_q)))
        # Tail: last `tail_len` positions interpolate peak → end.
        tail_start = max(self.ramp_len, n - self.tail_len)
        if i >= tail_start and self.tail_len > 0:
            into_tail = i - tail_start
            tail_size = n - tail_start
            if tail_size <= 1:
                return self.end_q
            frac = into_tail / max(tail_size - 1, 1)
            return int(round(self.peak_q + frac * (self.end_q - self.peak_q)))
        # Middle: hold peak.
        return self.peak_q


def resolve_quality_model(name: str, **kwargs) -> object:
    """Build a model from a string name + kwargs. Used by
    `SimulationResult.to_fastq`."""
    if name == "constant":
        return ConstantQualityModel(**kwargs)
    if name == "illumina":
        return IlluminaQualityModel(**kwargs)
    raise ValueError(
        f"unknown quality model {name!r}; expected 'constant' or 'illumina'"
    )


def phred_to_ascii(q_array: List[int]) -> str:
    """Encode a list of Phred scores as a Phred+33 ASCII string.

    Out-of-range values are clamped to [0, 93] (the printable
    range of Phred+33). This is defensive — the model classes
    enforce the range at construction.
    """
    return "".join(chr(33 + min(max(q, 0), 93)) for q in q_array)
