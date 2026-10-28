"""What a biological mechanism is, for the purposes of verification."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Mechanism:
    """One step GenAIRR simulates, and the truth it is responsible for recording.

    ``truth_fields`` is the point of the whole structure: it names the record
    fields this mechanism is answerable for. A test that verifies the mechanism
    is making a claim about those fields, and a mechanism with no verified
    claim is a field nobody has checked.
    """

    id: str
    stage: str
    summary: str
    truth_fields: tuple[str, ...]
    engine_passes: tuple[str, ...] = ()
    references: tuple[str, ...] = ()
    notes: str = ""
    #: Set when the mechanism's real truth is not exported at all. Without it,
    #: such a mechanism borrows whatever field happens to exist and reports as
    #: merely untested, which confuses "no one has written the test" with
    #: "there is nothing to test against".
    truth_gap: str = ""

    def __post_init__(self) -> None:
        if not self.id.startswith(f"{self.stage}."):
            raise ValueError(f"mechanism id {self.id!r} is not prefixed with stage {self.stage!r}")
