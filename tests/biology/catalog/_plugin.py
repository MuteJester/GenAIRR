"""The --biology-coverage report."""
from __future__ import annotations

from ._registry import coverage_rows

_REGISTRATIONS: dict[str, list[str]] = {}
_COVERED: dict[str, set[str]] = {}


def pytest_addoption(parser):
    parser.addoption("--biology-coverage", action="store_true",
                     help="report which biological mechanisms have verified truth")


def pytest_collection_modifyitems(items):
    _REGISTRATIONS.clear()
    _COVERED.clear()
    for item in items:
        registration = getattr(getattr(item, "function", None), "genairr_verifies", None)
        if not registration:
            continue
        mechanism_id, _claim, fields = registration
        _REGISTRATIONS.setdefault(mechanism_id, []).append(item.nodeid)
        _COVERED.setdefault(mechanism_id, set()).update(fields)


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    if not config.getoption("--biology-coverage", default=False):
        return
    tw = terminalreporter
    from ._rust import rust_coverage
    # Only registered Python tests count toward a field being covered. Rust
    # credit is shown in its own column and cannot move a status.
    rows = coverage_rows({k: tuple(v) for k, v in _REGISTRATIONS.items()},
                         rust_by_mechanism=rust_coverage(),
                         covered_by_mechanism={k: tuple(v) for k, v in _COVERED.items()})
    colour = {"verified": {"green": True}, "partial": {"yellow": True},
              "unverified": {"red": True}}
    tw.section("biology coverage", "=")
    stage = None
    for row in rows:
        if row.stage != stage:
            stage = row.stage
            tw.line("")
            tw.line(f"{stage}", bold=True)
        fields = f"{len(row.covered_fields)}/{len(row.truth_fields)} fields"
        asserted = f"rust-asserted:{len(row.rust_asserted)}" if row.rust_asserted else ""
        tw.write(f"  {row.mechanism_id:34s} ")
        tw.write(f"{row.status:11s}", **colour[row.status])
        tw.write(f" {fields:16s} ")
        tw.write(f"py:{len(row.python_tests)}".ljust(7))
        tw.line(f" {asserted}")
    counts = {s: sum(1 for r in rows if r.status == s)
              for s in ("verified", "partial", "unverified")}
    tw.line("")
    tw.line(
        f"  {counts['verified']} verified, {counts['partial']} partial, "
        f"{counts['unverified']} unverified, of {len(rows)} mechanisms")
    tw.line("")
    tw.line("  Status comes from registered Python tests only.", bold=True)
    tw.line("  rust-asserted counts engine_rs/tests/metadata_truth/ files recorded by")
    tw.line("  hand from their module docstrings. Those docstrings state intent, not")
    tw.line("  coverage, and spot-checks found credited fields the files never mention,")
    tw.line("  so they are shown but deliberately do not set a status.")
