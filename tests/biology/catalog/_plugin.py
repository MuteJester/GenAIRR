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
    from ._rust import rust_coverage, rust_covered_fields
    covered = {k: set(v) for k, v in rust_covered_fields().items()}
    for mechanism_id, fields in _COVERED.items():
        covered.setdefault(mechanism_id, set()).update(fields)
    rows = coverage_rows({k: tuple(v) for k, v in _REGISTRATIONS.items()},
                         rust_by_mechanism=rust_coverage(),
                         covered_by_mechanism={k: tuple(v) for k, v in covered.items()})
    colour = {"verified": {"green": True}, "partial": {"yellow": True},
              "unverified": {"red": True}}
    tw.section("biology coverage", "=")
    stage = None
    for row in rows:
        if row.stage != stage:
            stage = row.stage
            tw.line("")
            tw.line(f"{stage}", bold=True)
        sources = []
        if row.python_tests:
            sources.append(f"py:{len(row.python_tests)}")
        if row.rust_sources:
            sources.append(f"rs:{len(row.rust_sources)}")
        fields = f"{len(row.covered_fields)}/{len(row.truth_fields)} fields"
        tw.write(f"  {row.mechanism_id:34s} ")
        tw.write(f"{row.status:11s}", **colour[row.status])
        tw.line(f" {fields:16s} {' '.join(sources)}")
    counts = {s: sum(1 for r in rows if r.status == s)
              for s in ("verified", "partial", "unverified")}
    tw.line("")
    tw.line(
        f"  {counts['verified']} verified, {counts['partial']} partial, "
        f"{counts['unverified']} unverified, of {len(rows)} mechanisms")
