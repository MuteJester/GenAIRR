"""The --biology-coverage report."""
from __future__ import annotations

from ._registry import coverage_rows

#: Registrations seen at collection, before any test has run.
_COLLECTED: dict[str, list[tuple[str, tuple[str, ...]]]] = {}
#: Node ids the run observed PASS. A mechanism cannot be verified by a test
#: the same run watched fail, so credit is withheld until the outcome is known.
_PASSED: set[str] = set()
_REGISTRATIONS: dict[str, list[str]] = {}
_COVERED: dict[str, set[str]] = {}


def pytest_addoption(parser):
    parser.addoption("--biology-coverage", action="store_true",
                     help="report which biological mechanisms have verified truth")


#: Markers that mean a test does not establish its claim. A registration on
#: such a test must not credit the mechanism: known defects are pinned as
#: strict xfails, and counting those would be the same lie as a prose-sourced
#: green row, arriving by a different route.
NON_ESTABLISHING = ("xfail", "skip", "skipif")


def _establishes_its_claim(item) -> bool:
    return not any(mark.name in NON_ESTABLISHING for mark in item.iter_markers())


def pytest_collection_modifyitems(items):
    _COLLECTED.clear()
    _PASSED.clear()
    _REGISTRATIONS.clear()
    _COVERED.clear()
    for item in items:
        registration = getattr(getattr(item, "function", None), "genairr_verifies", None)
        if not registration or not _establishes_its_claim(item):
            continue
        mechanism_id, _claim, fields = registration
        _COLLECTED.setdefault(mechanism_id, []).append((item.nodeid, tuple(fields)))


def pytest_runtest_logreport(report):
    """Record which registered tests actually passed.

    Credit was granted at collection time, so a registered test whose body
    simply failed still reported its mechanism verified. The xfail route to a
    false green was closed by inspecting markers; this closes the plain-failure
    route, which no marker reveals.
    """
    if report.when == "call" and report.passed:
        _PASSED.add(report.nodeid)


def _resolve() -> None:
    """Fold the passing registrations into the reported coverage."""
    _REGISTRATIONS.clear()
    _COVERED.clear()
    for mechanism_id, entries in _COLLECTED.items():
        for node_id, fields in entries:
            if node_id not in _PASSED:
                continue
            _REGISTRATIONS.setdefault(mechanism_id, []).append(node_id)
            _COVERED.setdefault(mechanism_id, set()).update(fields)


def _merge(registrations: dict, covered: dict) -> None:
    for mechanism_id, nodes in registrations.items():
        existing = _REGISTRATIONS.setdefault(mechanism_id, [])
        existing.extend(n for n in nodes if n not in existing)
    for mechanism_id, fields in covered.items():
        _COVERED.setdefault(mechanism_id, set()).update(fields)


# ── crossing the xdist boundary ──────────────────────────────────────────────
# Registrations are gathered during collection, which under xdist happens in
# each worker, while the report is written by the controller. CI runs -n auto,
# so without shipping them across that boundary the report would read
# "unverified" the moment anyone annotated a test.

def pytest_sessionfinish(session):
    config = session.config
    output = getattr(config, "workeroutput", None)
    distributed = getattr(config.option, "dist", "no") != "no"
    # A worker resolves its own outcomes and ships them. A serial run resolves
    # here. The controller must NOT resolve: it has no outcomes of its own and
    # resolving would clear what the workers already sent.
    if output is not None or not distributed:
        _resolve()
    if output is None:
        return                                  # controller, or no xdist
    output["genairr_registrations"] = {k: list(v) for k, v in _REGISTRATIONS.items()}
    output["genairr_covered"] = {k: sorted(v) for k, v in _COVERED.items()}


def pytest_testnodedown(node, error):
    output = getattr(node, "workeroutput", None) or {}
    _merge(output.get("genairr_registrations", {}), output.get("genairr_covered", {}))


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    if not config.getoption("--biology-coverage", default=False):
        return
    if not _REGISTRATIONS and _COLLECTED \
            and getattr(config.option, "dist", "no") == "no":
        _resolve()                              # serial fallback
    tw = terminalreporter
    from ._rust import rust_coverage
    # Only registered Python tests count toward a field being covered. Rust
    # credit is shown in its own column and cannot move a status.
    rows = coverage_rows({k: tuple(v) for k, v in _REGISTRATIONS.items()},
                         rust_by_mechanism=rust_coverage(),
                         covered_by_mechanism={k: tuple(v) for k, v in _COVERED.items()})
    colour = {"verified": {"green": True}, "partial": {"yellow": True},
              "unverified": {"red": True}, "no-truth-exported": {"purple": True}}
    tw.section("biology coverage", "=")
    stage = None
    for row in rows:
        if row.stage != stage:
            stage = row.stage
            tw.line("")
            tw.line(f"{stage}", bold=True)
        fields = f"{len(row.covered_fields)}/{len(row.truth_fields)} fields"
        asserted = f"rust-asserted:{len(row.rust_asserted)}" if row.rust_asserted else ""
        # Shown even when the status is verified: a mechanism can be verified on
        # the one borrowed field it owns while the truth that actually matters
        # is still not exported, and "verified" alone would hide that.
        gap = "  (truth gap)" if row.truth_gap else ""
        tw.write(f"  {row.mechanism_id:34s} ")
        tw.write(f"{row.status:11s}", **colour[row.status])
        tw.write(f" {fields:16s} ")
        tw.write(f"py:{len(row.python_tests)}".ljust(7))
        tw.line(f" {asserted}{gap}")
    counts = {s: sum(1 for r in rows if r.status == s)
              for s in ("verified", "partial", "unverified", "no-truth-exported")}
    tw.line("")
    tw.line(
        f"  {counts['verified']} verified, {counts['partial']} partial, "
        f"{counts['unverified']} unverified, "
        f"{counts['no-truth-exported']} with no truth exported, "
        f"of {len(rows)} mechanisms")
    tw.line("")
    tw.line("  Status comes from registered Python tests only.", bold=True)
    tw.line("  no-truth-exported means the mechanism's truth is not in the record at")
    tw.line("  all, so no test could establish it. Writing a test does not fix those;")
    tw.line("  exporting the truth does. See each mechanism's truth_gap.")
    tw.line("  rust-asserted counts engine_rs/tests/metadata_truth/ files recorded by")
    tw.line("  hand from their module docstrings. Those docstrings state intent, not")
    tw.line("  coverage, and spot-checks found credited fields the files never mention,")
    tw.line("  so they are shown but deliberately do not set a status.")
