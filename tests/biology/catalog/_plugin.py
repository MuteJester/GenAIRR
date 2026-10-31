"""The --biology-coverage report."""
from __future__ import annotations

from ._registry import coverage_rows

#: Registrations seen at collection, before any test has run.
_COLLECTED: dict[str, list[tuple[str, tuple[str, ...], str]]] = {}
#: Node ids the run observed PASS. A mechanism cannot be verified by a test
#: the same run watched fail, so credit is withheld until the outcome is known.
_PASSED: set[str] = set()
_REGISTRATIONS: dict[str, list[str]] = {}
_COVERED: dict[str, set[str]] = {}
#: mechanism id -> {truth field: the strongest basis covering it}. Separate
#: from _COVERED because "which fields were checked" and "how well" are
#: different questions and the report answers both.
_BASIS_BY_FIELD: dict[str, dict[str, str]] = {}


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
    _BASIS_BY_FIELD.clear()
    for item in items:
        registration = getattr(getattr(item, "function", None), "genairr_verifies", None)
        if not registration or not _establishes_its_claim(item):
            continue
        mechanism_id, _claim, fields, basis = registration
        _COLLECTED.setdefault(mechanism_id, []).append(
            (item.nodeid, tuple(fields), basis))


def _basis_rank(basis: str) -> int:
    from ._registry import VERIFICATION_BASES
    return VERIFICATION_BASES.index(basis)


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
    _BASIS_BY_FIELD.clear()
    for mechanism_id, entries in _COLLECTED.items():
        for node_id, fields, basis in entries:
            if node_id not in _PASSED:
                continue
            _REGISTRATIONS.setdefault(mechanism_id, []).append(node_id)
            _COVERED.setdefault(mechanism_id, set()).update(fields)
            # Per field, the strongest basis any passing registration brings.
            # A field checked once against the reference and once by a
            # presence guard is established by the former; recording the
            # weakest would report a mechanism worse than it is.
            strengths = _BASIS_BY_FIELD.setdefault(mechanism_id, {})
            for field in fields:
                current = strengths.get(field)
                if current is None or _basis_rank(basis) < _basis_rank(current):
                    strengths[field] = basis


def _merge(registrations: dict, covered: dict, basis: dict | None = None) -> None:
    for mechanism_id, nodes in registrations.items():
        existing = _REGISTRATIONS.setdefault(mechanism_id, [])
        existing.extend(n for n in nodes if n not in existing)
    for mechanism_id, fields in covered.items():
        _COVERED.setdefault(mechanism_id, set()).update(fields)
    # The basis has to cross the boundary with the rest of it. Without this
    # the controller had registrations and covered fields but no basis at all,
    # so every covered field fell to the conservative default and the report
    # CI publishes read 0 verified, 25 consistent where the same tree read 20
    # verified serially.
    for mechanism_id, strengths in (basis or {}).items():
        merged = _BASIS_BY_FIELD.setdefault(mechanism_id, {})
        for field, value in strengths.items():
            current = merged.get(field)
            if current is None or _basis_rank(value) < _basis_rank(current):
                merged[field] = value


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
    output["genairr_basis"] = {k: dict(v) for k, v in _BASIS_BY_FIELD.items()}


def pytest_testnodedown(node, error):
    output = getattr(node, "workeroutput", None) or {}
    _merge(output.get("genairr_registrations", {}),
           output.get("genairr_covered", {}),
           output.get("genairr_basis", {}))


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
                         covered_by_mechanism={k: tuple(v) for k, v in _COVERED.items()},
                         basis_by_mechanism=_BASIS_BY_FIELD)
    colour = {"verified": {"green": True}, "consistent": {"cyan": True},
              "gap-limited": {"yellow": True},
              "partial": {"yellow": True}, "unverified": {"red": True},
              "no-truth-exported": {"purple": True}}
    tw.section("biology coverage", "=")
    stage = None
    for row in rows:
        if row.stage != stage:
            stage = row.stage
            tw.line("")
            tw.line(f"{stage}", bold=True)
        # Measured against the verifiable fields: a configuration echo listed
        # beside its mechanism is not something a test could establish, so
        # counting it would make its owner permanently short of its own total.
        fields = f"{len(row.covered_fields)}/{len(row.verifiable_fields)} fields"
        if len(row.verifiable_fields) != len(row.truth_fields):
            fields += f" of {len(row.truth_fields)}"
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
              for s in ("verified", "consistent", "gap-limited", "partial",
                        "unverified", "no-truth-exported")}
    tw.line("")
    tw.line(
        f"  {counts['verified']} verified, {counts['consistent']} consistent, "
        f"{counts['gap-limited']} gap-limited, {counts['partial']} partial, "
        f"{counts['unverified']} unverified, "
        f"{counts['no-truth-exported']} with no truth exported, "
        f"of {len(rows)} mechanisms")
    tw.line("")
    tw.line("  Status comes from registered Python tests only.", bold=True)
    tw.line("  verified means every exported field is covered by a test that can")
    tw.line("  establish its VALUE: against the reference data, against an independent")
    tw.line("  oracle over the recorded provenance, or against a second run of the")
    tw.line("  library. Each registration declares which.")
    tw.line("  consistent means every field is covered and at least one only by a")
    tw.line("  comparison between two exported views of the same number, or by an")
    tw.line("  assertion a wrong value would also satisfy. The fields are reached; what")
    tw.line("  is missing is a basis for believing them. Two mechanisms sat green under")
    tw.line("  the old rule while their counters were measurably wrong.")
    tw.line("  gap-limited means every exported field is covered on a strong basis and")
    tw.line("  some of what the mechanism means is still not exported, so it is not the")
    tw.line("  same as verified.")
    tw.line("  no-truth-exported means the mechanism's truth is not in the record at")
    tw.line("  all, so no test could establish it. Writing a test does not fix those;")
    tw.line("  exporting the truth does. See each mechanism's truth_gap.")
    tw.line("  rust-asserted counts engine_rs/tests/metadata_truth/ files recorded by")
    tw.line("  hand from their module docstrings. Those docstrings state intent, not")
    tw.line("  coverage, and spot-checks found credited fields the files never mention,")
    tw.line("  so they are shown but deliberately do not set a status.")
