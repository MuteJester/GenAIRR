"""Contract pins for the Docs Website Information
Architecture audit.

Companion to
[`docs/docs_website_audit.md`](../docs/docs_website_audit.md).

Pin set:

- ``pin_scaffold_*`` — current state inventory of the
  three parallel doc trees (`website/` live site,
  `docs/*.md` maintainer audits, `_old_docs/` abandoned
  earlier attempt), the deployment workflow, and the
  README's cross-reference shape.
- ``pin_present_*`` — broken-URL verification: the
  README references URLs that match the abandoned
  `_old_docs/` MkDocs scheme, NOT the live
  `website/*.html` scheme. The check is offline-safe —
  it parses the README and compares URL shapes against
  the live filesystem, no network required.
- ``pin_scaffold_*`` — content gaps: only 1 of 5
  estimators is mentioned in the live site;
  `website/reference.html` is hand-curated; the 35-symbol
  top-level export is not enumerated.
- ``pin_absence_*`` — framework configs: no `mkdocs.yml`
  / Sphinx `conf.py` / `docusaurus.config.js` / Jekyll
  `_config.yml` in the live tree today.

**Pre-flight verdict (audit §3): stop-and-report
triggered (mild).** A docs site exists; do not propose
blind replacement. The README's external URL scheme is
broken; that's the critical-path fix.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest


_REPO_ROOT = Path(__file__).resolve().parent.parent
_AUDIT_DOC = _REPO_ROOT / "docs" / "docs_website_audit.md"
_README = _REPO_ROOT / "README.md"
_WEBSITE = _REPO_ROOT / "website"
_DOCS = _REPO_ROOT / "docs"
_OLD_DOCS = _REPO_ROOT / "_old_docs"
_DEPLOY_WORKFLOW = _REPO_ROOT / ".github" / "workflows" / "deploy-docs.yml"


# ──────────────────────────────────────────────────────────────────
# A. pin_scaffold_* — current state inventory
# ──────────────────────────────────────────────────────────────────

# --- skip whole module if docs/ tree is absent ---
# The docs/ directory holds contributor-only audit + design markdowns
# and is gitignored. CI without docs/ skips this entire test file rather
# than report dozens of irrelevant failures.
_DOCS_DIR = _REPO_ROOT / "docs"
pytestmark = pytest.mark.skipif(
    not _DOCS_DIR.is_dir(),
    reason="docs/ is contributor-only and not present in this checkout",
)


def test_pin_scaffold_website_dir_exists_with_handwritten_html() -> None:
    """`website/` exists with handwritten HTML files. No
    build step intervenes — the deployment workflow
    uploads it verbatim. Pinned at filesystem level."""
    assert _WEBSITE.is_dir(), f"{_WEBSITE} missing"
    html_files = sorted(_WEBSITE.glob("*.html"))
    assert len(html_files) >= 25, (
        f"website/ has only {len(html_files)} HTML files; audit "
        f"§1.2 documented ~30 pages — content was deleted?"
    )
    # Hub pages present.
    for hub in ("index.html", "learn.html", "guides.html",
                "concepts.html", "reference.html"):
        assert (_WEBSITE / hub).is_file(), (
            f"website/{hub} missing — site navigation hub regressed"
        )
    # Lesson series present.
    for n in range(1, 6):
        assert (_WEBSITE / f"lesson-{n}.html").is_file(), (
            f"website/lesson-{n}.html missing — lesson series regressed"
        )


def test_pin_post_cutover_deploy_docs_workflow_builds_mkdocs_and_uploads_site() -> None:
    """Post-Phase-6-cutover: `.github/workflows/deploy-docs.yml`
    builds the MkDocs Material site under `site_docs/` with
    `make docs-build` (which runs `mkdocs build --strict`)
    and publishes the generated `site/` directory to GitHub
    Pages. The workflow's `path: site` line plus the
    `make docs-build` step are the load-bearing pins."""
    assert _DEPLOY_WORKFLOW.is_file(), (
        f"{_DEPLOY_WORKFLOW} missing — docs deployment workflow gone"
    )
    text = _DEPLOY_WORKFLOW.read_text(encoding="utf-8")
    # Build step present.
    assert "make docs-build" in text, (
        "deploy-docs.yml no longer builds the MkDocs site — cutover "
        "Phase 6's `make docs-build` step regressed"
    )
    # Upload target is the MkDocs build output, not website/.
    assert "path: site" in text, (
        "deploy-docs.yml no longer uploads `site/` — Phase 6 "
        "cutover state regressed"
    )
    # Sanity: no Sphinx / Docusaurus / other frameworks snuck in.
    assert "sphinx-build" not in text
    assert "docusaurus" not in text


def test_pin_scaffold_docs_dir_carries_audit_design_md_files() -> None:
    """`docs/` carries ≥ 35 markdown files (audits +
    designs + hubs). The audit catalogued 38 at write
    time; the pin is intentionally loose so adding new
    audits doesn't fail this test."""
    md_files = sorted(_DOCS.glob("*.md"))
    assert len(md_files) >= 35, (
        f"docs/ has only {len(md_files)} markdown files; audit §1.3 "
        f"documented 38 — audits/designs were deleted?"
    )
    # Spot-check canonical files that the validation matrix relies on.
    for name in (
        "validation_matrix.md",
        "engine_architecture.md",
        "adding_a_pass.md",
        "reference_cartridge.md",
        "airr_record_validator.md",
    ):
        assert (_DOCS / name).is_file(), (
            f"docs/{name} missing — contributor entry point regressed"
        )


def test_pin_scaffold_old_docs_dir_exists_as_abandoned_earlier_attempt() -> None:
    """`_old_docs/` exists as the deprecated earlier docs
    tree. Renamed with leading underscore signals
    deliberate abandonment. The audit recommends keeping
    it as a salvage source (tutorials + MkDocs
    requirements file)."""
    assert _OLD_DOCS.is_dir(), (
        f"{_OLD_DOCS} missing — abandoned-docs evidence trail gone"
    )
    # MkDocs Material breadcrumb.
    req_file = _OLD_DOCS / "doc_requirements.txt"
    assert req_file.is_file(), (
        "_old_docs/doc_requirements.txt missing — MkDocs Material "
        "abandoned-migration evidence regressed; audit §6.1's "
        "framework-history inference rests on this file"
    )
    req_content = req_file.read_text(encoding="utf-8")
    assert "mkdocs-material" in req_content
    assert "mkdocs-jupyter" in req_content
    # Sphinx breadcrumb (the Makefile still claims Sphinx).
    makefile = _OLD_DOCS / "Makefile"
    assert makefile.is_file()
    assert "sphinx-build" in makefile.read_text(encoding="utf-8")
    # 5 Jupyter tutorials.
    tutorials = sorted((_OLD_DOCS / "tutorials").glob("*.ipynb"))
    assert len(tutorials) >= 5, (
        f"_old_docs/tutorials/ has only {len(tutorials)} notebooks; "
        f"audit §1.4 documented 5"
    )


def test_pin_scaffold_docs_superpowers_subdir_holds_session_artifacts() -> None:
    """`docs/superpowers/plans/` carries Claude-session
    planning artefacts, NOT user-facing documentation.
    Pinned defensively so a future "is this a doc?" sweep
    treats it correctly."""
    superpowers = _DOCS / "superpowers"
    if not superpowers.is_dir():
        # If somebody moves it, that's fine — pin only the docs
        # subdirs that exist today.
        return
    plans = superpowers / "plans"
    if plans.is_dir():
        plans_files = sorted(plans.glob("*.md"))
        # The presence of a date-prefixed planning markdown is
        # the canonical evidence (e.g. 2026-05-18-mcp-redesign-v2.md).
        date_prefixed = [
            f for f in plans_files
            if re.match(r"^\d{4}-\d{2}-\d{2}-", f.name)
        ]
        assert date_prefixed, (
            "docs/superpowers/plans/ has no date-prefixed planning "
            "files — the session-artefact discipline broke; this "
            "directory's contents may have been recategorised as "
            "docs without an audit-doc update"
        )


def test_pin_scaffold_docs_build_subdir_holds_wheel_artefacts_not_docs() -> None:
    """`docs/build/` holds Python wheel build artefacts
    (a side-effect of `python -m build` choosing
    `docs/build/` for its temp output). NOT documentation.
    Pinned so a future audit doesn't mistake it for
    rendered docs output."""
    build_dir = _DOCS / "build"
    if not build_dir.is_dir():
        return  # build dir may not exist on a fresh clone
    # Wheel build leaves these distinctive subdirs.
    for marker in ("bdist.linux-x86_64", "lib.linux-x86_64-cpython-312",
                   "temp.linux-x86_64-cpython-312"):
        # Tolerate platform-specific names — just check the prefix.
        prefix = marker.split(".")[0]
        matches = list(build_dir.glob(f"{prefix}.*"))
        if matches:
            # At least one wheel build artefact present → confirmation.
            return
    # No wheel artefacts found — maybe docs/build/ is now used for
    # something else. Surface for review.
    pytest.skip(
        "docs/build/ exists but contains no recognised wheel "
        "artefacts; verify it's still a wheel build target"
    )


# ──────────────────────────────────────────────────────────────────
# B. pin_scaffold_* — README cross-reference shape
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_readme_links_many_md_files_in_docs_dir() -> None:
    """The README links to many `docs/*.md` files in
    normal prose. These render on GitHub (where the
    README is viewed pre-install) but are NOT part of
    the live `website/` site — the maintainer-leak
    problem the audit §2.2 names. The audit time count
    was 15 distinct `docs/*.md` markdown links."""
    readme_text = _README.read_text(encoding="utf-8")
    docs_links = re.findall(r"\(docs/[^)]+\.md[^)]*\)", readme_text)
    unique_links = sorted(set(docs_links))
    assert len(unique_links) >= 10, (
        f"README links to only {len(unique_links)} distinct "
        f"docs/*.md files; audit §2.2 documented 15 — verify "
        f"the README was not silently stripped of internal "
        f"references"
    )


def test_pin_scaffold_readme_advertises_mutejester_github_io_as_docs_home() -> None:
    """The README's "📖 Documentation" link points at
    `mutejester.github.io/GenAIRR/` — the GitHub Pages
    deployment target. The slugged URLs in the prose
    point at a DIFFERENT URL scheme (see the broken-URL
    pin)."""
    readme_text = _README.read_text(encoding="utf-8")
    assert "mutejester.github.io/GenAIRR" in readme_text, (
        "README no longer advertises the GitHub Pages URL — site "
        "homepage reference regressed"
    )


# ──────────────────────────────────────────────────────────────────
# C. pin_present_* — broken-URL verification (offline-safe)
# ──────────────────────────────────────────────────────────────────


def test_pin_present_readme_hosted_doc_urls_resolve_to_live_website_pages() -> None:
    """**Slice 1 (README Production URL Fix) — flipped
    from the prior broken-URL stop-and-report pin.**

    Every `mutejester.github.io/GenAIRR/<path>` URL in the
    README MUST correspond to either:

    - The site root (`mutejester.github.io/GenAIRR/` with
      no path component) — the homepage link.
    - An existing `website/<path>` file — the flat-HTML
      page the deployment workflow uploads.

    The check is offline-safe: parses the README,
    extracts every hosted-doc URL, and verifies the
    filesystem mapping. No network calls — the test runs
    against the source tree.

    Pre-slice this pin asserted the BROKEN state (slugged
    MkDocs-style URLs that 404'd in production); the
    Slice 1 implementation rewrote them all to flat
    `<page>.html` URLs. If a regression re-introduces a
    URL that doesn't resolve, this pin fires with the
    exact missing path."""
    readme_text = _README.read_text(encoding="utf-8")
    url_pattern = re.compile(
        r"https://mutejester\.github\.io/GenAIRR/[^)\s\"]*"
    )
    urls = sorted(set(url_pattern.findall(readme_text)))
    assert urls, (
        "README no longer contains any mutejester.github.io/GenAIRR/ "
        "URLs — verify the docs section was not silently stripped"
    )
    missing: list[tuple[str, str]] = []
    for url in urls:
        # Strip the hosting prefix to get the path component.
        path = url[len("https://mutejester.github.io/GenAIRR/"):]
        if not path:
            # Root URL — homepage link, always OK.
            continue
        # Strip query string / fragment if present.
        path = path.split("?", 1)[0].split("#", 1)[0]
        target = _WEBSITE / path
        if not target.is_file():
            missing.append((url, str(target)))
    assert not missing, (
        "README hosted-doc URLs do not resolve to live website "
        "pages:\n"
        + "\n".join(
            f"  {url} → {target} (missing)"
            for url, target in missing
        )
        + "\n\nEither fix the URL to point at an existing "
        "website/*.html page or use the site root "
        "(https://mutejester.github.io/GenAIRR/)."
    )


def test_pin_present_readme_contains_no_slugged_mkdocs_urls() -> None:
    """**Slice 1 regression guard.** The README MUST NOT
    re-introduce slugged URLs at
    `mutejester.github.io/GenAIRR/docs/...` (the
    abandoned MkDocs scheme that 404'd in production).
    Pinned defensively so a future doc edit that copies
    URLs from the abandoned `_old_docs/` migration
    surfaces here."""
    readme_text = _README.read_text(encoding="utf-8")
    slugged = re.findall(
        r"mutejester\.github\.io/GenAIRR/docs/[^)\s\"]+",
        readme_text,
    )
    assert not slugged, (
        f"README contains {len(slugged)} slugged "
        f"mutejester.github.io/GenAIRR/docs/... URLs — these "
        f"match the abandoned MkDocs scheme and 404 in "
        f"production. Either point at a flat <page>.html or use "
        f"the site root. First few: {slugged[:3]}"
    )


def test_pin_scaffold_live_website_uses_flat_html_url_scheme() -> None:
    """The live site uses flat `*.html` paths. Pinned so a
    structural migration that introduces subdirectories
    (like `lessons/01-intro.html`) surfaces here for
    explicit review of URL-stability impact."""
    html_files = sorted(_WEBSITE.glob("*.html"))
    assert html_files, "website/ has no HTML files at all"
    # No subdirectories with HTML files (modulo asset
    # subdirectories like static/, assets/, etc.).
    for sub_html in _WEBSITE.glob("**/*.html"):
        if sub_html.parent != _WEBSITE:
            pytest.fail(
                f"{sub_html} lives in a subdirectory — site URL "
                f"scheme has changed from flat to nested; audit "
                f"§3.1's URL-scheme analysis regressed"
            )


# ──────────────────────────────────────────────────────────────────
# D. pin_scaffold_* — content gaps (estimator + API reference)
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_only_one_of_five_estimators_referenced_in_live_site() -> None:
    """The 5 cartridge-builder estimators that landed in
    this session are:

    - `estimate_allele_usage`
    - `estimate_trim_distributions`
    - `estimate_np_length_distributions`
    - `estimate_np_base_model`
    - `estimate_p_nucleotide_lengths`

    Only `estimate_allele_usage` is mentioned in any
    `website/*.html` page today. The other 4 are
    completely undocumented user-side. Pinned so a
    future migration that closes the gap flips this
    pin (the assertion would then need updating to "5
    of 5")."""
    estimators = (
        "estimate_allele_usage",
        "estimate_trim_distributions",
        "estimate_np_length_distributions",
        "estimate_np_base_model",
        "estimate_p_nucleotide_lengths",
    )
    referenced = []
    for html_file in _WEBSITE.glob("*.html"):
        text = html_file.read_text(encoding="utf-8", errors="replace")
        for est in estimators:
            if est in text:
                referenced.append(est)
    referenced_unique = sorted(set(referenced))
    # Audit time: at most 1 estimator referenced (allele_usage). The
    # pin tolerates 0–2 to absorb minor doc updates without flapping;
    # crossing 3 indicates significant new content the audit didn't
    # anticipate — surface for explicit review.
    assert len(referenced_unique) <= 2, (
        f"{len(referenced_unique)} estimators now referenced in "
        f"website/*.html: {referenced_unique!r}. Audit §5.1 said "
        f"only `estimate_allele_usage` was documented; a closing of "
        f"the gap is welcome but should trigger an audit-doc update "
        f"so this pin can flip cleanly."
    )


def test_pin_scaffold_top_level_api_export_count_exceeds_thirty() -> None:
    """The top-level `GenAIRR.__all__` exports 35 symbols
    at audit time. The live site's `reference.html` does
    NOT enumerate them. Pin the count so the implementation
    slice that adds the auto-generated reference page knows
    how many symbols it needs to cover."""
    import GenAIRR as ga

    assert len(ga.__all__) >= 30, (
        f"GenAIRR.__all__ has only {len(ga.__all__)} symbols; "
        f"audit §1.4 documented 35 — verify the public API has not "
        f"contracted unexpectedly"
    )


def test_pin_scaffold_reference_html_is_hand_curated_not_auto_generated() -> None:
    """`website/reference.html` is hand-curated. There's
    no marker indicating it was auto-generated from
    `__all__` + docstrings. Pinned so a future
    auto-generation slice surfaces here for explicit
    review (the implementation should add a clear
    auto-generated marker)."""
    ref = _WEBSITE / "reference.html"
    if not ref.is_file():
        return
    text = ref.read_text(encoding="utf-8")
    # Heuristic: auto-generated pages typically carry a marker like
    # "GENERATED — DO NOT EDIT" or reference inspect.signature output.
    for marker in ("GENERATED — DO NOT EDIT",
                   "inspect.signature",
                   "auto-generated from __all__"):
        if marker in text:
            pytest.fail(
                f"reference.html now contains the auto-generated "
                f"marker {marker!r} — the API-reference gap is being "
                f"closed; audit §3.2 noted this would happen on "
                f"Path B / late Path A. Update the audit doc + flip "
                f"this pin."
            )


# ──────────────────────────────────────────────────────────────────
# E. pin_absence_* — framework configs
# ──────────────────────────────────────────────────────────────────


def test_pin_present_mkdocs_yml_landed_at_repo_root() -> None:
    """Post-Phase-0 — `mkdocs.yml` landed at the repo root
    as the MkDocs Material migration scaffold. Flipped
    from the prior absence pin in lockstep with the
    `docs_mkdocs_migration_plan.md` Phase 0 deliverable.

    The mkdocs.yml content + buildability is owned by the
    migration-contract test file
    (`tests/test_docs_mkdocs_migration_contract.py`);
    this pin is the lockstep flip on the docs-website-
    audit contract."""
    assert (_REPO_ROOT / "mkdocs.yml").is_file(), (
        "mkdocs.yml at repo root missing — Phase 0 scaffold "
        "regressed; see docs/docs_mkdocs_migration_plan.md"
    )
    # The migration plan's source-directory choice (site_docs/, not
    # docs/) means no `docs/mkdocs.yml` and no alternate yaml ext.
    for forbidden in (_REPO_ROOT / "mkdocs.yaml",
                      _DOCS / "mkdocs.yml"):
        assert not forbidden.is_file(), (
            f"{forbidden} exists — multiple mkdocs configs would "
            f"create framework ambiguity; consolidate on the repo-"
            f"root mkdocs.yml per the migration plan"
        )


def test_pin_absence_no_sphinx_conf_py_today() -> None:
    """No Sphinx `conf.py` in the live tree. The
    `_old_docs/Makefile` still claims Sphinx but the
    matching `source/conf.py` was deleted long ago."""
    for candidate in (_REPO_ROOT / "conf.py",
                      _DOCS / "conf.py",
                      _DOCS / "source" / "conf.py",
                      _OLD_DOCS / "source" / "conf.py",
                      _OLD_DOCS / "conf.py"):
        if candidate.is_file():
            pytest.fail(
                f"{candidate} exists — Sphinx is now configured; "
                f"audit §6.1's framework-history finding regressed"
            )


def test_pin_absence_no_docusaurus_or_jekyll_config_today() -> None:
    """No Docusaurus / Jekyll / Hugo / Eleventy config.
    Audit §6 considered only Sphinx and MkDocs as the
    realistic candidates given the existing breadcrumbs;
    pinning the absence of other frameworks frames the
    choice."""
    for candidate in (
        _REPO_ROOT / "docusaurus.config.js",
        _REPO_ROOT / "_config.yml",
        _REPO_ROOT / "config.toml",  # Hugo
        _REPO_ROOT / ".eleventy.js",
    ):
        assert not candidate.is_file(), (
            f"{candidate} exists — a third docs framework joined the "
            f"mix; audit needs an update to evaluate it as a Path C"
        )


# ──────────────────────────────────────────────────────────────────
# F. Doc anchor
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_audit_doc_exists_and_references_contract() -> None:
    """The audit doc exists and references the contract
    file by name; section structure intact."""
    if not _AUDIT_DOC.exists():
        import pytest
        pytest.skip("docs/ is contributor-only; audit doc not present in this checkout")
    doc = _AUDIT_DOC.read_text(encoding="utf-8")
    assert "test_docs_website_contract.py" in doc, (
        "audit doc no longer references the contract file"
    )
    for marker in (
        "## 1. Q1",
        "## 3. Q3",
        "## 6. Q6",
        "## 9. Implementation order",
        "## 10. Summary table",
    ):
        assert marker in doc, f"audit doc missing section marker {marker!r}"
