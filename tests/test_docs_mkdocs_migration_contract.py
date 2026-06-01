"""Contract pins for the MkDocs Migration Plan
(Phases 0 + 1, documented in
[`docs/docs_mkdocs_migration_plan.md`](../docs/docs_mkdocs_migration_plan.md)).

Pin set:

- ``pin_scaffold_*`` — Phase 0 scaffold files exist
  (`mkdocs.yml`, `site_docs/`, the starter pages).
- ``pin_scaffold_*`` — the migration plan doc carries
  the load-bearing sections referenced by other slices.
- ``pin_present_*`` — **redirect-map completeness**:
  every existing `website/*.html` page has a row in the
  migration plan's §3.1 redirect table, AND every
  hosted-doc URL the README still references appears in
  that same table. This is the load-bearing URL-
  preservation contract.
- ``pin_present_*`` — Phase 1 build tooling:
  `docs/requirements-docs.txt` exists with the four
  packages; plugins are active in `mkdocs.yml` (no
  longer commented out); every nav target under
  `site_docs/` exists as a real file; if the deps are
  installed, `mkdocs build --strict` succeeds.
- ``pin_absence_*`` — the cutover hasn't happened yet:
  `website/` content remains; `.github/workflows/deploy-docs.yml`
  still targets `website/`.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


_REPO_ROOT = Path(__file__).resolve().parent.parent
_MIGRATION_PLAN = _REPO_ROOT / "docs" / "docs_mkdocs_migration_plan.md"
_MKDOCS_YML = _REPO_ROOT / "mkdocs.yml"
_SITE_DOCS = _REPO_ROOT / "site_docs"
_WEBSITE = _REPO_ROOT / "website"
_README = _REPO_ROOT / "README.md"
_DEPLOY_WORKFLOW = _REPO_ROOT / ".github" / "workflows" / "deploy-docs.yml"


# ──────────────────────────────────────────────────────────────────
# A. pin_scaffold_* — scaffold files exist
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_mkdocs_yml_exists_at_repo_root() -> None:
    """`mkdocs.yml` exists at the repo root and carries
    the minimum-required Material theme + nav skeleton.
    Pinned shape, not content."""
    assert _MKDOCS_YML.is_file(), (
        f"{_MKDOCS_YML} missing — Phase 0 scaffold regressed"
    )
    text = _MKDOCS_YML.read_text(encoding="utf-8")
    for required in (
        "site_name: GenAIRR",
        "docs_dir: site_docs",
        "theme:",
        "name: material",
        "nav:",
    ):
        assert required in text, (
            f"mkdocs.yml missing required line {required!r} — "
            f"scaffold shape regressed"
        )


def test_pin_scaffold_site_docs_directory_exists_with_starter_files() -> None:
    """`site_docs/` exists with the 4 starter files per
    the user brief (Slice 3 scope)."""
    assert _SITE_DOCS.is_dir(), f"{_SITE_DOCS} missing"
    for required_page in (
        "index.md",
        "getting-started/quick-start.md",
        "reference/index.md",
        "architecture/index.md",
    ):
        target = _SITE_DOCS / required_page
        assert target.is_file(), (
            f"site_docs/{required_page} missing — Phase 0 scaffold "
            f"regressed (user brief required this starter file)"
        )


def test_pin_scaffold_site_docs_index_md_is_non_empty() -> None:
    """`site_docs/index.md` is the homepage — non-empty +
    carries the project name. Phase 0 also required a
    migration-status banner; Phase 1.5 replaced that banner
    with the polished landing-page IA (see
    `test_pin_present_homepage_no_longer_a_placeholder`)."""
    index = _SITE_DOCS / "index.md"
    text = index.read_text(encoding="utf-8")
    assert "GenAIRR" in text, "homepage doesn't mention the project name"
    # The body length is a coarse proxy for "non-placeholder".
    assert len(text) > 500, (
        f"homepage is only {len(text)} bytes — looks like a stub. "
        "Phase 1.5 expects a polished landing page."
    )


# ──────────────────────────────────────────────────────────────────
# B. pin_scaffold_* — migration plan carries load-bearing sections
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_migration_plan_exists_with_load_bearing_sections() -> None:
    """`docs/docs_mkdocs_migration_plan.md` exists and
    carries the section markers other contract pins
    reference (§3 redirect map, §4 ownership table)."""
    assert _MIGRATION_PLAN.is_file(), (
        f"{_MIGRATION_PLAN} missing — Phase 0 deliverable regressed"
    )
    text = _MIGRATION_PLAN.read_text(encoding="utf-8")
    for marker in (
        "## 1. Current-site inventory",
        "## 2. Target nav structure",
        "## 3. URL redirect map",
        "## 4. Content ownership table",
        "## 5. Phased migration plan",
        "## 6. Redirect strategy",
    ):
        assert marker in text, (
            f"migration plan missing section {marker!r} — other "
            f"contract pins depend on this section's content"
        )


# ──────────────────────────────────────────────────────────────────
# C. pin_present_* — redirect-map completeness
# ──────────────────────────────────────────────────────────────────


def _extract_redirect_pages_from_plan() -> set[str]:
    """Parse the migration plan §3.1 redirect table for
    every `website/*.html` filename referenced in the
    'Source `website/` file' column. Returns the set of
    filenames (just the basename, not the full path)."""
    text = _MIGRATION_PLAN.read_text(encoding="utf-8")
    # The redirect-table rows have the shape:
    #   | `/lesson-1.html` | `getting-started/...` | `lesson-1.html` |
    # We extract the rightmost backtick-wrapped HTML filename per row.
    pattern = re.compile(r"`([a-z0-9_\-]+\.html)`")
    return set(pattern.findall(text))


def test_pin_present_redirect_map_covers_every_website_html_page() -> None:
    """**Load-bearing URL-preservation contract.**

    Every `website/*.html` page MUST appear in the
    migration plan's §3.1 redirect map. The plan is the
    machine-verifiable target every subsequent migration
    phase uses to keep URLs alive.

    If a page is added to `website/` without a
    corresponding redirect-map row, this pin fires and
    the URL would silently break at Phase 6 cutover."""
    website_pages = {p.name for p in _WEBSITE.glob("*.html")}
    redirect_pages = _extract_redirect_pages_from_plan()
    missing = website_pages - redirect_pages
    assert not missing, (
        f"website/ pages missing from migration plan §3.1 redirect "
        f"map: {sorted(missing)}\n"
        f"Either add them to the plan, or document why they're "
        f"excluded (e.g. internal-only pages not to be migrated)."
    )


def test_pin_present_redirect_map_covers_every_readme_hosted_doc_url() -> None:
    """**README/redirect-map alignment.**

    Every `mutejester.github.io/GenAIRR/<page>.html` URL
    the README still references MUST also appear in the
    migration plan's §3.1 redirect map. Otherwise the
    README link would 404 after the Phase 6 cutover.

    The Slice 1 contract test
    (`tests/test_docs_website_contract.py::test_pin_present_readme_hosted_doc_urls_resolve_to_live_website_pages`)
    enforces filesystem-level coverage today; this pin
    extends the coverage to the migration plan."""
    readme_text = _README.read_text(encoding="utf-8")
    url_pattern = re.compile(
        r"https://mutejester\.github\.io/GenAIRR/([a-z0-9_\-]+\.html)"
    )
    readme_pages = set(url_pattern.findall(readme_text))
    redirect_pages = _extract_redirect_pages_from_plan()
    missing = readme_pages - redirect_pages
    assert not missing, (
        f"README hosted-doc URLs not covered by migration plan §3.1: "
        f"{sorted(missing)}\n"
        f"Every README-cited live page must have a redirect-map row "
        f"so Phase 6 cutover preserves the URL contract Slice 1 "
        f"committed to."
    )


# ──────────────────────────────────────────────────────────────────
# D. pin_scaffold_* — content-ownership table is complete
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_ownership_table_covers_every_docs_md_file() -> None:
    """Every `docs/*.md` file appears in the migration
    plan's §4 ownership table (under one of the four
    classification buckets). This is how Phase 2-4 know
    which files to migrate vs leave alone."""
    plan_text = _MIGRATION_PLAN.read_text(encoding="utf-8")
    docs_md_files = {p.name for p in (_REPO_ROOT / "docs").glob("*.md")}
    missing: list[str] = []
    for name in sorted(docs_md_files):
        if name not in plan_text:
            missing.append(name)
    assert not missing, (
        f"docs/*.md files missing from migration plan §4 "
        f"ownership table: {missing}\n"
        f"Every docs/ markdown file must be classified into one of "
        f"the 4 buckets (keep_as_user_guide / "
        f"rewrite_from_audit_into_user_guide / contributor_only / "
        f"archive) so phase ordering is deterministic."
    )


# ──────────────────────────────────────────────────────────────────
# E. pin_absence_* — cutover hasn't happened yet
# ──────────────────────────────────────────────────────────────────


def test_pin_post_cutover_website_dir_preserved_as_rollback_source() -> None:
    """Post-Phase-6-cutover: the hand-rolled `website/`
    directory stays in the repo as the rollback source and
    historical reference, but is no longer deployed by CI.
    The deploy workflow now builds and uploads `site/`
    (the MkDocs output); the `mkdocs-redirects` plugin
    preserves every live URL from the old site."""
    assert _WEBSITE.is_dir(), (
        "website/ was deleted post-cutover — Phase 6's "
        "rollback-source discipline regressed. The directory is "
        "preserved for emergency rollback and historical reference."
    )
    # And the deployment workflow now uploads the MkDocs build output.
    deploy_text = _DEPLOY_WORKFLOW.read_text(encoding="utf-8")
    assert "path: site" in deploy_text, (
        "deploy-docs.yml no longer uploads `site/` (the MkDocs "
        "build output) — Phase 6 cutover state regressed."
    )


# Phase 1 absence-pin removed: the absence-of-mkdocs check is
# replaced by the `pin_present_mkdocs_build_strict_succeeds` pin
# below, which skips cleanly when deps aren't installed.


# Phase 0 redirect-block-preview pin REMOVED: superseded by
# `pin_present_redirect_sources_map_to_existing_website_pages`
# (Phase 1 — the redirect block is now active, not commented).


# ──────────────────────────────────────────────────────────────────
# G. Doc anchor
# ──────────────────────────────────────────────────────────────────


def test_pin_scaffold_migration_plan_references_contract_file() -> None:
    """Cross-reference: the migration plan should mention
    the contract test file by name so a reader can find
    the machine-verifiable pins."""
    text = _MIGRATION_PLAN.read_text(encoding="utf-8")
    assert "test_docs_mkdocs_migration_contract.py" in text, (
        "migration plan no longer references the contract file"
    )


# ──────────────────────────────────────────────────────────────────
# H. pin_present_* — Phase 1 build tooling
# ──────────────────────────────────────────────────────────────────

_DOCS_REQ = _REPO_ROOT / "docs" / "requirements-docs.txt"
_REQUIRED_DOCS_PACKAGES = (
    "mkdocs-material",
    "mkdocs-redirects",
    "mkdocstrings",
    "mkdocs-jupyter",
)


def test_pin_present_docs_requirements_file_exists_with_four_packages() -> None:
    """Phase 1 — `docs/requirements-docs.txt` ships the
    four docs-build packages. Per the migration plan §6.2,
    these are docs-only deps (NOT runtime). The Phase 1
    contract is: `pip install -r docs/requirements-docs.txt`
    leaves the env with enough to run `mkdocs build --strict`."""
    assert _DOCS_REQ.is_file(), (
        f"{_DOCS_REQ} missing — Phase 1 build-tooling deliverable "
        f"regressed"
    )
    text = _DOCS_REQ.read_text(encoding="utf-8")
    for pkg in _REQUIRED_DOCS_PACKAGES:
        assert pkg in text, (
            f"docs/requirements-docs.txt missing {pkg!r} — Phase 1 "
            f"package list regressed (per migration plan §6.2)"
        )


def test_pin_present_mkdocs_yml_plugins_block_is_active_not_commented() -> None:
    """Phase 1 — the `plugins:` block in `mkdocs.yml` is
    actually active (top-level YAML key), NOT a commented
    preview. Pinned so a future edit that re-comments the
    block accidentally surfaces here."""
    yml_text = _MKDOCS_YML.read_text(encoding="utf-8")
    # Look for an actual top-level `plugins:` key (not "# plugins:").
    lines = yml_text.splitlines()
    has_active_plugins_block = any(
        line.rstrip() == "plugins:" for line in lines
    )
    assert has_active_plugins_block, (
        "mkdocs.yml has no active top-level `plugins:` key — the "
        "block is still commented out. Phase 1 requires the plugins "
        "to be enabled."
    )
    # The three Phase 1 plugin names must appear under the active
    # block (not inside a commented section).
    for plugin in ("- search", "- redirects:", "- mkdocstrings:"):
        # Find the first non-commented occurrence.
        found_active = False
        for line in lines:
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            if plugin in line:
                found_active = True
                break
        assert found_active, (
            f"mkdocs.yml plugins block missing active {plugin!r} — "
            f"Phase 1 plugin enablement regressed"
        )


def test_pin_present_redirect_sources_map_to_existing_website_pages() -> None:
    """Phase 1 — every LEFT-side entry in
    `mkdocs.yml`'s `redirect_maps:` block corresponds to
    an existing `website/<name>.html` file.

    The plugin uses `.md` on the LEFT side as a fake
    source path (because the plugin warns on non-`.md`
    extensions). With `use_directory_urls: false`, MkDocs
    translates each LEFT entry `<name>.md` to a URL
    `<name>.html` — which MUST correspond to a real
    `website/<name>.html` so the live-URL contract Slice 1
    committed to stays preserved."""
    yml_text = _MKDOCS_YML.read_text(encoding="utf-8")
    # Find every quoted-string-on-LEFT redirect entry.
    pattern = re.compile(r'^\s+"([a-z0-9_\-/]+\.md)":', re.MULTILINE)
    redirect_sources = pattern.findall(yml_text)
    assert redirect_sources, (
        "no redirect_maps entries found in mkdocs.yml — Phase 1 "
        "URL-preservation contract regressed"
    )
    missing: list[tuple[str, str]] = []
    for src_md in redirect_sources:
        # Translate `<name>.md` → `<name>.html` per the plugin's
        # URL-rendering logic with use_directory_urls: false.
        html_name = src_md[:-3] + ".html"
        target = _WEBSITE / html_name
        if not target.is_file():
            missing.append((src_md, str(target)))
    assert not missing, (
        "redirect_maps source paths don't correspond to existing "
        "website/*.html files:\n"
        + "\n".join(
            f"  {src_md} (would emit /{src_md[:-3]}.html) → {target} (missing)"
            for src_md, target in missing
        )
        + "\n\nEvery LEFT entry must map to an existing website page "
        "so the live-URL contract holds post-cutover."
    )


def test_pin_present_use_directory_urls_is_false_for_html_url_contract() -> None:
    """Phase 1 — `use_directory_urls: false` is required
    for the `<name>.html` URL contract. Otherwise MkDocs
    would emit `<name>/` clean-URLs and the
    `mkdocs-redirects` plugin would generate redirects at
    the wrong paths.

    Pinned so a future config edit that re-enables clean
    URLs surfaces here — that change would silently break
    every live README link."""
    yml_text = _MKDOCS_YML.read_text(encoding="utf-8")
    assert "use_directory_urls: false" in yml_text, (
        "mkdocs.yml no longer sets `use_directory_urls: false` — "
        "the URL contract requires `<name>.html` paths, not clean "
        "URLs. Re-enabling directory URLs would break every Slice 1 "
        "README link at cutover."
    )


def test_pin_present_every_nav_target_exists_under_site_docs() -> None:
    """Phase 1 — every file referenced in `mkdocs.yml`'s
    `nav:` block exists under `site_docs/`. Pinned so a
    future nav entry without a matching file surfaces
    here (and would fail `mkdocs build --strict`)."""
    yml_text = _MKDOCS_YML.read_text(encoding="utf-8")
    # Nav entries take the shape `<text>: path/to/file.md` or just
    # `- path/to/file.md`. Match any .md path appearing after a
    # colon or hyphen.
    pattern = re.compile(r"[:\-]\s+([a-z0-9_\-/]+\.md)\s*$", re.MULTILINE)
    candidates = pattern.findall(yml_text)
    # Filter to entries that look like nav targets (under site_docs/),
    # not redirect_maps entries (which are also .md but live under
    # the `redirect_maps:` block).
    nav_targets = set()
    in_nav = False
    in_redirects = False
    for line in yml_text.splitlines():
        stripped = line.strip()
        if stripped == "nav:":
            in_nav = True
            in_redirects = False
            continue
        if "redirect_maps:" in stripped:
            in_redirects = True
            in_nav = False
            continue
        # Top-level keys end both nav and redirects context.
        if line and not line.startswith(" ") and not line.startswith("#"):
            in_nav = False
            in_redirects = False
            continue
        if in_nav:
            # Extract any .md path from the line.
            for cand in re.findall(r"([a-z0-9_\-/]+\.md)", line):
                # Skip redirect entries (they live under
                # redirect_maps not nav).
                nav_targets.add(cand)
    missing = []
    for target in sorted(nav_targets):
        if not (_SITE_DOCS / target).is_file():
            missing.append(target)
    assert not missing, (
        "mkdocs.yml nav entries missing under site_docs/:\n"
        + "\n".join(f"  - {target}" for target in missing)
        + "\n\nEvery nav entry must exist as a real file — "
        "`mkdocs build --strict` would fail otherwise. Add "
        "placeholder stubs per the Phase 0+1 discipline."
    )


def test_pin_present_makefile_carries_docs_build_targets() -> None:
    """Phase 1 — the Makefile carries `docs-install`,
    `docs-build`, `docs-serve`, `docs-clean` targets so
    the contributor workflow is uniform with the other
    `make` targets (test, validate-fast, etc.)."""
    makefile = _REPO_ROOT / "Makefile"
    if not makefile.is_file():
        pytest.skip("no Makefile at repo root — docs targets optional")
    text = makefile.read_text(encoding="utf-8")
    for target in ("docs-install:", "docs-build:", "docs-serve:"):
        assert target in text, (
            f"Makefile missing target {target!r} — Phase 1 "
            f"contributor workflow regressed"
        )
    # The build target must use --strict.
    assert "mkdocs build --strict" in text, (
        "Makefile docs-build target doesn't use --strict — broken "
        "links / nav drift would silently pass CI"
    )


# ──────────────────────────────────────────────────────────────────
# I. pin_present_* — Phase 1.5 visual parity (GenAIRR theme layer)
# ──────────────────────────────────────────────────────────────────

_EXTRA_CSS = _SITE_DOCS / "stylesheets" / "extra.css"

# Canonical brand tokens that MUST appear in the theme layer. These
# are the visual signature — IBM Plex triad + three-color semantic
# system (green/cobalt/plum) — sourced from website/styles.css.
_REQUIRED_BRAND_TOKENS = {
    "fonts": ('"IBM Plex Sans"', '"IBM Plex Mono"', '"IBM Plex Serif"'),
    "colors": ("#fafaf7", "#0e0e10", "#1f8a4c", "#1d4ed8", "#6e3a8e"),
    "css_class_patterns": (
        ".eyebrow",
        ".feat-grid",
        ".cap-grid",
        ".cta-row",
        ".install-band",
        ".engine-strip",
        ".btn-primary",
        ".btn-ghost",
    ),
    "material_overrides": (
        "--md-primary-fg-color",
        "--md-accent-fg-color",
        "--md-typeset-color",
    ),
}


def test_pin_present_extra_css_exists_under_site_docs_stylesheets() -> None:
    """Phase 1.5 — `site_docs/stylesheets/extra.css`
    is the theme layer that prevents the site from
    looking like default MkDocs Material. Pinned at file
    level."""
    assert _EXTRA_CSS.is_file(), (
        f"{_EXTRA_CSS} missing — Phase 1.5 visual parity layer "
        f"regressed"
    )


def test_pin_present_mkdocs_yml_loads_extra_css() -> None:
    """Phase 1.5 — `mkdocs.yml`'s `extra_css:` block loads
    `stylesheets/extra.css`. Otherwise the theme tokens
    wouldn't reach the rendered pages."""
    yml_text = _MKDOCS_YML.read_text(encoding="utf-8")
    # `extra_css:` is a top-level key; look for the active block.
    assert "extra_css:" in yml_text, (
        "mkdocs.yml has no `extra_css:` block — the theme layer "
        "won't be loaded. Phase 1.5 regressed."
    )
    assert "stylesheets/extra.css" in yml_text, (
        "mkdocs.yml `extra_css:` block doesn't reference "
        "stylesheets/extra.css — Phase 1.5 wiring regressed."
    )


def test_pin_present_extra_css_carries_ibm_plex_font_stack() -> None:
    """Phase 1.5 — `extra.css` references the IBM Plex
    triad (Sans / Mono / Serif). Mono is the canonical
    visual signature of the design — it carries
    eyebrows, labels, code, and chrome typography."""
    text = _EXTRA_CSS.read_text(encoding="utf-8")
    for font in _REQUIRED_BRAND_TOKENS["fonts"]:
        assert font in text, (
            f"extra.css missing IBM Plex font reference {font!r} — "
            f"typography signature regressed"
        )


def test_pin_present_extra_css_carries_bench_science_color_palette() -> None:
    """Phase 1.5 — `extra.css` carries the bench-science
    palette: paper-white background, near-black ink, and
    the three-color semantic system (green accent, cobalt
    info, plum archive). Each is the canonical token from
    the website/styles.css source."""
    text = _EXTRA_CSS.read_text(encoding="utf-8")
    for color in _REQUIRED_BRAND_TOKENS["colors"]:
        assert color in text, (
            f"extra.css missing brand color {color!r} — the "
            f"bench-science palette signature regressed. Tokens: "
            f"#fafaf7=paper-white bg, #0e0e10=near-black ink, "
            f"#1f8a4c=GREEN accent, #1d4ed8=COBALT info, "
            f"#6e3a8e=PLUM archive."
        )


def test_pin_present_extra_css_carries_genairr_component_classes() -> None:
    """Phase 1.5 — `extra.css` defines the recurring
    component classes the homepage and hub pages use
    (`.eyebrow`, `.feat-grid`, `.cap-grid`, etc.). These
    are markdown-author-facing utilities that let content
    pages mirror the existing website's visual patterns."""
    text = _EXTRA_CSS.read_text(encoding="utf-8")
    for cls in _REQUIRED_BRAND_TOKENS["css_class_patterns"]:
        assert cls in text, (
            f"extra.css missing component class {cls!r} — content "
            f"pages that use this class would lose their styling"
        )


def test_pin_present_extra_css_overrides_material_default_palette() -> None:
    """Phase 1.5 — `extra.css` overrides Material's own
    palette CSS variables so the chrome (header, sidebar,
    accents) paints with GenAIRR's tokens. Without this,
    Material's indigo defaults would bleed through."""
    text = _EXTRA_CSS.read_text(encoding="utf-8")
    for override in _REQUIRED_BRAND_TOKENS["material_overrides"]:
        assert override in text, (
            f"extra.css missing Material variable override "
            f"{override!r} — Material's default palette would "
            f"bleed through the chrome (NOT acceptable per the "
            f"\"not generic MkDocs\" Phase 1.5 contract)"
        )


def test_pin_present_homepage_no_longer_a_placeholder() -> None:
    """Phase 1.5 — `site_docs/index.md` is rewritten as a
    polished landing page; no "Migration placeholder"
    text remains. The homepage is the first impression;
    it must carry GenAIRR's identity, not a stub."""
    text = (_SITE_DOCS / "index.md").read_text(encoding="utf-8")
    assert "Migration placeholder" not in text, (
        "site_docs/index.md still contains 'Migration placeholder' "
        "text — the Phase 1.5 homepage rewrite regressed"
    )
    # And the homepage must use the design-token utility classes
    # so the build actually applies the theme.
    for required_class in ("eyebrow", "cta-row", "feat-grid",
                           "cap-grid", "install-band"):
        assert required_class in text, (
            f"homepage missing the {required_class!r} class — the "
            f"polished landing-page IA regressed"
        )


def test_pin_present_homepage_carries_core_ctas() -> None:
    """Homepage exposes the four primary CTAs. Phase 2A
    repointed the Validation + Reference-cartridges CTAs
    at the real ported content (was: placeholder hubs)."""
    text = (_SITE_DOCS / "index.md").read_text(encoding="utf-8")
    cta_keywords = (
        ("First simulation", "getting-started/quick-start.md"),
        ("Reference cartridges", "concepts/reference-cartridge.md"),
        ("Validate AIRR records", "validation/validate-records.md"),
        ("API reference", "reference/index.md"),
    )
    for keyword, target_md in cta_keywords:
        assert keyword in text, (
            f"homepage missing CTA text {keyword!r} — user-brief "
            f"core CTA list regressed"
        )
        assert target_md in text, (
            f"homepage missing CTA target {target_md!r} — the link "
            f"that backs the {keyword!r} CTA is gone"
        )
        # Target must exist as a real file.
        assert (_SITE_DOCS / target_md).is_file(), (
            f"homepage CTA target {target_md!r} doesn't exist on "
            f"disk — `mkdocs build --strict` would fail"
        )


def test_pin_present_existing_website_directory_remains_untouched() -> None:
    """Phase 1.5 — the live `website/` site is NOT
    modified in this slice. The hand-rolled HTML
    continues to deploy until the Phase 6 cutover. Pinned
    so an accidental edit during visual-parity work
    surfaces here."""
    assert _WEBSITE.is_dir()
    # Spot-check that the load-bearing pages still exist with
    # non-trivial content (i.e. we didn't accidentally empty them).
    for page in ("index.html", "styles.css", "learn.html",
                 "guides.html", "concepts.html", "reference.html"):
        target = _WEBSITE / page
        assert target.is_file(), (
            f"{target} missing — website/ content modified during "
            f"Phase 1.5; deploy workflow would now break"
        )
        assert target.stat().st_size > 1000, (
            f"{target} suspiciously small (size {target.stat().st_size}) "
            f"— website/ may have been emptied"
        )


def test_pin_post_cutover_deploy_workflow_runs_mkdocs_build() -> None:
    """Post-Phase-6-cutover: the deploy workflow installs
    docs dependencies and runs `make docs-build` (which
    runs `mkdocs build --strict`) before uploading the
    generated `site/` to GitHub Pages."""
    deploy_text = _DEPLOY_WORKFLOW.read_text(encoding="utf-8")
    assert "make docs-build" in deploy_text, (
        "deploy-docs.yml no longer runs `make docs-build` — Phase 6 "
        "build step regressed; the upload would be stale or missing."
    )
    assert "docs/requirements-docs.txt" in deploy_text, (
        "deploy-docs.yml no longer installs docs dependencies from "
        "docs/requirements-docs.txt — the MkDocs build would fail "
        "on the runner."
    )
    assert "path: site" in deploy_text, (
        "deploy-docs.yml no longer uploads `site/` — the MkDocs "
        "build output is not the publish source."
    )


def test_pin_present_mkdocs_build_strict_succeeds() -> None:
    """Phase 1 — if the docs deps are installed in the
    test environment, `mkdocs build --strict` succeeds on
    the current scaffold. Otherwise the test skips
    cleanly with a helpful message.

    This is the load-bearing build-gate pin: every
    subsequent migration phase must keep the strict build
    passing. Once Phase 6 cutover lands, the
    `.github/workflows/deploy-docs.yml` workflow uses the
    same command."""
    try:
        import mkdocs  # noqa: F401
    except ImportError:
        pytest.skip(
            "docs deps not installed — run "
            "`pip install -r docs/requirements-docs.txt` (or "
            "`make docs-install`) to enable this pin. The deps are "
            "documented as optional per migration plan §6 — Phase 1 "
            "contributor workflow only."
        )
    # Run the build in a clean subprocess so it picks up the same
    # environment the user / CI would use.
    env = dict(os.environ)
    # Silence the properdocs nag in the subprocess output.
    env["DISABLE_MKDOCS_2_WARNING"] = "true"
    result = subprocess.run(
        [sys.executable, "-m", "mkdocs", "build", "--strict"],
        cwd=str(_REPO_ROOT),
        capture_output=True,
        text=True,
        env=env,
        timeout=120,
    )
    assert result.returncode == 0, (
        f"`mkdocs build --strict` failed (exit {result.returncode}).\n"
        f"--- stdout ---\n{result.stdout[-2000:]}\n"
        f"--- stderr ---\n{result.stderr[-2000:]}\n"
        "Either the build broke (fix the scaffold) or this pin "
        "needs adjusting (e.g. if the build target changed)."
    )
    # Clean up the build artifact so the test doesn't leave
    # `site/` behind on each run.
    site_dir = _REPO_ROOT / "site"
    if site_dir.is_dir():
        shutil.rmtree(site_dir)
