"""Every preserved URL must land on a page that says something.

`mkdocs.yml`'s `redirect_maps` keeps the pre-MkDocs URLs resolvable. A
redirect is only worth having if its destination answers the question the
old URL answered: a reader who follows a bookmark or a citation into a
placeholder is worse off than one who gets a 404, because the 404 at least
tells them to search.

This ran green while 22 of the 28 redirects pointed at pages whose entire
body was an admonition saying the page had not been written. It did not
catch that, because nothing checked it. These two tests do.

They deliberately do not depend on the contributor-only migration plan the
sibling `test_docs_mkdocs_migration_contract` module parses, which is
absent from every clone and skips that whole module in CI.
"""

from __future__ import annotations

import re

import pytest

from tests.support.paths import repo_root

_MKDOCS_YML = repo_root() / "mkdocs.yml"
_SITE_DOCS = repo_root() / "site_docs"

# Phrases a placeholder page uses to admit it is one. Matched
# case-insensitively against the whole page body.
_PLACEHOLDER_MARKERS = (
    "has not been written yet",
    "this page is a placeholder",
    "migration placeholder",
    "content lands in phase",
)

pytestmark = pytest.mark.pin


def _redirect_map() -> dict[str, str]:
    """`{fake source .md path: destination .md path}` from `redirect_maps`."""
    text = _MKDOCS_YML.read_text(encoding="utf-8")
    pairs = re.findall(r'^\s*"([^"]+\.md)":\s*"([^"]+\.md)"\s*$', text, re.M)
    assert pairs, "no redirect_maps entries found in mkdocs.yml"
    return dict(pairs)


def test_pin_present_every_redirect_target_is_a_real_file() -> None:
    """A redirect whose destination does not exist fails the build late.

    `mkdocs-redirects` emits the redirect file regardless, so the breakage
    shows up as a 404 for a reader rather than as an error for us.
    """
    missing = sorted(
        f"{source} -> {target}"
        for source, target in _redirect_map().items()
        if not (_SITE_DOCS / target).is_file()
    )
    assert not missing, (
        "redirect destinations that do not exist under site_docs/:\n  "
        + "\n  ".join(missing)
    )


def test_pin_present_no_redirect_lands_on_a_placeholder() -> None:
    """The preserved URL must answer the question the old URL answered.

    If a page has not been written, point its redirect at the nearest page
    that has been, and delete the placeholder. Keeping the placeholder as a
    destination spends the URL on an apology.
    """
    landing_on_placeholder = []
    for source, target in sorted(_redirect_map().items()):
        page = _SITE_DOCS / target
        if not page.is_file():
            continue  # the test above owns this failure
        body = page.read_text(encoding="utf-8").lower()
        for marker in _PLACEHOLDER_MARKERS:
            if marker in body:
                landing_on_placeholder.append(f"{source} -> {target} ({marker!r})")
                break
    assert not landing_on_placeholder, (
        "redirects pointing at placeholder pages:\n  "
        + "\n  ".join(landing_on_placeholder)
        + "\n\nRepoint each at the nearest page that has content, and delete "
        "the placeholder."
    )
