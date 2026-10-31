"""Inventory gate: prove a change to tests/ did not quietly weaken a test.

Records, for every test function, helper and module-level constant, a hash of
its normalized source, then compares a branch against its merge base. A moved
test is invisible; an edited one is a violation until the author names it.

Normalization drops the signature (so moving a function into a class, which
adds ``self``, is invisible) and every import statement (so rewriting an
import is invisible), then unparses the AST -- which discards comments,
blank lines and indentation.  Decorators are kept, because editing a
``parametrize`` table changes what is tested.

Nothing here knows which changes are legitimate. The author declares them,
in ``tools/inventory-declarations.json``, with a reason per entry:

    {"revised": {"test_x": "why its body changed"},
     "replaced": {"test_old": {"by": "test_new", "why": "..."}},
     "allow":    {"dead-code": "why a helper went away"}}

A declaration that no longer explains anything is reported stale, so the file
empties itself out as branches merge rather than growing into a standing
allowlist.

Usage:
    python tools/inventory_gate.py check --since origin/master
    python tools/inventory_gate.py check --baseline PATH [--revised NAME ...]
    python tools/inventory_gate.py capture --out PATH [--tests-dir DIR]
    python tools/inventory_gate.py duplicates
"""
from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ALLOWED_TRANSFORMS = frozenset({
    "coverage-registration",  # adding a @verifies decorator to an existing test
    "cross-import",        # rewriting an import to point at tests/support
    "path-resolution",     # replacing inline repo-root arithmetic with support.paths
    "audit-doc-filename",  # updating a quoted module filename after a rename
    "module-marker",       # adding a module-level pytestmark
    "module-split",        # a constant naming a split module, replaced per target
    "dead-code",           # a member removed because nothing referenced it
    "pin-resolved",        # a strict xfail removed because its defect is fixed
    "renamed",             # a test renamed, body unchanged -- checked, not trusted
})

OPTIONAL_PACKAGES = ("pandas", "numpy", "scipy", "fastmcp", "hypothesis", "xdist", "graphviz", "tqdm")

#: Sections a baseline must carry for the gate to run at full strength. Each
#: one gates a check that is skipped outright when the section is absent, so a
#: short baseline degrades the gate without saying so.
REQUIRED_BASELINE_SECTIONS = frozenset({
    "functions", "sources", "members", "member_sources",
    "constants", "constant_sources",
    # `collected` guards the node-id comparison -- the check whose own comment
    # says "a gain plus a loss that nets to zero was invisible, as was a
    # parametrize id changing". It was readable with `.get("collected", ())`
    # and skipped on a falsy result, so the first version of this requirement
    # left out the very check it was written for.
    "collected",
})


# ── hashing ──────────────────────────────────────────────────────────────

class _DropEmptyFPrefix(ast.NodeTransformer):
    """Rewrite ``f"no placeholders"`` to ``"no placeholders"``.

    ruff's F541 strips that prefix, and the result is the same string, so the
    difference is formatting of exactly the kind this normalizer exists to
    discard. Handled here rather than by an allow class: a class meaning "a
    linter touched it" would explain anything a linter can reach.

    Only a ``JoinedStr`` whose every part is a literal is rewritten. ``f"{x}"``
    keeps its prefix, so deleting a substitution -- which does change what is
    asserted -- still alters the hash.
    """

    def visit_JoinedStr(self, node: ast.JoinedStr):
        self.generic_visit(node)
        if all(isinstance(v, ast.Constant) and isinstance(v.value, str)
               for v in node.values):
            joined = "".join(v.value for v in node.values)
            return ast.copy_location(ast.Constant(value=joined), node)
        return node


def _unparse(node: ast.AST) -> str:
    """``ast.unparse`` with the formatting-only differences normalized away."""
    return ast.unparse(_DropEmptyFPrefix().visit(copy.deepcopy(node)))


def _record_lines(text: str) -> list[str]:
    """``text`` split into the lines the ``--allow`` token rules are matched on.

    ``ast.unparse`` renders a multi-line string -- a docstring, most often --
    as a single logical line with ``\n`` escapes in it. Splitting only on real
    newlines therefore made a whole docstring one entry, and the token rules
    require merely that *a* changed line mention a token the class is about: so
    one occurrence of ``tests/`` or ``.py`` anywhere in a docstring licensed
    rewriting every other sentence in it.

    Splitting the escapes too makes the rule mean what it says -- a token
    licenses the line it is on, and nothing else.
    """
    return [part for line in text.splitlines() for part in line.split("\\n")]


def _normalized_body(node: ast.AST) -> str:
    """The function's decorators and body, minus imports, as source text."""
    kept = [s for s in node.body if not isinstance(s, (ast.Import, ast.ImportFrom))]
    parts = [_unparse(d) for d in node.decorator_list]
    parts.extend(_unparse(s) for s in kept)
    return "\n".join(parts)


def normalized_sources(path: Path) -> dict[str, list[str]]:
    """Test function name -> its normalized body lines.

    Recorded so ``--allow`` can check a change against the class claimed for
    it rather than taking the class on trust.
    """
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    out: dict[str, list[str]] = defaultdict(list)
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_"):
            out[node.name].extend(_record_lines(_normalized_body(node)))
    return dict(out)


def function_hashes(path: Path) -> dict[str, list[str]]:
    """Map each test function name in ``path`` to the hashes of its occurrences.

    A name can occur more than once in one module (two classes, say), so the
    value is a list and comparisons are multiset comparisons.
    """
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    out: dict[str, list[str]] = defaultdict(list)
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_"):
            digest = hashlib.sha256(_normalized_body(node).encode("utf-8")).hexdigest()
            out[node.name].append(digest)
    return {name: sorted(hashes) for name, hashes in out.items()}


def member_sources(path: Path) -> dict[str, list[str]]:
    """Qualified name -> normalized body lines, for every def/class.

    Recorded for the same reason as test-function sources: a member change has
    to be checked against the transform class claimed for it, not suppressed
    because some class was named.
    """
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    out: dict[str, list[str]] = {}

    def walk(node, prefix: str) -> None:
        for child in getattr(node, "body", []):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                qualname = f"{prefix}{child.name}"
                out[qualname] = _record_lines(_normalized_body(child))
                if isinstance(child, ast.ClassDef):
                    walk(child, f"{qualname}.")

    walk(tree, "")
    return out


def member_hashes(path: Path) -> dict[str, list[str]]:
    """Qualified-name -> hashes of every def/class body in the module.

    Function-body hashing covered only names starting with ``test_``, which
    left the shared layer unaudited -- and that is where the risk is: a gutted
    assertion helper makes tests pass MORE easily, so "suite green + gate
    clean" was satisfiable by a silently weakened suite.

    Keyed by qualified name because a bare name is not an identity here: 22
    modules define a function of the same name, so a bare-name multiset could
    not tell two modules' bodies apart if they were swapped.
    """
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    out: dict[str, list[str]] = defaultdict(list)

    def walk(node, prefix: str) -> None:
        for child in getattr(node, "body", []):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                qualname = f"{prefix}{child.name}"
                digest = hashlib.sha256(_normalized_body(child).encode("utf-8")).hexdigest()
                out[qualname].append(digest)
                if isinstance(child, ast.ClassDef):
                    walk(child, f"{qualname}.")

    walk(tree, "")
    return {name: sorted(h) for name, h in out.items()}


def _module_functions(paths, tests_dir: str = "tests") -> dict[str, dict[str, list[str]]]:
    return {_as_suite_path(p, tests_dir): function_hashes(p) for p in paths}


# ── module-level constants ───────────────────────────────────────────────

def _assigned_names(node):
    """(name, value_node) for each top-level name a statement binds.

    Tuple unpacking binds several names at once; where the right-hand side is
    a matching tuple each name is paired with its own element, so changing one
    locus id does not blur across all three.
    """
    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
    value = node.value
    for target in targets:
        if isinstance(target, ast.Name):
            yield target.id, value
        elif isinstance(target, (ast.Tuple, ast.List)):
            elts = target.elts
            values = value.elts if isinstance(value, (ast.Tuple, ast.List)) \
                and len(value.elts) == len(elts) else [value] * len(elts)
            for elt, val in zip(elts, values):
                if isinstance(elt, ast.Name):
                    yield elt.id, val


def constant_sources(path: Path) -> dict[str, list[str]]:
    """Module-level bound name -> the source of its assigned value."""
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    out: dict[str, list[str]] = {}
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)) and node.value is not None:
            for name, value in _assigned_names(node):
                out[name] = _record_lines(_unparse(value))
    return out


def constant_hashes(path: Path) -> dict[str, list[str]]:
    """Map each module-level bound name to the hashes of its assigned values.

    Function-body hashing cannot see a changed constant, and this suite keeps
    load-bearing values -- locus identifiers, format tuples, data roots -- at
    module level.  Without this, a migration that retargeted one would pass a
    body-only gate clean.
    """
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    out: dict[str, list[str]] = defaultdict(list)
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)) and node.value is not None:
            for name, value in _assigned_names(node):
                digest = hashlib.sha256(ast.unparse(value).encode("utf-8")).hexdigest()
                out[name].append(digest)
    return {name: sorted(h) for name, h in out.items()}


# ── duplicates ───────────────────────────────────────────────────────────

def duplicate_groups(paths) -> list[list[tuple[str, str]]]:
    """Groups of (module, function) that share an identical normalized body.

    Candidates only.  Per design section 8.2 a candidate is a genuine
    duplicate only if its resolved inputs are identical too, which a human
    decides.  Occurrences inside a single module are never reported: a module
    that defines the same body twice is a different question from two modules
    that duplicate each other.
    """
    by_hash: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for path in paths:
        module = Path(path).name
        for name, hashes in function_hashes(path).items():
            for digest in hashes:
                by_hash[digest].append((module, name))
    groups = []
    for occurrences in by_hash.values():
        modules = {m for m, _ in occurrences}
        if len(modules) > 1:
            groups.append(sorted(occurrences))
    return sorted(groups)


# ── capture ──────────────────────────────────────────────────────────────

def _detected_extras() -> list[str]:
    import importlib.util
    return [name for name in OPTIONAL_PACKAGES if importlib.util.find_spec(name) is not None]


def _as_suite_path(path, tests_dir: str) -> str:
    """Key a module by where it sits in the suite, not by where it was read.

    A capture run with ``--tests-dir`` pointing at a git worktree wrote
    absolute keys. The member and constant checks fall back to comparing bare
    qualnames when a key's path differs, so such a baseline silently degraded
    both to name-only for the whole suite -- a helper moving between modules
    would have passed. Keyed relative to the suite root, a capture is the same
    wherever the suite is read from.
    """
    root = Path(tests_dir)
    return "/".join((root.name, *Path(path).relative_to(root).parts))


def _collected_node_ids(tests_dir: str) -> list[str]:
    proc = subprocess.run([sys.executable, "-m", "pytest", tests_dir, "--collect-only", "-q"],
                          capture_output=True, text=True)
    ids = [line.strip() for line in proc.stdout.splitlines() if "::" in line and not line.startswith(" ")]
    if not ids:
        raise SystemExit(f"collection produced no node ids; pytest said:\n{proc.stdout[-2000:]}")
    # pytest reports node ids under the path it was given, so a capture from
    # elsewhere would record ids no run from the repo root can produce.
    prefix = str(Path(tests_dir).parent).replace("\\", "/").rstrip("/") + "/"
    return [i[len(prefix):] if prefix != "./" and i.startswith(prefix) else i for i in ids]


def _scanned(tests_dir: str):
    """Every module under the suite, support packages included.

    Scanning only ``test_*.py`` would read a helper's move into
    ``tests/support`` as a deletion -- the gate would fire on exactly the
    refactor it exists to permit.
    """
    return sorted(p for p in Path(tests_dir).rglob("*.py") if "__pycache__" not in p.parts)


def capture(tests_dir: str = "tests") -> dict:
    paths = _scanned(tests_dir)
    modules = _module_functions(paths, tests_dir)
    functions: dict[str, list[str]] = defaultdict(list)
    locations: dict[str, list[str]] = defaultdict(list)
    for module, names in modules.items():
        for name, hashes in names.items():
            functions[name].extend(hashes)
            locations[name].extend([module] * len(hashes))
    # Keyed by path::name for the same reason members are: pooling by bare name
    # across 171 modules makes two modules swapping a constant invisible, and
    # makes a legitimate addition whose name exists elsewhere a false positive.
    constants: dict[str, list[str]] = {}
    constant_src: dict[str, list[str]] = {}
    for path in paths:
        relative = _as_suite_path(path, tests_dir)
        for name, hashes in constant_hashes(path).items():
            constants[f"{relative}::{name}"] = hashes
        for name, lines in constant_sources(path).items():
            constant_src[f"{relative}::{name}"] = lines
    # Members are keyed by path::qualname so module identity is kept.
    members: dict[str, list[str]] = {}
    member_src: dict[str, list[str]] = {}
    for path in paths:
        relative = _as_suite_path(path, tests_dir)
        for name, hashes in member_hashes(path).items():
            members[f"{relative}::{name}"] = hashes
        for name, lines in member_sources(path).items():
            member_src[f"{relative}::{name}"] = lines
    sources: dict[str, list[str]] = defaultdict(list)
    for path in paths:
        for name, lines in normalized_sources(path).items():
            sources[name].extend(lines)
    collected = _collected_node_ids(tests_dir)
    return {
        "env": {
            "python": ".".join(str(v) for v in sys.version_info[:3]),
            "extras": _detected_extras(),
            "collected_count": len(collected),
            "module_count": len(paths),
        },
        "collected": collected,
        "functions": {name: sorted(h) for name, h in functions.items()},
        "constants": constants,
        "constant_sources": constant_src,
        "members": members,
        "member_sources": member_src,
        "sources": {name: lines for name, lines in sources.items()},
        "modules": {name: sorted(set(loc)) for name, loc in locations.items()},
    }


# ── check ────────────────────────────────────────────────────────────────

def _env_line(baseline: dict, current: dict) -> str:
    be, ce = baseline.get("env", {}), current.get("env", {})
    return (f"baseline env: python={be.get('python')} extras={sorted(be.get('extras', []))} "
            f"count={be.get('collected_count')}\n"
            f"current  env: python={ce.get('python')} extras={sorted(ce.get('extras', []))} "
            f"count={ce.get('collected_count')}")


#: Classes that can explain a changed *body*, and the tokens every differing
#: line must mention for that class to account for it.
_BODY_TOKENS = {
    # Adding a registration decorator to a test that already existed. Its own
    # class because the alternative was cross-import explaining it: the token
    # "from" matched the word "from" inside the claim text, which is the same
    # prose-matching problem as the "/" token, found the same way.
    "coverage-registration": ("verifies(",),
    # "import" and "from" as whole statements, not as English words. A claim
    # string reading "...tell one case from another" was being accepted as an
    # import rewrite.
    "cross-import": ("import ", "from .", "from tests", "from GenAIRR"),
    # "/" used to be in this list and explained far too much: every changed
    # line containing a slash counted as path resolution, including a docstring
    # with a markdown link in it. The tokens here all name the mechanism of the
    # transform -- the support module, the attribute arithmetic it replaces --
    # rather than merely looking like a path.
    "path-resolution": ("paths.", "__file__", "Path(", "parent", "parents"),
    # ".md" and "audit-docs" belong here, not above: repointing a reference to
    # a design document is a documentation filename change, and giving it its
    # own token keeps it from borrowing path-resolution's.
    "audit-doc-filename": (".py", "tests/", ".md", "audit-docs"),
    "module-marker": ("pytestmark",),
}

#: Classes that explain a *removal* and nothing else. These had no token list,
#: and `if not tokens: return []` made them explain every changed body as well
#: -- so passing --allow dead-code, which the two largest commits on this
#: branch needed, switched the classification off entirely. A class listed here
#: can account for a vanished member or a removed constant, never an edit.
#: `pin-resolved` is for a strict xfail whose defect has been fixed and which
#: is replaced by the positive test. It is not dead-code -- nothing was
#: unreferenced -- and it is not a `--revised` body change, because the
#: function is gone, which `--revised` deliberately refuses to excuse. A
#: strict xfail that starts passing is a loud signal by design, so the commit
#: that removes one is exactly where it should have to be named.
#: `renamed` joins them: it cannot explain an edit either, and unlike the
#: others it carries its own evidence -- the vanished name is excused only
#: where its exact baseline body hash is found under another name.
_REMOVAL_ONLY = frozenset({"dead-code", "module-split", "pin-resolved", "renamed"})


def _body_tokens(allow) -> tuple[str, ...]:
    """Tokens from the body-explaining classes in ``allow``.

    Removal-only classes contribute nothing, so naming one cannot explain an
    edit; an unknown class is rejected in ``check`` before this is reached.
    """
    return tuple(tok for cls in allow for tok in _BODY_TOKENS.get(cls, ()))


def _note_revised(used, names) -> None:
    """Record that ``revised`` earned its keep on these names.

    Tracked so a declaration left behind after its branch merged is reported
    stale instead of quietly excusing the next edit to the same function.
    """
    if used is not None:
        used.setdefault("revised", set()).update(names)


def _unexplained(changed: list[str], baseline: dict, current: dict, allow,
                 revised=frozenset(), used=None) -> list[str]:
    """Which changed bodies the named transform classes do not account for.

    ``--allow`` used to validate a spelling and then permit an unlimited number
    of arbitrary edits, reported as a note. A class now has to explain the
    change: every differing line must mention a token that class is about, so
    naming ``cross-import`` cannot wave through a loosened comparison.
    """
    _note_revised(used, [name.split("::", 1)[-1] for name in changed
                         if name.split("::", 1)[-1] in revised])
    changed = [name for name in changed
               if name.split("::", 1)[-1] not in revised]
    if not allow:
        return list(changed)
    base_src, cur_src = baseline.get("sources"), current.get("sources")
    if not base_src or not cur_src:
        return list(changed)          # nothing to judge against: do not permit
    tokens = _body_tokens(allow)
    if not tokens:
        return list(changed)          # no body-explaining class was named
    out = []
    for name in changed:
        before = set(base_src.get(name, ()))
        after = set(cur_src.get(name, ()))
        differing = (before ^ after)
        if any(not any(tok in line for tok in tokens) for line in differing):
            out.append(name)
    return out


def _unexplained_constants(changed: list[str], baseline: dict, current: dict, allow,
                           revised=frozenset(), used=None) -> list[str]:
    """Constant changes the named body-explaining classes do not account for.

    `revised` names constants the same way it names functions: an enumeration,
    matched on the name after `::`. A constant that encodes a decision -- the
    mechanism catalog, a pinned distribution -- changes deliberately, and a
    text pattern is the wrong instrument for saying so.
    """
    _note_revised(used, [k.split("::", 1)[-1] for k in changed
                         if k.split("::", 1)[-1] in revised])
    changed = [key for key in changed if key.split("::", 1)[-1] not in revised]
    tokens = _body_tokens(allow)
    if not tokens:
        return list(changed)
    base_src, cur_src = baseline.get("constant_sources"), current.get("constant_sources")
    if not base_src or not cur_src:
        return list(changed)          # nothing to judge against: do not permit
    out = []
    for key in changed:
        differing = set(base_src.get(key, ())) ^ set(cur_src.get(key, ()))
        if any(not any(tok in line for tok in tokens) for line in differing):
            out.append(key)
    return out


def _unexplained_members(changed: list[str], baseline: dict, current: dict, allow,
                         revised=frozenset(), used=None) -> list[str]:
    """Member changes the named transform classes do not account for."""
    _note_revised(used, [k.split("::", 1)[-1] for k in changed
                         if k.split("::", 1)[-1] in revised])
    changed = [key for key in changed if key.split("::", 1)[-1] not in revised]
    tokens = _body_tokens(allow)
    if not tokens:
        return list(changed)
    base_src, cur_src = baseline.get("member_sources"), current.get("member_sources")
    if not base_src or not cur_src:
        return list(changed)          # nothing to judge against: do not permit          # no body-explaining class was named

    def candidates(mapping, key):
        """Bodies in ``mapping`` that could be this member.

        Exactly one when the member stayed put. When it relocated, every body
        sharing its qualified name is a candidate -- and they must be compared
        one at a time, not unioned: 22 modules define a function of the same
        name, and unioning their bodies makes the diff pick up lines from
        modules that have nothing to do with the change.
        """
        if key in mapping:
            return [set(mapping[key])]
        qualname = key.split("::", 1)[-1]
        return [set(value) for other, value in mapping.items()
                if other.split("::", 1)[-1] == qualname] or [set()]

    out = []
    for key in changed:
        before_options = candidates(base_src, key)
        after_options = candidates(cur_src, key)
        explained = any(
            all(any(tok in line for tok in tokens) for line in (before ^ after))
            for before in before_options for after in after_options)
        if not explained:
            out.append(key)
    return out


def check(baseline: dict, current: dict, allow=frozenset(), revised=frozenset(),
          replaced=None, used=None) -> list[str]:
    """Violations of the migration invariant.  Empty list means clean.

    ``revised`` names individual functions whose body change is accepted, by
    exact name. A transform class explains a change by matching text, which is
    how every leak in this gate has worked -- "/" matched any path-like line,
    "from" matched the word in a sentence. An enumeration cannot over-explain:
    it accepts the functions listed and nothing else, and the list is the
    audit trail. Use it for a change that genuinely alters an assertion, where
    a text pattern would be both too weak and too broad.

    ``replaced`` maps a retired test's name to the name that carries its
    guarantee now. It excuses that one removal, and only once the replacement
    is found among the current test functions and among the collected node
    ids: a class alone would be a rubber stamp, where this costs the author a
    name that has to exist and run.
    """
    replaced = dict(replaced or {})
    # ``used`` collects which declarations did any work, so main() can report
    # the ones that did none. Precise for `revised` and `replaced`, coarse for
    # `allow`: a class explains a change through shared token matching, and
    # which of several named classes supplied the matching token is not
    # recovered -- so the allow classes are reported used if any change was
    # attributed to the set.
    unknown = set(allow) - ALLOWED_TRANSFORMS
    if unknown:
        raise ValueError(f"unknown transform class(es): {sorted(unknown)}; "
                         f"known: {sorted(ALLOWED_TRANSFORMS)}")
    # Every allowed class must have a rule. An empty rule used to mean "explains
    # everything", which is the opposite of what a transform class is for.
    ruleless = sorted(set(allow) - set(_BODY_TOKENS) - _REMOVAL_ONLY)
    if ruleless:
        raise ValueError(f"transform class(es) with no rule: {ruleless}")

    # A baseline missing a section turns that section's check off. Every such
    # check used `if baseline.get(section):`, so an older or partially captured
    # baseline ran a fraction of the gate and printed `clean` -- the gate's
    # pre-remediation strength wearing the post-remediation name. This raises
    # rather than returning a violation: a violation reads as a finding about
    # the code, and a short baseline is a finding about the baseline.
    missing = sorted(section for section in REQUIRED_BASELINE_SECTIONS
                     if section not in baseline)
    if missing:
        raise ValueError(
            f"baseline is missing {missing}, so those checks would be skipped "
            f"silently. Required sections: {sorted(REQUIRED_BASELINE_SECTIONS)}. "
            "Re-capture it with `inventory_gate.py capture` from the branch "
            "point; a baseline captured by an older version of this script is "
            "not comparable.")

    violations: list[str] = []
    base_fns, cur_fns = baseline["functions"], current["functions"]

    if baseline.get("env", {}).get("collected_count") is None:
        raise ValueError(
            "baseline has no env.collected_count, so the collected-count check "
            "would not run. Re-capture it with `inventory_gate.py capture`.")
    base_count = baseline["env"]["collected_count"]
    cur_count = current.get("env", {}).get("collected_count")
    if cur_count is not None and cur_count < base_count:
        violations.append(
            f"collected count dropped {base_count} -> {cur_count}.\n{_env_line(baseline, current)}\n"
            "If an optional extra is merely absent this is an environment difference, "
            "not lost coverage -- compare against a baseline captured in the same environment.")

    base_ids, cur_ids = set(baseline["collected"]), set(current.get("collected", ()))
    #: bare test names the current tree actually collects, for --replaced.
    collected_names = {i.split("::")[-1].split("[", 1)[0] for i in cur_ids}
    if True:
        # Compare the sets, not the counts: a gain plus a loss that nets to
        # zero was invisible, as was a parametrize id changing.
        #
        # Unconditional: `if base_ids:` skipped the whole comparison for a
        # baseline whose section was present but empty, which is the same
        # silent degradation the section requirement exists to prevent. An
        # empty section now compares against nothing and finds nothing, which
        # is the honest result rather than a skip.
        gone_ids = sorted(base_ids - cur_ids)
        bare = {i.split("::", 1)[-1] for i in cur_ids}
        truly_gone = [i for i in gone_ids if i.split("::", 1)[-1] not in bare]
        # A function named in `revised` may also have gained or lost a
        # parametrize dimension, which changes its ids. That is part of the
        # revision being named, not a separate silent loss -- and the check
        # still protects every function not named.
        if revised or replaced:
            excused = set(revised) | set(replaced)
            truly_gone = [
                i for i in truly_gone
                if i.split("::", 1)[-1].split("[", 1)[0] not in excused]
        if truly_gone:
            violations.append(
                f"collected node ids no longer present ({len(truly_gone)}): "
                f"{truly_gone[:10]}{' ...' if len(truly_gone) > 10 else ''}")

    base_members = baseline["members"]
    if True:          # present-but-empty compares to nothing; absent raises above
        cur_members = current.get("members", {})
        gone = sorted(set(base_members) - set(cur_members))
        moved = {k.split("::", 1)[-1] for k in cur_members}
        vanished = []
        for key in gone:
            bare = key.split("::", 1)[-1]
            if bare in moved:
                continue
            if bare in replaced:
                if used is not None:
                    used.setdefault("replaced", set()).add(bare)
                continue
            vanished.append(key)
        if vanished and ({"dead-code", "pin-resolved"} & set(allow)) and used is not None:
            used.setdefault("allow", set()).update(
                {"dead-code", "pin-resolved"} & set(allow))
        if vanished and not ({"dead-code", "pin-resolved"} & set(allow)):
            violations.append(
                f"support/helper members no longer defined anywhere ({len(vanished)}): "
                f"{vanished[:10]}. If a member was removed because nothing "
                "referenced it, pass --allow dead-code and name it in the commit.")
        elif vanished:
            print(f"note: {len(vanished)} member(s) removed, attributed to dead-code: {vanished}")
        # Compare full keys where the member stayed put -- that is what catches
        # two modules' bodies being swapped, which a qualname-only comparison
        # cannot see. Only a member whose key is gone is matched by qualname,
        # which is how a genuine relocation is tolerated.
        shared_keys = set(base_members) & set(cur_members)
        changed_members = sorted(
            k for k in shared_keys
            if Counter(base_members[k]) != Counter(cur_members[k]))
        relocated = []
        for key in sorted(set(base_members) - set(cur_members)):
            qualname = key.split("::", 1)[-1]
            elsewhere = [v for k, v in cur_members.items()
                         if k.split("::", 1)[-1] == qualname]
            if not elsewhere:
                continue                     # already reported as vanished
            if not any(Counter(v) == Counter(base_members[key]) for v in elsewhere):
                relocated.append(key)
        changed_members.extend(relocated)

        if changed_members:
            unexplained = _unexplained_members(changed_members, baseline, current, allow, revised, used)
            if unexplained:
                violations.append(
                    f"helper bodies changed that {sorted(allow) or 'no transform class'} "
                    f"does not explain ({len(unexplained)}): {unexplained[:10]}. "
                    "A changed assertion helper can make tests pass MORE easily, so "
                    "the suite being green does not cover this.")
            explained = [m for m in changed_members if m not in set(unexplained)]
            if explained and used is not None and allow:
                used.setdefault("allow", set()).update(allow)
            if explained:
                print(f"note: {len(explained)} helper change(s) attributed to "
                      f"{sorted(allow)}: {explained[:10]}")

    missing = sorted(set(base_fns) - set(cur_fns))
    # A test function going missing is the gate's loudest finding. Only
    # `pin-resolved` excuses it, and deliberately NOT `dead-code`: dead-code
    # means an unreferenced helper, it is the flag the two largest commits on
    # the restructure branch passed, and letting it reach test removals would
    # be the same over-broad class the `_REMOVAL_ONLY` rule was written to
    # stop. The node-id check still requires the name, so removing a test
    # costs both a class and a name.
    if missing and "pin-resolved" not in allow:
        # `renamed` is narrower than a class: the vanished name is excused only
        # where its baseline body hash turns up under a different name in the
        # current capture, which is evidence that the test moved rather than
        # went away. A rename that also edited the body is not covered, and
        # should not be: that is two changes wearing one name.
        if "renamed" in allow:
            elsewhere = {h for name, hashes in cur_fns.items() for h in hashes}
            renamed = [n for n in missing if set(base_fns[n]) <= elsewhere]
            if renamed:
                print(f"note: {len(renamed)} test(s) attributed to renamed, body "
                      f"unchanged: {renamed}")
            missing = [n for n in missing if n not in set(renamed)]
        accounted = []
        for name in list(missing):
            successor = replaced.get(name)
            if successor is None:
                continue
            if successor not in cur_fns:
                violations.append(
                    f"{name} is declared replaced by {successor}, which is not "
                    "a test function in the current tree")
                continue
            if successor not in collected_names:
                violations.append(
                    f"{name} is declared replaced by {successor}, which exists "
                    "but is not collected, so nothing runs in its place")
                continue
            accounted.append(name)
        if accounted:
            print(f"note: {len(accounted)} test(s) replaced by a named, collected "
                  f"successor: {[f'{n} -> {replaced[n]}' for n in accounted]}")
            if used is not None:
                used.setdefault("replaced", set()).update(accounted)
        missing = [n for n in missing if n not in set(accounted)]
    if missing and "pin-resolved" not in allow:
        violations.append(f"test functions no longer present ({len(missing)}): {missing}")

    shrunk = sorted(n for n in set(base_fns) & set(cur_fns) if len(cur_fns[n]) < len(base_fns[n]))
    if shrunk:
        violations.append(f"test functions with fewer occurrences than baseline: {shrunk}")

    changed = sorted(n for n in set(base_fns) & set(cur_fns)
                     if Counter(base_fns[n]) != Counter(cur_fns[n]))
    if changed:
        unexplained = _unexplained(changed, baseline, current, allow, revised, used)
        if unexplained:
            violations.append(
                f"normalized body changed for {len(unexplained)} function(s) that "
                f"{sorted(allow) or 'no transform class'} does not explain: "
                f"{unexplained[:20]}{' ...' if len(unexplained) > 20 else ''}")
        explained = [n for n in changed if n not in set(unexplained)]
        if explained and used is not None and allow:
            used.setdefault("allow", set()).update(allow)
        if explained:
            print(f"note: {len(explained)} body change(s) attributed to {sorted(allow)}: "
                  f"{explained[:20]}{' ...' if len(explained) > 20 else ''}")

    # Module-level constants. A baseline that predates them is refused above,
    # so there is no "no section to compare" case left to tolerate here.
    base_consts = baseline["constants"]
    if True:          # present-but-empty compares to nothing; absent raises above
        cur_consts = current.get("constants", {})
        gone = sorted(set(base_consts) - set(cur_consts))
        if gone and "module-split" not in allow:
            violations.append(
                f"module-level constants no longer defined anywhere: {gone}. "
                "If a split replaced one with per-target constants, pass "
                "--allow module-split.")
        elif gone:
            print(f"note: {len(gone)} constant(s) removed, attributed to module-split: {gone}")
        retargeted = sorted(n for n in set(base_consts) & set(cur_consts)
                            if Counter(base_consts[n]) != Counter(cur_consts[n]))
        if retargeted:
            unexplained = _unexplained_constants(retargeted, baseline, current, allow, revised, used)
            if unexplained:
                violations.append(
                    f"module-level constant value(s) changed that "
                    f"{sorted(allow) or 'no transform class'} does not explain "
                    f"({len(unexplained)}): {unexplained[:10]}. A retargeted constant "
                    "is invisible to body hashing.")
            explained = [c for c in retargeted if c not in set(unexplained)]
            if explained and used is not None and allow:
                used.setdefault("allow", set()).update(allow)
            if explained:
                print(f"note: {len(explained)} constant change(s) attributed to "
                      f"{sorted(allow)}: {explained[:20]}"
                      f"{' ...' if len(explained) > 20 else ''}")
    return violations



# ── declarations ─────────────────────────────────────────────────────────

DECLARATIONS_DEFAULT = "tools/inventory-declarations.json"


def load_declarations(path) -> tuple[set, set, dict]:
    """``(allow, revised, replaced)`` from a declarations file.

    Every entry carries a reason, and an empty one is refused. The reason is
    not machine-checkable -- nothing here can tell a real explanation from
    "needed" -- but it is what a reviewer reads, and requiring it is why the
    file is a record rather than a list of names that accumulated.
    """
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    unknown_sections = sorted(set(data) - {"revised", "replaced", "allow"})
    if unknown_sections:
        raise ValueError(
            f"{path}: unknown section(s) {unknown_sections}; expected "
            "'revised', 'replaced' and 'allow'")

    revised = {}
    for name, why in (data.get("revised") or {}).items():
        if not isinstance(why, str) or not why.strip():
            raise ValueError(f"{path}: revised['{name}'] has no reason")
        revised[name] = why

    replaced = {}
    for name, entry in (data.get("replaced") or {}).items():
        if not isinstance(entry, dict) or not entry.get("by"):
            raise ValueError(f"{path}: replaced['{name}'] must name its successor as 'by'")
        if not str(entry.get("why", "")).strip():
            raise ValueError(f"{path}: replaced['{name}'] has no reason")
        replaced[name] = entry["by"]

    allow = {}
    for name, why in (data.get("allow") or {}).items():
        if not str(why).strip():
            raise ValueError(f"{path}: allow['{name}'] has no reason")
        allow[name] = why

    return set(allow), set(revised), replaced


def stale_declarations(path, used: dict) -> list[str]:
    """Declared entries that explained nothing on this run.

    A declaration outlives the branch that needed it: once the change is in
    the merge base there is nothing left to excuse, and the entry would sit
    there waving through the next edit to the same name. Reporting it is what
    keeps the file from becoming a standing allowlist.
    """
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    out = []
    for section in ("revised", "replaced", "allow"):
        names = set(data.get(section) or {})
        for name in sorted(names - set(used.get(section, ()))):
            out.append(f"{section}['{name}']")
    return out


def _capture_at(ref: str) -> dict:
    """Capture the suite as of ``ref``, via a throwaway worktree.

    A committed baseline was the alternative and is 8 MB of JSON per commit,
    which is why the gate reads its baseline out of git instead.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        tree = str(Path(tmp) / "base")
        add = subprocess.run(["git", "worktree", "add", "--detach", tree, ref],
                             capture_output=True, text=True)
        if add.returncode != 0:
            raise SystemExit(f"could not check out {ref}:\n{add.stderr}")
        try:
            return capture(str(Path(tree) / "tests"))
        finally:
            subprocess.run(["git", "worktree", "remove", "--force", tree],
                           capture_output=True, text=True)


# ── cli ──────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p_cap = sub.add_parser("capture")
    p_cap.add_argument("--out", required=True)
    p_cap.add_argument("--tests-dir", default="tests")

    p_chk = sub.add_parser("check")
    p_chk.add_argument("--baseline", help="a file written by `capture`")
    p_chk.add_argument("--since", metavar="REF",
                       help="capture the baseline from this git ref instead "
                            "(what CI uses: the branch's merge base)")
    p_chk.add_argument("--declarations", nargs="?", const=DECLARATIONS_DEFAULT,
                       help=f"read allow/revised/replaced from a file "
                            f"(default {DECLARATIONS_DEFAULT})")
    p_chk.add_argument("--tests-dir", default="tests")
    p_chk.add_argument("--revised", action="append", default=[],
                       help="a function whose body change is accepted, by exact "
                            "name. Enumerated rather than pattern-matched, so it "
                            "cannot explain anything it does not name.")
    p_chk.add_argument("--replaced", action="append", default=[], metavar="OLD=NEW",
                       help="a retired test and the test that carries its "
                            "guarantee now. NEW must exist and be collected.")
    p_chk.add_argument("--allow", action="append", default=[],
                       help=f"one of {sorted(ALLOWED_TRANSFORMS)}")

    p_dup = sub.add_parser("duplicates")
    p_dup.add_argument("--tests-dir", default="tests")

    args = parser.parse_args(argv)

    if args.command == "capture":
        data = capture(args.tests_dir)
        Path(args.out).write_text(json.dumps(data, indent=1, sort_keys=True), encoding="utf-8")
        print(f"captured {data['env']['collected_count']} tests "
              f"across {data['env']['module_count']} modules -> {args.out}")
        return 0

    if args.command == "check":
        if bool(args.baseline) == bool(args.since):
            parser.error("pass exactly one of --baseline and --since")
        allow, revised, pairs = set(args.allow), set(args.revised), {}
        for entry in args.replaced:
            if entry.count("=") != 1:
                parser.error(f"--replaced takes OLD=NEW, got {entry!r}")
            old_name, new_name = entry.split("=")
            pairs[old_name] = new_name
        if args.declarations:
            declared_allow, declared_revised, declared_pairs = \
                load_declarations(args.declarations)
            allow |= declared_allow
            revised |= declared_revised
            pairs.update(declared_pairs)
        baseline = (_capture_at(args.since) if args.since
                    else json.loads(Path(args.baseline).read_text(encoding="utf-8")))
        used: dict = {}
        violations = check(baseline, capture(args.tests_dir), allow=allow,
                           revised=revised, replaced=pairs, used=used)
        if args.declarations:
            stale = stale_declarations(args.declarations, used)
            if stale:
                violations.append(
                    f"{len(stale)} declaration(s) in {args.declarations} "
                    f"explained nothing on this run: {stale}. A declaration "
                    "whose change is already in the baseline excuses the next "
                    "edit to the same name -- delete them.")
        if not violations:
            print("clean")
            return 0
        for v in violations:
            print(f"VIOLATION: {v}")
        return 1

    groups = duplicate_groups(
        sorted(p for p in Path(args.tests_dir).rglob("test_*.py") if "__pycache__" not in p.parts))
    if not groups:
        print("no cross-module duplicate bodies")
        return 0
    print(f"{len(groups)} candidate group(s) -- identical bodies across modules.")
    print("A group is a genuine duplicate only if its resolved inputs match too.")
    for group in groups:
        print("  " + ", ".join(f"{m}::{n}" for m, n in group))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
