"""The named schemes a clonal run derives its sub-simulation seeds under.

The seeds in here are a wire value: a family-replay bundle records which
schedule drew it, so the first schedule's arithmetic has to stay exactly
what it was while the current one has to stop sharing streams between
runs.
"""

import pytest

from GenAIRR._seed_schedule import BY_LABEL, CURRENT, _Additive, _Keyed

pytestmark = pytest.mark.behavior


def _read_source(path):
    """Read a source file as UTF-8, whatever the platform would have picked.

    Both scans below walk `src/GenAIRR`, which contains at least one file
    Windows' default cp1252 cannot decode (`_s5f_loader.py`, byte 0x81 at
    offset 1529). An omitted `encoding` is therefore not a style point
    here: it makes the guards pass on Linux and raise `UnicodeDecodeError`
    on Windows, which is how both of them first reached CI red.
    """
    return path.read_text(encoding="utf-8")


def _relative_key(path, root):
    """Name a file under `root` the same way on every platform.

    `str(path.relative_to(root))` spells the result with the platform's
    separator, so one file is `_compiled/plain.py` here and
    `_compiled\\plain.py` on Windows -- and an allowlist keyed the first
    way matches nothing the second way, turning a sanctioned site into a
    reported violation.
    """
    return path.relative_to(root).as_posix()


# The first schedule's arithmetic, as `origin/master` writes it at
# `_compiled/repertoire.py:106`, `_compiled/lineage.py:164,178,298` and
# `_compiled/repertoire.py:134`. Quoted here rather than read back off
# the schedule object: a replay has to reproduce the seeds the capture
# used, so this is the contract, not an implementation detail.
def _master_clone(run_seed, clone_index):
    return int(run_seed) + clone_index * 1_000_000


def _master_attempt(clone_seed, attempt_index):
    return clone_seed + attempt_index * 7_919


def _master_read_or_node(clone_seed, index):
    return clone_seed + 1 + index


class TestTheScheduleThatBundlesRecord:
    def test_its_label_is_the_one_written_into_existing_bundles(self):
        assert _Additive.label == "clone_1000000_attempt_7919_node_plus_one_v1"

    def test_it_derives_exactly_what_it_always_derived(self):
        v1 = _Additive()
        for run_seed in (0, 1, 7, 2**40):
            for clone_index in (0, 1, 2, 37):
                clone_seed = v1.clone(run_seed, clone_index)
                assert clone_seed == _master_clone(run_seed, clone_index)
                for index in (0, 1, 49, 7918):
                    assert v1.attempt(clone_seed, index) == _master_attempt(clone_seed, index)
                    assert v1.read(clone_seed, index) == _master_read_or_node(clone_seed, index)
                    assert v1.node(clone_seed, index) == _master_read_or_node(clone_seed, index)

    def test_its_defects_are_still_present_because_a_replay_needs_them(self):
        """Not a schedule anyone should draw under -- a schedule a bundle may
        have been drawn under. The three collisions are asserted so that
        the contrast below is a measurement rather than a claim."""
        v1 = _Additive()
        assert v1.clone(0, 1) == v1.clone(1_000_000, 0)
        assert v1.read(4242, 7918) == v1.attempt(4242, 1)
        assert v1.node(4242, 7918) == v1.read(4242, 7918)


class TestTheScheduleRunsDrawUnderNow:
    def test_it_is_the_current_one_and_says_so(self):
        assert CURRENT is BY_LABEL[_Keyed.label]
        assert _Keyed.label.endswith("_v2")

    def test_a_different_run_seed_is_a_different_lattice(self):
        seen = {}
        for run_seed in (0, 1, 2, 1_000_000, 2_000_000, 2**63 - 1):
            for clone_index in range(48):
                word = CURRENT.clone(run_seed, clone_index)
                assert word not in seen, f"{(run_seed, clone_index)} reuses {seen[word]}"
                seen[word] = (run_seed, clone_index)

    def test_the_axes_under_one_clone_stay_apart(self):
        clone_seed = CURRENT.clone(0, 0)
        seen = {}
        for axis in (CURRENT.attempt, CURRENT.read, CURRENT.node):
            for index in list(range(128)) + [7918, 7919, 15837, 999_999]:
                word = axis(clone_seed, index)
                assert word not in seen, f"{(axis.__name__, index)} reuses {seen[word]}"
                seen[word] = (axis.__name__, index)

    def test_every_seed_it_derives_is_one_the_engine_can_be_given(self):
        """A seed is handed straight to the engine, which takes a u64."""
        words = [CURRENT.clone(2**64 - 1, 2**32 - 1)]
        words.append(CURRENT.read(words[0], 2**32 - 1))
        words.append(CURRENT.node(words[0], 0))
        words.append(CURRENT.attempt(words[0], 49))
        assert all(0 <= word <= 2**64 - 1 for word in words)

    def test_an_index_it_cannot_name_is_refused_not_folded(self):
        with pytest.raises((OverflowError, ValueError)):
            CURRENT.clone(0, 2**32)


def test_a_recorded_label_names_exactly_one_schedule():
    assert set(BY_LABEL) == {_Additive.label, _Keyed.label}
    assert all(label == type(schedule).label for label, schedule in BY_LABEL.items())


class TestNothingDerivesASeedByArithmetic:
    """The Python counterpart of the engine's `rng::discipline` guard.

    The schedule is only worth having if nothing bypasses it, and
    "someone adds an offset to a seed" is not a thing review catches
    twice. The package is walked for arithmetic on anything named like a
    seed, and the one place that is deliberate carries its reason.
    """

    #: Where arithmetic on a seed is the contract rather than a mistake,
    #: with how many sites to expect. A count, so that a second one in
    #: the same file still fails.
    SANCTIONED = {
        "_compiled/plain.py": (
            1,
            "`stream(n, seed)` yields the run seeded `seed + i` on purpose: a "
            "seed names a record on the plain path, which is what lets two "
            "batches stitch into one dataset (`n=N` covers `seed`..`seed+N-1`, "
            "so the next batch starts at `seed+N`). The engine's `run_batch` "
            "does the same thing and calls it seed-stitching. Changing it "
            "would be an API decision about that contract, not a bug fix.",
        ),
    }

    @staticmethod
    def _sites(path):
        import ast

        found = []
        for node in ast.walk(ast.parse(_read_source(path))):
            if not isinstance(node, ast.BinOp) or not isinstance(
                node.op, (ast.Add, ast.Sub, ast.Mult)
            ):
                continue
            named = [n.id for n in ast.walk(node) if isinstance(n, ast.Name)]
            named += [n.attr for n in ast.walk(node) if isinstance(n, ast.Attribute)]
            if any("seed" in name for name in named):
                found.append(node.lineno)
        return found

    def _scan(self):
        from tests.support import paths

        root = paths.repo_root() / "src" / "GenAIRR"
        assert (root / "_seed_schedule.py").exists(), f"premise: scanning the package at {root}"
        return {
            _relative_key(path, root): self._sites(path)
            for path in sorted(root.rglob("*.py"))
            # The schedule is where the arithmetic is allowed to live: it
            # has to keep deriving what the earlier scheme derived.
            if path.name != "_seed_schedule.py" and self._sites(path)
        }

    def test_the_scan_can_see_arithmetic_at_all(self):
        """Guard the premise: a scanner that finds nothing anywhere would
        pass the test below while checking nothing."""
        import ast
        import tempfile
        import pathlib

        with tempfile.TemporaryDirectory() as directory:
            planted = pathlib.Path(directory) / "planted.py"
            planted.write_text(
                "def f(seed, i):\n    return seed + i * 1_000_000\n",
                encoding="utf-8")
            assert self._sites(planted) == [2]
            innocent = pathlib.Path(directory) / "innocent.py"
            innocent.write_text(
                "def f(a, b):\n    return a + b * 1_000_000\n", encoding="utf-8")
            assert self._sites(innocent) == []
            assert isinstance(
                ast.parse(planted.read_text(encoding="utf-8")), ast.Module)

    def test_every_site_left_is_one_with_a_reason(self):
        found = self._scan()
        assert set(found) == set(self.SANCTIONED), (
            f"seed arithmetic outside the schedule: "
            f"{ {k: v for k, v in found.items() if k not in self.SANCTIONED} }"
        )
        for path, lines in found.items():
            expected, why = self.SANCTIONED[path]
            assert len(lines) == expected, f"{path} lines {lines}: {why}"


class TestEverySeedDerivedComesFromThePrimitive:
    """The other half of the discipline, for derivations that are calls.

    Arithmetic is one way to make a seed out of a seed; calling something
    is the other, and the scan above cannot see it. Two axes used to
    derive theirs from ``random.Random(seed).getrandbits(63)`` -- which
    mixes, so it was never the collision the first clonal schedule had,
    but which answered "where did this seed come from" with a second
    generator the engine's determinism version does not reach.

    So a name that looks like a seed may be assigned from a call only if
    the call is one of these. The list is short on purpose: a derivation
    worth adding is worth naming here.
    """

    #: Callables a seed may be assigned from. `sub_seed` and the schedule
    #: are the primitive; `_seed` reads a seed label off a replay bundle,
    #: which is parsing rather than deriving.
    DERIVERS = {"sub_seed", "_sub_seed", "substream_seed", "_seed"}

    #: Objects whose every method is a derivation (the named schedules).
    DERIVING_OBJECTS = {"_schedule", "schedule"}

    @classmethod
    def _sites(cls, path):
        import ast

        found = []
        for node in ast.walk(ast.parse(_read_source(path))):
            if isinstance(node, ast.Assign):
                targets = node.targets
            elif isinstance(node, ast.AnnAssign):
                targets = [node.target]
            else:
                continue
            named = [t.id for t in targets if isinstance(t, ast.Name)]
            if not any("seed" in name for name in named):
                continue
            if not isinstance(node.value, ast.Call):
                continue
            func = node.value.func
            if isinstance(func, ast.Name) and func.id in cls.DERIVERS:
                continue
            if isinstance(func, ast.Attribute):
                owner = func.value
                if func.attr in cls.DERIVERS:
                    continue
                if isinstance(owner, ast.Name) and owner.id in cls.DERIVING_OBJECTS:
                    continue
                if isinstance(owner, ast.Attribute) and owner.attr in cls.DERIVING_OBJECTS:
                    continue
            found.append((node.lineno, ast.unparse(node.value)))
        return found

    def test_the_scan_can_see_the_derivation_it_was_written_for(self):
        """Guard the premise with the code this guard replaced.

        The two superseded sites, as they stood, plus the two that
        replaced them -- so the scan is shown rejecting the one and
        accepting the other rather than merely finding nothing.
        """
        import pathlib
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            superseded = pathlib.Path(directory) / "superseded.py"
            superseded.write_text(
                "def f(seed, n):\n"
                "    base_rng = random.Random(seed)\n"
                "    for i in range(n):\n"
                "        sub_seed = base_rng.getrandbits(63)\n"
                "        yield sub_seed\n",
                encoding="utf-8",
            )
            assert [line for line, _ in self._sites(superseded)] == [4]

            current = pathlib.Path(directory) / "current.py"
            current.write_text(
                "def f(seed, n):\n"
                "    for i in range(n):\n"
                "        subject_seed = _sub_seed(seed, 'subject', i)\n"
                "        clone_seed = _schedule.clone(seed, i)\n"
                "        yield subject_seed, clone_seed\n",
                encoding="utf-8",
            )
            assert self._sites(current) == []

    def test_nothing_in_the_package_derives_a_seed_another_way(self):
        from tests.support import paths

        root = paths.repo_root() / "src" / "GenAIRR"
        assert (root / "_seed_schedule.py").exists(), f"premise: scanning the package at {root}"
        offenders = {
            _relative_key(path, root): sites
            for path in sorted(root.rglob("*.py"))
            if (sites := self._sites(path))
        }
        assert not offenders, f"a seed derived outside the primitive: {offenders}"


class TestTheseScansDeclareTheirEncoding:
    """Both guards above walk the package, so both read files they did not
    write. `Path.read_text()` with no `encoding` uses the platform default
    -- UTF-8 on Linux and macOS, cp1252 on Windows -- and `src/GenAIRR`
    contains a file cp1252 cannot decode. So the omission does not fail
    where it is written; it fails on one of three platforms, after review,
    which is exactly how both scans first went red in CI.

    A scan that reads source is the wrong place to learn that lesson
    twice, so this module checks its own text I/O.
    """

    @staticmethod
    def _undeclared(path):
        """Every line in `path` that reads or writes text without saying
        in what encoding. Binary modes are exempt: they decode nothing."""
        import ast

        found = []
        for node in ast.walk(ast.parse(_read_source(path))):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
            if name not in ("read_text", "write_text", "open"):
                continue
            if name == "open" and len(node.args) >= 2:
                mode = node.args[1]
                if isinstance(mode, ast.Constant) and "b" in str(mode.value):
                    continue
            if any(kw.arg == "encoding" for kw in node.keywords):
                continue
            found.append(node.lineno)
        return found

    def test_the_check_can_see_an_undeclared_read(self):
        """Guard the premise: plant both forms and watch them separate."""
        import pathlib
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            bad = pathlib.Path(directory) / "bad.py"
            bad.write_text("import pathlib\nx = pathlib.Path('a').read_text()\n",
                           encoding="utf-8")
            assert self._undeclared(bad) == [2]

            good = pathlib.Path(directory) / "good.py"
            good.write_text(
                "import pathlib\nx = pathlib.Path('a').read_text(encoding='utf-8')\n"
                "y = open('b', 'rb').read()\n",
                encoding="utf-8")
            assert self._undeclared(good) == []

    def test_this_module_reads_and_writes_nothing_by_platform_default(self):
        from tests.support import paths

        here = paths.repo_root() / "tests" / "components" / "clonal" / "test_seed_schedule.py"
        assert here.exists(), f"premise: checking this module at {here}"
        assert self._undeclared(here) == []

    def test_the_package_holds_a_file_the_windows_default_cannot_read(self):
        """The premise of all of this. If every file in the package were
        ASCII the omission would be harmless and this guard would be
        guarding nothing, so the file that makes it matter is named.
        """
        from tests.support import paths

        offender = paths.repo_root() / "src" / "GenAIRR" / "_s5f_loader.py"
        assert offender.exists(), f"premise: {offender} is still there"
        with pytest.raises(UnicodeDecodeError):
            offender.read_bytes().decode("cp1252")
        assert offender.read_text(encoding="utf-8"), "but UTF-8 reads it fine"


class TestTheseScansNameAFileTheSameWayEverywhere:
    """The second platform defect in these scans, after the encoding.

    A scan's allowlist is keyed by path, and `str(path.relative_to(root))`
    spells that path with the platform's separator -- so the entry
    `_compiled/plain.py`, written and reviewed on Linux, matched nothing
    on Windows, where the same file comes back as `_compiled\\plain.py`.
    The sanctioned site was therefore reported as an unsanctioned one on
    one of three platforms.

    Checked against a Windows path object rather than by running on
    Windows, so the guard holds wherever this suite runs.
    """

    def test_a_backslash_platform_produces_the_same_key(self):
        import pathlib

        root = pathlib.PureWindowsPath(r"C:\repo\src\GenAIRR")
        path = pathlib.PureWindowsPath(r"C:\repo\src\GenAIRR\_compiled\plain.py")
        assert _relative_key(path, root) == "_compiled/plain.py"

    def test_the_naive_spelling_is_the_one_that_broke(self):
        """Guard the premise: if `str()` already gave forward slashes
        there would be nothing here to fix, and the test above would pass
        against either implementation."""
        import pathlib

        root = pathlib.PureWindowsPath(r"C:\repo\src\GenAIRR")
        path = pathlib.PureWindowsPath(r"C:\repo\src\GenAIRR\_compiled\plain.py")
        assert str(path.relative_to(root)) == "_compiled\\plain.py"

    def test_every_sanctioned_key_is_spelled_the_way_a_scan_reports_one(self):
        """An entry no scan can ever match is an allowlist that silently
        does nothing, which is worse than no allowlist."""
        for key in TestNothingDerivesASeedByArithmetic.SANCTIONED:
            assert "\\" not in key, f"{key!r} cannot match a key the scan builds"
            assert key == key.replace("\\", "/")
