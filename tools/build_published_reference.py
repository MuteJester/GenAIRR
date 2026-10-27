#!/usr/bin/env python3
"""Generate the dossier of a reference build packaged with GenAIRR.

Maintainer tool (run from a source checkout; not part of the shipped wheel).  For one exact
build identifier it records the pinned provider exports under
``tests/data/reference_sources/ogrdb/`` as local files, together with the provider's data-terms
statement pages (the licence and the requested citation, recorded in the lock), builds the
reviewed recipe next to them with ``OgrdbSource.build``, writes the dossier into a new or empty
directory and checks it there with ``ReferenceBuild.read_dossier``.  It never uses the
network.  Acquisition receipts are machine-local records and are not written.

The dossier records the GenAIRR version that runs this script as its builder.  Regenerate the
packaged dossier whenever that version changes (after a release version bump) and pin the
printed manifest hash in ``src/GenAIRR/_reference_sources/bundled.py``.

Usage::

    # Write the starter into a scratch directory:
    python tools/build_published_reference.py --output ./human_igh_ogrdb_10_observed_v1

    # Regenerate the packaged starter (remove the old directory first):
    python tools/build_published_reference.py \\
        --output src/GenAIRR/data/published_references/human_igh_ogrdb_10_observed_v1
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import tempfile
from typing import NamedTuple

import GenAIRR
from GenAIRR._reference_sources import ogrdb
from GenAIRR._reference_sources.acquisition import SuppliedFile, ingest_files
from GenAIRR._reference_sources.digests import manifest_hash
from GenAIRR._reference_sources.lock import Terms
from GenAIRR.reference_sources import OgrdbSource, ReferenceBuild, ReferenceSourceError

REPO = Path(__file__).resolve().parent.parent
SOURCE_DATA = REPO / "tests" / "data" / "reference_sources"
EXPORT_FILES = {"airr": "airr.json", "airr_ex": "airr_ex.json", "gapped": "gapped.fasta",
                "gapped_ex": "gapped_ex.fasta"}
# OGRDB's data terms: each statement page kept as a fixture, and the section of the fixture
# manifest (sources.json) that records the URL it was retrieved from and what it states.
OGRDB_STATEMENTS = {"licence_statement": ("ogrdb/statements/licensing_statement.html", "data_licence"),
                    "citation_statement": ("ogrdb/statements/help.html", "citation")}


def _sources() -> dict:
    return json.loads((SOURCE_DATA / "sources.json").read_bytes().decode("utf-8"))
MANIFEST = "manifest.json"


class Published(NamedTuple):
    """A packaged build: the provider set release, its interpretation profile and its recipe."""

    set_id: str
    release: int
    profile: str
    recipe: str

    @property
    def exports(self) -> dict:
        directory = SOURCE_DATA / "ogrdb" / self.set_id / str(self.release)
        return {fmt: directory / name for fmt, name in EXPORT_FILES.items()}

    @property
    def terms(self) -> Terms:
        """The provider's data terms, as ``sources.json`` records them."""
        sources = _sources()
        return Terms(status="recorded", assets=tuple(OGRDB_STATEMENTS), licence=sources["data_licence"]["licence"],
                     citation=sources["citation"]["requested_by_provider"], reason=None)

    @property
    def recipe_path(self) -> Path:
        return SOURCE_DATA / self.recipe


BUILDS = {
    "human_igh_ogrdb_10_observed_v1": Published(
        set_id="9606.IGH_VDJ", release=10, profile="ogrdb_api_v2_human_ig_observed_v1",
        recipe="human_igh_ogrdb_10_observed_v1.recipe.json"),
}
DEFAULT_ID = "human_igh_ogrdb_10_observed_v1"


class GeneratorError(Exception):
    """A target the generator will not write into."""


def source(exact_id: str, *, cache) -> OgrdbSource:
    """The pinned exports of ``exact_id`` and the provider's statement pages, recorded as local files
    with the terms they state; ``cache`` is an empty scratch directory."""
    published = BUILDS[exact_id]
    files = {fmt: SuppliedFile(path) for fmt, path in published.exports.items()}
    sources = _sources()
    files.update({fmt: SuppliedFile(SOURCE_DATA / name, declared_url=sources[section]["statement_url"])
                  for fmt, (name, section) in OGRDB_STATEMENTS.items()})
    acquired = ingest_files(set_id=published.set_id, release=published.release, profile=published.profile,
                            files=files, cache=cache, identify=ogrdb.identify, terms=published.terms)
    return OgrdbSource._of(acquired)


def build(exact_id: str, *, cache) -> ReferenceBuild:
    """Build ``exact_id`` from its pinned exports, statements and recipe; ``cache`` is an empty scratch
    directory."""
    return source(exact_id, cache=cache).build(BUILDS[exact_id].recipe_path)


def check_target(target: Path) -> None:
    """``target`` must be absent with an existing parent, or an empty directory."""
    if target.is_symlink() or (target.exists() and not target.is_dir()):
        raise GeneratorError(f"{target} exists and is not a directory")
    if target.is_dir():
        if any(target.iterdir()):
            raise GeneratorError(f"{target} is not empty; the dossier is written into a new or empty directory")
    elif not target.parent.is_dir():
        raise GeneratorError(f"the parent directory of {target} does not exist")


def publish(result: ReferenceBuild, target: Path) -> ReferenceBuild:
    """Write the dossier members of ``result`` into ``target`` (new or empty), ``manifest.json`` last,
    and return the dossier checked there.  Receipts are not written."""
    check_target(target)
    names = [member.name for member in result.manifest.members] + [MANIFEST]
    with tempfile.TemporaryDirectory(prefix="genairr_dossier_") as scratch:
        staged = Path(scratch) / "dossier"
        result.write_dossier(staged)
        members = {name: staged.joinpath(*name.split("/")).read_bytes() for name in names}
    target.mkdir(exist_ok=True)
    for name in names:
        path = target.joinpath(*name.split("/"))
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "xb") as handle:
            handle.write(members[name])
    checked = ReferenceBuild.read_dossier(target)
    if checked.manifest != result.manifest:
        raise GeneratorError(f"the dossier written into {target} does not check as the build that was written")
    return checked


def summary(exact_id: str, result: ReferenceBuild) -> str:
    """The identifiers to pin: the typed manifest hash and the runtime reference hashes."""
    reference = result.reference
    return "\n".join([
        f"build {exact_id}",
        f"builder GenAIRR {result.manifest.builder.version}",
        f"manifest {manifest_hash(result.manifest)}",
        f"runtime execution {reference.execution_hash}",
        f"runtime catalogue {reference.catalogue_hash}",
        f"runtime models {reference.models_hash}",
        f"runtime artifact {reference.artifact_hash}",
    ])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--id", default=DEFAULT_ID, choices=sorted(BUILDS), help="exact build identifier")
    parser.add_argument("--output", required=True, type=Path, help="new or empty directory for the dossier")
    args = parser.parse_args(argv)
    try:
        check_target(args.output)
        with tempfile.TemporaryDirectory(prefix="genairr_sources_") as scratch:
            result = build(args.id, cache=Path(scratch) / "cache")
        checked = publish(result, args.output)
    except (GeneratorError, ReferenceSourceError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(summary(args.id, checked))
    print(f"checked with GenAIRR {GenAIRR.__version__}; pin the manifest hash above in "
          "src/GenAIRR/_reference_sources/bundled.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
