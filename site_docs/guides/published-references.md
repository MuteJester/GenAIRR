# Published germline sources and dossiers

<p class="lead">When a result has to be defensible a year from now, "we
used the current IMGT release" is not a reference. This page locks one
<em>named</em> published germline release, builds a runtime reference
from it, and writes a dossier that anyone can verify offline without
trusting you, the network, or the registry.</p>

Three names do all of it, from `GenAIRR.reference_sources`:

| Name | What it is for |
|---|---|
| `OgrdbSource` | Acquire and pin one explicit OGRDB germline-set release |
| `ReferenceBuild` | A runtime `Reference` plus the complete record of how it was made |
| `ReferenceSourceError` | Every failure here, with a stable `kind` and a `path` |

This is a different thing from
[Portable local references](portable-references.md). There you build a
reference from *your own* FASTA files. Here you pin a *published*
release and keep the evidence that you did.

## Start from the shipped build

One build ships inside the wheel, so the whole workflow is available
with no network at all:

```python
import GenAIRR as ga
from GenAIRR.reference_sources import ReferenceBuild

build = ReferenceBuild.bundled("human_igh_ogrdb_10_observed_v1")

result = (
    ga.Experiment.on(build.reference)
      .recombine()
      .run_records(n=1000, seed=0)
)
result.records[0]["v_call"]     # 'IGHV4-30-2*06' - IMGT names
```

Ask for a name that does not exist and the error lists the ones that
do, so the catalogue is always one call away:

```python
ReferenceBuild.bundled("nope")
# ReferenceSourceError(kind='missing_asset'):
#   no bundled reference build is named exactly
#   'human_igh_ogrdb_10_observed_v1'
```

!!! note "This is also the simplest route to IMGT allele names"
    The bundled build reports IMGT names (`IGHV4-30-2*06`), where
    `"human_igh"` reports ASC cluster names (`IGHVF6-G23*03`). If you
    are scoring a tool against an IMGT reference, that matters - see
    [Blinded handoff](../getting-started/export-results.md#blinded-handoff-strip-truth-columns).

## What a build carries

A `ReferenceBuild` is the reference *and* its provenance. Four
read-only parts:

| Member | Holds |
|---|---|
| `build.reference` | The runtime `Reference` you pass to `Experiment.on(...)` |
| `build.lock` | The pinned release: provider, `set_id`, every asset and its digest |
| `build.recipe` | The explicit projection - selection, decisions, defaults profile |
| `build.report` | What came out: counts, capabilities, claims, comparison, inventory |
| `build.manifest` | Member versions, hashes, and the generator that wrote it |

The manifest records the GenAIRR version that produced the build:

```python
build.manifest.builder
# Generator(name='GenAIRR', version='2.4.0')
```

That is deliberate. A reference built by a different version is a
different artifact even when the germline bytes are identical, because
the projection rules may have moved.

## Write a dossier, verify it offline

`write_dossier` takes a **directory**, not an archive:

```python
build.write_dossier("igh_dossier/")
```

It writes the named members plus the content-addressed blobs they
refer to:

```text
source.lock.json  recipe.json  report.json
reference.genairr-ref.json  manifest.json  annotations.json
v.fasta  d.fasta  j.fasta
<sha256>  <sha256>  ...
```

Anyone can check it with no network and no trust in you:

```python
from GenAIRR.reference_sources import ReferenceBuild

back = ReferenceBuild.read_dossier("igh_dossier/")
back.manifest.hashes == build.manifest.hashes   # True
```

`read_dossier` does not merely read the stored hashes. It regenerates
every derived file from the locked bytes and the recipe and compares
the result, so a tampered member is caught and located:

```python
# after changing one byte of reference.genairr-ref.json
ReferenceBuild.read_dossier("igh_dossier/")
# ReferenceSourceError(kind='digest_mismatch',
#                      path='/reference.genairr-ref.json')
```

This is the property worth having in a methods section: the dossier is
checkable by a reviewer who has only the directory.

## Pin a release yourself

`OgrdbSource.fetch` is the **only** operation in GenAIRR that uses the
network. It takes one explicit release - never "latest" - and writes a
lock that pins the exact bytes:

```python
from GenAIRR.reference_sources import OgrdbSource

source = OgrdbSource.fetch(
    set_id="Homo sapiens IGH",
    release=10,                                       # an int, never "latest"
    profile="ogrdb_api_v2_human_ig_observed_v1",
    cache="blobs/",                                   # content-addressed blob directory
)
source.write_lock("igh.lock.json")
```

`release` is an explicit positive integer; there is no "latest". A
`profile` names the interpretation rules applied to the release, and
`cache` is the content-addressed directory the downloaded bytes land
in. By default all four exports are acquired (`airr`, `airr_ex`,
`gapped`, `gapped_ex`) - a build compares both record exports, and V
junction anchors are accepted only when bound to the gapped ones.

Pass `expected_lock=` to reacquire a release and *require* its pinned
bytes: anything that changed under the same release is reported as
drift, and nothing is stored or replaced.

!!! warning "The download ignores proxy settings"
    `fetch` uses the standard-library HTTPS client, which does not
    honour `HTTPS_PROXY` and friends. Behind a proxy, download the
    exports with your own tools and record them with
    `OgrdbSource.from_files(...)`.

Three ways to get a source, depending on what you have:

| Call | Use when |
|---|---|
| `OgrdbSource.fetch(...)` | You are online and pinning a release for the first time |
| `OgrdbSource.from_files(...)` | You already hold the release's exports locally |
| `OgrdbSource.read(lock, cache=...)` | You have a lock and a populated cache - fully offline |

Then `inspect()` to see what the release holds, and `build(recipe)` to
project it into a `ReferenceBuild`:

```python
source = OgrdbSource.read("igh.lock.json", cache="blobs/")
source.inspect()                 # what is in the release
build = source.build(recipe)     # -> ReferenceBuild
```

`read` verifies every blob the lock needs against `cache` and never
touches the network; a missing or corrupt blob is an error naming the
asset and its pinned digest.

The recipe is explicit by design: which alleles are selected, which
decisions were taken, which defaults profile applies. Nothing about a
build is implied by the time you ran it.

## Errors

Everything on this surface raises `ReferenceSourceError` with two
fields worth branching on - a stable `kind` and a `path` locating the
problem:

```python
from GenAIRR.reference_sources import ReferenceSourceError

try:
    build = ReferenceBuild.read_dossier(path)
except ReferenceSourceError as exc:
    print(exc.kind)   # 'digest_mismatch', 'missing_asset', ...
    print(exc.path)   # '/reference.genairr-ref.json'
```

`kind` is part of the contract and safe to branch on; the message is
for people.
