# tools/

Maintainer tooling for GenAIRR. Not part of the shipped `GenAIRR` wheel —
these scripts are run from a source checkout.

## `build_imgt_configs.py`

Fetches IMGT V-QUEST germline FASTA and builds structural GenAIRR
cartridges (`<SPECIES>_<LOCUS>_IMGT.pkl`) via
`GenAIRR.ReferenceCartridgeBuilder`.

```bash
# Build one species into a scratch dir (never touches the bundled set):
python tools/build_imgt_configs.py --output-dir ./built_cartridges --species Mus_musculus

# Probe locus availability only:
python tools/build_imgt_configs.py --output-dir ./built_cartridges --dry-run
```

**Structural output.** Built cartridges carry the germline V/D/J
catalogue plus IMGT V-subregion (FWR/CDR) annotations, but **no**
data-derived empirical distributions (trim / NP length / NP base model /
gene usage) — those cannot be inferred from germline FASTA. At simulation
time the engine uses uniform defaults for those parameters. Fit real
distributions into the returned `DataConfig`, or estimate them with the
builder's `estimate_*` methods from your own AIRR data, if you need
empirically-grounded parameters. Only the bundled human IGH/IGK/IGL and
TCRB cartridges ship with real data-derived distributions.

`--output-dir` is **required**: a freshly-built structural cartridge
differs from the shipped `src/GenAIRR/data/builtin_dataconfigs/` set, so
the tool never overwrites it by default. Point `--output-dir` there
explicitly only when you deliberately intend to regenerate the bundled
cartridges (and expect the golden baselines to move).

## `build_published_reference.py`

Generates the packaged archive of a reference build shipped with GenAIRR
(currently `human_igh_ogrdb_10_observed_v1`). It records the pinned OGRDB
exports under `tests/data/reference_sources/ogrdb/` and OGRDB's licensing and
help pages (`ogrdb/statements/`, recorded in the lock as the data terms), builds
the reviewed recipe next to them (`tests/data/reference_sources/<id>.recipe.json`)
with `OgrdbSource.build`, writes the dossier into a private temporary directory,
checks it with `ReferenceBuild.read_dossier`, and writes its members as one
deterministic zip archive (entries sorted by name, stored uncompressed, fixed
times and modes; no acquisition receipts). It prints the hashes to pin. It never
uses the network.

```bash
# Write the starter archive into a scratch location:
python tools/build_published_reference.py --output ./human_igh_ogrdb_10_observed_v1.zip

# Regenerate the packaged archive (remove the old one first):
python tools/build_published_reference.py \
    --output src/GenAIRR/data/published_references/human_igh_ogrdb_10_observed_v1.zip
```

The dossier records the GenAIRR version that ran the script as its builder,
so its manifest hash and the archive digest change with the version. After a
release version bump, regenerate the packaged archive with the bumped version
installed, then pin the printed manifest hash in
`src/GenAIRR/_reference_sources/bundled.py` and in
`tests/test_published_reference_resources.py` (`MANIFEST_PIN`), and the printed
archive digest in the same test file (`ARCHIVE_PIN`). The runtime execution,
catalogue and models hashes do not change with the version.
