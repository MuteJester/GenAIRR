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
