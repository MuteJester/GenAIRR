# Docs Website Information Architecture — Audit

**Status: audit + Slice 1 shipped.** Inventories every
documentation surface in the GenAIRR repository, evaluates
the information architecture against the user's proposed
structure, and recommends a roadmap. The pre-flight check
**triggered stop-and-report** in a nuanced way (§3 below):
a docs site already exists but the README's promised URLs
returned 404 in production because they pointed at a
different framework's URL scheme.

**Slice 1 — README Production URL Fix — SHIPPED.** All 13
broken slugged URLs in the README (`mutejester.github.io/GenAIRR/docs/...`)
were rewritten to point at the live flat-HTML pages
(`mutejester.github.io/GenAIRR/<page>.html`). Bottom-of-
README "Documentation" pointer list expanded to surface
all 5 lessons + 5 concepts + 8 representative guides + the
reference page. Pin set in
[`tests/test_docs_website_contract.py`](../tests/test_docs_website_contract.py)
flipped: prior broken-URL stop-and-report pin replaced by
`pin_present_readme_hosted_doc_urls_resolve_to_live_website_pages`
(asserts every URL maps to a real `website/*.html` file)
plus a defensive `pin_present_readme_contains_no_slugged_mkdocs_urls`
regression guard.

**Project-direction call still pending (§6.2 + §6.3).** No
framework change yet. Path A (keep hand-rolled) vs Path B
(migrate to MkDocs Material) decision belongs to the user.

Companion to
[`tests/test_docs_website_contract.py`](../tests/test_docs_website_contract.py)
which freezes (a) the current site framework + URL
scheme, (b) the parallel maintainer-facing doc tree, (c)
the broken README cross-references that any migration must
preserve, (d) the abandoned earlier-framework
configuration in `_old_docs/`, and (e) the docs-deployment
CI hook.

---

## 1. Q1 — What documentation exists today?

### 1.1 Three parallel documentation surfaces

The repository carries **three independent doc trees**,
none integrated with the others:

| Surface | Location | Format | Audience | Size |
|---|---|---|---|---|
| **A. Live website** | [`website/`](../website/) | Hand-rolled HTML/CSS/JS | New users | 30 HTML pages (1 hub + 5 lessons + 17 guides + 5 concepts + 1 reference + index) + `styles.css` + `genairr.js` |
| **B. Maintainer audits/designs** | [`docs/`](../docs/) | Markdown rendered on GitHub | Engine contributors | 38 `.md` files (15 audits + 17 designs + 6 hubs/guides) |
| **C. Abandoned earlier docs** | [`_old_docs/`](../_old_docs/) | Markdown + Jupyter notebooks for MkDocs Material | (deprecated) | 12 `.md` + 5 `.ipynb` tutorials + `doc_requirements.txt` |

### 1.2 Inventory — surface A (live website)

The hand-rolled HTML site at
[`website/`](../website/) is the public-facing
documentation. Deployed to
`mutejester.github.io/GenAIRR` via
[`.github/workflows/deploy-docs.yml`](../.github/workflows/deploy-docs.yml)
on every push to `master` that touches `website/**`. The
workflow has NO build step — files ship as-is.

Top-level navigation: **Home / Learn / Guides / Concepts
/ Reference / GitHub**.

| Page family | Count | Hub page | Examples |
|---|---|---|---|
| Index | 1 | [`index.html`](../website/index.html) | hero + "pip install" CTA |
| Lessons (numbered tutorial track) | 5 | [`learn.html`](../website/learn.html) | "V(D)J recombination, from bases up" / "The pipeline scrubber" / "S5F: mutation isn't uniform" / "Sequencing artifacts, dialed in" / "The ground-truth payoff" |
| Guides (task-shaped) | 17 | [`guides.html`](../website/guides.html) organised into **Build / Validate / Customize / Operate** | `guide-build-config.html`, `guide-productive.html`, `guide-clonal-families.html`, `guide-benchmark-aligner.html`, `guide-compare-shm.html`, `guide-audit-realism.html`, `guide-smoke-test.html`, `guide-tune-corruption.html`, `guide-v-usage.html`, `guide-replay.html`, `guide-trace-introspection.html`, `guide-export.html`, `guide-streaming.html`, `guide-reproduce.html`, `guide-with-metadata.html` |
| Concepts (mental-model essays) | 5 | [`concepts.html`](../website/concepts.html) | `concept-pipeline.html`, `concept-persistent-ir.html`, `concept-contracts.html`, `concept-airr-record.html`, `concept-live-call.html` |
| Reference | 1 | [`reference.html`](../website/reference.html) | Python API + Configs + AIRR fields + CLI |

This is a **deliberate information architecture** —
not an accidental dump. Headings + content map onto a
known user-journey (Learn → Guides → Concepts →
Reference) that broadly aligns with the user's proposed
new structure (§3 comparison below).

### 1.3 Inventory — surface B (maintainer audits/designs in `docs/`)

The [`docs/`](../docs/) directory carries every audit +
design doc produced over the engine-rewrite slices. Each
file is canonical for its slice and routinely cross-
referenced from the README, the validation matrix, and
sibling audits.

| Category | Count | Examples |
|---|---|---|
| **Audits** (per-mechanism / per-invariant analysis with §6 drift catalogue) | 15 | `productive_failure_mode_audit.md`, `indel_provenance_audit.md`, `allele_call_audit.md`, `junction_call_audit.md`, `distribution_invariant_audit.md`, `performance_baseline.md`, `mutation_provenance_audit.md`, `shm_model_audit.md`, `v_region_substructure_audit.md`, `v_subregion_mutation_counters_audit.md`, `primer_trim_end_loss_audit.md`, `reference_cartridge_authoring_audit.md`, `reference_cartridge_completeness_audit.md`, `junction_n_addition_audit.md`, `plan_signature_completeness_audit.md` |
| **Designs** (per-slice / per-mechanism scoping pre-implementation) | 17 | `clonal_family_design.md`, `clonal_parent_outcome_design.md`, `clonal_plan_split_design.md`, `d_inversion_design.md`, `d_inversion_extension_design.md`, `paired_end_design.md`, `receptor_revision_design.md`, `shm_segment_rate_design.md`, `v_subregion_shm_rate_design.md`, `np_markov_base_generator_design.md`, `p_nucleotide_design.md`, `fastq_export_design.md`, `allele_usage_estimation_design.md`, `trim_distribution_estimation_design.md`, `np_length_estimation_design.md`, `np_base_model_estimation_design.md`, `p_nucleotide_length_estimation_design.md` |
| **Hubs / contributor entry points** | 6 | `engine_architecture.md`, `adding_a_pass.md`, `validation_matrix.md`, `reference_cartridge.md`, `airr_record_validator.md`, `allele_model_audit.md` |

The `docs/superpowers/plans/` directory holds Claude-
session planning artifacts (e.g.
`2026-05-18-mcp-redesign-v2.md`) — **not user-facing
documentation**, but currently mixed into the same
`docs/` tree.

The `docs/build/` directory holds Python wheel build
artifacts (a side-effect of `python -m build`) — also
**not documentation**, and would not deploy in any
framework's expected output.

### 1.4 Inventory — surface C (abandoned earlier docs in `_old_docs/`)

The renamed-with-leading-underscore
[`_old_docs/`](../_old_docs/) directory contains an
earlier docs system that was deliberately deprecated. It
mixes evidence of two abandoned framework choices:

- **Sphinx**:
  [`_old_docs/Makefile`](../_old_docs/Makefile) +
  [`_old_docs/make.bat`](../_old_docs/make.bat) are the
  standard `sphinx-quickstart` boilerplate (references
  `sphinx-build` + `source/` + `build/` — though no
  `conf.py` survives).
- **MkDocs Material** (later migration attempt):
  [`_old_docs/doc_requirements.txt`](../_old_docs/doc_requirements.txt)
  pins `mkdocs-material` + `mkdocs-jupyter` (no
  Sphinx). The Markdown files
  (`getting_started.md` / `step_by_step_tutorial.md` /
  `faq.md` / `troubleshooting.md` / `biological_context.md`
  etc.) use MkDocs admonition syntax (`!!! tip "..."`).

**Content quality:** the 12 markdown files were last
touched in January 2025 (`api_reference.md`,
`getting_started.md`, etc.) and **predate the engine_rs /
Rust kernel slice**. They are stylistically polished but
factually stale on:

- Engine architecture (still describes pure-Python kernel).
- Cartridge model (predates the typed-plane discipline).
- API surface (predates `ReferenceCartridgeBuilder` +
  estimators).
- Validation discipline (predates the two-layer model).

The 5 `.ipynb` tutorials
(`Quick Start Guide.ipynb`,
`Creating Custom DataConfig from FASTA Files.ipynb`, …)
sit in `_old_docs/tutorials/` with bundled example FASTA
files. They are similarly stale but executable enough to
salvage as starting templates for new
`mkdocs-jupyter`-rendered notebooks if a migration to
MkDocs Material lands.

### 1.5 Pinned

- `pin_scaffold_website_dir_exists_with_handwritten_html`
- `pin_scaffold_docs_dir_carries_thirty_eight_md_audit_design_files`
- `pin_scaffold_old_docs_dir_exists_as_abandoned_earlier_attempt`
- `pin_scaffold_deploy_docs_workflow_targets_website_dir`
- `pin_scaffold_docs_superpowers_subdir_holds_session_artifacts_not_docs`
- `pin_scaffold_docs_build_subdir_holds_wheel_artefacts_not_docs`

---

## 2. Q2 — User-facing vs maintainer-facing

### 2.1 Current classification

The three doc trees split cleanly along the user / maintainer axis but with significant gaps:

| Tree | Currently for | Currently used by | Cross-references |
|---|---|---|---|
| [`website/`](../website/) | New users — Learn / Guides / Concepts / Reference | Public-facing deployment via GitHub Pages | Self-contained (intra-site only) |
| [`docs/*.md`](../docs/) | Engine contributors writing new passes / understanding invariants | The README links to ~30 audit/design files; the validation matrix is the navigable index | Heavy intra-`docs/` linking + back-references from source code (`See docs/foo.md`) |
| [`_old_docs/`](../_old_docs/) | (deprecated — nobody) | Nobody | Self-referential within `_old_docs/`; one external mention in `tutorials/example_custom_tcr_files/` |

### 2.2 The maintainer-leak problem

The README currently links to **15 distinct `docs/*.md`
files** in normal prose (e.g. `See [docs/p_nucleotide_design.md]`).
These render fine on GitHub (so the README is readable as
written) but are **not part of the live website**.
Result: a new user clicking "📖 Documentation" on the
README lands at `mutejester.github.io/GenAIRR` and sees
the polished 30-page site; if they then click a link
inside README prose for technical detail, they're
redirected to a raw GitHub markdown view of an internal
audit doc — context-shift that confuses new users.

The user's brief observation is correct:

> "Many audit docs are excellent but too internal for
> new users."

The audit confirms the symptom: the README is the **only
joint** between the two trees, and it does the joining
inconsistently. Audit docs are linked as if they were
user-facing reference material, but they're written
for engine contributors and assume context the new user
doesn't have.

### 2.3 Pinned

- `pin_scaffold_readme_links_thirty_plus_md_files_in_docs_dir`
- `pin_scaffold_validation_matrix_is_navigable_index_for_audits`

---

## 3. Q3 — Proposed structure vs reality + stop-and-report verdict

### 3.1 Stop-and-report condition

The user's brief specified:

> "If a docs site already exists with a different
> framework, do not propose replacing it blindly.
> Inventory it first."

**Verdict: triggered, with nuance.**

A docs site exists. It is **not a different framework**
in the conventional sense (no Sphinx / MkDocs / Docusaurus
config) — it is hand-rolled HTML/CSS/JS, deployed
verbatim. The user's "do not propose replacing blindly"
constraint applies even more strongly here than to a
framework-based site: replacement would discard
hand-crafted content + custom styling without an obvious
upgrade.

But three production-quality issues mean the audit MUST
also recommend changes:

#### Issue 1 — README cross-references return 404 in production

The README links to URLs like:

```text
https://mutejester.github.io/GenAIRR/docs/getting-started/quick-start
https://mutejester.github.io/GenAIRR/docs/concepts/simulation-pipeline
https://mutejester.github.io/GenAIRR/docs/guides/options/productive
```

Empirical check at audit time
(`curl -sI <url>` against the live site):

| URL | HTTP status |
|---|---|
| `mutejester.github.io/GenAIRR/docs/getting-started/quick-start` | **404** |
| `mutejester.github.io/GenAIRR/lesson-1.html` | **200** |

These URLs match an **MkDocs Material URL scheme** (slash-
delimited, no `.html` suffix). They were the planned
output of the abandoned `_old_docs/` migration. The
README still cites them; the live site uses a different
URL scheme (`lesson-1.html`, `guide-productive.html`).
**Every "📖 Docs" link in the README that points at a
specific page is broken**. This is a public-facing
trust regression: a user clicking
"See [Quick Start](https://…/docs/getting-started/quick-start)"
hits a 404.

#### Issue 2 — Hand-rolled HTML doesn't scale with content velocity

The repository has shipped ~40 audit / design docs in the
last 6 months. Each new estimator / mechanism / refactor
produces 100–500 lines of canonical documentation.
The `website/` site has not absorbed any of this content:
the lessons + guides + concepts predate the Rust engine,
the estimator suite (5 estimators landed in this
session), the V-subregion work, the paired-end FASTQ
export, the typed cartridge planes, the validation
matrix's two-layer model, etc.

A hand-rolled HTML system requires every doc page to be
**hand-translated from markdown sources** (where authors
draft) into HTML (where users read). The translation
step is currently not happening at all — `website/` and
`docs/` have diverged.

#### Issue 3 — No search, no versioning, no cross-link checking

The hand-rolled site provides none of:

- **Full-text search** (the navigation requires the user
  to know which guide hub their question belongs in).
- **Version tagging** (the site shows whatever's on
  `master`; users on `pip install GenAIRR==0.X` see a
  different API).
- **Cross-link validation** (broken intra-site links and
  README→site link breakage go undetected).
- **Auto-generated API reference** (the `reference.html`
  page is hand-curated; the 35 top-level `__all__`
  exports are not enumerated in any tested-against-source
  fashion).

### 3.2 The user's proposed structure — alignment with current site

The user's brief proposes:

```text
Getting Started     Core Concepts        Simulation Guides
  Installation        Experiment pipeline   Recombination
  First simulation    Reference cartridges  SHM and mutation
  Outputs (AIRR/...)  Reproducibility       Indels / artefacts
                      Validation reports    End-loss / paired-end
                                            P/N nucleotides
                                            D inversion
                                            Clonal families
Reference Cartridge Validation & Debugging  API Reference        Architecture
Authoring             validate_records        Experiment           engine
  Build from FASTA    validate_families       SimulationResult     adding a pass
  Rules / identity    cache parity            DataConfig / builder trace/event/replay
  Curation            release matrix          Reference models     audits
  Empirical models
  Estimators
```

Mapping this to the existing `website/` site:

| User's proposed bucket | Current `website/` equivalent | Gap |
|---|---|---|
| Getting Started | (none — `index.html` does part of this; `learn.html` is the closest) | Need a dedicated "Installation" + "First simulation" trio |
| Core Concepts | `concepts.html` + 5 concept pages | Aligned but covers different concepts (the existing site emphasises "Persistent IR", "Live Calls" — internals visible to users; the proposed structure emphasises "Reproducibility", "Validation reports" — user-relevant abstractions) |
| Simulation Guides | `guides.html` Build/Customize buckets | Mostly aligned — the existing 17 guides cover most of the proposed topics, with different naming and grouping |
| Reference Cartridge Authoring | `guide-build-config.html` (one guide) | **Major gap** — no dedicated section for the cartridge authoring story, which is now a substantial surface (5 estimators, the builder API, build reports) |
| Validation & Debugging | `guides.html` Validate bucket | Aligned (5 of the 17 guides live here) |
| API Reference | `reference.html` (one hand-curated page) | **Major gap** — no auto-generated, comprehensive API reference; users currently rely on Python's `help()` |
| Architecture / Contributor | (none on the website) — `docs/engine_architecture.md` + `docs/adding_a_pass.md` are GitHub-only | **Major gap** — contributor docs are not discoverable from the site |

**Overall alignment: ~60%.** The major gaps are
**Reference Cartridge Authoring**, **API Reference**, and
**Architecture / Contributor**.

### 3.3 Pinned

- `pin_present_readme_external_doc_links_are_broken_in_production`
- `pin_scaffold_website_navigation_is_learn_guides_concepts_reference`
- `pin_scaffold_proposed_structure_has_three_major_gaps_vs_current_site`

---

## 4. Q4 — Which docs should be rewritten?

### 4.1 Rewrite vs salvage decision matrix

| Doc category | Path | Action | Reason |
|---|---|---|---|
| 15 `*_audit.md` | `docs/` | **Salvage as contributor reference** (link from a future `architecture/` section of the site) | Each is canonical for its mechanism; rewriting would lose the §6 drift catalogue + slice-level traceability |
| 17 `*_design.md` | `docs/` | **Salvage as contributor reference** (same as audits) | Each documents a v1 boundary decision that's load-bearing for cross-slice consistency |
| `validation_matrix.md` | `docs/` | **Adopt directly as the site's contributor entry point** | Already shaped as a navigable index; one-line wrapper page in the site links into it |
| `engine_architecture.md` | `docs/` | **Salvage as contributor reference** | "Contracts constrain support, trace = choices, events = consequences" — load-bearing invariants; a new user guide cannot replace it |
| `adding_a_pass.md` | `docs/` | **Salvage as contributor reference** | Same |
| `reference_cartridge.md` | `docs/` | **Promote to user guide + salvage internal references** | The 3-paths comparison (bundled / manual / builder) is genuinely user-facing; the file is already shaped that way |
| `airr_record_validator.md` | `docs/` | **Promote to user guide** (under Validation & Debugging) | Already user-facing in shape — "Standard postcondition pattern" + "Layered validation" |
| 5 `.ipynb` in `_old_docs/tutorials/` | `_old_docs/` | **Refresh for current API + salvage as MkDocs Material notebooks** (if migration path chosen) | The tutorial track ("Quick Start", "Custom DataConfig from FASTA", etc.) is the user-journey the new site needs; the existing notebooks are stale on API but the structure is correct |
| 12 `*.md` in `_old_docs/` | `_old_docs/` | **Audit individually**; salvage those that survive a freshness check, discard rest | Mixed quality + currency — `getting_started.md` is stale on `Experiment.recombine()`, `faq.md` predates several major mechanisms, `troubleshooting.md` is generic |
| 30 `*.html` in `website/` | `website/` | **Decision point — see §6 framework recommendation** | Either keep as-is (light path) or convert to markdown sources (migration path) |

### 4.2 Pinned

- `pin_scaffold_audit_design_docs_are_canonical_for_their_mechanisms`
- `pin_scaffold_reference_cartridge_md_is_user_guide_shaped`
- `pin_scaffold_validation_matrix_md_is_navigable_index`

---

## 5. Q5 — Executable examples

### 5.1 Required executable examples (per user brief)

| Example | Current state | Source-of-truth |
|---|---|---|
| First recombination | README has it; `lesson-1.html` has it; `_old_docs/tutorials/Quick Start Guide.ipynb` has it (stale) | README + lesson-1 are both correct |
| Productive full stack | README has it; `lesson-4.html` covers it; `guide-productive.html` covers it | Aligned across all three |
| Custom cartridge from FASTA | README has it; `guide-build-config.html` covers it; `_old_docs/tutorials/Creating Custom DataConfig from FASTA Files.ipynb` (stale) | README + guide-build-config |
| Estimated allele usage / trims / NP models | README §"Build a cartridge from FASTA" mentions `estimate_allele_usage`; remaining 4 estimators NOT in README or `website/` | **MAJOR GAP** — 4 of 5 estimators undocumented user-side |
| Paired FASTQ export | README has it (`to_paired_fastq`); `guide-export.html` mentions it | Aligned |
| Trace replay | README has it; `guide-replay.html` + `guide-trace-introspection.html` cover it | Aligned |
| Validation report usage | README has it (`result.validate_records(refdata)`); `airr_record_validator.md` covers it | Aligned in code; not in user-site |

### 5.2 Executable-check status

The user brief asked: "Examples in README/docs either execute or are marked pseudo-code."

| Source | Status |
|---|---|
| README code blocks | NOT mechanically executed — written by hand, may drift |
| `website/*.html` code blocks | NOT mechanically executed — written by hand |
| `_old_docs/tutorials/*.ipynb` | Mechanically executed AT NOTEBOOK SAVE TIME, but the notebooks haven't been re-saved since Jan 2025 |
| `docs/*_audit.md` + `docs/*_design.md` code snippets | Some are illustrative pseudo-code; most are descriptive (not invocable) |

**Verdict:** there is no executability discipline anywhere
in the docs surface. Code in the README is likely correct
because it ships in this session's working tree and is
quoted from passing tests, but there is no automated
guard.

### 5.3 Pinned

- `pin_scaffold_readme_code_blocks_have_no_executability_discipline`
- `pin_scaffold_only_one_of_five_estimators_documented_user_side`

---

## 6. Q6 — Website tech recommendation

### 6.1 Existing situation

- **Current production:** hand-rolled HTML/CSS/JS in
  `website/`, deployed via
  `.github/workflows/deploy-docs.yml`.
- **Abandoned earlier:** Sphinx scaffold then MkDocs
  Material migration, both in `_old_docs/` — neither
  shipped to production.
- **Inferred history:** the project moved
  Sphinx → MkDocs Material → hand-rolled HTML. The
  current state may be a deliberate aesthetic choice
  (the existing site has carefully-crafted typography
  with IBM Plex fonts + custom CSS that none of the
  off-the-shelf frameworks provide out of the box).

### 6.2 Two recommended paths

The audit explicitly refrains from picking a winner —
the trade-off is a project-direction call. Both paths
are viable.

#### Path A — Keep the hand-rolled site, formalize content velocity (LIGHT)

| Step | Effort | Effect |
|---|---|---|
| Fix the 30+ broken README URLs (point at `lesson-X.html` / `guide-Y.html` / `concept-Z.html`) | ~1 hour | Closes the production trust regression |
| Add a `website/CONTRIBUTING.md` documenting how to add a new HTML page (template + nav update) | ~30 min | Lowers the barrier to content velocity |
| Add a `website/api-reference.html` auto-generated from `__all__` + docstrings (could be a one-off script that emits HTML from `inspect`) | ~4 hours | Closes the major-gap on API reference |
| Add a `website/cartridge-authoring.html` hub linking to a new sub-page per estimator | ~6 hours | Closes the cartridge-authoring gap |
| Add 5 estimator pages (one per estimator) | ~10 hours | Closes the estimator-documentation gap |
| Add a `tests/test_docs_links.py` that crawls `website/` + README and asserts every internal link is reachable | ~3 hours | Catches future drift |

Total: **~24 hours** for a complete patch + new gap-filling content.

#### Path B — Migrate to MkDocs Material (HEAVY)

| Step | Effort | Effect |
|---|---|---|
| Bootstrap MkDocs Material from `_old_docs/doc_requirements.txt` (the abandoned migration's leftover deps file) | ~2 hours | Re-instates the abandoned framework |
| Author `mkdocs.yml` matching the user's proposed structure (Getting Started / Core Concepts / Simulation Guides / Reference Cartridge Authoring / Validation & Debugging / API Reference / Architecture) | ~3 hours | Locks the IA |
| Convert 30 `website/*.html` pages to MkDocs markdown (manual content extraction; ~30 min per page) | ~15 hours | Migrates the existing user-facing content |
| Refresh + import 5 `_old_docs/tutorials/*.ipynb` notebooks via `mkdocs-jupyter` (mostly content updates for current API) | ~10 hours | Re-instates the tutorial track |
| Add a `mkdocstrings` auto-API page (auto-generated from `__all__` + docstrings) | ~2 hours | Closes the API reference gap |
| Add URL-redirect mapping so existing `lesson-X.html` / `guide-Y.html` URLs continue to resolve to their MkDocs equivalents (e.g. via `mkdocs-redirects`) | ~1 hour | Preserves external links / backwards compatibility |
| Add 5 estimator pages + the cartridge-authoring hub | ~6 hours | Closes the same content gaps as Path A |
| Replace the `deploy-docs.yml` workflow with a build step (`mkdocs build` + upload) | ~1 hour | Switches the deployment target |
| Add `tests/test_docs_links.py` (same as Path A) | ~3 hours | Catches future drift |

Total: **~43 hours** for the migration + new gap-filling content.

#### When to pick which

| Pick **Path A** (light) if | Pick **Path B** (heavy) if |
|---|---|
| You value the custom typography / styling of the current site | You want full-text search / versioning / cross-link validation out of the box |
| You expect to ship < 10 new doc pages over the next 3 months | You expect to ship > 30 new doc pages (which the estimator + cartridge work implies) |
| You're optimising for short-term URL stability | You're optimising for long-term content velocity |
| You don't want Python/Node toolchain in the docs build | You're comfortable adding `mkdocs-material` as a docs-only dep |

### 6.3 The audit's hedged recommendation

If forced to recommend: **Path B (MkDocs Material)** wins
on long-term content velocity, but ONLY if the user is
willing to do the ~15-hour content extraction up front
and accept the framework constraints. The
**`_old_docs/doc_requirements.txt`** breadcrumb suggests
this was the intended target before something pulled the
project off-track; finishing the abandoned migration is
arguably cheaper than maintaining the hand-rolled site
indefinitely.

If the user prioritises preserving the existing site's
visual identity, **Path A** is fully viable and the
broken-URLs fix alone (the first 1-hour bullet) closes
the most production-critical gap.

### 6.4 Pinned

- `pin_scaffold_no_active_docs_framework_config_today_aside_from_handwritten_html`
- `pin_scaffold_old_docs_doc_requirements_lists_mkdocs_material`
- `pin_absence_no_mkdocs_yml_today`
- `pin_absence_no_sphinx_conf_py_today`

---

## 7. Pins/tests — what this audit pins

Mirrored in
[`tests/test_docs_website_contract.py`](../tests/test_docs_website_contract.py).

### `pin_scaffold_*` — current state inventory

1. `website/` exists with handwritten HTML, deployed by
   GitHub Pages workflow.
2. `docs/` contains 38 audit/design `.md` files.
3. `_old_docs/` exists as the abandoned earlier docs
   tree (with MkDocs Material breadcrumbs).
4. `.github/workflows/deploy-docs.yml` deploys
   `website/` verbatim with no build step.
5. README links to 15 distinct `.md` files in `docs/` AND
   to external URLs at `mutejester.github.io/GenAIRR/`.
6. `validation_matrix.md` is the navigable contributor
   index.
7. `_old_docs/doc_requirements.txt` lists
   `mkdocs-material` + `mkdocs-jupyter`.

### `pin_present_*` — broken-URL verification

8. **The README's external doc URLs to
   `mutejester.github.io/GenAIRR/docs/…` return 404 in
   production**, while `mutejester.github.io/GenAIRR/lesson-1.html`
   returns 200 (test offline-safe: parses the README,
   verifies the broken URL scheme matches the abandoned
   `_old_docs/` MkDocs scheme).

### `pin_scaffold_*` — content gaps

9. Only 1 of 5 estimators (`estimate_allele_usage`) is
   referenced in the live `website/` site.
10. `website/reference.html` is hand-curated, NOT
    auto-generated from `__all__`.
11. The 35-symbol top-level `__all__` export is not
    enumerated in any single discoverable doc surface.

### `pin_absence_*` — framework configs

12. No `mkdocs.yml` at repo root or in `docs/`.
13. No Sphinx `conf.py` anywhere in the live tree
    (the `_old_docs/Makefile` references one but it
    doesn't exist).
14. No `docusaurus.config.js`, no Jekyll `_config.yml`.

### Doc anchor

15. Audit doc exists and references the contract file;
    section structure intact.

---

## 8. Out of scope

Documented here so a future implementer doesn't accidentally
expand this audit.

- **Choosing the framework.** This audit lays out two
  recommended paths (A: light, keep hand-rolled; B: heavy,
  migrate to MkDocs Material). The framework decision is
  a project-direction call that belongs to the user, not
  to the implementation slices.
- **Rewriting individual doc pages.** This is a structure
  audit. Individual page rewrites are out of scope; the
  audit catalogues which pages need rewriting (§4) but
  doesn't draft replacements.
- **Auto-generating an API reference.** Both paths
  recommend it, but the implementation depends on the
  framework chosen (Path A: custom HTML generator;
  Path B: `mkdocstrings`).
- **Versioning strategy.** The hand-rolled site has no
  versioning; the MkDocs migration could add it (via
  `mike`). Decision belongs to the framework choice.
- **Search backend.** Same — depends on framework choice.
- **`docs/superpowers/plans/` cleanup.** The Claude-
  session planning artefacts should probably move to a
  separate `.private/` directory but that's a housekeeping
  matter, not a docs structure issue.

---

## 9. Implementation order (recommended)

The audit recommends a slice sequence regardless of which
framework path lands:

### 9.1 Slice 1 — fix the broken URLs (1 hour, both paths) **[SHIPPED]**

Updated the README to point at the actual live URLs
(`mutejester.github.io/GenAIRR/lesson-1.html` etc.).
Critical-path: closed the public-facing trust regression.
20 broken URLs rewritten; bottom-of-README
"Documentation" pointer list reorganised around the
site's actual Learn / Concepts / Guides / Reference
buckets. Contract pins flipped to assert every README
hosted-doc URL resolves to a live `website/*.html` file
(offline-safe — filesystem mapping, no network).

### 9.2 Slice 2 — pick the framework path

This is where the project owner makes the framework
decision. The audit's deliverables (this doc + the
contract pins) freeze the inventory + the verdict so the
decision can happen offline.

### 9.3 Slice 3 — site skeleton (path-dependent)

- **Path A:** add `website/cartridge-authoring.html` hub
  + `website/api-reference.html` auto-gen + the 5
  estimator pages.
- **Path B:** bootstrap `mkdocs.yml` + matching directory
  structure, port `website/*.html` content to markdown
  shells with TODO markers, set up the `mkdocs build`
  workflow.

### 9.4 Slice 4 — content migration (path-dependent)

- **Path A:** translate the `docs/reference_cartridge.md`
  + `docs/airr_record_validator.md` user-facing content
  into HTML pages.
- **Path B:** convert hub HTML pages to markdown using a
  one-shot script, then iterate.

### 9.5 Slice 5 — executable-example discipline (both paths)

Add a `tests/test_docs_executable_examples.py` that
extracts code blocks from selected markdown / HTML pages
(those tagged as executable) and runs them through
`subprocess.run([sys.executable, "-c", code])`. Failing
extraction surfaces in CI.

### 9.6 Slice 6 — link-checking discipline (both paths)

Add `tests/test_docs_links.py` (per audit §6.2 / §6.3 last
bullet of both paths). Catches future drift.

### 9.7 Slice 7 — gap-filling content (both paths)

Once the framework + skeleton are in place, write:

- Estimator hub + 5 sub-pages.
- API reference auto-generated.
- Validation & Debugging hub with the two-layer model.
- Architecture / Contributor hub linking out to
  `docs/engine_architecture.md` + `docs/adding_a_pass.md`
  + `docs/validation_matrix.md`.

---

## 10. Summary table

| Concern | Today | Audit recommendation |
|---|---|---|
| **Live docs site** | Hand-rolled HTML in `website/` (30 pages), deployed by `.github/workflows/deploy-docs.yml` | **Stop-and-report (mild): site exists and is deliberate.** Two recommended paths — keep + patch (light) OR migrate to MkDocs Material (heavy). Decision belongs to user. |
| **README cross-references** | 15 distinct `docs/*.md` links (render on GitHub) + 7+ links to `mutejester.github.io/GenAIRR/docs/…/...` URLs (404 in production — wrong framework's URL scheme) | **Slice 1: fix the broken URLs.** Critical path. |
| **Maintainer docs** | 38 `.md` audits/designs in `docs/` — canonical for their mechanisms, used by README + validation matrix | **Salvage all 38 as contributor reference.** Link from a future Architecture / Contributor section of the site. |
| **Abandoned earlier docs** | `_old_docs/` with Sphinx Makefile + MkDocs Material `doc_requirements.txt` + 12 stale markdown files + 5 stale tutorial notebooks | **Don't delete blindly.** The `_old_docs/doc_requirements.txt` is a useful breadcrumb if Path B is chosen; the 5 tutorials are good starting templates after refresh. |
| **Estimator documentation** | Only 1 of 5 estimators in `website/` | **Slice 7: add estimator hub + 5 sub-pages.** Largest content gap. |
| **API reference** | Hand-curated `website/reference.html` — does not enumerate the 35-symbol `__all__` export | **Slice 3 or 4: auto-generate from `__all__` + docstrings.** Path-dependent implementation. |
| **Executable examples** | No discipline anywhere | **Slice 5: tag-and-execute test.** Catches future drift. |
| **Link validation** | None | **Slice 6: link-check test.** Catches both README→site and intra-site drift. |
| **Contributor architecture docs** | `engine_architecture.md` + `adding_a_pass.md` + `validation_matrix.md` are GitHub-only — not discoverable from the site | **Slice 7: Architecture hub linking out to the canonical contributor entry points.** |
| **Pre-flight stop-and-report** | **Triggered (mild)** — site exists; do not propose blind replacement | Audit lays out the inventory + two paths; defers framework decision to user |
