"""The catalog must describe the engine that exists, not the one it remembers.

A catalog of biological mechanisms is only worth having if it cannot quietly
fall out of step with the simulation. These tests are what stop it becoming a
comfortable fiction: every pass the engine defines is accounted for exactly
once, every truth field named is a field a record actually carries, and a
missing registry fails loudly instead of skipping.
"""
from __future__ import annotations

import pytest

from GenAIRR._compiled import pass_registry
from tests.biology import catalog

pytestmark = pytest.mark.pin


def test_every_engine_pass_is_classified_exactly_once():
    """Biology or infrastructure, once each.

    Not every pass is biological: the engine documents EchoPass as its
    deterministic transform reference and SampleBasePass as its sampling
    reference. Forcing those into a biological mechanism would make the catalog
    assert something false, so infrastructure is a classification of its own
    rather than a gap.
    """
    claimed = catalog.classified_passes()
    duplicates = sorted({p for p in claimed if list(claimed).count(p) > 1})
    assert duplicates == [], f"passes classified more than once: {duplicates}"
    registry = set(pass_registry())
    assert set(claimed) == registry, (
        f"unclassified passes: {sorted(registry - set(claimed))}\n"
        f"classified but not in the engine: {sorted(set(claimed) - registry)}\n"
        "A new pass must be classified -- as a mechanism, or as infrastructure "
        "with a reason -- before this passes.")


def test_an_unavailable_registry_fails_rather_than_skips(monkeypatch):
    """Fail closed. A check that skips when it cannot see the engine reports the
    catalog complete on no evidence, which is worse than no check."""
    def boom():
        raise RuntimeError("registry unavailable")
    monkeypatch.setattr(catalog, "pass_registry", boom)
    with pytest.raises(RuntimeError):
        catalog.classified_passes()


def test_mechanism_ids_are_unique():
    ids = [m.id for m in catalog.all_mechanisms()]
    assert sorted(ids) == sorted(set(ids))


def test_every_mechanism_id_is_prefixed_with_its_stage():
    wrong = [m.id for m in catalog.all_mechanisms() if not m.id.startswith(f"{m.stage}.")]
    assert wrong == [], f"mechanism ids not prefixed with their stage: {wrong}"


def test_stages_are_exactly_the_declared_order():
    """Developmental order, then the two that are not developmental stages:
    assay is what the protocol did to the molecule afterwards, and projection is
    the surface a benchmark reads rather than a step anything performs."""
    assert catalog.STAGES == ("germline", "recombination", "junction", "expression",
                              "maturation", "clonality", "assay", "projection")
    assert {m.stage for m in catalog.all_mechanisms()} <= set(catalog.STAGES)


def test_every_mechanism_states_a_summary_and_owns_truth_fields():
    thin = [m.id for m in catalog.all_mechanisms() if not m.summary.strip() or not m.truth_fields]
    assert thin == [], f"mechanisms with no summary or no truth fields: {thin}"


def test_every_truth_field_is_a_field_some_mode_actually_produces():
    """Measured across every run mode, not checked against a hand-written list.

    The earlier version exempted a list of fields "not in a plain record", and
    that list named `generation` -- which exists nowhere in GenAIRR; the field
    is `lineage_generation`. An invented name in the one place the check stops
    looking is the exact failure this test exists to prevent, so the exemption
    is gone and the known set is taken from real records.
    """
    from tests.support.fixtures.record_modes import every_known_field
    known = every_known_field()
    unknown = sorted({f for m in catalog.all_mechanisms() for f in m.truth_fields
                      if f not in known})
    assert unknown == [], f"truth fields no run mode produces: {unknown}"


def test_the_lineage_truth_fields_have_owners():
    """These are what a clonal-inference tool is scored against, and a single
    clonality mechanism owning `clone_id` left all five unowned."""
    from tests.support.fixtures.record_modes import fields_by_mode
    lineage_fields = {f for f in fields_by_mode()["clonal_lineage"] if f.startswith("lineage_")}
    owned = {f for m in catalog.all_mechanisms() for f in m.truth_fields}
    assert lineage_fields <= owned, f"unowned lineage truth: {sorted(lineage_fields - owned)}"


def test_infrastructure_passes_each_state_why_they_are_not_biology():
    unexplained = [p for p, reason in catalog.INFRASTRUCTURE.items() if not reason.strip()]
    assert unexplained == [], f"infrastructure passes with no stated reason: {unexplained}"


def test_by_id_finds_a_mechanism_and_rejects_an_unknown_one():
    assert catalog.by_id("junction.n_addition").stage == "junction"
    with pytest.raises(KeyError):
        catalog.by_id("junction.n_additon")


# ── the projection surface, and truth that is not exported ───────────────────

def test_the_alignment_surface_has_an_owner():
    """These are the fields an aligner benchmark actually compares.

    They were deliberately unowned on the reasoning that no single mechanism
    produces them -- a trim, an N addition, a mutation and an indel all land in
    the same CIGAR. True, but the consequence was that the catalog could read
    fully green while every coordinate a benchmark consumes was wrong. The
    answer is a cross-cutting mechanism with its own claims, not an exemption.
    """
    owned = {f for m in catalog.all_mechanisms() for f in m.truth_fields}
    surface = {
        "v_cigar", "d_cigar", "j_cigar",
        "v_sequence_start", "v_sequence_end", "d_sequence_start", "d_sequence_end",
        "j_sequence_start", "j_sequence_end",
        "v_alignment_start", "v_alignment_end", "d_alignment_start", "d_alignment_end",
        "j_alignment_start", "j_alignment_end",
        "sequence_alignment", "germline_alignment", "germline_alignment_d_mask",
        "sequence_aa",
    }
    assert surface <= owned, f"alignment-surface fields owned by nobody: {sorted(surface - owned)}"


def test_no_exported_field_is_owned_by_nobody():
    """Every field a record carries is either claimed by a mechanism or named
    in the explicit unclaimed record, with a reason. A field nobody owns is a
    field nobody is answerable for."""
    from tests.support.fixtures.record_modes import every_known_field
    owned = {f for m in catalog.all_mechanisms() for f in m.truth_fields}
    unaccounted = sorted(every_known_field() - owned - set(catalog.UNCLAIMED_FIELDS))
    assert unaccounted == [], (
        f"exported fields with no owner and no recorded reason: {unaccounted}")


def test_every_unclaimed_field_states_why():
    assert catalog.UNCLAIMED_FIELDS
    thin = [f for f, reason in catalog.UNCLAIMED_FIELDS.items() if not reason.strip()]
    assert thin == [], f"unclaimed fields with no reason: {thin}"


def test_a_mechanism_whose_truth_is_not_exported_says_so():
    """`germline.genotype` owned only `locus`, and `assay.ambiguous_calls` and
    `assay.duplex_formation` owned only `sequence`, because their real truth is
    not in the record at all. Reporting those as merely untested confuses "no
    one has written the test" with "there is nothing to test against"."""
    gapped = [m for m in catalog.all_mechanisms() if m.truth_gap]
    assert gapped, "no mechanism declares a truth gap, which cannot be right"
    for mechanism in gapped:
        assert mechanism.truth_gap.strip(), mechanism.id


def test_a_truth_gap_reports_as_its_own_status_not_as_unverified():
    rows = {r.mechanism_id: r for r in catalog.coverage_rows({})}
    gapped = [m.id for m in catalog.all_mechanisms() if m.truth_gap]
    for mechanism_id in gapped:
        assert rows[mechanism_id].status == "no-truth-exported", (
            f"{mechanism_id} declares a truth gap but reports "
            f"{rows[mechanism_id].status!r}")


def test_a_truth_gap_stays_visible_even_when_the_mechanism_verifies():
    """A mechanism can be verified on the one borrowed field it owns while the
    truth that actually matters is still not exported. assay.ambiguous_calls is
    exactly that: the claim "an N is not a mutation" can be established against
    the sequence, but no N-count field exists, so "verified" alone would hide
    the gap rather than report it."""
    gapped = {m.id for m in catalog.all_mechanisms() if m.truth_gap}
    rows = {r.mechanism_id: r for r in catalog.coverage_rows(
        {mid: ("x::y",) for mid in gapped})}
    for mechanism_id in gapped:
        assert rows[mechanism_id].truth_gap, (
            f"{mechanism_id} lost its truth gap once a test was registered")


def test_a_field_the_catalog_calls_unclaimed_does_not_count_toward_coverage():
    """Thirteen fields are in both UNCLAIMED_FIELDS and some mechanism's truth
    fields, and the entries say why: a configuration echo, a harness label, a
    derived ratio, or an inference output of a search this tool does not run.
    Each is listed beside the mechanism it belongs to but is not something that
    mechanism can be held to as truth.

    Counting them in the denominator would make their owners permanently
    partial however thoroughly either was tested, so the denominator is the
    verifiable fields rather than every field listed.

    Two fields were wrongly in this set and are not any more. Both were the
    only reason their mechanism read verified, which is exactly the hazard the
    pin below exists to catch:

    - productivity_reasons is a per-record outcome -- which assessment failed
      -- and the suite asserts specific codes in 22 places. Calling it
      diagnostic text let expression.productivity read verified while a
      different tier silently covered the field.
    - insert_size is drawn from a distribution per record, and r2 is positioned
      by it. It looked like an echo because the fixture pinned it to a constant.
    """
    overlapping = {f for m in catalog.all_mechanisms() for f in m.truth_fields
                   if f in catalog.UNCLAIMED_FIELDS}
    assert overlapping, "the overlap is the thing under test; it cannot be empty"
    for mechanism in catalog.all_mechanisms():
        excluded = set(mechanism.truth_fields) & set(catalog.UNCLAIMED_FIELDS)
        assert not (excluded & set(mechanism.verifiable_fields)), (
            f"{mechanism.id} counts {sorted(excluded & set(mechanism.verifiable_fields))} "
            "toward its coverage, but the catalog says no mechanism is "
            "answerable for them")


def test_the_fields_excluded_from_coverage_are_pinned():
    """Shrinking the denominator is a way to manufacture green, so the exact
    set of owned-but-excluded fields is pinned rather than derived.

    Adding a field to UNCLAIMED_FIELDS raises its owner's coverage ratio
    without anyone testing anything. Each entry is a judgement -- that nothing
    a test could do would establish a configuration echo, a harness label or a
    derived ratio -- and a new one has to be argued for here, not slipped in.
    """
    excluded = {f for m in catalog.all_mechanisms() for f in m.truth_fields
                if f in catalog.UNCLAIMED_FIELDS}
    assert excluded == {
        "sequence_id",             # clonality.founder
        "locus",                   # germline.genotype
        "productivity_policy",     # expression.productivity
        "productivity_scope",      # expression.productivity
        "read_layout",             # assay.paired_end
        "mutation_rate",           # maturation.hypermutation
        "c_call",                  # maturation.class_switching
        # projection.alignment_scores. The one exclusion that is not a
        # modelling decision: a score and a support are what an alignment
        # SEARCH produces, and GenAIRR assembles from a known allele. The
        # columns are exported as null because the AIRR schema has them, which
        # test_the_score_columns_are_still_exported pins, and
        # test_score_and_support_are_null_on_every_record pins the premise --
        # if the engine ever populates one, the exclusion has to be argued
        # again rather than quietly outlived.
        "v_score", "d_score", "j_score",
        "v_support", "d_support", "j_support",
    }, (
        "the set of fields excluded from coverage changed. Each exclusion "
        "raises its mechanism's ratio without any test being written, so a "
        "change here needs its reason in UNCLAIMED_FIELDS and this pin updated "
        "deliberately.")


def test_every_mechanism_can_reach_verified_or_declares_why_not():
    """A mechanism whose every field is a configuration echo can never be
    verified however well it is tested, which must be stated rather than left
    to be inferred from a stuck report."""
    stuck = [m.id for m in catalog.all_mechanisms()
             if not m.verifiable_fields and not m.truth_gap]
    assert stuck == [], (
        f"mechanisms with no verifiable field and no declared truth gap: {stuck}")


def test_a_truth_field_has_one_owner_unless_the_sharing_is_declared():
    """Two mechanisms quietly answerable for the same field means neither is.

    `sequence` and `original_v_call` are shared on purpose and say so. Any
    other overlap is an accident, and an accident here makes a field look
    covered because a different mechanism's test touched it.
    """
    from collections import Counter
    counted = Counter(f for m in catalog.all_mechanisms() for f in m.truth_fields)
    shared = {f for f, n in counted.items() if n > 1}
    assert shared == set(catalog.SHARED_FIELDS), (
        f"undeclared shared fields: {sorted(shared - set(catalog.SHARED_FIELDS))}\n"
        f"declared but no longer shared: {sorted(set(catalog.SHARED_FIELDS) - shared)}")
    thin = [f for f, reason in catalog.SHARED_FIELDS.items() if not reason.strip()]
    assert thin == [], f"shared fields with no reason: {thin}"


def test_a_mechanism_with_a_declared_truth_gap_cannot_read_verified():
    """The status column is the thing that gets quoted, so a declared gap must
    be visible there and not only in a neighbouring field.

    assay.ambiguous_calls is the case: the claim "an N is not a mutation" can
    be established against the sequence, but no N-count field exists. Reporting
    that as `verified` says the mechanism is covered when a third of what it
    means is not exported at all.
    """
    gapped = {m.id for m in catalog.all_mechanisms() if m.truth_gap}
    covered = {m.id: m.verifiable_fields for m in catalog.all_mechanisms()}
    rows = {r.mechanism_id: r for r in catalog.coverage_rows(
        {mid: ("x::y",) for mid in gapped},
        covered_by_mechanism={mid: covered[mid] for mid in gapped})}
    for mechanism_id in gapped:
        assert rows[mechanism_id].status != "verified", (
            f"{mechanism_id} declares a truth gap but its status reads verified")
        assert rows[mechanism_id].truth_gap, f"{mechanism_id} lost its truth gap"


def test_every_claimed_field_is_read_by_the_test_that_claims_it():
    """A registration must claim only fields its test actually inspects.

    This is the failure the whole scheme exists to prevent, and it happened
    anyway: seven registrations claimed eleven fields their bodies never read,
    including end_loss_3_length (zero in every record of its own fixture) and
    the six germline-coordinate record fields, whose values come from the
    column walker rather than the IR track the test inspects. Each produced a
    green or partial row on nothing.

    The decorator is excluded from the searched text, since it names every
    claimed field by construction -- which is why an earlier version of this
    analysis reported no offenders.
    """
    import ast
    import re

    from tests.support import paths

    offenders = {}
    for module in sorted(paths.tests_root().rglob("test_*.py")):
        if "__pycache__" in module.parts:
            continue
        tree = ast.parse(module.read_text(encoding="utf-8"))
        helpers = {node.name: "\n".join(ast.unparse(s) for s in node.body)
                   for node in tree.body
                   if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if not node.name.startswith("test_"):
                continue            # a nested sample used to prove a rejection
            for decorator in node.decorator_list:
                if not (isinstance(decorator, ast.Call)
                        and getattr(decorator.func, "id", "") == "verifies"):
                    continue
                claimed = next((ast.literal_eval(k.value) for k in decorator.keywords
                                if k.arg == "fields"), ())
                body = "\n".join(ast.unparse(s) for s in node.body)
                for name, text in helpers.items():
                    if name != node.name and re.search(rf"\b{name}\b", body):
                        body += "\n" + text
                missing = [f for f in claimed if f not in body]
                if missing:
                    offenders[f"{module.name}::{node.name}"] = missing
    assert offenders == {}, (
        "registrations claiming fields their test never reads:\n"
        + "\n".join(f"  {where}: {fields}" for where, fields in sorted(offenders.items()))
        + "\n\nClaim only what the body inspects. Read fields by name rather than "
          "through an f-string, so this check can see them.")


def _has_non_establishing_marker(node) -> bool:
    """Whether a test carries xfail, skip or skipif, by any spelling.

    Matched on the attribute path's last two names, so ``pytest.mark.xfail``,
    ``mark.xfail`` and a bare decorator or a call all read the same.
    """
    import ast

    from tests.biology.catalog._plugin import NON_ESTABLISHING

    for decorator in node.decorator_list:
        target = decorator.func if isinstance(decorator, ast.Call) else decorator
        names = []
        while isinstance(target, ast.Attribute):
            names.append(target.attr)
            target = target.value
        if isinstance(target, ast.Name):
            names.append(target.id)
        if any(name in NON_ESTABLISHING for name in names):
            return True
    return False


def _rank(basis: str) -> int:
    from tests.biology.catalog._registry import VERIFICATION_BASES
    return VERIFICATION_BASES.index(basis)


def _declared_registrations():
    """Every ``@verifies`` in the suite, read from the source.

    Static rather than live: the plugin's numbers depend on which tests ran,
    and a pin that moves with the selection pins nothing. What this measures is
    the coverage the tree *declares*, which is an upper bound on the report --
    the report additionally requires each registered test to have been observed
    passing, so it can be equal or worse, never better.
    """
    import ast
    from collections import defaultdict

    from tests.support import paths

    tests, covered = defaultdict(list), defaultdict(set)
    bases = defaultdict(dict)
    for module in sorted(paths.tests_root().rglob("test_*.py")):
        if "__pycache__" in module.parts:
            continue
        tree = ast.parse(module.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if not node.name.startswith("test_"):
                # A nested sample used to prove a rejection -- test_registration
                # registers a deliberately typo'd id inside pytest.raises. It is
                # not coverage, and counting it would make the typo the scan is
                # supposed to catch indistinguishable from a real one.
                continue
            # A registration on an xfailed or skipped test does not establish
            # its fields, and the plugin that prints the report discounts it
            # (catalog/_plugin.NON_ESTABLISHING). Ignoring markers here made
            # this scan disagree with the artifact it pins: the scan said 16
            # verified while the report printed 15.
            if _has_non_establishing_marker(node):
                continue
            for decorator in node.decorator_list:
                if not (isinstance(decorator, ast.Call)
                        and getattr(decorator.func, "id", "") == "verifies"):
                    continue
                mechanism_id = ast.literal_eval(decorator.args[0])
                fields = next((ast.literal_eval(k.value) for k in decorator.keywords
                               if k.arg == "fields"), ())
                basis = next((ast.literal_eval(k.value) for k in decorator.keywords
                              if k.arg == "basis"), None)
                assert basis is not None, (
                    f"{module.name}::{node.name} declares no basis; the kwarg "
                    "is required, so this means the scan is reading a call "
                    "shape it does not understand")
                tests[mechanism_id].append(f"{module.name}::{node.name}")
                covered[mechanism_id].update(fields)
                strengths = bases[mechanism_id]
                for one in fields:
                    current = strengths.get(one)
                    if current is None or _rank(basis) < _rank(current):
                        strengths[one] = basis
    return ({k: tuple(v) for k, v in tests.items()},
            {k: tuple(v) for k, v in covered.items()},
            {k: dict(v) for k, v in bases.items()})


def test_the_coverage_report_numbers_are_pinned():
    """The report is the artifact anyone quotes, and nothing pinned it.

    Without this, coverage can fall silently: delete a registration, or widen a
    mechanism's truth_fields, and the suite stays green while the report says
    something weaker. Both are legitimate changes -- so this is a pin, not a
    floor, and moving it is how you record that the number moved on purpose.

    It last moved when registrations began declaring HOW they establish their
    claim. The distribution was 16 verified, 1 gap-limited, 2 partial, 9
    unverified, 3 no-truth. Six of those sixteen were reached only by
    comparisons between two exported views of the same number, or by
    assertions a wrong value would also satisfy, and are now reported
    ``consistent``; the previously gap-limited row joined them for the same
    reason. Nothing about the engine changed when this number fell -- what
    changed is that the report stopped counting the weaker basis as proof.

    It then rose by two: junction.assembly and junction.coding_end_processing
    went from no registered tests at all to verified on a reference basis --
    against the anchor coordinate and the conserved residue the reference
    declares, rather than against another exported field.

    The clonality block then moved four rows at once, without a single new
    engine field. Two gained a strong basis and a declared gap
    (clonality.founder on clone_id, clonality.sampling on both its fields);
    two were reclassified to no-truth-exported, because measurement showed
    there is nothing in the record to test -- lineage_parent_id names no
    exported record, lineage_node_id is unique only within a clone, and
    lineage_affinity is the 0.0 sentinel at the default selection strength.
    "unverified" fell from 7 to 3 and "no truth exported" rose from 3 to 5,
    which is the honest reading: those two were never waiting for a test.

    germline.allele_choice then went from no coverage to verified on all
    seven of its fields, which is the one that matters most: those calls are
    what an allele caller is scored against. The check that was missing is
    that the recorded bases come from the allele the call NAMES -- the
    engine's own comparison recomputes the call from the same assignment and
    agrees with itself either way, and the only external check on v_call
    compared the named allele's LENGTH, which cannot separate 195 of the 198
    IGH V alleles.

    Then every remaining `consistent` row was given a basis that establishes
    its values, which is what that status existed to make visible. Six
    mechanisms moved and the count is now zero: paired_end reads its reads
    against the molecule rather than against their own length, contamination
    against the divergence of the observed bases from the germline the record
    names, ambiguous_calls against the same seed without the pass,
    segment_coordinates against the span its bases occupy, receptor_revision
    against an unrevised run of the same seed, and p_nucleotide against the
    tip palindrome at all four coding ends rather than one.

    assay.polymerase_indels joined them from unverified: n_indels is the net
    length the molecule changed by against the same seed without the pass,
    120 of 120 in each direction, so it counts applied indels and not drawn
    ones -- worth establishing because the two error counters beside it did
    report the drawn number. Its ceiling is gap-limited because the three
    per-segment counters do not partition it: an indel in an NP region
    belongs to none of them, on 18 of 200 records.

    That leaves one unverified mechanism, assay.duplex_formation, an
    experimental surface.

    Both partial rows then closed. junction.n_addition reached its gapped
    ceiling on np1_aa and np2_aa, the whole codons of the NP interval in the
    record's own frame. maturation.hypermutation reached verified on all
    eighteen: the four per-segment and six per-subregion EVENT counters are
    now compared against what the read shows, bucketed by segment and by the
    allele's IMGT table, which is an equality at one mutation per record and
    an inequality in a known direction above it. The partition those counters
    had before is satisfied by a mutation routed to the wrong bucket; this is
    not.

    So every mechanism is now at verified, at a declared gapped ceiling, or
    carrying no exported truth at all -- except assay.duplex_formation.

    Then an audit of the basis column against the bodies that declare it found
    six mechanisms claiming a basis stronger than their test supported on at
    least one field, and the headline fell from 20 to 19. Three were
    strengthened rather than relabelled, because the stronger check was
    available and nobody had written it: assay.paired_end now reads its
    windows against the declared layout and its reads against the engine's own
    per-base record, junction.assembly derives the junction window from the
    reference anchor codons instead of slicing at the bounds the record
    reports, and projection.germline_alignment compares the alignment to the
    recorded molecule and the D mask to the span the provenance assigns to D.
    Two were label corrections between two strong bases, so no status moved.
    One, recombination.assembly, is now reported `consistent`: sequence_length
    against len(sequence) is two views of one exported string and there is no
    stronger comparison available for it, which is the status existing to be
    visible rather than a regression in the engine.

    The last two rows then closed, and no status is left that means "nobody
    looked". projection.alignment_scores reached verified on the three
    identity fields -- each segment's fraction of matching templated positions,
    counted off the provenance, exact across mutation, reverse complement and
    an inverted D -- once score and support were recorded as unclaimable
    rather than as missing truth: they are an alignment search's outputs, and
    this tool assembles from a known allele. assay.duplex_formation reached
    verified against repair scripts authored in the test, including the one
    measurement that justifies its whole plane -- a base resected and
    refilled produces a molecule byte-identical to a blunt ligation's, so the
    read cannot tell them apart and duplex_repair is the only record that a
    base was cut and resynthesised.

    So the surface is now: 21 verified, 5 at a declared gapped ceiling, 4 with
    no exported truth at all, and 1 consistent. `unverified` is empty.
    """
    tests, covered, bases = _declared_registrations()
    rows = catalog.coverage_rows(tests, covered_by_mechanism=covered,
                                 basis_by_mechanism=bases)
    assert len(rows) == 31, f"the catalog has {len(rows)} mechanisms, not 31"
    # Pinned per mechanism, not as a distribution. The distribution was the
    # original pin and it cannot see a swap: one mechanism falling out of
    # verified while another rises into it leaves every count untouched, and
    # the row that fell is exactly the thing a pin is for. Six of these rows
    # moved individually over this branch's work and the headline count moved
    # three times, which is how the blind spot was noticed.
    status = {row.mechanism_id: row.status for row in rows}
    expected = {
        "germline.allele_choice":         "verified",
        "germline.genotype":              "no-truth-exported",
        "recombination.assembly":         "consistent",
        "recombination.d_inversion":      "verified",
        "recombination.trimming":         "verified",
        "junction.assembly":              "verified",
        "junction.coding_end_processing": "verified",
        "junction.n_addition":            "gap-limited",
        "junction.p_nucleotide":          "verified",
        "expression.productivity":        "verified",
        "expression.receptor_revision":   "verified",
        "maturation.class_switching":     "no-truth-exported",
        "maturation.hypermutation":       "verified",
        "clonality.founder":              "gap-limited",
        "clonality.sampling":             "gap-limited",
        "clonality.selection":            "no-truth-exported",
        "clonality.tree_topology":        "no-truth-exported",
        "assay.ambiguous_calls":          "gap-limited",
        "assay.contamination":            "verified",
        "assay.duplex_formation":         "verified",
        "assay.end_loss":                 "verified",
        "assay.orientation":              "verified",
        "assay.paired_end":               "verified",
        "assay.pcr_error":                "verified",
        "assay.polymerase_indels":        "gap-limited",
        "assay.quality_error":            "verified",
        "projection.alignment_scores":    "verified",
        "projection.cigar":               "verified",
        "projection.germline_alignment":  "verified",
        "projection.segment_coordinates": "verified",
        "projection.translation":         "verified",
    }
    assert status == expected, (
        "the declared coverage changed, per mechanism:\n"
        + "\n".join(f"    {mechanism}: {was} -> {status.get(mechanism, 'absent')}"
                     for mechanism, was in expected.items()
                     if status.get(mechanism) != was)
        + "\n".join(f"    {mechanism}: not pinned, reports {value}"
                     for mechanism, value in status.items()
                     if mechanism not in expected)
        + "\nIf that was deliberate, update this pin in the same commit and "
        "say why in the message. If it was not, a registration or a truth "
        "field moved without anyone noticing.")
    counts = {name: sum(1 for value in status.values() if value == name)
              for name in sorted(set(status.values()))}
    assert counts == {"consistent": 1, "gap-limited": 5, "no-truth-exported": 4,
                      "verified": 21}, counts


def test_the_pinned_report_is_not_a_fiction_about_which_tests_exist():
    """The pin above is computed from the source, so it would survive every
    registered test being deleted from the run. This ties it to reality."""
    tests, _, _ = _declared_registrations()
    found = sum(len(v) for v in tests.values())
    assert found >= 20, (
        f"only {found} registrations found; the scan is not reading the suite. "
        "This is a floor on the scanner, not a coverage target -- the "
        "distribution is pinned above.")
    unknown = sorted(set(tests) - {m.id for m in catalog.all_mechanisms()})
    assert unknown == [], f"registrations naming no known mechanism: {unknown}"


def test_the_published_report_is_collected_over_the_whole_suite():
    """Coverage is claimed by registration, not by location -- so the command
    that publishes the report has to collect everywhere registrations live.

    The first version of the CI step ran `pytest tests/biology
    --biology-coverage` and published 15 verified where the suite establishes
    16, because one mechanism's registration lives under components/. That is
    the one thing "registration, not location" exists to allow, under-counted
    by the artifact meant to report it.
    """
    from tests.support import paths

    workflow = (paths.repo_root() / ".github" / "workflows" / "test.yml")
    assert workflow.exists(), "the test workflow is missing"
    lines = [line.strip() for line in workflow.read_text(encoding="utf-8").splitlines()
             if "--biology-coverage" in line]
    assert lines, "no CI step publishes the coverage report"
    for line in lines:
        assert " tests/ " in line or " tests " in line, (
            f"the report is collected from a subdirectory, so it under-counts "
            f"any mechanism registered outside it: {line!r}")


def test_registrations_live_outside_the_biology_directory():
    """The guard above only matters while that is true, and it is the design.

    If every registration moved under biology/, the location-independence
    would be untested and the guard would be protecting nothing.
    """
    from tests.support import paths

    tests, _, _ = _declared_registrations()
    outside = {module for module in
               (m.name for m in sorted(paths.tests_root().rglob("test_*.py"))
                if "__pycache__" not in m.parts and "biology" not in m.parts)
               if any(module in entry for entries in tests.values() for entry in entries)}
    assert outside, (
        "every registration lives under biology/, so coverage-by-registration "
        "is not actually exercised and the workflow guard protects nothing")


def _loop_products(function):
    """Names a loop in ``function`` writes to: counters, lists, sets."""
    import ast
    produced = set()
    loops = [node for node in ast.walk(function)
             if isinstance(node, (ast.For, ast.While))]
    for loop in loops:
        for node in ast.walk(loop):
            if isinstance(node, ast.AugAssign) and isinstance(node.target, ast.Name):
                produced.add(node.target.id)
            elif isinstance(node, ast.Assign):
                produced.update(target.id for target in node.targets
                                if isinstance(target, ast.Name))
            elif (isinstance(node, ast.Call)
                  and isinstance(node.func, ast.Attribute)
                  and isinstance(node.func.value, ast.Name)
                  and node.func.attr in {"append", "add", "update", "extend"}):
                produced.add(node.func.value.id)
    return produced


def test_a_registered_test_that_can_skip_its_assertions_proves_it_did_not():
    """The vacuous-pass guard.

    A registered test that iterates records and `continue`s past the
    uninteresting ones can compare nothing at all and still report green. That
    is not hypothetical: the NP-translation test keyed both spans off the D
    coordinates, so on a VJ locus every record was skipped and it passed over
    zero comparisons. The oracle probe that bucketed D and J coordinates
    against the V subregion table reported 3 of 60 agreement for the same
    reason -- it was comparing the wrong thing, and nothing said how much.

    So: a registered test containing a `continue` must, after its loops,
    assert something that refers to a name the loop writes to -- a counter it
    incremented, or a collection it filled. "I compared 40 of them" is such an
    assertion; "the experiment exists" is not.

    Two limits, stated rather than papered over. It cannot tell a meaningful
    floor from a trivial one: `assert checked >= 0` would satisfy it. And it
    only fires on `continue`, so a test that skips work another way -- an `if`
    wrapping its assertions -- is not caught. It targets the shape that has
    actually failed here, twice.
    """
    import ast

    from tests.support import paths

    offenders = []
    for module in sorted(paths.tests_root().rglob("test_*.py")):
        source = module.read_text(encoding="utf-8")
        if "@verifies(" not in source:
            continue
        for function in ast.walk(ast.parse(source)):
            if not isinstance(function, ast.FunctionDef):
                continue
            if not any(isinstance(decorator, ast.Call)
                       and getattr(decorator.func, "id", "") == "verifies"
                       for decorator in function.decorator_list):
                continue
            if not any(isinstance(node, ast.Continue)
                       for node in ast.walk(function)):
                continue
            loops = [node for node in ast.walk(function)
                     if isinstance(node, (ast.For, ast.While))]
            inside = {id(node) for loop in loops for node in ast.walk(loop)
                      if isinstance(node, ast.Assert)}
            after = [node for node in ast.walk(function)
                     if isinstance(node, ast.Assert) and id(node) not in inside]
            produced = _loop_products(function)
            speaks = any(isinstance(name, ast.Name) and name.id in produced
                         for node in after for name in ast.walk(node.test))
            if not speaks:
                offenders.append(f"{module.name}::{function.name}")
    assert offenders == [], (
        "these registered tests can skip every comparison and still pass, "
        "because nothing after their loops refers to what the loop produced:\n"
        + "\n".join(f"  {name}" for name in offenders)
        + "\n\nAdd an assertion on the count compared, or on the collection "
          "filled, so a run that compared nothing fails.")
