"""T2-13: AIRR-spec conformance for D-gene inversion.

Background:
  When V(D)J recombination uses inversional joining (RAG cleavage at
  matching-RSS-orientation 12-bp spacers on D), the genomic DNA at the
  rearranged locus carries the REVERSE-COMPLEMENT of the canonical D
  allele on the sense strand. AID/SHM then acts on this rearranged
  template — the "wild-type" reference for any downstream mutation
  IS the inverted base, not the canonical D allele base.

  AIRR Rearrangement Schema requires `germline_alignment` to be column-
  wise consistent with `sequence_alignment`. Since the read carries
  inverted bases at inverted-D positions, the germline column must
  also carry the inverted base. Every mainstream AIRR-seq tool
  (IgBLAST, partis, MiXCR) does this.

What this file locks in:
  * d_inverted records: sequence == germline_alignment over the D
    region when there is no SHM, no corruption, no indels. Both
    reflect the rearranged template.
  * d_inverted records: germline_alignment over the D region equals
    revcomp(canonical_D[d_germline_start:d_germline_end]).
    The germline reference at inverted-D columns is the revcomp.
  * d_inverted + SHM: every mutation entry whose position falls in
    the D region shows the inverted base as the "from" side, not
    the canonical D base. Column-wise mutation calling is correct.
  * NOT-inverted records (d_inverted=False) MUST behave normally:
    germline_alignment over D equals the canonical D allele bases
    directly (no revcomp).

Pre-fix (per audit T2-13's literal recommendation), `n->germline`
would have stayed canonical — which would have made every inverted-D
column look like an SHM mutation against the canonical reference.
This test file would have caught that immediately.
"""
from __future__ import annotations

from GenAIRR import Experiment
from GenAIRR.ops import with_d_inversion, rate
from GenAIRR.protocol import _resolve_config


_COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    'N': 'N', 'n': 'n',
}


def _revcomp(seq: str) -> str:
    return "".join(_COMPLEMENT.get(b, 'N') for b in reversed(seq))


def _flatten_d_alleles(dc):
    out = {}
    if dc.d_alleles is None:
        return out
    for alleles in dc.d_alleles.values():
        for allele in alleles:
            out[allele.name] = allele.ungapped_seq
    return out


def _parse_mutations(mutation_str):
    """Parse 'pos:X>Y,pos:X>Y,...' into list of (pos, from_base, to_base)."""
    out = []
    if not mutation_str:
        return out
    for entry in mutation_str.split(","):
        entry = entry.strip()
        if not entry or ":" not in entry:
            continue
        pos_str, change = entry.split(":", 1)
        if ">" not in change:
            continue
        from_b, to_b = change.split(">", 1)
        try:
            out.append((int(pos_str), from_b.strip(), to_b.strip()))
        except ValueError:
            continue
    return out


N_SEQUENCES = 80
SEED = 4213


class TestSequenceMatchesGermlineNoSHM:
    """With inversion + no SHM, sequence and germline_alignment are
    identical over the D region (both encode the rearranged template)."""

    def test_d_region_columns_match(self):
        records = list(Experiment.on("human_igh")
                       .recombine(with_d_inversion(1.0))
                       .run(n=N_SEQUENCES, seed=SEED))

        inverted = [r for r in records if r["d_inverted"]]
        assert inverted, (
            "Expected at least one d_inverted record with prob=1.0")

        for r in inverted:
            d_start = r["d_sequence_start"]
            d_end = r["d_sequence_end"]
            if d_end <= d_start:
                continue
            seq_d = r["sequence"][d_start:d_end]
            germ_d = r["germline_alignment"][d_start:d_end]
            assert seq_d == germ_d, (
                f"d_inverted record but sequence != germline_alignment "
                f"over D region (no SHM):\n"
                f"  sequence  D: {seq_d!r}\n"
                f"  germline  D: {germ_d!r}\n"
                f"  d_call: {r.get('d_call')}")


class TestGermlineIsRevcompOfCanonical:
    """With inversion, germline_alignment over D equals revcomp of the
    canonical D allele (sliced at d_germline_start..d_germline_end).

    Note on coordinates: under inversion, d_germline_start..d_germline_end
    in the AIRR record refer to coordinates ON THE INVERTED template.
    The CANONICAL allele bases at the corresponding "post-trim" window,
    when read 3'→5' and complemented (i.e. reverse-complemented), give
    the inverted template bases — which is what the germline_alignment
    column should show. We verify that equivalence here.
    """

    def test_germline_d_matches_revcomp_canonical(self):
        dc = _resolve_config("human_igh")
        d_alleles = _flatten_d_alleles(dc)

        records = list(Experiment.on("human_igh")
                       .recombine(with_d_inversion(1.0))
                       .run(n=N_SEQUENCES, seed=SEED))

        inverted = [r for r in records if r["d_inverted"]]
        assert inverted

        checked = 0
        for r in inverted:
            d_call = r.get("d_call", "")
            # d_call may be an ambiguous comma-separated list; for the
            # invariant we only need the picked allele (d_call_true).
            d_call_true = r.get("d_call_true") or d_call.split(",")[0]
            if d_call_true not in d_alleles:
                continue
            allele_seq = d_alleles[d_call_true]

            d_seq_start = r["d_sequence_start"]
            d_seq_end = r["d_sequence_end"]
            d_g_start = r["d_germline_start"]
            d_g_end = r["d_germline_end"]
            if d_seq_end <= d_seq_start or d_g_end <= d_g_start:
                continue
            if d_g_end > len(allele_seq):
                continue

            # The germline_alignment D-region columns must equal the
            # revcomp of the canonical D bases at those germline
            # coordinates. Because we set both `n->germline` and
            # `n->current` to revcomp during inversion, sequence and
            # germline_alignment should both equal this revcomp.
            canonical_window = allele_seq[d_g_start:d_g_end]
            expected_revcomp = _revcomp(canonical_window)
            actual_germ = r["germline_alignment"][d_seq_start:d_seq_end]

            # Length match is required for a meaningful comparison; D
            # boundary ambiguity / trims may already be reflected in
            # both sides (d_sequence_* and d_germline_*).
            if len(actual_germ) != len(expected_revcomp):
                continue

            assert actual_germ == expected_revcomp, (
                f"germline_alignment D region != revcomp(canonical D):\n"
                f"  expected (revcomp): {expected_revcomp!r}\n"
                f"  got     (germ):     {actual_germ!r}\n"
                f"  d_call_true: {d_call_true}\n"
                f"  d_germline_start..end: {d_g_start}..{d_g_end}")
            checked += 1

        assert checked > 0, (
            "No d_inverted records had a resolvable d_call_true mapping "
            "to verify revcomp invariant — this test verified nothing.")


class TestMutationsUseInvertedBase:
    """With inversion + SHM, every mutation entry whose position lies
    inside the D region shows the inverted (revcomp) base on the FROM
    side of the 'pos:X>Y' annotation — never the canonical D base.

    This is the column-wise consistency invariant: the germline column
    at position `pos` is `germline_alignment[pos]`, which equals the
    inverted base after D inversion. A mutation entry must be
    `germline_alignment[pos] > sequence[pos]`.
    """

    def test_mutation_from_base_matches_germline_alignment(self):
        records = list(Experiment.on("human_igh")
                       .recombine(with_d_inversion(1.0))
                       .mutate(rate(0.05, 0.15))
                       .run(n=N_SEQUENCES, seed=SEED))

        inverted_with_d_mutations = 0
        for r in records:
            if not r["d_inverted"]:
                continue
            d_start = r["d_sequence_start"]
            d_end = r["d_sequence_end"]
            germ = r["germline_alignment"]
            seq = r["sequence"]

            for pos, from_b, to_b in _parse_mutations(r.get("mutations", "")):
                if not (d_start <= pos < d_end):
                    continue  # not in D region
                if pos >= len(germ) or pos >= len(seq):
                    continue
                # Column-wise consistency.
                assert from_b == germ[pos], (
                    f"Mutation 'from' base disagrees with germline_alignment "
                    f"column at inverted-D position {pos}:\n"
                    f"  mutation says: {from_b}>{to_b}\n"
                    f"  germline col:  {germ[pos]}\n"
                    f"  sequence col:  {seq[pos]}\n"
                    f"  This is the T2-13 audit-fix bug: germline must "
                    f"track the rearranged (inverted) template.")
                assert to_b == seq[pos], (
                    f"Mutation 'to' base disagrees with sequence column "
                    f"at position {pos}: mutation says >{to_b}, "
                    f"sequence has {seq[pos]}")
                inverted_with_d_mutations += 1

        assert inverted_with_d_mutations > 0, (
            "No mutations landed inside any inverted-D region across "
            f"{N_SEQUENCES} sequences with rate(0.05, 0.15) — test did "
            "not exercise the invariant. Bump rate or N if this becomes "
            "flaky.")


class TestNonInvertedRecordsUnchanged:
    """The fix must not affect records that did NOT undergo inversion.
    For those, germline_alignment over D continues to equal the
    canonical D allele bases directly (no revcomp)."""

    def test_non_inverted_germline_matches_canonical(self):
        dc = _resolve_config("human_igh")
        d_alleles = _flatten_d_alleles(dc)

        # Use prob=0.0 so we know none can be inverted; this also
        # acts as a regression test that the per-position flag does
        # not leak when the step doesn't fire.
        records = list(Experiment.on("human_igh")
                       .recombine(with_d_inversion(0.0))
                       .run(n=N_SEQUENCES, seed=SEED))

        assert all(not r["d_inverted"] for r in records)

        checked = 0
        for r in records:
            d_call_true = r.get("d_call_true") or r.get("d_call", "").split(",")[0]
            if d_call_true not in d_alleles:
                continue
            allele_seq = d_alleles[d_call_true]
            d_seq_start = r["d_sequence_start"]
            d_seq_end = r["d_sequence_end"]
            d_g_start = r["d_germline_start"]
            d_g_end = r["d_germline_end"]
            if d_seq_end <= d_seq_start or d_g_end <= d_g_start:
                continue
            if d_g_end > len(allele_seq):
                continue

            canonical_window = allele_seq[d_g_start:d_g_end]
            actual_germ = r["germline_alignment"][d_seq_start:d_seq_end]
            if len(actual_germ) != len(canonical_window):
                continue

            assert actual_germ == canonical_window, (
                f"Non-inverted record but germline_alignment D region "
                f"!= canonical D bases:\n"
                f"  expected (canonical): {canonical_window!r}\n"
                f"  got      (germ):      {actual_germ!r}\n"
                f"  d_call_true: {d_call_true}")
            checked += 1

        assert checked > 0, "No non-inverted records were verifiable"
