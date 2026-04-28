"""T2-12: ``with_ns`` no longer flips ``productive`` when the user
asked for ``productive=True``.

Background — the bug:
  AIRR ``productive`` is recomputed from the FINAL post-corruption
  sequence (``airr.c::derive_final_productivity``). If
  ``step_insert_ns`` lands an N on the V or J anchor codon, the codon
  retranslates to ``?``, the conserved-residue rule fails, and
  ``productive`` silently flips to ``false`` — even though the
  rearrangement that came out of the assembly stage was productive.

The fix (Option B from the design discussion):
  When the user requested ``PRODUCTIVITY_PRODUCTIVE_ONLY``,
  ``step_insert_ns`` protects the V and J anchor codons. Other modes
  (``MIXED``, ``NON_PRODUCTIVE_ONLY``) keep the full noise model so
  the realistic AIRR-seq behavior of "Ns can land anywhere" is
  preserved when the user has no productive contract to honor.

What this file tests:
  * ``productive=True`` + heavy ``with_ns`` → essentially zero flip
    rate (the protection works).
  * ``productive=True`` + no corruption → 100% productive (sanity:
    we didn't accidentally degrade the productive contract elsewhere).
  * ``productive=True`` + heavy ``with_ns`` still inserts Ns at NON-
    anchor positions (we're not over-protecting).
"""
from __future__ import annotations

from GenAIRR import Experiment
from GenAIRR.ops import with_ns


# Use a fairly aggressive N rate. At 5% per base over a ~350-nt sequence,
# expectation is ~17 N's per sequence — enough to almost-certainly hit
# an anchor codon by chance pre-fix, but per the new logic the anchor
# codons are guarded entirely.
N_PROB = 0.05
N_SEQUENCES = 200
SEED = 42


class TestProductiveContractHonored:
    """productive=True must mean productive=True even when N corruption
    is enabled."""

    def test_no_flips_under_productive_with_ns(self):
        records = list(Experiment.on("human_igh")
                       .observe(with_ns(prob=N_PROB))
                       .run(n=N_SEQUENCES, seed=SEED, productive=True))

        flipped = [r for r in records if not r["productive"]]
        # Pre-fix this was ~25%. Post-fix should be 0 (or very near 0
        # — there's no other corruption path that touches anchors).
        flip_rate = len(flipped) / len(records)
        assert flip_rate == 0.0, (
            f"productive=True + with_ns({N_PROB}) flipped "
            f"{len(flipped)}/{len(records)} ({100*flip_rate:.0f}%) — "
            f"anchor protection is broken")

    def test_productive_without_corruption_near_100pct(self):
        """Sanity: the protection mechanism shouldn't affect runs
        that don't use with_ns at all. The flip rate without N
        corruption should be near zero — there's a tiny known
        caveat (T2-16: SHM can flip productive after the retry
        boundary) so we allow ≤2% slack rather than asserting 100%."""
        records = list(Experiment.on("human_igh")
                       .run(n=N_SEQUENCES, seed=SEED, productive=True))
        n_productive = sum(1 for r in records if r["productive"])
        rate = n_productive / len(records)
        assert rate >= 0.98, (
            f"productive=True without corruption produced "
            f"{n_productive}/{len(records)} ({100*rate:.1f}%) productive "
            f"— broken (T2-16 baseline allows ≥98%)")


class TestNonAnchorPositionsStillCorrupted:
    """Anchor protection shouldn't suppress N-insertion at OTHER
    positions. Real AIRR-seq noise everywhere except the anchors."""

    def test_ns_still_inserted_at_non_anchor_positions(self):
        records = list(Experiment.on("human_igh")
                       .observe(with_ns(prob=N_PROB))
                       .run(n=N_SEQUENCES, seed=SEED, productive=True))

        # With prob=0.05 across ~350 nt sequences, expectation is
        # ~17 N's per sequence. Even with anchor codons (~6 nt total)
        # excluded, the bulk of the sequence still gets the noise.
        seqs_with_ns = sum(1 for r in records if "N" in r["sequence"])
        assert seqs_with_ns >= int(0.95 * len(records)), (
            f"Only {seqs_with_ns}/{len(records)} sequences contain N — "
            f"anchor protection is over-broad and suppressing too much")

        # Also: average N count should be close to expected (~17),
        # confirming we're not skipping huge swaths.
        avg_ns = sum(r["sequence"].count("N") for r in records) / len(records)
        assert avg_ns >= 5.0, (
            f"avg N count is {avg_ns:.1f} — anchor protection is "
            f"over-broad")


class TestMixedModeUnaffected:
    """In MIXED mode there's no productive contract to enforce, so the
    protection deliberately doesn't apply — N can land anywhere
    including anchors. We verify MIXED still works without erroring,
    AND that the protection mechanism doesn't accidentally leak into
    MIXED mode."""

    def test_mixed_mode_still_runs_with_ns(self):
        # productive default = MIXED.
        records = list(Experiment.on("human_igh")
                       .observe(with_ns(prob=N_PROB))
                       .run(n=N_SEQUENCES, seed=SEED))
        # Some productive, some not — that's the whole point of MIXED.
        # Just check we got records back.
        assert len(records) == N_SEQUENCES

    def test_mixed_anchor_can_be_corrupted(self):
        """Empirical sanity check that the protection is conditional
        on PRODUCTIVE_ONLY: in MIXED mode at high N rate, at least
        some V or J anchor codons MUST end up with N — otherwise
        the protection is leaking into MIXED."""
        # Use a much higher rate so the test isn't flaky on small N.
        records = list(Experiment.on("human_igh")
                       .observe(with_ns(prob=0.20))
                       .run(n=N_SEQUENCES, seed=SEED))

        # Count records where the V anchor codon contains an N.
        # The V anchor codon sits at the LAST 3 nts of V.
        n_with_anchor_n = 0
        for r in records:
            v_end = r.get("v_sequence_end", 0)
            if v_end >= 3:
                v_anchor_codon = r["sequence"][v_end-3:v_end].upper()
                if "N" in v_anchor_codon:
                    n_with_anchor_n += 1
        # At 20% per nt, expectation is each codon has 1 - 0.8^3 ≈ 49%
        # chance of containing N. Across 200 sequences, ~98 expected.
        # Assert at least some — confirms protection is OFF in MIXED.
        assert n_with_anchor_n >= 10, (
            f"Only {n_with_anchor_n}/{N_SEQUENCES} MIXED-mode sequences "
            f"have N at V anchor codon — protection is leaking out of "
            f"PRODUCTIVE_ONLY mode")
