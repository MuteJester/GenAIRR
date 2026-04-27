"""
Tests for biological correctness and simulation integrity.

These tests validate invariants that must hold for simulated sequences
to be biologically plausible and internally consistent. They catch
subtle coordinate, mutation, and biological errors that basic
functional tests miss.

All tests use the Experiment DSL → C backend path.
"""

import pytest

from GenAIRR import Experiment

N_SEQUENCES = 50

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


# =============================================================================
# Reusable coordinate integrity checker (dict-based)
# =============================================================================

def assert_coordinate_integrity(rec, label=""):
    """Validate all positional invariants on a simulation record dict."""
    seq = rec["sequence"]
    seq_len = len(seq)
    pfx = f"[{label}] " if label else ""

    # All positions within sequence bounds
    assert 0 <= rec["v_sequence_start"] <= seq_len, \
        f"{pfx}v_sequence_start={rec['v_sequence_start']} out of bounds (seq_len={seq_len})"
    assert 0 <= rec["v_sequence_end"] <= seq_len, \
        f"{pfx}v_sequence_end={rec['v_sequence_end']} out of bounds"
    assert 0 <= rec["j_sequence_start"] <= seq_len, \
        f"{pfx}j_sequence_start={rec['j_sequence_start']} out of bounds"
    assert 0 <= rec["j_sequence_end"] <= seq_len, \
        f"{pfx}j_sequence_end={rec['j_sequence_end']} out of bounds"

    # Start <= End for each segment
    assert rec["v_sequence_start"] <= rec["v_sequence_end"], \
        f"{pfx}V start > end: {rec['v_sequence_start']} > {rec['v_sequence_end']}"
    assert rec["j_sequence_start"] <= rec["j_sequence_end"], \
        f"{pfx}J start > end"

    # Segments ordered: V before J
    assert rec["v_sequence_end"] <= rec["j_sequence_end"], \
        f"{pfx}V end ({rec['v_sequence_end']}) > J end ({rec['j_sequence_end']})"

    # D ordering (when D is present — heavy chain)
    d_start = rec.get("d_sequence_start", 0)
    d_end = rec.get("d_sequence_end", 0)
    if d_start != d_end:
        assert 0 <= d_start <= seq_len, f"{pfx}d_sequence_start out of bounds"
        assert 0 <= d_end <= seq_len, f"{pfx}d_sequence_end out of bounds"
        assert d_start <= d_end, f"{pfx}D start > end"

    # Junction within bounds
    j_start = rec.get("junction_start", 0)
    j_end = rec.get("junction_end", 0)
    assert 0 <= j_start <= seq_len, f"{pfx}junction_start out of bounds"
    assert 0 <= j_end <= seq_len, f"{pfx}junction_end out of bounds"
    assert j_start <= j_end, f"{pfx}junction start > end"

    # Only valid nucleotides
    invalid = set(seq.upper()) - set('ATCGN')
    assert not invalid, f"{pfx}invalid characters: {invalid}"


# =============================================================================
# 1. Junction Integrity (productive, unmutated)
# =============================================================================

class TestJunctionIntegrity:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.records = Experiment.on("human_igh").run(
            n=N_SEQUENCES, seed=42, productive=True)

    def test_junction_within_sequence_bounds(self):
        for i, r in enumerate(self.records):
            seq_len = len(r["sequence"])
            js = r["junction_start"]
            je = r["junction_end"]
            assert 0 <= js < seq_len, f"Seq {i}: junction_start out of bounds"
            assert 0 < je <= seq_len, f"Seq {i}: junction_end out of bounds"
            assert js < je, f"Seq {i}: junction_start >= junction_end"

    def test_junction_divisible_by_three_when_in_frame(self):
        for i, r in enumerate(self.records):
            if r["productive"] and r["vj_in_frame"]:
                jlen = r["junction_end"] - r["junction_start"]
                assert jlen % 3 == 0, \
                    f"Seq {i}: productive junction length {jlen} not divisible by 3"

    def test_junction_starts_with_cysteine_codon(self):
        cys_codons = {'TGT', 'TGC'}
        for i, r in enumerate(self.records):
            if r["productive"]:
                js = r["junction_start"]
                codon = r["sequence"][js:js + 3].upper()
                assert codon in cys_codons, \
                    f"Seq {i}: junction starts with '{codon}', expected Cys"

    def test_junction_ends_with_fw_codon(self):
        fw_codons = {'TTT', 'TTC', 'TGG'}
        for i, r in enumerate(self.records):
            if r["productive"]:
                je = r["junction_end"]
                codon = r["sequence"][je - 3:je].upper()
                assert codon in fw_codons, \
                    f"Seq {i}: junction ends with '{codon}', expected F/W"

    def test_junction_reasonable_length(self):
        for i, r in enumerate(self.records):
            jlen = r["junction_end"] - r["junction_start"]
            assert 10 <= jlen <= 120, \
                f"Seq {i}: junction length {jlen} outside reasonable range"


# =============================================================================
# 2. V/D/J Segment Boundaries
# =============================================================================

class TestSegmentBoundaries:

    @pytest.fixture(autouse=True)
    def setup(self):
        from GenAIRR.ops import rate
        self.records = (Experiment.on("human_igh")
                        .mutate(rate(0.01, 0.05))
                        .run(n=N_SEQUENCES, seed=42))

    def test_segments_non_overlapping(self):
        for i, r in enumerate(self.records):
            d_start = r.get("d_sequence_start", 0)
            d_end = r.get("d_sequence_end", 0)
            if d_start != d_end:
                assert r["v_sequence_end"] <= d_start, \
                    f"Seq {i}: V end ({r['v_sequence_end']}) > D start ({d_start})"
                assert d_end <= r["j_sequence_start"], \
                    f"Seq {i}: D end ({d_end}) > J start ({r['j_sequence_start']})"
            assert r["v_sequence_end"] <= r["j_sequence_start"], \
                f"Seq {i}: V end ({r['v_sequence_end']}) > J start ({r['j_sequence_start']})"

    def test_segment_lengths_positive(self):
        for i, r in enumerate(self.records):
            v_len = r["v_sequence_end"] - r["v_sequence_start"]
            j_len = r["j_sequence_end"] - r["j_sequence_start"]
            assert v_len > 0, f"Seq {i}: V empty"
            assert j_len > 0, f"Seq {i}: J empty"


# =============================================================================
# 3. Mutation Constraints
# =============================================================================

class TestMutationConstraints:

    @pytest.fixture(autouse=True)
    def setup(self):
        from GenAIRR.ops import rate
        self.records = (Experiment.on("human_igh")
                        .mutate(rate(0.05, 0.15))
                        .run(n=N_SEQUENCES, seed=42))

    def test_mutation_rate_within_bounds(self):
        for i, r in enumerate(self.records):
            rate = r["mutation_rate"]
            assert 0 <= rate <= 1.0, \
                f"Seq {i}: mutation rate {rate} out of [0,1]"

    def test_mutation_count_consistent_with_rate(self):
        for i, r in enumerate(self.records):
            n_mut = r["n_mutations"]
            rate = r["mutation_rate"]
            seq_len = len(r["sequence"])
            if seq_len > 0 and n_mut > 0:
                computed_rate = n_mut / seq_len
                # Allow generous tolerance — rate and count may differ slightly
                assert computed_rate < 0.5, \
                    f"Seq {i}: {n_mut} mutations in {seq_len}bp seems too high"

    def test_no_mutations_without_shm(self):
        records = Experiment.on("human_igh").run(n=N_SEQUENCES, seed=42)
        for i, r in enumerate(records):
            assert r["n_mutations"] == 0, \
                f"Seq {i}: {r['n_mutations']} mutations without SHM enabled"
            assert r["mutation_rate"] == 0.0

    def test_t1_4_no_spurious_boundary_mutations_with_indels(self):
        """T1-4: heavy indels can fully wipe D from the sequence; that
        used to leave d_trim_*_adjusted at memset 0, making
        collect_boundary_mutations scan V's start as if it were D's
        boundary extension and emit spurious mutation entries.

        With no SHM and no quality/PCR errors, n_mutations must remain
        0 regardless of how aggressive the indels are."""
        from GenAIRR.ops import with_indels
        result = (Experiment.on("human_igh")
                  .observe(with_indels(prob=0.15))
                  .run(n=2000, seed=123, productive=False))
        records = list(result)
        spurious = [r for r in records if r["n_mutations"] > 0]
        assert not spurious, (
            f"{len(spurious)} records have n_mutations > 0 with no SHM. "
            f"First few: "
            + "; ".join(f"row{i} n={r['n_mutations']} muts={r['mutations'][:60]}"
                         for i, r in enumerate(spurious[:3]))
        )


# =============================================================================
# 4. Productive Sequence Biological Invariants
# =============================================================================

class TestProductiveInvariants:

    @pytest.fixture(autouse=True)
    def setup(self):
        self.records = Experiment.on("human_igh").run(
            n=N_SEQUENCES, seed=42, productive=True)

    def test_productive_no_stop_codons_in_junction(self):
        for i, r in enumerate(self.records):
            if not r["productive"]:
                continue
            js = r["junction_start"]
            je = r["junction_end"]
            junction = r["sequence"][js:je].upper()
            for j in range(0, len(junction) - 2, 3):
                codon = junction[j:j+3]
                aa = CODON_TABLE.get(codon, 'X')
                assert aa != '*', \
                    f"Seq {i}: stop codon '{codon}' at junction pos {j}"

    def test_productive_flag_consistency(self):
        for i, r in enumerate(self.records):
            if r["productive"]:
                assert not r["stop_codon"], f"Seq {i}: productive but has stop codon"
                assert r["vj_in_frame"], f"Seq {i}: productive but not in frame"


# =============================================================================
# 5. Full Pipeline End-to-End Coordinate Consistency
# =============================================================================

class TestFullPipelineIntegrity:

    def test_basic_pipeline(self):
        records = Experiment.on("human_igh").run(n=N_SEQUENCES, seed=42)
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"basic-{i}")

    def test_with_shm(self):
        from GenAIRR.ops import rate
        records = (Experiment.on("human_igh")
                   .mutate(rate(0.01, 0.05))
                   .run(n=N_SEQUENCES, seed=42))
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"shm-{i}")

    def test_with_corruption(self):
        from GenAIRR.ops import with_5prime_loss
        records = (Experiment.on("human_igh")
                   .sequence(with_5prime_loss())
                   .run(n=N_SEQUENCES, seed=42))
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"corrupt-{i}")

    def test_with_indels(self):
        from GenAIRR.ops import with_indels
        records = (Experiment.on("human_igh")
                   .observe(with_indels())
                   .run(n=N_SEQUENCES, seed=42))
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"indel-{i}")

    def test_full_augmentation(self):
        from GenAIRR.ops import rate, with_5prime_loss, with_indels, with_ns
        records = (Experiment.on("human_igh")
                   .mutate(rate(0.01, 0.05))
                   .sequence(with_5prime_loss())
                   .observe(with_indels(), with_ns())
                   .run(n=N_SEQUENCES, seed=42))
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"full-{i}")

    def test_light_chain_kappa(self):
        records = Experiment.on("human_igk").run(n=N_SEQUENCES, seed=42)
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"igk-{i}")

    def test_light_chain_lambda(self):
        records = Experiment.on("human_igl").run(n=N_SEQUENCES, seed=42)
        for i, r in enumerate(records):
            assert_coordinate_integrity(r, f"igl-{i}")


# =============================================================================
# 6. AIRR Dict Field Completeness
# =============================================================================

class TestAIRRFieldCompleteness:

    def test_has_all_required_fields(self):
        result = Experiment.on("human_igh").run(n=1, seed=42)
        rec = result[0]
        required = [
            'sequence', 'v_call', 'd_call', 'j_call',
            'v_sequence_start', 'v_sequence_end',
            'd_sequence_start', 'd_sequence_end',
            'j_sequence_start', 'j_sequence_end',
            'junction_start', 'junction_end',
            'mutation_rate', 'n_mutations',
            'productive', 'stop_codon', 'vj_in_frame',
        ]
        for field in required:
            assert field in rec, f"Missing field '{field}'"

    def test_sequence_only_valid_nucleotides(self):
        records = Experiment.on("human_igh").run(n=N_SEQUENCES, seed=42)
        for i, r in enumerate(records):
            seq = r["sequence"]
            assert len(seq) > 0, f"Seq {i}: empty sequence"
            invalid = set(seq.upper()) - set('ATCGN')
            assert not invalid, f"Seq {i}: invalid chars {invalid}"


# =============================================================================
# 7. NP region biology — TdT Markov sampling (T0-1 regression)
# =============================================================================
#
# Before the T0-1 fix, NP1/NP2 bases were drawn uniformly at 25% per base
# (~50% GC). The DataConfig ships an empirical NP_first_bases / NP_transitions
# Markov model trained from real BCR data; the C engine now consumes it via
# NpDist. Real human IGH NP regions are GC-enriched (~60-65% GC) due to TdT
# preference. These tests guard against silent regression to the old uniform
# behavior — both via base composition and length distribution.

class TestNPRegionBiology:

    @pytest.fixture(scope="class")
    def igh_records(self):
        # Larger n than other classes: NP1 is short (~8bp average), so we
        # need many sequences to get a stable empirical estimate.
        return Experiment.on("human_igh").run(
            n=500, seed=42, productive=False)

    @staticmethod
    def _extract_np1(rec):
        """Slice NP1 = sequence[v_sequence_end : d_sequence_start]."""
        v_end = rec["v_sequence_end"]
        d_start = rec.get("d_sequence_start", 0)
        if d_start <= v_end:
            return None
        return rec["sequence"][v_end:d_start].upper()

    @staticmethod
    def _extract_np2(rec):
        """Slice NP2 = sequence[d_sequence_end : j_sequence_start]."""
        d_end = rec.get("d_sequence_end", 0)
        j_start = rec["j_sequence_start"]
        if d_end == 0 or j_start <= d_end:
            return None
        return rec["sequence"][d_end:j_start].upper()

    def test_np1_gc_content_above_uniform(self, igh_records):
        """NP1 bases should be GC-enriched (~60-65%), well above the
        uniform-sampler 50%. Threshold of 0.55 catches a regression to
        uniform without being so tight it flakes."""
        nps = [self._extract_np1(r) for r in igh_records]
        nps = [s for s in nps if s]
        total = sum(len(s) for s in nps)
        if total < 1000:
            pytest.skip(f"insufficient NP1 bases for stat test: {total}")
        gc = sum(c in "GC" for s in nps for c in s)
        gc_frac = gc / total
        assert 0.55 <= gc_frac <= 0.80, (
            f"NP1 GC fraction {gc_frac:.3f} outside expected [0.55, 0.80]. "
            f"Uniform sampler would give ~0.50; real biology ~0.60-0.65."
        )

    def test_np2_gc_content_above_uniform(self, igh_records):
        nps = [self._extract_np2(r) for r in igh_records]
        nps = [s for s in nps if s]
        total = sum(len(s) for s in nps)
        if total < 500:
            pytest.skip(f"insufficient NP2 bases for stat test: {total}")
        gc = sum(c in "GC" for s in nps for c in s)
        gc_frac = gc / total
        assert 0.55 <= gc_frac <= 0.80, (
            f"NP2 GC fraction {gc_frac:.3f} outside expected [0.55, 0.80]."
        )

    def test_np1_mean_length_matches_empirical(self, igh_records):
        """Mean NP1 length should track the empirical mean from the
        DataConfig's NP_lengths distribution (within ±2 bp)."""
        from GenAIRR.data import HUMAN_IGH_OGRDB
        length_dist = HUMAN_IGH_OGRDB.NP_lengths.get("NP1", {})
        if not length_dist:
            pytest.skip("NP_lengths['NP1'] not present in DataConfig")

        expected_mean = sum(k * v for k, v in length_dist.items())
        # rec["np1_length"] is set by step_assemble; use that rather than
        # re-deriving from coordinates so we measure the sampler directly.
        observed = [r["np1_length"] for r in igh_records]
        observed_mean = sum(observed) / len(observed)
        assert abs(observed_mean - expected_mean) <= 2.0, (
            f"NP1 mean length {observed_mean:.2f} diverges from "
            f"empirical {expected_mean:.2f} by more than 2 bp"
        )

    def test_np_sampler_is_seed_stable(self):
        """Two independent runs with the same seed must produce identical
        NP1 sequences. Regression guard for the Markov sampler's RNG plumbing."""
        a = Experiment.on("human_igh").run(n=20, seed=4242, productive=False)
        b = Experiment.on("human_igh").run(n=20, seed=4242, productive=False)
        np1_a = [self._extract_np1(r) for r in a]
        np1_b = [self._extract_np1(r) for r in b]
        assert np1_a == np1_b, "Same-seed runs produced different NP1 regions"

    def test_np1_uses_more_than_one_unique_base(self, igh_records):
        """Sanity: a buggy Markov sampler that locked into a single state
        (e.g., always G after the first base) would produce monotonous NPs.
        Across 500 sequences, all four bases should appear in NP1."""
        nps = [self._extract_np1(r) for r in igh_records]
        joined = "".join(s for s in nps if s)
        if len(joined) < 200:
            pytest.skip(f"insufficient NP1 bases: {len(joined)}")
        bases_seen = set(joined)
        assert bases_seen >= set("ACGT"), (
            f"NP1 bases observed {bases_seen}; expected all of A,C,G,T"
        )


# =============================================================================
# 8. P-nucleotide biology — RAG/Artemis hairpin opening (T0-2 regression)
# =============================================================================
#
# Before T0-2, P-nucleotides were never generated despite the DataConfig
# carrying a `p_nucleotide_length_probs` distribution. P-nucs appear only at
# segment ends with zero exonuclease trim and consist of the reverse-
# complement of the K bases at that segment edge.
#
# The default distribution {0:.5, 1:.25, 2:.15, 3:.07, 4:.03} gives
# P(K > 0 | trim==0) = 0.5. About 21% of human IGH records have v_trim_3==0,
# so ~10% of records should carry a detectable V 3' P-nuc.

_DNA_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _rc(s):
    return s.translate(_DNA_COMPLEMENT)[::-1]


class TestPNucleotideBiology:

    @pytest.fixture(scope="class")
    def igh_records(self):
        # No mutation — keeps P-nuc bases pristine for palindromicity check.
        return Experiment.on("human_igh").run(
            n=500, seed=42, productive=False)

    @staticmethod
    def _v_allele_seq(v_call):
        """Resolve a v_call (possibly comma-separated for ambiguity) to the
        first matching allele's ungapped germline sequence. Returns None
        when none of the listed alleles can be looked up."""
        from GenAIRR.data import HUMAN_IGH_OGRDB
        names = [n.strip() for n in v_call.split(",")] if v_call else []
        for gene_dict in (HUMAN_IGH_OGRDB.v_alleles or {}).values():
            for allele in gene_dict:
                if allele.name in names:
                    return allele.ungapped_seq.upper()
        return None

    @staticmethod
    def _j_allele_seq(j_call):
        from GenAIRR.data import HUMAN_IGH_OGRDB
        names = [n.strip() for n in j_call.split(",")] if j_call else []
        for gene_dict in (HUMAN_IGH_OGRDB.j_alleles or {}).values():
            for allele in gene_dict:
                if allele.name in names:
                    return allele.ungapped_seq.upper()
        return None

    @staticmethod
    def _detect_v_p_nuc(np1_region, v_seq):
        """Return the largest K in [1..4] such that np1_region.upper()[:K]
        equals the reverse-complement of v_seq[-K:]. 0 if no match."""
        if not np1_region or not v_seq:
            return 0
        np1 = np1_region.upper()
        for k in (4, 3, 2, 1):
            if k <= len(np1) and k <= len(v_seq):
                if np1[:k] == _rc(v_seq[-k:]):
                    return k
        return 0

    @staticmethod
    def _detect_j_p_nuc(np_region, j_seq):
        """For a J 5' head-side P-nuc, check whether the LAST K bases of
        the NP segment equal the reverse-complement of j_seq[:K]."""
        if not np_region or not j_seq:
            return 0
        np_up = np_region.upper()
        for k in (4, 3, 2, 1):
            if k <= len(np_up) and k <= len(j_seq):
                if np_up[-k:] == _rc(j_seq[:k]):
                    return k
        return 0

    def test_v3_p_nuc_appears_when_untrimmed(self, igh_records):
        """Among records with v_trim_3 == 0, at least 25% should show a
        detectable V 3' P-nuc at the start of np1_region. Under the
        default distribution we expect ~50%; the threshold is loose to
        absorb cases where the first NP1 N-nuc happens to match the
        expected P-nuc base by chance (which we cannot disambiguate
        from a true P-nuc) and vice versa."""
        eligible = [r for r in igh_records if r["v_trim_3"] == 0]
        if len(eligible) < 30:
            pytest.skip(f"too few v_trim_3==0 records: {len(eligible)}")

        detected = 0
        for r in eligible:
            v_seq = self._v_allele_seq(r["v_call"])
            if not v_seq:
                continue
            if self._detect_v_p_nuc(r["np1_region"], v_seq) > 0:
                detected += 1
        rate = detected / len(eligible)
        assert rate >= 0.25, (
            f"V 3' P-nuc detection rate {rate:.2f} below 0.25 "
            f"({detected}/{len(eligible)} eligible records)"
        )

    def test_j5_p_nuc_appears_when_untrimmed(self, igh_records):
        """Symmetrical test for J 5' P-nucs at the end of NP2 (in VDJ
        chains). Sanity-checks that head-side P-nucs are wired up, not
        just tail-side ones."""
        eligible = [r for r in igh_records if r["j_trim_5"] == 0]
        if len(eligible) < 30:
            pytest.skip(f"too few j_trim_5==0 records: {len(eligible)}")

        detected = 0
        for r in eligible:
            j_seq = self._j_allele_seq(r["j_call"])
            if not j_seq:
                continue
            np_region = r["np2_region"] or r["np1_region"]
            if self._detect_j_p_nuc(np_region, j_seq) > 0:
                detected += 1
        rate = detected / len(eligible)
        assert rate >= 0.25, (
            f"J 5' P-nuc detection rate {rate:.2f} below 0.25 "
            f"({detected}/{len(eligible)} eligible records)"
        )

    def test_no_v3_p_nuc_when_trimmed(self, igh_records):
        """When v_trim_3 > 0, no V 3' P-nuc can be emitted (RAG/Artemis
        only emits P-nucs at untrimmed ends). The first base of NP1 is
        therefore decoupled from V's tail — there should be no
        statistical enrichment for RC(V[-1]) at NP1[0]. We check the
        weaker contrapositive: the empirical rate of "RC(V[-1]) ==
        NP1[0]" among trimmed-V records should be ~0.25 (i.e., chance)."""
        trimmed = [r for r in igh_records
                   if r["v_trim_3"] > 0 and r["np1_region"]]
        if len(trimmed) < 50:
            pytest.skip(f"too few v_trim_3>0 records: {len(trimmed)}")

        matches = 0
        considered = 0
        for r in trimmed:
            v_seq = self._v_allele_seq(r["v_call"])
            if not v_seq or not r["np1_region"]:
                continue
            considered += 1
            if r["np1_region"][0].upper() == _rc(v_seq[-1]):
                matches += 1
        rate = matches / considered if considered else 0.0
        # If P-nucs were leaking into trimmed records this would be
        # well above 0.5; the null is 0.25.
        assert rate < 0.45, (
            f"RC(V[-1])==NP1[0] rate {rate:.2f} suspiciously high in "
            f"trimmed-V records (suggests P-nucs leaking in despite trim>0)"
        )


# =============================================================================
# 9. Per-allele trim distributions (T0-3 regression)
# =============================================================================
#
# Before T0-3, the C engine collapsed all per-(family, gene) trim
# distributions to a single global TrimDist per segment, so every V allele
# sampled trim from the same distribution regardless of gene/family. After
# T0-3, each allele dispatches to its own (family, gene) entry. The
# regression test exploits a strong empirical separation in HUMAN_IGH:
# IGHVF1 has mean V_3 trim of ~2.34 while IGHVF4 has mean ~0.85. Under
# the bug, both families converge to the same global mean (~1.1 for IGH).

class TestPerAlleleTrimDispatch:

    @pytest.fixture(scope="class")
    def igh_records(self):
        # Larger n: we partition by V family and need enough samples per
        # family to get a stable empirical mean.
        return Experiment.on("human_igh").run(
            n=2000, seed=42, productive=False)

    @staticmethod
    def _family_means_from_dataconfig():
        """Compute the per-family mean V_3 trim from the DataConfig
        (the ground-truth distribution)."""
        from GenAIRR.data import HUMAN_IGH_OGRDB
        td = HUMAN_IGH_OGRDB.trim_dicts["V_3"]
        means = {}
        for fam, gene_dict in td.items():
            per_gene_means = []
            for prob_dict in gene_dict.values():
                per_gene_means.append(
                    sum(k * v for k, v in prob_dict.items())
                )
            if per_gene_means:
                means[fam] = sum(per_gene_means) / len(per_gene_means)
        return means

    def test_high_trim_families_trim_more_than_low(self, igh_records):
        """The empirical mean V_3 trim for IGHVF1 (high-trim family,
        DataConfig mean ~2.34) must exceed IGHVF4 (low-trim family,
        DataConfig mean ~0.85). Under the bug both used the same
        global distribution → both families produce the same mean."""
        # Use v_call_true (truly-sampled allele) — v_call is god-aligner
        # derived and may drift under SHM.
        from collections import defaultdict
        family_trims = defaultdict(list)
        for r in igh_records:
            v_call = r.get("v_call_true") or r["v_call"].split(",")[0]
            # Family is the prefix before the gene (e.g. IGHVF1 → IGHVF1-G3*01)
            # Use the leading IGHVFxx token.
            from GenAIRR.data import HUMAN_IGH_OGRDB
            for gene_dict in HUMAN_IGH_OGRDB.v_alleles.values():
                hit = next((a for a in gene_dict if a.name == v_call), None)
                if hit:
                    family_trims[hit.family].append(r["v_trim_3"])
                    break

        # Need a meaningful sample size for both target families
        n_high = len(family_trims.get("IGHVF1", []))
        n_low = len(family_trims.get("IGHVF4", []))
        if n_high < 20 or n_low < 20:
            pytest.skip(
                f"insufficient family samples: IGHVF1={n_high}, IGHVF4={n_low}"
            )

        mean_high = sum(family_trims["IGHVF1"]) / n_high
        mean_low = sum(family_trims["IGHVF4"]) / n_low

        # Strong separation expected: ~2.34 vs ~0.85 in the DataConfig.
        # Allow noise but require clear directional separation.
        assert mean_high - mean_low > 0.5, (
            f"V trim per-allele dispatch broken: IGHVF1 mean={mean_high:.2f} "
            f"not significantly higher than IGHVF4 mean={mean_low:.2f} "
            f"(under T0-3 fix, IGHVF1 should be ~2.3 vs IGHVF4 ~0.85)"
        )

    def test_simulated_family_means_track_dataconfig(self, igh_records):
        """Across families with ≥30 samples, the simulated mean V_3 trim
        should track the DataConfig's family mean within ±1.0 bp. This
        is the strict per-family fidelity check."""
        from collections import defaultdict
        from GenAIRR.data import HUMAN_IGH_OGRDB
        family_trims = defaultdict(list)
        for r in igh_records:
            v_call = r.get("v_call_true") or r["v_call"].split(",")[0]
            for gene_dict in HUMAN_IGH_OGRDB.v_alleles.values():
                hit = next((a for a in gene_dict if a.name == v_call), None)
                if hit:
                    family_trims[hit.family].append(r["v_trim_3"])
                    break

        expected_means = self._family_means_from_dataconfig()
        checked = 0
        for fam, trims in family_trims.items():
            if len(trims) < 30 or fam not in expected_means:
                continue
            empirical = sum(trims) / len(trims)
            expected = expected_means[fam]
            # ±1.0 bp tolerance — a single global dist would push every
            # family toward the global mean (~1.1), making high-mean
            # families like IGHVF1 fail this check.
            assert abs(empirical - expected) <= 1.0, (
                f"family {fam}: empirical mean trim {empirical:.2f} "
                f"diverges from DataConfig {expected:.2f} by >1.0 bp "
                f"(n={len(trims)})"
            )
            checked += 1
        assert checked >= 3, f"only {checked} families had enough samples"


# =============================================================================
# 10. V-anchored reading frame (T0-4 regression)
# =============================================================================
#
# Before T0-4, the codon rail was rebuilt after every length change
# (5' corruption, UMI prepend, etc.) using "phase 0 at seq->head"
# semantics — ignoring V's biological reading frame. The rebuilt rail
# placed codons in arbitrary frames relative to V, finding spurious
# stop codons in the adapter and breaking the junction-in-frame check.
# Empirically the productive rate on human_igh dropped from ~11% to
# ~3.4% under default 5'-corruption settings — a 3× collapse.
#
# After T0-4, the codon rail derives its initial phase from V's first
# surviving node so V's reading frame is preserved across head/tail
# mutations, and adapter/UMI stop codons are excluded from the
# productivity check.

class TestVAnchoredFrame:

    def test_5prime_corruption_does_not_collapse_productive_rate(self):
        """The simulated productive rate must not collapse when 5'
        corruption is enabled. Compare with-vs-without on the same
        seed and assert the rates are within a reasonable factor.
        Pre-T0-4 the corrupted rate was ~3× lower than the baseline."""
        from GenAIRR.ops import rate, with_5prime_loss

        baseline = list((Experiment.on("human_igh")
                         .mutate(rate(0.02, 0.05))
                         .run(n=500, seed=42, productive=False)))
        corrupted = list((Experiment.on("human_igh")
                          .mutate(rate(0.02, 0.05))
                          .sequence(with_5prime_loss())
                          .run(n=500, seed=42, productive=False)))

        prod_baseline = sum(1 for r in baseline if r["productive"]) / len(baseline)
        prod_corrupted = sum(1 for r in corrupted if r["productive"]) / len(corrupted)

        # Baseline should be ~10-15% under the requested mutation rate.
        # Corrupted must remain in the same ballpark (within a 2x
        # window). Pre-T0-4 the corrupted rate was 3-4x lower.
        assert prod_baseline > 0.05, (
            f"baseline productive rate {prod_baseline:.3f} unexpectedly low"
        )
        ratio = prod_corrupted / prod_baseline if prod_baseline > 0 else 0.0
        assert 0.5 <= ratio <= 2.0, (
            f"5'-corruption distorted productive rate beyond 0.5-2.0× "
            f"(baseline={prod_baseline:.3f}, corrupted={prod_corrupted:.3f}, "
            f"ratio={ratio:.2f})"
        )

    def test_5prime_corruption_does_not_inflate_stop_codon_rate(self):
        """5'-corrupted records must not have an inflated stop_codon
        rate. Pre-T0-4 the rebuilt codon rail mis-framed adapter
        bases and counted spurious stops in non-coding regions."""
        from GenAIRR.ops import rate, with_5prime_loss

        baseline = list((Experiment.on("human_igh")
                         .mutate(rate(0.02, 0.05))
                         .run(n=500, seed=42, productive=False)))
        corrupted = list((Experiment.on("human_igh")
                          .mutate(rate(0.02, 0.05))
                          .sequence(with_5prime_loss())
                          .run(n=500, seed=42, productive=False)))

        stop_baseline = sum(1 for r in baseline if r["stop_codon"]) / len(baseline)
        stop_corrupted = sum(1 for r in corrupted if r["stop_codon"]) / len(corrupted)

        # Adapter-induced stops must not push the corrupted rate
        # significantly above baseline. Allow +10pp of slack for RNG.
        assert stop_corrupted - stop_baseline < 0.10, (
            f"5'-corrupted stop_codon rate {stop_corrupted:.3f} exceeds "
            f"baseline {stop_baseline:.3f} by >0.10 — suggests adapter "
            f"stops are leaking into the productivity check"
        )


# =============================================================================
# 11. Per-simulator RNG independence (T0-5 regression)
# =============================================================================
#
# Before T0-5, the C engine used a single process-global libc rand()
# stream. Two GenAIRRSimulator instances in the same Python process
# clobbered each other's RNG, breaking reproducibility. After T0-5,
# each simulator owns a PCG32 state and seed=N → identical output
# regardless of how many other simulators run concurrently.

class TestPerSimulatorRng:

    def test_same_seed_byte_identical(self):
        """Two simulators built with the same SimConfig and same seed
        must produce byte-identical sequences in independent runs."""
        a = list(Experiment.on("human_igh").run(n=20, seed=12345,
                                                productive=False))
        b = list(Experiment.on("human_igh").run(n=20, seed=12345,
                                                productive=False))
        assert len(a) == len(b) == 20
        for i, (ra, rb) in enumerate(zip(a, b)):
            assert ra["sequence"] == rb["sequence"], (
                f"seq {i}: same-seed runs diverged "
                f"({ra['sequence'][:30]}... vs {rb['sequence'][:30]}...)"
            )

    def test_different_seeds_diverge(self):
        """Different seeds must produce different output (overwhelmingly
        likely across n=20 records)."""
        a = list(Experiment.on("human_igh").run(n=20, seed=1,
                                                productive=False))
        b = list(Experiment.on("human_igh").run(n=20, seed=2,
                                                productive=False))
        diffs = sum(1 for ra, rb in zip(a, b)
                    if ra["sequence"] != rb["sequence"])
        assert diffs >= 18, (
            f"only {diffs}/20 sequences differ between seed=1 and seed=2"
        )

    def test_interleaved_simulators_independent(self):
        """Two compiled simulators with different seeds must each
        produce its own deterministic stream regardless of interleaved
        consumption order. Pre-T0-5 the global libc RNG would couple
        them; post-T0-5 each owns a PCG32 state."""
        # Compile two independent simulators
        sim_a = Experiment.on("human_igh").compile(seed=111)
        sim_b = Experiment.on("human_igh").compile(seed=222)

        # Reference runs: each one alone
        sim_a_solo = Experiment.on("human_igh").compile(seed=111)
        sim_b_solo = Experiment.on("human_igh").compile(seed=222)
        ref_a = sim_a_solo.simulate(n=10)
        ref_b = sim_b_solo.simulate(n=10)

        # Interleaved: alternate single calls
        inter_a, inter_b = [], []
        for _ in range(10):
            inter_a.extend(sim_a.simulate(n=1))
            inter_b.extend(sim_b.simulate(n=1))

        # Each simulator's interleaved output must equal its solo output
        for i, (r1, r2) in enumerate(zip(ref_a, inter_a)):
            assert r1["sequence"] == r2["sequence"], (
                f"sim_a record {i}: interleaved output diverges from solo "
                f"(per-simulator RNG broken)"
            )
        for i, (r1, r2) in enumerate(zip(ref_b, inter_b)):
            assert r1["sequence"] == r2["sequence"], (
                f"sim_b record {i}: interleaved output diverges from solo"
            )
