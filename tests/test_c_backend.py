"""
Integration tests for the C-backed simulation pipeline.

All tests use the Experiment DSL — every simulation explicitly
declares what happens, no hidden convenience wrappers.
"""

import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def igh_config():
    from GenAIRR.data import HUMAN_IGH_OGRDB
    return HUMAN_IGH_OGRDB


@pytest.fixture
def igk_config():
    from GenAIRR.data import HUMAN_IGK_OGRDB
    return HUMAN_IGK_OGRDB


@pytest.fixture
def tcrb_config():
    from GenAIRR.data import HUMAN_TCRB_IMGT
    return HUMAN_TCRB_IMGT


# ---------------------------------------------------------------------------
# Basic simulation
# ---------------------------------------------------------------------------

class TestBasicSimulation:
    """Test basic rearrangement and output format."""

    def test_rearrangement_returns_result(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=5, seed=42)
        assert len(result) == 5

    def test_result_to_dataframe(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=10, seed=42)
        df = result.to_dataframe()
        assert df.shape[0] == 10
        assert "sequence" in df.columns
        assert "v_call" in df.columns

    def test_airr_fields_present(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=1, seed=42)
        rec = result[0]
        required = [
            "sequence", "v_call", "d_call", "j_call",
            "v_sequence_start", "v_sequence_end",
            "d_sequence_start", "d_sequence_end",
            "j_sequence_start", "j_sequence_end",
            "junction_start", "junction_end",
            "mutation_rate", "productive",
        ]
        for field in required:
            assert field in rec, f"Missing field: {field}"

    def test_sequence_is_nucleotide(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=10, seed=42)
        for rec in result:
            seq = rec["sequence"]
            assert len(seq) > 0
            assert all(c in "ACGTNacgtn" for c in seq), f"Invalid: {seq[:50]}"

    def test_unmutated_rearrangement_only(self):
        from GenAIRR import Experiment
        # No .somatic_hypermutation() → no mutations
        result = Experiment.on("human_igh").run(n=10, seed=42)
        for rec in result:
            assert rec["mutation_rate"] == 0.0
            assert rec["n_mutations"] == 0


# ---------------------------------------------------------------------------
# Productive mode
# ---------------------------------------------------------------------------

class TestProductiveMode:

    def test_mostly_productive(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=100, seed=42, productive=True)
        n_prod = sum(1 for r in result if r["productive"])
        assert n_prod >= 95, f"Only {n_prod}/100 productive"

    def test_productive_no_stop_codons(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=100, seed=42, productive=True)
        productive_recs = [r for r in result if r["productive"]]
        for rec in productive_recs:
            assert rec["stop_codon"] is False

    def test_productive_in_frame(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=100, seed=42, productive=True)
        productive_recs = [r for r in result if r["productive"]]
        for rec in productive_recs:
            assert rec["vj_in_frame"] is True


# ---------------------------------------------------------------------------
# Somatic Hypermutation
# ---------------------------------------------------------------------------

class TestSomaticHypermutation:

    def test_mutation_rate_in_range(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate
        result = (Experiment.on("human_igh")
                  .mutate(rate(0.02, 0.08))
                  .run(n=100, seed=42))
        for rec in result:
            rate = rec["mutation_rate"]
            assert rate >= 0.0
            assert rate <= 0.20  # generous upper bound

    def test_mutations_present(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate
        result = (Experiment.on("human_igh")
                  .mutate(rate(0.01, 0.05))
                  .run(n=50, seed=42))
        n_with_mutations = sum(1 for r in result if r["n_mutations"] > 0)
        assert n_with_mutations > 0

    def test_no_mutations_without_shm(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=50, seed=42)
        for rec in result:
            assert rec["n_mutations"] == 0


# ---------------------------------------------------------------------------
# Segment boundaries
# ---------------------------------------------------------------------------

class TestSegmentBoundaries:

    def test_v_before_d_before_j(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=50, seed=42)
        for rec in result:
            assert rec["v_sequence_start"] <= rec["v_sequence_end"]
            assert rec["j_sequence_start"] <= rec["j_sequence_end"]
            assert rec["v_sequence_end"] <= rec["j_sequence_end"]

    def test_junction_within_sequence(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=50, seed=42)
        for rec in result:
            seq_len = len(rec["sequence"])
            assert rec["junction_start"] >= 0
            assert rec["junction_end"] <= seq_len

    def test_segment_lengths_positive(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igh").run(n=50, seed=42)
        for rec in result:
            v_len = rec["v_sequence_end"] - rec["v_sequence_start"]
            j_len = rec["j_sequence_end"] - rec["j_sequence_start"]
            assert v_len > 0, f"V length {v_len} <= 0"
            assert j_len > 0, f"J length {j_len} <= 0"


# ---------------------------------------------------------------------------
# Experiment compile/reuse API
# ---------------------------------------------------------------------------

class TestExperimentCompileAPI:

    def test_compile_and_simulate(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate
        sim = (Experiment.on("human_igh")
               .mutate(rate(0.01, 0.05))
               .compile(seed=42))
        result = sim.simulate(n=5)
        assert len(result) == 5

    def test_compile_productive(self):
        from GenAIRR import Experiment
        sim = (Experiment.on("human_igh")
               .compile(productive=True, seed=42))
        result = sim.simulate(n=20)
        for rec in result:
            assert rec["productive"] is True

    def test_compiled_simulator_reuse(self):
        from GenAIRR import Experiment
        sim = Experiment.on("human_igh").compile(seed=42)
        r1 = sim.simulate(n=5)
        r2 = sim.simulate(n=5)
        assert len(r1) == 5
        assert len(r2) == 5


# ---------------------------------------------------------------------------
# Corruption features
# ---------------------------------------------------------------------------

class TestCorruptionFeatures:

    def test_flip_strand(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import with_reverse_complement
        result = (Experiment.on("human_igh")
                  .sequence(with_reverse_complement(0.5))
                  .run(n=100, seed=42))
        n_rc = sum(1 for r in result if r["is_reverse_complement"])
        assert n_rc > 20
        assert n_rc < 80

    def test_contaminants(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import with_contaminants
        result = (Experiment.on("human_igh")
                  .observe(with_contaminants(0.1))
                  .run(n=200, seed=42))
        n_contam = sum(1 for r in result if r["is_contaminant"])
        assert n_contam > 0

    def test_indels(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import with_indels
        result = (Experiment.on("human_igh")
                  .observe(with_indels())
                  .run(n=50, seed=42))
        assert len(result) == 50

    def test_pcr_errors(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import with_pcr
        result = (Experiment.on("human_igh")
                  .prepare(with_pcr(error_rate=1e-4, cycles=30))
                  .run(n=50, seed=42))
        n_pcr = sum(1 for r in result if r["n_pcr_errors"] > 0)
        assert n_pcr > 0

    def test_corrupt_5prime(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import with_5prime_loss
        result = (Experiment.on("human_igh")
                  .sequence(with_5prime_loss())
                  .run(n=50, seed=42))
        assert len(result) == 50

    def test_corrupt_3prime(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import with_3prime_loss
        result = (Experiment.on("human_igh")
                  .sequence(with_3prime_loss())
                  .run(n=50, seed=42))
        assert len(result) == 50


# ---------------------------------------------------------------------------
# Multi-species
# ---------------------------------------------------------------------------

class TestMultiSpecies:

    def test_mouse_igh(self):
        from GenAIRR import Experiment
        result = Experiment.on("mouse_igh").run(n=5, seed=42)
        assert len(result) == 5
        assert result[0]["v_call"].startswith("IGHV")

    def test_human_igk(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_igk").run(n=5, seed=42)
        assert len(result) == 5
        assert result[0]["v_call"].startswith("IGK")

    def test_human_tcrb(self):
        from GenAIRR import Experiment
        result = Experiment.on("human_tcrb").run(n=5, seed=42)
        assert len(result) == 5
        assert result[0]["v_call"].startswith("TRBV")


# ---------------------------------------------------------------------------
# Seed reproducibility
# ---------------------------------------------------------------------------

class TestReproducibility:

    def test_same_seed_same_results(self):
        from GenAIRR import Experiment
        r1 = Experiment.on("human_igh").run(n=10, seed=42)
        r2 = Experiment.on("human_igh").run(n=10, seed=42)
        for a, b in zip(r1, r2):
            assert a["sequence"] == b["sequence"]
            assert a["v_call"] == b["v_call"]

    def test_different_seeds_different_results(self):
        from GenAIRR import Experiment
        r1 = Experiment.on("human_igh").run(n=10, seed=42)
        r2 = Experiment.on("human_igh").run(n=10, seed=99)
        n_same = sum(1 for a, b in zip(r1, r2) if a["sequence"] == b["sequence"])
        assert n_same < 10


# ---------------------------------------------------------------------------
# CSR (Class Switch Recombination)
# ---------------------------------------------------------------------------

class TestCSR:

    def test_csr_produces_higher_mutation_rates(self):
        """CSR adjusts mutation rates per isotype — IgG/IgA/IgE get higher rates."""
        from GenAIRR import Experiment
        from GenAIRR.ops import rate, with_isotype_rates
        # Base SHM: low rates
        base_result = (Experiment.on("human_igh")
                       .mutate(rate(0.001, 0.03))
                       .run(n=200, seed=42))
        # CSR: isotype-adjusted rates (IgG1 gets 0.05-0.12, etc.)
        csr_result = (Experiment.on("human_igh")
                      .mutate(rate(0.001, 0.03), with_isotype_rates())
                      .run(n=200, seed=42))

        base_avg = sum(r["mutation_rate"] for r in base_result) / len(base_result)
        csr_avg = sum(r["mutation_rate"] for r in csr_result) / len(csr_result)
        # CSR should produce higher average rates since most isotypes
        # (IgG, IgA, IgE) have higher rates than the base IgM range
        assert csr_avg > base_avg, f"CSR avg {csr_avg:.4f} should exceed base avg {base_avg:.4f}"

    def test_csr_deterministic(self):
        from GenAIRR import Experiment
        from GenAIRR.ops import rate, with_isotype_rates
        def _run():
            return (Experiment.on("human_igh")
                    .mutate(rate(0.001, 0.03), with_isotype_rates())
                    .run(n=50, seed=42))
        r1 = _run()
        r2 = _run()
        rates1 = [r["mutation_rate"] for r in r1]
        rates2 = [r["mutation_rate"] for r in r2]
        assert rates1 == rates2

    def test_csr_c_allele_present(self):
        """CSR sequences should have c_call set."""
        from GenAIRR import Experiment
        from GenAIRR.ops import rate, with_isotype_rates
        result = (Experiment.on("human_igh")
                  .mutate(rate(0.01, 0.05), with_isotype_rates())
                  .run(n=20, seed=42))
        for rec in result:
            assert rec.get("c_call", "") != "", "c_call should be set with CSR"


# ---------------------------------------------------------------------------
# CSimulator direct API
# ---------------------------------------------------------------------------

class TestCSimulatorDirect:

    def test_create_from_gdc_bytes(self):
        from GenAIRR._native import CSimulator, get_version
        from GenAIRR.protocol import _resolve_config
        from GenAIRR.dataconfig.gdc_io import to_gdc_bytes

        config = _resolve_config("human_igh")
        gdc = to_gdc_bytes(config)
        sim = CSimulator(gdc_bytes=gdc)
        assert sim is not None
        sim.destroy()

    def test_version(self):
        from GenAIRR._native import get_version
        v = get_version()
        assert v == "1.0.0"

    def test_simulate_to_file(self, tmp_path):
        from GenAIRR import Experiment
        sim = Experiment.on("human_igh").compile(seed=42)
        path = str(tmp_path / "output.tsv")
        n = sim.simulate_to_file(10, path)
        assert n == 10

        import csv
        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 10
