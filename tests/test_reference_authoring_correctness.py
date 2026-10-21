"""C0 regressions for reference authoring.

See audit-docs/reference_curation_pipeline_design_2026_10_21.md, §2 and §12 C0.

Three defects in the pre-catalogue authoring path:

1. The DataConfig bridge assigned fixed IMGT subregions to any V whose
   ``gapped_seq`` was non-empty, including ungapped FASTA with no declared
   numbering alignment.
2. ``infer_identity()`` overwrote loader/header provenance on every allele
   and stamped wall-clock dates, so identical inputs built different
   cartridges on different days.
3. The PyO3 bridge narrowed sequence lengths to ``u16`` without a check.
"""
from datetime import date, datetime

import pytest

from GenAIRR import _engine
from GenAIRR.alleles import VAllele
from GenAIRR.cartridge_builder import ReferenceCartridgeBuilder
from GenAIRR._refdata_resolver import dataconfig_to_refdata
from GenAIRR.utilities.imgt_regions import compute_v_region_boundaries

# 105 codons with no stop: an ordinary ungapped V-length sequence.
BARE_V_315 = "GCT" * 104 + "TGT"
J_FASTA = ">TESTJ1*01\nTGGTTCGGCGGAGGGACCAAGCTGACCGTCCTAG\n"
# IMGT-gapped V: three gap dots inside FWR1; Cys104 at gapped 309..312.
GAPPED_V_315 = "AAA" + "..." + "AAA" * 101 + "TGT" + "AAA"
IMGT_V_HEADER = ">ACC1|TRAV1*01|Mus musculus|P|V-REGION|"
FIXED_DATE = date(2026, 1, 2)
NATIVE_MAX_ALLELE_LEN = 65_535


def _bare_builder(chain_type="BCR_LIGHT_KAPPA"):
    return ReferenceCartridgeBuilder.from_fasta(
        v_fasta=">TESTV1*01\n" + BARE_V_315 + "\n", j_fasta=J_FASTA, chain_type=chain_type
    )


def _imgt_header_builder():
    return ReferenceCartridgeBuilder.from_fasta(
        v_fasta=IMGT_V_HEADER + "\n" + GAPPED_V_315 + "\n",
        j_fasta=">ACC2|TRAJ1*01|Mus musculus|F|J-REGION|\nTTTGGGAAAGGG\n",
        chain_type="TCR_ALPHA",
    )


def _only_v(cfg):
    (alleles,) = cfg.v_alleles.values()
    (allele,) = alleles
    return allele


# ── 1. Undeclared IMGT subregion inference ─────────────────────────────


def test_bare_fasta_v_receives_no_imgt_subregions_in_native_reference():
    assert len(BARE_V_315) == 315 and "." not in BARE_V_315
    native = dataconfig_to_refdata(_bare_builder().build())
    assert native.v_allele(0).subregions == []


def test_declared_imgt_gapped_v_keeps_derived_subregions():
    # Hand-derived: IMGT gapped bounds 78/114/165/195/312 shifted by the
    # three FWR1 gap dots; FWR1 starts at ungapped 0.
    native = dataconfig_to_refdata(_imgt_header_builder().build())
    assert native.v_allele(0).subregions == [
        ("FWR1", 0, 75), ("CDR1", 75, 111), ("FWR2", 111, 162),
        ("CDR2", 162, 192), ("FWR3", 192, 309),
    ]


def test_explicit_subregions_on_bare_v_survive_the_bridge():
    cfg = _bare_builder().build()
    _only_v(cfg).subregions = {"FWR1": (0, 60), "CDR1": (60, 90)}
    native = dataconfig_to_refdata(cfg)
    assert native.v_allele(0).subregions == [("FWR1", 0, 60), ("CDR1", 60, 90)]


def test_imgt_boundary_helper_refuses_sequence_without_declared_gaps():
    with pytest.raises(ValueError, match="IMGT"):
        compute_v_region_boundaries(VAllele("TESTV1*01", BARE_V_315, len(BARE_V_315)))


# ── 2. Provenance and deterministic build metadata ─────────────────────


def test_infer_identity_preserves_header_provenance():
    cfg = _imgt_header_builder().infer_identity(species="HUMAN", source="lab-build").build()
    v = _only_v(cfg)
    assert v.source == "ACC1"
    assert v.species == "Mus musculus"
    assert cfg.metadata.species.name == "HUMAN"


def test_infer_identity_reports_preserved_species_labels():
    # Header labels and the Species enum use different vocabularies, so C0
    # reports what it preserved rather than guessing a disagreement.
    cfg = _imgt_header_builder().infer_identity(species="HUMAN").build()
    stage = next(s for s in cfg.build_report.stages if s["stage"] == "infer_identity")
    assert stage["inferred"]["preserved_allele_species"] == {"Mus musculus": 2}


def test_infer_identity_fills_missing_provenance():
    cfg = _bare_builder().infer_identity(species="HUMAN", source="lab-build").build()
    v = _only_v(cfg)
    assert (v.source, v.species, v.locus) == ("lab-build", "Human", "IGK")


def test_repeated_infer_identity_updates_only_builder_filled_provenance():
    bare = _bare_builder().infer_identity(species="HUMAN", source="first")
    v = _only_v(bare.infer_identity(species="MOUSE", source="second").build())
    assert (v.source, v.species) == ("second", "Mouse")
    header = _imgt_header_builder().infer_identity(species="HUMAN", source="first")
    v = _only_v(header.infer_identity(species="MOUSE", source="second").build())
    assert (v.source, v.species) == ("ACC1", "Mus musculus")


@pytest.mark.parametrize("with_identity", [True, False])
def test_fixed_build_date_sets_cartridge_date(with_identity):
    builder = _bare_builder()
    if with_identity:
        builder = builder.infer_identity(species="HUMAN")
    cfg = builder.with_build_date(FIXED_DATE).build()
    assert cfg.metadata.last_updated == FIXED_DATE


def test_build_date_applies_regardless_of_call_order():
    cfg = _bare_builder().with_build_date(FIXED_DATE).infer_identity(species="HUMAN").build()
    assert cfg.metadata.last_updated == FIXED_DATE


def test_fixed_build_date_reproduces_cartridge_checksum():
    first = _bare_builder().infer_identity(species="HUMAN").with_build_date(FIXED_DATE).build()
    second = _bare_builder().infer_identity(species="HUMAN").with_build_date(FIXED_DATE).build()
    assert first.compute_checksum() == second.compute_checksum()


@pytest.mark.parametrize("value", ["2026-01-02", datetime(2026, 1, 2, 3, 4), None])
def test_build_date_must_be_a_calendar_date(value):
    with pytest.raises(TypeError):
        _bare_builder().with_build_date(value)


# ── 3. Native coordinate narrowing ─────────────────────────────────────


@pytest.mark.parametrize("segment", ["v", "d", "j"])
def test_native_bridge_accepts_longest_representable_allele(segment):
    refdata = _engine.RefDataConfig("vdj")
    getattr(refdata, f"add_{segment}_allele")("X*01", "X", b"A" * NATIVE_MAX_ALLELE_LEN)
    assert getattr(refdata, f"{segment}_pool_size")() == 1


@pytest.mark.parametrize("segment", ["v", "d", "j"])
def test_native_bridge_rejects_oversized_allele_before_pool_mutation(segment):
    refdata = _engine.RefDataConfig("vdj")
    with pytest.raises(ValueError, match="65535"):
        getattr(refdata, f"add_{segment}_allele")("X*01", "X", b"A" * (NATIVE_MAX_ALLELE_LEN + 1))
    assert getattr(refdata, f"{segment}_pool_size")() == 0


def test_oversized_v_with_subregions_reports_length_limit_not_wrapped_length():
    refdata = _engine.RefDataConfig("vdj")
    with pytest.raises(ValueError, match="65535"):
        refdata.add_v_allele("X*01", "X", b"A" * (NATIVE_MAX_ALLELE_LEN + 1),
                             subregions=[("FWR1", 0, 10)])
    assert refdata.v_pool_size() == 0


def test_subregion_may_end_at_the_last_representable_base():
    refdata = _engine.RefDataConfig("vdj")
    refdata.add_v_allele("X*01", "X", b"A" * NATIVE_MAX_ALLELE_LEN,
                         subregions=[("FWR3", 100, NATIVE_MAX_ALLELE_LEN)])
    assert refdata.v_allele(0).subregions == [("FWR3", 100, NATIVE_MAX_ALLELE_LEN)]
