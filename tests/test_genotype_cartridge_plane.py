"""Cartridge genotype plane: PopulationGenotypeModel + Genotype.sample consumption."""
import pickle

import pytest

import GenAIRR as ga
import GenAIRR.data as gdata
from GenAIRR.genotype import Genotype
from GenAIRR.genotype_priors import PopulationGenotypeModel, PopulationNovelAllele


def _cfg():
    return gdata.HUMAN_IGH_OGRDB


def _vg_two_alleles(cfg):
    """A V gene with >= 2 alleles, plus its first two allele names."""
    vg = next(g for g, al in cfg.v_alleles.items() if len(al) >= 2)
    a0, a1 = (a.name for a in cfg.v_alleles[vg][:2])
    return vg, a0, a1


def test_model_requires_identity():
    m = PopulationGenotypeModel(model_id="", source="x")
    with pytest.raises(ValueError, match="model_id"):
        m.validate()
    m = PopulationGenotypeModel(model_id="x", source="")
    with pytest.raises(ValueError, match="source"):
        m.validate()


def test_model_validates_shapes():
    cfg = _cfg()
    vg, a0, a1 = _vg_two_alleles(cfg)
    ok = PopulationGenotypeModel(
        model_id="toy", source="unit-test",
        allele_frequencies={"V": {vg: {a0: 2.0, a1: 1.0}}},
        haplotype_deletion_prob={"V": {vg: 0.1}},
        chromosome_weights=(0.5, 0.5),
    )
    ok.validate(chain_type="vdj")  # no raise

    with pytest.raises(ValueError, match="finite"):
        PopulationGenotypeModel(model_id="t", source="s",
            allele_frequencies={"V": {vg: {a0: float("nan")}}}).validate()
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        PopulationGenotypeModel(model_id="t", source="s",
            haplotype_deletion_prob={"V": {vg: 1.5}}).validate()
    with pytest.raises(ValueError, match="non-negative"):
        PopulationGenotypeModel(model_id="t", source="s",
            chromosome_weights=(-1.0, 1.0)).validate()
    with pytest.raises(ValueError, match="at least one"):
        PopulationGenotypeModel(model_id="t", source="s",
            allele_frequencies={"V": {vg: {a0: 0.0, a1: 0.0}}}).validate()


def test_novel_shape_validation():
    base = PopulationNovelAllele(name="IGHV1-2*99", segment="V",
        base_allele="IGHV1-2*02", sequence="ACGT", frequency=1.0)
    PopulationGenotypeModel(model_id="t", source="s", novel_alleles=[base]).validate()

    with pytest.raises(ValueError, match="must be > 0"):
        PopulationGenotypeModel(model_id="t", source="s", novel_alleles=[
            PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHV1-2*02", sequence="ACGT", frequency=0.0)]).validate()
    with pytest.raises(ValueError, match="A/C/G/T|DNA"):
        PopulationGenotypeModel(model_id="t", source="s", novel_alleles=[
            PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHV1-2*02", sequence="ACGX", frequency=1.0)]).validate()
    with pytest.raises(ValueError, match="duplicate|unique"):
        PopulationGenotypeModel(model_id="t", source="s", novel_alleles=[
            PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHV1-2*02", sequence="ACGT", frequency=1.0),
            PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHV1-2*02", sequence="ACGA", frequency=1.0)]).validate()


def test_d_on_vj_rejected():
    m = PopulationGenotypeModel(model_id="t", source="s",
        allele_frequencies={"D": {"IGHD1-1": {"IGHD1-1*01": 1.0}}})
    with pytest.raises(ValueError, match="VJ|D-segment|D segment"):
        m.validate(chain_type="vj")


def test_d_on_vj_rejected_via_chaintype_object():
    class _FakeChainType:  # ChainType-like: exposes has_d
        has_d = False
    m = PopulationGenotypeModel(model_id="t", source="s",
        haplotype_deletion_prob={"D": {"IGHD1-1": 0.5}})
    with pytest.raises(ValueError, match="VJ|D-segment|D segment"):
        m.validate(chain_type=_FakeChainType())


def test_gene_keys_must_be_strings():
    with pytest.raises(ValueError, match="gene"):
        PopulationGenotypeModel(model_id="t", source="s",
            allele_frequencies={"V": {123: {"IGHV1-2*02": 1.0}}}).validate()
    with pytest.raises(ValueError, match="gene"):
        PopulationGenotypeModel(model_id="t", source="s",
            haplotype_deletion_prob={"V": {123: 0.1}}).validate()


def test_allow_nonfunctional_must_be_bool():
    with pytest.raises(ValueError, match="allow_nonfunctional"):
        PopulationGenotypeModel(model_id="t", source="s", novel_alleles=[
            PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHV1-2*02", sequence="ACGT", frequency=1.0,
                allow_nonfunctional="yes")]).validate()


def test_content_checksum_is_canonical():
    vg = "IGHV1-2"
    m1 = PopulationGenotypeModel(model_id="m", source="s",
        allele_frequencies={"V": {vg: {"IGHV1-2*02": 2, "IGHV1-2*04": 1}}},
        novel_alleles=[PopulationNovelAllele(name="IGHV1-2*99", segment="V",
            base_allele="IGHV1-2*02", sequence="acgt", frequency=1.0)])
    # same content, different dict insertion order + int-vs-float + DNA case
    m2 = PopulationGenotypeModel(model_id="m", source="s",
        allele_frequencies={"V": {vg: {"IGHV1-2*04": 1.0, "IGHV1-2*02": 2.0}}},
        novel_alleles=[PopulationNovelAllele(name="IGHV1-2*99", segment="V",
            base_allele="IGHV1-2*02", sequence="ACGT", frequency=1.0)])
    assert m1.content_checksum() == m2.content_checksum()

    m3 = PopulationGenotypeModel(model_id="m", source="s",
        allele_frequencies={"V": {vg: {"IGHV1-2*02": 3.0, "IGHV1-2*04": 1.0}}})
    assert m3.content_checksum() != m1.content_checksum()


def test_genotype_priors_field_default_and_checksum_invariant():
    cfg = _cfg()
    # default-new and bundled both have no plane
    assert getattr(cfg, "genotype_priors", "MISSING") is None
    base_checksum = cfg.compute_checksum()

    import copy
    cfg2 = copy.deepcopy(cfg)
    cfg2.genotype_priors = None  # explicitly None must not change the checksum
    assert cfg2.compute_checksum() == base_checksum

    cfg3 = copy.deepcopy(cfg)
    cfg3.genotype_priors = PopulationGenotypeModel(model_id="m", source="s")
    assert cfg3.compute_checksum() != base_checksum  # a real plane is cartridge identity


def test_genotype_priors_pickle_round_trip():
    import copy
    cfg = copy.deepcopy(_cfg())
    m = PopulationGenotypeModel(model_id="m", source="s",
        haplotype_deletion_prob={"V": {next(iter(cfg.v_alleles)): 0.2}})
    cfg.genotype_priors = m
    back = pickle.loads(pickle.dumps(cfg, protocol=4))
    assert back.genotype_priors.model_id == "m"
    assert back.genotype_priors.content_checksum() == m.content_checksum()


def test_manifest_genotype_priors_block():
    import copy
    cfg = copy.deepcopy(_cfg())
    block = cfg.cartridge_manifest()["models"]["genotype_priors"]
    assert block["available"] is False

    vg = next(iter(cfg.v_alleles))
    cfg.genotype_priors = PopulationGenotypeModel(
        model_id="m1", source="VDJbase-toy", version="1",
        allele_frequencies={"V": {vg: {cfg.v_alleles[vg][0].name: 1.0}}},
        haplotype_deletion_prob={"V": {vg: 0.1}},
        novel_alleles=[PopulationNovelAllele(name="IGHV1-2*99", segment="V",
            base_allele="IGHV1-2*02", sequence="ACGT", frequency=1.0)],
    )
    block = cfg.cartridge_manifest()["models"]["genotype_priors"]
    assert block["available"] is True
    assert block["model_id"] == "m1"
    assert block["source"] == "VDJbase-toy"
    assert block["model_checksum"] == cfg.genotype_priors.content_checksum()
    assert block["freq_gene_counts"]["V"] == 1
    assert block["deletion_gene_counts"]["V"] == 1
    assert block["novel_allele_count"] == 1
    assert block["chromosome_weights"] == [0.5, 0.5]
    assert block["source_field"] == "DataConfig.genotype_priors"


# ── Task 5: provenance scaffolding ───────────────────────────────


def test_prior_provenance_default_and_metadata():
    cfg = _cfg()
    vg, a0, _a1 = _vg_two_alleles(cfg)
    g = Genotype.from_dataconfig(cfg).homozygous(vg, a0).with_subject("S1")
    # builder-path genotype: every source non-cartridge, no model id
    prov = g.prior_provenance
    assert prov["allele_frequencies"] == "manual"
    assert prov["model_id"] is None
    md = g.to_metadata()
    assert md["subject_id"] == "S1"
    assert md["prior_provenance"] == prov
    assert "source_refdata_hash" in md


def test_prior_provenance_survives_snapshot():
    cfg = _cfg()
    g = Genotype.sample(cfg, seed=1)
    snap = g._snapshot()
    assert snap.prior_provenance == g.prior_provenance


# ── Task 6: sample plane consumption (catalogue alleles) ─────────


def _planed_cfg():
    import copy
    cfg = copy.deepcopy(_cfg())
    vg, a0, a1 = _vg_two_alleles(cfg)
    cfg.genotype_priors = PopulationGenotypeModel(
        model_id="m1", source="toy",
        allele_frequencies={"V": {vg: {a0: 100.0, a1: 1.0}}},
        haplotype_deletion_prob={"V": {vg: 0.0}},
        chromosome_weights=(0.5, 0.5),
    )
    return cfg, vg, a0, a1


def test_sample_auto_uses_plane_with_provenance():
    cfg, vg, a0, a1 = _planed_cfg()
    homo_a0 = sum(1 for s in range(60)
                  if Genotype.sample(cfg, seed=s).carried_alleles("V", vg) == {a0})
    assert homo_a0 > 40  # dominant plane allele -> usually homozygous-common
    g = Genotype.sample(cfg, seed=1)
    assert g.prior_provenance["allele_frequencies"] == "cartridge"
    assert g.prior_provenance["haplotype_deletion_prob"] == "cartridge"
    assert g.prior_provenance["chromosome_weights"] == "cartridge"
    assert g.prior_provenance["model_id"] == "m1"
    assert g.prior_provenance["model_checksum"] == cfg.genotype_priors.content_checksum()


def test_sample_opt_out_is_uniform():
    cfg, vg, a0, a1 = _planed_cfg()
    g = Genotype.sample(cfg, seed=1, use_cartridge_priors=False)
    assert g.prior_provenance["allele_frequencies"] == "uniform"
    assert g.prior_provenance["haplotype_deletion_prob"] == "default"
    assert g.prior_provenance["chromosome_weights"] == "default"
    assert g.prior_provenance["model_id"] is None
    # uniform: a1 should appear materially more than under the biased plane
    a1_seen = sum(1 for s in range(60)
                  if a1 in Genotype.sample(cfg, seed=s, use_cartridge_priors=False)
                  .carried_alleles("V", vg))
    assert a1_seen > 5


def test_sample_explicit_overrides_plane():
    cfg, vg, a0, a1 = _planed_cfg()
    g = Genotype.sample(cfg, seed=1, allele_frequencies={"V": {vg: {a1: 1.0}}})
    assert g.prior_provenance["allele_frequencies"] == "explicit"
    assert g.carried_alleles("V", vg) <= {a1}


def test_sample_mixed_sourcing():
    cfg, vg, a0, a1 = _planed_cfg()
    g = Genotype.sample(cfg, seed=2, allele_frequencies={"V": {vg: {a0: 1.0}}})
    # explicit freq, cartridge deletion, cartridge chromosome weights
    assert g.prior_provenance["allele_frequencies"] == "explicit"
    assert g.prior_provenance["haplotype_deletion_prob"] == "cartridge"
    assert g.prior_provenance["chromosome_weights"] == "cartridge"


def test_include_cartridge_novel_alleles_validated():
    cfg, vg, a0, a1 = _planed_cfg()
    with pytest.raises(ValueError, match="include_cartridge_novel_alleles"):
        Genotype.sample(cfg, seed=1, include_cartridge_novel_alleles="yes")


def test_sample_no_plane_unchanged():
    cfg = _cfg()  # no plane
    g = Genotype.sample(cfg, seed=3)
    assert g.prior_provenance["allele_frequencies"] == "uniform"
    assert g.prior_provenance["haplotype_deletion_prob"] == "default"


# ── Task 7: candidate-vs-carried novel injection ─────────────────


def _functional_novel_seq(cfg, base_name):
    """Find a single-base substitution off ``base_name`` that Genotype accepts
    as functional (avoids stop codons / anchor breakage), returning the seq."""
    base_allele = next(a for g in cfg.v_alleles.values() for a in g if a.name == base_name)
    base_seq = base_allele.ungapped_seq.upper()
    for pos in range(len(base_seq)):
        for nt in "ACGT":
            if nt == base_seq[pos]:
                continue
            cand = base_seq[:pos] + nt + base_seq[pos + 1:]
            try:
                (Genotype.from_dataconfig(cfg)
                 .add_novel_allele(f"{base_name.split('*')[0]}*97", base=base_name,
                                   sequence=cand, segment="V"))
                return cand
            except ValueError:
                continue
    raise AssertionError("no functional single-base novel found")


def _planed_cfg_with_novel(freq_novel=1000.0):
    import copy
    cfg = copy.deepcopy(_cfg())
    vg = next(g for g, al in cfg.v_alleles.items() if len(al) >= 1)
    base = cfg.v_alleles[vg][0].name
    novel_seq = _functional_novel_seq(cfg, base)
    novel = f"{vg}*97"
    cfg.genotype_priors = PopulationGenotypeModel(
        model_id="mN", source="toy",
        allele_frequencies={"V": {vg: {base: 1.0}}},
        haplotype_deletion_prob={"V": {vg: 0.0}},
        novel_alleles=[PopulationNovelAllele(name=novel, segment="V",
            base_allele=base, sequence=novel_seq, frequency=freq_novel)],
    )
    return cfg, vg, base, novel


def test_plane_novel_is_drawn_and_carried():
    cfg, vg, base, novel = _planed_cfg_with_novel()
    g = Genotype.sample(cfg, seed=1)  # auto -> cartridge freqs -> novels injected
    assert g.prior_provenance["novel_alleles"] == "cartridge"
    assert novel in g.carried_alleles("V", vg)  # dominant novel frequency
    eff = g.effective_dataconfig()
    assert any(a.name == novel for a in eff.v_alleles[vg])


def test_uncarried_plane_novel_absent_from_export():
    # rare novel: base dominates -> novel essentially never carried
    cfg, vg, base, novel = _planed_cfg_with_novel(freq_novel=1e-9)
    g = Genotype.sample(cfg, seed=7)
    assert novel not in g.carried_alleles("V", vg)
    eff = g.effective_dataconfig()
    assert all(a.name != novel for a in eff.v_alleles[vg])


def test_novel_skipped_with_explicit_freqs_auto():
    cfg, vg, base, novel = _planed_cfg_with_novel()
    g = Genotype.sample(cfg, seed=1, allele_frequencies={"V": {vg: {base: 1.0}}})
    assert g.prior_provenance["novel_alleles"] == "none"
    assert novel not in g.carried_alleles("V", vg)


def test_novel_forced_with_explicit_freqs_true():
    cfg, vg, base, novel = _planed_cfg_with_novel()
    g = Genotype.sample(cfg, seed=1, allele_frequencies={"V": {vg: {base: 1.0}}},
                        include_cartridge_novel_alleles=True)
    assert g.prior_provenance["novel_alleles"] == "cartridge"
    # synthesized table: explicit base at 1.0 + novel at huge frequency -> novel wins
    assert novel in g.carried_alleles("V", vg)


def test_novel_disabled():
    cfg, vg, base, novel = _planed_cfg_with_novel()
    g = Genotype.sample(cfg, seed=1, include_cartridge_novel_alleles=False)
    assert g.prior_provenance["novel_alleles"] == "none"
    assert novel not in g.carried_alleles("V", vg)
