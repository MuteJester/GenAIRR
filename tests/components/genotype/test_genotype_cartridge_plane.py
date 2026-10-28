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
    # integers must NOT slip through (1 == True / 0 == False in Python)
    for bad in (0, 1, 2):
        with pytest.raises(ValueError, match="include_cartridge_novel_alleles"):
            Genotype.sample(cfg, seed=1, include_cartridge_novel_alleles=bad)


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


# ── Task 8: from_genotypes pure estimator ────────────────────────


def test_from_genotypes_counts_chromosomes_and_deletions():
    cfg = _cfg()
    vg, a0, a1 = _vg_two_alleles(cfg)
    # subject A: homozygous a0 (2 chromosomes of a0)
    gA = (Genotype.from_dataconfig(cfg).homozygous(vg, a0)
          .complete_from_reference().with_subject("A"))
    # subject B: heterozygous a0/a1 (1 each)
    gB = (Genotype.from_dataconfig(cfg).heterozygous(vg, a0, a1)
          .complete_from_reference().with_subject("B"))
    # subject C: a0 on hap0, deleted on hap1 (hemizygous) -> 1 deleted haplotype
    gC = (Genotype.from_dataconfig(cfg).homozygous(vg, a0)
          .delete_gene(vg, haplotype=1).complete_from_reference().with_subject("C"))

    m = PopulationGenotypeModel.from_genotypes([gA, gB, gC], cfg=cfg,
                                               model_id="est", source="unit")
    # a0 chromosomes: A=2, B=1, C=1 => 4 ; a1: B=1 => 1
    assert m.allele_frequencies["V"][vg][a0] == 4.0
    assert m.allele_frequencies["V"][vg][a1] == 1.0
    # deletions for vg: only C hap1 => 1 / (2*3) subjects
    assert m.haplotype_deletion_prob["V"][vg] == pytest.approx(1.0 / 6.0)
    m.validate(chain_type="vdj")  # estimator output is shape-valid


def test_from_genotypes_duplicate_subject_id_rejected():
    cfg = _cfg()
    g1 = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("X")
    g2 = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("X")
    with pytest.raises(ValueError, match="duplicate subject"):
        PopulationGenotypeModel.from_genotypes([g1, g2], cfg=cfg,
                                               model_id="e", source="u")


def test_from_genotypes_min_subjects():
    cfg = _cfg()
    g1 = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("X")
    with pytest.raises(ValueError, match="min_subjects"):
        PopulationGenotypeModel.from_genotypes([g1], cfg=cfg, min_subjects=5,
                                               model_id="e", source="u")


def test_from_genotypes_rejects_duplicated_gene():
    cfg = _cfg()
    vg, a0, a1 = _vg_two_alleles(cfg)
    g = (Genotype.from_dataconfig(cfg).homozygous(vg, a0)
         .duplicate_gene(vg, [a0, a1], haplotype=0)  # copy-number > 1 on a slot
         .complete_from_reference().with_subject("D"))
    with pytest.raises(ValueError, match="copy-number|duplicated gene"):
        PopulationGenotypeModel.from_genotypes([g], cfg=cfg, model_id="e", source="u")


def test_from_genotypes_pseudocount_scope():
    cfg = _cfg()
    vg, a0, a1 = _vg_two_alleles(cfg)
    gA = (Genotype.from_dataconfig(cfg).homozygous(vg, a0)
          .complete_from_reference().with_subject("A"))
    m = PopulationGenotypeModel.from_genotypes([gA], cfg=cfg, pseudocount=0.5,
                                               subject_id_policy="allow_duplicates",
                                               model_id="e", source="u")
    # a0 observed twice + 0.5 ; a1 unobserved but catalogue -> 0.5 (pseudocount only)
    assert m.allele_frequencies["V"][vg][a0] == pytest.approx(2.5)
    assert m.allele_frequencies["V"][vg][a1] == pytest.approx(0.5)


# ── Task 9: builder set/estimate_genotype_priors ─────────────────


def _builder_from_bundled():
    """A ReferenceCartridgeBuilder seeded from the bundled IGH catalogue."""
    from GenAIRR.cartridge_builder import ReferenceCartridgeBuilder
    from GenAIRR.dataconfig.enums import ChainType
    cfg = _cfg()
    b = ReferenceCartridgeBuilder(ChainType.BCR_HEAVY)
    b._v_alleles = {g: list(a) for g, a in cfg.v_alleles.items()}
    b._d_alleles = {g: list(a) for g, a in (cfg.d_alleles or {}).items()}
    b._j_alleles = {g: list(a) for g, a in cfg.j_alleles.items()}
    b._metadata = cfg.metadata
    return b, cfg


def test_set_genotype_priors_catalogue_aware():
    b, cfg = _builder_from_bundled()
    vg = next(iter(cfg.v_alleles))
    good = PopulationGenotypeModel(model_id="m", source="s",
        allele_frequencies={"V": {vg: {cfg.v_alleles[vg][0].name: 1.0}}})
    assert b.set_genotype_priors(good) is b  # chainable
    assert b._genotype_priors is good

    with pytest.raises(ValueError, match="not in the cartridge|unknown"):
        b.set_genotype_priors(PopulationGenotypeModel(model_id="m", source="s",
            allele_frequencies={"V": {"NOPE-GENE": {"NOPE*01": 1.0}}}))
    with pytest.raises(ValueError, match="not a known allele|allele"):
        b.set_genotype_priors(PopulationGenotypeModel(model_id="m", source="s",
            allele_frequencies={"V": {vg: {"IGHVNOPE*99": 1.0}}}))


def test_set_genotype_priors_validates_novel():
    b, cfg = _builder_from_bundled()
    with pytest.raises(ValueError, match="base allele|not found"):
        b.set_genotype_priors(PopulationGenotypeModel(model_id="m", source="s",
            novel_alleles=[PopulationNovelAllele(name="IGHV1-2*99", segment="V",
                base_allele="IGHVNOPE*01", sequence="ACGT", frequency=1.0)]))


def test_estimate_genotype_priors_chainable_and_attaches():
    b, cfg = _builder_from_bundled()
    vg, a0, a1 = _vg_two_alleles(cfg)
    gA = Genotype.from_dataconfig(cfg).homozygous(vg, a0).complete_from_reference().with_subject("A")
    gB = Genotype.from_dataconfig(cfg).heterozygous(vg, a0, a1).complete_from_reference().with_subject("B")
    out = b.estimate_genotype_priors([gA, gB], model_id="est", source="cohort")
    assert out is b
    assert b._genotype_priors.allele_frequencies["V"][vg][a0] == 3.0


def test_from_genotypes_requires_complete_genotypes():
    # A genotype missing a catalogue gene (no complete_from_reference) must NOT
    # be silently scored as a double deletion — it should raise.
    cfg = _cfg()
    vg, a0, _a1 = _vg_two_alleles(cfg)
    g = Genotype.from_dataconfig(cfg).homozygous(vg, a0).with_subject("P")  # partial
    with pytest.raises(ValueError, match="does not specify|complete_from_reference"):
        PopulationGenotypeModel.from_genotypes([g], cfg=cfg, model_id="x", source="y")


def test_from_genotypes_all_none_subject_ids_rejected_under_require_unique():
    cfg = _cfg()
    g1 = Genotype.from_dataconfig(cfg).complete_from_reference()  # no subject id
    g2 = Genotype.from_dataconfig(cfg).complete_from_reference()
    with pytest.raises(ValueError, match="subject_id"):
        PopulationGenotypeModel.from_genotypes([g1, g2], cfg=cfg,
                                               model_id="x", source="y")


def test_estimate_with_novel_round_trips_and_samples():
    cfg, vg, base, novel = _planed_cfg_with_novel()  # gives us a functional novel seq
    novel_seq = next(a for a in cfg.genotype_priors.novel_alleles).sequence
    base_cfg = _cfg()
    gn = (Genotype.from_dataconfig(base_cfg)
          .add_novel_allele(novel, base=base, sequence=novel_seq, segment="V")
          .homozygous(vg, novel).complete_from_reference().with_subject("N"))
    from GenAIRR.cartridge_builder import ReferenceCartridgeBuilder
    from GenAIRR.dataconfig.enums import ChainType
    b = ReferenceCartridgeBuilder(ChainType.BCR_HEAVY)
    b._v_alleles = {g: list(a) for g, a in base_cfg.v_alleles.items()}
    b._d_alleles = {g: list(a) for g, a in (base_cfg.d_alleles or {}).items()}
    b._j_alleles = {g: list(a) for g, a in base_cfg.j_alleles.items()}
    b._metadata = base_cfg.metadata
    # round-trip must not raise (novel lives in novel_alleles, not allele_frequencies)
    b.estimate_genotype_priors([gn], model_id="e", source="c",
                               subject_id_policy="allow_duplicates")
    model = b._genotype_priors
    assert any(nv.name == novel for nv in model.novel_alleles)
    assert novel not in model.allele_frequencies.get("V", {}).get(vg, {})


def test_from_genotypes_pseudocount_validation():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("A")
    for bad in (-1.0, float("nan"), float("inf"), True):
        with pytest.raises(ValueError, match="pseudocount"):
            PopulationGenotypeModel.from_genotypes([g], cfg=cfg, pseudocount=bad,
                                                   model_id="x", source="y")


def test_from_genotypes_rejects_unknown_carried_allele():
    cfg = _cfg()
    vg, a0, _a1 = _vg_two_alleles(cfg)
    g = Genotype.from_dataconfig(cfg).homozygous(vg, a0).complete_from_reference().with_subject("A")
    # forge an unregistered, non-catalogue allele directly into a slot
    g._slots["V"][vg] = [[("BOGUS*01", 1, 1.0)], [(a0, 1, 1.0)]]
    with pytest.raises(ValueError, match="neither a catalogue allele nor a registered novel"):
        PopulationGenotypeModel.from_genotypes([g], cfg=cfg, model_id="x", source="y")


def test_from_genotypes_rejects_copy_count_gt_one():
    cfg = _cfg()
    vg, a0, _a1 = _vg_two_alleles(cfg)
    g = Genotype.from_dataconfig(cfg).homozygous(vg, a0).complete_from_reference().with_subject("A")
    g._slots["V"][vg] = [[(a0, 2, 1.0)], [(a0, 1, 1.0)]]  # single entry, copies=2
    with pytest.raises(ValueError, match="copy-number|duplicated"):
        PopulationGenotypeModel.from_genotypes([g], cfg=cfg, model_id="x", source="y")


def test_from_genotypes_validates_identity_before_return():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("A")
    with pytest.raises(ValueError, match="model_id"):
        PopulationGenotypeModel.from_genotypes([g], cfg=cfg, model_id="", source="y")


def test_sample_validates_directly_attached_plane():
    import copy
    cfg = copy.deepcopy(_cfg())
    vg, a0, _a1 = _vg_two_alleles(cfg)
    # invalid plane attached directly (bypassing the builder)
    cfg.genotype_priors = PopulationGenotypeModel(model_id="", source="",
        allele_frequencies={"V": {vg: {a0: 1.0}}})
    with pytest.raises(ValueError, match="model_id|source"):
        Genotype.sample(cfg, seed=1)


def test_manifest_invalid_plane_marked_not_valid():
    import copy, math
    cfg = copy.deepcopy(_cfg())
    cfg.genotype_priors = PopulationGenotypeModel(model_id="m", source="s",
        chromosome_weights=(float("nan"), 1.0))
    block = cfg.cartridge_manifest()["models"]["genotype_priors"]
    assert block["available"] is True
    assert block["valid"] is False
    # no non-JSON-clean (NaN) numerics leak through
    cw = block["chromosome_weights"]
    assert cw is None or all(math.isfinite(x) for x in cw)


def test_to_metadata_effective_refdata_hash_for_novel():
    cfg, vg, base, novel = _planed_cfg_with_novel()
    g = Genotype.sample(cfg, seed=1)
    assert novel in g.carried_alleles("V", vg)
    md = g.to_metadata()
    assert "effective_refdata_hash" in md
    assert md["effective_refdata_hash"] != md["source_refdata_hash"]
    # a no-novel genotype: effective == source
    g2 = Genotype.from_dataconfig(_cfg()).complete_from_reference()
    md2 = g2.to_metadata()
    assert md2["effective_refdata_hash"] == md2["source_refdata_hash"]


def test_manifest_never_crashes_on_garbage_plane():
    import copy, math
    cfg = copy.deepcopy(_cfg())
    # non-finite frequency weight + bad chromosome weight type -> content_checksum
    # and float() would raise; the manifest must stay JSON-clean and not crash.
    cfg.genotype_priors = PopulationGenotypeModel(model_id="m", source="s",
        chromosome_weights=("x", 1.0))
    block = cfg.cartridge_manifest()["models"]["genotype_priors"]
    assert block["available"] is True and block["valid"] is False
    assert block["chromosome_weights"] is None
    assert block["model_checksum"] is None or isinstance(block["model_checksum"], str)

    # genotype_priors set to a non-model object must not crash the manifest either
    cfg.genotype_priors = object()
    block = cfg.cartridge_manifest()["models"]["genotype_priors"]
    assert block["available"] is True and block["valid"] is False


def test_use_cartridge_priors_must_be_bool():
    cfg, vg, a0, a1 = _planed_cfg()
    for bad in ("False", 1, 0, None):
        with pytest.raises(ValueError, match="use_cartridge_priors"):
            Genotype.sample(cfg, seed=1, use_cartridge_priors=bad)


def test_from_genotypes_include_novel_must_be_bool():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("A")
    for bad in ("yes", 1, None):
        with pytest.raises(ValueError, match="include_novel"):
            PopulationGenotypeModel.from_genotypes([g], cfg=cfg, include_novel=bad,
                                                   model_id="x", source="y")


def test_from_genotypes_min_subjects_and_segments_validation():
    cfg = _cfg()
    g = Genotype.from_dataconfig(cfg).complete_from_reference().with_subject("A")
    for bad in (-1, 0, True, "2"):
        with pytest.raises(ValueError, match="min_subjects"):
            PopulationGenotypeModel.from_genotypes([g], cfg=cfg, min_subjects=bad,
                                                   model_id="x", source="y")
    with pytest.raises(ValueError, match="segment"):
        PopulationGenotypeModel.from_genotypes([g], cfg=cfg, segments=("X",),
                                               model_id="x", source="y")
    with pytest.raises(ValueError, match="segment"):
        PopulationGenotypeModel.from_genotypes([g], cfg=cfg, segments="V",
                                               model_id="x", source="y")
    with pytest.raises(ValueError, match="segment"):
        PopulationGenotypeModel.from_genotypes([g], cfg=cfg, segments=("V", "V"),
                                               model_id="x", source="y")


def test_sample_draw_independent_of_freq_dict_order():
    import copy
    cfg = _cfg()
    vg, a0, a1 = _vg_two_alleles(cfg)
    c1 = copy.deepcopy(cfg)
    c1.genotype_priors = PopulationGenotypeModel(model_id="m", source="s",
        allele_frequencies={"V": {vg: {a0: 2.0, a1: 3.0}}})
    c2 = copy.deepcopy(cfg)
    c2.genotype_priors = PopulationGenotypeModel(model_id="m", source="s",
        allele_frequencies={"V": {vg: {a1: 3.0, a0: 2.0}}})  # reversed insertion order
    assert c1.genotype_priors.content_checksum() == c2.genotype_priors.content_checksum()
    for s in range(15):
        assert (Genotype.sample(c1, seed=s).to_table()
                == Genotype.sample(c2, seed=s).to_table())


# ── Task 10: end-to-end ──────────────────────────────────────────


def test_end_to_end_planed_cartridge_truth_calls_carried():
    cfg, vg, a0, a1 = _planed_cfg()
    g = Genotype.sample(cfg, seed=11)
    assert g.prior_provenance["allele_frequencies"] == "cartridge"
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(
        n=200, seed=2, expose_provenance=True)
    assert len(res) == 200
    for r in res:
        for seg, col in (("V", "truth_v_call"), ("D", "truth_d_call"), ("J", "truth_j_call")):
            call = r[col]
            if not call:
                continue
            gene = call.split("*")[0]
            assert call in g.carried_alleles(seg, gene), (seg, call)


def test_end_to_end_planed_novel_flows_to_records():
    cfg, vg, base, novel = _planed_cfg_with_novel()
    g = Genotype.sample(cfg, seed=1)
    assert novel in g.carried_alleles("V", vg)
    res = ga.Experiment.on(cfg).with_genotype(g).recombine().run_records(
        n=50, seed=3, expose_provenance=True)
    # the novel can legitimately appear as a V truth call
    assert any(r["truth_v_call"] == novel for r in res)
