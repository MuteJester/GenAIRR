"""T2-8 Phase 2 + 3a Python tests for the multi-format reference loaders.

Covers:
  * ReferenceLoader.open_plain / open_airrc / open_ogrdb / open_igblast / open_auto
  * format auto-detection rules (.json, '|'-headers vs plain, sibling JSON)
  * RandomDataConfigBuilder.make_from_reference end-to-end
  * BuildReport accumulation (accepted + rejected + summary)
"""
from __future__ import annotations

import json
import textwrap

import pytest

from GenAIRR._native._anchor import (
    AnchorConfidence,
    Locus,
    ReferenceLoader,
    Segment,
)


# ── Plain FASTA loader ─────────────────────────────────────────────


class TestPlainFastaLoader:

    def test_basic_round_trip(self, tmp_path):
        f = tmp_path / "ref.fa"
        f.write_text(">IGHV1-2*01\nacgtacgtacgt\n>IGHJ4*02\ntggggcaaaggc\n")
        with ReferenceLoader.open_plain(str(f), locus_hint=Locus.IGH) as L:
            recs = list(L)
        assert len(recs) == 2
        assert recs[0].name == "IGHV1-2*01"
        assert recs[0].source == "plain-fasta"
        assert recs[0].gap_convention_imgt is False
        assert recs[0].gapped_sequence is None      # plain has no gapped form
        assert recs[1].name == "IGHJ4*02"

    def test_locus_hint_used_when_name_unknown(self, tmp_path):
        f = tmp_path / "novel.fa"
        f.write_text(">novel_X\nacgtacgt\n")
        with ReferenceLoader.open_plain(str(f), locus_hint=Locus.TRG,
                                         segment_hint=Segment.V) as L:
            rec = next(iter(L))
        assert rec.locus is Locus.TRG
        assert rec.segment is Segment.V


# ── AIRR-C JSON loader ─────────────────────────────────────────────


_AIRRC_SAMPLE = {
    "GermlineSet": [{
        "germline_set_name": "test",
        "species": {"label": "Homo sapiens"},
        "allele_descriptions": [
            {
                "label": "IGHV1-2*02",
                "locus": "IGH",
                "sequence": "acgtacgtacgttgtaaaaaaaaaaaaaaaaa",
                "functional": True,
                "v_gene_delineations": [
                    {"delineation_scheme": "IMGT", "cdr3_start": 13},
                ],
            },
            {
                "label": "IGHJ4*01",
                "locus": "IGH",
                "sequence": "tggggcaaaggcaaaaaaaa",
                "functional": True,
                "j_cdr3_end": 3,
            },
        ],
    }]
}


class TestAirrcLoader:

    def test_round_trip_with_explicit_anchors(self, tmp_path):
        f = tmp_path / "germline.json"
        f.write_text(json.dumps(_AIRRC_SAMPLE))
        with ReferenceLoader.open_airrc(str(f)) as L:
            recs = list(L)
        assert len(recs) == 2
        v, j = recs
        assert v.name == "IGHV1-2*02"
        assert v.locus is Locus.IGH
        assert v.segment is Segment.V
        assert v.species == "Homo sapiens"
        assert v.source == "airrc-germline-set"
        # cdr3_start (1-based) = 13 → explicit_anchor (0-based) = 12
        assert v.explicit_anchor == 12
        # j_cdr3_end (1-based) = 3 → conserved-codon start (0-based) = 0
        assert j.explicit_anchor == 0


# ── OGRDB combined loader ──────────────────────────────────────────


class TestOgrdbLoader:

    def test_with_sidecar_pulls_explicit_anchor(self, tmp_path):
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">IGHV1-2*02\nacgtacgtacgttgtaaa\n")
        sidecar = tmp_path / "ref.json"
        sidecar.write_text(json.dumps({
            "GermlineSet": [{
                "allele_descriptions": [{
                    "label": "IGHV1-2*02",
                    "functional": True,
                    "v_gene_delineations": [
                        {"delineation_scheme": "IMGT", "cdr3_start": 13},
                    ],
                }],
            }],
        }))
        with ReferenceLoader.open_ogrdb(str(fasta), str(sidecar),
                                          locus_hint=Locus.IGH,
                                          segment_hint=Segment.V) as L:
            recs = list(L)
        assert len(recs) == 1
        assert recs[0].name == "IGHV1-2*02"
        assert recs[0].explicit_anchor == 12
        assert recs[0].source == "ogrdb"

    def test_without_sidecar_falls_back(self, tmp_path):
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">IGHV1-2*02\nacgtacgt\n")
        with ReferenceLoader.open_ogrdb(str(fasta), None,
                                          locus_hint=Locus.IGH) as L:
            rec = next(iter(L))
        # No sidecar → no explicit anchor; resolver will fall back
        # to motif on the resulting record.
        assert rec.explicit_anchor == -1


# ── IgBLAST bundle loader ──────────────────────────────────────────


class TestIgBlastLoader:

    def test_aux_j_anchor_wired(self, tmp_path):
        fasta = tmp_path / "j.fasta"
        fasta.write_text(">IGHJ4*01\ntggggcaaaggc\n")
        # cdr3_end is 0-based; for "tggggcaaaggc" the conserved Trp ends
        # at position 2 → cdr3_end=2 → anchor = 0.
        aux = tmp_path / "j.aux"
        aux.write_text("IGHJ4*01\t0\tVH\t2\n")
        with ReferenceLoader.open_igblast(str(fasta), str(aux),
                                            locus_hint=Locus.IGH,
                                            segment_hint=Segment.J) as L:
            rec = next(iter(L))
        assert rec.explicit_anchor == 0
        assert rec.source == "igblast-bundle"


# ── Auto-dispatcher ────────────────────────────────────────────────


class TestAutoDispatch:

    def test_json_extension_routes_to_airrc(self, tmp_path):
        f = tmp_path / "g.json"
        f.write_text("[]")
        with ReferenceLoader.open_auto(str(f)) as L:
            recs = list(L)
        assert recs == []  # empty array → no records but no error

    def test_pipe_header_routes_to_imgt(self, tmp_path):
        f = tmp_path / "ig.fasta"
        f.write_text(">X1|IGHV1-2*02|Homo sapiens|F|V-REGION|...\nacgtacgt\n")
        with ReferenceLoader.open_auto(str(f)) as L:
            rec = next(iter(L))
        assert rec.source == "imgt-vquest"

    def test_sibling_json_routes_to_ogrdb(self, tmp_path):
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">IGHV1-2*02\nacgtacgt\n")
        sidecar = tmp_path / "ref.json"
        sidecar.write_text(json.dumps({"GermlineSet": [{
            "allele_descriptions": [{"label": "IGHV1-2*02", "functional": True}]
        }]}))
        with ReferenceLoader.open_auto(str(fasta), locus_hint=Locus.IGH) as L:
            rec = next(iter(L))
        assert rec.source == "ogrdb"

    def test_plain_fasta_fallback(self, tmp_path):
        f = tmp_path / "plain.fa"
        f.write_text(">novel_X\nacgtacgt\n")
        with ReferenceLoader.open_auto(str(f), locus_hint=Locus.IGH) as L:
            rec = next(iter(L))
        assert rec.source == "plain-fasta"


# ── make_from_reference end-to-end ─────────────────────────────────


class TestMakeFromReference:

    def test_imgt_format_explicit(self, tmp_path):
        """Build a tiny but valid IMGT-style DataConfig from synthetic
        FASTAs."""
        from GenAIRR.dataconfig.make.random import RandomDataConfigBuilder
        from GenAIRR.dataconfig.enums import Species, ChainType

        # IMGT-shaped V FASTA: gapped sequence with Cys at IMGT pos 309.
        v_fasta = tmp_path / "v.fasta"
        v_seq = "a" * 309 + "tgt" + "a" * 30      # ungapped 342 nt
        v_fasta.write_text(
            f">X1|IGHV1-1*01|Homo sapiens|F|V-REGION|...\n{v_seq}\n"
        )
        # J FASTA with WGXG motif at start.
        j_fasta = tmp_path / "j.fasta"
        j_fasta.write_text(
            ">X2|IGHJ1*01|Homo sapiens|F|J-REGION|...\n"
            "tggggcaaaggcacc\n"
        )
        # Minimal D for IGH.
        d_fasta = tmp_path / "d.fasta"
        d_fasta.write_text(">X3|IGHD1*01|Homo sapiens|F|D-REGION|...\nacgtacgt\n")

        builder = RandomDataConfigBuilder(
            species=Species.HUMAN,
            chain_type=ChainType.BCR_HEAVY,
        )
        cfg, report = builder.make_from_reference(
            v_reference=str(v_fasta),
            j_reference=str(j_fasta),
            d_reference=str(d_fasta),
            v_format="imgt",
            j_format="imgt",
        )
        assert cfg is not None
        assert report.accepted >= 3       # 1 V + 1 D + 1 J = 3
        assert report.total_seen >= 3
        # V allele anchor wired through.
        v_alleles = list(cfg.v_alleles.values())[0]
        assert v_alleles[0].anchor == 309
        assert "summary" in report.summary().lower() or "accepted" in report.summary().lower()

    def test_airrc_format(self, tmp_path):
        """Build from AIRR-C JSON references."""
        from GenAIRR.dataconfig.make.random import RandomDataConfigBuilder
        from GenAIRR.dataconfig.enums import Species, ChainType

        v_seq = "a" * 12 + "tgt" + "a" * 30      # cdr3_start=13 (1-based)
        v_path = tmp_path / "v.json"
        v_path.write_text(json.dumps({
            "GermlineSet": [{
                "species": {"label": "Homo sapiens"},
                "allele_descriptions": [{
                    "label": "IGHV1-1*01",
                    "locus": "IGH",
                    "sequence": v_seq,
                    "functional": True,
                    "v_gene_delineations": [
                        {"delineation_scheme": "IMGT", "cdr3_start": 13},
                    ],
                }],
            }],
        }))
        j_path = tmp_path / "j.json"
        j_path.write_text(json.dumps({
            "GermlineSet": [{
                "allele_descriptions": [{
                    "label": "IGHJ1*01",
                    "locus": "IGH",
                    "sequence": "tggggcaaaggcaaa",
                    "functional": True,
                    "j_cdr3_end": 3,
                }],
            }],
        }))
        d_path = tmp_path / "d.json"
        d_path.write_text(json.dumps({
            "GermlineSet": [{
                "allele_descriptions": [{
                    "label": "IGHD1*01",
                    "locus": "IGH",
                    "sequence": "acgtacgt",
                    "functional": True,
                }],
            }],
        }))
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN,
            chain_type=ChainType.BCR_HEAVY,
        )
        cfg, report = builder.make_from_reference(
            v_reference=str(v_path),
            j_reference=str(j_path),
            d_reference=str(d_path),
            v_format="airrc", j_format="airrc",
        )
        assert cfg is not None
        v_alleles = list(cfg.v_alleles.values())[0]
        assert v_alleles[0].anchor == 12       # 1-based 13 → 0-based 12
        assert report.accepted >= 3

    def test_empty_v_raises(self, tmp_path):
        from GenAIRR.dataconfig.make.random import RandomDataConfigBuilder
        from GenAIRR.dataconfig.enums import Species, ChainType
        empty = tmp_path / "empty.fasta"
        empty.write_text("")
        j = tmp_path / "j.fasta"
        j.write_text(">IGHJ1*01\ntggggcaaaggc\n")
        builder = RandomDataConfigBuilder(
            species=Species.HUMAN,
            chain_type=ChainType.BCR_HEAVY,
        )
        with pytest.raises(ValueError, match="empty"):
            builder.make_from_reference(
                v_reference=str(empty), j_reference=str(j),
                v_format="plain", j_format="plain",
                v_locus_hint=Locus.IGH, j_locus_hint=Locus.IGH,
            )
