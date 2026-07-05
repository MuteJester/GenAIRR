#!/usr/bin/env python3
"""Fetch IMGT V-QUEST germline FASTA and build GenAIRR cartridges.

Maintainer tool (not part of the shipped ``GenAIRR`` wheel). It downloads
the IMGT V-QUEST reference directory FASTA for each species/locus and
builds a structural :class:`GenAIRR.DataConfig` cartridge per
species+chain via :class:`GenAIRR.ReferenceCartridgeBuilder`.

**Structural output.** Cartridges carry the germline V/D/J allele
catalogue plus IMGT-derived V-subregion (FWR/CDR) annotations. They do
**not** carry data-derived empirical distributions (trim lengths, NP
lengths, NP base model, gene usage): those are not inferable from
germline FASTA alone. At simulation time the engine falls back to its
uniform defaults for those parameters — i.e. a uniform representation of
the scenarios a repertoire can express. Fit real distributions into the
returned ``DataConfig`` (or estimate them with the builder's
``estimate_*`` methods from your own AIRR data) if you need
empirically-grounded parameters.

Only the bundled ``HUMAN_IGH`` / ``HUMAN_IGK`` / ``HUMAN_IGL`` /
``HUMAN_TCRB`` cartridges ship with real data-derived distributions;
everything this tool builds is uniform-by-construction.

Usage::

    # Build one species into a scratch dir (never touches the bundled set):
    python tools/build_imgt_configs.py --output-dir ./built_cartridges \\
        --species Mus_musculus

    # Probe availability only, no build:
    python tools/build_imgt_configs.py --output-dir ./built_cartridges --dry-run

``--output-dir`` is required: writing a freshly-built (structural)
cartridge over the shipped ``src/GenAIRR/data/builtin_dataconfigs/`` set
would change simulation behaviour and break the golden tests, so it is
never the default. Point ``--output-dir`` there explicitly only when you
intend to regenerate the bundled set.
"""
from __future__ import annotations

import argparse
import logging
import os
import pickle
import sys
import time
import urllib.error
import urllib.request

from GenAIRR import DataConfig, ReferenceCartridgeBuilder
from GenAIRR.dataconfig.enums import ChainType, Species

logging.basicConfig(level=logging.INFO, format="%(levelname)-8s %(message)s")
logger = logging.getLogger("build_imgt_configs")

# ─────────────────────────────────────────────────────────────
# IMGT V-QUEST reference directory layout
# ─────────────────────────────────────────────────────────────

IMGT_BASE = (
    "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory"
)

# IMGT V-QUEST species directories carried by GenAIRR's bundled set.
IMGT_SPECIES = [
    "Aotus_nancymaae",
    "Bos_taurus",
    "Camelus_dromedarius",
    "Canis_lupus_familiaris",
    "Capra_hircus",
    "Danio_rerio",
    "Equus_caballus",
    "Felis_catus",
    "Gallus_gallus",
    "Gorilla_gorilla_gorilla",
    "Homo_sapiens",
    "Macaca_fascicularis",
    "Macaca_mulatta",
    "Mus_musculus",
    "Mus_musculus_C57BL6J",
    "Mustela_putorius_furo",
    "Oncorhynchus_mykiss",
    "Ornithorhynchus_anatinus",
    "Oryctolagus_cuniculus",
    "Ovis_aries",
    "Rattus_norvegicus",
    "Salmo_salar",
    "Sus_scrofa",
    "Vicugna_pacos",
]

# IMGT directory name → (Species enum, short label used in the cartridge name).
SPECIES_MAP = {
    "Homo_sapiens":             (Species.HUMAN,              "HUMAN"),
    "Mus_musculus":             (Species.MOUSE,              "MOUSE"),
    "Mus_musculus_C57BL6J":     (Species.MOUSE_C57BL6J,      "MOUSE_C57BL6J"),
    "Rattus_norvegicus":        (Species.RAT,                "RAT"),
    "Oryctolagus_cuniculus":    (Species.RABBIT,             "RABBIT"),
    "Macaca_mulatta":           (Species.RHESUS_MACAQUE,     "RHESUS"),
    "Macaca_fascicularis":      (Species.CYNOMOLGUS_MACAQUE, "CYNOMOLGUS"),
    "Gorilla_gorilla_gorilla":  (Species.GORILLA,            "GORILLA"),
    "Aotus_nancymaae":          (Species.OWL_MONKEY,         "OWL_MONKEY"),
    "Bos_taurus":               (Species.COW,                "COW"),
    "Ovis_aries":               (Species.SHEEP,              "SHEEP"),
    "Capra_hircus":             (Species.GOAT,               "GOAT"),
    "Sus_scrofa":               (Species.PIG,                "PIG"),
    "Equus_caballus":           (Species.HORSE,              "HORSE"),
    "Canis_lupus_familiaris":   (Species.DOG,                "DOG"),
    "Felis_catus":              (Species.CAT,                "CAT"),
    "Camelus_dromedarius":      (Species.DROMEDARY_CAMEL,    "DROMEDARY"),
    "Vicugna_pacos":            (Species.ALPACA,             "ALPACA"),
    "Mustela_putorius_furo":    (Species.FERRET,             "FERRET"),
    "Ornithorhynchus_anatinus": (Species.PLATYPUS,           "PLATYPUS"),
    "Gallus_gallus":            (Species.CHICKEN,            "CHICKEN"),
    "Danio_rerio":              (Species.ZEBRAFISH,          "ZEBRAFISH"),
    "Oncorhynchus_mykiss":      (Species.TROUT,              "TROUT"),
    "Salmo_salar":              (Species.SALMON,             "SALMON"),
}

# IMGT locus → (ChainType, chain label, V file, D file or None, J file).
LOCUS_DEFS = {
    "IGH": (ChainType.BCR_HEAVY,        "IGH",  "IGHV.fasta", "IGHD.fasta", "IGHJ.fasta"),
    "IGK": (ChainType.BCR_LIGHT_KAPPA,  "IGK",  "IGKV.fasta", None,          "IGKJ.fasta"),
    "IGL": (ChainType.BCR_LIGHT_LAMBDA, "IGL",  "IGLV.fasta", None,          "IGLJ.fasta"),
    "TRA": (ChainType.TCR_ALPHA,        "TCRA", "TRAV.fasta", None,          "TRAJ.fasta"),
    "TRB": (ChainType.TCR_BETA,         "TCRB", "TRBV.fasta", "TRBD.fasta", "TRBJ.fasta"),
    "TRG": (ChainType.TCR_GAMMA,        "TCRG", "TRGV.fasta", None,          "TRGJ.fasta"),
    "TRD": (ChainType.TCR_DELTA,        "TCRD", "TRDV.fasta", "TRDD.fasta", "TRDJ.fasta"),
}

# ─────────────────────────────────────────────────────────────
# Download / probe helpers
# ─────────────────────────────────────────────────────────────


def download_file(url: str, dest: str, retries: int = 3, delay: float = 1.0) -> bool:
    """Download ``url`` to ``dest`` with retries. Returns ``True`` on success,
    ``False`` on a 404 (locus absent for the species) or exhausted retries."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "GenAIRR-DataConfig-Builder/1.0"}
            )
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with open(dest, "wb") as fh:
                fh.write(data)
            return True
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return False
            logger.warning("HTTP %d for %s (attempt %d/%d)", exc.code, url, attempt + 1, retries)
        except Exception as exc:  # network hiccup — retry
            logger.warning("Error downloading %s: %s (attempt %d/%d)", url, exc, attempt + 1, retries)
        time.sleep(delay * (attempt + 1))
    return False


def probe_loci(imgt_species: str, cache_dir: str) -> dict:
    """Probe which loci exist for a species and cache their FASTA locally.

    Returns ``{locus_name: {"V": path, "J": path, "D": path?}}`` for every
    locus with (at least) a V and J file. A locus with no V or no J is
    skipped entirely; a D-bearing locus without a D file keeps V/J only.
    """
    results: dict = {}
    for locus_name, (_chain, _label, v_file, d_file, j_file) in LOCUS_DEFS.items():
        receptor = "IG" if locus_name.startswith("IG") else "TR"
        base_url = f"{IMGT_BASE}/{imgt_species}/{receptor}"
        species_dir = os.path.join(cache_dir, imgt_species, receptor)

        segments: dict = {}

        v_dest = os.path.join(species_dir, v_file)
        if os.path.exists(v_dest) or download_file(f"{base_url}/{v_file}", v_dest):
            segments["V"] = v_dest
        else:
            continue  # no V → skip the locus

        j_dest = os.path.join(species_dir, j_file)
        if os.path.exists(j_dest) or download_file(f"{base_url}/{j_file}", j_dest):
            segments["J"] = j_dest
        else:
            continue  # no J → skip

        if d_file is not None:
            d_dest = os.path.join(species_dir, d_file)
            if os.path.exists(d_dest) or download_file(f"{base_url}/{d_file}", d_dest):
                segments["D"] = d_dest

        results[locus_name] = segments
    return results


def fasta_has_sequences(path: str) -> bool:
    """True if ``path`` exists and contains at least one FASTA record."""
    if not os.path.exists(path):
        return False
    with open(path) as fh:
        return any(line.startswith(">") for line in fh)


# ─────────────────────────────────────────────────────────────
# Build (ported to ReferenceCartridgeBuilder)
# ─────────────────────────────────────────────────────────────


def build_cartridge(
    species_enum: Species,
    chain_type: ChainType,
    chain_label: str,
    species_label: str,
    v_path: str,
    j_path: str,
    d_path: str | None = None,
    reference_set: str = "IMGT",
) -> DataConfig | None:
    """Build one structural cartridge from IMGT FASTA.

    Runs the ``ReferenceCartridgeBuilder`` chain
    (``from_fasta → infer_identity → infer_v_subregions → build``) and
    returns the resulting :class:`GenAIRR.DataConfig`, or ``None`` if the
    build fails (e.g. a D-bearing locus whose D FASTA was unavailable).
    """
    if chain_type.has_d and d_path is None:
        logger.warning(
            "%s_%s: chain has a D segment but no D FASTA — skipping",
            species_label, chain_label,
        )
        return None
    try:
        return (
            ReferenceCartridgeBuilder.from_fasta(
                v_fasta=v_path,
                j_fasta=j_path,
                d_fasta=d_path,
                chain_type=chain_type,
            )
            .infer_identity(
                species=species_enum,
                locus=chain_label,
                reference_set=reference_set,
                name=f"{species_label}_{chain_label}_{reference_set}",
                source="IMGT",
            )
            .infer_v_subregions()
            .build()
        )
    except Exception as exc:  # malformed FASTA, anchor issues, etc.
        logger.error("Failed to build %s_%s: %s", species_label, chain_label, exc)
        return None


# ─────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build structural GenAIRR cartridges from IMGT V-QUEST FASTA.",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory to write built <SPECIES>_<LOCUS>_IMGT.pkl cartridges. "
             "REQUIRED (never defaults to the shipped builtin_dataconfigs set).",
    )
    parser.add_argument(
        "--cache-dir", default="/tmp/imgt_cache",
        help="Directory to cache downloaded IMGT FASTA (default: /tmp/imgt_cache).",
    )
    parser.add_argument(
        "--species", nargs="*", default=None,
        help="Only process these IMGT species directory names (default: all).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Probe and report available loci; do not build or write anything.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.cache_dir, exist_ok=True)

    species_list = args.species if args.species else IMGT_SPECIES
    built: list[str] = []
    skipped: list[tuple[str, str]] = []
    failed: list[str] = []

    for imgt_species in species_list:
        if imgt_species not in SPECIES_MAP:
            logger.warning("No Species enum mapping for %s — skipping", imgt_species)
            skipped.append((imgt_species, "no enum mapping"))
            continue

        species_enum, species_label = SPECIES_MAP[imgt_species]
        logger.info("=== %s (%s) ===", imgt_species, species_label)

        loci = probe_loci(imgt_species, args.cache_dir)
        if not loci:
            logger.info("  no loci found")
            skipped.append((imgt_species, "no loci"))
            continue

        for locus_name, segments in sorted(loci.items()):
            chain_type, chain_label, _, _, _ = LOCUS_DEFS[locus_name]
            config_name = f"{species_label}_{chain_label}_IMGT"

            if not fasta_has_sequences(segments["V"]):
                skipped.append((config_name, "empty V"))
                continue
            if not fasta_has_sequences(segments["J"]):
                skipped.append((config_name, "empty J"))
                continue

            d_path = segments.get("D")
            logger.info(
                "  %s: V=%s J=%s D=%s", locus_name,
                os.path.basename(segments["V"]),
                os.path.basename(segments["J"]),
                os.path.basename(d_path) if d_path else "none",
            )
            if args.dry_run:
                built.append(config_name)
                continue

            config = build_cartridge(
                species_enum, chain_type, chain_label, species_label,
                v_path=segments["V"], j_path=segments["J"], d_path=d_path,
            )
            if config is None:
                failed.append(config_name)
                continue

            n_v = sum(len(a) for a in config.v_alleles.values()) if config.v_alleles else 0
            n_d = sum(len(a) for a in config.d_alleles.values()) if config.d_alleles else 0
            n_j = sum(len(a) for a in config.j_alleles.values()) if config.j_alleles else 0

            pkl_path = os.path.join(args.output_dir, f"{config_name}.pkl")
            with open(pkl_path, "wb") as fh:
                pickle.dump(config, fh, protocol=pickle.HIGHEST_PROTOCOL)
            logger.info(
                "  built %s: %d V, %d D, %d J alleles → %s (%.1f KB)",
                config_name, n_v, n_d, n_j, pkl_path,
                os.path.getsize(pkl_path) / 1024,
            )
            built.append(config_name)

    print("\n" + "=" * 60)
    print(f"BUILT: {len(built)}")
    for name in sorted(built):
        print(f"  {name}")
    if skipped:
        print(f"\nSKIPPED: {len(skipped)}")
        for name, reason in skipped:
            print(f"  {name}: {reason}")
    if failed:
        print(f"\nFAILED: {len(failed)}")
        for name in failed:
            print(f"  {name}")
    print("=" * 60)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
