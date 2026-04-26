"""Fetch the canonical TP53 CDS (NM_000546.6) from NCBI E-utilities.

Outputs:
    data/raw/{accession}.fasta              - raw FASTA as returned by NCBI
    data/raw/{accession}.cds.txt            - the CDS as a single line of bases
    data/raw/{accession}.meta.json          - URL, timestamp, status, MD5
    data/processed/tp53_codon_index.csv     - residue -> CDS coords + codon
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Tuple

import requests

from .config import Config
from .utils import (
    codon_index_for_residue,
    ensure_dir,
    md5_str,
    project_root,
    translate,
    utc_timestamp,
    write_json,
)

log = logging.getLogger(__name__)


def fetch_fasta(cfg: Config, *, timeout: int = 30) -> Tuple[str, dict]:
    """Hit NCBI efetch and return (fasta_text, metadata)."""
    params = {
        "db": "nuccore",
        "id": cfg.refseq_accession,
        "rettype": "fasta_cds_na",
        "retmode": "text",
    }
    log.info("Fetching %s from NCBI E-utilities ...", cfg.refseq_accession)
    response = requests.get(cfg.ncbi_efetch_url, params=params, timeout=timeout)
    response.raise_for_status()
    fasta_text = response.text
    meta = {
        "accession": cfg.refseq_accession,
        "url": response.url,
        "http_status": response.status_code,
        "fetched_at_utc": utc_timestamp(),
        "content_length": len(fasta_text),
    }
    return fasta_text, meta


def parse_fasta_to_cds(fasta_text: str) -> str:
    """Return the concatenated bases of the first FASTA record."""
    bases = []
    seen_header = False
    for line in fasta_text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if seen_header:
                # Stop at the second record (we only want the first CDS).
                break
            seen_header = True
            continue
        bases.append(line.strip().upper())
    cds = "".join(bases)
    if not cds:
        raise ValueError("Empty CDS extracted from FASTA response.")
    return cds


def validate_cds(cds: str, cfg: Config) -> None:
    """Sanity-check the fetched CDS against expectations."""
    if len(cds) != cfg.expected_cds_length:
        raise ValueError(
            f"CDS length {len(cds)} != expected {cfg.expected_cds_length} "
            f"for {cfg.refseq_accession}."
        )
    if not cds.startswith("ATG"):
        raise ValueError(f"CDS does not start with ATG: {cds[:6]!r}")
    stop = cds[-3:]
    if stop not in {"TAA", "TAG", "TGA"}:
        raise ValueError(f"CDS does not end with a stop codon: {stop!r}")
    protein = translate(cds)
    if protein.count("*") != 1 or not protein.endswith("*"):
        raise ValueError("CDS contains an internal stop codon.")
    # Spot-check the canonical TP53 hotspot codon translations.
    expected_aa = {175: "R", 220: "Y", 245: "G", 248: "R",
                   249: "R", 273: "R", 282: "R"}
    for residue, aa in expected_aa.items():
        start, end = codon_index_for_residue(residue)
        codon = cds[start:end]
        observed = translate(codon)
        if observed != aa:
            raise ValueError(
                f"Residue {residue} codon {codon!r} translates to {observed!r}, "
                f"expected {aa!r}."
            )


def write_codon_index(cds: str, out_path: Path) -> None:
    """Write a residue -> CDS coords table for downstream modules."""
    ensure_dir(out_path.parent)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["residue_1based", "amino_acid", "cds_start_0based",
             "cds_end_0based_exclusive", "codon"]
        )
        # 393 amino acids + stop codon; iterate only over coding residues.
        for residue in range(1, len(cds) // 3):
            start, end = codon_index_for_residue(residue)
            codon = cds[start:end]
            aa = translate(codon)
            writer.writerow([residue, aa, start, end, codon])


def run(cfg: Config | None = None) -> Path:
    """Fetch + validate + persist the TP53 CDS. Return path to cds.txt."""
    cfg = cfg or Config.load()
    raw_dir = ensure_dir(project_root() / "data" / "raw")
    proc_dir = ensure_dir(project_root() / "data" / "processed")

    fasta_text, meta = fetch_fasta(cfg)
    cds = parse_fasta_to_cds(fasta_text)
    validate_cds(cds, cfg)

    fasta_path = raw_dir / f"{cfg.refseq_accession}.fasta"
    cds_path = raw_dir / f"{cfg.refseq_accession}.cds.txt"
    meta_path = raw_dir / f"{cfg.refseq_accession}.meta.json"
    codon_index_path = proc_dir / "tp53_codon_index.csv"

    fasta_path.write_text(fasta_text)
    cds_path.write_text(cds + "\n")
    meta["cds_md5"] = md5_str(cds)
    meta["cds_length"] = len(cds)
    write_json(meta_path, meta)
    write_codon_index(cds, codon_index_path)

    log.info(
        "Saved %s (length=%d, md5=%s)",
        cfg.refseq_accession, len(cds), meta["cds_md5"]
    )
    return cds_path


if __name__ == "__main__":  # pragma: no cover
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run()
