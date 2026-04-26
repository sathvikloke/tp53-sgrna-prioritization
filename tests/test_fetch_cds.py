"""Tests for src/fetch_cds.py (offline, no NCBI calls)."""

from __future__ import annotations

import pytest

from src.fetch_cds import parse_fasta_to_cds, validate_cds


def _make_min_valid_cds():
    """Return a length-1182 CDS with hotspot codons matching the real TP53.

    For testing purposes we only need:
        - starts with ATG
        - ends with a stop codon
        - exactly 1182 nt
        - residues 175,220,245,248,249,273,282 translate to R,Y,G,R,R,R,R
        - no internal stop codons
    """
    # 394 codons total (including stop). We fill with "GCC" (Ala) which has
    # no stop codons, then patch the required codons.
    codons = ["GCC"] * 394
    codons[0] = "ATG"          # M start
    codons[174] = "CGT"        # R175
    codons[219] = "TAT"        # Y220
    codons[244] = "GGT"        # G245
    codons[247] = "CGT"        # R248
    codons[248] = "CGT"        # R249
    codons[272] = "CGT"        # R273
    codons[281] = "CGT"        # R282
    codons[393] = "TAA"        # stop
    return "".join(codons)


def test_parse_fasta_extracts_cds():
    fasta = ">NM_000546.6 cds\nATGGCC\nGCCTAA\n"
    assert parse_fasta_to_cds(fasta) == "ATGGCCGCCTAA"


def test_parse_fasta_only_first_record():
    fasta = ">first\nATGAAA\n>second\nGGGCCC\n"
    assert parse_fasta_to_cds(fasta) == "ATGAAA"


def test_parse_fasta_empty_raises():
    with pytest.raises(ValueError):
        parse_fasta_to_cds(">only header\n")


def test_validate_cds_accepts_canonical(monkeypatch):
    from src.config import Config, TP53Locus

    cfg = Config(
        raw={"run": {"version": "test"}},
        cas_enzyme="SpCas9", pam="NGG", spacer_length=20,
        cut_offset_from_pam=3, pam_rescue={}, edit_window_nt=10,
        framing_mode="hdr",
        hotspots=[], hotspot_variants=[], cohorts=[],
        refseq_accession="NM_000546.6",
        expected_cds_length=1182,
        ncbi_efetch_url="x",
        tp53_locus=TP53Locus(
            chromosome="chr17", start=7668421, end=7687490, strand="-",
            exons_grch38=[],
        ),
        cbioportal_base_url="x",
        tp53_entrez_id=7157,
        weights={"on_target": 0.5, "off_target": 0.3, "frequency": 0.2},
        off_target={},
    )
    cds = _make_min_valid_cds()
    assert len(cds) == 1182
    validate_cds(cds, cfg)   # must not raise
