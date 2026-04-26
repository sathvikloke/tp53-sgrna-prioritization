"""Tests for src/enumerate_sgrnas.py."""

from __future__ import annotations

from src.config import Hotspot
from src.enumerate_sgrnas import (
    candidates_for_hotspot,
    _iter_antisense_pams,
    _iter_sense_pams,
)


def test_iter_sense_pams_finds_ngg():
    seq = "AAAAGGAAAA"   # NGG starts at index 3 (A,G,G) and 4? let's check
    # cds[3]='A', cds[4]='G', cds[5]='G' -> NGG at p=3
    pams = list(_iter_sense_pams(seq))
    assert 3 in pams
    # cds[4]='G', cds[5]='G', cds[6]='A' is GGA, not NGG (needs G at p+1, G at p+2).
    # Our test: only p=3 should qualify.
    assert pams == [3]


def test_iter_antisense_pams_finds_ccn():
    seq = "AAACCAAAA"      # CCN at q=3 (C,C,A)
    pams = list(_iter_antisense_pams(seq))
    assert 3 in pams


def test_candidates_cut_site_in_window():
    """Build a synthetic CDS where exactly ONE sense sgRNA should pass.

    Codon 11 (residue 11) sits at CDS[30:33]. If we place an NGG at
    CDS[33:36], the cut site lands at CDS[30] (offset 0). The protospacer
    is CDS[13:33] (20 nt) and the 30-mer context is CDS[9:39].
    """
    cds = (
        "ATGAAAAAAAAA"           # 12 nt: 4 codons (residues 1-4)
        "AAAAAAAAAAAAAAAAAA"     # 18 nt fill (residues 5-10)
        "CGT"                    # codon 11 = R (CGT) at CDS[30:33]
        "AGG"                    # NGG at CDS[33:36]; cut at CDS[30]
        "AAAAAAAAA"              # 9 nt right context
    )
    assert len(cds) == 12 + 18 + 3 + 3 + 9
    h = Hotspot(residue=11, ref_aa="R", label="R11")
    cands = candidates_for_hotspot(cds, h, edit_window_nt=10, spacer_length=20)
    sense = [c for c in cands if c["strand"] == "+"]
    assert len(sense) >= 1
    sense_top = sense[0]
    assert sense_top["pam"] == "AGG"
    assert sense_top["cut_position_cds_0based"] == 30
    assert sense_top["codon_start_cds_0based"] == 30
    assert sense_top["cut_offset_from_codon_start_nt"] == 0


def test_candidates_outside_window_excluded():
    """An NGG far from the codon must NOT yield a candidate."""
    cds = (
        "ATG" + "A" * 30                      # residue 1 codon + 30 A's
        + "CGT"                               # residue ~12 codon (not target)
        + "A" * 100                           # 100 nt buffer
        + "AGG"                               # an NGG far from residue 1
        + "A" * 30
    )
    h = Hotspot(residue=1, ref_aa="M", label="M1")
    cands = candidates_for_hotspot(cds, h, edit_window_nt=10, spacer_length=20)
    # No NGG within +/- 10 of CDS[0] in this sequence (the only NGG is far).
    assert all(abs(c["cut_offset_from_codon_start_nt"]) <= 10 for c in cands)


def test_antisense_candidate_well_formed():
    """Verify that an antisense candidate produces a 20-nt revcomp spacer."""
    # Construct a CDS where a CCN is positioned so the cut lands on residue 5's codon.
    # Residue 5 codon: CDS[12:15]. Antisense cut at CDS[q+6]; we want q+6 = 12 -> q = 6.
    # CDS[6:9] must be "CCN". Protospacer on sense at CDS[9:29].
    cds = (
        "ATGAAA"                              # 0..5: residues 1..2
        + "CCA"                               # 6..8: CCN (q=6)
        + "T" * 20                            # 9..28: 20-nt protospacer (sense)
        + "ACGT" * 5                          # padding for context
    )
    # residue 5 codon at CDS[12:15] (within the T's)
    h = Hotspot(residue=5, ref_aa="L", label="L5")
    cands = candidates_for_hotspot(cds, h, edit_window_nt=10, spacer_length=20)
    antisense = [c for c in cands if c["strand"] == "-"]
    assert any(c["cut_position_cds_0based"] == 12 for c in antisense)
    for c in antisense:
        assert len(c["spacer_5p_to_3p"]) == 20
        assert len(c["context_30mer"]) == 30
