"""Tests for src/utils.py."""

from __future__ import annotations

import pytest

from src.utils import (
    codon_index_for_residue,
    gc_content,
    matches_pam,
    reverse_complement,
    translate,
)


def test_reverse_complement_basic():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAA") == "TTT"
    assert reverse_complement("ATCGN") == "NCGAT"


def test_reverse_complement_double_inversion():
    seq = "GATTACA"
    assert reverse_complement(reverse_complement(seq)) == seq


def test_translate_known_codons():
    assert translate("ATG") == "M"
    assert translate("TAA") == "*"
    assert translate("ATGGCC") == "MA"


def test_translate_invalid_length():
    with pytest.raises(ValueError):
        translate("ATGA")


def test_gc_content_extremes():
    assert gc_content("") == 0.0
    assert gc_content("AAAA") == 0.0
    assert gc_content("GGGG") == 1.0
    assert gc_content("ATGC") == 0.5


def test_codon_index_for_residue():
    assert codon_index_for_residue(1) == (0, 3)
    assert codon_index_for_residue(2) == (3, 6)
    assert codon_index_for_residue(175) == (522, 525)
    with pytest.raises(ValueError):
        codon_index_for_residue(0)


def test_matches_pam_iupac():
    assert matches_pam("AGG", "NGG") is True
    assert matches_pam("CGG", "NGG") is True
    assert matches_pam("AAG", "NGG") is False
    assert matches_pam("AAG", "NRG") is True   # R = A/G
    assert matches_pam("AGT", "NGG") is False
    assert matches_pam("AG",  "NGG") is False  # length mismatch
