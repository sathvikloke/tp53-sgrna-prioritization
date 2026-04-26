"""Tests for src/score_on_target.py (feature backend)."""

from __future__ import annotations

import pytest

from src.score_on_target import score_30mer_features


def test_score_in_unit_interval():
    seq = "AAAA" + "ACGTACGTACGTACGTACGT" + "AGG" + "AAA"   # 4+20+3+3 = 30
    s = score_30mer_features(seq)
    assert 0.0 <= s <= 1.0


def test_score_rejects_wrong_length():
    with pytest.raises(ValueError):
        score_30mer_features("ACGT")


def test_high_gc_balanced_better_than_extreme():
    """A balanced spacer should outscore an all-A extreme."""
    balanced = "ACGT" + "ACGTACGTACGTACGTACGT" + "AGG" + "AAA"
    poor    = "AAAA" + "AAAAAAAAAAAAAAAAAAAA" + "AGG" + "AAA"
    assert score_30mer_features(balanced) > score_30mer_features(poor)


def test_pam_n_preference():
    spacer = "ACGT" + "ACGTACGTACGTACGTACGT"
    after_pam = "AAA"
    s_g = score_30mer_features(spacer + "GGG" + after_pam)
    s_a = score_30mer_features(spacer + "AGG" + after_pam)
    s_c = score_30mer_features(spacer + "CGG" + after_pam)
    # G > A > C bonus
    assert s_g >= s_a >= s_c
