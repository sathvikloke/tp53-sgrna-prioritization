"""Tests for src/cfd.py."""

from __future__ import annotations

import pytest

from src.cfd import aggregate_score, cfd


SPACER = "ACGTACGTACGTACGTACGT"  # 20 nt
assert len(SPACER) == 20


def test_perfect_match_ngg_returns_one():
    # NGG (TGG) PAM and perfect match -> score = 1.0
    assert cfd(SPACER, SPACER, "TGG") == pytest.approx(1.0)


def test_unknown_pam_returns_zero():
    assert cfd(SPACER, SPACER, "AAA") == 0.0


def test_seed_mismatch_penalized_more_than_distal():
    distal_off = "T" + SPACER[1:]      # mismatch at position 1 (distal)
    seed_off = SPACER[:-1] + "T" if SPACER[-1] != "T" else SPACER[:-1] + "A"
    score_distal = cfd(distal_off, SPACER, "AGG")
    score_seed = cfd(seed_off, SPACER, "AGG")
    assert score_distal > score_seed


def test_aggregate_in_unit_interval():
    assert aggregate_score([]) == pytest.approx(1.0)
    s = aggregate_score([0.1, 0.2, 0.5])
    assert 0.0 < s <= 1.0


def test_aggregate_decreases_with_more_off_targets():
    s1 = aggregate_score([0.5])
    s2 = aggregate_score([0.5, 0.5])
    assert s1 > s2


def test_cfd_requires_20mer():
    with pytest.raises(ValueError):
        cfd("ACGT", SPACER, "AGG")
