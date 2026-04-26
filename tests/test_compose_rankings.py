"""Tests for src/compose_rankings.py."""

from __future__ import annotations

from src.compose_rankings import (
    _build_residue_freq,
    _compose_for_cohort,
    _normalize_minmax,
    _tier,
)


def test_normalize_minmax_basic():
    out = _normalize_minmax([1.0, 2.0, 3.0])
    assert out == [0.0, 0.5, 1.0]


def test_normalize_constant_returns_midpoint():
    assert _normalize_minmax([0.5, 0.5, 0.5]) == [0.5, 0.5, 0.5]


def test_normalize_empty():
    assert _normalize_minmax([]) == []


def test_build_residue_freq_sums_per_residue():
    rows = [
        {"cohort_id": "A", "residue": "248", "frequency": "0.10"},
        {"cohort_id": "A", "residue": "248", "frequency": "0.15"},
        {"cohort_id": "A", "residue": "175", "frequency": "0.20"},
        {"cohort_id": "B", "residue": "175", "frequency": "0.05"},
    ]
    out = _build_residue_freq(rows)
    assert out["A"][248] == 0.25
    assert out["A"][175] == 0.20
    assert out["B"][175] == 0.05


def test_tier_assignment():
    assert _tier(0) == "tier_1"
    assert _tier(1) == "tier_2"
    assert _tier(2) == "tier_2"
    assert _tier(3) == "tier_3"


def test_compose_for_cohort_orders_correctly():
    candidates = [
        {"sgRNA_id": "g1", "residue": "175", "on_target_score": "0.9",
         "off_target_score": "0.9"},
        {"sgRNA_id": "g2", "residue": "175", "on_target_score": "0.5",
         "off_target_score": "0.5"},
        {"sgRNA_id": "g3", "residue": "248", "on_target_score": "0.7",
         "off_target_score": "0.7"},
    ]
    weights = {"on_target": 0.5, "off_target": 0.3, "frequency": 0.2}
    residue_freq = {175: 0.30, 248: 0.10}
    ranked = _compose_for_cohort(candidates, residue_freq, weights)
    # g1 should rank above g2 within residue 175
    by_id = {r["sgRNA_id"]: r for r in ranked}
    assert float(by_id["g1"]["composite_score"]) > float(by_id["g2"]["composite_score"])
    # g1 (residue 175, freq 0.30) should rank above g3 (residue 248, freq 0.10) overall
    g1_rank = next(i for i, r in enumerate(ranked) if r["sgRNA_id"] == "g1")
    g3_rank = next(i for i, r in enumerate(ranked) if r["sgRNA_id"] == "g3")
    # We sort by (residue, -composite); g1 is residue 175 (sorts before 248).
    assert g1_rank < g3_rank
    # g2 should be tier_2 (second within residue 175).
    assert by_id["g2"]["tier"] == "tier_2"
    assert by_id["g1"]["tier"] == "tier_1"
