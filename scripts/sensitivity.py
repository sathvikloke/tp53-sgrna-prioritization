#!/usr/bin/env python3
"""Weight sensitivity analysis (Flag #6).

The composite score uses three weights — on_target, off_target,
frequency — that the user fixed in config/params.yaml. Their default
values (0.5, 0.3, 0.2) reflect a defensible prioritization but are
ultimately a choice. This script re-runs the composite step over a
small grid of weight settings and reports how much the top-tier
selections move.

The intent is NOT to retune weights for any specific outcome; it is
to provide a transparent picture of which sgRNAs are robust to
weight perturbations and which only rank because of one particular
choice.

Inputs:
    data/processed/candidates.scored.csv       (on + off scores)
    data/processed/cohort_frequencies.csv      (per-cohort frequencies)

Outputs:
    results/tables/sensitivity_summary.csv
        one row per (cohort, sgRNA) summarizing rank stability
    results/tables/sensitivity_grid.csv
        full per-(cohort, weight-setting, sgRNA) ranks

Run:
    python -m scripts.sensitivity
"""
from __future__ import annotations

import csv
import logging
from itertools import product
from pathlib import Path
from typing import Any, Dict, List, Tuple

from src.compose_rankings import _build_residue_freq, _compose_for_cohort, _read_csv
from src.config import Config
from src.utils import ensure_dir, project_root


log = logging.getLogger("sensitivity")


def _grid() -> List[Dict[str, float]]:
    """Generate a small simplex of (on, off, freq) weight settings.

    All non-negative, summing to 1, with each weight in {0.2, 0.3, ..., 0.8}.
    """
    out: List[Dict[str, float]] = []
    steps = [0.1 * i for i in range(2, 9)]      # 0.2 ... 0.8
    for on, off in product(steps, steps):
        freq = round(1.0 - on - off, 6)
        if freq < 0.0 or freq > 0.8 + 1e-9:
            continue
        out.append({
            "on_target": round(on, 6),
            "off_target": round(off, 6),
            "frequency": round(freq, 6),
        })
    # Always include the configured default at the front.
    return out


def _rank_within_cohort(rows: List[Dict[str, str]]) -> Dict[str, int]:
    """Return {sgRNA_id: 1-based rank by composite_score (desc), residue-grouped}."""
    ranks: Dict[str, int] = {}
    by_res: Dict[int, List[Dict[str, str]]] = {}
    for r in rows:
        by_res.setdefault(int(r["residue"]), []).append(r)
    for residue, group in by_res.items():
        group.sort(key=lambda r: -float(r["composite_score"]))
        for i, r in enumerate(group, start=1):
            ranks[r["sgRNA_id"]] = i
    return ranks


def run() -> Tuple[Path, Path]:
    cfg = Config.load()
    proc = project_root() / "data" / "processed"
    cands = _read_csv(proc / "candidates.scored.csv")
    freqs = _read_csv(proc / "cohort_frequencies.csv")
    residue_freq = _build_residue_freq(freqs)

    grid = _grid()
    log.info("Sensitivity grid: %d weight settings.", len(grid))

    # full grid: (cohort, weights, sgRNA) -> rank
    grid_rows: List[Dict[str, Any]] = []
    summary: Dict[Tuple[str, str], Dict[str, Any]] = {}

    for cohort in cfg.cohorts:
        cohort_freq = residue_freq.get(cohort.id, {})
        for w in grid:
            ranked = _compose_for_cohort(cands, cohort_freq, w)
            ranks = _rank_within_cohort(ranked)
            for r in ranked:
                grid_rows.append({
                    "cohort_id": cohort.id,
                    "w_on": w["on_target"],
                    "w_off": w["off_target"],
                    "w_freq": w["frequency"],
                    "sgRNA_id": r["sgRNA_id"],
                    "residue": r["residue"],
                    "composite_score": r["composite_score"],
                    "composite_score_absolute": r.get("composite_score_absolute", ""),
                    "rank_within_residue": ranks[r["sgRNA_id"]],
                })
                key = (cohort.id, r["sgRNA_id"])
                rec = summary.setdefault(key, {
                    "cohort_id": cohort.id,
                    "sgRNA_id": r["sgRNA_id"],
                    "residue": r["residue"],
                    "n_settings": 0,
                    "ranks": [],
                })
                rec["n_settings"] += 1
                rec["ranks"].append(ranks[r["sgRNA_id"]])

    out_dir = ensure_dir(project_root() / "results" / "tables")
    grid_path = out_dir / "sensitivity_grid.csv"
    if grid_rows:
        with grid_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(grid_rows[0].keys()))
            w.writeheader()
            w.writerows(grid_rows)
    log.info("Wrote %s (n=%d)", grid_path, len(grid_rows))

    summary_rows: List[Dict[str, Any]] = []
    for (cohort_id, sgRNA_id), rec in summary.items():
        ranks = rec["ranks"]
        summary_rows.append({
            "cohort_id": cohort_id,
            "sgRNA_id": sgRNA_id,
            "residue": rec["residue"],
            "n_settings": rec["n_settings"],
            "rank_min": min(ranks),
            "rank_max": max(ranks),
            "rank_median": sorted(ranks)[len(ranks) // 2],
            "frac_top_in_residue": (
                sum(1 for r in ranks if r == 1) / len(ranks)
            ),
        })
    summary_rows.sort(key=lambda r: (r["cohort_id"], int(r["residue"]),
                                     r["rank_median"]))

    summary_path = out_dir / "sensitivity_summary.csv"
    if summary_rows:
        with summary_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
            w.writeheader()
            w.writerows(summary_rows)
    log.info("Wrote %s (n=%d)", summary_path, len(summary_rows))
    return grid_path, summary_path


if __name__ == "__main__":  # pragma: no cover
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run()
