"""Create summary outputs: results/summary.json + figures."""

from __future__ import annotations

import csv
import json
import logging
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Dict, List

from .config import Config
from .utils import (
    ensure_dir,
    md5_str,
    project_root,
    read_json,
    utc_timestamp,
    write_json,
)

log = logging.getLogger(__name__)


def _safe_version(pkg: str) -> str:
    try:
        return version(pkg)
    except PackageNotFoundError:
        return "n/a"


def _read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh))


def _build_summary(cfg: Config) -> Dict[str, Any]:
    root = project_root()
    cds_path = root / "data" / "raw" / f"{cfg.refseq_accession}.cds.txt"
    cds_meta_path = root / "data" / "raw" / f"{cfg.refseq_accession}.meta.json"
    candidates_path = root / "data" / "processed" / "candidates.scored.csv"
    freqs_path = root / "data" / "processed" / "cohort_frequencies.csv"

    summary: Dict[str, Any] = {
        "run": {
            "version": cfg.raw["run"]["version"],
            "generated_at_utc": utc_timestamp(),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "package_versions": {
                pkg: _safe_version(pkg)
                for pkg in ("requests", "pyyaml", "pandas",
                            "numpy", "matplotlib")
            },
        },
        "config": {
            "weights": cfg.weights,
            "edit_window_nt": cfg.edit_window_nt,
            "cas_enzyme": cfg.cas_enzyme,
            "pam": cfg.pam,
            "hotspots": [h.label for h in cfg.hotspots],
            "cohorts": [c.id for c in cfg.cohorts],
        },
        "inputs": {},
        "outputs": {},
    }

    if cds_meta_path.exists():
        summary["inputs"]["tp53_cds"] = read_json(cds_meta_path)
    if cds_path.exists():
        cds = cds_path.read_text().strip()
        summary["inputs"]["tp53_cds_md5"] = md5_str(cds)
        summary["inputs"]["tp53_cds_length"] = len(cds)

    summary["outputs"]["n_candidates"] = len(_read_csv(candidates_path))
    summary["outputs"]["n_cohort_freq_rows"] = len(_read_csv(freqs_path))
    return summary


def _make_figures(cfg: Config) -> List[Path]:
    """Generate figures if matplotlib is available; else skip silently."""
    try:
        import matplotlib  # noqa: F401
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        log.warning("matplotlib unavailable; skipping figures (%s)", exc)
        return []

    figs_dir = ensure_dir(project_root() / "results" / "figures")
    out_paths: List[Path] = []

    freqs_path = project_root() / "data" / "processed" / "cohort_frequencies.csv"
    rows = _read_csv(freqs_path)
    if rows:
        # Hotspot frequency heatmap (cohort x variant).
        cohorts = sorted({r["cohort_id"] for r in rows})
        variants = cfg.hotspot_variants
        matrix = [
            [
                next(
                    (float(r["frequency"]) for r in rows
                     if r["cohort_id"] == c and r["variant"] == v),
                    0.0
                )
                for v in variants
            ]
            for c in cohorts
        ]
        fig, ax = plt.subplots(figsize=(8, 4))
        im = ax.imshow(matrix, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(variants)))
        ax.set_xticklabels(variants, rotation=45, ha="right")
        ax.set_yticks(range(len(cohorts)))
        ax.set_yticklabels(cohorts)
        ax.set_title("TP53 hotspot frequency by cohort")
        fig.colorbar(im, ax=ax, label="frequency")
        fig.tight_layout()
        path = figs_dir / "hotspot_freq_heatmap.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        out_paths.append(path)

    # Composite score distribution per cohort.
    tables_dir = project_root() / "results" / "tables"
    composites: Dict[str, List[float]] = {}
    for cohort in cfg.cohorts:
        f = tables_dir / f"{cohort.id}_ranked_sgRNAs.csv"
        rows2 = _read_csv(f)
        composites[cohort.label] = [
            float(r["composite_score"]) for r in rows2 if "composite_score" in r
        ]
    if any(v for v in composites.values()):
        fig, ax = plt.subplots(figsize=(8, 4))
        labels = list(composites.keys())
        ax.boxplot(
            [composites[l] for l in labels], labels=labels, vert=True,
        )
        ax.set_ylabel("Composite score")
        ax.set_title("Composite sgRNA score by cohort")
        fig.tight_layout()
        path = figs_dir / "composite_score_boxplot.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        out_paths.append(path)

    return out_paths


def run(cfg: Config | None = None) -> Path:
    cfg = cfg or Config.load()
    summary = _build_summary(cfg)
    fig_paths = _make_figures(cfg)
    summary["outputs"]["figures"] = [str(p.relative_to(project_root()))
                                     for p in fig_paths]
    out = project_root() / "results" / "summary.json"
    write_json(out, summary)
    log.info("Wrote %s", out)
    return out


if __name__ == "__main__":  # pragma: no cover
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run()
