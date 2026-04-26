"""Per-cohort composite scoring + ranking.

Inputs:
    data/processed/candidates.scored.csv         (sgRNA + on/off scores)
    data/processed/cohort_frequencies.csv        (per-cohort variant counts)

Outputs:
    results/tables/{cohort_id}_ranked_sgRNAs.csv
    results/tables/all_cohorts_top3.csv
    results/tables/pam_deserts.csv               (residues with no candidates)
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Dict, List

from .config import Config
from .utils import ensure_dir, project_root

log = logging.getLogger(__name__)


def _read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh))


def _normalize_minmax(values: List[float]) -> List[float]:
    if not values:
        return []
    lo, hi = min(values), max(values)
    if hi <= lo:
        return [0.5 for _ in values]
    return [(v - lo) / (hi - lo) for v in values]


def _build_residue_freq(freq_rows: List[Dict[str, str]]
                        ) -> Dict[str, Dict[int, float]]:
    """Return {cohort_id: {residue: total_variant_frequency_at_that_residue}}.

    For residues with multiple hotspot variants (e.g. R248Q + R248W), we
    use the SUM of variant frequencies as the cohort-relevance signal at
    that residue, because either variant is a valid editing target and a
    single sgRNA may correct both.
    """
    out: Dict[str, Dict[int, float]] = {}
    for row in freq_rows:
        cohort = row["cohort_id"]
        residue = int(row["residue"])
        freq = float(row["frequency"])
        out.setdefault(cohort, {})
        out[cohort][residue] = out[cohort].get(residue, 0.0) + freq
    return out


def _tier(idx_in_residue: int) -> str:
    if idx_in_residue == 0:
        return "tier_1"
    if idx_in_residue <= 2:
        return "tier_2"
    return "tier_3"


def _compose_for_cohort(candidates: List[Dict[str, str]],
                        residue_freq: Dict[int, float],
                        weights: Dict[str, float]) -> List[Dict[str, str]]:
    """Compute composite scores for one cohort. Returns sorted, tiered rows.

    Two composite scores are emitted (Flag #5 fix):

      * `composite_score` — within-residue min-max normalization of on/off,
        with frequency normalized across residues. This scale is **only**
        meaningful for ranking spacers WITHIN the same residue, because the
        within-residue normalization erases absolute differences in spacer
        quality across residues. Use this for tier_1/2/3 assignment.

      * `composite_score_absolute` — uses the raw on/off values directly
        (already in [0,1] by construction: on_target_score is the Doench
        feature score / Azimuth probability; off_target_score is the
        CRISPOR aggregate specificity 100/(100+sum CFDs) rescaled to [0,1]).
        Frequency is min-max normalized across residues. This scale IS
        comparable across residues and is the recommended sort key for the
        cross-residue master table.
    """
    if not candidates:
        return []

    # --- step 1: compute frequency-normalized values across all residues ---
    # We min-max normalize *across* the residues that appear in this run so
    # that frequency scales the same way as on/off target.
    freqs = [residue_freq.get(int(c["residue"]), 0.0) for c in candidates]
    freq_norm = _normalize_minmax(freqs)

    # --- step 2: normalize on/off WITHIN residue (for tier ranking) ---
    by_res: Dict[int, List[int]] = {}
    for i, c in enumerate(candidates):
        by_res.setdefault(int(c["residue"]), []).append(i)

    on_norm = [0.0] * len(candidates)
    off_norm = [0.0] * len(candidates)
    for indices in by_res.values():
        on_vals = [float(candidates[i]["on_target_score"]) for i in indices]
        off_vals = [float(candidates[i]["off_target_score"]) for i in indices]
        on_n = _normalize_minmax(on_vals)
        off_n = _normalize_minmax(off_vals)
        for idx, on_v, off_v in zip(indices, on_n, off_n):
            on_norm[idx] = on_v
            off_norm[idx] = off_v

    # --- step 3: composite scores (within-residue + absolute) ---
    out: List[Dict[str, str]] = []
    for i, c in enumerate(candidates):
        on_raw = float(c["on_target_score"])
        off_raw = float(c["off_target_score"])
        composite = (
            weights["on_target"]  * on_norm[i] +
            weights["off_target"] * off_norm[i] +
            weights["frequency"]  * freq_norm[i]
        )
        composite_abs = (
            weights["on_target"]  * on_raw +
            weights["off_target"] * off_raw +
            weights["frequency"]  * freq_norm[i]
        )
        row = dict(c)
        row["on_target_norm"] = f"{on_norm[i]:.6f}"
        row["off_target_norm"] = f"{off_norm[i]:.6f}"
        row["frequency_norm"] = f"{freq_norm[i]:.6f}"
        row["cohort_residue_frequency"] = f"{residue_freq.get(int(c['residue']), 0.0):.6f}"
        row["composite_score"] = f"{composite:.6f}"
        row["composite_score_absolute"] = f"{composite_abs:.6f}"
        out.append(row)

    # --- step 4: sort and assign tiers within residue (uses within-residue
    # composite, since tier semantics are "best within this residue"). ---
    out.sort(key=lambda r: (int(r["residue"]),
                            -float(r["composite_score"])))
    res_seen: Dict[int, int] = {}
    for r in out:
        res = int(r["residue"])
        idx = res_seen.get(res, 0)
        r["tier"] = _tier(idx)
        res_seen[res] = idx + 1
    return out


def run(cfg: Config | None = None) -> Path:
    cfg = cfg or Config.load()
    proc = project_root() / "data" / "processed"
    cands_path = proc / "candidates.scored.csv"
    freqs_path = proc / "cohort_frequencies.csv"
    if not cands_path.exists():
        raise FileNotFoundError(f"{cands_path} missing; run off-target step first.")
    if not freqs_path.exists():
        raise FileNotFoundError(f"{freqs_path} missing; run cohort-fetch first.")

    candidates = _read_csv(cands_path)
    freqs = _read_csv(freqs_path)
    residue_freq = _build_residue_freq(freqs)

    tables_dir = ensure_dir(project_root() / "results" / "tables")
    all_top3: List[Dict[str, str]] = []

    for cohort in cfg.cohorts:
        cohort_residue_freq = residue_freq.get(cohort.id, {})
        ranked = _compose_for_cohort(candidates, cohort_residue_freq, cfg.weights)
        for r in ranked:
            r["cohort_id"] = cohort.id
            r["cohort_label"] = cohort.label
        out_path = tables_dir / f"{cohort.id}_ranked_sgRNAs.csv"
        if ranked:
            fieldnames = list(ranked[0].keys())
            with out_path.open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fieldnames)
                w.writeheader()
                w.writerows(ranked)
        else:
            out_path.write_text("# no candidates passed enumeration\n")
        log.info("Wrote %s (n=%d)", out_path, len(ranked))

        # Top-3 per residue per cohort for the master table.
        seen_per_res: Dict[int, int] = {}
        for r in ranked:
            res = int(r["residue"])
            seen_per_res[res] = seen_per_res.get(res, 0) + 1
            if seen_per_res[res] <= 3:
                all_top3.append(r)

    master_path = tables_dir / "all_cohorts_top3.csv"
    if all_top3:
        # Sort the master table by the ABSOLUTE composite (Flag #5): the
        # within-residue normalized score is not meaningful across residues.
        all_top3.sort(key=lambda r: -float(r.get("composite_score_absolute",
                                                 r["composite_score"])))
        fieldnames = list(all_top3[0].keys())
        with master_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(all_top3)
    else:
        master_path.write_text("# no top-3 rows\n")
    log.info("Wrote %s (n=%d)", master_path, len(all_top3))

    # PAM-desert table: residues with zero candidates.
    desert_path = tables_dir / "pam_deserts.csv"
    seen_residues = {int(c["residue"]) for c in candidates}
    desert_rows: List[Dict[str, str]] = []
    for h in cfg.hotspots:
        if h.residue not in seen_residues:
            desert_rows.append({
                "residue": str(h.residue),
                "ref_aa": h.ref_aa,
                "label": h.label,
                "spcas9_NGG_status": "PAM_DESERT",
                "rescue_options": ", ".join(
                    f"{cas}({pam})" for cas, pam in cfg.pam_rescue.items()
                ),
            })
    with desert_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["residue", "ref_aa", "label",
                    "spcas9_NGG_status", "rescue_options"])
        for r in desert_rows:
            w.writerow([r["residue"], r["ref_aa"], r["label"],
                        r["spcas9_NGG_status"], r["rescue_options"]])
    log.info("Wrote %s (n=%d residues with no SpCas9 candidates)",
             desert_path, len(desert_rows))
    return master_path


if __name__ == "__main__":  # pragma: no cover
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run()
