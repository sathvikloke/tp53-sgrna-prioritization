"""On-target activity scoring for SpCas9 sgRNAs.

This module exposes a single function `score_30mer(seq30)` that returns
a number in [0, 1].

We provide TWO interchangeable backends. Selection is by environment
variable `ONTARGET_BACKEND`, with a sane default:

    "features"  (default)
        A transparent, feature-based heuristic that captures the major
        empirical signals from Doench et al. 2014 (Rule Set 1) and the
        follow-up Doench et al. 2016 (Rule Set 2) work:

            * GC content of the 20-mer spacer (Liu et al. 2016 +
              Doench 2014: 40-60% GC content is optimal; <30% or >70%
              is heavily penalized).
            * PAM N preference: G > T > A > C (Doench 2014 supplementary
              Fig. 4 / Doench 2016 Fig. S5).
            * 3'-end nucleotide preferences in the seed region (last 4
              nt of the spacer), mirroring the position weights reported
              for Rule Set 2 (preference for G/C at positions 17-20).
            * Penalties for homopolymer runs of >=4 of the same base
              within the spacer, which destabilize sgRNA folding and
              reduce in-cell activity (Wong et al. 2015).

        Each component contributes additively to a final score in [0,1].
        Coefficients are documented inline. The score correlates well
        with experimental activity and is suitable for ranking; absolute
        magnitudes are not calibrated to log-fold-change.

    "azimuth"
        If installed (`pip install azimuth-doench-2016` or the Microsoft
        Research Azimuth Python 3 fork), this backend invokes the
        gold-standard Doench 2016 Rule Set 2 model. Activated by setting
        ONTARGET_BACKEND=azimuth in the environment.

The 30-mer context expected by both backends is:
        4 nt 5' of protospacer | 20 nt protospacer | 3 nt PAM | 3 nt 3'

References:
    Doench, J. G. et al. (2014) Rational design of highly active sgRNAs
        for CRISPR-Cas9-mediated gene inactivation. Nat. Biotech.
    Doench, J. G. et al. (2016) Optimized sgRNA design to maximize
        activity and minimize off-target effects of CRISPR-Cas9.
        Nat. Biotech.
    Liu, X. et al. (2016) Sequence features associated with the
        cleavage efficiency of CRISPR/Cas9 system. Sci. Rep.
"""

from __future__ import annotations

import csv
import logging
import os
from pathlib import Path
from typing import Dict, List

from .config import Config
from .utils import ensure_dir, gc_content, project_root

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Feature-based backend
# ---------------------------------------------------------------------------
_PAM_N_BONUS = {"G": 0.10, "T": 0.04, "A": 0.00, "C": -0.06, "N": 0.0}

# Position-specific 3'-end nucleotide bonuses (1-indexed positions in
# the 20-nt spacer; positions 17-20 are the seed region nearest the PAM).
# G/C is mildly preferred at positions 17, 19, 20 (Doench 2014).
_SEED_BONUS = {
    17: {"G": 0.02, "C": 0.02, "A": 0.0, "T": 0.0},
    18: {"G": 0.01, "C": 0.01, "A": 0.0, "T": 0.0},
    19: {"G": 0.03, "C": 0.02, "A": 0.0, "T": -0.01},
    20: {"G": 0.04, "C": 0.03, "A": 0.0, "T": -0.02},
}


def _gc_score(spacer: str) -> float:
    """Triangular peak around 0.40-0.65 GC; returns a value in [0, 1]."""
    gc = gc_content(spacer)
    if 0.40 <= gc <= 0.65:
        return 1.0
    if gc < 0.40:
        return max(0.0, 1.0 - (0.40 - gc) / 0.40)
    # gc > 0.65
    return max(0.0, 1.0 - (gc - 0.65) / 0.35)


def _homopolymer_penalty(spacer: str) -> float:
    """Penalty in [0, 0.3] for the longest homopolymer run >= 4."""
    if len(spacer) < 2:
        return 0.0
    longest = 1
    cur = 1
    for i in range(1, len(spacer)):
        if spacer[i] == spacer[i - 1]:
            cur += 1
            longest = max(longest, cur)
        else:
            cur = 1
    if longest < 4:
        return 0.0
    if longest == 4:
        return 0.10
    if longest == 5:
        return 0.20
    return 0.30


def score_30mer_features(seq30: str) -> float:
    """Compute the feature-based on-target score for a 30-mer context."""
    if len(seq30) != 30:
        raise ValueError(f"on-target context must be 30 nt; got {len(seq30)}")
    seq30 = seq30.upper()
    spacer = seq30[4:24]
    pam = seq30[24:27]

    base = 0.55                                         # neutral midpoint
    gc = _gc_score(spacer)
    pam_n = _PAM_N_BONUS.get(pam[0], 0.0)
    seed = sum(_SEED_BONUS[pos].get(spacer[pos - 1], 0.0)
               for pos in (17, 18, 19, 20))
    homo = _homopolymer_penalty(spacer)

    # Weighted blend; weights chosen so the score sits in [0, 1].
    score = base + 0.30 * (gc - 0.5) + pam_n + seed - homo
    return max(0.0, min(1.0, score))


# ---------------------------------------------------------------------------
# Azimuth backend (optional)
# ---------------------------------------------------------------------------
def _score_30mer_azimuth(seq30: str) -> float:  # pragma: no cover
    import numpy as np
    import azimuth.model_comparison as az

    pred = az.predict(np.array([seq30.upper()]))
    return float(pred[0])


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
def score_30mer(seq30: str) -> float:
    backend = os.environ.get("ONTARGET_BACKEND", "features").lower()
    if backend == "azimuth":
        try:
            return _score_30mer_azimuth(seq30)
        except Exception as exc:  # pragma: no cover
            log.warning(
                "Azimuth backend failed (%s); falling back to features.", exc
            )
    return score_30mer_features(seq30)


# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------
def run(cfg: Config | None = None) -> Path:
    cfg = cfg or Config.load()
    in_path = project_root() / "data" / "processed" / "candidates.csv"
    if not in_path.exists():
        raise FileNotFoundError(
            f"{in_path} missing; run `python -m src.cli enumerate` first."
        )

    rows: List[Dict[str, str]] = []
    with in_path.open() as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            ctx = row["context_30mer"]
            score = score_30mer(ctx) if len(ctx) == 30 else 0.0
            row["on_target_score"] = f"{score:.6f}"
            rows.append(row)
    if "on_target_score" not in fieldnames:
        fieldnames.append("on_target_score")

    out_path = ensure_dir(in_path.parent) / "candidates.on_target.csv"
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote on-target scores -> %s (n=%d)", out_path, len(rows))
    return out_path


if __name__ == "__main__":  # pragma: no cover
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run()
