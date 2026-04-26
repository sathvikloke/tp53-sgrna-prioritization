"""CFD-style off-target scoring.

Two scoring modes are available:

  1) DOENCH-2016 mode (preferred). When `data/raw/cfd/mismatch_scores.csv`
     and `data/raw/cfd/pam_scores.csv` are present, this module loads the
     full Doench-2016 (CFD) mismatch and PAM matrices via
     `cfd_loader.load_weights()`. Use `scripts/download_cfd_weights.py`
     to fetch and convert them from the canonical CRISPOR distribution.

  2) APPROX mode (fallback). If the CSVs are missing/invalid, this module
     falls back to a transparent position-only approximation:
       * mismatches in the PAM-proximal seed region (positions 13-20)
         are heavily penalized;
       * mismatches in the PAM-distal region (positions 1-12) carry a
         small CFD score so off-targets there still cleave;
       * non-NGG PAMs incur an additional penalty (NAG ~ 0.3 of NGG;
         other NPNs much lower).
     The approximation is informative for sanity checking but should not
     be used for publishable specificity numbers — it lacks the per-
     nucleotide-pair (ref->alt) resolution of the real CFD matrix.

The scalar `cfd(off_target_seq, on_target_seq, pam)` returns a value in
[0, 1] where 1 == perfect match. The aggregate specificity score for a
spacer is computed by `aggregate_score()` using the CRISPOR-style
formula  100 / (100 + sum_off_target_CFDs). Call `cfd_mode()` to check
which mode is active.
"""

from __future__ import annotations

import logging
from typing import Iterable, Optional

from .cfd_loader import MismatchMatrix, PamMatrix, load_weights

log = logging.getLogger(__name__)

# Position-weighted penalty: weight at spacer position i (1-indexed
# from the 5' end of the protospacer). Higher weight => stronger cleavage
# at a mismatch (i.e., LESS penalty), so we use it as the CFD value
# directly when there is a mismatch at that position.
_POS_CFD = {
    1:  0.85, 2:  0.80, 3:  0.75, 4:  0.70,
    5:  0.65, 6:  0.60, 7:  0.55, 8:  0.50,
    9:  0.45, 10: 0.40, 11: 0.35, 12: 0.30,
    13: 0.25, 14: 0.20, 15: 0.15, 16: 0.12,
    17: 0.08, 18: 0.05, 19: 0.03, 20: 0.02,
}

# Non-NGG PAM penalty (multiplicative against the per-position product).
_PAM_PENALTY = {
    "AGG": 1.00, "CGG": 1.00, "GGG": 1.00, "TGG": 1.00,    # NGG = canonical
    "AAG": 0.30, "CAG": 0.30, "GAG": 0.30, "TAG": 0.30,    # NAG
    "ACG": 0.05, "CCG": 0.05, "GCG": 0.05, "TCG": 0.05,    # NCG
    "ATG": 0.05, "CTG": 0.05, "GTG": 0.05, "TTG": 0.05,    # NTG
}


def _pam_penalty(pam: str) -> float:
    # Only the first 3 nt are the canonical SpCas9 PAM; ignore any extras.
    return _PAM_PENALTY.get(pam.upper()[:3], 0.0)


# Try to load the published Doench-2016 weights at import time. Falls
# back to the embedded approximation if files are missing or invalid.
_MISMATCH_MATRIX: Optional[MismatchMatrix]
_PAM_MATRIX: Optional[PamMatrix]
_MISMATCH_MATRIX, _PAM_MATRIX = load_weights()
_MODE = "doench2016" if (_MISMATCH_MATRIX and _PAM_MATRIX) else "approx"
if _MODE == "approx":
    log.warning(
        "Doench-2016 CFD weights not found at data/raw/cfd/*.csv; "
        "using the position-only approximation. Run "
        "scripts/download_cfd_weights.py to upgrade."
    )
else:
    log.info("CFD scoring uses Doench-2016 weights "
             "(%d mismatch entries, %d PAM entries).",
             len(_MISMATCH_MATRIX), len(_PAM_MATRIX))


def cfd_mode() -> str:
    """Return 'doench2016' or 'approx' depending on which weights are loaded."""
    return _MODE


def cfd(off_target: str, on_target: str, pam: str) -> float:
    """Compute the CFD score for one off-target site.

    Parameters
    ----------
    off_target, on_target : str
        20-nt protospacer sequences (5'->3') for the off-target site
        and the perfect on-target spacer.
    pam : str
        The 3-nt PAM observed at the off-target site (e.g. "AGG").

    Returns
    -------
    float in [0, 1]
    """
    if len(off_target) != 20 or len(on_target) != 20:
        raise ValueError("CFD requires 20-nt protospacers")

    pam3 = pam.upper()[:3]
    if _MODE == "doench2016":
        score = _PAM_MATRIX.get(pam3, 0.0)
    else:
        score = _pam_penalty(pam3)
    if score == 0.0:
        return 0.0

    on = on_target.upper()
    off = off_target.upper()
    for i in range(20):
        ref = on[i]
        alt = off[i]
        if ref == alt or "N" in (ref, alt):
            continue
        if _MODE == "doench2016":
            # Doench 2016 keys mismatches by (position, ref, alt). If a
            # specific pair is absent (rare), fall back to the position
            # approximation rather than zeroing the score.
            w = _MISMATCH_MATRIX.get((i + 1, ref, alt))
            if w is None:
                w = _POS_CFD[i + 1]
        else:
            w = _POS_CFD[i + 1]
        score *= w
        if score == 0.0:
            return 0.0
    return score


def aggregate_score(off_target_cfds: Iterable[float]) -> float:
    """CRISPOR-style specificity: 100 / (100 + sum(CFDs)).

    Returns a value in [0, 1] (technically (0, 1]); higher = safer.
    The on-target hit itself MUST be excluded by the caller.
    """
    total = float(sum(off_target_cfds))
    return 100.0 / (100.0 + total)         # in (0, 1]


__all__ = ["cfd", "aggregate_score", "cfd_mode"]
