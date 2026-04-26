"""Load Doench-2016 CFD weights from local CSV files (Flag #2).

The CFD off-target score is parameterized by:
    * a "mismatch matrix" of CFD penalties keyed on
      (position 1..20, ref_nt, alt_nt)
    * a "PAM matrix" of CFD penalties keyed on the trinucleotide PAM

This module reads the matrices from:
    data/raw/cfd/mismatch_scores.csv
    data/raw/cfd/pam_scores.csv

If both files are present and parse cleanly, `load_weights()` returns the
matrices. If either is missing or invalid, it returns (None, None) and
the caller (cfd.py) falls back to the embedded position-only
approximation. This keeps the pipeline runnable in environments without
network access while letting users opt into the published weights.

CSV formats
-----------
mismatch_scores.csv — header: position,ref,alt,score
    one row per (position, ref, alt) triple where ref != alt; ref/alt are
    one of A/C/G/T; position is 1-indexed (1..20); score is in [0, 1].

pam_scores.csv — header: pam,score
    one row per 3-nt PAM (uppercase). Only PAMs present in the file are
    accepted; missing PAMs default to 0.0 in the consumer.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple

from .utils import project_root

log = logging.getLogger(__name__)


MismatchKey = Tuple[int, str, str]
MismatchMatrix = Dict[MismatchKey, float]
PamMatrix = Dict[str, float]


def _cfd_dir() -> Path:
    return project_root() / "data" / "raw" / "cfd"


def _load_mismatch(path: Path) -> Optional[MismatchMatrix]:
    out: MismatchMatrix = {}
    valid_nt = set("ACGT")
    with path.open() as fh:
        reader = csv.DictReader(fh)
        required = {"position", "ref", "alt", "score"}
        if not required.issubset(reader.fieldnames or []):
            log.warning("mismatch_scores.csv missing columns; got %s",
                        reader.fieldnames)
            return None
        for row in reader:
            try:
                pos = int(row["position"])
                ref = row["ref"].strip().upper()
                alt = row["alt"].strip().upper()
                score = float(row["score"])
            except (ValueError, KeyError) as exc:
                log.warning("mismatch_scores.csv row %r invalid: %s", row, exc)
                return None
            if pos < 1 or pos > 20:
                log.warning("mismatch_scores.csv: position %d out of [1,20]", pos)
                return None
            if ref not in valid_nt or alt not in valid_nt:
                log.warning("mismatch_scores.csv: bad nt ref=%r alt=%r",
                            ref, alt)
                return None
            if not (0.0 <= score <= 1.0):
                log.warning("mismatch_scores.csv: score %f out of [0,1]", score)
                return None
            out[(pos, ref, alt)] = score
    if not out:
        return None
    return out


def _load_pam(path: Path) -> Optional[PamMatrix]:
    out: PamMatrix = {}
    with path.open() as fh:
        reader = csv.DictReader(fh)
        required = {"pam", "score"}
        if not required.issubset(reader.fieldnames or []):
            log.warning("pam_scores.csv missing columns; got %s",
                        reader.fieldnames)
            return None
        for row in reader:
            try:
                pam = row["pam"].strip().upper()
                score = float(row["score"])
            except (ValueError, KeyError) as exc:
                log.warning("pam_scores.csv row %r invalid: %s", row, exc)
                return None
            if len(pam) != 3 or any(c not in "ACGT" for c in pam):
                log.warning("pam_scores.csv: invalid PAM %r", pam)
                return None
            if not (0.0 <= score <= 1.0):
                log.warning("pam_scores.csv: score %f out of [0,1]", score)
                return None
            out[pam] = score
    if not out:
        return None
    return out


def load_weights() -> Tuple[Optional[MismatchMatrix], Optional[PamMatrix]]:
    """Return (mismatch_matrix, pam_matrix) or (None, None) if unavailable."""
    d = _cfd_dir()
    mm_path = d / "mismatch_scores.csv"
    pam_path = d / "pam_scores.csv"
    if not mm_path.exists() or not pam_path.exists():
        return None, None
    try:
        mm = _load_mismatch(mm_path)
        pm = _load_pam(pam_path)
    except Exception as exc:                            # pragma: no cover
        log.warning("Failed to load CFD weights: %s", exc)
        return None, None
    if mm is None or pm is None:
        return None, None
    return mm, pm


__all__ = ["load_weights", "MismatchMatrix", "PamMatrix"]
