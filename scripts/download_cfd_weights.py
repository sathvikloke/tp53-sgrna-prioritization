#!/usr/bin/env python3
"""Download Doench-2016 CFD weights and convert them to CSV (Flag #2).

Source
------
The canonical distribution of the Doench-2016 CFD scoring tables is the
CRISPOR repository (https://github.com/maximilianh/crisporWebsite). It
ships two Python pickle files under cfd-scoring/:

    mismatch_score.pkl   keys: "rNxN" (e.g. "rGxA"), values: per-position
                         lists of length 20 of CFD penalties in [0, 1].
    pam_scores.pkl       keys: 3-nt PAMs, values: CFD penalty in [0, 1].

This script downloads the pickles from a pinned commit, verifies their
SHA256 hashes against a list the user controls, loads them, and writes
clean CSVs to data/raw/cfd/{mismatch_scores.csv, pam_scores.csv}. The
CSVs are what cfd_loader.load_weights() consumes at runtime.

WHY a pinned commit + hash check?
---------------------------------
We're loading pickles, which is unsafe for untrusted data. By pinning to
a specific commit and verifying the SHA256 BEFORE loading, the user is
explicitly attesting to which file they trust. If the upstream file
changes, the hash check fails and the pickle is never loaded.

Usage
-----
    python scripts/download_cfd_weights.py

To pin a different commit, edit the COMMIT and EXPECTED_SHA256 constants
below and re-run. To use a local mirror, pass `--mirror /path/to/dir/`.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import pickle  # nosec: B403 — guarded by SHA256 verification below
import sys
import urllib.request
from pathlib import Path
from typing import Any, Dict


# Pinned CRISPOR commit. Update this AND the hashes together.
COMMIT = "61bf3ff6cb9d4c6c8b75d9e5e44e95e62f0f5eb6"
BASE_URL = (
    "https://raw.githubusercontent.com/maximilianh/crisporWebsite/"
    f"{COMMIT}/cfd-scoring"
)
FILES = {
    "mismatch_score.pkl": "mismatch_score.pkl",
    "pam_scores.pkl":     "pam_scores.pkl",
}

# Expected SHA256 hashes. Populate these from a trusted source the first
# time you run the script (the script will print the hashes it sees if
# the entry is None, allowing you to record them after manual review).
EXPECTED_SHA256: Dict[str, str | None] = {
    "mismatch_score.pkl": None,   # paste the SHA256 you trust here
    "pam_scores.pkl":     None,
}


log = logging.getLogger("download_cfd_weights")


def _fetch(url: str) -> bytes:
    log.info("GET %s", url)
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read()


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _verify_or_print(name: str, data: bytes) -> bool:
    h = _sha256(data)
    expected = EXPECTED_SHA256.get(name)
    if expected is None:
        print(
            f"[WARN] No expected SHA256 pinned for {name}.\n"
            f"       Observed: {h}\n"
            f"       If you trust this commit, paste the hash into "
            f"EXPECTED_SHA256 in this script and re-run.",
            file=sys.stderr,
        )
        return False
    if h != expected:
        print(
            f"[ERROR] SHA256 mismatch for {name}.\n"
            f"        Expected: {expected}\n"
            f"        Observed: {h}\n"
            f"        Refusing to load this file as a pickle.",
            file=sys.stderr,
        )
        return False
    return True


def _convert_mismatch(obj: Any, out_path: Path) -> None:
    """Convert CRISPOR mismatch_score.pkl -> our mismatch_scores.csv.

    CRISPOR format (v1): dict[str, list[float]] where key is "rNxN"
    (e.g. "rGxA") meaning a guide nt = N1 paired against an off-target
    DNA nt = N2. The list has length 20 with the per-position penalty.
    """
    if not isinstance(obj, dict):
        raise ValueError("mismatch_score.pkl: expected dict at top level")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["position", "ref", "alt", "score"])
        n_rows = 0
        for key, vals in obj.items():
            # key examples: "rGxA" (RNA G against DNA A), "rTxC" (rare;
            # CRISPOR uses uracil-equivalent encoding). We flatten to
            # ref=DNA-of-RNA, alt=off-target-DNA.
            if not isinstance(key, str) or len(key) < 4 or "x" not in key:
                continue
            try:
                rna_part, dna_part = key.split("x")
                ref_rna = rna_part[-1].upper()
                alt_dna = dna_part[-1].upper()
            except Exception:
                continue
            ref = "T" if ref_rna == "U" else ref_rna
            for pos, score in enumerate(vals, start=1):
                if pos < 1 or pos > 20:
                    continue
                w.writerow([pos, ref, alt_dna, f"{float(score):.6f}"])
                n_rows += 1
        log.info("Wrote %d mismatch entries -> %s", n_rows, out_path)


def _convert_pam(obj: Any, out_path: Path) -> None:
    if not isinstance(obj, dict):
        raise ValueError("pam_scores.pkl: expected dict at top level")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pam", "score"])
        n_rows = 0
        for pam, score in obj.items():
            if not isinstance(pam, str):
                continue
            pam = pam.upper().strip()
            if len(pam) != 3 or any(c not in "ACGT" for c in pam):
                continue
            w.writerow([pam, f"{float(score):.6f}"])
            n_rows += 1
        log.info("Wrote %d PAM entries -> %s", n_rows, out_path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mirror", type=Path,
                        help="Local directory containing the two .pkl files; "
                             "if set, files are read from disk instead of HTTP.")
    parser.add_argument("--out", type=Path,
                        default=Path("data/raw/cfd"),
                        help="Output directory (default: data/raw/cfd).")
    parser.add_argument("--allow-unpinned", action="store_true",
                        help="Skip SHA256 enforcement (NOT recommended).")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    payloads: Dict[str, bytes] = {}
    for fname in FILES:
        if args.mirror:
            data = (args.mirror / fname).read_bytes()
        else:
            data = _fetch(f"{BASE_URL}/{fname}")
        ok = _verify_or_print(fname, data)
        if not ok and not args.allow_unpinned:
            log.error("Aborting; pin a SHA256 for %s before loading the pickle.", fname)
            return 2
        payloads[fname] = data

    for fname, data in payloads.items():
        obj = pickle.loads(data)  # nosec: B301 — guarded by SHA256 verification
        if fname == "mismatch_score.pkl":
            _convert_mismatch(obj, out_dir / "mismatch_scores.csv")
        elif fname == "pam_scores.pkl":
            _convert_pam(obj, out_dir / "pam_scores.csv")
    log.info("Done. CFD weights written to %s/.", out_dir)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
