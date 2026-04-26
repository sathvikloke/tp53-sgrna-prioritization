"""Shared utilities: sequence helpers, codon mapping, I/O helpers."""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Tuple

# IUPAC complements (incl. ambiguity).
_COMPLEMENT: Mapping[str, str] = {
    "A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D",
}

# Standard genetic code (DNA).
_CODON_TABLE: Mapping[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------
def reverse_complement(seq: str) -> str:
    """Return the reverse complement of `seq` (uppercased)."""
    seq = seq.upper()
    try:
        return "".join(_COMPLEMENT[b] for b in reversed(seq))
    except KeyError as exc:  # pragma: no cover - defensive
        raise ValueError(f"Unknown base in sequence: {exc.args[0]!r}") from exc


def translate(cds: str) -> str:
    """Translate a CDS (must be a multiple of 3) into a protein string."""
    cds = cds.upper()
    if len(cds) % 3 != 0:
        raise ValueError(
            f"CDS length {len(cds)} is not a multiple of 3; cannot translate."
        )
    protein: List[str] = []
    for i in range(0, len(cds), 3):
        codon = cds[i : i + 3]
        protein.append(_CODON_TABLE.get(codon, "X"))
    return "".join(protein)


def gc_content(seq: str) -> float:
    """Return the GC fraction of `seq` (in [0, 1])."""
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = sum(1 for b in seq if b in ("G", "C"))
    return gc / len(seq)


def codon_index_for_residue(residue_1based: int) -> Tuple[int, int]:
    """Return the 0-indexed CDS [start, end) of the codon for residue N."""
    if residue_1based < 1:
        raise ValueError("residue must be 1-indexed and >= 1")
    start = (residue_1based - 1) * 3
    return start, start + 3


def matches_pam(triplet: str, pam_pattern: str) -> bool:
    """Test whether `triplet` matches the IUPAC pattern `pam_pattern`.

    The pattern uses standard IUPAC codes: A/C/G/T literal, N = any,
    R = A/G, Y = C/T, S = G/C, W = A/T, K = G/T, M = A/C, V = ACG,
    H = ACT, D = AGT, B = CGT.
    """
    if len(triplet) != len(pam_pattern):
        return False
    iupac: Mapping[str, set] = {
        "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
        "N": {"A", "C", "G", "T"},
        "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
        "K": {"G", "T"}, "M": {"A", "C"},
        "V": {"A", "C", "G"}, "H": {"A", "C", "T"},
        "D": {"A", "G", "T"}, "B": {"C", "G", "T"},
    }
    triplet = triplet.upper()
    pam_pattern = pam_pattern.upper()
    for actual, allowed_code in zip(triplet, pam_pattern):
        if actual not in iupac.get(allowed_code, set()):
            return False
    return True


# ---------------------------------------------------------------------------
# I/O + provenance helpers
# ---------------------------------------------------------------------------
def utc_timestamp() -> str:
    """ISO-8601 UTC timestamp with second resolution and explicit `Z`."""
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def md5_str(text: str) -> str:
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def write_json(path: Path, data: Mapping) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True))


def read_json(path: Path):
    return json.loads(Path(path).read_text())


def project_root() -> Path:
    """Return the project root (one level above this `src/` package)."""
    return Path(__file__).resolve().parent.parent


def ensure_dir(path: Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def lazy_iter_kmers(seq: str, k: int) -> Iterable[Tuple[int, str]]:
    """Yield `(start_index, kmer)` pairs of every k-mer in `seq`."""
    n = len(seq)
    for i in range(0, n - k + 1):
        yield i, seq[i : i + k]


__all__ = [
    "reverse_complement",
    "translate",
    "gc_content",
    "codon_index_for_residue",
    "matches_pam",
    "utc_timestamp",
    "md5_str",
    "write_json",
    "read_json",
    "project_root",
    "ensure_dir",
    "lazy_iter_kmers",
]
