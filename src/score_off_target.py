"""Off-target scoring against GRCh38 via Cas-OFFinder + CFD.

Inputs:
    data/processed/candidates.on_target.csv

Outputs:
    data/processed/cas_offinder_input.txt        (intermediate)
    data/raw/cas_offinder/{sgRNA}.tsv             (per-spacer off-target raw)
    data/processed/candidates.scored.csv          (on + off scores)
    data/processed/off_target_top5.csv            (top-5 off-targets per spacer)

The user must:
    1) Have built a GRCh38 directory (see scripts/setup_grch38.sh).
    2) Have `cas-offinder` on PATH (https://github.com/snugel/cas-offinder).
    3) Set the `GRCH38_DIR` environment variable to point at that dir.
"""

from __future__ import annotations

import csv
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

from .cfd import aggregate_score, cfd
from .config import Config, TP53Locus
from .utils import ensure_dir, project_root

log = logging.getLogger(__name__)


def _check_cas_offinder() -> str:
    path = shutil.which("cas-offinder")
    if not path:
        raise FileNotFoundError(
            "`cas-offinder` not found on PATH. Install it from "
            "https://github.com/snugel/cas-offinder and ensure the binary "
            "is reachable, or run with the test-mode override."
        )
    return path


def _resolve_genome_dir(cfg: Config) -> Path:
    raw = cfg.off_target.get("genome_dir", "")
    if raw and not raw.startswith("${"):
        gpath = Path(raw)
    else:
        env = os.environ.get("GRCH38_DIR")
        if not env:
            raise EnvironmentError(
                "GRCH38_DIR environment variable is not set. "
                "Run scripts/setup_grch38.sh and set GRCH38_DIR=<that dir>."
            )
        gpath = Path(env)
    if not gpath.is_dir():
        raise NotADirectoryError(f"GRCh38 directory not found: {gpath}")
    return gpath


def _build_input_file(spacers: List[str], pam_pattern: str,
                      genome_dir: Path, max_mm: int, out: Path) -> None:
    """Write the Cas-OFFinder input file.

    Format (one line each):
        line 1: directory containing chromosome FASTA files
        line 2: search pattern (e.g. NNNNNNNNNNNNNNNNNNNNNGG)
        line 3+: <spacer><PAM-placeholder> <max_mismatches>
    The PAM placeholder is N's of length equal to the PAM portion of the
    pattern (the suffix after the trailing run of N's).
    """
    ensure_dir(out.parent)
    pam_pattern = pam_pattern.upper()
    # Strip the leading run of N's representing the protospacer to recover
    # the PAM portion (e.g. NNN...NNNNGG -> "NGG"; we then need 3 placeholder
    # N's in the query).
    n_proto = 0
    for ch in pam_pattern:
        if ch == "N":
            n_proto += 1
        else:
            break
    pam_len = len(pam_pattern) - n_proto
    if pam_len <= 0:
        raise ValueError(f"PAM portion is empty in pattern {pam_pattern!r}")
    placeholder = "N" * pam_len
    lines = [str(genome_dir), pam_pattern]
    for s in spacers:
        s = s.upper()
        if len(s) != 20:
            raise ValueError(f"spacer must be 20 nt; got {len(s)}: {s!r}")
        lines.append(f"{s}{placeholder} {max_mm}")
    out.write_text("\n".join(lines) + "\n")


def _run_cas_offinder(input_file: Path, out_file: Path,
                      device: str = "C") -> None:
    """Run cas-offinder. `device` is C (CPU), G (GPU), or A (accelerator)."""
    binary = _check_cas_offinder()
    cmd = [binary, str(input_file), device, str(out_file)]
    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def _parse_cas_offinder_output(path: Path) -> List[Dict[str, str]]:
    """Return one dict per off-target hit. Tolerant to format variants."""
    out: List[Dict[str, str]] = []
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("Bulge"):
                continue
            parts = line.split("\t")
            # Modern cas-offinder w/o bulges: query, chrom, pos, target, strand, mm
            # With bulges (8 columns), the leading two columns are bulge sizes.
            if len(parts) >= 8 and parts[0].lstrip("-").isdigit():
                _, _, query, chrom, pos, target, strand, mm = parts[:8]
            elif len(parts) >= 6:
                query, chrom, pos, target, strand, mm = parts[:6]
            else:
                continue
            out.append({
                "query": query.upper(),
                "chrom": chrom,
                "position": pos,
                "target": target.upper(),
                "strand": strand,
                "mismatches": mm,
            })
    return out


def _query_to_spacer(query_with_pam_placeholder: str,
                     spacer_length: int = 20) -> str:
    """Take the first `spacer_length` bases of a Cas-OFFinder query line.

    The Cas-OFFinder query format is `<spacer><PAM-placeholder>`; everything
    past `spacer_length` is the placeholder, regardless of its length.
    """
    return query_with_pam_placeholder.upper()[:spacer_length]


def _hit_pam(target_with_pam: str, spacer_length: int = 20) -> str:
    """The PAM at the off-target site (everything past `spacer_length`)."""
    return target_with_pam[spacer_length:].upper()


def _hit_protospacer(target_with_pam: str, spacer_length: int = 20) -> str:
    return target_with_pam[:spacer_length].upper()


def _normalize_chrom(name: str) -> str:
    """Normalize chromosome names so 'chr17' and '17' compare equal."""
    if not name:
        return ""
    n = name.strip()
    return n[3:] if n.lower().startswith("chr") else n


def _is_on_target_locus(hit: Dict[str, str], locus: TP53Locus) -> bool:
    """True iff a Cas-OFFinder hit falls inside the TP53 genomic locus.

    Coordinate-based exclusion (Flag #3 fix). Previously we excluded any hit
    with zero mismatches and an NGG PAM, which would have wrongly excluded
    matches in the TP53P1 pseudogene on chromosome 1 and would have failed
    to exclude on-target hits with non-canonical PAMs (e.g. NAG bystanders).
    """
    if _normalize_chrom(hit.get("chrom", "")) != _normalize_chrom(locus.chromosome):
        return False
    try:
        pos = int(hit.get("position", "-1"))
    except (TypeError, ValueError):
        return False
    return locus.start <= pos <= locus.end


def _score_spacer(spacer: str, hits: List[Dict[str, str]],
                  tp53_locus: TP53Locus
                  ) -> Tuple[float, List[Dict[str, str]]]:
    """Return (aggregate score in [0,1], top-5 detailed off-targets).

    `tp53_locus` is used to exclude on-target hits by genomic coordinate
    rather than by mismatch count, so that hits in the TP53P1 pseudogene on
    chromosome 1 are correctly counted as off-targets, not silently dropped.
    """
    cfds: List[Tuple[float, Dict[str, str]]] = []
    for hit in hits:
        target = hit["target"].upper()
        if len(target) < 23:                       # 20 nt + >=3 nt PAM
            continue
        # Coordinate-based on-target exclusion: skip iff the hit lies inside
        # TP53 itself. Pseudogene hits (TP53P1 on chr1) are kept and scored.
        if _is_on_target_locus(hit, tp53_locus):
            continue
        proto = _hit_protospacer(target)
        pam = _hit_pam(target)
        score = cfd(proto, spacer, pam)
        if score > 0:
            cfds.append((score, {**hit, "cfd": f"{score:.6f}"}))
    cfds.sort(key=lambda x: x[0], reverse=True)
    aggregate = aggregate_score(s for s, _ in cfds)
    top5 = [h for _, h in cfds[:5]]
    return aggregate, top5


def run(cfg: Config | None = None, *, device: str = "C") -> Path:
    cfg = cfg or Config.load()
    in_path = project_root() / "data" / "processed" / "candidates.on_target.csv"
    if not in_path.exists():
        raise FileNotFoundError(
            f"{in_path} missing; run on-target scoring first."
        )

    rows: List[Dict[str, str]] = []
    with in_path.open() as fh:
        rows = list(csv.DictReader(fh))

    if not rows:
        raise RuntimeError("No candidate sgRNAs to score.")

    genome_dir = _resolve_genome_dir(cfg)
    pam_pattern = cfg.off_target["pam_pattern"]
    max_mm = int(cfg.off_target["max_mismatches"])

    work_dir = ensure_dir(project_root() / "data" / "processed")
    raw_dir = ensure_dir(project_root() / "data" / "raw" / "cas_offinder")
    in_file = work_dir / "cas_offinder_input.txt"
    out_file = raw_dir / "all_hits.tsv"

    spacers = [r["spacer_5p_to_3p"] for r in rows]
    unique_spacers = sorted(set(spacers))
    _build_input_file(unique_spacers, pam_pattern, genome_dir, max_mm, in_file)
    _run_cas_offinder(in_file, out_file, device=device)
    hits = _parse_cas_offinder_output(out_file)

    # Index hits by spacer (the query field with the trailing NNN stripped).
    by_spacer: Dict[str, List[Dict[str, str]]] = {}
    for h in hits:
        s = _query_to_spacer(h["query"])
        by_spacer.setdefault(s, []).append(h)

    # Score each candidate.
    fieldnames = list(rows[0].keys())
    if "off_target_score" not in fieldnames:
        fieldnames += ["off_target_score", "n_off_targets_total"]

    top5_records: List[Dict[str, str]] = []
    for r in rows:
        spacer = r["spacer_5p_to_3p"].upper()
        spacer_hits = by_spacer.get(spacer, [])
        agg, top5 = _score_spacer(spacer, spacer_hits, cfg.tp53_locus)
        # Off-targets = all Cas-OFFinder hits that did NOT land inside the
        # TP53 locus (Flag #3 fix). This handles pseudogene hits correctly.
        n_off = sum(1 for h in spacer_hits
                    if not _is_on_target_locus(h, cfg.tp53_locus))
        r["off_target_score"] = f"{agg:.6f}"
        r["n_off_targets_total"] = str(n_off)
        for h in top5:
            top5_records.append({"sgRNA_id": r["sgRNA_id"], **h})

    out_csv = work_dir / "candidates.scored.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote off-target scores -> %s (n=%d)", out_csv, len(rows))

    top5_csv = work_dir / "off_target_top5.csv"
    with top5_csv.open("w", newline="") as fh:
        if top5_records:
            w = csv.DictWriter(fh, fieldnames=list(top5_records[0].keys()))
            w.writeheader()
            w.writerows(top5_records)
        else:
            fh.write("sgRNA_id,query,chrom,position,target,strand,mismatches,cfd\n")
    log.info("Wrote top-5 off-targets per spacer -> %s", top5_csv)
    return out_csv


if __name__ == "__main__":  # pragma: no cover
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run()
