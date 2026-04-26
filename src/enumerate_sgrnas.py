"""Enumerate candidate sgRNAs around each TP53 hotspot residue.

Cut-site convention (SpCas9):
  - For a sense-strand sgRNA whose PAM (NGG) starts at CDS position p
    (0-indexed), the protospacer is CDS[p-20 : p] and the cut occurs
    between CDS positions p-4 and p-3 (i.e., 3 nt 5' of the PAM,
    between bp 17 and 18 of the protospacer counted 1-indexed from
    the 5' end). We report `cut_position = p - 3`.

  - For an antisense-strand sgRNA whose PAM is "NGG" on the reverse-
    complement, equivalently a "CCN" pattern at CDS position q on the
    sense CDS:
        protospacer (5'->3') = revcomp(CDS[q+3 : q+23])
        cut occurs between CDS positions q+5 and q+6;
        we report `cut_position = q + 6`.

The edit-window check is:
        |cut_position - codon_start| <= edit_window_nt

where `codon_start` is the 0-indexed first base of the hotspot codon.

CAVEAT: PAMs and protospacers are searched within the CDS only. Hotspot
codons of interest (residues 175, 220, 245, 248, 249, 273, 282) sit
several codons inside their respective TP53 exons, so the +/- 10 nt
window does not cross exon boundaries; the CDS is therefore a faithful
proxy of genomic DNA at those windows. For other residues this assumption
must be re-checked against the genomic context.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from .config import Config, Hotspot
from .utils import (
    codon_index_for_residue,
    ensure_dir,
    project_root,
    reverse_complement,
)

log = logging.getLogger(__name__)


CONTEXT_5 = 4   # nt 5' of protospacer needed for the 30-mer on-target context
CONTEXT_3 = 3   # nt 3' of PAM for the 30-mer (4 + 20 + 3 + 3 = 30)


def _iter_sense_pams(cds: str) -> Iterable[int]:
    """Yield CDS positions p where CDS[p..p+2] matches NGG (sense strand)."""
    for p in range(len(cds) - 2):
        if cds[p + 1] == "G" and cds[p + 2] == "G":
            yield p


def _iter_antisense_pams(cds: str) -> Iterable[int]:
    """Yield CDS positions q where CDS[q..q+2] = CCN (antisense PAM)."""
    for q in range(len(cds) - 2):
        if cds[q] == "C" and cds[q + 1] == "C":
            yield q


def _build_record(*, residue: int, ref_aa: str, label: str, strand: str,
                  spacer: str, pam: str, context30: str,
                  cut_position: int, codon_start: int,
                  pam_position_in_cds: int) -> Dict[str, object]:
    return {
        "residue": residue,
        "ref_aa": ref_aa,
        "hotspot_label": label,
        "strand": strand,
        "spacer_5p_to_3p": spacer,
        "pam": pam,
        "context_30mer": context30,
        "cut_position_cds_0based": cut_position,
        "codon_start_cds_0based": codon_start,
        "cut_offset_from_codon_start_nt": cut_position - codon_start,
        "pam_position_cds_0based": pam_position_in_cds,
    }


def _crosses_exon_boundary(start: int, end: int,
                           boundaries: List[int]) -> bool:
    """True iff any boundary in `boundaries` lies inside (start, end).

    `start` and `end` are 0-based CDS positions; the half-open interval
    [start, end) is what's checked. Empty `boundaries` -> always False
    (the check is opt-in via config).
    """
    for b in boundaries:
        if start < b < end:
            return True
    return False


def candidates_for_hotspot(cds: str, hotspot: Hotspot,
                           edit_window_nt: int,
                           spacer_length: int,
                           cds_exon_boundaries: List[int] | None = None,
                           ) -> List[Dict[str, object]]:
    """Return all sense + antisense candidates whose cut lands in window.

    If `cds_exon_boundaries` is provided (non-empty), each candidate is
    tagged with `crosses_exon_boundary` (True/False) reflecting whether
    its protospacer + 30-mer context spans a CDS exon boundary. Such
    candidates are still emitted but should be reviewed before scoring
    (Flag #8 safety check).
    """
    boundaries = cds_exon_boundaries or []
    codon_start, _ = codon_index_for_residue(hotspot.residue)
    candidates: List[Dict[str, object]] = []

    # ---- sense strand: NGG PAM at p ----
    for p in _iter_sense_pams(cds):
        cut_position = p - 3                          # 0-indexed
        if abs(cut_position - codon_start) > edit_window_nt:
            continue
        spacer_start = p - spacer_length
        if spacer_start < 0:
            continue
        ctx_start = spacer_start - CONTEXT_5
        ctx_end = p + 3 + CONTEXT_3                  # PAM is 3 nt
        if ctx_start < 0 or ctx_end > len(cds):
            continue
        spacer = cds[spacer_start:p]
        pam = cds[p:p + 3]
        context30 = cds[ctx_start:ctx_end]
        if len(context30) != 30:
            continue
        rec = _build_record(
            residue=hotspot.residue, ref_aa=hotspot.ref_aa,
            label=hotspot.label, strand="+",
            spacer=spacer, pam=pam, context30=context30,
            cut_position=cut_position, codon_start=codon_start,
            pam_position_in_cds=p,
        )
        rec["crosses_exon_boundary"] = (
            _crosses_exon_boundary(ctx_start, ctx_end, boundaries)
        )
        candidates.append(rec)

    # ---- antisense strand: CCN PAM at q on CDS ----
    for q in _iter_antisense_pams(cds):
        cut_position = q + 6
        if abs(cut_position - codon_start) > edit_window_nt:
            continue
        proto_start = q + 3
        proto_end = q + 3 + spacer_length
        if proto_end > len(cds):
            continue
        # 30-mer context on sense strand spans 3 nt left of CCN through
        # 4 nt right of protospacer end, then reverse-complemented.
        ctx_start = q - CONTEXT_3                     # 3 nt 5' of CCN on sense
        ctx_end = proto_end + CONTEXT_5
        if ctx_start < 0 or ctx_end > len(cds):
            continue
        sense_protospacer = cds[proto_start:proto_end]
        spacer = reverse_complement(sense_protospacer)
        # PAM on the antisense strand (5'->3'): revcomp of CDS[q:q+3].
        pam = reverse_complement(cds[q:q + 3])
        context30 = reverse_complement(cds[ctx_start:ctx_end])
        if len(context30) != 30:
            continue
        rec = _build_record(
            residue=hotspot.residue, ref_aa=hotspot.ref_aa,
            label=hotspot.label, strand="-",
            spacer=spacer, pam=pam, context30=context30,
            cut_position=cut_position, codon_start=codon_start,
            pam_position_in_cds=q,
        )
        rec["crosses_exon_boundary"] = (
            _crosses_exon_boundary(ctx_start, ctx_end, boundaries)
        )
        candidates.append(rec)

    return candidates


def assign_ids(candidates: List[Dict[str, object]]) -> None:
    """Assign deterministic sgRNA IDs (in-place).

    ID format: TP53-<residue>-<strand>-<spacer-md5-prefix>
    """
    from .utils import md5_str
    for c in candidates:
        token = f"{c['spacer_5p_to_3p']}|{c['strand']}"
        c["sgRNA_id"] = (
            f"TP53-{c['residue']}-{'P' if c['strand']=='+' else 'M'}-"
            f"{md5_str(token)[:8]}"
        )


def run(cfg: Config | None = None) -> Path:
    cfg = cfg or Config.load()
    cds_path = project_root() / "data" / "raw" / f"{cfg.refseq_accession}.cds.txt"
    if not cds_path.exists():
        raise FileNotFoundError(
            f"{cds_path} missing; run `python -m src.cli fetch-cds` first."
        )
    cds = cds_path.read_text().strip().upper()

    boundaries = cfg.tp53_locus.cds_exon_boundaries_0based
    all_rows: List[Dict[str, object]] = []
    for h in cfg.hotspots:
        rows = candidates_for_hotspot(
            cds, h,
            edit_window_nt=cfg.edit_window_nt,
            spacer_length=cfg.spacer_length,
            cds_exon_boundaries=boundaries,
        )
        n_boundary = sum(1 for r in rows if r.get("crosses_exon_boundary"))
        if n_boundary:
            log.warning(
                "residue %d (%s): %d/%d candidates straddle a CDS exon "
                "boundary; review their context_30mer before scoring "
                "(genomic context may not equal CDS context).",
                h.residue, h.label, n_boundary, len(rows),
            )
        log.info("residue %d (%s): %d candidates", h.residue, h.label, len(rows))
        all_rows.extend(rows)

    assign_ids(all_rows)

    out_dir = ensure_dir(project_root() / "data" / "processed")
    out_path = out_dir / "candidates.csv"
    fieldnames = [
        "sgRNA_id", "residue", "ref_aa", "hotspot_label", "strand",
        "spacer_5p_to_3p", "pam", "context_30mer",
        "pam_position_cds_0based", "cut_position_cds_0based",
        "codon_start_cds_0based", "cut_offset_from_codon_start_nt",
        "crosses_exon_boundary",
    ]
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in all_rows:
            w.writerow({k: r[k] for k in fieldnames})
    log.info("Wrote %d candidates -> %s", len(all_rows), out_path)
    return out_path


if __name__ == "__main__":  # pragma: no cover
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run()
