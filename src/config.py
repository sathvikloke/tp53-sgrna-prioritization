"""Load + validate `config/params.yaml`."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping

import yaml

from .utils import project_root


_ENV_PATTERN = re.compile(r"\$\{([A-Z_][A-Z0-9_]*)\}")


def _expand_env(value: Any) -> Any:
    """Recursively expand ${VAR} substrings in strings using os.environ."""
    if isinstance(value, str):
        def repl(m: "re.Match[str]") -> str:
            var = m.group(1)
            return os.environ.get(var, m.group(0))
        return _ENV_PATTERN.sub(repl, value)
    if isinstance(value, list):
        return [_expand_env(v) for v in value]
    if isinstance(value, dict):
        return {k: _expand_env(v) for k, v in value.items()}
    return value


@dataclass
class Hotspot:
    residue: int
    ref_aa: str
    label: str


@dataclass
class Cohort:
    id: str
    label: str


@dataclass
class TP53Locus:
    chromosome: str
    start: int                  # 1-based inclusive
    end: int                    # 1-based inclusive
    strand: str
    exons_grch38: List[Dict[str, int]]
    # CDS-coordinate (0-based) positions where one CDS exon ends and the
    # next begins. Used by enumerate_sgrnas to flag any candidate whose
    # edit window straddles an exon boundary (Flag #8). Empty list disables
    # the check (no false alarms when boundaries are not configured).
    cds_exon_boundaries_0based: List[int] = field(default_factory=list)


@dataclass
class Config:
    raw: Dict[str, Any]
    cas_enzyme: str
    pam: str
    spacer_length: int
    cut_offset_from_pam: int
    pam_rescue: Dict[str, str]
    edit_window_nt: int
    framing_mode: str
    hotspots: List[Hotspot]
    hotspot_variants: List[str]
    cohorts: List[Cohort]
    refseq_accession: str
    expected_cds_length: int
    ncbi_efetch_url: str
    tp53_locus: TP53Locus
    cbioportal_base_url: str
    tp53_entrez_id: int
    weights: Dict[str, float]
    off_target: Dict[str, Any]

    @classmethod
    def load(cls, path: Path | None = None) -> "Config":
        path = Path(path) if path else project_root() / "config" / "params.yaml"
        if not path.exists():
            raise FileNotFoundError(f"Config not found: {path}")
        with path.open() as f:
            raw = yaml.safe_load(f)
        raw = _expand_env(raw)

        cas = raw["cas"]
        ref = raw["reference"]
        cbio = raw["cbioportal"]
        scoring = raw["scoring"]
        framing = raw.get("framing", {"mode": "hdr", "edit_window_nt":
                                       raw.get("edit_window_nt", 10)})
        locus_raw = raw["tp53_locus"]

        weights = {k: float(v) for k, v in scoring["weights"].items()}
        wsum = sum(weights.values())
        if abs(wsum - 1.0) > 1e-6:
            raise ValueError(
                f"Composite weights must sum to 1.0; got {wsum:.6f} ({weights})"
            )

        if framing["mode"] not in {"hdr"}:
            raise ValueError(
                f"Unsupported framing mode {framing['mode']!r}. "
                "This pipeline currently supports only 'hdr'; "
                "knockout / base_editing variants would require additional "
                "constraints not yet implemented."
            )

        edit_window = int(framing.get("edit_window_nt",
                                       raw.get("edit_window_nt", 10)))
        tp53_locus = TP53Locus(
            chromosome=str(locus_raw["chromosome"]),
            start=int(locus_raw["start"]),
            end=int(locus_raw["end"]),
            strand=str(locus_raw["strand"]),
            exons_grch38=[
                {"number": int(e["number"]), "start": int(e["start"]),
                 "end": int(e["end"])}
                for e in locus_raw["exons_grch38"]
            ],
            cds_exon_boundaries_0based=[
                int(b) for b in locus_raw.get("cds_exon_boundaries_0based", [])
            ],
        )

        return cls(
            raw=raw,
            cas_enzyme=cas["enzyme"],
            pam=cas["pam"],
            spacer_length=int(cas["spacer_length"]),
            cut_offset_from_pam=int(cas["cut_offset_from_pam"]),
            pam_rescue=dict(cas.get("pam_rescue", {})),
            edit_window_nt=edit_window,
            framing_mode=framing["mode"],
            hotspots=[Hotspot(**h) for h in raw["hotspots"]],
            hotspot_variants=list(raw["hotspot_variants"]),
            cohorts=[Cohort(**c) for c in raw["cohorts"]],
            refseq_accession=ref["refseq_accession"],
            expected_cds_length=int(ref["expected_cds_length"]),
            ncbi_efetch_url=ref["ncbi_efetch_url"],
            tp53_locus=tp53_locus,
            cbioportal_base_url=cbio["base_url"],
            tp53_entrez_id=int(cbio["tp53_entrez_id"]),
            weights=weights,
            off_target=raw["off_target"],
        )


__all__ = ["Config", "Hotspot", "Cohort", "TP53Locus"]
