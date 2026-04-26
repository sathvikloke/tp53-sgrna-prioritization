"""Fetch TP53 mutation calls per cohort from cBioPortal.

For each study we:
    1) GET /molecular-profiles/{study}_mutations               (verify exists)
    2) GET /sample-lists/{study}_sequenced                     (verify exists)
    3) POST /molecular-profiles/{study}_mutations/mutations/fetch
       body: {"entrezGeneIds": [7157], "sampleListId": "{study}_sequenced"}
    4) Parse mutations and count hotspot variants.

Outputs:
    data/raw/cbioportal/{study_id}.mutations.json     (raw response)
    data/raw/cbioportal/{study_id}.meta.json          (provenance)
    data/processed/cohort_frequencies.csv             (residue-level counts)
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping

import requests

from .config import Config
from .utils import ensure_dir, project_root, utc_timestamp, write_json

log = logging.getLogger(__name__)


# Published TP53 mutation prevalence floors per cohort (Flag #9 sanity check).
# These are deliberately conservative (well below the canonical mid-point so we
# only fire on serious data drift, not on minor methodological differences).
#   HGSOC: typically 90-96% (Cancer Genome Atlas Research Network, Nature 2011);
#          floor 0.85 catches dataset replacement / endpoint drift.
#   PDAC:  typically 60-75% (PAAD TCGA, ICGC); floor 0.50.
#   CRC:   typically 50-60% (TCGA COAD/READ); floor 0.40.
# If a cohort here is missing, the check is skipped for that cohort with INFO.
TP53_PREVALENCE_SENTINELS: Dict[str, float] = {
    "ov_tcga_pan_can_atlas_2018":       0.85,
    "paad_tcga_pan_can_atlas_2018":     0.50,
    "coadread_tcga_pan_can_atlas_2018": 0.40,
}


def _distinct_tp53_mutated_samples(mutations: List[Dict[str, Any]]) -> int:
    """Count distinct samples carrying ANY TP53 mutation in the response."""
    seen: set = set()
    for m in mutations:
        sid = m.get("sampleId") or m.get("uniqueSampleKey")
        if sid:
            seen.add(sid)
    return len(seen)


def _post_mutations(base_url: str, study_id: str, entrez_id: int,
                    *, timeout: int = 60) -> List[Dict[str, Any]]:
    """POST to /molecular-profiles/.../mutations/fetch and return parsed JSON."""
    url = f"{base_url}/molecular-profiles/{study_id}_mutations/mutations/fetch"
    payload = {
        "entrezGeneIds": [entrez_id],
        "sampleListId": f"{study_id}_sequenced",
    }
    params = {"projection": "DETAILED"}
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    log.info("POST %s (%s)", url, study_id)
    resp = requests.post(url, json=payload, params=params,
                         headers=headers, timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def _get_n_sequenced(base_url: str, study_id: str, *, timeout: int = 30) -> int:
    """Return the number of sequenced samples in the study."""
    url = f"{base_url}/sample-lists/{study_id}_sequenced"
    resp = requests.get(url, timeout=timeout,
                        headers={"Accept": "application/json"})
    resp.raise_for_status()
    body = resp.json()
    n = body.get("sampleCount")
    if n is None and "sampleIds" in body:
        n = len(body["sampleIds"])
    if not isinstance(n, int):
        raise ValueError(
            f"Could not determine sample count for {study_id} from {url}"
        )
    return n


def _count_hotspots(mutations: List[Dict[str, Any]],
                    hotspot_variants: List[str]) -> Dict[str, set]:
    """Return {variant_label: set(sample_ids)} restricted to the panel."""
    panel = set(hotspot_variants)
    out: Dict[str, set] = defaultdict(set)
    for m in mutations:
        variant = (m.get("proteinChange") or "").strip()
        # cBioPortal sometimes returns "p.R175H"; strip the "p." prefix.
        if variant.startswith("p."):
            variant = variant[2:]
        if variant in panel:
            sample_id = m.get("sampleId") or m.get("uniqueSampleKey")
            if sample_id:
                out[variant].add(sample_id)
    return out


def fetch_one_cohort(study_id: str, cfg: Config) -> Dict[str, Any]:
    """Fetch + summarize one cohort. Persist raw + meta. Return summary dict."""
    raw_dir = ensure_dir(project_root() / "data" / "raw" / "cbioportal")

    # Sample count (denominator) + raw mutations.
    n_sequenced = _get_n_sequenced(cfg.cbioportal_base_url, study_id)
    mutations = _post_mutations(
        cfg.cbioportal_base_url, study_id, cfg.tp53_entrez_id
    )

    raw_path = raw_dir / f"{study_id}.mutations.json"
    write_json(raw_path, {"mutations": mutations})

    meta_path = raw_dir / f"{study_id}.meta.json"
    write_json(meta_path, {
        "study_id": study_id,
        "entrez_gene_id": cfg.tp53_entrez_id,
        "fetched_at_utc": utc_timestamp(),
        "n_mutations_returned": len(mutations),
        "n_sequenced_samples": n_sequenced,
        "endpoint": (
            f"{cfg.cbioportal_base_url}/molecular-profiles/"
            f"{study_id}_mutations/mutations/fetch"
        ),
    })

    counts = _count_hotspots(mutations, cfg.hotspot_variants)

    # Sanity-check (Flag #9): TP53 prevalence floor per cohort.
    n_tp53_mut = _distinct_tp53_mutated_samples(mutations)
    overall_prev = (n_tp53_mut / n_sequenced) if n_sequenced else 0.0
    floor = TP53_PREVALENCE_SENTINELS.get(study_id)
    if floor is None:
        log.info(
            "[%s] no TP53 prevalence sentinel configured; skipping sanity check.",
            study_id,
        )
    elif overall_prev < floor:
        log.warning(
            "[%s] TP53 mutation prevalence %.1f%% (%d/%d) is below the "
            "published floor of %.0f%%. Possible API drift, study-id "
            "rename, or endpoint change. Inspect the raw response before "
            "trusting cohort_frequencies.csv.",
            study_id, overall_prev * 100, n_tp53_mut, n_sequenced, floor * 100,
        )
    else:
        log.info(
            "[%s] TP53 prevalence %.1f%% (%d/%d) >= floor %.0f%% — OK.",
            study_id, overall_prev * 100, n_tp53_mut, n_sequenced, floor * 100,
        )

    summary: Dict[str, Any] = {
        "study_id": study_id,
        "n_sequenced": n_sequenced,
        "n_tp53_mutated_samples": n_tp53_mut,
        "tp53_prevalence": overall_prev,
        "tp53_prevalence_floor": floor,
        "by_variant": {v: len(counts.get(v, set()))
                       for v in cfg.hotspot_variants},
    }
    return summary


def write_frequency_table(summaries: List[Dict[str, Any]], out_path: Path,
                          cfg: Config) -> None:
    ensure_dir(out_path.parent)
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "cohort_id", "variant", "residue", "n_patients",
            "n_sequenced", "frequency",
        ])
        # Map variant label like "R175H" -> residue number 175.
        def variant_to_residue(v: str) -> int:
            digits = "".join(ch for ch in v[1:-1] if ch.isdigit())
            return int(digits)

        for s in summaries:
            n = s["n_sequenced"]
            for variant in cfg.hotspot_variants:
                k = s["by_variant"].get(variant, 0)
                freq = (k / n) if n else 0.0
                w.writerow([
                    s["study_id"], variant, variant_to_residue(variant),
                    k, n, f"{freq:.6f}",
                ])


def run(cfg: Config | None = None) -> Path:
    cfg = cfg or Config.load()
    summaries = [fetch_one_cohort(c.id, cfg) for c in cfg.cohorts]
    out = project_root() / "data" / "processed" / "cohort_frequencies.csv"
    write_frequency_table(summaries, out, cfg)
    log.info("Wrote %s", out)
    return out


if __name__ == "__main__":  # pragma: no cover
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run()
