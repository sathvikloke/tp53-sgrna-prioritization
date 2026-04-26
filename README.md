# TP53 sgRNA prioritization

Cancer-type-stratified prioritization of CRISPR-Cas9 sgRNAs targeting TP53
missense hotspots, using live data from NCBI (the canonical TP53 CDS) and
cBioPortal (per-cohort variant frequencies in TCGA Pan-Cancer Atlas
ovarian, pancreatic, and colorectal cohorts).

The pipeline outputs three ranked sgRNA tables (one per cancer type), a
PAM-desert supplement listing hotspot codons that lack a usable SpCas9
NGG-PAM guide, and a reproducibility bundle (run metadata, MD5s,
versioned config).

## Therapeutic framing and scope

The framing is locked at the configuration level: this pipeline assumes
**HDR-mediated correction** of a TP53 hotspot codon back to the wild-type
amino acid via a homology-directed-repair donor template. The pipeline
**only** selects, scores, and ranks the sgRNAs that make HDR feasible
near each hotspot. The following are explicitly **out of scope** in this
version:

* Donor template (HDR repair template) design.
* Knockout / NHEJ-based strategies (different cut-site constraints).
* Base-editing or prime-editing strategies (different scoring needs).
* In silico delivery vehicle selection (LNP / AAV / RNP).
* Any wet-lab validation; the outputs are computational candidates, not
  therapies.

The framing is enforced by `config/params.yaml::framing.mode = "hdr"`;
loading any other value raises `ValueError`. If you need a different
framing, that is a future-work item, not a config switch.

## Related work

This pipeline is small and deliberately transparent. It does not compete
with the major sgRNA design tools — it leans on them.

* **CRISPOR** (Concordet & Haeussler, 2018) is the reference SpCas9
  on/off-target design tool. Our off-target aggregation uses the
  CRISPOR-style `100/(100+sum CFD)` specificity score and our CFD
  weights are loaded from the canonical Doench-2016 distribution that
  CRISPOR ships (`scripts/download_cfd_weights.py`).
* **CHOPCHOP** (Labun et al., 2019) and **Benchling** offer per-region
  guide selection with similar scoring families.
* **Doench 2016 / Azimuth** (Rule Set 2) is the published gold-standard
  on-target scorer. We support it as an opt-in backend
  (`ONTARGET_BACKEND=azimuth`); the default is a transparent
  feature-based scorer used as a fallback when Azimuth is not
  installed.
* **CRISPRme** (Cancellieri et al., 2023) extends off-target search to
  variant-aware genomes; we currently search GRCh38 reference only.

The contribution of this pipeline is not a new scorer. It is the
**cancer-type stratification**: re-ranking the same hotspot-targeting
sgRNAs by per-cohort variant prevalence (HGSOC vs PDAC vs CRC) and
surfacing PAM-desert hotspots — i.e. residues with no usable SpCas9
guide — as a first-class output instead of silently dropping them.

## Quick start

```bash
# 1. Create the conda environment (installs Python deps + cas-offinder).
make env
conda activate tp53-sgrna

# 2. Download GRCh38 once and point at it.
scripts/setup_grch38.sh /path/to/grch38_dir
export GRCH38_DIR=/path/to/grch38_dir

# 3. Run the full pipeline.
make all

# 4. Inspect outputs.
ls results/tables/
cat results/summary.json
```

Per-step targets: `make fetch-cds`, `make fetch-cohorts`, `make enumerate`,
`make score-on-target`, `make score-off-target`, `make rank`, `make tables`.

Tests: `make test`.

## What the pipeline does

1. **Fetch the TP53 CDS** (NM_000546.6, 1,182 nt) from NCBI E-utilities and
   build a residue -> codon index. Pinned MD5 in `data/raw/.meta.json`.
2. **Enumerate candidate sgRNAs** within +/- 10 nt of every configured
   hotspot codon, on both strands, with explicit cut-site bookkeeping.
3. **Score on-target activity.** Two backends:
   * **Default (feature-based).** A transparent, lightweight scorer over
     GC content, PAM-proximal nucleotide preference, seed-region
     position weights, and homopolymer penalty. The features are
     motivated by Doench 2014/2016 but the scorer is **not** Doench
     Rule Set 2 — it is a much smaller approximation and should be
     treated as a sanity check, not a publication-grade activity score.
   * **Azimuth backend (opt-in).** The Doench 2016 Rule Set 2 (Azimuth)
     gradient-boosted model. Activate with `ONTARGET_BACKEND=azimuth`
     after `pip install azimuth`. This is the recommended setting for
     any results that will leave the lab.
4. **Score off-target activity** by running Cas-OFFinder against GRCh38
   (up to 4 mismatches) and computing a CRISPOR-style aggregate
   specificity score from a CFD matrix. By default the **published
   Doench-2016 CFD weights** (mismatch + PAM matrices) are loaded from
   `data/raw/cfd/*.csv` after running
   `scripts/download_cfd_weights.py`. If those CSVs are absent the
   pipeline falls back to a transparent position-only approximation
   and emits a warning. The active mode is reported by
   `cfd.cfd_mode()`. The on-target hit is excluded by **genomic
   coordinate** (the configured TP53 locus on chr17), not by a
   zero-mismatch heuristic — this correctly counts the TP53P1
   pseudogene on chr1 as off-target rather than hiding it.
5. **Fetch per-cohort hotspot frequencies** from cBioPortal for
   ov_tcga_pan_can_atlas_2018, paad_tcga_pan_can_atlas_2018,
   coadread_tcga_pan_can_atlas_2018.
6. **Compose a per-cohort composite score**, in two flavors:

   * `composite_score` — within-residue min-max normalization of on/off
     target scores, with frequency normalized across residues. This
     score is meaningful **only** for ranking candidates targeting the
     **same** residue (which is what tier_1/2/3 labels are for); it is
     **not** meaningful across residues, because the within-residue
     normalization erases absolute differences in spacer quality.
   * `composite_score_absolute` — uses the raw on/off scores directly
     (already in [0,1]) plus the population-normalized frequency. This
     **is** comparable across residues and is the recommended sort key
     for any cross-residue master table. The default weights are 0.5
     on-target / 0.3 off-target / 0.2 frequency; use
     `scripts/sensitivity.py` to inspect rank stability under other
     weight settings.

7. **Surface PAM deserts** as a first-class negative finding in
   `results/tables/pam_deserts.csv`, with rescue options listed for
   SpCas9-NG (NGN), SpCas9-VQR (NGA), and Cas12a (TTTV).

## Project layout

```
config/params.yaml                     # all locked parameters
src/                                   # pipeline modules
  utils.py                             # sequence helpers, IO, provenance
  config.py                            # YAML loader + dataclass
  fetch_cds.py                         # NM_000546.6 fetch + validate
  fetch_cohort_mutations.py            # cBioPortal fetch + frequencies
  enumerate_sgrnas.py                  # PAM/cut-site enumeration
  score_on_target.py                   # on-target scoring (feature/Azimuth)
  cfd.py                               # CFD-style off-target scoring core
  score_off_target.py                  # Cas-OFFinder driver + aggregation
  compose_rankings.py                  # composite scoring + tiers
  make_tables.py                       # summary JSON + figures
  cli.py                               # unified `python -m src.cli ...`
tests/                                 # pytest unit + smoke tests
scripts/setup_grch38.sh                # one-time GRCh38 download/prep
data/                                  # raw + processed (ignored by git)
results/                               # tables + figures + summary.json
```

## Configuration

All parameters live in `config/params.yaml`. The file is loaded once at
the top of every script; do not edit during a run. Composite weights
must sum to 1.0 (validated at load time).

To reproduce a published run, copy the matching `params.yaml` from the
release tag.

## Caveats and limitations

* **Computational only.** No wet-lab validation; the pipeline produces
  prioritization candidates, not therapies.
* **CDS-based PAM search.** Hotspot residues 175, 220, 245, 248, 249,
  273, 282 sit several codons inside their respective TP53 exons, so a
  +/- 10 nt edit window stays within a single exon. For other residues,
  the user must re-validate against genomic context.
* **PAM deserts are real.** Some hotspot codons have no NGG-PAM'd
  20-mer with a cut site in the +/- 10 nt window. They are reported as
  `pam_deserts.csv` with explicit Cas-variant rescue suggestions.
* **Off-target search reflects the reference genome.** Patient-specific
  variants in off-target sites are not modeled.
* **Cohort sizes are modest** (PAAD n ~107). Frequencies are
  TCGA-specific and may differ in other populations.
* **On-target backend.** The default feature-based scorer is a
  lightweight approximation, not the Doench 2016 Rule Set 2 model. It
  is useful for sanity checking and for reproducible rankings without
  external dependencies, but it should not be reported as Rule Set 2.
  For any publication-grade rankings, install Azimuth and set
  `ONTARGET_BACKEND=azimuth`.
* **CFD weights.** Run `scripts/download_cfd_weights.py` to fetch the
  published Doench-2016 CFD matrices (CRISPOR-distributed pickles,
  pinned to a specific commit and verified by SHA256 before any
  pickle is loaded) and convert them to CSV. Without this step, the
  pipeline runs with a clearly-labeled position-only approximation
  and logs a warning at startup. Off-target numbers from approx mode
  should not be used for publication.
* **Cross-residue comparison.** Use `composite_score_absolute`, not
  the within-residue normalized `composite_score`, when comparing
  spacers across different residues.
* **Weight choices.** The 0.5 / 0.3 / 0.2 default is a defensible
  prior, not the Bayes-optimal choice. `scripts/sensitivity.py`
  re-ranks across a small simplex of weight settings so reviewers can
  see which spacers are robust and which only top the table at one
  particular weighting.

## Reproducibility

Every result is traceable. `results/summary.json` records:

* the run version, UTC timestamp, Python and platform versions;
* the URL, HTTP status, fetch time, and MD5 of the TP53 CDS;
* the cBioPortal endpoint hit per cohort with a fetched-at timestamp;
* package versions (`requests`, `pyyaml`, `pandas`, `numpy`, `matplotlib`);
* counts (number of candidates, number of frequency rows).

For a deposited run, tag a release in git and push to Zenodo via the
GitHub-Zenodo integration.

## License

MIT (see `LICENSE`).

## Citation

If you use this pipeline, please cite the repository (see `pyproject.toml`
for version metadata) and the underlying data sources:

* NCBI RefSeq (NM_000546.6).
* TCGA Pan-Cancer Atlas / cBioPortal.
* GRCh38 / UCSC.
* Doench, J. G. et al. (2014, 2016). On-target / off-target sgRNA models.
* Bae, S., Park, J., Kim, J.-S. (2014). Cas-OFFinder. *Bioinformatics*.
