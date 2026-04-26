#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Download GRCh38 primary-assembly chromosomes and lay them out for Cas-OFFinder.
#
# Usage:
#     scripts/setup_grch38.sh /path/to/grch38_dir
#
# After the script finishes, set:
#     export GRCH38_DIR=/path/to/grch38_dir
#
# The script downloads chr1..22, chrX, chrY, chrM (mitochondrial) from UCSC
# in 2bit format and converts them to FASTA. Total disk usage ~3.1 GB.
#
# Requirements:
#     - curl
#     - twoBitToFa (from UCSC; install via `conda install -c bioconda ucsc-twobittofa`)
# -----------------------------------------------------------------------------
set -euo pipefail

DEST=${1:-grch38}
mkdir -p "$DEST"
cd "$DEST"

if ! command -v twoBitToFa >/dev/null 2>&1; then
    echo "ERROR: twoBitToFa is required. Install with:"
    echo "    conda install -c bioconda ucsc-twobittofa"
    exit 1
fi

# Single 2bit file for the whole genome:
URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit"
if [[ ! -f hg38.2bit ]]; then
    echo "Downloading $URL ..."
    curl -L -o hg38.2bit "$URL"
fi

# Split into per-chromosome FASTAs (Cas-OFFinder accepts a directory of FASTAs).
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
        chr21 chr22 chrX chrY chrM)

for c in "${CHROMS[@]}"; do
    if [[ ! -f "${c}.fa" ]]; then
        echo "Extracting ${c} ..."
        twoBitToFa -seq="$c" hg38.2bit "${c}.fa"
    fi
done

echo
echo "GRCh38 ready in: $(pwd)"
echo "Now run:    export GRCH38_DIR=$(pwd)"
