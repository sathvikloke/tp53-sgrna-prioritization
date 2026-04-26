PY ?= python
SRC := src
PKG := src

.PHONY: all env fetch-cds fetch-cohorts enumerate score-on-target score-off-target rank tables clean test lint help

help:
	@echo "Targets:"
	@echo "  make env               - create conda env from environment.yml"
	@echo "  make fetch-cds         - download NM_000546.6 from NCBI"
	@echo "  make fetch-cohorts     - pull cBioPortal mutation calls"
	@echo "  make enumerate         - enumerate candidate sgRNAs"
	@echo "  make score-on-target   - on-target scoring"
	@echo "  make score-off-target  - off-target scoring (requires GRCH38_DIR)"
	@echo "  make rank              - compose per-cohort rankings"
	@echo "  make tables            - summary JSON + figures"
	@echo "  make all               - run every step in order"
	@echo "  make test              - run pytest"
	@echo "  make clean             - remove generated artifacts"

env:
	conda env create -f environment.yml || conda env update -f environment.yml --prune

fetch-cds:
	$(PY) -m $(PKG).cli fetch-cds

fetch-cohorts:
	$(PY) -m $(PKG).cli fetch-cohorts

enumerate: fetch-cds
	$(PY) -m $(PKG).cli enumerate

score-on-target: enumerate
	$(PY) -m $(PKG).cli score-on-target

score-off-target: score-on-target
	$(PY) -m $(PKG).cli score-off-target

rank: score-off-target fetch-cohorts
	$(PY) -m $(PKG).cli rank

tables: rank
	$(PY) -m $(PKG).cli tables

all:
	$(PY) -m $(PKG).cli all

test:
	$(PY) -m pytest -q tests

lint:
	$(PY) -m compileall -q $(SRC) tests

clean:
	rm -rf data/raw/* data/processed/* results/tables/* results/figures/* results/summary.json
