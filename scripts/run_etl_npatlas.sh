#!/usr/bin/env bash
set -euo pipefail
source .env || true
python -m bionexus.cli upgrade head
python -m bionexus.cli load-npatlas --path "${NPATLAS_FILE:-/data/npatlas.csv}" --format auto --dedupe-on inchikey --chunk-size 20000
python -m bionexus.cli compute-fp --batch 4000
python -m bionexus.cli dump-db --out "dumps/bionexus_$(date +%Y%m%d).dump"
