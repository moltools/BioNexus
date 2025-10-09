# BioNexus

## Quick start

1) Copy `.env.example` -> `.env` and adjust

2) Start Postgres + Adminer:

    ```bash
    docker-compose up -d db adminer
    ```

3) (Optional) Create local env, using `venv`:

    ```bash
    python -m venv .venv
    source .venv/bin/activate
    pip install poetry
    ```

    or using `conda`:

    ```bash
    conda create -n bionexus python=3.10
    conda activate bionexus
    pip install poetry
    ```

4) Install deps:

    ```bash
    poetry install --all-extras
    ```

5) Create/upgrade schema:

    ```bash
    bionexus upgrade head
    ```

6) Load NPAtlas file (CSV/TSV/JSONL accepted; autodetecs headers):

    ```bash
    bionexus load-npatlas
    ```

7) (Optional) Compute fingerprints (requires RDKit):

    ```bash
    bionexus compute-fp
    ```

8) Dump for deployoment:

    ```bash
    bionexus dump-db --out dumps/bionexus_$(date +%Y%m%d).dump
    ```

Adminer: http://localhost:8080 (server: db, user: bionexus, db: bionexus)

## Notes
* The NPAtlas loader is header-flexible. It tries common key names (e.g., `npa_id`, `smiles`, `inchikey`). Update `HEADER_MAP` if export differs.
* Dedupe key defaults to `inchikey`; switch to `ext_id` or `name` if needed.
* `pgvector` is enabled for `compound.fp`. Heavy chemistry remains in Python.