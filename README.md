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

6) Load NPAtlas:

    ```bash
    bionexus load-npatlas
    ```

7) (Optional) Compute fingerprints (requires RDKit):

    Make sure to have installed the `chem` extras:

    ```bash
    poetry install --extras chem
    ```

    Then run:

    ```bash
    bionexus compute-fp
    ```

8) Dump for deployoment:

    ```bash
    bionexus dump-db --out dumps/bionexus_$(date +%Y%m%d).dump
    ```

Adminer: http://localhost:8080 (server: db, user: bionexus, db: bionexus)

## Compound search

Make sure to have installed the `chem` extras:

```bash
poetry install --extras chem
```

Exact search:

```bash
bionexus search-jaccard --smiles "CCO" --top-k 5
```

Exact search returns `top-k` results with lowest Jaccard distance (highest similarity). Hybrid search speeds up the process by using the `pgvector` extensions and the `<#>` operator (inner-product distance) to quickly find the 2000 nearest neighbors and then re-ranks these using the exact Jaccard distance.

Hybrid search:
```bash
bionexus search-jaccard --smiles "CCO" --top-k 5 --hybrid   
```
