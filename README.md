# BioNexus

## Setting up database

1) Copy `.env.example` -> `.env` and adjust

2) Start Postgres + Adminer:

    ```bash
    docker-compose up -d db adminer
    ```

    Adminer available at http://localhost:8080 (server: db, user: bionexus, db: bionexus)

3) (Optional) Create local env using `conda`:

    ```bash
    conda create -n bionexus python=3.10
    conda activate bionexus
    pip install -e .[dev]
    ```

4) Create/upgrade schema:

    ```bash
    bionexus upgrade head
    ```

5) Load data:
    
    RetroMol and BioCracker results can be loaded with the load command. For more information, run:

    ```bash 
    bionexus load --help
    ```

    Examples of scripts that extract annotations and references from RetroMol and BioCracker can be found in the `scripts/` folder.

6) Dump for deployoment:

    Make sure the database is running in Docker, then:

    ```bash
    bionexus dump --out ~/Downloads/dumps/bionexus_$(date +%Y%m%d).dump
    ```
