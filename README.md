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

6) Load data:

    * Load NPAtlas compound structures, their names and properties, and organism annotations:

        ```bash
        bionexus load-npatlas
        ```

    * Load MIBiG:

        Make sure to have installed the `bio` extras:

        ```bash
        poetry install --extras bio
        ```

        Then run:

        ```bash
        bionexus load-mibig
        ```

    * (Optional) Compute compound fingerprints:

        Make sure to have installed the `chem` extras:

        ```bash
        poetry install --extras chem
        ```

        Then run:

        ```bash
        bionexus compute-fp-morgan
        ```

    * (Optional) Supplement NPClassifier predictions:

        ```bash
        bionexus annotate-npc
        ```

    * (Optional) Parse compounds with RetroMol:

        Make sure to have installed the `retromol` extras:

        ```bash
        poetry install --extras retromol
        ```

        Then run:

        ```bash
        bionexus parse-compounds
        ```
    

7) Dump for deployoment:

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

<!-- ## Retrobiosynthetic search

Not yet implemented -->

## Database schema & migrations

The BioNexus database schema is defined in `bionexus/db/models.py` using SQLAlchemy ORM classes.  
Schema changes (new tables, columns, constraints, etc.) are tracked with **Alembic migrations** so all environments stay in sync.

### Typical workflow

1. **Edit the models**  
   Add or modify table definitions in `bionexus/db/models.py`.

   ```python
   class MyNewTable(Base):
       __tablename__ = "my_new_table"
       id = mapped_column(BigInteger, primary_key=True)
       name = mapped_column(String(128), nullable=False)
   ```

2. **Generate a new migration**  
   Autogenerate a migration from model diffs:

   ```bash
   bionexus revision -m "add MyNewTable"
   ```

   A file appears under `alembic/versions/` with the DDL operations.

3. **Review and tweak**  
   Open the new file and verify `upgrade()` / `downgrade()`.  
   You can add triggers, extra indexes, or server defaults as needed.

4. **Apply the migration**

   ```bash
   bionexus upgrade
   ```

   (Runs `alembic upgrade head`.)

5. **Rollback if needed**

   ```bash
   bionexus downgrade -1
   ```

---

### Tips

- Ensure **all** models are imported by `bionexus/db/models.py` (or submodules imported there) so Alembic can detect them.  
- Alembic compares `Base.metadata` to the live DB. Confirm your `.env` has:

  ```bash
  BIONEXUS_DB_URL=postgresql+psycopg://user:pass@host:5432/dbname
  ```

- If generated revisions are empty, check your template `alembic/script.py.mako` includes:

  ```mako
  ${upgrades if upgrades else "pass"}
  ${downgrades if downgrades else "pass"}
  ```

- Handy commands:

  ```bash
  bionexus history   # list migrations
  bionexus current   # show current DB revision
  bionexus heads     # show head revision(s)
  ```

