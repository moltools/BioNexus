"""Module for parsing compounds using RetroMol."""

from __future__ import annotations
from typing import List, Dict, Any
import os
import hashlib
import logging
import shutil
import subprocess
from pathlib import Path
from importlib.resources import files

from tqdm import tqdm
from sqlalchemy import select, and_, literal, exists, literal_column
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.dialects.postgresql import insert

from bionexus.db.models import Compound, Ruleset, RetroMolCompound
from bionexus.db.engine import SessionLocal


logger = logging.getLogger(__name__)


ID_COL_NAME = "compound_id"
SMILES_COL_NAME = "smiles"


def _sha256_hex(s: str) -> str:
    """
    Compute the SHA256 hex digest of a given string.

    :param s: input string
    :return: SHA256 hex digest
    """
    return hashlib.sha256((s or "").encode("utf-8")).hexdigest()


def parse_compounds_with_retromol(
    cache_dir: Path, recompute: bool = False, chunk_size: int = 2000, workers: int = 1
) -> int:
    """
    Parse compounds in the database using RetroMol and store results.

    :param cache_dir: directory to use for caching input/output files
    :param recompute: if True, recompute results for all compounds
    :param chunk_size: number of results to upsert in each database transaction
    :param workers: number of parallel workers to use for RetroMol processing
    :return: number of new results inserted
    """
    try:
        # Import retromol package
        import retromol
        import retromol.data
        from retromol.helpers import iter_json

        # Check if retromol is an installed command line tool
        if shutil.which("retromol") is None:
            raise ImportError("retromol command line tool not found in PATH")
    except ImportError:
        logger.error(
            "retromol package or command line tool not found. Please install retromol."
        )
        return -1

    # Limit number of workers to at least one, and 1 below number of CPU cores
    max_workers = max(1, os.cpu_count() - 1)
    workers = min(max(1, workers), max_workers)
    logger.info(
        f"Using {workers} parallel workers for RetroMol processing (max available: {max_workers})"
    )

    # Ensure cache_dir exists
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Get sha256 of default rulesets
    path_matching = files(retromol.data).joinpath("default_matching_rules.yml")
    path_reaction = files(retromol.data).joinpath("default_reaction_rules.yml")
    matching_rules_yaml = path_matching.read_text()
    reaction_rules_yaml = path_reaction.read_text()
    matching_sha = _sha256_hex(matching_rules_yaml)
    reaction_sha = _sha256_hex(reaction_rules_yaml)

    # Determine number of compounds to process
    with SessionLocal() as s:
        # Look up by BOTH section SHAs (each is unique; together uniquely identify row)
        ruleset = s.execute(
            select(Ruleset).where(
                and_(
                    Ruleset.matching_rules_sha256 == matching_sha,
                    Ruleset.reaction_rules_sha256 == reaction_sha,
                )
            )
        ).scalar_one_or_none()

        if ruleset is None:
            ruleset = Ruleset(
                matching_rules_yaml=matching_rules_yaml,
                matching_rules_sha256=matching_sha,
                reaction_rules_yaml=reaction_rules_yaml,
                reaction_rules_sha256=reaction_sha,
                # DO NOT set ruleset_sha256; DB trigger computes it
            )
            s.add(ruleset)
            s.flush()  # obtain PK before commit if needed downstream
            s.commit()
            s.refresh(ruleset)  # get trigger-computed fields (version, ruleset_sha256)

        ruleset_id = ruleset.id

        # Build query of compounds to process
        base_q = select(Compound.id, Compound.smiles).where(
            Compound.smiles.is_not(None)
        )

        if recompute:
            q = base_q
        else:
            # Anti-join: compounds WITHOUT a retromol result for this ruleset
            exists_subq = (
                select(literal(1))
                .select_from(RetroMolCompound)
                .where(
                    and_(
                        RetroMolCompound.compound_id == Compound.id,
                        RetroMolCompound.ruleset_id == ruleset.id,
                    )
                )
                .limit(1)
            )
            q = base_q.where(~exists(exists_subq))

        rows = s.execute(q).all()
        logger.info(f"{len(rows)} compounds to process with RetroMol")

    # Write out rows to TSV in cache_dir, save path
    smiles_path = cache_dir.joinpath("retromol_input.tsv")
    with smiles_path.open("w") as f:
        # Write header
        f.write(f"{ID_COL_NAME}\t{SMILES_COL_NAME}\n")
        for compound_id, smiles in rows:
            f.write(f"{compound_id}\t{smiles}\n")
    logger.info(f"Wrote input TSV for RetroMol to {smiles_path}")

    # Ouputs are written to cache_dir/results.jsonl; delete if exists and recompute
    results_path = cache_dir.joinpath("results.jsonl")
    if recompute and results_path.exists():
        results_path.unlink()
        logger.info(
            f"Deleted existing results file at {results_path} due to recompute=True"
        )

    if not results_path.exists():
        # Run RetroMol on compounds
        subprocess.run(
            [
                "retromol",
                "-o",
                str(cache_dir),
                "batch",
                "--results",
                "jsonl",
                "--jsonl-path",
                str(results_path),
                "--table",
                str(smiles_path),
                "--separator",
                "tab",
                "--id-col",
                ID_COL_NAME,
                "--smiles-col",
                SMILES_COL_NAME,
                "--workers",
                str(workers),
            ],
            check=True,
        )

    # Stream results.jsonl and upset in chunks
    to_upsert = []
    n_seen = 0
    n_bad = 0

    n_inserted_total = 0
    n_updated_total = 0

    def flush_chunk(session, payload: List[Dict[str, Any]]) -> tuple[int, int]:
        """
        Flush a chunk of results to the database using an upsert operation.

        :param session: SQLAlchemy session
        :param payload: list of result dicts to insert/update
        :return: tuple of (n_inserted, n_updated)
        """
        if not payload:
            return (0, 0)
        stmt = insert(RetroMolCompound).values(payload)
        stmt = stmt.on_conflict_do_update(
            constraint="uq_retromol_compound_compound_ruleset",
            set_={
                "result_json": stmt.excluded.result_json,
            },
        )
        stmt = stmt.returning(literal_column("xmax = 0").label("inserted"))
        res = session.execute(stmt)
        session.commit()
        rows = res.fetchall()
        ins = sum(1 for r in rows if r.inserted)
        upd = len(rows) - ins
        return (ins, upd)

    with SessionLocal() as s:
        for rec in tqdm(
            iter_json(results_path, jsonl=True), desc="Loading RetroMol results"
        ):
            n_seen += 1

            # Basic validation
            result_json = rec.get("result")
            if isinstance(result_json, dict):
                payload_json = result_json
            else:
                n_bad += 1
                continue

            # Extract compound_id
            cid = result_json.get("input_id")
            if cid is None:
                logger.warning("Missing compound_id in result; skipping")
                n_bad += 1
                continue

            row = {
                "compound_id": int(cid),
                "ruleset_id": int(ruleset_id),
                "result_json": payload_json,  # plain dict, JSON-safe
            }
            to_upsert.append(row)

            if len(to_upsert) >= chunk_size:
                try:
                    ins, upd = flush_chunk(s, to_upsert)
                    n_inserted_total += ins
                    n_updated_total += upd
                except SQLAlchemyError as e:
                    s.rollback()
                    logger.error(f"Database error during upsert: {e}")
                    raise
                finally:
                    to_upsert.clear()

        # Flush any remaining
        if to_upsert:
            try:
                ins, upd = flush_chunk(s, to_upsert)
                n_inserted_total += ins
                n_updated_total += upd
            except SQLAlchemyError as e:
                s.rollback()
                logger.error(f"Database error during final upsert: {e}")
                raise
            finally:
                to_upsert.clear()

    logger.info(
        f"Inserted {n_inserted_total} new RetroMol results, updated {n_updated_total} existing results"
    )

    return n_inserted_total
