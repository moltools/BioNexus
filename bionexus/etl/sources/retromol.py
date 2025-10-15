from __future__ import annotations
from typing import List, Tuple, Dict, Any, Optional
import os, hashlib, logging, shutil, subprocess, json
from pathlib import Path
from requests import Session
from tqdm import tqdm
from sqlalchemy import select, delete, and_, literal, exists
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.dialects.postgresql import insert
from bionexus.db.models import Compound, Ruleset, RetroMolCompound
from bionexus.db.engine import SessionLocal
from importlib.resources import files

logger = logging.getLogger(__name__)

ID_COL_NAME = "compound_id"
SMILES_COL_NAME = "smiles"

def _sha256_hex(s: str) -> str:
    return hashlib.sha256((s or "").encode("utf-8")).hexdigest()

def parse_compounds_with_retromol(
    cache_dir: Path,
    recompute: bool = False,
    chunk_size: int = 2000,
    workers: int = 1
) -> int:
    try:
        # import retromol package
        import retromol
        import retromol.data
        from retromol import streaming

        # check if retromol is an installed command line tool
        if shutil.which("retromol") is None:
            raise ImportError("retromol command line tool not found in PATH")
    except ImportError as e:
        logger.error("retromol package or command line tool not found. Please install retromol.")
        return -1
    
    # limit number of workers to at least one, and 1 below number of CPU cores
    max_workers = max(1, os.cpu_count() - 1)
    if workers < 1:
        workers = 1
    elif workers > max_workers:
        workers = max_workers
    logger.info(f"Using {workers} parallel workers for RetroMol processing (max available: {max_workers})")

    # ensure cache_dir exists
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # get sha256 of default rulesets
    path_matching = files(retromol.data).joinpath("default_matching_rules.yml")
    path_reaction = files(retromol.data).joinpath("default_reaction_rules.yml")
    matching_rules_yaml = path_matching.read_text()
    reaction_rules_yaml = path_reaction.read_text()
    matching_sha = _sha256_hex(matching_rules_yaml)
    reaction_sha = _sha256_hex(reaction_rules_yaml)
    ruleset_sha  = _sha256_hex((matching_rules_yaml or "") + "||" + (reaction_rules_yaml or ""))

    # determine number of compounds to process
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
            s.flush()   # obtain PK before commit if needed downstream
            s.commit()
            s.refresh(ruleset)  # get trigger-computed fields (version, ruleset_sha256)

        # Build query of compounds to process
        base_q = select(Compound.id, Compound.smiles).where(Compound.smiles.is_not(None))

        if recompute:
            q = base_q
        else:
            # anti-join: compounds WITHOUT a retromol result for this ruleset
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

    # write out rows to TSV in cache_dir, save path
    smiles_path = cache_dir.joinpath("retromol_input.tsv")
    with smiles_path.open("w") as f:
        # write header
        f.write(f"{ID_COL_NAME}\t{SMILES_COL_NAME}\n")
        for compound_id, smiles in rows:
            f.write(f"{compound_id}\t{smiles}\n")
    logger.info(f"Wrote input TSV for RetroMol to {smiles_path}")

    # TODO: get hashes of most recent ruleset
    
    # TODO: retrieve all compounds with unique id that do not have a retromol with most recent ruleset yet, or all if recompute

    # TODO: use multiprocessing command line interface of retromol to compute retromol results for all compounds, write results to cache_dir

    return 0

def parse_bgcs_with_retromol(recompute: bool = False, chunk_size: int = 2000) -> int:
    return 0
