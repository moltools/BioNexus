"""Module for RetroMol-related tasks."""

import logging
from pathlib import Path
from typing import Generator

from sqlalchemy import delete, select
from sqlalchemy.dialects.postgresql import insert
from tqdm import tqdm

from bionexus.config import default_cache_dir
from bionexus.db.engine import SessionLocal
from bionexus.db.models import RetroMolCompound, Ruleset, RetroFingerprint

logger = logging.getLogger(__name__)


def _retro_bits_and_vec(
    result_json: dict[str, dict],
    generator,
) -> Generator[tuple[str, int, list[float], list[float] | None], None, None]:
    """
    Compute RetroMol biosynthetic fingerprint bits, popcount, and ANN vector from result JSON.

    :param result_json: the RetroMol result JSON
    :param generator: the RetroMol FingerprintGenerator instance
    :return: yields tuples of (bit string, popcount, ANN vector binary, ANN vector counted), or None on error
    """
    num_bits = 512

    try:
        from retromol.io import Result as RetroMolResult
        from retromol.fingerprint import FingerprintGenerator
        import numpy as np
        from numpy.typing import NDArray

        result_obj = RetroMolResult.from_serialized(result_json)
        generator: FingerprintGenerator = generator
        # Get fingerprint in shape (N, 512); could return multiple fingerprints if there are multiple optimal mappings
        # fps_counted: NDArray[np.int8] | None = generator.fingerprint_from_result(result_obj, num_bits=num_bits, counted=True, kmer_sizes=[1, 2, 3], kmer_weights={1: 2, 2: 4, 3: 8})
        fps_counted: NDArray[np.int8] | None = generator.fingerprint_from_result(result_obj, num_bits=num_bits, counted=True, strict=False)
        if fps_counted is None:
            # Can happen when coverage is 0.0
            yield None
            return
        for fp_counted in fps_counted:
            fp_binary = (fp_counted > 0).astype(np.int8)
            # Get fp as bit string (e.g. '101010...')
            bitstr = "".join(str(b) for b in fp_binary.flatten().tolist())
            pop = bitstr.count("1")
            # ANN vector binary: 0/1 floats length 512
            vec_binary = [float(b) for b in fp_binary.flatten().tolist()]
            vec_counted = [float(c) for c in fp_counted.flatten().tolist()]
            assert (
                len(bitstr) == num_bits
                and len(vec_binary) == num_bits
                and len(vec_counted) == num_bits
            ), "Fingerprint lengths do not match expected size"
            yield bitstr, pop, vec_binary, vec_counted
    except Exception as e:
        logger.warning(f"Error computing RetroMol fingerprint: {e}")
        yield None


def backfill_retro_fingerprints(
    cache_dir: Path | str, 
    batch: int = 1000,
    recompute: bool = False,
) -> int:
    """
    Compute RetroMol biosynthetic fingerptins for RetroMol results.

    :param cache_dir: directory to use for caching input/output files
    :param batch: numbe rof compounds to process per transaction
    :param recompute: whether to recompute fingerprints even if they exist
    :return: total number of results updated
    """
    done = 0
    last_id = 0

    # Make Path out of cache_dir if it's a string
    if cache_dir is None:
        cache_dir = default_cache_dir()

    if isinstance(cache_dir, str):
        cache_dir = Path(cache_dir)

    # make sure cache_dir exists
    cache_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Import relevant retromol modules
        from retromol.fingerprint import (
            FingerprintGenerator,
            NameSimilarityConfig,
            polyketide_family_of,
            polyketide_ancestors_of,
        )
        from retromol.rules import get_path_default_matching_rules
    except ImportError:
        logger.error("retromol package not found. Please install retromol.")
        return done
    
    COLLAPSE_BY_NAME = {
        "glycosylation": ["glycosyltransferase"],
        "methylation": ["methyltransferase"],
    }
        
    def _setup_fingerprint_generator(yaml_path: str) -> FingerprintGenerator:
        """
        Setup and return a FingerprintGenerator instance.

        :return: FingerprintGenerator instance
        """
        collapse_by_name: list[str] = list(COLLAPSE_BY_NAME.keys())
        cfg = NameSimilarityConfig(
            # family_of=polyketide_family_of,
            # family_repeat_scale=1,
            ancestors_of=polyketide_ancestors_of,
            ancestor_repeat_scale=1,
            symmetric=True,
        )
        generator = FingerprintGenerator(
            matching_rules_yaml=yaml_path,
            collapse_by_name=collapse_by_name,
            name_similarity=cfg
        )
        return generator
    
    # Retrieve all matching rules versions from database and save to yaml files in cache_dir
    with SessionLocal() as s:
        # Get all unique ruleset_id values from RetroMolCompound
        ruleset_ids = s.scalars(select(RetroMolCompound.ruleset_id).distinct()).all()
        
        # Get matching_rules_yaml values for each ruleset_id and save to files
        ruleset_files = {}
        for rid in ruleset_ids:
            rules_q = select(Ruleset).where(Ruleset.id == rid).limit(1)
            ruleset = s.scalars(rules_q).first()
            if ruleset:
                yaml_path = cache_dir / f"matching_rules_{rid}.yaml"
                with open(yaml_path, "w") as f:
                    f.write(ruleset.matching_rules_yaml)
                ruleset_files[rid] = yaml_path

    # Create separate fingerprint generator for each matching ruleset
    generators = {}
    for rid, yaml_path in ruleset_files.items():
        # generators[rid] = _setup_fingerprint_generator(str(yaml_path))
        generators[rid] = _setup_fingerprint_generator(get_path_default_matching_rules())

    # Calculate fingerprints
    total_inserted = 0
    with SessionLocal() as s:
        while True:
            q = select(RetroMolCompound).where(RetroMolCompound.result_json.is_not(None))
            if not recompute:
                q = q.where(~RetroMolCompound.retrofingerprints.any())
            q = q.where(RetroMolCompound.id > last_id).order_by(RetroMolCompound.id).limit(batch)
            
            rows = s.scalars(q).all()
            if not rows:
                break

            # If recomputing, clear existing fingerprints for this batch of compounds
            if recompute:
                s.execute(
                    delete(RetroFingerprint).where(
                        RetroFingerprint.retromol_compound_id.in_([r.id for r in rows])
                    )
                )

            mappings: list[dict] = []

            for c in tqdm(rows, desc="Computering RetroMol fingerprints"):
                # Get appropriate generator for this compound
                gen = generators.get(c.ruleset_id)
                if gen is None:
                    logger.warning(f"No fingerprint generator found for ruleset_id {c.ruleset_id}, skipping compound ID {c.id}")
                    continue

                for fp_result in _retro_bits_and_vec(c.result_json, gen):
                    if fp_result is None:
                        continue
                    bitstr, pop, vec_binary, vec_counted = fp_result

                    mappings.append({
                        "retromol_compound_id": c.id,
                        "biocracker_genbank_id": None,
                        "fp_retro_b512_bit": bitstr,
                        "fp_retro_b512_pop": int(pop),
                        "fp_retro_b512_vec_binary": [float(x) for x in vec_binary],
                        "fp_retro_b512_vec_counted": [float(x) for x in vec_counted],
                    })

                    # Flush in chunks to keep memory + round-trips low
                    if len(mappings) >= batch:
                        s.execute(insert(RetroFingerprint), mappings)
                        total_inserted += len(mappings)
                        mappings.clear()

                done += 1
            
            # Insert any remaning from this page of results
            if mappings:
                s.execute(insert(RetroFingerprint), mappings)
                total_inserted += len(mappings)
                mappings.clear()

            s.commit()
            last_id = rows[-1].id
            logger.info(f"[batch] committed {done} total updated")

    logger.info(f"Fingerprint backfill complete: {done} RetroMol results updated")

    return done
