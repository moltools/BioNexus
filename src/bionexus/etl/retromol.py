"""Module for RetroMol-related tasks."""

import logging
import tempfile
import uuid
from pathlib import Path
from typing import Any, Generator, Optional

import joblib
import numpy as np
from sqlalchemy import delete, select
from sqlalchemy.dialects.postgresql import insert
from tqdm import tqdm

from bionexus.config import default_cache_dir
from bionexus.db.engine import SessionLocal
from bionexus.db.models import RetroMolCompound, Ruleset, RetroFingerprint

logger = logging.getLogger(__name__)


_model_cache: dict[str, object | None] = {}


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


    
COLLAPSE_BY_NAME = {
    "glycosylation": ["glycosyltransferase"],
    "methylation": ["methyltransferase"],
}


def _setup_fingerprint_generator(yaml_path: str):
    """
    Setup and return a FingerprintGenerator instance.

    :return: FingerprintGenerator instance
    """
    try:
        # Import relevant retromol modules
        from retromol.fingerprint import (
            FingerprintGenerator,
            NameSimilarityConfig,
            polyketide_family_of,
            polyketide_ancestors_of,
        )
    except ImportError:
        logger.error("retromol package not found. Please install retromol.")
        raise

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
        from retromol.rules import get_path_default_matching_rules
    except ImportError:
        logger.error("retromol package not found. Please install retromol.")
        return done
    
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




def get_paras_model(paras_model_path: Optional[str] = None) -> object | None:
    """
    Load and return the PARAS model from disk, caching it in memory.
    
    :return: the loaded PARAS model, or None if not found
    """
    # Check if model is already cached
    if "paras" in _model_cache:
        return _model_cache["paras"]

    # Check if model path is defined
    if paras_model_path:
        # Model path is defined; attempt to load the model
        path = Path(paras_model_path)
        if path.is_file():
            _model_cache["paras"] = joblib.load(path)
        else:
            _model_cache["paras"] = None
        return _model_cache["paras"]
    else:
        # Model path is not defined
        return None
    

def get_unique_identifier() -> str:
    """
    Generate a unique identifier string.
    
    :return: unique identifier as a string
    """
    return str(uuid.uuid4())


def kmerize_sequence(sequence: list[Any], k: int) -> list[list[Any]]:
    """
    Generate k-mers from a given sequence (forward and backward).
    
    :param sequence: list of elements (e.g., amino acids)
    :param k: length of each k-mer
    :return: list of k-mer strings
    """
    kmers = []
    seq_length = len(sequence)
    
    # Forward k-mers
    for i in range(seq_length - k + 1):
        kmer = sequence[i:i + k]
        kmers.append(kmer)
    
    # Backward k-mers
    for i in range(seq_length - k, -1, -1):
        kmer = sequence[i:i + k]
        kmers.append(kmer)
    
    return kmers


def _compute_gene_cluster(
    cache_dir: Path | str,
    fp_generator,
    gbk_str,
    paras_model_path: Path | str | None = None,
    top_level: str = "cand_cluster",
):
    try:
        from biocracker.antismash import parse_region_gbk_file
        from biocracker.readout import (
            NRPSModuleReadout,
            PKSModuleReadout,
            linear_readouts as biocracker_linear_readouts,
        )
        from biocracker.text_mining import get_default_tokenspecs, mine_virtual_tokens
        from retromol.chem import smiles_to_mol, mol_to_fpr
    except ImportError:
        logger.error("biocracker package not found")
        return

    with tempfile.NamedTemporaryFile(delete=True, suffix=".gbk") as temp_gbk_file:
        temp_gbk_file.write(gbk_str.encode("utf-8"))
        temp_gbk_file.flush()
        gbk_path = temp_gbk_file.name
        tokenspecs = get_default_tokenspecs()
        targets = parse_region_gbk_file(gbk_path, top_level=top_level)  # region or cand_cluster

    level = "gene"
    avg_pred_vals, fps, linear_readouts = [], [], []
    for target in targets:
        pred_vals, raw_kmers = [], []

        # Mine for tokenspecs (i.e., family tokens)
        for mined_tokenspec in mine_virtual_tokens(target, tokenspecs):
            if token_spec := mined_tokenspec.get("token"):
                for token_name, values in COLLAPSE_BY_NAME.items():
                    if token_spec in values:
                        raw_kmers.append([(token_name, None)])
        
        # Optionally load PARAS model
        paras_model = get_paras_model(paras_model_path)

        # Exctract module kmers
        for readout in biocracker_linear_readouts(
            target,
            model=paras_model,
            cache_dir_override=cache_dir,
            level="gene",
            pred_threshold=0.1,
        ):
            kmer, linear_readout = [], []
            for module in readout["readout"]:
                match module:
                    case PKSModuleReadout(module_type="PKS_A") as m:
                        kmer.append(("A", None))
                        pred_vals.append(1.0)
                        linear_readout.append({
                            "id": get_unique_identifier(),
                            "name": "A",
                            "displayName": None,
                            "tags": [],
                            "smiles": None,
                            "morganfingerprint2048r2": None,
                        })
                    case PKSModuleReadout(module_type="PKS_B") as m:
                        kmer.append(("B", None))
                        pred_vals.append(1.0)
                        linear_readout.append({
                            "id": get_unique_identifier(),
                            "name": "B",
                            "displayName": None,
                            "tags": [],
                            "smiles": None,
                            "morganfingerprint2048r2": None,
                        })
                    case PKSModuleReadout(module_type="PKS_C") as m:
                        kmer.append(("C", None))
                        pred_vals.append(1.0)
                        linear_readout.append({
                            "id": get_unique_identifier(),
                            "name": "C",
                            "displayName": None,
                            "tags": [],
                            "smiles": None,
                            "morganfingerprint2048r2": None,
                        })
                    case PKSModuleReadout(module_type="PKS_D") as m:
                        kmer.append(("D", None))
                        pred_vals.append(1.0)
                        linear_readout.append({
                            "id": get_unique_identifier(),
                            "name": "D",
                            "displayName": None,
                            "tags": [],
                            "smiles": None,
                            "morganfingerprint2048r2": None,
                        })
                    case NRPSModuleReadout() as m:
                        substrate_name = m.get("substrate_name", None)
                        substrate_smiles = m.get("substrate_smiles", None)
                        substrate_score = m.get("score", 0.0)
                        kmer.append((substrate_name, substrate_smiles))
                        pred_vals.append(substrate_score)

                        # Calculate fingerprint for SMILES if present
                        if substrate_smiles:
                            clean_mol = smiles_to_mol(substrate_smiles)
                            fp_bits = mol_to_fpr(clean_mol, rad=2, nbs=2048).reshape(-1).astype(np.uint8)
                        else:
                            fp_bits = None

                        linear_readout.append({
                            "id": get_unique_identifier(),
                            "name": substrate_name or "unknown",
                            "displayName": None,
                            "tags": [],
                            "smiles": substrate_smiles,
                            "morganfingerprint2048r2": fp_bits
                        })
                    case _: raise ValueError("Unknown module readout type")

            if len(kmer) > 0:
                raw_kmers.append(kmer)

            if len(linear_readout) >= 2:  # skip too short
                linear_readouts.append({
                    "id": get_unique_identifier(),
                    "name": f"readout_{len(linear_readouts)+1}",
                    "parentSmilesTagged": None,
                    "sequence": linear_readout,
                })

        # Mine for kmers of lengths 1 to 3
        kmers = []
        kmer_lengths = [1, 2, 3]
        for k in kmer_lengths:
            for raw_kmer in raw_kmers:
                kmers.extend(kmerize_sequence(raw_kmer, k))

        # Generate fingerprint
        fp: np.ndarray = fp_generator.fingerprint_from_kmers(kmers, num_bits=512, counted=False)

        # Calculate average prediction value
        avg_pred_val = float(np.mean(pred_vals)) if len(pred_vals) > 0 else 0.0
        
        avg_pred_vals.append(avg_pred_val)
        fps.append(fp)

    return avg_pred_vals, fps, linear_readouts


def backfill_retro_fingerprints_gbk(
    cache_dir: Path | str, 
    batch: int = 1000,
    recompute: bool = False,
    paras_model_path: Optional[Path | str] = None,
) -> int:
    """
    Compute RetroMol biosynthetic fingerptins for GenBank records.

    :param cache_dir: directory to use for caching input/output files
    :param batch: numbe rof compounds to process per transaction
    :param recompute: whether to recompute fingerprints even if they exist
    :param paras_model_path: optional path to the PARAS model file
    :return: total number of results updated
    """
    done = 0
    last_id = 0

    if cache_dir is None:
        cache_dir = default_cache_dir()

    if isinstance(cache_dir, str):
        cache_dir = Path(cache_dir)

    try:
        from retromol.rules import get_path_default_matching_rules
    except ImportError:
        logger.error("retromol and/or biocracker package(s) not found.")

    fp_generator = _setup_fingerprint_generator(get_path_default_matching_rules())
    print(fp_generator)

    total_inserted = 0
    with SessionLocal() as s:
        # while True:
        #     q = select
        pass

    # TODO: update so we save intermediate results in biocracker_genbank (e.g., paras predictions)
    # TODO: read out JSONs biocracker_genbank readouts and parse into fingerprinta for retrofingerprint and link biocracker_genbank ID

    return done
