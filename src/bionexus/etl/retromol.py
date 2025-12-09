"""Module for RetroMol-related tasks."""

import logging
import tempfile
import uuid
from pathlib import Path
from typing import Any, Generator, Optional
import re

import joblib
import numpy as np
from sqlalchemy import delete, select, func, insert
from sqlalchemy.dialects.postgresql import insert
from tqdm import tqdm

from bionexus.config import default_cache_dir
from bionexus.db.engine import SessionLocal
from bionexus.db.models import (
    BioCrackerGenBank,
    RetroMolCompound,
    Ruleset,
    RetroFingerprint,
)

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


TOKENSPEC_EPOXIDATION: dict = {
    "any": {
        # Direct epoxidation-related names
        "epoxidase",
        "epoxygenase",
        "epoxidation",
        "epoxide-forming monooxygenase",
        "epoxide hydrolase",  # often annotated in same neighbourhood
    },
    "rx": [
        # epoxidase / epoxygenase variants
        re.compile(r"\b(epoxidase|epoxygenase)\b", re.I),
        re.compile(r"\bepoxidation\b", re.I),
        # epoxy* with oxygenase/monooxygenase
        re.compile(r"\bepoxy\w*(monooxygenase|oxygenase)\b", re.I),
        # Some classic style annotations ("styrene monooxygenase", etc.) with epoxide context
        re.compile(r"\bstyrene\s+monooxygenase\b", re.I),
    ],
    "bonus_if": {
        # Typical cofactors / tailoring context
        "FAD",
        "FMN",
        "flavin",
        "monooxygenase",
        "P450",
        "cytochrome P450",
        "oxidoreductase",
        "NRPS",
        "PKS",
        "T1PKS",
        "trans-AT",
    },
    "weight": 1.0,
    "bonus_weight": 0.3,
    # Require at least a reasonably strong signal; most true epoxidases
    # will have 'epoxidase'/'epoxidation' or similar in the annotation.
    "min_score": 2.5,
}


TOKENSPEC_AMINATION: dict = {
    "any": {
        "aminotransferase",
        "amine transferase",
        "transaminase",
        "amination",
        "aminating enzyme",
        "aminase",
    },
    "rx": [
        re.compile(r"\baminotransferase\b", re.I),
        re.compile(r"\btransaminase\b", re.I),
        re.compile(r"\bomega[-\s]?transaminase\b", re.I),
        # PLP-dependent is a strong hint when combined with above
        re.compile(r"\bPLP[-\s]?dependent\b", re.I),
    ],
    "bonus_if": {
        "PLP",
        "pyridoxal phosphate",
        "pyridoxamine",
        "beta-amino",
        "amino transfer",
        "amino group transfer",
        "NRPS",
        "PKS",
        "T1PKS",
    },
    "weight": 1.0,
    "bonus_weight": 0.3,
    # One strong direct hit (e.g. 'aminotransferase') + maybe a bonus
    # should be enough to fire, so keep threshold modest.
    "min_score": 2.0,
}


# --- Halogenation split by halogen ---

TOKENSPEC_FLUORINATION: dict = {
    "any": {
        "fluorinase",
        "fluorination",
        "fluorinating enzyme",
        "fluorinated",
    },
    "rx": [
        re.compile(r"\bfluorinase\b", re.I),
        re.compile(r"\bfluorinat(e|ing|ion|ed)\b", re.I),
        # specific halogenases with fluoro-, but exclude dehalogenases
        re.compile(r"(?<!de)\bfluoro\w*halogenase\b", re.I),
    ],
    "bonus_if": {
        "FAD",
        "FMN",
        "flavin",
        "monooxygenase",
        "oxidoreductase",
        "fluorination",
        "PKS",
        "NRPS",
        "trans-AT",
        "tAT",
    },
    "weight": 1.0,
    "bonus_weight": 0.4,
    "min_score": 3.0,
}


TOKENSPEC_CHLORINATION: dict = {
    "any": {
        "chlorinase",
        "chlorination",
        "chlorinated",
    },
    "rx": [
        re.compile(r"\bchlorinase\b", re.I),
        re.compile(r"\bchlorinat(e|ing|ion|ed)\b", re.I),
        # halogenase explicitly tied to chloro- context, avoid dehalogenase
        re.compile(r"(?<!de)\bchloro\w*halogenase\b", re.I),
        # classic tryptophan chlorinase style annotations
        re.compile(r"\btryptophan\s*\d?-?(chlorinase|chlorinase-like)\b", re.I),
    ],
    "bonus_if": {
        "FAD",
        "FMN",
        "flavin",
        "monooxygenase",
        "oxidoreductase",
        "chlorination",
        "PKS",
        "NRPS",
        "trans-AT",
        "tAT",
    },
    "weight": 1.0,
    "bonus_weight": 0.4,
    "min_score": 3.0,
}


TOKENSPEC_BROMINATION: dict = {
    "any": {
        "brominase",
        "bromination",
        "brominated",
    },
    "rx": [
        re.compile(r"\bbrominase\b", re.I),
        re.compile(r"\bbrominat(e|ing|ion|ed)\b", re.I),
        re.compile(r"(?<!de)\bbromo\w*halogenase\b", re.I),
        re.compile(r"\btryptophan\s*\d?-?brominase\b", re.I),
    ],
    "bonus_if": {
        "FAD",
        "FMN",
        "flavin",
        "monooxygenase",
        "oxidoreductase",
        "bromination",
        "PKS",
        "NRPS",
        "trans-AT",
        "tAT",
    },
    "weight": 1.0,
    "bonus_weight": 0.4,
    "min_score": 3.0,
}


TOKENSPEC_IODINATION: dict = {
    "any": {
        "iodinase",
        "iodination",
        "iodinated",
    },
    "rx": [
        re.compile(r"\biodinase\b", re.I),
        re.compile(r"\biodinat(e|ing|ion|ed)\b", re.I),
        re.compile(r"(?<!de)\biodo\w*halogenase\b", re.I),
    ],
    "bonus_if": {
        "FAD",
        "FMN",
        "flavin",
        "monooxygenase",
        "oxidoreductase",
        "iodination",
        "PKS",
        "NRPS",
        "trans-AT",
        "tAT",
    },
    "weight": 1.0,
    "bonus_weight": 0.4,
    "min_score": 3.0,
}


def _compute_gene_cluster(
    cache_dir: Path | str,
    fp_generator,
    targets,
    biocracker_genbank_id: int,
):
    try:
        from biocracker.readout import (
            NRPSModuleReadout,
            PKSModuleReadout,
            linear_readouts as biocracker_linear_readouts,
        )
        from biocracker.text_mining import get_default_tokenspecs, mine_virtual_tokens
    except ImportError:
        logger.error("biocracker package not found")
        return
    

    def get_tokenspecs():
        default = get_default_tokenspecs()
        default["epoxidation"] = TOKENSPEC_EPOXIDATION
        default["amination"] = TOKENSPEC_AMINATION
        default["fluorination"] = TOKENSPEC_FLUORINATION
        default["chlorination"] = TOKENSPEC_CHLORINATION
        default["bromination"] = TOKENSPEC_BROMINATION
        default["iodination"] = TOKENSPEC_IODINATION
        return default

    tokenspecs = get_tokenspecs()
    mappings: list[dict[str, Any]] = []

    for target in targets:
        raw_kmers: list[list[Any]] = []

        if not target:
            continue

        # Mine for tokenspecs (i.e., family tokens)
        if not (isinstance(target, list) and target and isinstance(target[0], dict)):
            for mined_tokenspec in mine_virtual_tokens(target, tokenspecs):
                if token_spec := mined_tokenspec.get("token"):
                    for token_name, values in COLLAPSE_BY_NAME.items():
                        if token_spec in values:
                            raw_kmers.append([(token_name, None)])

            # Extract module kmers directly from BioCracker readouts
            for readout in biocracker_linear_readouts(
                target,
                model=None,
                cache_dir_override=cache_dir,
                level="gene",
                pred_threshold=0.1,
            ):
                kmer = []
                for module in readout["readout"]:
                    match module:
                        case PKSModuleReadout(module_type="PKS_A") as m: kmer.append(("A", None))
                        case PKSModuleReadout(module_type="PKS_B") as m: kmer.append(("B", None))
                        case PKSModuleReadout(module_type="PKS_C") as m: kmer.append(("C", None))
                        case PKSModuleReadout(module_type="PKS_D") as m: kmer.append(("D", None))
                        case NRPSModuleReadout() as m:
                            substrate_name = m.get("substrate_name", None)
                            substrate_smiles = m.get("substrate_smiles", None)
                            kmer.append((substrate_name, substrate_smiles))
                        case _: raise ValueError("Unknown module readout type")

                if len(kmer) > 0:
                    raw_kmers.append(kmer)
        else:
            # Use cached BioCracker readouts stored on the BioCrackerGenBank record
            for readout in target:
                kmer = []
                for module in readout.get("readout", []):
                    module_type = module.get("module_type")
                    match module_type:
                        case "PKS_A": kmer.append(("A", None))
                        case "PKS_B": kmer.append(("B", None))
                        case "PKS_C": kmer.append(("C", None))
                        case "PKS_D": kmer.append(("D", None))
                        # case PKSModuleReadout(kind="PKS_module", module_type="UNCLASSIFIED") as m: kmer.append(("A", None))
                        case "UNCLASSIFIED": kmer.append(("A", None))
                        case _:
                            substrate_name = module.get("substrate_name", None)
                            substrate_smiles = module.get("substrate_smiles", None)
                            # if substrate_name is None and substrate_smiles is None:
                            #     raise ValueError("Unknown module readout type")
                            if substrate_smiles == "O=NN(O)CCC[C@H](N)(C(=O)O":
                                substrate_smiles = "OC(C(CCCN(N=O)O)N)=O"
                            kmer.append((substrate_name, substrate_smiles))

                if len(kmer) > 0:
                    raw_kmers.append(kmer)

        if not raw_kmers:
            continue

        # Mine for kmers of lengths 1 to 3
        kmers = []
        kmer_lengths = [1, 2, 3]
        for k in kmer_lengths:
            for raw_kmer in raw_kmers:
                kmers.extend(kmerize_sequence(raw_kmer, k))

        if not kmers:
            continue

        # Generate fingerprint
        fps_counted: np.ndarray | None = fp_generator.fingerprint_from_kmers(
            kmers,
            num_bits=512,
            counted=True,
        )
        if fps_counted is None:
            continue
        fps_counted = np.atleast_2d(fps_counted)

        for fp_counted in fps_counted:
            fp_binary = (fp_counted > 0).astype(np.int8)
            bitstr = "".join(str(int(b)) for b in fp_binary.flatten().tolist())
            pop = int(fp_binary.sum())

            mappings.append({
                "retromol_compound_id": None,
                "biocracker_genbank_id": biocracker_genbank_id,
                "fp_retro_b512_bit": bitstr,
                "fp_retro_b512_pop": pop,
                "fp_retro_b512_vec_binary": [float(x) for x in fp_binary.flatten().tolist()],
                "fp_retro_b512_vec_counted": [float(x) for x in fp_counted.flatten().tolist()],
            })

    return mappings


def backfill_retro_fingerprints_gbk(
    cache_dir: Path | str, 
    batch: int = 1000,
    recompute: bool = False,
) -> int:
    """
    Compute RetroMol biosynthetic fingerptins for GenBank records.

    :param cache_dir: directory to use for caching input/output files
    :param batch: numbe rof compounds to process per transaction
    :param recompute: whether to recompute fingerprints even if they exist
    :return: total number of results updated
    """
    done = 0
    total_inserted = 0
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

    with SessionLocal() as s:

        # compute total for tqdm
        count_q = select(func.count()).select_from(BioCrackerGenBank)
        if not recompute:
            count_q = count_q.where(~BioCrackerGenBank.retrofingerprints.any())
        total_to_process = s.scalars(count_q).first() or 0
        pbar = tqdm(total=total_to_process, desc="Computing RetroMol fingerprints for BioCrackerGenBank")

        while True:
            q = (
                select(BioCrackerGenBank)
                .where(BioCrackerGenBank.id > last_id)
                .order_by(BioCrackerGenBank.id)
                .limit(batch)
            )
            if not recompute:
                q = q.where(~BioCrackerGenBank.retrofingerprints.any())

            rows = s.scalars(q).all()
            if not rows:
                break

            deleted_existing = False
            batch_ids = [r.id for r in rows]
            if recompute and batch_ids:
                s.execute(
                    delete(RetroFingerprint).where(
                        RetroFingerprint.biocracker_genbank_id.in_(batch_ids)
                    )
                )
                deleted_existing = True

            batch_inserted = 0

            for gbk in rows:
                targets = gbk.result_json.get("targets") or []
                if targets:
                    mappings = _compute_gene_cluster(cache_dir, fp_generator, targets, gbk.id)
                    if mappings:
                        s.execute(insert(RetroFingerprint), mappings)
                        inserted = len(mappings)
                        batch_inserted += inserted
                        total_inserted += inserted
                        done += 1

                pbar.update(1)

            if deleted_existing or batch_inserted:
                s.commit()
            last_id = rows[-1].id

        pbar.close()

    logger.info(f"Fingerprint backfill complete: {done} BioCrackerGenBank entries updated, {total_inserted} fingerprints inserted")

    return done
