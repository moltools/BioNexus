"""Module for chemistry-related tasks, such as fingerprint computation and InChIKey conversion."""

from __future__ import annotations

import logging

from sqlalchemy import or_, select
from tqdm import tqdm

from bionexus.db.engine import SessionLocal
from bionexus.db.models import Compound

logger = logging.getLogger(__name__)


def get_atom_counts(smiles: str) -> dict[str, int] | None:
    """
    Get atom counts for common elements from a SMILES string using RDKit.

    :param smiles: the SMILES string to analyze
    :return: a dictionary with atom counts, or None if parsing fails
    """
    try:
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")  # suppress warnings
        m = Chem.MolFromSmiles(smiles)
        if not m:
            return None
        atom_counts = {el: 0 for el in ["C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I"]}
        for atom in m.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in atom_counts:
                atom_counts[symbol] += 1
        return {
            "c_count": atom_counts["C"],
            "h_count": atom_counts["H"],
            "n_count": atom_counts["N"],
            "o_count": atom_counts["O"],
            "s_count": atom_counts["S"],
            "p_count": atom_counts["P"],
            "f_count": atom_counts["F"],
            "cl_count": atom_counts["Cl"],
            "br_count": atom_counts["Br"],
            "i_count": atom_counts["I"],
        }
    except Exception:
        return None


def smiles_to_inchikey(smiles: str) -> str | None:
    """
    Convert a SMILES string to an InChIKey using RDKit.

    :param smiles: the SMILES string to convert
    :return: the InChIKey string, or None if conversion fails
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem

        RDLogger.DisableLog("rdApp.*")  # suppress warnings
        m = Chem.MolFromSmiles(smiles)
        if not m:
            return None
        inchi = Chem.MolToInchi(m)
        inchikey = AllChem.InchiToInchiKey(inchi)
        return inchikey
    except Exception:
        return None


def _morgan_bits_and_vec(
    smiles: str, radius: int = 2, nbits: int = 2048
) -> tuple[str | None, int | None, list[float] | None]:
    """
    Compute Morgan fingerprint bitstring, population count, and vector from SMILES.

    :param smiles: the SMILES string
    :param radius: the radius for the Morgan fingerprint
    :param nbits: the number of bits for the fingerprint
    :return: a tuple of (bitstring, population count, vector), or (None, None, None) if computation fails
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

        RDLogger.DisableLog("rdApp.*")  # suppress warnings
        m = Chem.MolFromSmiles(smiles)
        if not m:
            return None, None, None
        gen = GetMorganGenerator(radius=radius, fpSize=nbits, includeChirality=False)
        bv = gen.GetFingerprint(m)
        bitstr = bv.ToBitString()  # "010101..."
        pop = bitstr.count("1")
        # ANN vector: 0/1 floats length 2048
        vec = [1.0 if ch == "1" else 0.0 for ch in bitstr]
        return bitstr, pop, vec
    except Exception:
        return None, None, None


def backfill_fingerprints(batch: int = 1000, recompute: bool = False, radius: int = 2, nbits: int = 2048) -> int:
    """
    Compute Morgan(2048, r=2) fingerprints for compounds with SMILES.
    If recompute=True, overwrite all fingerprints.
    Otherwise, only fill in rows where any fingerprint column is missing.
    Commits per chunk (no server cursor).

    :param batch: number of compounds to process per transaction
    :param recompute: whether to recompute fingerprints even if they exist
    :param radius: radius for Morgan fingerprint
    :param nbits: number of bits for Morgan fingerprint
    :return: total number of compounds updated
    """
    done = 0
    last_id = 0

    with SessionLocal() as s:
        while True:
            q = select(Compound).where(Compound.smiles.is_not(None))
            if not recompute:
                q = q.where(
                    or_(
                        Compound.fp_morgan_b2048_r2_bit.is_(None),
                        Compound.fp_morgan_b2048_r2_pop.is_(None),
                        Compound.fp_morgan_b2048_r2_vec.is_(None),
                    )
                )
            q = q.where(Compound.id > last_id).order_by(Compound.id).limit(batch)
            rows = s.scalars(q).all()

            if not rows:
                break

            for c in tqdm(rows, desc="Computing fingerprints"):
                bits, pop, vec = _morgan_bits_and_vec(c.smiles, radius=radius, nbits=nbits)
                c.fp_morgan_b2048_r2_bit = bits
                c.fp_morgan_b2048_r2_pop = pop
                c.fp_morgan_b2048_r2_vec = vec
                done += 1

            s.commit()
            last_id = rows[-1].id
            logger.info(f"[batch] committed {done} total updated")

    logger.info(f"Fingerprint backfill complete: {done} compounds updated")

    return done
