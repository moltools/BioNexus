from __future__ import annotations
import logging
from sqlalchemy import select, or_
from bionexus.db.engine import SessionLocal
from bionexus.db.models import Compound
from tqdm import tqdm

logger = logging.getLogger(__name__)

def smiles_to_inchikey(smiles: str) -> str | None:
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

def _morgan_bits_and_vec(smiles: str, radius: int = 2, nbits: int = 2048) -> tuple[str | None, int | None, list[float] | None]:
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