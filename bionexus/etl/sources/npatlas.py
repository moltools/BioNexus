from __future__ import annotations
from sqlalchemy import select
from tqdm import tqdm
from bionexus.utils.io import iter_json
from bionexus.db.engine import SessionLocal
from bionexus.db.models import Compound

def load_npatlas_file(path: str, chunk_size: int = 10000) -> int:
    inserted = 0

    with SessionLocal() as s:
        batch: list[Compound] = []
        seen: set[str] = set()  # batch level

        for d in tqdm(iter_json(path)):
            # unpack fields
            npaid = d["npaid"]
            name = d["original_name"]
            synonyms = [x["name"] for x in d["synonyms"]]
            mol_formula = d["mol_formula"]
            mol_weight = d["mol_weight"]
            exact_mass = d["exact_mass"]
            inchi = d["inchi"]
            inchikey = d["inchikey"]
            smiles = d["smiles"]
            m_plus_h = d["m_plus_h"]
            m_plus_na = d["m_plus_na"]
            
            # dedupe on batch level
            if inchikey in seen:
                continue
            seen.add(inchikey)

            # avoid duplicates already in DB
            exists = s.execute(select(Compound.id).where(Compound.inchikey == inchikey)).first()
            if exists:
                continue

            c = Compound(
                source="npatlas",
                ext_id=npaid,
                name=name,
                synonyms=synonyms,
                mol_formula=mol_formula,
                mol_weight=mol_weight,
                exact_mass=exact_mass,
                smiles=smiles,
                inchi=inchi,
                inchikey=inchikey,
                m_plus_h=m_plus_h,
                m_plus_na=m_plus_na,
            )
            batch.append(c)

            if len(batch) >= chunk_size:
                s.add_all(batch)
                s.commit()
                inserted += len(batch)
                batch.clear()

        if batch:
            s.add_all(batch)
            s.commit()
            inserted += len(batch)

    return inserted
