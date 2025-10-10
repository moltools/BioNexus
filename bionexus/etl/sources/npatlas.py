from __future__ import annotations
from sqlalchemy import select
from tqdm import tqdm
from bionexus.utils.io import iter_json
from bionexus.db.engine import SessionLocal
from bionexus.db.models import Compound, CompoundRecord

def load_npatlas_file(path: str, chunk_size: int = 10000) -> int:
    new_compounds = 0
    batch_compounds: list[Compound] = []
    batch_records: list[CompoundRecord] = []

    # de-dedupe within this load (by inchikey and by (source, ext_id))
    seen_inchikey: set[str] = set()
    seen_record_key: set[tuple[str, str]] = set()

    with SessionLocal() as s:
        for d in tqdm(iter_json(path), desc="Loading NPAtlas"):
            # unpack fields
            npaid = d["npaid"]
            inchikey = d["inchikey"]
            inchi = d["inchi"]
            smiles = d["smiles"]
            name = d["original_name"]
            synonyms = [x["name"] for x in d["synonyms"]]

            mol_formula = d["mol_formula"]
            mol_weight = d["mol_weight"]
            exact_mass = d["exact_mass"]
            m_plus_h = d["m_plus_h"]
            m_plus_na = d["m_plus_na"]

            # must have an inchikey to unify; skip if missing
            if not inchikey:
                continue

            # batch level de-dupe of records (source, ext_id)
            rec_key = ("npatlas", npaid)
            if rec_key in seen_record_key:
                continue
            seen_record_key.add(rec_key)

            # 1) find or create canonical compound by inchikey
            comp = s.scalars(select(Compound).where(Compound.inchikey == inchikey)).first()
            if not comp:
                # also avoid creating the same compound twice in one batch before flush
                if inchikey in seen_inchikey:
                    # if we have already staged it in this batch, fit it from batch list
                    comp = next((c for c in batch_compounds if c.inchikey == inchikey), None)
                if not comp:
                    comp = Compound(
                        inchikey=inchikey,
                        smiles=smiles,
                        inchi=inchi,
                        mol_formula=mol_formula,
                        mol_weight=mol_weight,
                        exact_mass=exact_mass,
                        m_plus_h=m_plus_h,
                        m_plus_na=m_plus_na,
                    )
                    batch_compounds.append(comp)
                    seen_inchikey.add(inchikey)
                    new_compounds += 1
            
            else:
                # optionally update canonical props to fill in gaps only
                if comp.smiles is None and smiles: comp.smiles = smiles
                if comp.inchi is None and inchi: comp.inchi = inchi
                if comp.mol_formula is None and mol_formula: comp.mol_formula = mol_formula
                if comp.mol_weight is None and mol_weight: comp.mol_weight = mol_weight
                if comp.exact_mass is None and exact_mass: comp.exact_mass = exact_mass
                if comp.m_plus_h is None and m_plus_h: comp.m_plus_h = m_plus_h
                if comp.m_plus_na is None and m_plus_na: comp.m_plus_na = m_plus_na

            # 2) ensure a CompoundRecord (unique on source+ext_id)
            rec_exists = s.scalars(
                select(CompoundRecord).where(
                    CompoundRecord.source == "npatlas",
                    CompoundRecord.ext_id == npaid
                )
            ).first()
            if not rec_exists:
                # if comp is new and not flushed yet, SQLA will link after add_all
                rec = CompoundRecord(
                    compound=comp,
                    source="npatlas",
                    ext_id=npaid,
                    name=name,
                    synonyms=synonyms or None,
                )
                batch_records.append(rec)

            # commit in chunks
            if len(batch_compounds) + len(batch_records) >= chunk_size:
                if batch_compounds: s.add_all(batch_compounds)
                if batch_records: s.add_all(batch_records)
                s.commit()
                batch_compounds.clear()
                batch_records.clear()
                seen_inchikey.clear()  # safe to clear after flush
                seen_record_key.clear()  # safe to clear after flush

        # final flush
        if batch_compounds: s.add_all(batch_compounds)
        if batch_records: s.add_all(batch_records)
        s.commit()

    return new_compounds
