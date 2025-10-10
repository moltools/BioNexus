from __future__ import annotations
import json
from sqlalchemy import select
from tqdm import tqdm
from bionexus.etl.chemistry import smiles_to_inchikey
from bionexus.db.engine import SessionLocal
from bionexus.db.models import Compound, CompoundRecord

def load_mibig_files(json_paths: list[str] | None, gbk_paths: str | list[str], chunk_size: int = 1000) -> tuple[int, int]:
    new_compounds = 0
    new_bgcs = 0
    batch_compounds: list[Compound] = []
    batch_records: list[CompoundRecord] = []

    # de-dupe within this load (by inchikey and by (source, ext_id))
    seen_inchikey: set[str] = set()
    seen_record_key: set[tuple[int, str, str]] = set()

    with SessionLocal() as s:
        for path in tqdm(json_paths or [], desc="Loading MIBiG JSON"):
            data = json.load(open(path))
            
            # unpack fields
            accession = data["accession"]
            version = data["version"]
            quality = data["quality"]
            status = data["status"]
            completeness = data["completeness"]
            compounds = data.get("compounds", [])

            ext_id = f"{accession}.{version}"

            for c in compounds:
                name = c.get("name")
                smiles = c.get("structure")
                mol_formula = c.get("formula")
                exact_mass = c.get("mass")

                if smiles:
                    inchikey: str | None = smiles_to_inchikey(smiles)

                    # must have an inchikey to unify; skip if missing
                    if not inchikey:
                        continue
                    
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
                                mol_formula=mol_formula,
                                exact_mass=exact_mass,
                            )
                            batch_compounds.append(comp)
                            seen_inchikey.add(inchikey)
                            new_compounds += 1
                    
                    else:
                        # optionally update canonical props to fill in gaps only
                        if comp.smiles is None and smiles: comp.smiles = smiles
                        if comp.mol_formula is None and mol_formula: comp.mol_formula = mol_formula
                        if comp.exact_mass is None and exact_mass: comp.exact_mass = exact_mass

                    # batch level de-dupe of records (compound_id, source, ext_id)
                    record_key = (comp.id, "mibig", ext_id)
                    if record_key in seen_record_key:
                        continue
                    seen_record_key.add(record_key)

                    # 2) ensure a CompoundRecord (unique on (compound_id, source, ext_id))
                    rec_exists = s.scalars(
                        select(CompoundRecord).where(
                            CompoundRecord.compound_id == comp.id,
                            CompoundRecord.source == "mibig",
                            CompoundRecord.ext_id == ext_id
                        )
                    ).first()
                    if not rec_exists:
                        # if comp is new and not flushed yet, SQLA will link after add_all
                        rec = CompoundRecord(
                            compound=comp,
                            source="mibig",
                            ext_id=ext_id,
                            name=name,
                        )
                        batch_records.append(rec)

                    # commit in chunks
                    if len(batch_compounds) + len(batch_records) >= chunk_size:
                        s.add_all(batch_compounds)
                        s.add_all(batch_records)
                        s.commit()
                        batch_compounds.clear()
                        batch_records.clear()
                        seen_inchikey.clear()  # safe to clear after flush
                        seen_record_key.clear()  # safe to clear after flush

        # final flush      
        if batch_compounds: s.add_all(batch_compounds)
        if batch_records: s.add_all(batch_records)
        s.commit()

    return new_compounds, new_bgcs