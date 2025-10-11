from __future__ import annotations
import json, hashlib, logging
from sqlalchemy import select
from sqlalchemy.exc import SQLAlchemyError
from tqdm import tqdm
from bionexus.config import DEFAULT_MAX_BYTES_GBK
from bionexus.etl.chemistry import smiles_to_inchikey
from bionexus.db.engine import SessionLocal
from bionexus.db.models import Compound, CompoundRecord, GenBankRegion

logger = logging.getLogger(__name__)

def _read_gbk_text(path: str) -> str:
    # allow .gbk, .gb, .gbff etc.
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        return fh.read()

def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()

def load_mibig_files(json_paths: list[str] | None, gbk_paths: str | list[str], chunk_size: int = 1000) -> tuple[int, int, int, int]:
    source = "mibig"

    new_compounds = 0
    new_records = 0
    new_regions = 0
    new_annotations = 0

    batch_compounds: list[Compound] = []
    batch_records: list[CompoundRecord] = []
    batch_regions: list[GenBankRegion] = []

    # collect annotations as dicts and upsert in bulk per chunk
    batch_ann_dicts: list[dict] = []  # each has {"_comp_obj"/"_gbk_obj", "scheme", "key", "value"}
    pending_ann = 0  # count staged annotations to trigger chunk commit even if only ann present

    # de-dupe within this load (by inchikey and by (inchikey, source, ext_id) pre-flush)
    seen_inchikey: set[str] = set()
    seen_record_key: set[tuple[str, str, str]] = set()  # (inchikey, source, ext_id)
    seen_region_key: set[str] = set()  # (sha256)
    seen_annotation_key: set[tuple[str, str, str, str]] = set()  # (inchikey/sha256, scheme, key, value)

    with SessionLocal() as s:
        # --- compounds from JSON files ----------------------------------------
        for path in tqdm(json_paths or [], desc="Loading MIBiG JSONs"):
            data = json.load(open(path))
            ext_id = data["accession"]

            for c in data.get("compounds", []):
                name = c.get("name")
                smiles = c.get("structure")
                mol_formula = c.get("formula")
                exact_mass = c.get("mass")

                if not smiles:
                    continue

                if smiles:
                    inchikey: str | None = smiles_to_inchikey(smiles)

                    # must have an inchikey to unify; skip if missing
                    if not inchikey:
                        continue
                    
                    # 1) find or stage canonical compound by inchikey
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

                    # 2) ensure CompoundRecord, batch de-dupe by inchikey+source+ext_id (pre-flush)
                    rk = (inchikey, source, ext_id)
                    if rk in seen_record_key:
                        pass
                    else:
                        seen_record_key.add(rk)
                        rec_exists = s.scalars(
                            select(CompoundRecord).where(
                                CompoundRecord.compound_id == comp.id,  # None for new; query will skip
                                CompoundRecord.source == source,
                                CompoundRecord.ext_id == ext_id
                            )
                        ).first()

                        if not rec_exists:
                            batch_records.append(
                                CompoundRecord(
                                    compound=comp,  # safe even if comp not flushed
                                    source=source,
                                    ext_id=ext_id,
                                    name=name,
                                )
                            )
                            new_records += 1

                    # flush if batch full
                    if (len(batch_compounds) + len(batch_records) + len(batch_regions) + pending_ann) >= chunk_size:
                        try:
                            if batch_compounds: s.add_all(batch_compounds)
                            if batch_records: s.add_all(batch_records)
                            if batch_regions: s.add_all(batch_regions)
                            s.flush()  # assign ids to new objects for FKs in annotations

                            # bulk upsert annotations
                            # TODO

                            s.commit()

                        except SQLAlchemyError as e:
                            s.rollback()
                            logger.error(f"Error during MIBiG load chunk commit: {e}")
                            raise

                        batch_compounds.clear()
                        batch_records.clear()
                        batch_regions.clear()
                        batch_ann_dicts.clear()
                        seen_inchikey.clear()
                        seen_record_key.clear()
                        seen_region_key.clear()
                        seen_annotation_key.clear()
                        pending_ann = 0

        # --- compounds from GenBank files -------------------------------------
        for path in tqdm(gbk_paths or [], desc="Loading MIBiG GenBank files"):
            ext_id = path.split("/")[-1].replace(".gbk", "").replace(".gb", "")

            gbk_text = _read_gbk_text(path)
            enc = gbk_text.encode("utf-8", errors="replace")
            size_bytes = len(enc)
            sha256 = _sha256_bytes(enc)
            
            # skip if already seen in this batch
            if sha256 in seen_region_key:
                continue
            seen_region_key.add(sha256)

            # skip if too large
            if size_bytes > DEFAULT_MAX_BYTES_GBK:
                logger.warning(f"Skipping {ext_id} ({size_bytes // 1_000_000:.2f} MB > {DEFAULT_MAX_BYTES_GBK // 1_000_000:.2f} MB max)")
                continue

            # prefer dedup by sha256 when present
            by_hash = s.scalars(select(GenBankRegion).where(GenBankRegion.sha256 == sha256)).first()
            if by_hash:
                continue
            
            # also check uniqueness on (source, ext_id)
            by_key = s.scalars(
                select(GenBankRegion).where(
                    GenBankRegion.source == source,
                    GenBankRegion.ext_id == ext_id,
                )
            ).first()
            if by_key:
                continue
            
            batch_regions.append(
                GenBankRegion(
                    source=source,
                    ext_id=ext_id,
                    gbk_text=gbk_text,
                    size_bytes=size_bytes,
                    sha256=sha256,
                )
            )
            new_regions += 1

            if (len(batch_compounds) + len(batch_records) + len(batch_regions) + pending_ann) >= chunk_size:
                try:
                    if batch_compounds: s.add_all(batch_compounds)
                    if batch_records: s.add_all(batch_records)
                    if batch_regions: s.add_all(batch_regions)
                    s.flush()  # assign ids to new objects for FKs in annotations

                    # bulk upsert annotations
                    # TODO

                    s.commit()

                except SQLAlchemyError as e:
                    s.rollback()
                    logger.error(f"Error during MIBiG load chunk commit: {e}")
                    raise

                batch_compounds.clear()
                batch_records.clear()
                batch_regions.clear()
                batch_ann_dicts.clear()
                seen_inchikey.clear()
                seen_record_key.clear()
                seen_region_key.clear()
                seen_annotation_key.clear()
                pending_ann = 0

        # final flush
        if batch_compounds or batch_records or batch_regions or batch_ann_dicts:
            try:
                if batch_compounds: s.add_all(batch_compounds)
                if batch_records: s.add_all(batch_records)
                if batch_regions: s.add_all(batch_regions)
                s.flush()  # assign ids to new objects for FKs in annotations

                # bulk upsert annotations
                # TODO

                s.commit()

            except SQLAlchemyError as e:
                s.rollback()
                logger.error(f"Error during MIBiG final load commit: {e}")
                raise

    return new_compounds, new_records, new_regions, new_annotations