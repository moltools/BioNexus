"""Load NPAtlas JSON data into the BioNexus database."""

from __future__ import annotations
import logging

from sqlalchemy import select
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.dialects.postgresql import insert
from tqdm import tqdm

from bionexus.utils.io import iter_json
from bionexus.etl.chemistry import get_atom_counts
from bionexus.db.engine import SessionLocal
from bionexus.db.models import Annotation, Compound, CompoundRecord


logger = logging.getLogger(__name__)


def load_npatlas_file(path: str, chunk_size: int = 10000) -> tuple[int, int]:
    """
    Load NPAtlas data from a JSON file into the database.

    :param path: path to the NPAtlas JSON file
    :param chunk_size: number of records to process per database commit
    :return: tuple of (new_compounds, new_records, new_annotations) counts
    """
    new_compounds = 0
    new_records = 0
    new_annotations = 0

    batch_compounds: list[Compound] = []
    batch_records: list[CompoundRecord] = []

    # Collect annotations as dicts and upsert in bulk per chunk
    batch_ann_dicts: list[dict] = []  # each has {"_comp_obj", "scheme", "key", "value"}
    pending_ann = (
        0  # count staged annotations to trigger chunk commit even if only ann present
    )

    # De-dedupe within this load (by inchikey and by (source, ext_id))
    seen_inchikey: set[str] = set()
    seen_record_key: set[tuple[str, str]] = set()
    seen_annotation_key: set[tuple[str, str, str, str]] = (
        set()
    )  # (inchikey, scheme, key, value)

    with SessionLocal() as s:
        for d in tqdm(iter_json(path), desc="Loading NPAtlas"):
            # Unpack fields
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

            tax = d["origin_organism"]
            organism_type = tax["type"]
            organism_genus = tax["genus"]
            organism_species_name = tax["species"]
            organism_species = (
                f"{organism_genus} {organism_species_name}"
                if organism_genus and organism_species_name
                else None
            )

            npclassifier = d.get("npclassifier", None)
            if npclassifier:
                npclassifier_isglycoside: bool | None = npclassifier.get(
                    "isglycoside", None
                )
                npclassifier_class: list[str] = npclassifier.get("class_results", [])
                npclassifier_pathway: list[str] = npclassifier.get(
                    "pathway_results", []
                )
                npclassifier_superclass: list[str] = npclassifier.get(
                    "superclass_results", []
                )

            # Must have an inchikey to unify; skip if missing
            if not inchikey:
                continue

            # Batch level de-dupe of records (source, ext_id)
            rec_key = ("npatlas", npaid)
            if rec_key in seen_record_key:
                continue
            seen_record_key.add(rec_key)

            # Compute atom counts
            atom_counts = get_atom_counts(smiles)

            # 1) Find or create canonical compound by inchikey
            comp = s.scalars(
                select(Compound).where(Compound.inchikey == inchikey)
            ).first()
            if not comp:
                # Also avoid creating the same compound twice in one batch before flush
                if inchikey in seen_inchikey:
                    # If we have already staged it in this batch, fetch it from batch list
                    comp = next(
                        (c for c in batch_compounds if c.inchikey == inchikey), None
                    )
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
                        **(atom_counts or {}),
                    )
                    batch_compounds.append(comp)
                    seen_inchikey.add(inchikey)
                    new_compounds += 1

            else:
                # Optionally update canonical props to fill in gaps only
                if comp.smiles is None and smiles:
                    comp.smiles = smiles
                if comp.inchi is None and inchi:
                    comp.inchi = inchi
                if comp.mol_formula is None and mol_formula:
                    comp.mol_formula = mol_formula
                if comp.mol_weight is None and mol_weight:
                    comp.mol_weight = mol_weight
                if comp.exact_mass is None and exact_mass:
                    comp.exact_mass = exact_mass
                if comp.m_plus_h is None and m_plus_h:
                    comp.m_plus_h = m_plus_h
                if comp.m_plus_na is None and m_plus_na:
                    comp.m_plus_na = m_plus_na
                # update atom counts if missing
                if atom_counts:
                    for k, v in atom_counts.items():
                        if getattr(comp, k) is None and v is not None:
                            setattr(comp, k, v)

            # 2) Ensure a CompoundRecord (unique on source+ext_id)
            rec_exists = s.scalars(
                select(CompoundRecord).where(
                    CompoundRecord.source == "npatlas", CompoundRecord.ext_id == npaid
                )
            ).first()
            if not rec_exists:
                # If comp is new and not flushed yet, SQLA will link after add_all
                rec = CompoundRecord(
                    compound=comp,
                    source="npatlas",
                    ext_id=npaid,
                    name=name,
                    synonyms=synonyms or None,
                )
                batch_records.append(rec)
                new_records += 1

            # 3) Collect annotations
            annotations = (
                [
                    # (scheme, key, value)
                    ("taxonomy", "type", organism_type),
                    ("taxonomy", "genus", organism_genus),
                    ("taxonomy", "species", organism_species),
                ]
                + [
                    (
                        "npclassifier",
                        "isglycoside",
                        str(npclassifier_isglycoside).lower(),
                    )
                    if npclassifier_isglycoside is not None
                    else None,
                ]
                + [("npclassifier", "class", cls) for cls in npclassifier_class]
                + [("npclassifier", "pathway", pw) for pw in npclassifier_pathway]
                + [("npclassifier", "superclass", sc) for sc in npclassifier_superclass]
            )

            for scheme, key, value in annotations:
                if not value:
                    continue

                # Avoid annotation dupes in this batch
                ann_key = (inchikey, scheme, key, value)
                if ann_key in seen_annotation_key:
                    continue
                seen_annotation_key.add(ann_key)

                batch_ann_dicts.append(
                    {
                        "_comp_obj": comp,  # link to compound object to resolve after flush
                        "scheme": scheme,
                        "key": key,
                        "value": value,
                    }
                )
                pending_ann += 1

            # Commit in chunks
            if (len(batch_compounds) + len(batch_records) + pending_ann) >= chunk_size:
                try:
                    if batch_compounds:
                        s.add_all(batch_compounds)
                    if batch_records:
                        s.add_all(batch_records)
                    s.flush()  # assign comp.id (PKs) before inserting annotations

                    # Bulk upsert annotations
                    if batch_ann_dicts:
                        ann_values = [
                            {
                                "compound_id": d["_comp_obj"].id,
                                "scheme": d["scheme"],
                                "key": d["key"],
                                "value": d["value"],
                                "metadata_json": None,
                            }
                            for d in batch_ann_dicts
                        ]
                        if ann_values:
                            stmt = (
                                insert(Annotation)
                                .values(ann_values)
                                .on_conflict_do_nothing(
                                    constraint="uq_annotation_target_scheme_key_value"
                                )
                                .returning(Annotation.id)
                            )
                            inserted_ids = s.execute(stmt).scalars().all()
                            new_annotations += len(inserted_ids)

                    s.commit()

                except SQLAlchemyError as e:
                    s.rollback()
                    logger.error(f"Error during NPAtlas load chunk commit: {e}")
                    raise

                batch_compounds.clear()
                batch_records.clear()
                batch_ann_dicts.clear()
                seen_inchikey.clear()
                seen_record_key.clear()
                seen_annotation_key.clear()
                pending_ann = 0

        # Final flush
        if batch_compounds or batch_records or batch_ann_dicts:
            try:
                if batch_compounds:
                    s.add_all(batch_compounds)
                if batch_records:
                    s.add_all(batch_records)
                s.flush()  # assign comp.id (PKs) before inserting annotations

                # Bulk upsert annotations
                if batch_ann_dicts:
                    ann_values = [
                        {
                            "compound_id": d["_comp_obj"].id,
                            "scheme": d["scheme"],
                            "key": d["key"],
                            "value": d["value"],
                            "metadata_json": None,
                        }
                        for d in batch_ann_dicts
                    ]
                    if ann_values:
                        stmt = (
                            insert(Annotation)
                            .values(ann_values)
                            .on_conflict_do_nothing(
                                constraint="uq_annotation_target_scheme_key_value"
                            )
                            .returning(Annotation.id)
                        )
                        inserted_ids = s.execute(stmt).scalars().all()
                        new_annotations += len(inserted_ids)

                s.commit()

            except SQLAlchemyError as e:
                s.rollback()
                logger.error(f"Error during NPAtlas load final commit: {e}")
                raise

    return new_compounds, new_records, new_annotations
