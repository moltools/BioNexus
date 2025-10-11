from __future__ import annotations
import logging
import requests
from tqdm import tqdm
from sqlalchemy import select
from sqlalchemy.exc import SQLAlchemyError
from bionexus.db.models import Compound, Annotation
from bionexus.db.engine import SessionLocal

logger = logging.getLogger(__name__)

API_URL: str = "https://npclassifier.gnps2.org/classify?smiles={}"

def query_npclassifier(smiles: str) -> dict | None:
    """Query the NPClassifier API with a SMILES string."""
    try:
        response = requests.get(API_URL.format(smiles), timeout=10)
        response.raise_for_status()
        annotations = response.json()
        return annotations
    except (requests.RequestException, ValueError) as e:
        logger.error(f"Error querying NPClassifier for SMILES: {smiles}\n{e}")
        return None

def annotate_with_npclassifier(recompute: bool, chunk_size: int = 10000) -> int:
    num_annotated = 0

    with SessionLocal() as s:
        # 1) determine which compounds need annotation
        if not recompute:
            # subquery: does an npclassifier annotation exist for this compound?
            npclass_exists = (
                select(1)
                .where(
                    (Annotation.compound_id == Compound.id) &
                    (Annotation.scheme == "npclassifier")
                ).exists()
            )
            # main query: compounds for which no npclassifier annotation exists
            q = (
                s.scalars(
                    select(Compound.id)
                    .where(~npclass_exists)
                )
            )
            compound_ids = q.all()
            logger.info(f"Found {len(compound_ids)} compounds without NPClassifier annotations")
        else:
            # get all compound IDs
            compound_ids = s.scalars(select(Compound.id)).all()
            logger.info(f"Recomputing NPClassifier annotations for all {len(compound_ids)} compounds")
    
        # 2) retrieve annotations using NPClassifier API, and store in DB
        for i in tqdm(range(0, len(compound_ids), chunk_size), desc="Annotating chunks"):
            try:
                chunk = compound_ids[i:i + chunk_size]
                compounds = s.scalars(select(Compound).where(Compound.id.in_(chunk))).all()
                for compound in tqdm(compounds, desc="Processing compounds in chunk", leave=False):
                    if not compound.smiles:
                        continue
                    annotations = query_npclassifier(compound.smiles)
                    if not annotations:
                        continue
                    
                    # remove existing NPClassifier annotations if recompute is True
                    if recompute:
                        s.query(Annotation).filter(
                            (Annotation.compound_id == compound.id) &
                            (Annotation.scheme == "npclassifier")
                        ).delete(synchronize_session="fetch")
                    
                    # add new annotations
                    if annotations:
                        npclassifier_isglycoside: bool | None = annotations.get("isglycoside", None)
                        npclassifier_class: list[str] = annotations.get("class_results", [])
                        npclassifier_pathway: list[str] = annotations.get("pathway_results", [])
                        npclassifier_superclass: list[str] = annotations.get("superclass_results", [])

                    # [(scheme, key, value), ...]
                    to_add = [
                        ("npclassifier", "isglycoside", str(npclassifier_isglycoside).lower()) if npclassifier_isglycoside is not None else None,
                    ] + [
                        ("npclassifier", "class", cls) for cls in npclassifier_class
                    ] + [
                        ("npclassifier", "pathway", pw) for pw in npclassifier_pathway
                    ] + [
                        ("npclassifier", "superclass", sc) for sc in npclassifier_superclass
                    ]

                    for scheme, key, value in to_add:
                        ann = Annotation(
                            compound_id=compound.id,
                            scheme=scheme,
                            key=key,
                            value=value
                        )
                        s.add(ann)
                        num_annotated += 1
                
                s.commit()
                logger.info(f"Processed chunk {i // chunk_size + 1}, total annotated so far: {num_annotated}")

            except SQLAlchemyError as e:
                s.rollback()
                logger.error(f"Database error while processing chunk starting at index {i}: {e}")
                raise

    return num_annotated
