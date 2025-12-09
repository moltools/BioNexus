"""BioCracker parsing of GBK files."""

import logging
import tempfile
import json
from pathlib import Path
from dataclasses import asdict
from typing import Any

from tqdm import tqdm

from bionexus.config import default_cache_dir
from bionexus.db.engine import SessionLocal
from bionexus.db.models import GenBankRegion, BioCrackerGenBank

logger = logging.getLogger(__name__)


def to_jsonable(obj: dict) -> dict:
    """
    Convert a BioCracker readout object to a JSON-serializable dictionary.
    
    :param obj: BioCracker readout object
    :return: JSON-serializable dictionary
    """
    rec = obj.get("rec", None)
    readout = obj.get("readout", None)

    if rec is None or readout is None:
        raise ValueError("Invalid BioCracker readout object")
    
    rec_as_dict = rec.to_dict()
    readout_as_dict = [asdict(x) for x in readout]

    return {
        "rec": rec_as_dict,
        "readout": readout_as_dict,
    }


def parse_gbks_with_biocracker(
    cache_dir: Path | str,
    recompute: bool = False,
    chunk_size: int = 1000,
    workers: int = 1,
    to_jsonl: bool = False,
) -> int:
    """
    Parse GBK files using BioCracker and store results.

    :param cache_dir: directory to use for caching input/output files
    :param recompute: if True, recompute results for all compounds
    :param chunk_size: number of results to upsert in each database transaction
    :param workers: number of parallel workers to use for RetroMol processing
    :param to_jsonl: if True, output results to JSONL format into cache_dir
    """
    if cache_dir is None:
        cache_dir = default_cache_dir()

    if isinstance(cache_dir, str):
        cache_dir = Path(cache_dir)

    # Import BioCracker components
    try:
        from biocracker.antismash import parse_region_gbk_file
        from biocracker.readout import linear_readouts as biocracker_linear_readouts
    except ImportError:
        logger.error("biocracker package not found")
        return
    
    jsonl_file = None
    if to_jsonl:
        cache_dir.mkdir(parents=True, exist_ok=True)
        jsonl_file = (cache_dir / "biocracker_readouts.jsonl").open("w" , encoding="utf-8")

    processed = 0
    
    with SessionLocal() as s:
        # base query for counting and id pagination
        base_query = s.query(GenBankRegion.id)

        if not recompute:
            base_query = base_query.outerjoin(
                BioCrackerGenBank,
                GenBankRegion.id == BioCrackerGenBank.genbank_region_id,
            ).filter(BioCrackerGenBank.id.is_(None))

        # Report on total items to process
        total_to_process = base_query.count()
        logger.info(f"Processing {total_to_process} GenBankRegion entries with BioCracker")

        if total_to_process == 0:
            if jsonl_file is not None:
                jsonl_file.close()
            return 0
        
        pbar = tqdm(total=total_to_process, desc="BioCracker GBK Parsing")

        last_id: int | None = None

        while True:

            # keyset paginate over IDs: id > last_id
            q_ids = base_query
            if last_id is not None:
                q_ids = q_ids.filter(GenBankRegion.id > last_id)

            id_rows = (
                q_ids.order_by(GenBankRegion.id)
                .limit(chunk_size)
                .all()
            )
            ids = [r.id for r in id_rows]

            if not ids:
                break # no more IDs to process

            # load GenBankRegion objects for this chunk
            regions = (
                s.query(GenBankRegion)
                .filter(GenBankRegion.id.in_(ids))
                .order_by(GenBankRegion.id)
                .all()
            )

            for region in regions:
                with tempfile.NamedTemporaryFile(delete=True, suffix=".gbk") as temp_gbk:
                    temp_gbk.write(region.gbk_text.encode("utf-8"))
                    temp_gbk.flush()
                    gbk_path = Path(temp_gbk.name)
                    targets = parse_region_gbk_file(gbk_path, top_level="cand_cluster")

                # Build JSONable structure per region
                region_result: dict[str, Any] = {
                    "genbank_region": {
                        "id": region.id,
                        "source": region.source,
                        "ext_id": region.ext_id,
                        "size_bytes": region.size_bytes,
                        "sha256": region.sha256,
                    },
                    "targets": [],
                }
                
                for target in targets:
                    gene_readouts = []
                    for readout in biocracker_linear_readouts(
                        target,
                        model=None,
                        cache_dir_override=cache_dir,
                        level="gene",
                        pred_threshold=0.1
                    ):
                        # JSON serialize the readout
                        readout_as_dict = to_jsonable(readout)
                        gene_readouts.append(readout_as_dict)
                    region_result["targets"].append(gene_readouts)

                # Upsert into BioCrackerGenBank
                existing = (
                    s.query(BioCrackerGenBank)
                    .filter(BioCrackerGenBank.genbank_region_id == region.id)
                    .one_or_none()
                )

                if existing is None:
                    obj = BioCrackerGenBank(
                        genbank_region_id=region.id,
                        result_json=region_result,
                    )
                    s.add(obj)
                else:
                    if recompute:
                        existing.result_json = region_result

                # Optional JSONL write
                if jsonl_file is not None:
                    jsonl_file.write(json.dumps(region_result) + "\n")

                processed += 1
                pbar.update(1)

            # Commit once per chunk
            s.commit()
            last_id = ids[-1]

        pbar.close()

    if jsonl_file is not None:
        jsonl_file.close()

    logger.info(f"Processed {processed} GenBankRegion entries with BioCracker")

    return processed


