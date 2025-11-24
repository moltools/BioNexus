"""Module for performing Jaccard similarity searches on chemical compound fingerprints."""

from __future__ import annotations
import logging
from pathlib import Path
from typing import Literal

import numpy as np
from pgvector.sqlalchemy import Vector
from sqlalchemy import Float, bindparam, text
from sqlalchemy.dialects.postgresql import BIT

from bionexus.db.engine import SessionLocal


logger = logging.getLogger(__name__)


def jaccard_search_exact(bitstr: str, top_k: int = 50) -> list[dict]:
    """
    Perform an exact Jaccard similarity search on compound fingerprints.

    :param bitstr: the query fingerprint as a bit string
    :param top_k: the number of top results to return
    :return: a list of dictionaries containing compound information and Jaccard similarity scores
    """
    # limit top-k to 500 for performance reasons
    if top_k > 500:
        logger.warning("top_k limited to 500 for performance reasons")
    top_k = min(top_k, 500)

    sql = text("""
    WITH top_compounds AS (
        SELECT 
            c.id,
            c.smiles,
            CASE
                WHEN length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')) = 0
                THEN 1.0
                ELSE length(replace((c.fp_morgan_b2048_r2_bit & :qb)::text, '0',''))::float
                    / NULLIF(length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')), 0)
            END AS jacc
        FROM compound c
        WHERE c.fp_morgan_b2048_r2_bit IS NOT NULL
        ORDER BY jacc DESC
        LIMIT :k
    )
    SELECT
        t.id,
        r.source,
        r.name,
        t.jacc,
        t.smiles
    FROM top_compounds t
    LEFT JOIN compound_record r
        ON r.compound_id = t.id
    WHERE r.name IS NOT NULL
    ORDER BY t.jacc DESC, r.source, r.name, t.id
    """).bindparams(bindparam("qb", type_=BIT(2048)), bindparam("k"))

    with SessionLocal() as s:
        rows = s.execute(sql, {"qb": bitstr, "k": top_k}).mappings().all()
        return rows


def jaccard_search_hybrid(bitstr: str, qvec: list[float], top_k: int = 50, cand_k: int = 2000) -> list[dict]:
    """
    Perform a hybrid Jaccard similarity search on compound fingerprints using an initial
    candidate selection based on vector similarity.

    :param bitstr: the query fingerprint as a bit string
    :param qvec: the query fingerprint as a vector of floats
    :param top_k: the number of top results to return
    :param cand_k: the number of candidate compounds to consider
    :return: a list of dictionaries containing compound information and Jaccard similarity scores
    """
    # limit top-k to 500 for performance reasons
    if top_k > 500:
        logger.warning("top_k limited to 500 for performance reasons")
    top_k = min(top_k, 500)

    sql = text("""
    WITH cand AS (
      SELECT c.id, c.smiles, c.fp_morgan_b2048_r2_bit
      FROM compound c
      WHERE c.fp_morgan_b2048_r2_bit IS NOT NULL
      ORDER BY c.fp_morgan_b2048_r2_vec <#> :qv   -- inner-product distance
      LIMIT :cand_k
    ),
    scored AS (
      SELECT
        c.id,
        c.smiles,
        CASE
          WHEN length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')) = 0
          THEN 1.0
          ELSE length(replace((c.fp_morgan_b2048_r2_bit & :qb)::text, '0',''))::float
               / NULLIF(length(replace((c.fp_morgan_b2048_r2_bit | :qb)::text, '0','')), 0)
        END AS jacc
      FROM cand c
    ),
    top_compounds AS (
      SELECT id, smiles, jacc
      FROM scored
      ORDER BY jacc DESC
      LIMIT :k
    )
    SELECT
      t.id,
      t.smiles,
      r.source,
      r.name,
      t.jacc
    FROM top_compounds t
    LEFT JOIN compound_record r
      ON r.compound_id = t.id
    WHERE r.name IS NOT NULL
    ORDER BY t.jacc DESC, r.source, r.ext_id
    """).bindparams(
        bindparam("qb", type_=BIT(2048)),
        bindparam("qv", type_=Vector(2048)),
        bindparam("k"),
        bindparam("cand_k"),
    )

    with SessionLocal() as s:
        return s.execute(sql, {"qb": bitstr, "qv": qvec, "k": top_k, "cand_k": cand_k}).mappings().all()


def retro_search_compound(
    smiles: str,
    top_k: int = 50,
    counted: bool = False
) -> list[dict]:
    """
    Perform a RetroMol fingerprint similarity search.

    :param smiles: the query compound SMILES
    :param top_k: the number of top results to return
    :param counted: whether to use counted vector similarity or binary vector similarity
    :return: a list of dictionaries containing compound information and similarity scores
    """
    # Import RetroMol methods
    try:
        from retromol.api import run_retromol_with_timeout
        from retromol.io import Input as RetroMolInput
        from retromol.fingerprint import (
            FingerprintGenerator,
            NameSimilarityConfig,
            polyketide_family_of
        )
        from retromol.rules import get_path_default_matching_rules
    except ImportError:
        logger.error("RetroMol library is not installed. Please install it to use retro_search.")
        return []
    
    raise RuntimeError("Please update generator setup to match recent changes in retromol ETL code; centralize generator setup?")

    # limit top-k to 500 for performance reasons
    if top_k > 500:
        logger.warning("top_k limited to 500 for performance reasons")
    top_k = min(top_k, 500)

    # Parse input compound
    input_cmp = RetroMolInput(cid="input", repr=smiles)
    result = run_retromol_with_timeout(input_cmp)
    coverage = result.best_total_coverage()
    logger.info(f"RetroMol coverage for input compound: {coverage:.2%}")

    # Setup generator
    collapse_by_name = ["glycosylation", "methylation"]
    cfg = NameSimilarityConfig(family_of=polyketide_family_of, symmetric=True, family_repeat_scale=1)
    generator = FingerprintGenerator(
        matching_rules_yaml=get_path_default_matching_rules(),
        collapse_by_name=collapse_by_name,
        name_similarity=cfg
    )
    logger.info(f"Initialized RetroMol FingerprintGenerator: {generator}")

    # Calculate fingerprints of shape (N, 512)
    fps = generator.fingerprint_from_result(result, num_bits=512, counted=counted)
    if fps is None or fps.shape[0] == 0:
        logger.warning("No RetroMol fingerprints could be generated for the input compound.")
        return []
    
    vec_col = "fp_retro_b512_vec_counted" if counted else "fp_retro_b512_vec_binary"

    # Cosine distance operator in pgvector is '<=>'; cosine similarity = 1 - distance
    sql = text(f"""
        SELECT
            rf.id AS id,
            cr.source AS source,
            cr.name AS name,
            (1.0 - (rf.{vec_col} <=> :qv)) AS cosine,
            c.smiles AS smiles
        FROM retrofingerprint rf
        JOIN retromol_compound rmc
          ON rmc.id = rf.retromol_compound_id
        JOIN compound c
          ON c.id = rmc.compound_id
        LEFT JOIN compound_record cr
          ON cr.compound_id = c.id
        ORDER BY rf.{vec_col} <=> :qv, rf.id
        LIMIT :k
    """).bindparams(
        bindparam("qv", type_=Vector(512)),
        bindparam("k")
    )

    # Collect top_k result for every individiual fingerprint, then merge and sort and return top_k overall
    best_by_id: dict[int, dict] = {}
    with SessionLocal() as s:
        for fp in fps:
            qv = [float(x) for x in fp]
            rows = s.execute(sql, {"qv": qv, "k": top_k}).mappings().all()
            for r in rows:
                rid = r["id"]
                source = r["source"]
                cos = float(r["cosine"]) if r["cosine"] is not None else None
                if rid not in best_by_id or (cos is not None and cos > best_by_id[rid]["cosine"]):
                    best_by_id[(rid, source)] = {
                        "id": r["id"],
                        "source": r["source"],
                        "name": r["name"],
                        "cosine": cos,
                        "smiles": r["smiles"],
                    }

    results = sorted(
        best_by_id.values(),
        key=lambda d: (d["cosine"] is not None, d["cosine"]),
        reverse=True
    )[:top_k]

    return results


def retro_search_gbk(
    path: Path | str,
    top_k: int = 50,
    readout_toplevel: Literal["region", "cand"] = "region",
    readout_sublevel: Literal["rec", "gene"] = "rec",
    counted: bool = False,
    cache_dir: Path | str | None = None,
    kmer_sizes: list[int] | None = None,
    metric: Literal["cosine", "tanimoto"] = "tanimoto",
) -> list[dict]:
    """
    Perform a RetroMol fingerprint similarity search on compounds in a GenBank file.

    :param path: path to the GenBank file
    :param top_k: the number of top results to return per compound
    :param readout_toplevel: whether to read out fingerprints at 'region' or 'cand' level
    :param readout_sublevel: whether to read out fingerprints at 'rec' or 'gene' level
    :param counted: whether to use counted vector similarity or binary vector similarity
    :param cache_dir: optional path to a cache directory for intermediate files
    :param kmer_sizes: optional list of k-mer sizes to use for fingerprint generation
    :param metric: similarity metric to use ('cosine' or 'tanimoto')
    :return: a list of dictionaries containing compound information and similarity scores
    """
    try:
        from biocracker.antismash import parse_region_gbk_file
        from biocracker.config import LOGGER_LEVEL, LOGGER_NAME
        from biocracker.readout import NRPSModuleReadout, PKSModuleReadout, linear_readouts
        from biocracker.text_mining import get_default_tokenspecs, mine_virtual_tokens
        from retromol.fingerprint import (
            FingerprintGenerator,
            NameSimilarityConfig,
            get_kmers,
            polyketide_family_of,
            # polyketide_ancestors_of
        )
        from retromol.rules import get_path_default_matching_rules
    except ImportError as e:
        logger.error("biocracker and/or retromol library is not installed. Please install it to use retro_search_gbk.")
        return []

    if isinstance(path, str):
        path = Path(path)

    if cache_dir is not None:
        cache_dir_override = Path(cache_dir)
    else:
        cache_dir_override = None

    if kmer_sizes is None:
        kmer_sizes = [1, 2, 3]

    # TODO: actually retrieve correction version of matching rules; only compare to fingerprints parsed with the same rules version

    # Setup generator
    collapse_by_name = [
        "glycosylation", "methylation",
        *[f"{x}{i}" for x in "ABCDEF" for i in range(1,16)]
    ]
    cfg = NameSimilarityConfig(family_of=polyketide_family_of, symmetric=True, family_repeat_scale=1)
    generator = FingerprintGenerator(
        matching_rules_yaml=get_path_default_matching_rules(),
        collapse_by_name=collapse_by_name,
        name_similarity=cfg
    )
    logger.info(f"Initialized RetroMol FingerprintGenerator: {generator}")

    if metric == "tanimoto":
        bit_col = "fp_retro_b512_bit"
    else:
        vec_col = "fp_retro_b512_vec_counted" if counted else "fp_retro_b512_vec_binary"

    # Build SQL query
    if metric == "cosine":
        sql = text(f"""
            SELECT
                rf.id AS id,
                cr.source AS source,
                cr.name AS name,
                (1.0 - (rf.{vec_col} <=> :qv)) AS score,
                (1.0 - (rf.{vec_col} <=> :qv)) AS cosine,
                c.smiles AS smiles
            FROM retrofingerprint rf
            JOIN retromol_compound rmc
            ON rmc.id = rf.retromol_compound_id
            JOIN compound c
            ON c.id = rmc.compound_id
            LEFT JOIN compound_record cr
            ON cr.compound_id = c.id
            ORDER BY rf.{vec_col} <=> :qv, rf.id
            LIMIT :k
        """).bindparams(
            bindparam("qv", type_=Vector(512)),
            bindparam("k")
        )

    elif metric == "tanimoto":
        # Exact Jaccard over BIT(512), same form as your jaccard_search_exact
        sql = text(f"""
            WITH top_hits AS (
                SELECT
                    rf.id,
                    cr.source,
                    cr.name,
                    c.smiles,
                    CASE
                    WHEN length(replace((rf.{bit_col} | :qb)::text, '0','')) = 0
                    THEN 1.0
                    ELSE length(replace((rf.{bit_col} & :qb)::text, '0',''))::float
                        / NULLIF(length(replace((rf.{bit_col} | :qb)::text, '0','')), 0)
                    END AS score
                FROM retrofingerprint rf
                JOIN retromol_compound rmc ON rmc.id = rf.retromol_compound_id
                JOIN compound c ON c.id = rmc.compound_id
                LEFT JOIN compound_record cr ON cr.compound_id = c.id
                WHERE rf.{bit_col} IS NOT NULL
                ORDER BY score DESC
                LIMIT :k
            )
            SELECT id, source, name, score, smiles
            FROM top_hits
            WHERE name IS NOT NULL
            ORDER BY score DESC, source, name, id
        """).bindparams(
            bindparam("qb", type_=BIT(512)),
            bindparam("k"),
        )
    else:
        raise ValueError(f"Unknown metric: {metric}")

    targets = parse_region_gbk_file(path, top_level=readout_toplevel)
    fp_labels = []
    fps = []
    for target in targets:
        kmers = []

        # Get tokenspects, currently only using it to mine for glycosylation
        tokenspecs = get_default_tokenspecs()
        mined_tokenspecs = mine_virtual_tokens(target, tokenspecs)
        found_glycosylation = any(ts["token"] == "glycosyltransferase" for ts in mined_tokenspecs)
        found_methylation = any(ts["token"] == "methyltransferase" for ts in mined_tokenspecs)
        if found_glycosylation:
            kmers.append((("glycosylation", None),))
        if found_methylation:
            kmers.append((("methylation", None),))

        # NOTE: readout for polyketides doesn't include structures
        fp_labels.append(f"{target.record_id}_{readout_toplevel}_{target.accession}")
        for readout in linear_readouts(
            target,
            cache_dir_override=cache_dir_override,
            level=readout_sublevel,
            pred_threshold=0.1
        ):
            seq = []
            for module in readout["readout"]:
                if isinstance(module, PKSModuleReadout):
                    seq.append((module.module_type.split("_")[1] + "0", None))
                elif isinstance(module, NRPSModuleReadout):
                    seq.append((module.substrate_name, module.substrate_smiles))
                else:
                    raise ValueError(f"Unsupported module type: {type(module)}")

            for k in kmer_sizes:
                kmers.extend(get_kmers(seq, k=k))

        fp = generator.fingerprint_from_kmers(kmers, num_bits=512, counted=counted, kmer_weights={1: 2, 2: 4, 3: 8})
        if fp is not None:
            fps.append(fp)

    # Stack fingerprints
    if len(fps) == 0:
        logger.warning("No RetroMol fingerprints could be generated for the input GenBank file.")
        return []
    else:
        fps = np.vstack(fps)

    # return results
    best_by_id: dict[tuple[int, str | None], dict] = {}
    with SessionLocal() as s:
        for fp_idx, fp in enumerate(fps):
            if metric == "cosine":
                qv = [float(x) for x in fp]
                rows = s.execute(sql, {"qv": qv, "k": top_k}).mappings().all()
            else:  # tanimoto: build 512-bit string (same idea as your exact Jaccard)
                qb = "".join("1" if float(x) > 0 else "0" for x in fp)  # length 512
                rows = s.execute(sql, {"qb": qb, "k": top_k}).mappings().all()

            for r in rows:
                key = (r["id"], r["source"])
                score = float(r["score"]) if r["score"] is not None else None
                prev = best_by_id.get(key)
                if prev is None or (score is not None and score > (prev["score"] or float("-inf"))):
                    best_by_id[key] = {
                        "record": fp_labels[fp_idx],
                        "id": r["id"],
                        "source": r["source"],
                        "name": r["name"],
                        "score": score,
                        "metric": metric,
                        "cosine": float(r["cosine"]) if metric == "cosine" and r.get("cosine") is not None else None,
                        "smiles": r["smiles"],
                    }

    results = sorted(best_by_id.values(),
                    key=lambda d: (d["score"] is not None, d["score"]),
                    reverse=True)[:top_k]
    return results
