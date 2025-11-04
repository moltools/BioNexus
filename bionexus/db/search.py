"""Module for performing Jaccard similarity searches on chemical compound fingerprints."""

from __future__ import annotations

from sqlalchemy import text, bindparam
from sqlalchemy.dialects.postgresql import BIT
from pgvector.sqlalchemy import Vector

from bionexus.db.engine import SessionLocal


def jaccard_search_exact(bitstr: str, top_k: int = 50) -> list[dict]:
    """
    Perform an exact Jaccard similarity search on compound fingerprints.

    :param bitstr: the query fingerprint as a bit string
    :param top_k: the number of top results to return
    :return: a list of dictionaries containing compound information and Jaccard similarity scores
    """
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


def jaccard_search_hybrid(
    bitstr: str, qvec: list[float], top_k: int = 50, cand_k: int = 2000
) -> list[dict]:
    """
    Perform a hybrid Jaccard similarity search on compound fingerprints using an initial
    candidate selection based on vector similarity.

    :param bitstr: the query fingerprint as a bit string
    :param qvec: the query fingerprint as a vector of floats
    :param top_k: the number of top results to return
    :param cand_k: the number of candidate compounds to consider
    :return: a list of dictionaries containing compound information and Jaccard similarity scores
    """
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
        return (
            s.execute(sql, {"qb": bitstr, "qv": qvec, "k": top_k, "cand_k": cand_k})
            .mappings()
            .all()
        )
