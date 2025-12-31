#!/usr/bin/env python3
"""Find compounds that are biosynthetically similar but structurally different."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations, islice
from typing import Iterable

from sqlalchemy import select
from tqdm import tqdm

from bionexus.db.engine import SessionLocal
from bionexus.db.models import (
    Annotation,
    Compound,
    CompoundRecord,
    RetroFingerprint,
    RetroMolCompound,
)


@dataclass(slots=True)
class CompoundFP:
    """Container for the fingerprints we need to compare."""

    cid: int
    name: str | None
    source: str | None
    smiles: str | None
    morgan_bits: int
    morgan_pop: int
    retro_bits: list[tuple[int, int]]

    def label(self) -> str:
        if self.name and self.source:
            return f"{self.name} [{self.source}]"
        if self.name:
            return self.name
        return f"compound_{self.cid}"


def _bits_to_int(val: object) -> int | None:
    """Convert a BIT column payload into an int for fast popcount operations."""
    if val is None:
        return None
    if isinstance(val, str):
        val = val.strip()
        return int(val, 2) if val else None
    if isinstance(val, memoryview):
        val = val.tobytes()
    if isinstance(val, (bytes, bytearray)):
        return int.from_bytes(val, "big") if val else None
    if isinstance(val, int):
        return val
    raise TypeError(f"Unsupported BIT payload type: {type(val)!r}")


def _tanimoto(a_bits: int, a_pop: int, b_bits: int, b_pop: int) -> float:
    """Fast Tanimoto over integer bitmasks."""
    inter = (a_bits & b_bits).bit_count()
    union = a_pop + b_pop - inter
    return inter / union if union else 0.0


def _best_retro_similarity(
    lhs: list[tuple[int, int]],
    rhs: list[tuple[int, int]],
    early_stop: float,
) -> float:
    """Return the max Tanimoto across two sets of biosynthetic fingerprints."""
    best = 0.0
    for a_bits, a_pop in lhs:
        for b_bits, b_pop in rhs:
            score = _tanimoto(a_bits, a_pop, b_bits, b_pop)
            if score > best:
                best = score
                if best > early_stop:
                    return best
    return best


def _load_compounds(
    min_coverage: float,
    limit: int | None,
    np_classes: list[str] | None,
) -> list[CompoundFP]:
    """
    Load compounds with both Morgan and RetroMol fingerprints.

    Only RetroMol results with coverage >= min_coverage are considered.
    Optionally restrict to compounds carrying an npclassifier class annotation.
    """
    if np_classes:
        np_classes = [c.strip() for c in np_classes if c and c.strip()]
        if not np_classes:
            np_classes = None

    stmt = (
        select(
            Compound.id.label("cid"),
            Compound.smiles,
            Compound.fp_morgan_b2048_r2_bit.label("morgan_bit"),
            Compound.fp_morgan_b2048_r2_pop.label("morgan_pop"),
            RetroFingerprint.fp_retro_b512_bit.label("retro_bit"),
            RetroFingerprint.fp_retro_b512_pop.label("retro_pop"),
            CompoundRecord.name,
            CompoundRecord.source,
        )
        .join(RetroMolCompound, RetroMolCompound.compound_id == Compound.id)
        .join(RetroFingerprint, RetroFingerprint.retromol_compound_id == RetroMolCompound.id)
        .outerjoin(CompoundRecord, CompoundRecord.compound_id == Compound.id)
        .where(
            RetroMolCompound.coverage.is_not(None),
            RetroMolCompound.coverage >= min_coverage,
            Compound.fp_morgan_b2048_r2_bit.is_not(None),
            RetroFingerprint.fp_retro_b512_bit.is_not(None),
        )
        .order_by(Compound.id)
    )

    if np_classes:
        stmt = stmt.join(
            Annotation,
            (
                (Annotation.compound_id == Compound.id)
                & (Annotation.scheme == "npclassifier")
                & (Annotation.key == "class")
                & (Annotation.value.in_(np_classes))
            ),
        )

    compounds: dict[int, CompoundFP] = {}
    with SessionLocal() as session:
        rows: Iterable[dict] = session.execute(stmt).mappings()
        if limit:
            rows = islice(rows, limit)

        for row in rows:
            cid = row["cid"]
            morgan_int = _bits_to_int(row["morgan_bit"])
            retro_int = _bits_to_int(row["retro_bit"])

            if morgan_int is None or retro_int is None:
                continue

            comp = compounds.get(cid)
            if comp is None:
                morgan_pop = int(row["morgan_pop"]) if row["morgan_pop"] is not None else morgan_int.bit_count()
                comp = CompoundFP(
                    cid=cid,
                    name=row.get("name"),
                    source=row.get("source"),
                    smiles=row.get("smiles"),
                    morgan_bits=morgan_int,
                    morgan_pop=morgan_pop,
                    retro_bits=[],
                )
                compounds[cid] = comp
            elif not comp.name and row.get("name"):
                comp.name = row["name"]
                comp.source = row.get("source")

            retro_pop = int(row["retro_pop"]) if row["retro_pop"] is not None else retro_int.bit_count()
            comp.retro_bits.append((retro_int, retro_pop))

    return list(compounds.values())


def _find_congener_pairs(
    compounds: list[CompoundFP],
    max_struct_sim: float,
    min_biosyn_sim: float,
) -> list[tuple[float, float, CompoundFP, CompoundFP]]:
    """Pairwise scan for biosynthetically similar yet structurally distant compounds."""
    results: list[tuple[float, float, CompoundFP, CompoundFP]] = []
    for lhs, rhs in tqdm(combinations(compounds, 2)):
        struct_sim = _tanimoto(lhs.morgan_bits, lhs.morgan_pop, rhs.morgan_bits, rhs.morgan_pop)
        if struct_sim >= max_struct_sim:
            continue

        biosyn_sim = _best_retro_similarity(lhs.retro_bits, rhs.retro_bits, min_biosyn_sim)
        if biosyn_sim > min_biosyn_sim:
            results.append((biosyn_sim, struct_sim, lhs, rhs))

    results.sort(key=lambda item: (item[0], -item[1]), reverse=True)
    return results


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-coverage", type=float, default=0.9, help="Minimum RetroMol coverage to include.")
    parser.add_argument("--max-struct-sim", type=float, default=0.675, help="Upper bound on structural Tanimoto.")
    parser.add_argument("--min-biosyn-sim", type=float, default=0.6, help="Lower bound on biosynthetic Tanimoto.")
    parser.add_argument(
        "--npclassifier-class",
        action="append",
        dest="np_classes",
        help="Restrict to compounds annotated with npclassifier class (repeatable).",
    )
    parser.add_argument(
        "--max-compounds",
        type=int,
        default=None,
        help="Optional cap on compounds loaded (first N by ID) for quick scans.",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=50,
        help="Only print the top N pairs by biosynthetic similarity.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    compounds = _load_compounds(args.min_coverage, args.max_compounds, args.np_classes)

    if not compounds:
        raise SystemExit("No compounds found with the requested coverage and fingerprints.")

    pairs = _find_congener_pairs(compounds, args.max_struct_sim, args.min_biosyn_sim)
    if not pairs:
        print("No candidate congeners found with the given thresholds.")
        return

    print(
        f"Loaded {len(compounds)} compounds with coverage >= {args.min_coverage:.2f}; "
        f"found {len(pairs)} candidate pairs."
    )
    for idx, (bio, struct, lhs, rhs) in enumerate(islice(pairs, args.max_results), start=1):
        print(
            f"{idx:>3}: {lhs.label()} (id={lhs.cid}) vs {rhs.label()} (id={rhs.cid}) | "
            f"biosyn={bio:.3f}, struct={struct:.3f}"
        )


if __name__ == "__main__":
    main()
