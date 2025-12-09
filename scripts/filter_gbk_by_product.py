#!/usr/bin/env python3
"""Filter GenBank files by antiSMASH cand_cluster/protocluster product, EXACT product match, and copy matches."""

from __future__ import annotations

import argparse
import os
import re
import shutil
from pathlib import Path
from typing import Iterable, Sequence

from tqdm import tqdm


# Regex to capture all product qualifiers
PRODUCT_RE = re.compile(r'/product="([^"]*)"', re.IGNORECASE)


def _iter_gbk_files(root: Path, exts: Sequence[str]) -> Iterable[Path]:
    """Yield all files with specified extensions."""
    exts = {e.lower() for e in exts}
    for p in root.rglob("*"):
        if p.is_file() and p.suffix.lower() in exts:
            yield p


def _file_matches(
    gbk_path: str,
    products: Sequence[str],
) -> dict | None:
    """
    Return match info if the GBK contains cand_cluster/protocluster and has a product EXACTLY equal
    to one of the given search terms (case-insensitive).
    """
    path = Path(gbk_path)

    try:
        text = path.read_text(errors="ignore")
    except Exception as exc:
        return {"error": f"read failed: {exc}", "file": gbk_path}

    lower = text.lower()

    # quick skip
    if "cand_cluster" not in lower and "protocluster" not in lower:
        return None

    # extract all /product="..."`
    found_products = PRODUCT_RE.findall(text)
    found_products_norm = {fp.strip().lower() for fp in found_products}

    # find exact matches
    matched = sorted(fp for fp in found_products_norm if fp in products)

    if not matched:
        return None

    return {
        "file": path.name,
        "path": str(path),
        "matched_products": matched,
    }


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Directory containing GenBank files")
    parser.add_argument("output_dir", type=Path, help="Directory to copy matching GBKs into")
    parser.add_argument(
        "-p",
        "--products",
        nargs="+",
        required=True,
        help='Exact product terms to match (case-insensitive)',
    )
    parser.add_argument(
        "--exts",
        nargs="+",
        default=[".gbk", ".gbff", ".gb"],
        help="File extensions to include",
    )
    return parser.parse_args()


def main() -> None:
    args = cli()
    input_dir: Path = args.input_dir
    output_dir: Path = args.output_dir

    # normalize to lowercase for exact matching
    products = [p.lower() for p in args.products]

    if not input_dir.exists():
        raise SystemExit(f"Input directory does not exist: {input_dir}")

    files = list(_iter_gbk_files(input_dir, args.exts))
    if not files:
        raise SystemExit("No GenBank files found.")

    output_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    errors = 0

    with tqdm(total=len(files), desc="Scanning GBKs", unit="file") as pbar:
        for path in files:
            result = _file_matches(str(path), products)
            pbar.update(1)

            if not result:
                continue

            if "error" in result:
                errors += 1
                pbar.set_postfix(copied=copied, errors=errors)
                continue

            try:
                shutil.copy2(result["path"], output_dir / result["file"])
                copied += 1
            except Exception as exc:
                errors += 1

            pbar.set_postfix(copied=copied, errors=errors)

    print(f"\nDone. Copied: {copied}, errors: {errors}, total scanned: {len(files)}")


if __name__ == "__main__":
    main()
