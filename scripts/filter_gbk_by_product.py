#!/usr/bin/env python3
"""Filter GenBank files by antiSMASH cand_cluster/protocluster product, EXACT product match, and copy matches."""

from __future__ import annotations

import argparse
import os
import re
import shutil
from pathlib import Path
from typing import Sequence

from tqdm import tqdm

# Regex to capture all product qualifiers on a line
PRODUCT_RE = re.compile(r'/product="([^"]*)"', re.IGNORECASE)


def _iter_gbk_files(root: Path, exts: Sequence[str]):
    """Greedily walk the tree and yield matching files as we see them."""
    exts = tuple(ext.lower() for ext in exts)
    root = str(root)
    for dirpath, _dirnames, filenames in os.walk(root):
        for name in filenames:
            # simple suffix match, case-insensitive
            lower_name = name.lower()
            if lower_name.endswith(exts):
                yield Path(dirpath) / name


def _file_matches(gbk_path: Path, products: Sequence[str]) -> dict | None:
    """
    Return match info if the GBK contains cand_cluster/protocluster and has a product
    EXACTLY equal to one of the given search terms (case-insensitive).
    Stream + early exit.
    """
    has_cluster = False
    matched: set[str] = set()

    try:
        with gbk_path.open("r", errors="ignore") as fh:
            for line in fh:
                lower = line.lower()

                # detect cluster keywords
                if not has_cluster and ("cand_cluster" in lower or "protocluster" in lower):
                    has_cluster = True
                    if matched:
                        return {
                            "file": gbk_path.name,
                            "path": str(gbk_path),
                            "matched_products": sorted(matched),
                        }

                # detect product qualifiers
                if '/product="' in line:
                    for m in PRODUCT_RE.finditer(line):
                        prod_val = m.group(1).strip().lower()
                        if prod_val in products:
                            matched.add(prod_val)
                            if has_cluster:
                                return {
                                    "file": gbk_path.name,
                                    "path": str(gbk_path),
                                    "matched_products": sorted(matched),
                                }

        if has_cluster and matched:
            return {
                "file": gbk_path.name,
                "path": str(gbk_path),
                "matched_products": sorted(matched),
            }

        return None

    except Exception as exc:
        return {"error": f"read failed: {exc}", "file": str(gbk_path), "detail": str(exc)}


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Directory containing GenBank files")
    parser.add_argument("output_dir", type=Path, help="Directory to copy matching GBKs into")
    parser.add_argument(
        "-p",
        "--products",
        nargs="+",
        required=True,
        help="Exact product terms to match (case-insensitive)",
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

    products = [p.lower() for p in args.products]

    if not input_dir.exists():
        raise SystemExit(f"Input directory does not exist: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    errors = 0
    scanned = 0

    # No total: just chew through files as we discover them
    with tqdm(desc="Scanning GBKs", unit="file") as pbar:
        for path in _iter_gbk_files(input_dir, args.exts):
            scanned += 1
            result = _file_matches(path, products)

            if result:
                if "error" in result:
                    errors += 1
                else:
                    try:
                        shutil.copy2(result["path"], output_dir / result["file"])
                        copied += 1
                    except Exception:
                        errors += 1

            pbar.update(1)
            pbar.set_postfix(copied=copied, errors=errors)

    if scanned == 0:
        print("No GenBank files found with the given extensions.")
    else:
        print(f"\nDone. Copied: {copied}, errors: {errors}, total scanned: {scanned}")


if __name__ == "__main__":
    main()
