"""Filter GenBank files by antiSMASH cand_cluster/protocluster product and copy matches."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Iterable, Sequence

from tqdm import tqdm


def _iter_gbk_files(root: Path, exts: Sequence[str]) -> Iterable[Path]:
    for ext in exts:
        yield from root.rglob(f"*{ext}")


def _file_matches(
    gbk_path: str,
    products: Sequence[str],
) -> dict | None:
    """
    Return match info if the GBK contains cand_cluster/protocluster with desired product.
    Uses lightweight text search to avoid full parsing overhead.
    """
    path = Path(gbk_path)
    try:
        text = path.read_text(errors="ignore")
    except Exception as exc:  # noqa: BLE001
        return {"error": f"read failed: {exc}", "file": gbk_path}

    lower = text.lower()
    if "cand_cluster" not in lower and "protocluster" not in lower:
        return None

    matched: set[str] = set()
    for prod in products:
        pat = re.compile(r'/product="[^"]*' + re.escape(prod.lower()) + r'[^"]*"', re.IGNORECASE)
        if pat.search(text):
            matched.add(prod)

    if not matched:
        return None

    return {
        "file": path.name,
        "path": str(path),
        "matched_products": sorted(matched),
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
        help='Product substrings to match (e.g., "NRPS", "T1PKS")',
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=os.cpu_count() or 1,
        help="Number of worker processes",
    )
    parser.add_argument(
        "--exts",
        nargs="+",
        default=[".gbk", ".gbff", ".gb"],
        help="File extensions to include",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=16,
        help="Chunksize for process pool mapping",
    )
    return parser.parse_args()


def main() -> None:
    args = cli()
    input_dir: Path = args.input_dir
    output_dir: Path = args.output_dir
    products = [p.lower() for p in args.products]

    if not input_dir.exists():
        raise SystemExit(f"Input directory does not exist: {input_dir}")

    files = list(_iter_gbk_files(input_dir, args.exts))
    if not files:
        raise SystemExit("No GenBank files found.")

    output_dir.mkdir(parents=True, exist_ok=True)

    total = len(files)
    copied = 0
    errors = 0

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        with tqdm(total=total, desc="Scanning GBKs", unit="file") as pbar:
            for result in executor.map(
                _file_matches,
                (str(p) for p in files),
                [products] * len(files),
                chunksize=args.chunksize,
            ):
                pbar.update(1)

                if not result:
                    continue

                if isinstance(result, dict) and "error" in result:
                    errors += 1
                    pbar.set_postfix(copied=copied, errors=errors)
                    continue

                src_path = Path(result["path"])
                dst_path = output_dir / src_path.name
                try:
                    shutil.copy2(src_path, dst_path)
                    copied += 1
                except Exception as exc:  # noqa: BLE001
                    errors += 1
                    result = {"error": f"copy failed: {exc}", "file": str(src_path)}

                pbar.set_postfix(copied=copied, errors=errors)


if __name__ == "__main__":
    main()
