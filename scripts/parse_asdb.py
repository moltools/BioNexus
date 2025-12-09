"""Parse a directory of GenBank files with BioCracker and dump readouts to JSONL."""

from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Iterable, Sequence

from tqdm import tqdm

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))


def _iter_gbk_files(root: Path, exts: Sequence[str]) -> Iterable[Path]:
    """Yield GenBank files under root matching extensions."""
    for ext in exts:
        yield from root.rglob(f"*{ext}")


def _process_file(
    gbk_path: str,
    cache_dir: str | None,
    pred_threshold: float,
) -> dict | None:
    """Worker: parse one GBK file and return JSON-serializable payload."""
    from bionexus.etl.sources.biocracker import to_jsonable

    try:
        from biocracker.antismash import parse_region_gbk_file
        from biocracker.readout import linear_readouts as biocracker_linear_readouts
    except ImportError as exc:
        return {"error": f"biocracker import failed: {exc}", "file": gbk_path}

    path = Path(gbk_path)
    try:
        targets = parse_region_gbk_file(path, top_level="cand_cluster")
    except Exception as exc:  # noqa: BLE001
        return {"error": f"parse failed: {exc}", "file": gbk_path}

    all_targets: list[list[dict]] = []
    accessions: set[str] = set()

    for target in targets:
        readouts_for_target: list[dict] = []
        try:
            readouts = biocracker_linear_readouts(
                target,
                model=None,
                cache_dir_override=cache_dir,
                level="gene",
                pred_threshold=pred_threshold,
            )
        except Exception as exc:  # noqa: BLE001
            return {"error": f"readout failed: {exc}", "file": gbk_path}

        for readout in readouts:
            try:
                readout_dict = to_jsonable(readout)
            except Exception as exc:  # noqa: BLE001
                return {"error": f"serialization failed: {exc}", "file": gbk_path}

            rec = readout_dict.get("rec", {})
            acc = rec.get("accession") or rec.get("acc") or rec.get("id")
            if acc:
                accessions.add(acc)

            readouts_for_target.append(readout_dict)

        if readouts_for_target:
            all_targets.append(readouts_for_target)

    if not all_targets:
        return None

    return {
        "file": path.name,
        "path": str(path),
        "accessions": sorted(accessions),
        "targets": all_targets,
    }


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Directory containing GenBank files")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("asdb_readouts.jsonl"),
        help="Output JSONL file (one record per GenBank file)",
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=os.cpu_count() or 1,
        help="Number of worker processes",
    )
    parser.add_argument(
        "--pred-threshold",
        type=float,
        default=0.1,
        help="Prediction threshold passed to BioCracker linear readouts",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Optional cache dir override for BioCracker",
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
        default=8,
        help="Chunksize for process pool mapping",
    )
    return parser.parse_args()


def main() -> None:
    args = cli()
    input_dir: Path = args.input_dir
    if not input_dir.exists():
        raise SystemExit(f"Input directory does not exist: {input_dir}")

    files = list(_iter_gbk_files(input_dir, args.exts))
    if not files:
        raise SystemExit("No GenBank files found.")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    total = len(files)
    processed = 0
    successes = 0
    errors = 0

    with args.output.open("w", encoding="utf-8") as outf:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            with tqdm(total=total, desc="Parsing GBKs", unit="file") as pbar:
                for result in executor.map(
                    _process_file,
                    (str(p) for p in files),
                    [str(args.cache_dir) if args.cache_dir else None] * len(files),
                    [args.pred_threshold] * len(files),
                    chunksize=args.chunksize,
                ):
                    processed += 1
                    pbar.update(1)

                    if not result:
                        continue

                    if isinstance(result, dict) and "error" in result:
                        errors += 1
                    else:
                        successes += 1

                    outf.write(json.dumps(result, ensure_ascii=False) + "\n")
                    pbar.set_postfix(ok=successes, errors=errors)


if __name__ == "__main__":
    main()
