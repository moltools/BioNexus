from __future__ import annotations
import os, sys, argparse, subprocess, logging, json
import pandas as pd
from pathlib import Path
from rich.table import Table
from rich.console import Console
from bionexus.config import DEFAULT_NPATLAS_URL, DEFAULT_MIBIG_JSON_URL
from bionexus.utils.logging import setup_logging
from bionexus.utils.net import ensure_download, extract_if_needed
from bionexus.db.models import Base
from bionexus.db.engine import engine
from sqlalchemy import text

console = Console()

setup_logging()
logger = logging.getLogger(__name__)

# ---------- commands ----------

def cmd_init_db(args: argparse.Namespace) -> None:
    if args.drop:
        with engine.begin() as conn:
            conn.execute(text("DROP SCHEMA public CASCADE; CREATE SCHEMA public;"))
            conn.execute(text("CREATE EXTENSION IF NOT EXISTS pgvector;"))
    Base.metadata.create_all(engine)
    console.print("[green]Schema ready.[/]")

def cmd_upgrade(args: argparse.Namespace) -> None:
    rev = args.rev or "head"
    rc = subprocess.call([sys.executable, "-m", "alembic", "upgrade", rev])
    sys.exit(rc)

def cmd_load_npatlas(args: argparse.Namespace) -> None:
    from bionexus.etl.sources.npatlas import load_npatlas_file

    # download and extract if needed
    url = args.url or DEFAULT_NPATLAS_URL
    if not url:
        console.print("[red]No NPAtlas URL set. Put NPATLAS_URL in .env or pass --url[/]")
        sys.exit(2)
    local_path = ensure_download(url, args.cache_dir, args.force)
    extracted = extract_if_needed(local_path, args.cache_dir)

    def pick_data(files):  # choose a data file to load
        cands = [f for f in files if f.lower().endswith((".json"))]
        return cands[0] if cands else local_path
    
    data_path = pick_data(extracted) if extracted else local_path
    n_compounds, n_records, n_annotations = load_npatlas_file(path=data_path, chunk_size=args.chunk_size)
    console.print(f"Loaded {n_compounds} NPAtlas compounds, {n_records} associated compound records, and {n_annotations} associated annotations")

def cmd_load_mibig(args: argparse.Namespace) -> None:
    from bionexus.etl.sources.mibig import load_mibig_files, download_mibig_gbks

    # download and extract if needed
    url_json = args.url_json or DEFAULT_MIBIG_JSON_URL
    if not url_json:
        console.print("[red]No MIBiG URL set. Put MIBIG_JSON_URL and MIBIG_GBK_URL in .env or pass --url-json and --url-gbk[/]")
        sys.exit(2)
    local_path_json = ensure_download(url_json, args.cache_dir, args.force)
    extracted_jsons = extract_if_needed(local_path_json, args.cache_dir)

    # quickly loop over JSONS and extract accessions with their versions so we can download up-to-date GBKs with product info
    # gbk_url we are using: https://mibig.secondarymetabolites.org/repository/{accession}.{version}/generated/{accession}.gbk
    versioned_accessions: tuple[str, str] = []
    for json_path in (extracted_jsons if extracted_jsons else [local_path_json]):
        with open(json_path) as f:
            data = json.load(f)
            accession = data.get("accession", None)
            version = data.get("version", None)
            if accession and version:
                versioned_accessions.append((accession, version))
    console.print(f"Found {len(versioned_accessions)} MIBiG JSON files with accessions and versions")      
    gbk_paths: list[Path] = download_mibig_gbks(accessions=versioned_accessions, outdir=os.path.join(args.cache_dir, "mibig_gbks"))
    console.print(f"Downloaded {len(gbk_paths)} MIBiG GenBank files")

    n_compounds, n_records, n_regions, n_annotations = load_mibig_files(
        json_paths=extracted_jsons,
        gbk_paths=[str(p) for p in gbk_paths],
        chunk_size=args.chunk_size
    )
    console.print(f"Loaded {n_regions} MIBiG regions, {n_compounds} associated compound structures, {n_records} associated compound records, and {n_annotations} associated annotations")

def cmd_annot_npc(args: argparse.Namespace) -> None:
    from bionexus.etl.sources.npclassifier import annotate_with_npclassifier
    n_annotated = annotate_with_npclassifier(args.recompute, chunk_size=args.batch)
    console.print(f"Used NPClassifier to add {n_annotated} annotations for compounds")

def cmd_compute_fp_morgan(args: argparse.Namespace) -> None:
    from bionexus.etl.chemistry import backfill_fingerprints
    done = backfill_fingerprints(batch=args.batch, recompute=args.recompute)
    console.print(f"Computed fingerprints for {done} compounds (batch={args.batch})")

def cmd_compute_fp_retro(args: argparse.Namespace) -> None:
    console.print("[yellow]Biosynthetic fingerprint computation is not yet implemented[/]")
    exit(1)

def cmd_dump_db(args: argparse.Namespace) -> None:
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    url = os.getenv("BIONEXUS_DB_URL")
    if not url:
        console.print("[red]BIONEXUS_DB_URL is not set[/]")
        sys.exit(2)
    url_plain = url.replace("+psycopg", "")
    cmd = [
        "pg_dump", "--format=custom", "--no-owner", "--no-privileges",
        f"--dbname={url_plain}", "-f", args.out
    ]
    sys.exit(subprocess.call(cmd))

def cmd_restore_db(args: argparse.Namespace) -> None:
    url = os.getenv("BIONEXUS_DB_URL")
    if not url:
        console.print("[red]BIONEXUS_DB_URL is not set[/]")
        sys.exit(2)
    url_plain = url.replace("+psycopg", "")
    cmd = [
        "pg_restore", "--clean", "--if-exists", "--no-owner", "--no-privileges",
        f"--dbname={url_plain}", args.dump
    ]
    sys.exit(subprocess.call(cmd))

def cmd_search_morgan(args):
    from bionexus.etl.chemistry import _morgan_bits_and_vec
    from bionexus.db.search import jaccard_search_exact, jaccard_search_hybrid

    bits, _pop, vec = _morgan_bits_and_vec(args.smiles)
    if not bits:
        logger.warning("Could not compute fingerprint for SMILES")  
        return

    rows = jaccard_search_hybrid(bits, vec, args.top_k) \
           if getattr(args, "hybrid", False) \
           else jaccard_search_exact(bits, args.top_k)

    df = pd.DataFrame(rows)
    if not args.out:
        # console.print(df[["source", "name", "jacc"]])
        table = Table(show_header=True, header_style="bold magenta")
        table.add_column("id", style="dim", width=6)
        table.add_column("source", style="dim", width=10)
        table.add_column("name", style="white", width=30)
        table.add_column("jacc", justify="right")
        for _, row in df.iterrows():
            table.add_row(f"{row['id']}", str(row["source"]), f"[cyan]{row['name']}", f"{row['jacc']:.3f}")
        console.print(table)
    else:
        sep = "," if args.out.lower().endswith(".csv") else "\t"
        df.to_csv(args.out, index=False, sep=sep)
        console.print(f"Wrote {len(df)} results to [green]{args.out}[/]")

def cmd_search_retro(args):
    console.print("[yellow]Searching by biosynthetic fingerprint is not yet implemented[/]")
    exit(1)

# ---------- parser ----------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="bionexus", description="Bionexus CLI (argparse)")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_init = sub.add_parser("init-db", help="Create schema (optionally drop first)")
    p_init.add_argument("--drop", action="store_true")
    p_init.set_defaults(func=cmd_init_db)

    p_up = sub.add_parser("upgrade", help="Run Alembic upgrade")
    p_up.add_argument("rev", nargs="?", default="head")
    p_up.set_defaults(func=cmd_upgrade)

    p_np = sub.add_parser("load-npatlas", help="Download (if needed) and load NPAtlas")
    p_np.add_argument("--url", default=None, help="Override NPAtlas download URL")
    p_np.add_argument("--cache-dir", default=None, help="Cache/work dir")
    p_np.add_argument("--force", action="store_true", help="Redownload even if present")
    p_np.add_argument("--chunk-size", type=int, default=10000)
    p_np.set_defaults(func=cmd_load_npatlas)

    p_mibig = sub.add_parser("load-mibig", help="Download (if needed) and load MIBiG")
    p_mibig.add_argument("--url-json", default=None, help="Override MIBiG JSON download URL")
    p_mibig.add_argument("--url-gbk", default=None, help="Override MIBiG GenBank download URL")
    p_mibig.add_argument("--cache-dir", default=None, help="Cache/work dir")
    p_mibig.add_argument("--force", action="store_true", help="Redownload even if present")
    p_mibig.add_argument("--chunk-size", type=int, default=1000)
    p_mibig.set_defaults(func=cmd_load_mibig)

    p_annot_npc = sub.add_parser("annotate-npc", help="Annotate NPClassifier classes for compounds without one")
    p_annot_npc.add_argument("--batch", type=int, default=2000)
    p_annot_npc.add_argument("--recompute", action="store_true", help="Force recomputation for all compounds")
    p_annot_npc.set_defaults(func=cmd_annot_npc)

    p_fp = sub.add_parser("compute-fp-morgan", help="Compute fingerprints for compounds with SMILES")
    p_fp.add_argument("--batch", type=int, default=2000)
    p_fp.add_argument("--recompute", action="store_true", help="Force recomputation for all compounds")
    p_fp.set_defaults(func=cmd_compute_fp_morgan)

    p_do = sub.add_parser("compute-fp-retro", help="Compute biosynthetic fingerprints for compounds and/or GenBank records")
    p_do.add_argument("--for", choices=["compounds", "gbks", "both"], default="both", help="What to compute fingerprints for")
    p_do.add_argument("--batch", type=int, default=2000)
    p_do.add_argument("--recompute", action="store_true", help="Force recomputation for all compounds/records")
    p_do.set_defaults(func=cmd_compute_fp_retro)

    p_dump = sub.add_parser("dump-db", help="Write pg_dump custom format")
    p_dump.add_argument("--out", default="dumps/bionexus.dump")
    p_dump.set_defaults(func=cmd_dump_db)

    p_res = sub.add_parser("restore-db", help="Restore from pg_dump file")
    p_res.add_argument("--dump", required=True)
    p_res.set_defaults(func=cmd_restore_db)

    p_search_m = sub.add_parser("search-morgan", help="Search compounds by Jaccard to a SMILES")
    p_search_m.add_argument("--smiles", required=True)
    p_search_m.add_argument("--top-k", type=int, default=20)
    p_search_m.add_argument("--hybrid", action="store_true", help="Use pgvector candidate search")
    p_search_m.add_argument("--out", default=None, help="Optional output file (TSV/CSV)")
    p_search_m.set_defaults(func=cmd_search_morgan)

    p_search_r = sub.add_parser("search-retro", help="Search compoundsa and BGCs to a biosynthetic fingerprint")
    p_search_r.add_argument("--for", required=True, choices=["compound", "gbk"], help="Input type")
    p_search_r.add_argument("--input", required=True, type=str, help="SMILES for 'compound' input type or path to GBK region file for 'gbk' input type")
    p_search_r.add_argument("--out", default=None, help="Optional output file (TSV/CSV)")
    p_search_r.set_defaults(func=cmd_search_retro)

    return p

def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)

if __name__ == "__main__":
    main()