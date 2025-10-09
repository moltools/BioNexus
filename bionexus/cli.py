from __future__ import annotations
from re import sub
import os, sys, argparse, subprocess
from rich.console import Console
from bionexus.config import DEFAULT_NPATLAS_URL
from bionexus.utils.logging import setup_logging
from bionexus.utils.net import ensure_download, extract_if_needed
from bionexus.db.models import Base
from bionexus.db.engine import engine
from sqlalchemy import text

console = Console()

setup_logging()

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

    url = args.url or DEFAULT_NPATLAS_URL
    if not url:
        console.print("[red]No NPAtlas URL set. Put NPATLAS_URL in .env or pass --url[/]")
        sys.exit(2)
    local_path = ensure_download(url, args.cache_dir, args.filename, args.force)
    extracted = extract_if_needed(local_path, args.cache_dir)

    def pick_data(files):  # choose a data file to load
        cands = [f for f in files if f.lower().endswith((".json"))]
        return cands[0] if cands else local_path
    
    data_path = pick_data(extracted) if extracted else local_path
    n = load_npatlas_file(path=data_path, chunk_size=args.chunk_size)
    console.print(f"Loaded {n} NPAtlas compounds")

def cmd_compute_fp(args: argparse.Namespace) -> None:
    from bionexus.etl.chemistry import backfill_fingerprints
    done = backfill_fingerprints(batch=args.batch)
    console.print(f"Computed fingerprints for {done} compounds (batch={args.batch})")

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
    p_np.add_argument("--filename", default=None, help="Save-as filename")
    p_np.add_argument("--force", action="store_true", help="Redownload even if present")
    p_np.add_argument("--chunk-size", type=int, default=10000)
    p_np.set_defaults(func=cmd_load_npatlas)

    p_fp = sub.add_parser("compute-fp", help="Compute fingerprints for compounds with SMILES")
    p_fp.add_argument("--batch", type=int, default=2000)
    p_fp.set_defaults(func=cmd_compute_fp)

    p_dump = sub.add_parser("dump-db", help="Write pg_dump custom format")
    p_dump.add_argument("--out", default="dumps/bionexus.dump")
    p_dump.set_defaults(func=cmd_dump_db)

    p_res = sub.add_parser("restore-db", help="Restore from pg_dump file")
    p_res.add_argument("--dump", required=True)
    p_res.set_defaults(func=cmd_restore_db)

    return p

def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)

if __name__ == "__main__":
    main()