"""BioNexus command-line interface (CLI) using argparse."""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from rich.console import Console
from rich.table import Table
from sqlalchemy import text

from bionexus.config import DEFAULT_MIBIG_JSON_URL, DEFAULT_NPATLAS_URL
from bionexus.db.engine import engine
from bionexus.db.models import Base
from bionexus.utils.logging import setup_logging
from bionexus.utils.net import ensure_download, extract_if_needed

console = Console()


setup_logging()
logger = logging.getLogger(__name__)


def _run_alembic(*args: str) -> int:
    """
    Run an Alembic command via subprocess.

    :param args: arguments to pass to Alembic
    :return: return code from the Alembic command
    """
    cmd = [sys.executable, "-m", "alembic", *args]
    logger.debug(f"Running Alembic command: {' '.join(cmd)}")
    return subprocess.call(cmd)


def cmd_init_db(args: argparse.Namespace) -> None:
    """
    Initialize the database schema, optionally dropping existing schema first.

    :param args: command-line arguments
    """
    if args.drop:
        with engine.begin() as conn:
            conn.execute(text("DROP SCHEMA public CASCADE; CREATE SCHEMA public;"))
            conn.execute(text("CREATE EXTENSION IF NOT EXISTS pgvector;"))
    Base.metadata.create_all(engine)
    console.print("[green]Schema ready.[/]")


def cmd_upgrade(args: argparse.Namespace) -> None:
    """
    Upgrade the database schema to a specified revision.

    :param args: command-line arguments
    """
    rev = args.rev or "head"
    rc = _run_alembic("upgrade", rev)
    sys.exit(rc)


def cmd_downgrade(args: argparse.Namespace) -> None:
    """
    Downgrade the database schema to a specified revision.

    :param args: command-line arguments
    """
    rev = args.rev
    if not rev:
        console.print("[red]You must provide a revision to downgrade to (e.g., -1, base, or a rev id)[/]")
        sys.exit(2)
    sys.exit(_run_alembic("downgrade", rev))


def cmd_revision(args: argparse.Namespace) -> None:
    """
    Create a new Alembic revision with autogeneration.

    :param args: command-line arguments
    """
    if not args.message:
        console.print("[red]You must provide a message for the revision[/]")
        sys.exit(2)

    cmd = ["revision", "--autogenerate", "-m", args.message]
    if args.head:
        cmd += ["--head", args.head]
    if args.splice:
        cmd += ["--splice"]
    if args.rev_id:
        cmd += ["--rev-id", args.rev_id]
    for dep in args.depends_on or []:
        cmd += ["--depends-on", dep]
    for label in args.branch_label or []:
        cmd += ["--branch-label", label]

    sys.exit(_run_alembic(*cmd))


def cmd_history(args: argparse.Namespace) -> None:
    """
    Show Alembic migration history.

    :param args: command-line arguments
    """
    sys.exit(_run_alembic("history", "--verbose"))


def cmd_current(args: argparse.Namespace) -> None:
    """
    Show the current Alembic revision.

    :param args: command-line arguments
    """
    sys.exit(_run_alembic("current", "--verbose"))


def cmd_heads(args: argparse.Namespace) -> None:
    """
    Show the current Alembic heads.

    :param args: command-line arguments
    """
    sys.exit(_run_alembic("heads", "--verbose"))


def cmd_stamp(args: argparse.Namespace) -> None:
    """
    Stamp the database with a specified Alembic revision without running migrations.

    :param args: command-line arguments
    """
    if not args.rev:
        console.print("[red]Provide a revision to stamp (e.g. head or a specific rev id)[/]")
        sys.exit(2)
    sys.exit(_run_alembic("stamp", args.rev))


def cmd_load_npatlas(args: argparse.Namespace) -> None:
    """
    Download (if needed) and load NPAtlas data into the database.

    :param args: command-line arguments
    """
    from bionexus.etl.sources.npatlas import load_npatlas_file

    # Download and extract if needed
    url = args.url or DEFAULT_NPATLAS_URL
    if not url:
        console.print("[red]No NPAtlas URL set. Put NPATLAS_URL in .env or pass --url[/]")
        sys.exit(2)
    local_path = ensure_download(url, args.cache_dir, args.force)
    extracted = extract_if_needed(local_path, args.cache_dir)

    def pick_data(files):  # choose a data file to load
        cands = [f for f in files if f.lower().endswith(".json")]
        return cands[0] if cands else local_path

    data_path = pick_data(extracted) if extracted else local_path
    n_compounds, n_records, n_annotations = load_npatlas_file(path=data_path, chunk_size=args.chunk_size)
    console.print(
        f"Loaded {n_compounds} NPAtlas compounds, "
        f"{n_records} associated compound records, "
        f"and {n_annotations} associated annotations"
    )


def cmd_load_mibig(args: argparse.Namespace) -> None:
    """
    Download (if needed) and load MIBiG data into the database.

    :param args: command-line arguments
    """
    from bionexus.etl.sources.mibig import download_mibig_gbks, load_mibig_files

    # Download and extract if needed
    url_json = args.url_json or DEFAULT_MIBIG_JSON_URL
    if not url_json:
        console.print(
            "[red]No MIBiG URL set. Put MIBIG_JSON_URL and MIBIG_GBK_URL in .env or pass --url-json and --url-gbk[/]"
        )
        sys.exit(2)
    local_path_json = ensure_download(url_json, args.cache_dir, args.force)
    extracted_jsons = extract_if_needed(local_path_json, args.cache_dir)

    # Loop over JSONS and extract accessions with their versions so we can download up-to-date GBKs with product info
    # gbk_url we are using:
    #   https://mibig.secondarymetabolites.org/repository/{accession}.{version}/generated/{accession}.gbk
    versioned_accessions: tuple[str, str] = []
    for json_path in extracted_jsons if extracted_jsons else [local_path_json]:
        with open(json_path) as f:
            data = json.load(f)
            accession = data.get("accession", None)
            version = data.get("version", None)
            if accession and version:
                versioned_accessions.append((accession, version))
    console.print(f"Found {len(versioned_accessions)} MIBiG JSON files with accessions and versions")
    gbk_paths: list[Path] = download_mibig_gbks(
        accessions=versioned_accessions,
        outdir=os.path.join(args.cache_dir, "mibig_gbks"),
    )
    console.print(f"Downloaded {len(gbk_paths)} MIBiG GenBank files")

    n_compounds, n_records, n_regions, n_annotations = load_mibig_files(
        json_paths=extracted_jsons,
        gbk_paths=[str(p) for p in gbk_paths],
        chunk_size=args.chunk_size,
    )
    console.print(
        f"Loaded {n_regions} MIBiG regions, "
        f"{n_compounds} associated compound structures, "
        f"{n_records} associated compound records, "
        f"and {n_annotations} associated annotations"
    )


def cmd_annot_npc(args: argparse.Namespace) -> None:
    """
    Annotate compounds without NPClassifier classes.

    :param args: command-line arguments
    """
    from bionexus.etl.sources.npclassifier import annotate_with_npclassifier

    n_annotated = annotate_with_npclassifier(args.recompute, chunk_size=args.batch)
    console.print(f"Used NPClassifier to add {n_annotated} annotations for compounds")


def cmd_parse_compounds(args: argparse.Namespace) -> None:
    """
    Parse compounds from the database using RetroMol.

    :param args: command-line arguments
    """
    from bionexus.etl.sources.retromol import parse_compounds_with_retromol

    cache_dir = Path(args.cache_dir) if args.cache_dir else Path("./retromol_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    n_parsed = parse_compounds_with_retromol(
        cache_dir=cache_dir,
        recompute=args.recompute,
        chunk_size=args.batch,
        workers=args.workers,
    )
    console.print(f"Used RetroMol to parse {n_parsed} compounds")


def cmd_compute_fp_morgan(args: argparse.Namespace) -> None:
    """
    Compute Morgan fingerprints for compounds with SMILES.

    :param args: command-line arguments
    """
    from bionexus.etl.chemistry import backfill_fingerprints

    done = backfill_fingerprints(batch=args.batch, recompute=args.recompute)
    console.print(f"Computed fingerprints for {done} compounds (batch={args.batch})")


def cmd_compute_fp_retro_compound(args: argparse.Namespace) -> None:
    """
    Compute biosynthetic fingerprints for compounds and/or GenBank records.

    :param args: command-line arguments
    """
    from bionexus.etl.retromol import backfill_retro_fingerprints

    done = backfill_retro_fingerprints(cache_dir=args.cache_dir, batch=args.batch, recompute=args.recompute)
    console.print(f"Computed RetroMol fingerprints for {done} entries (batch={args.batch})")


def cmd_dump_db(args: argparse.Namespace) -> None:
    """
    Dump the database to a pg_dump custom format file.

    :param args: command-line arguments
    """
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    url = os.getenv("BIONEXUS_DB_URL")
    if not url:
        console.print("[red]BIONEXUS_DB_URL is not set[/]")
        sys.exit(2)
    url_plain = url.replace("+psycopg", "")
    cmd = [
        "pg_dump",
        "--format=custom",
        "--no-owner",
        "--no-privileges",
        f"--dbname={url_plain}",
        "-f",
        args.out,
    ]
    sys.exit(subprocess.call(cmd))


def cmd_restore_db(args: argparse.Namespace) -> None:
    """
    Restore the database from a pg_dump custom format file.

    :param args: command-line arguments
    """
    url = os.getenv("BIONEXUS_DB_URL")
    if not url:
        console.print("[red]BIONEXUS_DB_URL is not set[/]")
        sys.exit(2)
    url_plain = url.replace("+psycopg", "")
    cmd = [
        "pg_restore",
        "--clean",
        "--if-exists",
        "--no-owner",
        "--no-privileges",
        f"--dbname={url_plain}",
        args.dump,
    ]
    sys.exit(subprocess.call(cmd))


def cmd_search_morgan(args: argparse.Namespace) -> None:
    """
    Search compounds by Jaccard similarity to a given SMILES using Morgan fingerprints.

    :param args: command-line arguments
    """
    from bionexus.db.search import jaccard_search_exact, jaccard_search_hybrid
    from bionexus.etl.chemistry import _morgan_bits_and_vec

    bits, _pop, vec = _morgan_bits_and_vec(args.smiles)
    if not bits:
        logger.warning("Could not compute fingerprint for SMILES")
        return

    rows = (
        jaccard_search_hybrid(bits, vec, args.top_k)
        if getattr(args, "hybrid", False)
        else jaccard_search_exact(bits, args.top_k)
    )

    df = pd.DataFrame(rows)
    if not args.out:
        # console.print(df[["source", "name", "jacc"]])
        table = Table(show_header=True, header_style="bold magenta")
        table.add_column("compound_id", style="dim", width=6)
        table.add_column("source", style="dim", width=10)
        table.add_column("name", style="white", width=30)
        table.add_column("jacc", justify="right")
        for _, row in df.iterrows():
            table.add_row(
                f"{row['id']}",
                str(row["source"]),
                f"[cyan]{row['name']}",
                f"{row['jacc']:.3f}",
            )
        console.print(table)
    else:
        sep = "," if args.out.lower().endswith(".csv") else "\t"
        df.to_csv(args.out, index=False, sep=sep)
        console.print(f"Wrote {len(df)} results to [green]{args.out}[/]")


def cmd_search_retro_compound(args: argparse.Namespace) -> None:
    """
    Search compounds by biosynthetic similarity to a given RetroMol fingerprint.

    :param args: command-line arguments
    """
    from bionexus.db.search import retro_search_compound

    rows = retro_search_compound(
        smiles=args.smiles,
        top_k=args.top_k,
        counted=getattr(args, "counted", False),
    )
    df = pd.DataFrame(rows)

    if not args.out:
        # console.print(df[["source", "name", "jacc"]])
        table = Table(show_header=True, header_style="bold magenta")
        table.add_column("compound_id", style="dim", width=6)
        table.add_column("source", style="dim", width=10)
        table.add_column("name", style="white", width=30)
        table.add_column("cosine", justify="right")
        for _, row in df.iterrows():
            table.add_row(
                f"{row['id']}",
                str(row["source"]),
                f"[cyan]{row['name']}",
                f"{row['cosine']:.3f}",
            )
        console.print(table)
    else:
        sep = "," if args.out.lower().endswith(".csv") else "\t"
        df.to_csv(args.out, index=False, sep=sep)
        console.print(f"Wrote {len(df)} results to [green]{args.out}[/]")


def cmd_search_retro_gbk(args: argparse.Namespace) -> None:
    """
    Search compounds by biosynthetic similarity to a GenBank file.

    :param args: command-line arguments
    """
    from bionexus.db.search import retro_search_gbk

    rows = retro_search_gbk(
        path=args.path,
        top_k=args.top_k,
        readout_toplevel=args.readout_toplevel,
        readout_sublevel=args.readout_sublevel,
        counted=getattr(args, "counted", False),
        cache_dir=getattr(args, "cache_dir", None),
    )
    df = pd.DataFrame(rows)

    if not args.out:
        # console.print(df[["source", "name", "jacc"]])
        table = Table(show_header=True, header_style="bold magenta")
        table.add_column("record", style="dim", width=30)
        table.add_column("compound_id", style="dim", width=6)
        table.add_column("source", style="dim", width=10)
        table.add_column("name", style="white", width=30)
        table.add_column("cosine", justify="right")
        for _, row in df.iterrows():
            table.add_row(
                f"{row['record']}",
                f"{row['id']}",
                str(row["source"]),
                f"[cyan]{row['name']}",
                f"{row['cosine']:.3f}",
            )
        console.print(table)
    else:
        sep = "," if args.out.lower().endswith(".csv") else "\t"
        df.to_csv(args.out, index=False, sep=sep)
        console.print(f"Wrote {len(df)} results to [green]{args.out}[/]")


def build_parser() -> argparse.ArgumentParser:
    """
    Build the argument parser for the Bionexus CLI.

    :return: configured ArgumentParser instance
    """
    p = argparse.ArgumentParser(prog="bionexus", description="Bionexus CLI (argparse)")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_init = sub.add_parser("init-db", help="Create schema (optionally drop first)")
    p_init.add_argument("--drop", action="store_true")
    p_init.set_defaults(func=cmd_init_db)

    p_up = sub.add_parser("upgrade", help="Run Alembic upgrade")
    p_up.add_argument("rev", nargs="?", default="head")
    p_up.set_defaults(func=cmd_upgrade)

    p_rev = sub.add_parser("revision", help="Create Alembic revision with --autogenerate")
    p_rev.add_argument("-m", "--message", required=True, help="Revision message")
    p_rev.add_argument("--head", default="head", help="Base head for the new revision (default: head)")
    p_rev.add_argument("--splice", action="store_true", help="Create a new branch from the given head")
    p_rev.add_argument("--rev-id", default=None, help="Hardcode a specific revision id")
    p_rev.add_argument(
        "--depends-on",
        action="append",
        default=[],
        help="Add depends_on revision id (repeatable)",
    )
    p_rev.add_argument(
        "--branch-label",
        action="append",
        default=[],
        help="Add branch label (repeatable)",
    )
    p_rev.set_defaults(func=cmd_revision)

    p_hist = sub.add_parser("history", help="Show Alembic history")
    p_hist.set_defaults(func=cmd_history)

    p_curr = sub.add_parser("current", help="Show current DB revision(s)")
    p_curr.set_defaults(func=cmd_current)

    p_heads = sub.add_parser("heads", help="Show head revision(s)")
    p_heads.set_defaults(func=cmd_heads)

    p_stamp_cmd = sub.add_parser("stamp", help="Stamp the database with a revision without running migrations")
    p_stamp_cmd.add_argument("rev", help="Revision to stamp (e.g. head, base, or a rev id)")
    p_stamp_cmd.set_defaults(func=cmd_stamp)

    p_down = sub.add_parser("downgrade", help="Run Alembic downgrade")
    p_down.add_argument("rev", help="Revision to downgrade to (e.g. -1, base, or a rev id)")
    p_down.set_defaults(func=cmd_downgrade)

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

    p_parse_compounds = sub.add_parser("parse-compounds", help="Parse compounds from database with RetroMol")
    p_parse_compounds.add_argument("--cache-dir", default=None, help="Cache/work dir for RetroMol")
    p_parse_compounds.add_argument("--batch", type=int, default=2000)
    p_parse_compounds.add_argument("--recompute", action="store_true", help="Force recomputation for all compounds")
    p_parse_compounds.add_argument("--workers", type=int, default=1, help="Number of parallel workers to use")
    p_parse_compounds.set_defaults(func=cmd_parse_compounds)

    p_fp = sub.add_parser("compute-fp-morgan", help="Compute fingerprints for compounds with SMILES")
    p_fp.add_argument("--batch", type=int, default=2000)
    p_fp.add_argument("--recompute", action="store_true", help="Force recomputation for all compounds")
    p_fp.set_defaults(func=cmd_compute_fp_morgan)

    p_do = sub.add_parser("compute-fp-retro-compound",help="Compute biosynthetic fingerprints for compounds")
    p_do.add_argument("--cache-dir", default=None, help="Cache/work dir for RetroMol")
    p_do.add_argument("--batch", type=int, default=2000)
    p_do.add_argument("--recompute", action="store_true", help="Force recomputation")
    p_do.set_defaults(func=cmd_compute_fp_retro_compound)

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

    p_search_r = sub.add_parser("search-retro-compound", help="Search compounds by biosynthetic similarity to a RetroMol fingerprint")
    p_search_r.add_argument("--smiles", required=True)
    p_search_r.add_argument("--top-k", type=int, default=20)
    p_search_r.add_argument("--out", default=None, help="Optional output file (TSV/CSV)")
    p_search_r.add_argument("--counted", action="store_true", help="Use counted fingerprint for search")
    p_search_r.set_defaults(func=cmd_search_retro_compound)

    p_search_r_gbk = sub.add_parser("search-retro-gbk", help="Search compounds by biosynthetic similarity to a GenBank file")
    p_search_r_gbk.add_argument("--path", required=True, help="Path to GenBank file")
    p_search_r_gbk.add_argument("--top-k", type=int, default=50)
    p_search_r_gbk.add_argument("--readout-toplevel", choices=["region", "cand_cluster"], default="region",
                              help="Read out fingerprints at 'region' or 'cand_cluster' level")
    p_search_r_gbk.add_argument("--readout-sublevel", choices=["rec", "gene"], default="rec",
                              help="Read out fingerprints at 'rec' or 'gene' level")
    p_search_r_gbk.add_argument("--out", default=None, help="Optional output file (TSV/CSV)")
    p_search_r_gbk.add_argument("--counted", action="store_true", help="Use counted fingerprint for search")
    p_search_r_gbk.add_argument("--cache-dir", default=None, help="Cache/work dir for RetroMol")
    p_search_r_gbk.set_defaults(func=cmd_search_retro_gbk)

    return p


def main(argv: list[str] | None = None) -> None:
    """
    Main entry point for the Bionexus CLI.

    :param argv: list of command-line arguments (defaults to sys.argv)
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
