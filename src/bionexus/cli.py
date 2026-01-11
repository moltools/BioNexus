"""Command line interface for interacting with the BioNexus database."""

from dotenv import load_dotenv
from pathlib import Path

try:
    load_dotenv(Path(__file__).parent.parent.parent / ".env")
except Exception:
    exit("could not load .env file")

import argparse
import logging
import subprocess

from alembic import command
from alembic.config import Config

from bionexus.version import __version__
from bionexus.utils.logging import setup_logging
from bionexus.etl.annotation import load_annotations
from bionexus.etl.compound import load_compounds
from bionexus.etl.cluster import load_clusters


log = logging.getLogger(__name__)


def get_alembic_config() -> Config:
    """
    Get the Alembic configuration for database migrations.

    :return: Alembic configuration object
    """
    root = Path(__file__).resolve().parents[2]  # BioNexus/
    cfg = Config(str(root / "alembic.ini"))

    # Make script_location absolute so Alembic never has to guess where it is
    cfg.set_main_option("script_location", str(root / "alembic"))

    return cfg


def cmd_upgrade(args: argparse.Namespace) -> None:
    """
    Upgrade the database to the specified revision.
    """
    command.upgrade(get_alembic_config(), args.revision)


def cmd_downgrade(args: argparse.Namespace) -> None:
    """
    Downgrade the database to the specified revision.
    """
    command.downgrade(get_alembic_config(), args.revision)


def cmd_dump(args: argparse.Namespace) -> None:
    """
    Dump the database to a pg_dump custom-format file.
    """
    out = Path(args.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)

    container = "bionexus-db-1"
    tmp_path = "/tmp/bionexus.dump"

    subprocess.run(
        [
            "docker", "exec", "-i", container,
            "pg_dump",
            "-U", "bionexus",
            "-Fc",
            "-d", "bionexus",
            "-f", tmp_path,
        ],
        check=True,
    )

    subprocess.run(
        ["docker", "cp", f"{container}:{tmp_path}", str(out)],
        check=True,
    )

def cmd_load_annotations(args: argparse.Namespace) -> None:
    """
    Load annotations into database from a specified file.
    """
    filepath = Path(args.file).expanduser()
    load_annotations(filepath=filepath)


def cmd_load_compounds(args: argparse.Namespace) -> None:
    """
    Load compounds into database from a specifified file.
    """
    jsonl_path = Path(args.jsonl).expanduser()
    load_compounds(jsonl=jsonl_path)


def cmd_load_clusters(args: argparse.Namespace) -> None:
    """
    Load candidate clusters into database from a specified file.
    """
    jsonl_path = Path(args.jsonl).expanduser()
    load_clusters(jsonl=jsonl_path)


def cli() -> argparse.ArgumentParser:
    """
    Parse command line arguments for the BioNexus CLI.

    :return: parsed command line arguments
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(title="command", required=True)

    # Alembic commands
    alembic_parser = subparsers.add_parser("alembic", help="database migration commands")
    alembic_sub = alembic_parser.add_subparsers(dest="action", required=True)

    upgrade_parser = alembic_sub.add_parser("upgrade", help="upgrade the database")
    upgrade_parser.add_argument("--revision", nargs="?", default="head", help="target revision (default: head)")
    upgrade_parser.set_defaults(func=cmd_upgrade)

    downgrade_parser = alembic_sub.add_parser("downgrade", help="downgrade the database")
    downgrade_parser.add_argument("--revision", nargs="?", default="-1", help="target revision (default: -1)")
    downgrade_parser.set_defaults(func=cmd_downgrade)

    # Dump database command
    dump_parser = subparsers.add_parser("dump", help="dump the database using pg_dump")
    dump_parser.add_argument("--out", help="output file path for the database dump")
    dump_parser.set_defaults(func=cmd_dump)

    # Load data
    load_parser = subparsers.add_parser("load", help="load data into the database")
    load_sub = load_parser.add_subparsers(title="data type", required=True)

    # Load compounds command
    load_compounds_parser = load_sub.add_parser("compounds", help="load compounds into the database")
    load_compounds_parser.add_argument("--jsonl", help="file containing compounds to load")
    load_compounds_parser.set_defaults(func=cmd_load_compounds)

    # Load candidate clusters command
    load_clusters_parser = load_sub.add_parser("clusters", help="load candidate clusters into the database")
    load_clusters_parser.add_argument("--jsonl", help="file containing candidate clusters to load")
    load_clusters_parser.set_defaults(func=cmd_load_clusters)

    # Load annotations command
    load_annotations_parser = load_sub.add_parser("annotations", help="load annotations into the database")
    load_annotations_parser.add_argument("--file", help="file containing annotations to load")
    load_annotations_parser.set_defaults(func=cmd_load_annotations)

    return parser


def main(argv: list[str] | None = None) -> None:
    """
    Entry point for the BioNexus CLI.

    :param argv: command line arguments
    """
    setup_logging(level="DEBUG")
    parser = cli()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
