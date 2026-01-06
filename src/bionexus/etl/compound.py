"""ETL for compound data."""

from pathlib import Path


def load_compounds(jsonl: Path | str) -> None:
    """
    Load compounds from a JSONL file into the database.

    :param jsonl: path to the JSONL file containing compound data
    """
    raise NotImplementedError("compound loading not yet implemented")
