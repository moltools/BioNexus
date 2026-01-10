"""ETL for annotation data."""

from pathlib import Path


def load_annotations(jsonl: Path | str) -> None:
    """
    Load annotations from a JSONL file into the database.

    :param jsonl: path to the JSONL file containing annotation data
    """
    if isinstance(jsonl, str):
        jsonl = Path(jsonl)

    raise NotImplementedError("annotation loading not yet implemented")
