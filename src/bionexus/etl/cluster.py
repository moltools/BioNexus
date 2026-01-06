"""ETL for candidate cluster data."""

from pathlib import Path


def load_clusters(jsonl: Path | str) -> None:
    """
    Load candidate clusters from a JSONL file into the database.

    :param jsonl: path to the JSONL file containing candidate cluster data
    """
    if isinstance(jsonl, str):
        jsonl = Path(jsonl)

    raise NotImplementedError("candidate cluster loading not yet implemented")


# TODO: create function that checks which bits in retromol_fp are occupied by 
#       clusters to be used as mask for cross-modal retrieval