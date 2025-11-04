"""Utility functions for file I/O operations, including handling gzip files and iterating over JSON items."""

from __future__ import annotations

import gzip
from collections.abc import Generator
from typing import IO

import ijson


def open_maybe_gzip(path: str, mode: str = "rt") -> IO:
    """
    Open a file that may be gzip-compressed.

    :param path: path to the file
    :param mode: mode to open the file
    :return: file object
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode)

    return open(path, mode)


def iter_json(path: str) -> Generator[dict, None, None]:
    """
    Iterate over items in a JSON file.

    :param path: path to the JSON file
    :yield: items from the JSON file
    """
    with open(path, "rb") as f:
        for item in ijson.items(f, "item"):
            yield from item
