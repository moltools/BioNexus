from __future__ import annotations
import gzip
import ijson
from statistics import mode
from typing import IO

def open_maybe_gzip(path: str, mode: str = "rt") -> IO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def iter_json(path: str):
    with open(path, "rb") as f:
        for item in ijson.items(f, "item"):
            yield item
