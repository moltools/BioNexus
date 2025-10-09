from __future__ import annotations
import logging
from bionexus.config import DEFAULT_LOGGING_LVL

def setup_logging(level: int = DEFAULT_LOGGING_LVL) -> None:
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
