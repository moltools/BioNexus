"""Utility functions for setting up logging in BioNexus."""

from __future__ import annotations

import logging

from bionexus.config import DEFAULT_LOGGING_LVL


def setup_logging(level: int = DEFAULT_LOGGING_LVL) -> None:
    """
    Set up logging configuration for the BioNexus application.

    :param level: Logging level (default is DEFAULT_LOGGING_LVL)
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
