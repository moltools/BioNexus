"""Configuration for Bionexus package."""

import os
from pathlib import Path
from platform import system

# Load .env from repo root
try:
    from dotenv import load_dotenv

    load_dotenv(Path(__file__).resolve().parent.parent.parent / ".env")
except Exception:
    pass

DEFAULT_NPATLAS_URL = os.getenv("NPATLAS_URL", "").strip()
DEFAULT_MIBIG_JSON_URL = os.getenv("MIBIG_JSON_URL", "").strip()
DEFAULT_MIBIG_GBK_URL = os.getenv("MIBIG_GBK_URL", "").strip()

DEFAULT_MAX_BYTES_GBK = int(os.getenv("MAX_BYTES_GBK", "1048576").strip())  # default 1 MiB

DEFAULT_LOGGING_LVL = os.getenv("BIONEXUS_LOGGING_LVL", "INFO").strip().upper()


def default_cache_dir() -> Path:
    """
    Get the default cache directory for Bionexus.

    :return: Path to the cache directory.
    """
    root = os.getenv("BIONEXUS_CACHE_DIR")
    if root:
        return Path(root)
    s = system()
    if s == "Windows":
        base = Path(os.getenv("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
        return base / "Bionexus" / "Cache"
    elif s == "Darwin":
        return Path.home() / "Library" / "Caches" / "Bionexus"
    else:
        return Path(os.getenv("XDG_CACHE_HOME", Path.home() / ".cache")) / "bionexus"
