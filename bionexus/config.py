import os
from pathlib import Path
from platform import system

# load .env from repo root
try:
    from dotenv import load_dotenv
    load_dotenv(Path(__file__).resolve().parent.parent / ".env")
except Exception:
    pass

DEFAULT_NPATLAS_URL = os.getenv("NPATLAS_URL", "").strip()
DEFAULT_LOGGING_LVL = os.getenv("BIONEXUS_LOGGING_LVL", "INFO").strip().upper()

def default_cache_dir() -> Path:
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
