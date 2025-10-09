from __future__ import annotations
import sys, tarfile, zipfile, gzip, shutil
import logging
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen, Request
from bionexus.config import default_cache_dir

logger = logging.getLogger(__name__)

def ensure_download(url: str, cache_dir: str | None = None, filename: str | None = None, force: bool = False) -> str:
    cache = Path(cache_dir) if cache_dir else default_cache_dir()
    cache.mkdir(parents=True, exist_ok=True)
    name = filename or Path(urlparse(url).path).name or "downloaded.file"
    dst = cache / name

    if dst.exists() and not force:
        logger.info(f"Using cached file: {dst}")
        return str(dst)

    label = f"Downloading {name}"
    logger.info(label)

    req = Request(url, headers={"User-Agent": "bionexus/1.0"})
    with urlopen(req) as r, open(dst, "wb") as f:
        total = r.length or r.headers.get("Content-Length")
        total = int(total) if total else None
        downloaded = 0
        chunk_size = 8192

        while True:
            chunk = r.read(chunk_size)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)

            if total:
                done = int(40 * downloaded / total)
                percent = downloaded / total * 100
                bar = "=" * done + " " * (40 - done)
                sys.stdout.write(f"\r[{bar}] {percent:5.1f}%")
                sys.stdout.flush()
        sys.stdout.write("\n")

    logger.info(f"Finished downloading: {dst}")
    return str(dst)

def extract_if_needed(path: str, cache_dir: str | None = None) -> list[str]:
    p = Path(path)
    outdir = (Path(cache_dir) if cache_dir else default_cache_dir()) / (p.stem + "_extracted")
    files: list[str] = []

    if zipfile.is_zipfile(p):
        outdir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(p) as z:
            z.extractall(outdir)
        files = [str(outdir / m) for m in z.namelist() if not m.endswith("/")]
    elif tarfile.is_tarfile(p):
        outdir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(p) as t:
            t.extractall(outdir)
        files = [str(outdir / m.name) for m in t.getmembers() if m.isfile()]
    elif p.suffix == ".gz" and not p.name.endswith(".tar.gz"):
        out = outdir.with_suffix("")  # pointless but keeps interface
        outdir.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(p, "rb") as src, open(out, "wb") as dst:
            shutil.copyfileobj(src, dst)
        files = [str(out)]
    else:
        files = []

    return files
