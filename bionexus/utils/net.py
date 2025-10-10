from __future__ import annotations
import re, sys, tarfile, zipfile, gzip, shutil
import logging
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen, Request
from bionexus.config import default_cache_dir

logger = logging.getLogger(__name__)

def ensure_download(url: str, cache_dir: str | None = None, filename: str | None = None, force: bool = False) -> str:
    cache = Path(cache_dir) if cache_dir else default_cache_dir()
    cache.mkdir(parents=True, exist_ok=True)
    name = filename or Path(urlparse(url).path).name
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

def _default_cache_dir() -> Path:
    # use your existing default_cache_dir() if you have one
    from bionexus.utils.net import default_cache_dir  # or adjust import
    return Path(default_cache_dir())

def _safe_join(base: Path, *paths: str) -> Path:
    p = (base / Path(*paths)).resolve()
    base_resolved = base.resolve()
    if not str(p).startswith(str(base_resolved)):
        raise ValueError(f"Unsafe path detected: {p}")
    return p

def _safe_extract_tar(t: tarfile.TarFile, outdir: Path) -> None:
    for m in t.getmembers():
        # skip absolute paths and parent traversal
        if m.name.startswith("/") or ".." in Path(m.name).parts:
            continue
        target = _safe_join(outdir, m.name)
        if m.isdir():
            target.mkdir(parents=True, exist_ok=True)
        elif m.issym() or m.islnk():
            # ignore links for safety; adapt if you need them
            continue
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            with t.extractfile(m) as src, open(target, "wb") as dst:
                if src is not None:
                    shutil.copyfileobj(src, dst)

def _safe_extract_zip(z: zipfile.ZipFile, outdir: Path) -> None:
    for name in z.namelist():
        if name.endswith("/"):
            # directory entry
            d = _safe_join(outdir, name)
            d.mkdir(parents=True, exist_ok=True)
            continue
        if name.startswith("/") or ".." in Path(name).parts:
            continue
        target = _safe_join(outdir, name)
        target.parent.mkdir(parents=True, exist_ok=True)
        with z.open(name) as src, open(target, "wb") as dst:
            shutil.copyfileobj(src, dst)

def _base_for_outdir(p: Path) -> str:
    name = p.name
    # strip double/compound suffixes first
    for pat in (r"\.tar\.(gz|bz2|xz)$", r"\.t(gz|bz2|xz)$"):
        if re.search(pat, name, flags=re.IGNORECASE):
            return re.sub(pat, "", name, flags=re.IGNORECASE)
    # simple single-suffix cases
    for suf in (".zip", ".gz", ".bz2", ".xz"):
        if name.lower().endswith(suf):
            return name[: -len(suf)]
    return p.stem

def extract_if_needed(path: str, cache_dir: str | None = None) -> list[str]:
    p = Path(path)
    outroot = Path(cache_dir) if cache_dir else _default_cache_dir()
    outdir = outroot / (f"{_base_for_outdir(p)}_extracted")
    outdir.mkdir(parents=True, exist_ok=True)

    files: list[str] = []

    # ZIP
    if zipfile.is_zipfile(p):
        with zipfile.ZipFile(p) as z:
            members = [n for n in z.namelist() if not n.endswith("/")]
            _safe_extract_zip(z, outdir)
        files = [str(outdir / m) for m in members]
        return files

    # TAR family
    if tarfile.is_tarfile(p):
        # open with r:* to auto-detect compression
        with tarfile.open(p, mode="r:*") as t:
            members = [m.name for m in t.getmembers() if m.isfile()]
            _safe_extract_tar(t, outdir)
        files = [str(outdir / m) for m in members]
        return files

    # Plain single-file .gz (NOT .tar.gz)
    name_lower = p.name.lower()
    if name_lower.endswith(".gz") and not (name_lower.endswith(".tar.gz") or name_lower.endswith(".tgz")):
        # write decompressed file into the outdir with the .gz stripped
        out_file = outdir / _base_for_outdir(p)
        out_file.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(p, "rb") as src, open(out_file, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return [str(out_file)]

    # Not an archive we handle
    return []
