#!/usr/bin/env python3

import sys
from pathlib import Path

# Ensure repo root is importable when executing this script directly
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import argparse
import gzip
import math
import random
import shutil
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Iterable, Optional, Tuple

from scripts.utils.logger import TEPLogger, print_status, set_step_logger


DEFAULT_DRPVER = "v3_1_1"
DEFAULT_DAPVER = "3.1.0"
DEFAULT_DAPTYPE = "HYB10-MILESHC-MASTARSSP"
DEFAULT_BASE_URL = "https://data.sdss.org/sas/dr17/manga/spectro/analysis"
DEFAULT_DRPALL_BASE_URL = "https://data.sdss.org/sas/dr17/manga/spectro/redux"
DEFAULT_CMB_RA_DEG = 168.0
DEFAULT_CMB_DEC_DEG = -7.0


def _format_bytes(n: Optional[int]) -> str:
    if n is None:
        return "unknown"
    x = float(n)
    for unit in ["B", "KiB", "MiB", "GiB", "TiB"]:
        if x < 1024.0 or unit == "TiB":
            if unit == "B":
                return f"{int(x)} {unit}"
            return f"{x:.2f} {unit}"
        x /= 1024.0
    return f"{n} B"


def _head_content_length(url: str, timeout_s: float = 20.0) -> Optional[int]:
    import urllib.request

    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as r:
            v = r.headers.get("Content-Length")
            if v is None:
                return None
            return int(v)
    except Exception:
        return None


def _ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _http_to_rsync_url(url: str, rsync_host: str = "dtn.sdss.org") -> str:
    from urllib.parse import urlparse

    u = urlparse(url)
    path = u.path
    if path.startswith("/sas/"):
        path = path[len("/sas") :]
    return f"rsync://{rsync_host}{path}"


def _rsync_to_path(rsync_url: str, dest_path: Path) -> None:
    _ensure_parent_dir(dest_path)
    # rsync writes the file into the destination directory with the remote basename
    cmd = ["rsync", "-avz", "--no-motd", rsync_url, str(dest_path.parent) + "/"]
    print_status(f"rsync: {rsync_url}", "PROCESS")
    print_status(f"To: {dest_path.parent}", "PROCESS")
    import subprocess

    r = subprocess.run(cmd, check=True, capture_output=True, text=True)
    if not dest_path.exists():
        raise FileNotFoundError(f"rsync completed but file not found: {dest_path}")
    try:
        size_n = dest_path.stat().st_size
    except Exception:
        size_n = None
    if r.stdout:
        for line in r.stdout.strip().splitlines()[-2:]:
            print_status(f"rsync: {line}", "DEBUG")
    print_status(f"rsync complete: {dest_path.name} | size={_format_bytes(size_n)}", "SUCCESS")


def _download_to_path(
    url: str,
    dest_path: Path,
    *,
    progress_interval_s: float = 2.0,
    chunk_bytes: int = 1024 * 1024,
    retries: int = 3,
    retry_backoff_s: float = 2.0,
) -> None:
    import urllib.request

    _ensure_parent_dir(dest_path)

    tmp_path = dest_path.with_suffix(dest_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()

    print_status(f"Downloading: {url}", "PROCESS")
    print_status(f"To: {dest_path}", "PROCESS")

    last_err: Optional[BaseException] = None
    for attempt in range(1, max(retries, 1) + 1):
        try:
            with urllib.request.urlopen(url) as r, open(tmp_path, "wb") as f:
                expected = r.headers.get("Content-Length")
                expected_n = int(expected) if expected is not None else None
                if expected_n is not None:
                    print_status(f"Expected size: {_format_bytes(expected_n)}", "PROCESS")

                t0 = time.time()
                last_report = t0
                downloaded = 0

                while True:
                    chunk = r.read(chunk_bytes)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)

                    now = time.time()
                    if now - last_report >= progress_interval_s:
                        dt = max(now - t0, 1e-6)
                        rate = downloaded / dt
                        if expected_n is not None and expected_n > 0:
                            frac = min(max(downloaded / expected_n, 0.0), 1.0)
                            remain = max(expected_n - downloaded, 0)
                            eta_s = remain / max(rate, 1e-6)
                            print_status(
                                f"Progress: {frac*100:5.1f}% ({_format_bytes(downloaded)} / {_format_bytes(expected_n)}) | "
                                f"{_format_bytes(int(rate))}/s | ETA {eta_s:,.0f}s",
                                "PROCESS",
                            )
                        else:
                            print_status(
                                f"Progress: {_format_bytes(downloaded)} | {_format_bytes(int(rate))}/s",
                                "PROCESS",
                            )
                        last_report = now

                total_dt = max(time.time() - t0, 1e-6)
                rate = downloaded / total_dt
                print_status(
                    f"Downloaded: {_format_bytes(downloaded)} in {total_dt:,.1f}s ({_format_bytes(int(rate))}/s)",
                    "SUCCESS",
                )
            last_err = None
            break
        except Exception as e:
            last_err = e
            if tmp_path.exists():
                tmp_path.unlink()
            if attempt >= max(retries, 1):
                break
            sleep_s = retry_backoff_s * (2 ** (attempt - 1))
            print_status(f"Download failed (attempt {attempt}/{retries}): {e} | retrying in {sleep_s:.1f}s", "WARNING")
            time.sleep(sleep_s)

    if last_err is not None:
        raise last_err

    tmp_path.replace(dest_path)
    print_status(f"Download complete: {dest_path.name}", "SUCCESS")


def _maybe_gunzip_inplace(path: Path) -> Path:
    if not path.name.endswith(".gz"):
        return path

    out_path = path.with_suffix("")
    if out_path.exists():
        return out_path

    print_status(f"Decompressing: {path.name}", "PROCESS")
    with gzip.open(path, "rb") as fin, open(out_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)

    print_status(f"Decompressed to: {out_path.name}", "SUCCESS")
    return out_path


def construct_dapall_url(drpver: str, dapver: str, base_url: str = DEFAULT_BASE_URL) -> str:
    # DR17 SAS hosts dapall as an uncompressed FITS file (no .gz)
    return f"{base_url}/{drpver}/{dapver}/dapall-{drpver}-{dapver}.fits"


def construct_drpall_url(drpver: str, base_url: str = DEFAULT_DRPALL_BASE_URL) -> str:
    return f"{base_url}/{drpver}/drpall-{drpver}.fits"


def construct_maps_url(
    plate: int,
    ifudesign: int,
    daptype: str,
    drpver: str,
    dapver: str,
    base_url: str = DEFAULT_BASE_URL,
) -> str:
    plate_str = str(int(plate))
    ifu_str = str(int(ifudesign))
    maps_name = f"manga-{plate_str}-{ifu_str}-MAPS-{daptype}.fits.gz"
    return f"{base_url}/{drpver}/{dapver}/{daptype}/{plate_str}/{ifu_str}/{maps_name}"


def parse_plateifu(plateifu: str) -> Tuple[int, int]:
    parts = plateifu.strip().split("-")
    if len(parts) != 2:
        raise ValueError(f"Invalid PLATEIFU: {plateifu}")
    return int(parts[0]), int(parts[1])


def _load_plateifu_list(path: Path) -> Iterable[str]:
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            continue
        yield s


def _vec_from_radec_deg(ra_deg: float, dec_deg: float) -> tuple[float, float, float]:
    ra = math.radians(float(ra_deg))
    dec = math.radians(float(dec_deg))
    return (
        math.cos(dec) * math.cos(ra),
        math.cos(dec) * math.sin(ra),
        math.sin(dec),
    )


def _dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return float(a[0] * b[0] + a[1] * b[1] + a[2] * b[2])


def select_plateifus_from_dapall(
    dapall_fits_path: Path,
    daptype: str,
    n: int,
    seed: int,
) -> list[str]:
    try:
        from astropy.io import fits
    except Exception as e:
        raise RuntimeError(
            "Astropy is required to select plateifus from dapall. "
            "Either install astropy or pass --plateifu-list."
        ) from e

    with fits.open(dapall_fits_path, memmap=True) as hdul:
        if daptype not in hdul:
            available = [hdu.name for hdu in hdul]
            raise KeyError(f"DAPTYPE extension '{daptype}' not found in dapall. Available: {available}")

        data = hdul[daptype].data
        if data is None:
            raise RuntimeError(f"No table data found in dapall extension '{daptype}'")

        if "DAPDONE" in data.names:
            good = data[data["DAPDONE"] == 1]
        else:
            good = data

        if "PLATEIFU" not in good.names:
            raise KeyError(f"PLATEIFU column not found in dapall extension '{daptype}'")

        plateifus = [str(x).strip() for x in good["PLATEIFU"]]

    if not plateifus:
        raise RuntimeError("No plateifus available after selection")

    rng = random.Random(seed)
    rng.shuffle(plateifus)
    return plateifus[:n]


def select_plateifus_from_dapall_stratified_xcmb(
    dapall_fits_path: Path,
    daptype: str,
    n: int,
    seed: int,
    n_bins: int,
    cmb_ra_deg: float,
    cmb_dec_deg: float,
) -> list[str]:
    try:
        from astropy.io import fits
    except Exception as e:
        raise RuntimeError(
            "Astropy is required to select plateifus from dapall. "
            "Either install astropy or pass --plateifu-list."
        ) from e

    if n_bins < 2:
        n_bins = 2
    if n <= 0:
        return []

    cmb_vec = _vec_from_radec_deg(cmb_ra_deg, cmb_dec_deg)

    with fits.open(dapall_fits_path, memmap=True) as hdul:
        if daptype not in hdul:
            available = [hdu.name for hdu in hdul]
            raise KeyError(f"DAPTYPE extension '{daptype}' not found in dapall. Available: {available}")

        data = hdul[daptype].data
        if data is None:
            raise RuntimeError(f"No table data found in dapall extension '{daptype}'")

        if "DAPDONE" in data.names:
            good = data[data["DAPDONE"] == 1]
        else:
            good = data

        if "PLATEIFU" not in good.names:
            raise KeyError(f"PLATEIFU column not found in dapall extension '{daptype}'")

        ra_col = "IFURA" if "IFURA" in good.names else ("OBJRA" if "OBJRA" in good.names else None)
        dec_col = "IFUDEC" if "IFUDEC" in good.names else ("OBJDEC" if "OBJDEC" in good.names else None)
        if ra_col is None or dec_col is None:
            raise KeyError("Could not find RA/DEC columns (IFURA/IFUDEC or OBJRA/OBJDEC) in dapall")

        plateifus_all: list[str] = []
        x_all: list[float] = []
        for i in range(len(good)):
            pl = str(good["PLATEIFU"][i]).strip()
            try:
                ra = float(good[ra_col][i])
                dec = float(good[dec_col][i])
            except Exception:
                continue
            if not (math.isfinite(ra) and math.isfinite(dec)):
                continue
            x = _dot(_vec_from_radec_deg(ra, dec), cmb_vec)
            if not math.isfinite(x):
                continue
            plateifus_all.append(pl)
            x_all.append(float(x))

    if not plateifus_all:
        raise RuntimeError("No plateifus available after RA/DEC filtering")

    rng = random.Random(seed)

    edges = [-1.0 + 2.0 * k / float(n_bins) for k in range(n_bins + 1)]
    bins: list[list[int]] = [[] for _ in range(n_bins)]
    for idx, x in enumerate(x_all):
        if x <= -1.0:
            j = 0
        elif x >= 1.0:
            j = n_bins - 1
        else:
            j = int((x + 1.0) / 2.0 * n_bins)
            j = min(max(j, 0), n_bins - 1)
        bins[j].append(idx)

    for b in bins:
        rng.shuffle(b)

    target = int(n)
    per = target // n_bins
    rem = target - per * n_bins
    want = [per + (1 if k < rem else 0) for k in range(n_bins)]

    chosen: list[int] = []
    leftover: list[int] = []
    for j in range(n_bins):
        take = min(want[j], len(bins[j]))
        chosen.extend(bins[j][:take])
        leftover.extend(bins[j][take:])

    if len(chosen) < target:
        rng.shuffle(leftover)
        need = target - len(chosen)
        chosen.extend(leftover[:need])

    rng.shuffle(chosen)
    chosen = chosen[:target]
    return [plateifus_all[i] for i in chosen]


def main() -> None:
    parser = argparse.ArgumentParser(description="TEP-COS Step 1.0 - MaNGA data acquisition")

    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--drpver", default=DEFAULT_DRPVER)
    parser.add_argument("--dapver", default=DEFAULT_DAPVER)
    parser.add_argument("--daptype", default=DEFAULT_DAPTYPE)

    parser.add_argument(
        "--dapall-url",
        default=None,
        help="Optional explicit URL for dapall-*.fits(.gz). If omitted, a standard DR17 SAS URL is constructed.",
    )

    parser.add_argument(
        "--drpall-url",
        default=None,
        help="Optional explicit URL for drpall-*.fits. If omitted, a standard DR17 SAS URL is constructed.",
    )

    parser.add_argument(
        "--download-dapall",
        action="store_true",
        help="Download dapall summary table (recommended).",
    )

    parser.add_argument(
        "--download-drpall",
        action="store_true",
        help="Download drpall summary table (needed for distance-residual dipole analysis).",
    )

    parser.add_argument(
        "--plateifu-list",
        type=str,
        default=None,
        help="Path to a newline-delimited PLATEIFU list (e.g., 8485-1901). If provided, MAPS download will use this list.",
    )

    parser.add_argument(
        "--select-from-dapall",
        action="store_true",
        help="Select a random subset of PLATEIFU values from dapall (requires astropy).",
    )

    parser.add_argument(
        "--stratify-xcmb",
        action="store_true",
        help="When selecting from dapall, stratify the sample to balance x_CMB = n·n_CMB across [-1,1] bins (improves dipole sensitivity at fixed N).",
    )
    parser.add_argument("--xcmb-bins", type=int, default=10)
    parser.add_argument("--cmb-ra-deg", type=float, default=DEFAULT_CMB_RA_DEG)
    parser.add_argument("--cmb-dec-deg", type=float, default=DEFAULT_CMB_DEC_DEG)

    parser.add_argument("--sample-n", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)

    parser.add_argument(
        "--download-maps",
        action="store_true",
        help="Download MaNGA DAP MAPS files for the selected PLATEIFU subset.",
    )

    parser.add_argument(
        "--no-size-estimate",
        action="store_true",
        help="Disable HEAD-based size estimation (otherwise prints total download size when possible).",
    )
    parser.add_argument("--progress-interval", type=float, default=2.0)
    parser.add_argument("--chunk-mb", type=float, default=1.0)

    parser.add_argument(
        "--download-method",
        choices=["http", "rsync"],
        default="http",
        help="Download method. SDSS recommends rsync for bulk transfers.",
    )
    parser.add_argument(
        "--rsync-host",
        default="dtn.sdss.org",
        help="Rsync host to use (SDSS recommends dtn.sdss.org for bulk transfers).",
    )
    parser.add_argument("--max-workers", type=int, default=2)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-backoff", type=float, default=2.0)

    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    logger = TEPLogger(
        "step_1_0_data_acquisition",
        log_file_path=PROJECT_ROOT / "logs" / "step_1_0_data_acquisition.log",
    )
    set_step_logger(logger)

    data_dir = PROJECT_ROOT / "data"
    dapall_dir = data_dir / "dapall"
    drpall_dir = data_dir / "drpall"
    maps_dir = data_dir / "maps"
    outputs_dir = PROJECT_ROOT / "results" / "outputs"

    dapall_dir.mkdir(parents=True, exist_ok=True)
    drpall_dir.mkdir(parents=True, exist_ok=True)
    maps_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    dapall_url = args.dapall_url or construct_dapall_url(args.drpver, args.dapver, base_url=args.base_url)
    dapall_path = dapall_dir / Path(dapall_url).name

    drpall_url = args.drpall_url or construct_drpall_url(args.drpver)
    drpall_path = drpall_dir / Path(drpall_url).name

    dapall_fits_path: Optional[Path] = None

    estimate_sizes = not args.no_size_estimate
    chunk_bytes = max(int(args.chunk_mb * 1024 * 1024), 1024)
    progress_interval_s = max(float(args.progress_interval), 0.1)
    max_workers = max(int(args.max_workers), 1)

    if args.download_dapall:
        if dapall_path.exists() and not args.overwrite:
            print_status(f"dapall already exists, skipping: {dapall_path}", "WARNING")
        else:
            if estimate_sizes:
                size_n = _head_content_length(dapall_url)
                print_status(f"Planned download: dapall | size={_format_bytes(size_n)} | {dapall_path}", "PROCESS")
            if args.download_method == "rsync":
                rsync_url = _http_to_rsync_url(dapall_url, rsync_host=args.rsync_host)
                _rsync_to_path(rsync_url, dapall_path)
            else:
                _download_to_path(
                    dapall_url,
                    dapall_path,
                    progress_interval_s=progress_interval_s,
                    chunk_bytes=chunk_bytes,
                    retries=args.retries,
                    retry_backoff_s=args.retry_backoff,
                )

        dapall_fits_path = _maybe_gunzip_inplace(dapall_path)

    if args.download_drpall:
        if drpall_path.exists() and not args.overwrite:
            print_status(f"drpall already exists, skipping: {drpall_path}", "WARNING")
        else:
            if estimate_sizes:
                size_n = _head_content_length(drpall_url)
                print_status(f"Planned download: drpall | size={_format_bytes(size_n)} | {drpall_path}", "PROCESS")
            if args.download_method == "rsync":
                rsync_url = _http_to_rsync_url(drpall_url, rsync_host=args.rsync_host)
                _rsync_to_path(rsync_url, drpall_path)
            else:
                _download_to_path(
                    drpall_url,
                    drpall_path,
                    progress_interval_s=progress_interval_s,
                    chunk_bytes=chunk_bytes,
                    retries=args.retries,
                    retry_backoff_s=args.retry_backoff,
                )

    plateifus: list[str] = []

    if args.plateifu_list:
        plateifus = list(_load_plateifu_list(Path(args.plateifu_list)))
        print_status(f"Loaded PLATEIFU list: N={len(plateifus)}", "SUCCESS")

    if args.select_from_dapall:
        if dapall_fits_path is None:
            dapall_fits_path = _maybe_gunzip_inplace(dapall_path)
            if not dapall_fits_path.exists():
                raise FileNotFoundError(
                    "dapall FITS not found. Run with --download-dapall first or provide --dapall-url and --download-dapall."
                )

        if args.stratify_xcmb:
            plateifus = select_plateifus_from_dapall_stratified_xcmb(
                dapall_fits_path=dapall_fits_path,
                daptype=args.daptype,
                n=args.sample_n,
                seed=args.seed,
                n_bins=args.xcmb_bins,
                cmb_ra_deg=args.cmb_ra_deg,
                cmb_dec_deg=args.cmb_dec_deg,
            )
        else:
            plateifus = select_plateifus_from_dapall(
                dapall_fits_path=dapall_fits_path,
                daptype=args.daptype,
                n=args.sample_n,
                seed=args.seed,
            )
        print_status(f"Selected PLATEIFU subset from dapall: N={len(plateifus)}", "SUCCESS")

    if plateifus:
        selection_path = outputs_dir / "step_1_0_plateifu_selection.txt"
        selection_path.write_text("\n".join(plateifus) + "\n", encoding="utf-8")
        print_status(f"Saved PLATEIFU selection: {selection_path}", "SUCCESS")

    if args.download_maps:
        if not plateifus:
            raise RuntimeError(
                "No PLATEIFU values available. Provide --plateifu-list or use --select-from-dapall."
            )

        planned: list[tuple[str, Path, str]] = []
        skipped = 0
        for plateifu in plateifus:
            plate, ifu = parse_plateifu(plateifu)
            url = construct_maps_url(
                plate=plate,
                ifudesign=ifu,
                daptype=args.daptype,
                drpver=args.drpver,
                dapver=args.dapver,
                base_url=args.base_url,
            )
            out_path = maps_dir / args.daptype / str(plate) / str(ifu) / Path(url).name
            if out_path.exists() and not args.overwrite:
                skipped += 1
                continue
            planned.append((url, out_path, plateifu))

        print_status(
            f"MAPS download plan: total={len(plateifus)} | to_download={len(planned)} | skipped_existing={skipped}",
            "TITLE",
        )
        print_status(f"MAPS destination root: {maps_dir / args.daptype}", "PROCESS")
        if estimate_sizes and planned:
            total_n: int = 0
            known = 0
            for i, (url, out_path, plateifu) in enumerate(planned, start=1):
                n = _head_content_length(url)
                if n is not None:
                    total_n += n
                    known += 1
                print_status(
                    f"[{i}/{len(planned)}] {plateifu} | size={_format_bytes(n)} | {out_path.name}",
                    "PROCESS",
                )
            if known:
                print_status(
                    f"Estimated MAPS total: {_format_bytes(total_n)} (known={known}/{len(planned)})",
                    "SUCCESS",
                )
            else:
                print_status("Estimated MAPS total: unknown (no Content-Length available)", "WARNING")

        for i, (url, out_path, plateifu) in enumerate(planned, start=1):
            print_status(f"[{i}/{len(planned)}] Queued MAPS for {plateifu}", "PROCESS")

        def _fetch_one(item: tuple[str, Path, str]) -> str:
            url, out_path, plateifu = item
            if args.download_method == "rsync":
                rsync_url = _http_to_rsync_url(url, rsync_host=args.rsync_host)
                _rsync_to_path(rsync_url, out_path)
            else:
                _download_to_path(
                    url,
                    out_path,
                    progress_interval_s=progress_interval_s,
                    chunk_bytes=chunk_bytes,
                    retries=args.retries,
                    retry_backoff_s=args.retry_backoff,
                )
            return plateifu

        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(_fetch_one, item): item for item in planned}
            for j, fut in enumerate(as_completed(futures), start=1):
                url, out_path, plateifu = futures[fut]
                try:
                    fut.result()
                    print_status(
                        f"Completed MAPS [{j}/{len(planned)}]: {plateifu} -> {out_path.name}",
                        "SUCCESS",
                    )
                except Exception as e:
                    print_status(
                        f"Failed MAPS [{j}/{len(planned)}]: {plateifu} | {e}",
                        "ERROR",
                    )
                    raise

    print_status("Step 1.0 complete", "SUCCESS")


if __name__ == "__main__":
    main()
