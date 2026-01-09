#!/usr/bin/env python3
"""
Step 3.16 (self-extraction): PG1115 Pan-STARRS1 multi-band light curves

Purpose
- Retrieve PS1 single-epoch warp images (g/r/i/z/y) around PG1115+080.
- Perform forced photometry at fixed image positions (A/B/C/D) from priors.
- Output per-band, per-image light curves suitable for temporal-shear / chromaticity tests.

Inputs required
- Priors file (JSON) with sky positions (ICRS degrees) for images A,B,C,D.
  Example (save as data/priors/pg1115_hst_priors.json):
  {
    "A": {"ra": 169.57083, "dec": 7.76614},
    "B": {"ra": 169.56944, "dec": 7.76642},
    "C": {"ra": 169.57023, "dec": 7.76525},
    "D": {"ra": 169.57150, "dec": 7.76574}
  }
  (These are placeholders; replace with HST-fit coordinates you trust.)

Outputs
- data/selfextract/pg1115_ps1/lightcurves_pg1115_ps1_{band}.csv
  columns: band, image, mjd, flux, flux_err, mag, mag_err, zp, airmass, exposure, fname

Notes
- Photometry: simple circular aperture on background-subtracted image.
  (First-pass; swap in PSF/deconvolution later if needed.)
- Zeropoint: uses header MAGZP/MAGZERO/MAGZERO_RP; if missing, stores NaN mag.
- This script prioritizes reproducibility over absolute calibration; QC and
  PSF-model refinement can be added iteratively.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List
from urllib.request import urlretrieve

import numpy as np
from astropy.io import fits
from astropy.table import unique
from astropy.time import Time
from astropy.wcs import WCS
from astroquery.mast import Observations
from photutils.aperture import ApertureStats, CircularAnnulus, CircularAperture, aperture_photometry
from photutils.detection import DAOStarFinder

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("pg1115_ps1")


def load_priors(path: Path) -> Dict[str, Dict[str, float]]:
    with open(path, "r") as f:
        priors = json.load(f)
    # For PG1115, ground-based monitoring typically uses A (A1+A2 blended), B, C.
    # We allow an optional D for other systems.
    required = {"A", "B", "C"}
    missing = required - set(priors.keys())
    if missing:
        raise ValueError(f"Priors missing keys: {missing}")
    return priors


def save_priors(path: Path, priors: Dict[str, Dict[str, float]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(priors, f, indent=2, sort_keys=True)
        f.write("\n")


def download_hst_cutout(url: str, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    local = outdir / Path(url).name
    if local.exists():
        return local
    log.info(f"Downloading HST cutout: {url}")
    urlretrieve(url, local)
    return local


def robust_bkg_rms(data: np.ndarray) -> float:
    finite = np.isfinite(data)
    if finite.sum() < 50:
        return float("nan")
    vals = np.asarray(data[finite], dtype=float)
    med = np.median(vals)
    mad = np.median(np.abs(vals - med))
    sigma = 1.4826 * mad
    if np.isfinite(sigma) and sigma > 0:
        return float(sigma)

    # Fallback: percentile-based robust sigma (works when many values are exactly 0).
    p16, p84 = np.percentile(vals, [16.0, 84.0])
    sigma = 0.5 * (p84 - p16)
    if np.isfinite(sigma) and sigma > 0:
        return float(sigma)

    # Final fallback: clipped standard deviation.
    lo, hi = np.percentile(vals, [5.0, 95.0])
    clipped = vals[(vals >= lo) & (vals <= hi)]
    if clipped.size < 50:
        return float("nan")
    sigma = float(np.std(clipped))
    return sigma if np.isfinite(sigma) and sigma > 0 else float("nan")


def derive_priors_from_hst_cutout(
    hst_fits: Path,
    n_images: int = 4,
    fwhm_pix: float = 2.5,
    nsigma: float = 10.0,
) -> Dict[str, Dict[str, float]]:
    """Detect the brightest point sources in an HST cutout and return A/B/C/D priors.

    The A/B/C/D labels are assigned deterministically by sorting the detected sources
    by (Dec descending, RA ascending). This is a reproducible but not literature-anchored
    labeling.
    """
    data, hdr = read_first_image_hdu(hst_fits)
    wcs = WCS(hdr)
    if not wcs.has_celestial:
        raise RuntimeError("HST image does not contain a celestial WCS; cannot derive RA/Dec priors.")

    bkg = np.nanmedian(data)
    img = data - bkg
    rms = robust_bkg_rms(img)
    if not np.isfinite(rms) or rms <= 0:
        raise RuntimeError("Could not estimate background RMS from HST cutout.")

    finder = DAOStarFinder(fwhm=fwhm_pix, threshold=nsigma * rms)
    sources = finder(img)
    if sources is None or len(sources) < n_images:
        raise RuntimeError(
            f"DAOStarFinder found {0 if sources is None else len(sources)} sources; need >= {n_images}."
        )

    # Take the n brightest sources.
    sources.sort("flux")
    sources = sources[::-1][:n_images]

    ras, decs = wcs.all_pix2world(sources["xcentroid"], sources["ycentroid"], 0)

    det = []
    for ra, dec in zip(ras, decs):
        det.append((float(ra), float(dec)))

    # For PG1115, A1/A2 are a close pair that are not resolved in typical ground-based data.
    # We therefore merge the closest pair into a single blended component "A".
    det_arr = np.array(det, dtype=float)  # shape (4,2): ra, dec
    # Pairwise angular distance proxy in degrees (small-angle Euclidean on sphere ok here).
    dmin = None
    imin = jmin = None
    for i in range(len(det_arr)):
        for j in range(i + 1, len(det_arr)):
            dra = (det_arr[i, 0] - det_arr[j, 0]) * np.cos(np.deg2rad(0.5 * (det_arr[i, 1] + det_arr[j, 1])))
            ddec = det_arr[i, 1] - det_arr[j, 1]
            dij = np.hypot(dra, ddec)
            if dmin is None or dij < dmin:
                dmin = dij
                imin, jmin = i, j

    if imin is None or jmin is None:
        raise RuntimeError("Could not identify closest pair for PG1115 A1/A2 merge.")

    a_ra = float(np.mean([det_arr[imin, 0], det_arr[jmin, 0]]))
    a_dec = float(np.mean([det_arr[imin, 1], det_arr[jmin, 1]]))
    remaining = [k for k in range(len(det_arr)) if k not in (imin, jmin)]
    rem = [(float(det_arr[k, 0]), float(det_arr[k, 1])) for k in remaining]
    # Deterministic B/C labeling: sort by Dec desc then RA asc.
    rem.sort(key=lambda t: (-t[1], t[0]))
    priors = {
        "A": {"ra": a_ra, "dec": a_dec},
        "B": {"ra": rem[0][0], "dec": rem[0][1]},
        "C": {"ra": rem[1][0], "dec": rem[1][1]},
    }
    return priors


def download_hst_drz_from_mast(
    ra: float,
    dec: float,
    radius_arcsec: float,
    outdir: Path,
    max_products: int = 4,
) -> Path:
    """Download an HST drizzled image with celestial WCS from MAST and return local path."""
    outdir.mkdir(parents=True, exist_ok=True)
    obs = Observations.query_region(f"{ra} {dec}", radius=f"{radius_arcsec} arcsec")
    if obs is None or len(obs) == 0:
        raise RuntimeError("No MAST observations found for HST query_region.")
    if "obs_collection" in obs.colnames:
        keep = np.array([str(x).lower() == "hst" for x in obs["obs_collection"]], dtype=bool)
        obs = obs[keep]
    if obs is None or len(obs) == 0:
        raise RuntimeError("No HST observations found in MAST region.")

    products = Observations.get_product_list(obs)
    if products is None or len(products) == 0:
        raise RuntimeError("No HST products found for observations.")

    # Prefer drizzled science images.
    if "productFilename" in products.colnames:
        fn = np.array([str(x).lower() for x in products["productFilename"]], dtype=str)
        is_fits = np.char.endswith(fn, ".fits")
        is_drz = np.array([("_drz" in x) or ("_drc" in x) for x in fn], dtype=bool)
        products = products[is_fits & is_drz]

    if products is None or len(products) == 0:
        raise RuntimeError("No HST drizzled FITS products (_drz/_drc) found in MAST products.")

    if len(products) > max_products:
        products = products[:max_products]

    manifest = Observations.download_products(
        products,
        mrp_only=False,
        download_dir=str(outdir),
        cache=True,
        curl_flag=False,
    )
    if manifest is None or len(manifest) == 0 or "Local Path" not in manifest.colnames:
        raise RuntimeError("HST MAST download produced an empty manifest.")

    # Return the first FITS path.
    for p in manifest["Local Path"]:
        if p and str(p).lower().endswith(".fits"):
            return Path(str(p))
    raise RuntimeError("HST MAST download did not yield any FITS local paths.")


def priors_are_sane(priors: Dict[str, Dict[str, float]]) -> bool:
    try:
        for k in ["A", "B", "C"]:
            ra = float(priors[k]["ra"])
            dec = float(priors[k]["dec"])
            if not np.isfinite(ra) or not np.isfinite(dec):
                return False
            if ra < 0.0 or ra >= 360.0:
                return False
            if dec < -90.0 or dec > 90.0:
                return False
        return True
    except Exception:
        return False


def derive_priors_for_pg1115(args, outdir: Path) -> Dict[str, Dict[str, float]]:
    """Derive PG1115 priors from an HST image with celestial WCS.

    Preference order:
    1) User-provided --hst-cutout-local
    2) Download from --hst-cutout-url (may fail if file lacks WCS)
    3) Download an HST drizzled image from MAST
    """
    if args.hst_cutout_local is not None:
        hst_fits = args.hst_cutout_local
        return derive_priors_from_hst_cutout(hst_fits)

    # Try URL first (historical behavior)
    try:
        hst_fits = download_hst_cutout(args.hst_cutout_url, outdir / "priors")
        return derive_priors_from_hst_cutout(hst_fits)
    except Exception as e:
        log.warning("Could not derive priors from HST URL cutout (%s): %s", args.hst_cutout_url, e)

    # Fallback: download a drizzled HST science image with WCS via MAST
    hst_fits = download_hst_drz_from_mast(args.ra, args.dec, radius_arcsec=60.0, outdir=outdir / "priors")
    return derive_priors_from_hst_cutout(hst_fits)


def search_ps1(ra: float, dec: float, radius_arcmin: float, filters: List[str], max_products: int):
    # MAST schema and accepted query_criteria keys can vary between releases.
    # Start with a strict query_criteria attempt; if it yields no results,
    # fall back to a broad query_region + filtering.
    obs = Observations.query_criteria(
        coordinates=f"{ra} {dec}",
        radius=f"{radius_arcmin} arcmin",
        obs_collection="PS1",
        dataproduct_type="image",
        calib_level=[1, 2],
    )

    if obs is None or len(obs) == 0:
        radius_deg = radius_arcmin / 60.0
        log.warning(
            "No results for strict PS1 query_criteria; trying query_region fallback (radius=%.4f deg).",
            radius_deg,
        )
        obs = Observations.query_region(f"{ra} {dec}", radius=f"{radius_deg} deg")
        if obs is None or len(obs) == 0:
            raise RuntimeError(
                "No MAST observations found in region; network/MAST issue or target/radius problem. "
                "Try increasing --radius or check MAST availability."
            )

        if "obs_collection" in obs.colnames:
            keep = []
            for v in obs["obs_collection"]:
                s = str(v).lower()
                keep.append(("ps1" in s) or ("panstarr" in s) or ("pan-starr" in s))
            obs = obs[np.array(keep, dtype=bool)]

        if "dataproduct_type" in obs.colnames:
            obs = obs[obs["dataproduct_type"] == "image"]

    if obs is None or len(obs) == 0:
        raise RuntimeError("No PS1 observations found after fallback filtering.")

    # Filter on requested bandpasses if possible.
    if filters and "filters" in obs.colnames:
        obs_keep = []
        want = {f.strip() for f in filters}
        for v in obs["filters"]:
            obs_keep.append(str(v).strip() in want)
        obs = obs[np.array(obs_keep, dtype=bool)]

    if obs is None or len(obs) == 0:
        raise RuntimeError(
            "MAST returned PS1 observations, but none match requested filters. "
            "Try --filters r,i,z (or omit --filters)."
        )

    products = Observations.get_product_list(obs)
    # Keep warp products only
    if "productSubGroupDescription" in products.colnames:
        mask = np.array([str(x).lower() == "warp" for x in products["productSubGroupDescription"]], dtype=bool)
        products = products[mask]
    if len(products) == 0:
        raise RuntimeError(
            "No PS1 warp products found (productSubGroupDescription != 'warp'). "
            "MAST may not be exposing warp products for this target; consider alternative PS1 endpoints."
        )

    # Prefer science image products and exclude common ancillary products.
    # The PS1 warp bundles often include *mask* and *wt* images.
    if "productFilename" in products.colnames:
        fn = np.array([str(x) for x in products["productFilename"]], dtype=str)
        # Keep only FITS images (exclude metadata products like .cmf/.mdc).
        is_fits = np.char.endswith(np.char.lower(fn), ".fits")
        bad = (
            np.char.endswith(fn, ".mask.fits")
            | np.char.endswith(fn, ".wt.fits")
            | np.char.endswith(fn, ".exp.fits")
            | np.char.endswith(fn, ".invvar.fits")
        )
        products = products[is_fits & (~bad)]
    if len(products) == 0:
        raise RuntimeError(
            "No PS1 science warp FITS images remained after filtering (only ancillary/metadata products were available)."
        )

    # Deduplicate per observation / filter. MAST column names can vary.
    obsid_key = None
    for k in ["obsid", "obsID", "obs_id", "parent_obsid", "parent_obsID"]:
        if k in products.colnames:
            obsid_key = k
            break

    filt_key = None
    for k in ["filters", "filter", "FILTER"]:
        if k in products.colnames:
            filt_key = k
            break

    keys = []
    if obsid_key is not None:
        keys.append(obsid_key)
    if filt_key is not None:
        keys.append(filt_key)
    if not keys and "productFilename" in products.colnames:
        keys = ["productFilename"]

    if keys:
        products = unique(products, keys=keys, keep="first")
    else:
        log.warning("Could not find suitable de-duplication keys in MAST products: %s", products.colnames)
    # Limit size
    if len(products) > max_products:
        products = products[:max_products]
    log.info(f"Selected {len(products)} PS1 warp products.")
    return products


def download_products(products, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = Observations.download_products(
        products,
        mrp_only=False,
        download_dir=str(outdir),
        cache=True,
        curl_flag=False,
    )
    # Keep local file paths only
    paths: List[Path] = []
    if manifest is None or len(manifest) == 0:
        return paths
    if "Local Path" not in manifest.colnames:
        raise RuntimeError("MAST download manifest missing 'Local Path' column.")

    for p in manifest["Local Path"]:
        if p and str(p).endswith(".fits"):
            paths.append(Path(str(p)))
    log.info(f"Downloaded {len(paths)} FITS files.")
    return paths


def read_first_image_hdu(fname: Path):
    """Return (data, header) from the first HDU that contains 2D image data."""
    with fits.open(fname) as hdul:
        for hdu in hdul:
            if hdu.data is None:
                continue
            data = np.asarray(hdu.data)
            if data.ndim == 2:
                return data.astype(float), hdu.header
    raise RuntimeError("No 2D image HDU found.")


def estimate_pixel_scale(wcs: WCS) -> float:
    cd = wcs.pixel_scale_matrix
    # arcsec/pixel from matrix diagonal
    scale_deg = np.sqrt(np.abs(cd[0, 0] * cd[1, 1]))
    return scale_deg * 3600.0


def photometer_file(fname: Path, priors: Dict[str, Dict[str, float]], aperture_arcsec: float = 1.2):
    # Skip ancillary PS1 products defensively.
    s = fname.name.lower()
    if s.endswith(".mask.fits") or s.endswith(".wt.fits") or s.endswith(".invvar.fits"):
        raise RuntimeError("Ancillary PS1 product (mask/weight/invvar), not a science image.")

    data, hdr = read_first_image_hdu(fname)
    wcs = WCS(hdr)

    # NaNs are common in warps; mask them so photutils can ignore rather than propagate.
    mask = ~np.isfinite(data)
    bkg_rms_global = robust_bkg_rms(data)
    global_med = float(np.nanmedian(data[~mask])) if np.any(~mask) else np.nan

    # Pixel scale and aperture
    pixscale = estimate_pixel_scale(wcs)
    r_pix = aperture_arcsec / pixscale
    # Background annulus (local sky)
    r_in = 2.0 * r_pix
    r_out = 3.5 * r_pix

    labels = sorted(priors.keys())

    # Preserve conventional ordering if present.
    preferred = [k for k in ["A", "B", "C", "D"] if k in priors]
    remaining = [k for k in labels if k not in preferred]
    labels = preferred + remaining

    # Build sky coords
    ras = [priors[k]["ra"] for k in labels]
    decs = [priors[k]["dec"] for k in labels]
    x, y = wcs.all_world2pix(ras, decs, 0)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if (not np.isfinite(x).all()) or (not np.isfinite(y).all()):
        raise RuntimeError("Non-finite WCS projection for priors (target off-chip or invalid WCS/priors).")

    ny, nx = data.shape
    if (x < 0).any() or (x >= nx).any() or (y < 0).any() or (y >= ny).any():
        raise RuntimeError("Projected priors fall outside this image footprint (off-chip).")
    positions = list(zip(x, y))
    aperture = CircularAperture(positions, r=r_pix)
    annulus = CircularAnnulus(positions, r_in=r_in, r_out=r_out)

    phot_table = aperture_photometry(data, aperture, mask=mask, method="subpixel", subpixels=5)
    ann_stats = ApertureStats(data, annulus, mask=mask)
    # Zeropoint
    zp_keys = ["MAGZP", "MAGZERO", "MAGZERO_RP", "PHOTZP"]
    zp = None
    for k in zp_keys:
        if k in hdr:
            zp = float(hdr[k])
            break

    mjd = None
    if "MJD-OBS" in hdr:
        mjd = float(hdr["MJD-OBS"])
    elif "DATE-OBS" in hdr:
        mjd = Time(hdr["DATE-OBS"]).mjd

    filt = hdr.get("FILTER", "").strip()
    if not filt:
        # PS1 warps commonly store filter in a hierarch-style keyword.
        filt = str(hdr.get("HIERARCH FPA.FILTER", hdr.get("FPA.FILTER", ""))).strip()
    airmass = hdr.get("AIRMASS", np.nan)
    exptime = hdr.get("EXPTIME", np.nan)
    obstag = hdr.get("OBSTYPE", "").strip()

    area = float(getattr(aperture, "area", np.pi * r_pix * r_pix))

    rows = []
    for idx, img in enumerate(labels):
        aper_sum = float(phot_table["aperture_sum"][idx])

        # Local background estimate from annulus
        bkg_med = float(ann_stats.median[idx]) if np.isfinite(ann_stats.median[idx]) else np.nan
        bkg_std = float(ann_stats.std[idx]) if np.isfinite(ann_stats.std[idx]) else np.nan
        if not np.isfinite(bkg_std):
            bkg_std = float(bkg_rms_global)

        # If the annulus is fully masked, fall back to a global sky estimate.
        if not np.isfinite(bkg_med):
            bkg_med = global_med

        if np.isfinite(aper_sum) and np.isfinite(bkg_med):
            flux = aper_sum - bkg_med * area
        else:
            flux = np.nan

        # crude Poisson + background variance estimate
        if np.isfinite(aper_sum) and np.isfinite(bkg_std):
            flux_err = float(np.sqrt(np.abs(aper_sum) + (bkg_std * np.sqrt(area)) ** 2 + 1e-8))
        else:
            flux_err = np.nan

        if zp is not None and np.isfinite(flux) and np.isfinite(flux_err) and flux > 0:
            mag = -2.5 * np.log10(flux) + zp
            mag_err = 1.0857 * flux_err / flux
        else:
            mag = np.nan
            mag_err = np.nan
        rows.append(
            {
                "band": filt,
                "image": img,
                "mjd": mjd,
                "flux": flux,
                "flux_err": flux_err,
                "mag": mag,
                "mag_err": mag_err,
                "zp": zp if zp is not None else np.nan,
                "airmass": airmass,
                "exptime": exptime,
                "obstype": obstag,
                "fname": fname.name,
            }
        )
    return rows


def process_all(fits_files: List[Path], priors: Dict[str, Dict[str, float]], outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    per_band: Dict[str, List[Dict]] = {}
    for fname in fits_files:
        try:
            rows = photometer_file(fname, priors)
            for r in rows:
                band = r["band"]
                per_band.setdefault(band, []).append(r)
        except Exception as e:
            log.warning(f"Failed on {fname}: {e}")

    # write CSVs
    import csv

    for band, rows in per_band.items():
        if not band:
            band = "unknown"
        out_csv = outdir / f"lightcurves_pg1115_ps1_{band}.csv"
        fieldnames = [
            "band",
            "image",
            "mjd",
            "flux",
            "flux_err",
            "mag",
            "mag_err",
            "zp",
            "airmass",
            "exptime",
            "obstype",
            "fname",
        ]
        with open(out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in rows:
                writer.writerow(r)
        log.info(f"Wrote {len(rows)} rows -> {out_csv}")


def main():
    parser = argparse.ArgumentParser(description="Self-extract PG1115 PS1 multi-band light curves.")
    parser.add_argument("--ra", type=float, default=169.57035, help="ICRS RA (deg) of lens center (PG1115)")
    parser.add_argument("--dec", type=float, default=7.76627, help="ICRS Dec (deg) of lens center (PG1115)")
    parser.add_argument("--radius", type=float, default=4.0, help="Search radius (arcmin)")
    parser.add_argument("--filters", type=str, default="g,r,i,z,y", help="Comma-separated PS1 filters")
    parser.add_argument("--max-products", type=int, default=120, help="Max warp products to download")
    parser.add_argument("--priors", type=Path, default=Path("data/priors/pg1115_hst_priors.json"), help="JSON priors path")
    parser.add_argument(
        "--hst-cutout-url",
        type=str,
        default="https://www.ast.cam.ac.uk/ioa/research/lensedquasars/HST_data/PG1115+080_ibgw22010_FILTER_F218W_EXPTIME2608.0_drz_cutout.fits",
        help="Optional HST cutout URL used to auto-derive priors if priors JSON is missing.",
    )
    parser.add_argument(
        "--hst-cutout-local",
        type=Path,
        default=None,
        help="Optional local HST FITS path (overrides URL download) used to auto-derive priors.",
    )
    parser.add_argument("--outdir", type=Path, default=Path("data/selfextract/pg1115_ps1"), help="Output directory")
    args = parser.parse_args()

    filt_list = [x.strip() for x in args.filters.split(",") if x.strip()]

    log.info("Querying PS1...")
    products = search_ps1(args.ra, args.dec, args.radius, filt_list, args.max_products)

    log.info("Downloading products...")
    raw_dir = args.outdir / "raw"
    fits_files = download_products(products, raw_dir)

    # Priors handling:
    # - If priors file exists and is sane, use it.
    # - Otherwise derive from an HST image with celestial WCS (MAST-backed).
    priors = None
    if args.priors.exists():
        try:
            priors_loaded = load_priors(args.priors)
            if priors_are_sane(priors_loaded):
                priors = priors_loaded
            else:
                log.warning("Existing priors file is not sane; will regenerate from HST: %s", args.priors)
        except Exception as e:
            log.warning("Failed to read existing priors file; will regenerate from HST: %s (%s)", args.priors, e)

    if priors is None:
        priors = derive_priors_for_pg1115(args, args.outdir)
        if not priors_are_sane(priors):
            raise RuntimeError("Derived priors were not sane (RA/Dec bounds).")
        save_priors(args.priors, priors)
        log.info(f"Wrote HST-derived priors -> {args.priors}")

    log.info("Photometering downloaded FITS...")
    process_all(fits_files, priors, args.outdir)

    log.info("Done.")


if __name__ == "__main__":
    main()
