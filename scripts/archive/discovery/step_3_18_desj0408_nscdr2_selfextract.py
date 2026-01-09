#!/usr/bin/env python3
"""Step 3.18: DESJ0408 multi-epoch light curves via NOIRLab Data Lab NSC DR2

This step builds multi-band, multi-epoch light curves for DESJ0408-5354 by querying
NOIRLab Data Lab's NOIRLab Source Catalog (NSC) DR2 measurement-level tables via TAP.

Rationale
- DES single-epoch calibrated images are publicly described, but end-to-end image
  retrieval + forced photometry is operationally heavy.
- NSC DR2 already provides per-exposure photometry and timestamps across DECam
  and other instruments, accessible through the public Data Lab TAP service.
- This script is a reproducible first-pass multi-band, multi-epoch light-curve
  extractor for chromaticity tests.

Outputs
- data/selfextract/desj0408_nscdr2/lightcurves_desj0408_nscdr2_<band>.csv
  columns: band,image,mjd,mag,mag_err,ra,dec,sep_arcsec

Notes
- Astrometric priors are derived from an HST drizzled image with celestial WCS
  downloaded from MAST. Images are labeled A/B/C/D by deterministic sorting.
- Each NSC measurement is assigned to the nearest image prior if within
  --max-sep-arcsec.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u

from astroquery.mast import Observations

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("desj0408_nscdr2")


def read_first_image_hdu(fname: Path) -> Tuple[np.ndarray, fits.Header]:
    with fits.open(fname) as hdul:
        for hdu in hdul:
            if hdu.data is None:
                continue
            data = np.asarray(hdu.data)
            if data.ndim == 2:
                return data.astype(float), hdu.header
    raise RuntimeError("No 2D image HDU found.")


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
    p16, p84 = np.percentile(vals, [16.0, 84.0])
    sigma = 0.5 * (p84 - p16)
    if np.isfinite(sigma) and sigma > 0:
        return float(sigma)
    lo, hi = np.percentile(vals, [5.0, 95.0])
    clipped = vals[(vals >= lo) & (vals <= hi)]
    if clipped.size < 50:
        return float("nan")
    sigma = float(np.std(clipped))
    return sigma if np.isfinite(sigma) and sigma > 0 else float("nan")


def download_hst_drz_from_mast(ra: float, dec: float, radius_arcsec: float, outdir: Path) -> Path:
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

    if "productFilename" in products.colnames:
        fn = np.array([str(x).lower() for x in products["productFilename"]], dtype=str)
        is_fits = np.char.endswith(fn, ".fits")
        is_drz = np.array([("_drz" in x) or ("_drc" in x) for x in fn], dtype=bool)
        products = products[is_fits & is_drz]

    if products is None or len(products) == 0:
        raise RuntimeError("No HST drizzled FITS products (_drz/_drc) found in MAST products.")

    # Download just one
    products = products[:1]
    manifest = Observations.download_products(
        products,
        mrp_only=False,
        download_dir=str(outdir),
        cache=True,
        curl_flag=False,
    )
    if manifest is None or len(manifest) == 0 or "Local Path" not in manifest.colnames:
        raise RuntimeError("HST MAST download produced an empty manifest.")

    for p in manifest["Local Path"]:
        if p and str(p).lower().endswith(".fits"):
            return Path(str(p))
    raise RuntimeError("HST MAST download did not yield any FITS local paths.")


def derive_priors_from_hst(
    hst_fits: Path,
    ra_center: float,
    dec_center: float,
    box_half_size_pix: int = 160,
    nsigma: float = 10.0,
    fwhm_pix: float = 2.5,
) -> Dict[str, Dict[str, float]]:
    from photutils.detection import DAOStarFinder

    data, hdr = read_first_image_hdu(hst_fits)
    wcs = WCS(hdr)
    if not wcs.has_celestial:
        raise RuntimeError("HST image does not contain a celestial WCS; cannot derive priors.")

    # Restrict detection to a cutout around the lens center to avoid selecting unrelated field sources.
    x0, y0 = wcs.all_world2pix([ra_center], [dec_center], 0)
    x0, y0 = float(x0[0]), float(y0[0])
    if not np.isfinite(x0) or not np.isfinite(y0):
        raise RuntimeError("Lens center could not be projected into HST pixel coordinates.")

    ny, nx = data.shape
    if x0 < 0 or x0 >= nx or y0 < 0 or y0 >= ny:
        raise RuntimeError("Lens center projects outside HST image bounds.")
    x0i, y0i = int(round(x0)), int(round(y0))
    x1 = max(0, x0i - box_half_size_pix)
    x2 = min(nx, x0i + box_half_size_pix)
    y1 = max(0, y0i - box_half_size_pix)
    y2 = min(ny, y0i + box_half_size_pix)
    cut = data[y1:y2, x1:x2]
    if cut.size == 0:
        raise RuntimeError("Empty HST cutout.")

    bkg = np.nanmedian(cut)
    img = cut - bkg
    rms = robust_bkg_rms(img)
    if not np.isfinite(rms) or rms <= 0:
        raise RuntimeError("Could not estimate background RMS from HST image.")

    finder = DAOStarFinder(fwhm=fwhm_pix, threshold=nsigma * rms)
    sources = finder(img)
    if sources is None or len(sources) < 4:
        raise RuntimeError(f"DAOStarFinder found {0 if sources is None else len(sources)} sources; need >=4")

    # Select the 4 sources closest to the lens center (not simply the brightest).
    xs_all = np.array(sources["xcentroid"], dtype=float) + x1
    ys_all = np.array(sources["ycentroid"], dtype=float) + y1
    ras_all, decs_all = wcs.all_pix2world(xs_all, ys_all, 0)
    sc_all = SkyCoord(ra=ras_all * u.deg, dec=decs_all * u.deg)
    center = SkyCoord(ra_center * u.deg, dec_center * u.deg)
    seps = sc_all.separation(center).to(u.arcsec).value
    order = np.argsort(seps)
    order = order[:4]
    det = [(float(ras_all[i]), float(decs_all[i])) for i in order]

    # Deterministic label: sort by Dec desc then RA asc
    det.sort(key=lambda t: (-t[1], t[0]))
    labels = ["A", "B", "C", "D"]
    priors = {lab: {"ra": ra, "dec": dec} for lab, (ra, dec) in zip(labels, det)}
    return priors


def priors_near_center(priors: Dict[str, Dict[str, float]], ra_center: float, dec_center: float, max_sep_arcsec: float = 10.0) -> bool:
    try:
        center = SkyCoord(ra_center * u.deg, dec_center * u.deg)
        for k in ["A", "B", "C", "D"]:
            sc = SkyCoord(float(priors[k]["ra"]) * u.deg, float(priors[k]["dec"]) * u.deg)
            if center.separation(sc).to(u.arcsec).value > max_sep_arcsec:
                return False
        return True
    except Exception:
        return False


def discover_nsc_meas_table_and_cols(tap_url: str = "https://datalab.noirlab.edu/tap"):
    from pyvo.dal import TAPService

    svc = TAPService(tap_url)

    # List candidate tables
    q_tables = """
    SELECT table_name
    FROM TAP_SCHEMA.tables
    WHERE table_name LIKE 'nsc_dr2.%'
    """.strip()
    tables = svc.search(q_tables).to_table()["table_name"].tolist()

    candidates = [t for t in tables if ("meas" in t.lower()) or ("measurement" in t.lower())]
    # Put the most likely first
    candidates.sort(key=lambda s: (0 if "meas" in s.lower() else 1, len(s)))
    if not candidates:
        candidates = tables

    # Column name preferences
    ra_names = ["ra", "ra_meas", "ra_psf", "ra_icrs"]
    dec_names = ["dec", "dec_meas", "dec_psf", "dec_icrs"]
    time_names = ["mjd", "mjd_obs", "mjd_mid", "jd"]
    band_names = ["filter", "filt", "band", "f"]
    mag_names = ["mag", "mag_psf", "mag_auto", "psfmag", "cmag"]
    magerr_names = [
        "magerr_auto",
        "magerr",
        "e_mag",
        "mag_err",
        "e_mag_psf",
        "psfmagerr",
    ]

    def get_cols(table_name: str) -> List[str]:
        q_cols = f"SELECT column_name FROM TAP_SCHEMA.columns WHERE table_name='{table_name}'"
        return svc.search(q_cols).to_table()["column_name"].tolist()

    def pick(cols: List[str], options: List[str]) -> Optional[str]:
        cols_l = {c.lower(): c for c in cols}
        for o in options:
            if o.lower() in cols_l:
                return cols_l[o.lower()]
        return None

    best = None
    best_score = -1
    best_map = None

    for t in candidates[:25]:
        cols = get_cols(t)
        ra = pick(cols, ra_names)
        dec = pick(cols, dec_names)
        tim = pick(cols, time_names)
        band = pick(cols, band_names)
        mag = pick(cols, mag_names)
        magerr = pick(cols, magerr_names)

        score = sum([ra is not None, dec is not None, tim is not None, band is not None, mag is not None, magerr is not None])
        if score > best_score:
            best_score = score
            best = t
            best_map = {"ra": ra, "dec": dec, "time": tim, "band": band, "mag": mag, "magerr": magerr}
        if score == 6:
            break

    if best is None or best_map is None or best_score < 5:
        raise RuntimeError(f"Could not identify a usable NSC DR2 measurement table (best={best}, score={best_score}).")

    if best_map.get("magerr") is None:
        raise RuntimeError(f"Could not identify magnitude error column for {best}. Found map={best_map}")

    return svc, best, best_map


def query_nsc_measurements(
    ra: float,
    dec: float,
    radius_arcsec: float,
    max_rows: int,
    bands: List[str],
    tap_url: str,
):
    svc, table, col = discover_nsc_meas_table_and_cols(tap_url)

    radius_deg = radius_arcsec / 3600.0
    # RA half-width adjusted by cos(dec) so the box fully contains the cone.
    dra = radius_deg / max(1e-6, np.cos(np.deg2rad(dec)))
    ra_min, ra_max = ra - dra, ra + dra
    dec_min, dec_max = dec - radius_deg, dec + radius_deg
    band_list = ",".join([f"'{b}'" for b in bands])

    # Compose ADQL. Prefer q3c_radial_query if supported; otherwise fall back to a bounding box.
    q_q3c = f"""
    SELECT TOP {int(max_rows)}
      objectid AS objectid,
      {col['ra']} AS ra,
      {col['dec']} AS dec,
      {col['time']} AS mjd,
      {col['band']} AS band,
      {col['mag']} AS mag,
      {col['magerr']} AS mag_err
    FROM {table}
    WHERE q3c_radial_query({col['ra']}, {col['dec']}, {ra}, {dec}, {radius_deg})
      AND {col['band']} IN ({band_list})
    """.strip()

    q_box = f"""
    SELECT TOP {int(max_rows)}
      objectid AS objectid,
      {col['ra']} AS ra,
      {col['dec']} AS dec,
      {col['time']} AS mjd,
      {col['band']} AS band,
      {col['mag']} AS mag,
      {col['magerr']} AS mag_err
    FROM {table}
    WHERE {col['ra']} BETWEEN {ra_min} AND {ra_max}
      AND {col['dec']} BETWEEN {dec_min} AND {dec_max}
      AND {col['band']} IN ({band_list})
    """.strip()

    log.info(f"NSC query table: {table}")
    log.info(f"NSC query columns: {col}")

    try:
        tab = svc.search(q_q3c).to_table()
        return tab
    except Exception as e:
        log.warning(f"q3c_radial_query failed; falling back to RA/Dec bounding box: {e}")
        tab = svc.search(q_box).to_table()
        return tab


def assign_to_images(tab, priors: Dict[str, Dict[str, float]], max_sep_arcsec: float):
    prior_labels = ["A", "B", "C", "D"]
    pri = SkyCoord(
        ra=[priors[k]["ra"] for k in prior_labels] * u.deg,
        dec=[priors[k]["dec"] for k in prior_labels] * u.deg,
    )

    # TAP/pyvo sometimes attaches non-standard unit metadata to columns; cast to float arrays.
    ra_det = np.array(tab["ra"], dtype=float)
    dec_det = np.array(tab["dec"], dtype=float)
    det = SkyCoord(ra=ra_det * u.deg, dec=dec_det * u.deg)
    idx, sep2d, _ = det.match_to_catalog_sky(pri)

    labels = np.array([prior_labels[i] for i in idx], dtype=object)
    sep_arcsec = sep2d.to(u.arcsec).value
    ok = sep_arcsec <= max_sep_arcsec

    return labels, sep_arcsec, ok


def pick_objectids_by_image(
    tab,
    labels: np.ndarray,
    sep_arcsec: np.ndarray,
    ok_mask: np.ndarray,
    min_count: int = 3,
) -> Dict[str, str]:
    """Pick a single NSC objectid per image label.

    Strategy: among rows assigned to a given label, select the objectid with the
    smallest median separation to the corresponding prior (and require >= min_count
    measurements). This is more robust than picking the mode when the assignment
    radius is permissive and nearby sources contaminate the label.
    """
    if "objectid" not in tab.colnames:
        raise RuntimeError("NSC table does not include objectid; cannot refine.")

    objectids = np.array([str(x) for x in tab["objectid"]], dtype=object)
    out: Dict[str, str] = {}
    for lab in ["A", "B", "C", "D"]:
        m = ok_mask & (labels == lab)
        if not np.any(m):
            continue
        vals = objectids[m]
        seps = np.array(sep_arcsec[m], dtype=float)

        uniq = np.unique(vals)
        best_id = None
        best_med = None
        best_n = 0
        for oid in uniq:
            mm = vals == oid
            if int(mm.sum()) < min_count:
                continue
            med = float(np.nanmedian(seps[mm]))
            if not np.isfinite(med):
                continue
            if best_med is None or med < best_med or (med == best_med and int(mm.sum()) > best_n):
                best_med = med
                best_id = str(oid)
                best_n = int(mm.sum())

        if best_id is not None:
            out[lab] = best_id
    return out


def query_nsc_measurements_for_objectids(
    svc,
    table: str,
    col: Dict[str, str],
    objectids: Dict[str, int],
    bands: List[str],
    max_rows: int,
):
    # Pull all measurements for selected objectids.
    obj_list = ",".join([f"'{str(v)}'" for v in sorted(set(objectids.values()))])
    band_list = ",".join([f"'{b}'" for b in bands])
    q = f"""
    SELECT TOP {int(max_rows)}
      objectid AS objectid,
      {col['ra']} AS ra,
      {col['dec']} AS dec,
      {col['time']} AS mjd,
      {col['band']} AS band,
      {col['mag']} AS mag,
      {col['magerr']} AS mag_err
    FROM {table}
    WHERE objectid IN ({obj_list})
      AND {col['band']} IN ({band_list})
    """.strip()
    return svc.search(q).to_table()


def write_lightcurves(outdir: Path, system_tag: str, tab, labels, sep_arcsec, ok_mask):
    outdir.mkdir(parents=True, exist_ok=True)

    rows_by_band: Dict[str, List[Dict]] = {}
    for i in range(len(tab)):
        if not bool(ok_mask[i]):
            continue
        band = str(tab["band"][i]).strip()
        rows_by_band.setdefault(band, []).append(
            {
                "band": band,
                "image": str(labels[i]),
                "mjd": float(tab["mjd"][i]) if np.isfinite(float(tab["mjd"][i])) else np.nan,
                "mag": float(tab["mag"][i]) if np.isfinite(float(tab["mag"][i])) else np.nan,
                "mag_err": float(tab["mag_err"][i]) if np.isfinite(float(tab["mag_err"][i])) else np.nan,
                "ra": float(tab["ra"][i]),
                "dec": float(tab["dec"][i]),
                "sep_arcsec": float(sep_arcsec[i]),
                "objectid": str(tab["objectid"][i]) if "objectid" in tab.colnames else "",
            }
        )

    for band, rows in rows_by_band.items():
        out_csv = outdir / f"lightcurves_{system_tag}_{band}.csv"
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(
                f,
                fieldnames=["band", "image", "mjd", "mag", "mag_err", "ra", "dec", "sep_arcsec", "objectid"],
            )
            w.writeheader()
            for r in rows:
                w.writerow(r)
        log.info(f"Wrote {len(rows)} rows -> {out_csv}")


def main():
    parser = argparse.ArgumentParser(description="DESJ0408 NSC DR2 multi-epoch light-curve extraction")
    parser.add_argument("--ra", type=float, default=62.0905, help="ICRS RA deg")
    parser.add_argument("--dec", type=float, default=-53.8999, help="ICRS Dec deg")
    parser.add_argument("--radius-arcsec", type=float, default=10.0, help="Cone search radius")
    parser.add_argument("--max-rows", type=int, default=50000, help="Max rows to fetch")
    parser.add_argument("--bands", type=str, default="g,r,i,z,Y", help="Comma-separated bands")
    parser.add_argument("--max-sep-arcsec", type=float, default=2.0, help="Max separation to assign to an image prior")
    parser.add_argument("--tap-url", type=str, default="https://datalab.noirlab.edu/tap")
    parser.add_argument(
        "--refine-by-objectid",
        action="store_true",
        help="Refine assignments by selecting one NSC objectid per image, then re-querying all meas for those objectids.",
    )
    parser.add_argument("--outdir", type=Path, default=Path("data/selfextract/desj0408_nscdr2"))
    parser.add_argument("--priors-json", type=Path, default=Path("data/priors/desj0408_hst_priors.json"))
    args = parser.parse_args()

    # Priors
    priors = None
    if args.priors_json.exists():
        try:
            priors_loaded = json.loads(args.priors_json.read_text())
            if priors_near_center(priors_loaded, args.ra, args.dec, max_sep_arcsec=10.0):
                priors = priors_loaded
            else:
                log.warning("Existing priors are far from lens center; regenerating: %s", args.priors_json)
        except Exception as e:
            log.warning("Failed to read existing priors; regenerating: %s (%s)", args.priors_json, e)

    if priors is None:
        hst_dir = args.outdir / "priors"
        hst_fits = download_hst_drz_from_mast(args.ra, args.dec, radius_arcsec=60.0, outdir=hst_dir)
        priors = derive_priors_from_hst(hst_fits, ra_center=args.ra, dec_center=args.dec)
        args.priors_json.parent.mkdir(parents=True, exist_ok=True)
        args.priors_json.write_text(json.dumps(priors, indent=2, sort_keys=True) + "\n")
        log.info(f"Wrote priors -> {args.priors_json}")

    bands = [b.strip() for b in args.bands.split(",") if b.strip()]

    tab = query_nsc_measurements(
        ra=args.ra,
        dec=args.dec,
        radius_arcsec=args.radius_arcsec,
        max_rows=args.max_rows,
        bands=bands,
        tap_url=args.tap_url,
    )

    log.info(f"NSC rows returned: {len(tab)}")
    if len(tab) == 0:
        raise RuntimeError("NSC query returned zero rows.")

    labels, sep_arcsec, ok = assign_to_images(tab, priors, args.max_sep_arcsec)
    n_assigned = int(np.sum(ok))
    log.info(f"Assigned to priors within {args.max_sep_arcsec:.2f} arcsec: {n_assigned} / {len(tab)}")
    if n_assigned == 0:
        raise RuntimeError(
            "No NSC measurements were within max separation of the HST-derived priors. "
            "Try increasing --max-sep-arcsec or inspect data/priors/desj0408_hst_priors.json."
        )
    if args.refine_by_objectid:
        svc, table, col = discover_nsc_meas_table_and_cols(args.tap_url)
        obj_map = pick_objectids_by_image(tab, labels, sep_arcsec, ok, min_count=3)
        log.info(f"Selected objectids per image: {obj_map}")
        if len(obj_map) < 2:
            raise RuntimeError("ObjectID refinement found <2 images with objectids; cannot build pairwise light curves.")
        tab2 = query_nsc_measurements_for_objectids(svc, table, col, obj_map, bands=bands, max_rows=args.max_rows)
        log.info(f"Re-query by objectid returned rows: {len(tab2)}")

        # Label rows by objectid
        objectid_arr = np.array([str(x) for x in tab2["objectid"]], dtype=object)
        inv = {str(v): k for k, v in obj_map.items()}
        labels2 = np.array([inv.get(str(oid), "") for oid in objectid_arr], dtype=object)
        ok2 = labels2 != ""

        # Separation to priors for QC
        pri_sc = SkyCoord(
            [priors[k]["ra"] for k in ["A", "B", "C", "D"]] * u.deg,
            [priors[k]["dec"] for k in ["A", "B", "C", "D"]] * u.deg,
        )
        det_sc = SkyCoord(np.array(tab2["ra"], dtype=float) * u.deg, np.array(tab2["dec"], dtype=float) * u.deg)
        # Map label->index
        lab_to_idx = {"A": 0, "B": 1, "C": 2, "D": 3}
        sep2 = np.full(len(tab2), np.nan, dtype=float)
        for i in range(len(tab2)):
            lab = labels2[i]
            if lab in lab_to_idx:
                sep2[i] = det_sc[i].separation(pri_sc[lab_to_idx[lab]]).to(u.arcsec).value

        # Attach objectid to output rows
        write_lightcurves(args.outdir, "desj0408_nscdr2", tab2, labels2, sep2, ok2)
    else:
        # Attach objectid column for provenance
        write_lightcurves(args.outdir, "desj0408_nscdr2", tab, labels, sep_arcsec, ok)


if __name__ == "__main__":
    main()
