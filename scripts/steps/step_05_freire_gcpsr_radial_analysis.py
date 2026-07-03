#!/usr/bin/env python3

import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import time
from urllib.request import urlopen, HTTPError, URLError
from socket import timeout as SocketTimeout

import numpy as np
from scipy import stats


FREIRE_GCPSR_URL = "https://www3.mpifr-bonn.mpg.de/staff/pfreire/GCpsr.txt"

RESULTS_DIR = Path(__file__).resolve().parents[2] / "results" / "outputs"
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
RAW_OUT_PATH = DATA_DIR / "freire_GCpsr.txt"
SUMMARY_CSV_PATH = RESULTS_DIR / "step_05_freire_gcpsr_radial_summary.csv"
SUMMARY_JSON_PATH = RESULTS_DIR / "step_05_freire_gcpsr_radial_summary.json"
SUMMARY_MD_PATH = RESULTS_DIR / "step_05_freire_gcpsr_radial_summary.md"


_NUM_RE = re.compile(r"^([+-]?(?:\d+\.?\d*|\d*\.?\d+))(?:\(.*\))?$")


def _parse_num(tok: str):
    tok = tok.strip()
    if tok in ("*", "i", ""):
        return None
    m = _NUM_RE.match(tok)
    if not m:
        return None
    return float(m.group(1))


def download_freire_catalog(max_retries: int = 3, timeout: int = 60):
    """Download Freire GC PSR catalog with retry logic and error handling."""
    for attempt in range(max_retries):
        try:
            print(f"Downloading Freire GC PSR catalog (attempt {attempt + 1}/{max_retries})...")
            raw = urlopen(FREIRE_GCPSR_URL, timeout=timeout).read()
            print(f"Downloaded {len(raw)} bytes")
            sha256 = hashlib.sha256(raw).hexdigest()
            RESULTS_DIR.mkdir(parents=True, exist_ok=True)
            RAW_OUT_PATH.write_bytes(raw)
            return raw.decode("utf-8", errors="replace"), sha256, len(raw)
        except HTTPError as e:
            print(f"HTTP Error {e.code}: {e.reason}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"Retrying in {wait_time}s...")
                time.sleep(wait_time)
            else:
                raise RuntimeError(f"Failed to download Freire catalog: HTTP {e.code}") from e
        except (URLError, SocketTimeout) as e:
            print(f"Network error: {e}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"Retrying in {wait_time}s...")
                time.sleep(wait_time)
            else:
                raise RuntimeError(f"Failed to download Freire catalog: Network error") from e
        except Exception as e:
            print(f"Unexpected error: {e}")
            raise RuntimeError(f"Failed to download Freire catalog: {e}") from e
    
    raise RuntimeError(f"Failed to download Freire catalog after {max_retries} attempts")


def parse_freire_gcpsr(text: str):
    """Parse Freire GCpsr.txt.

    Column meanings per header line:
    - r: projected offset from cluster center in arcminutes
    - P: period in milliseconds
    - Pdot: in units of 1e-20 (s/s)

    Returns list of dict rows with keys:
    cluster, pulsar, r_arcmin, p_ms, pdot_1e20
    """

    rows = []
    cluster = None

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if not line.strip():
            continue

        if "\t" not in line and not line.startswith(("J", "B")):
            cluster = line.strip()
            continue

        if not line.startswith(("J", "B")):
            continue

        parts = [p for p in line.split("\t") if p != ""]
        if len(parts) < 4:
            parts = line.split()
        if len(parts) < 4:
            continue

        pulsar = parts[0].strip()
        r_arcmin = _parse_num(parts[1])
        p_ms = _parse_num(parts[2])
        pdot_1e20 = _parse_num(parts[3])

        rows.append(
            {
                "cluster": cluster,
                "pulsar": pulsar,
                "r_arcmin": r_arcmin,
                "p_ms": p_ms,
                "pdot_1e20": pdot_1e20,
            }
        )

    return rows


def compute_cluster_correlations(rows, min_n=5):
    clusters = defaultdict(list)
    for r in rows:
        if r["cluster"] is None:
            continue
        if r["r_arcmin"] is None or r["pdot_1e20"] is None:
            continue
        clusters[r["cluster"]].append(r)

    results = []
    for cluster, rs in clusters.items():
        if len(rs) < min_n:
            continue

        r_arcmin = np.array([x["r_arcmin"] for x in rs], dtype=float)
        pdot_abs = np.array([abs(x["pdot_1e20"]) for x in rs], dtype=float)
        pdot_signed = np.array([x["pdot_1e20"] for x in rs], dtype=float)

        mask_abs = pdot_abs > 0
        r_arcmin_abs = r_arcmin[mask_abs]
        log10_abs = np.log10(pdot_abs[mask_abs])

        if len(r_arcmin_abs) < min_n:
            continue

        pear_r_abs, pear_p_abs = stats.pearsonr(r_arcmin_abs, log10_abs)
        spear_rho_abs, spear_p_abs = stats.spearmanr(r_arcmin_abs, log10_abs)

        mask_signed = pdot_signed != 0
        r_arcmin_signed = r_arcmin[mask_signed]
        pdot_signed_use = pdot_signed[mask_signed]

        pear_r_signed = None
        pear_p_signed = None
        spear_rho_signed = None
        spear_p_signed = None
        if len(r_arcmin_signed) >= min_n:
            pear_r_signed, pear_p_signed = stats.pearsonr(r_arcmin_signed, pdot_signed_use)
            spear_rho_signed, spear_p_signed = stats.spearmanr(r_arcmin_signed, pdot_signed_use)

        span_arcmin = float(r_arcmin_abs.max() - r_arcmin_abs.min())
        span_arcsec = span_arcmin * 60.0

        results.append(
            {
                "cluster": cluster,
                "n": int(len(r_arcmin_abs)),
                "r_span_arcmin": span_arcmin,
                "r_span_arcsec": span_arcsec,
                "pearson_r_logabs": float(pear_r_abs),
                "pearson_p_logabs": float(pear_p_abs),
                "spearman_rho_logabs": float(spear_rho_abs),
                "spearman_p_logabs": float(spear_p_abs),
                "pearson_r_signed": None if pear_r_signed is None else float(pear_r_signed),
                "pearson_p_signed": None if pear_p_signed is None else float(pear_p_signed),
                "spearman_rho_signed": None if spear_rho_signed is None else float(spear_rho_signed),
                "spearman_p_signed": None if spear_p_signed is None else float(spear_p_signed),
            }
        )

    results.sort(key=lambda x: (-x["n"], -x["r_span_arcsec"]))
    return results


def write_outputs(meta, results):
    with SUMMARY_CSV_PATH.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "cluster",
                "n",
                "r_span_arcsec",
                "pearson_r_logabs",
                "pearson_p_logabs",
                "spearman_rho_logabs",
                "spearman_p_logabs",
                "pearson_r_signed",
                "pearson_p_signed",
                "spearman_rho_signed",
                "spearman_p_signed",
            ]
        )
        for r in results:
            w.writerow(
                [
                    r["cluster"],
                    r["n"],
                    f"{r['r_span_arcsec']:.3f}",
                    f"{r['pearson_r_logabs']:.6f}",
                    f"{r['pearson_p_logabs']:.6g}",
                    f"{r['spearman_rho_logabs']:.6f}",
                    f"{r['spearman_p_logabs']:.6g}",
                    "" if r["pearson_r_signed"] is None else f"{r['pearson_r_signed']:.6f}",
                    "" if r["pearson_p_signed"] is None else f"{r['pearson_p_signed']:.6g}",
                    "" if r["spearman_rho_signed"] is None else f"{r['spearman_rho_signed']:.6f}",
                    "" if r["spearman_p_signed"] is None else f"{r['spearman_p_signed']:.6g}",
                ]
            )

    out = {"meta": meta, "results": results}
    SUMMARY_JSON_PATH.write_text(json.dumps(out, indent=2))

    lines = []
    lines.append(f"# Freire GCpsr Radial Correlation Summary\n")
    lines.append(f"**Source URL:** {meta['source_url']}\n")
    lines.append(f"**Downloaded:** {meta['downloaded_at']}\n")
    lines.append(f"**SHA256:** `{meta['sha256']}`\n")
    lines.append(f"**Bytes:** {meta['bytes']}\n")
    lines.append("\n")
    lines.append(
        "This uses the Freire catalog columns exactly as documented: `r` is projected offset in **arcminutes** and `Pdot` is in units of `1e-20 s/s`.\n"
    )
    lines.append("We report correlations of `r` vs `log10(|Pdot|)` (magnitude) and also `r` vs `Pdot` (signed).\n")
    lines.append("\n")
    lines.append("| Cluster | N | span(\") | Pearson r (log|Pdot|) | p | Spearman rho | p |\n")
    lines.append("|---|---:|---:|---:|---:|---:|---:|\n")
    for r in results:
        lines.append(
            f"| {r['cluster']} | {r['n']} | {r['r_span_arcsec']:.1f} | {r['pearson_r_logabs']:+.3f} | {r['pearson_p_logabs']:.4g} | {r['spearman_rho_logabs']:+.3f} | {r['spearman_p_logabs']:.4g} |\n"
        )

    SUMMARY_MD_PATH.write_text("".join(lines))


def main():
    """Download and analyze Freire GC pulsar catalog for radial correlations.
    
    Downloads the Freire GCPSR catalog, parses pulsar data, computes
    radial correlations for each cluster, and writes output files.
    """
    text, sha256, nbytes = download_freire_catalog()
    rows = parse_freire_gcpsr(text)
    results = compute_cluster_correlations(rows, min_n=5)

    meta = {
        "source_url": FREIRE_GCPSR_URL,
        "downloaded_at": datetime.now(timezone.utc).isoformat() + "Z",
        "sha256": sha256,
        "bytes": nbytes,
        "min_n": 5,
        "analysis": {
            "x": "r_arcmin (projected offset)",
            "y1": "log10(|Pdot_1e20|)",
            "y2": "Pdot_1e20 (signed)",
        },
    }

    write_outputs(meta, results)

    print(f"Wrote: {RAW_OUT_PATH}")
    print(f"Wrote: {SUMMARY_CSV_PATH}")
    print(f"Wrote: {SUMMARY_JSON_PATH}")
    print(f"Wrote: {SUMMARY_MD_PATH}")


if __name__ == "__main__":
    main()
