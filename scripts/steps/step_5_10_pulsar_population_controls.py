#!/usr/bin/env python3

import csv
import hashlib
import io
import json
import math
import os
import re
import tarfile
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import time
from urllib.request import urlopen, HTTPError, URLError
from socket import timeout as SocketTimeout

import numpy as np
from scipy import stats

# Constants
MSP_P0_CUT_SECONDS = 0.03  # MSP period threshold
BOOTSTRAP_ITERATIONS = 2000
RANDOM_SEED = 42


FREIRE_GCPSR_URL = "https://www3.mpifr-bonn.mpg.de/staff/pfreire/GCpsr.txt"
ATNF_PSRCAT_TGZ_URL = "https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz"

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "outputs"
DATA_DIR = REPO_ROOT / "data"

OUT_JSON = RESULTS_DIR / "step_5_10_pulsar_population_controls.json"
OUT_CSV = RESULTS_DIR / "step_5_10_pulsar_population_controls.csv"
OUT_MD = RESULTS_DIR / "step_5_10_pulsar_population_controls.md"

RAW_FREIRE_PATH = RESULTS_DIR / "freire_GCpsr.txt"
RAW_ATNF_TGZ_PATH = RESULTS_DIR / "atnf_psrcat_pkg.tar.gz"
RAW_ATNF_DB_PATH = RESULTS_DIR / "atnf_psrcat.db"


_NUM_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\d*\.?\d+)(?:[eE][+-]?\d+)?$")


def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _download_with_retry(url: str, max_retries: int = 3, timeout: int = 60) -> tuple[bytes, str]:
    """Download data with retry logic and error handling."""
    for attempt in range(max_retries):
        try:
            print(f"  Downloading (attempt {attempt + 1}/{max_retries}): {url[:60]}...")
            raw = urlopen(url, timeout=timeout).read()
            print(f"  Downloaded {len(raw)} bytes")
            return raw, _sha256_bytes(raw)
        except HTTPError as e:
            print(f"  HTTP Error {e.code}: {e.reason}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt  # Exponential backoff
                print(f"  Retrying in {wait_time}s...")
                time.sleep(wait_time)
            else:
                raise RuntimeError(f"Failed to download {url}: HTTP {e.code}") from e
        except (URLError, SocketTimeout) as e:
            print(f"  Network error: {e}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"  Retrying in {wait_time}s...")
                time.sleep(wait_time)
            else:
                raise RuntimeError(f"Failed to download {url}: Network error") from e
        except Exception as e:
            print(f"  Unexpected error: {e}")
            raise RuntimeError(f"Failed to download {url}: {e}") from e
    
    raise RuntimeError(f"Failed to download {url} after {max_retries} attempts")


def _parse_freire_gcpsr(text: str) -> list[dict]:
    """Parse Freire GCpsr.

    Columns per file header:
      Pulsar, r(arcmin), P(ms), Pdot(1e-20), DM, ...

    We extract:
      cluster, name, P_ms, P1_e20 (signed)
    """

    rows = []
    cluster = None

    for raw_line in text.splitlines():
        line = raw_line.strip("\n")
        if not line.strip():
            continue

        # Cluster header lines have no leading whitespace and don't start with pulsar names
        if not line.startswith((" ", "\t", "J", "B")) and not line.startswith("#"):
            cluster = line.strip()
            continue

        # Pulsar lines start with J or B (after stripping leading whitespace)
        if not line.lstrip().startswith(("J", "B")):
            continue
        
        # Use stripped line for parsing
        line_stripped = line.lstrip()
        parts = [p for p in line_stripped.split("\t") if p != ""]
        if len(parts) < 4:
            parts = line.split()
        if len(parts) < 4:
            continue

        name = parts[0].strip()
        p_ms = parts[2].strip()
        p1_e20 = parts[3].strip()

        if p_ms in ("*", "i") or p1_e20 in ("*", "i"):
            continue

        # strip uncertainties like -4.9850(6)
        def parse_num(tok: str):
            tok = tok.strip()
            if tok in ("*", "i", ""):
                return None
            tok = re.sub(r"\(.*\)$", "", tok)
            if not _NUM_RE.match(tok):
                return None
            return float(tok)

        p_ms_f = parse_num(p_ms)
        p1_e20_f = parse_num(p1_e20)
        if p_ms_f is None or p1_e20_f is None:
            continue

        rows.append(
            {
                "source": "freire",
                "environment": "globular_cluster",
                "cluster": cluster,
                "name": name,
                "P0_s": p_ms_f / 1000.0,
                "P_ms": p_ms_f,
                "P1_sps": p1_e20_f * 1e-20,
                "P1_e20": p1_e20_f,
                "assoc": cluster,
            }
        )

    return rows


def _extract_psrcat_db_from_tgz(tgz_bytes: bytes) -> bytes:
    with tarfile.open(fileobj=io.BytesIO(tgz_bytes), mode="r:gz") as tf:
        # Try common names
        for member in tf.getmembers():
            base = Path(member.name).name
            if base.lower() in {"psrcat.db", "psrcat.db.txt"}:
                f = tf.extractfile(member)
                if f is None:
                    continue
                return f.read()

        # Fallback: first file with .db extension
        for member in tf.getmembers():
            if member.name.lower().endswith(".db"):
                f = tf.extractfile(member)
                if f is None:
                    continue
                return f.read()

    raise RuntimeError("Could not locate psrcat.db inside psrcat_pkg.tar.gz")


def _parse_atnf_psrcat_db(db_text: str) -> list[dict]:
    """Parse ATNF psrcat.db (ASCII block format).

    We extract:
      PSRJ (or PSRB), P0, P1, ASSOC

    Records are separated by lines starting with '@'.
    """

    rows = []
    current = {}

    def flush():
        nonlocal current
        if not current:
            return

        name = current.get("PSRJ") or current.get("PSRB")

        # ATNF entries may provide P0/P1 directly, or provide F0/F1.
        # Convert F0/F1 -> P0/P1 if needed:
        #   P0 = 1/F0
        #   P1 = -F1 / F0^2
        p0 = current.get("P0")
        p1 = current.get("P1")
        f0 = current.get("F0")
        f1 = current.get("F1")
        if p0 is None and f0 is not None and f0 != 0:
            p0 = 1.0 / f0
        if p1 is None and f0 is not None and f1 is not None and f0 != 0:
            p1 = -f1 / (f0 * f0)

        if name is None or p0 is None or p1 is None:
            current = {}
            return

        rows.append(
            {
                "source": "atnf",
                "environment": None,
                "cluster": None,
                "name": name,
                "P0_s": p0,
                "P_ms": p0 * 1000.0,
                "P1_sps": p1,
                "P1_e20": p1 / 1e-20,
                "assoc": current.get("ASSOC", ""),
            }
        )
        current = {}

    for raw_line in db_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#"):
            continue
        if line.startswith("@"):  # record separator
            flush()
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        key = parts[0].strip()
        # ATNF format: KEY VALUE [uncertainty] [ref]
        # We parse the first token as VALUE for numeric keys.
        val0 = parts[1].strip()
        val_rest = " ".join(parts[1:]).strip()

        # Numeric parsing for P0 and P1.
        if key in {"P0", "P1", "F0", "F1"}:
            # Values can be like '1.234E-20' or '-0.14e-16'
            if _NUM_RE.match(val0):
                current[key] = float(val0)
            continue

        if key in {"PSRJ", "PSRB", "ASSOC"}:
            if key in {"PSRJ", "PSRB"}:
                current[key] = val0
            else:
                current[key] = val_rest

        if key == "TYPE":
            current[key] = val_rest

    flush()
    return rows


def _looks_like_gc_assoc(assoc: str, freire_cluster_tokens: set[str]) -> bool:
    if not assoc:
        return False
    a = assoc.lower()
    if "gc" in a:
        return True
    # Common globular cluster tokens
    for tok in freire_cluster_tokens:
        if tok and tok.lower() in a:
            return True
    return False


def _compute_proxies(rows: list[dict]) -> list[dict]:
    out = []
    for r in rows:
        p0 = r["P0_s"]
        p1 = r["P1_sps"]
        p1_abs = abs(p1)
        if p0 <= 0 or p1_abs <= 0:
            continue

        # Proxies commonly used (with caveat: for GC pulsars observed Pdot is contaminated by acceleration)
        tau_c_s = p0 / (2.0 * p1_abs)
        b_proxy = math.sqrt(p0 * p1_abs)  # proportional to Bsurf up to constant factor

        r2 = dict(r)
        r2["logP"] = float(math.log10(p0))
        r2["logPdot_abs"] = float(math.log10(p1_abs))
        r2["tau_c_s"] = float(tau_c_s)
        r2["log_tau_c"] = float(math.log10(tau_c_s))
        r2["b_proxy"] = float(b_proxy)
        r2["log_b_proxy"] = float(math.log10(b_proxy))
        out.append(r2)

    return out


def _ttest_logpdot(gc: np.ndarray, field: np.ndarray) -> dict:
    t_stat, p_value = stats.ttest_ind(gc, field, equal_var=False)
    mw_u, mw_p = stats.mannwhitneyu(gc, field, alternative="two-sided")
    return {
        "gc_mean": float(np.mean(gc)),
        "field_mean": float(np.mean(field)),
        "diff_dex": float(np.mean(gc) - np.mean(field)),
        "t_stat": float(t_stat),
        "t_p": float(p_value),
        "mw_u": float(mw_u),
        "mw_p": float(mw_p),
        "gc_n": int(len(gc)),
        "field_n": int(len(field)),
    }


def _period_matched_bootstrap(gc_rows: list[dict], field_rows: list[dict], n_boot=2000, seed=42) -> dict:
    """Bootstrap period-matched comparison WITHOUT replacement.

    For each GC pulsar, select the nearest field pulsar in logP.
    Each field pulsar is used at most once to avoid bias from overmatching.
    """

    rng = np.random.default_rng(seed)
    gc_logp = np.array([r["logP"] for r in gc_rows])
    gc_logpdot = np.array([r["logPdot_abs"] for r in gc_rows])

    field_logp = np.array([r["logP"] for r in field_rows])
    field_logpdot = np.array([r["logPdot_abs"] for r in field_rows])

    # Pre-compute all pairwise distances
    n_gc = len(gc_rows)
    n_field = len(field_rows)
    distances = np.zeros((n_gc, n_field))
    for i in range(n_gc):
        distances[i, :] = np.abs(field_logp - gc_logp[i])

    diffs = []
    for _ in range(n_boot):
        # Resample GC pulsars with replacement (bootstrap)
        idx_gc = rng.integers(0, n_gc, size=n_gc)
        
        # Match WITHOUT replacement
        used_field = set()
        f_sel = []
        
        # Randomize order
        order = rng.permutation(len(idx_gc))
        
        for idx in order:
            i = idx_gc[idx]
            # Find nearest unused field pulsar
            sorted_indices = np.argsort(distances[i, :])
            for j in sorted_indices:
                if j not in used_field:
                    used_field.add(j)
                    f_sel.append(field_logpdot[j])
                    break
        
        f_sel = np.array(f_sel)
        diffs.append(float(np.mean(gc_logpdot[idx_gc]) - np.mean(f_sel)))

    diffs = np.array(diffs)
    return {
        "n_boot": int(n_boot),
        "diff_mean": float(np.mean(diffs)),
        "diff_ci16": float(np.quantile(diffs, 0.16)),
        "diff_ci84": float(np.quantile(diffs, 0.84)),
        "p_two_sided": float(2 * min(np.mean(diffs <= 0), np.mean(diffs >= 0))),
    }


def _two_dim_match_bootstrap(gc_rows: list[dict], field_rows: list[dict], n_boot=2000, seed=42) -> dict:
    """Bootstrap matching in (logP, log_b_proxy) WITHOUT replacement.

    Each GC pulsar is matched to a unique field pulsar (no reuse) to avoid
    bias from overmatching. Uses greedy nearest neighbor with random order
    to ensure fair matching.
    """

    rng = np.random.default_rng(seed)
    gc_x = np.array([[r["logP"], r["log_b_proxy"]] for r in gc_rows])
    gc_y = np.array([r["logPdot_abs"] for r in gc_rows])

    field_x = np.array([[r["logP"], r["log_b_proxy"]] for r in field_rows])
    field_y = np.array([r["logPdot_abs"] for r in field_rows])

    # Pre-compute all pairwise distances
    n_gc = len(gc_rows)
    n_field = len(field_rows)
    distances = np.zeros((n_gc, n_field))
    for i in range(n_gc):
        dx = field_x[:, 0] - gc_x[i, 0]
        dy = field_x[:, 1] - gc_x[i, 1]
        distances[i, :] = np.sqrt(dx*dx + dy*dy)

    diffs = []
    for _ in range(n_boot):
        # Resample GC pulsars with replacement (bootstrap)
        idx_gc = rng.integers(0, n_gc, size=n_gc)
        
        # Match WITHOUT replacement: each field pulsar used at most once
        used_field = set()
        f_sel = []
        
        # Randomize order to avoid systematic bias
        order = rng.permutation(len(idx_gc))
        
        for idx in order:
            i = idx_gc[idx]
            # Find nearest unused field pulsar
            sorted_indices = np.argsort(distances[i, :])
            for j in sorted_indices:
                if j not in used_field:
                    used_field.add(j)
                    f_sel.append(field_y[j])
                    break
        
        f_sel = np.array(f_sel)
        diffs.append(float(np.mean(gc_y[idx_gc]) - np.mean(f_sel)))

    diffs = np.array(diffs)
    return {
        "n_boot": int(n_boot),
        "diff_mean": float(np.mean(diffs)),
        "diff_ci16": float(np.quantile(diffs, 0.16)),
        "diff_ci84": float(np.quantile(diffs, 0.84)),
        "p_two_sided": float(2 * min(np.mean(diffs <= 0), np.mean(diffs >= 0))),
    }


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # --- Freire ---
    freire_raw, freire_sha = _download_with_retry(FREIRE_GCPSR_URL)
    RAW_FREIRE_PATH.write_bytes(freire_raw)
    freire_rows = _parse_freire_gcpsr(freire_raw.decode("utf-8", errors="replace"))

    # Cluster tokens for conservative GC exclusion on ATNF side
    freire_cluster_tokens = set()
    for r in freire_rows:
        if r.get("cluster"):
            for tok in re.split(r"[^A-Za-z0-9]+", r["cluster"]):
                tok = tok.strip()
                if tok:
                    freire_cluster_tokens.add(tok)

    # --- ATNF psrcat package ---
    atnf_tgz_raw, atnf_sha = _download_with_retry(ATNF_PSRCAT_TGZ_URL)
    RAW_ATNF_TGZ_PATH.write_bytes(atnf_tgz_raw)
    db_raw = _extract_psrcat_db_from_tgz(atnf_tgz_raw)
    RAW_ATNF_DB_PATH.write_bytes(db_raw)
    atnf_rows_raw = _parse_atnf_psrcat_db(db_raw.decode("utf-8", errors="replace"))

    # Define samples
    # MSP cut
    def is_msp(r):
        return r["P0_s"] is not None and r["P0_s"] < 0.030

    gc_rows = [r for r in freire_rows if is_msp(r)]

    # Field MSPs from ATNF: MSP cut + P1 present + exclude anything that looks like a GC association
    field_rows = []
    for r in atnf_rows_raw:
        if not is_msp(r):
            continue
        if _looks_like_gc_assoc(r.get("assoc", ""), freire_cluster_tokens):
            continue
        if r.get("TYPE") and "gc" in str(r.get("TYPE")).lower():
            continue
        if any(fr["name"] == r.get("name") for fr in freire_rows):
            continue
        r2 = dict(r)
        r2["environment"] = "field"
        field_rows.append(r2)

    # Compute proxies
    gc_rows_p = _compute_proxies(gc_rows)
    field_rows_p = _compute_proxies(field_rows)

    if len(field_rows_p) == 0:
        raise RuntimeError(
            "ATNF parsing produced zero field MSPs after exclusions. "
            "This indicates a parsing/filtering error or a connectivity issue with ATNF download."
        )

    # Base comparison
    gc_logpdot = np.array([r["logPdot_abs"] for r in gc_rows_p])
    field_logpdot = np.array([r["logPdot_abs"] for r in field_rows_p])
    base = _ttest_logpdot(gc_logpdot, field_logpdot)

    # Controls
    period_match = _period_matched_bootstrap(gc_rows_p, field_rows_p)
    two_dim_match = _two_dim_match_bootstrap(gc_rows_p, field_rows_p)

    meta = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "freire_gcpsr": {"url": FREIRE_GCPSR_URL, "sha256": freire_sha, "bytes": len(freire_raw)},
            "atnf_psrcat_pkg": {"url": ATNF_PSRCAT_TGZ_URL, "sha256": atnf_sha, "bytes": len(atnf_tgz_raw)},
        },
        "selection": {
            "msp_cut": "P0 < 0.03 s",
            "field_exclusion": "Exclude ATNF rows with ASSOC containing 'GC' or tokens matching Freire cluster headers",
            "note": "GC observed Pdot can be contaminated by line-of-sight acceleration; tau_c and B proxies are treated as proxies, not intrinsic ages/fields.",
        },
        "counts": {
            "gc_msp": len(gc_rows_p),
            "field_msp": len(field_rows_p),
        },
    }

    out = {
        "meta": meta,
        "base_log10_abs_pdot": base,
        "controls": {
            "period_matched": period_match,
            "period_and_bproxy_matched": two_dim_match,
        },
    }

    OUT_JSON.write_text(json.dumps(out, indent=2))

    # Flat CSV for downstream inspection
    with OUT_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "source",
                "environment",
                "cluster",
                "name",
                "P0_s",
                "P_ms",
                "P1_sps",
                "assoc",
                "logP",
                "logPdot_abs",
                "log_b_proxy",
                "log_tau_c",
            ]
        )
        for r in gc_rows_p + field_rows_p:
            w.writerow(
                [
                    r.get("source"),
                    r.get("environment"),
                    r.get("cluster"),
                    r.get("name"),
                    f"{r.get('P0_s'):.12g}",
                    f"{r.get('P_ms'):.6g}",
                    f"{r.get('P1_sps'):.12g}",
                    r.get("assoc", ""),
                    f"{r.get('logP'):.8g}",
                    f"{r.get('logPdot_abs'):.8g}",
                    f"{r.get('log_b_proxy'):.8g}",
                    f"{r.get('log_tau_c'):.8g}",
                ]
            )

    md = []
    md.append("# Pulsar Population Controls (Freire + ATNF)\n")
    md.append(f"**Freire GCpsr URL:** {FREIRE_GCPSR_URL}\\\n")
    md.append(f"**Freire SHA256:** `{freire_sha}`\\\n")
    md.append(f"**ATNF psrcat_pkg URL:** {ATNF_PSRCAT_TGZ_URL}\\\n")
    md.append(f"**ATNF SHA256:** `{atnf_sha}`\\\n")
    md.append("\n")
    md.append("## Sample sizes\n")
    md.append(f"- **GC MSPs (Freire, P<30 ms, measured Pdot):** {len(gc_rows_p)}\n")
    md.append(f"- **Field MSPs (ATNF, P<30 ms, Pdot present, non-GC ASSOC):** {len(field_rows_p)}\n")
    md.append("\n")
    md.append("## Base test (log10|Pdot|)\n")
    md.append(f"- **GC mean:** {base['gc_mean']:.3f}\\\n")
    md.append(f"- **Field mean:** {base['field_mean']:.3f}\\\n")
    md.append(f"- **Difference (GC-Field):** {base['diff_dex']:.3f} dex\\\n")
    md.append(f"- **Welch t-test p:** {base['t_p']:.3g}\\\n")
    md.append(f"- **Mann-Whitney p:** {base['mw_p']:.3g}\\\n")
    md.append("\n")
    md.append("## Controls\n")
    md.append("### Period-matched bootstrap\n")
    md.append(
        f"- **Mean diff:** {period_match['diff_mean']:.3f} dex (16–84%: {period_match['diff_ci16']:.3f} to {period_match['diff_ci84']:.3f})\\\n"
    )
    md.append(f"- **Two-sided p:** {period_match['p_two_sided']:.3g}\\\n")
    md.append("\n")
    md.append("### Period + B-proxy matched bootstrap\n")
    md.append(
        f"- **Mean diff:** {two_dim_match['diff_mean']:.3f} dex (16–84%: {two_dim_match['diff_ci16']:.3f} to {two_dim_match['diff_ci84']:.3f})\\\n"
    )
    md.append(f"- **Two-sided p:** {two_dim_match['p_two_sided']:.3g}\\\n")

    OUT_MD.write_text("".join(md))

    print(f"Wrote: {OUT_JSON}")
    print(f"Wrote: {OUT_CSV}")
    print(f"Wrote: {OUT_MD}")


if __name__ == "__main__":
    # Setup file logging when run manually
    import sys
    from pathlib import Path
    # Add repo root and scripts directory to path
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root))
    sys.path.insert(0, str(repo_root / "scripts"))
    from utils.logger import setup_step_logger
    logger = setup_step_logger("step_5_10_pulsar_population_controls")
    main()
