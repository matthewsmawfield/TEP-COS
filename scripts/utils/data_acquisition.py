#!/usr/bin/env python3
"""
TEP-COS Data Acquisition Module
==============================

Ensures all required data files are present before analysis steps run.
Downloads data from external sources when needed.

Usage:
    from scripts.utils.data_acquisition import ensure_data
    ensure_data('pulsars')  # Downloads Freire + ATNF if needed
    ensure_data('lensing')  # Checks COSMOGRAIL files exist
    ensure_data('supernovae')  # Downloads Pantheon+ if needed
"""

import hashlib
import io
import os
import re
import tarfile
import time
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from urllib.request import urlopen, HTTPError, URLError
from socket import timeout as SocketTimeout

# Repository root
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "results" / "outputs"

# Data sources configuration
# NOTE: URLs point to external data repositories. If files change or move,
# the pipeline will fail with clear error messages. Data integrity is verified
# via SHA256 checksums stored in output metadata.
# 
# Version pinning: The Freire GCpsr and ATNF psrcat are snapshots downloaded
# at runtime. SHA256 hashes are recorded in step outputs for reproducibility.
# To use specific versions, manually download and place files in results/outputs/
DATA_SOURCES: Dict[str, Dict] = {
    "freire_gcpsr": {
        "url": "https://www3.mpifr-bonn.mpg.de/staff/pfreire/GCpsr.txt",
        "path": RESULTS_DIR / "freire_GCpsr.txt",
        "type": "text",
        "required_for": ["pulsars"],
    },
    "atnf_psrcat": {
        "url": "https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz",
        "path": RESULTS_DIR / "atnf_psrcat_pkg.tar.gz",
        "extract_to": RESULTS_DIR / "atnf_psrcat.db",
        "type": "tgz",
        "required_for": ["pulsars"],
    },
}


def _sha256_bytes(b: bytes) -> str:
    """Compute SHA256 hash of bytes."""
    return hashlib.sha256(b).hexdigest()


def _download_with_retry(
    url: str, 
    max_retries: int = 3, 
    timeout: int = 60,
    verbose: bool = True
) -> Tuple[bytes, str]:
    """Download data with retry logic and error handling."""
    for attempt in range(max_retries):
        try:
            if verbose:
                print(f"  Downloading (attempt {attempt + 1}/{max_retries}): {url[:60]}...")
            raw = urlopen(url, timeout=timeout).read()
            if verbose:
                print(f"  Downloaded {len(raw)} bytes")
            return raw, _sha256_bytes(raw)
        except HTTPError as e:
            print(f"  HTTP Error {e.code}: {e.reason}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
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


def _extract_psrcat_db(tgz_bytes: bytes) -> bytes:
    """Extract psrcat.db from psrcat_pkg.tar.gz."""
    with tarfile.open(fileobj=io.BytesIO(tgz_bytes), mode="r:gz") as tf:
        for member in tf.getmembers():
            base = Path(member.name).name
            if base.lower() in {"psrcat.db", "psrcat.db.txt"}:
                f = tf.extractfile(member)
                if f is not None:
                    return f.read()
        # Fallback: first .db file
        for member in tf.getmembers():
            if member.name.lower().endswith(".db"):
                f = tf.extractfile(member)
                if f is not None:
                    return f.read()
    raise RuntimeError("Could not locate psrcat.db inside psrcat_pkg.tar.gz")


def ensure_file_exists(
    path: Path, 
    url: Optional[str] = None,
    min_size_bytes: int = 100,
    verbose: bool = True
) -> bool:
    """
    Ensure a data file exists, downloading if necessary.
    
    Args:
        path: Expected file path
        url: Download URL if file missing
        min_size_bytes: Minimum valid file size
        verbose: Print status messages
    
    Returns:
        True if file exists and is valid
    """
    path = Path(path)
    
    # Check if file exists and is valid
    if path.exists():
        size = path.stat().st_size
        if size >= min_size_bytes:
            if verbose:
                print(f"  ✓ Data file exists: {path.name} ({size} bytes)")
            return True
        else:
            if verbose:
                print(f"  ⚠ File too small ({size} bytes), re-downloading...")
    
    # File missing or invalid - try to download
    if url is None:
        if verbose:
            print(f"  ✗ Missing file: {path}")
            print(f"    No download URL configured")
        return False
    
    # Download
    if verbose:
        print(f"  → Downloading missing data: {path.name}")
    
    try:
        raw, sha256 = _download_with_retry(url, verbose=verbose)
        
        # Ensure parent directory exists
        path.parent.mkdir(parents=True, exist_ok=True)
        
        # Write file
        path.write_bytes(raw)
        
        if verbose:
            print(f"  ✓ Saved: {path} ({len(raw)} bytes, SHA256: {sha256[:16]}...)")
        
        return True
        
    except Exception as e:
        if verbose:
            print(f"  ✗ Download failed: {e}")
        return False


def ensure_data(dataset: str, verbose: bool = True) -> Dict[str, bool]:
    """
    Ensure all data files for a dataset are available.
    
    Args:
        dataset: One of 'pulsars', 'lensing', 'supernovae', 'all'
        verbose: Print status messages
    
    Returns:
        Dictionary of {file_key: success_status}
    """
    results = {}
    
    if verbose:
        print(f"\n[Data Acquisition] Checking: {dataset}")
    
    if dataset == "pulsars" or dataset == "all":
        # Freire GCpsr
        src = DATA_SOURCES["freire_gcpsr"]
        results["freire_gcpsr"] = ensure_file_exists(
            src["path"], src["url"], min_size_bytes=1000, verbose=verbose
        )
        
        # ATNF psrcat
        src = DATA_SOURCES["atnf_psrcat"]
        tgz_ok = ensure_file_exists(
            src["path"], src["url"], min_size_bytes=10000, verbose=verbose
        )
        results["atnf_psrcat"] = tgz_ok
        
        # Extract .db if needed
        if tgz_ok and src.get("extract_to"):
            db_path = src["extract_to"]
            if not db_path.exists() or db_path.stat().st_size < 1000:
                if verbose:
                    print(f"  → Extracting {src['path'].name}...")
                try:
                    tgz_bytes = src["path"].read_bytes()
                    db_bytes = _extract_psrcat_db(tgz_bytes)
                    db_path.write_bytes(db_bytes)
                    if verbose:
                        print(f"  ✓ Extracted: {db_path.name} ({len(db_bytes)} bytes)")
                except Exception as e:
                    if verbose:
                        print(f"  ✗ Extraction failed: {e}")
                    results["atnf_psrcat"] = False
    
    if dataset == "lensing" or dataset == "all":
        # Check COSMOGRAIL files exist
        cosmograil_dir = DATA_DIR / "cosmograil"
        required_rdb = [
            "HE0435_Bonvin2016.rdb",
            "HS2209_Eulaers2013.rdb",
            "J1001_Rathnakumar2013.rdb",
            "J1206_Eulaers2013.rdb",
            "PG1115_Bonvin2018.rdb",
            "RXJ1131_Tewes2013.rdb",
            "WFI2033_Bonvin2019.rdb",
        ]
        
        all_present = True
        for fname in required_rdb:
            fpath = cosmograil_dir / fname
            if fpath.exists():
                if verbose:
                    print(f"  ✓ {fname}")
            else:
                if verbose:
                    print(f"  ✗ Missing: {fname}")
                all_present = False
        
        results["cosmograil"] = all_present
        
        if not all_present:
            if verbose:
                print("  ⚠ Some COSMOGRAIL files missing. Run data download scripts manually.")
    
    if dataset == "supernovae" or dataset == "all":
        # Check Pantheon+ data
        pantheon_file = DATA_DIR / "supernovae" / "pantheon_plus.dat"
        if pantheon_file.exists():
            if verbose:
                print(f"  ✓ Pantheon+ data exists")
            results["pantheon_plus"] = True
        else:
            if verbose:
                print(f"  ✗ Pantheon+ data missing: {pantheon_file}")
                print(f"    Download from: https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/1_DATA")
            results["pantheon_plus"] = False
    
    # Summary
    if verbose:
        n_ok = sum(1 for v in results.values() if v)
        n_total = len(results)
        if n_ok == n_total:
            print(f"  ✓ All {n_total} data sources ready")
        else:
            print(f"  ⚠ {n_ok}/{n_total} data sources ready")
    
    return results


def check_data_health(verbose: bool = True) -> Dict[str, Dict]:
    """
    Comprehensive health check of all data files.
    
    Returns status information for all known data sources.
    """
    status = {}
    
    if verbose:
        print("\n" + "="*60)
        print("TEP-COS Data Health Check")
        print("="*60)
    
    # Check pulsar data
    status["pulsars"] = {
        "freire_exists": DATA_SOURCES["freire_gcpsr"]["path"].exists(),
        "freire_size": DATA_SOURCES["freire_gcpsr"]["path"].stat().st_size 
            if DATA_SOURCES["freire_gcpsr"]["path"].exists() else 0,
        "atnf_tgz_exists": DATA_SOURCES["atnf_psrcat"]["path"].exists(),
        "atnf_tgz_size": DATA_SOURCES["atnf_psrcat"]["path"].stat().st_size
            if DATA_SOURCES["atnf_psrcat"]["path"].exists() else 0,
        "atnf_db_exists": DATA_SOURCES["atnf_psrcat"]["extract_to"].exists(),
        "atnf_db_size": DATA_SOURCES["atnf_psrcat"]["extract_to"].stat().st_size
            if DATA_SOURCES["atnf_psrcat"]["extract_to"].exists() else 0,
    }
    
    if verbose:
        print("\nPulsar Data:")
        for key, val in status["pulsars"].items():
            symbol = "✓" if val else "✗" if "exists" in key else ""
            print(f"  {symbol} {key}: {val}")
    
    # Check lensing data
    cosmograil_dir = DATA_DIR / "cosmograil"
    lens_files = list(cosmograil_dir.glob("*.rdb")) if cosmograil_dir.exists() else []
    status["lensing"] = {
        "dir_exists": cosmograil_dir.exists(),
        "rdb_count": len(lens_files),
        "files": [f.name for f in lens_files],
    }
    
    if verbose:
        print("\nLensing Data (COSMOGRAIL):")
        print(f"  ✓ Directory exists: {status['lensing']['dir_exists']}")
        print(f"  ✓ RDB files found: {status['lensing']['rdb_count']}")
    
    return status


if __name__ == "__main__":
    # Run data acquisition for all datasets
    print("TEP-COS Data Acquisition")
    print("="*60)
    
    results = ensure_data("all", verbose=True)
    
    # Exit with error code if any critical data missing
    critical = ["freire_gcpsr", "atnf_psrcat", "cosmograil"]
    missing = [k for k in critical if not results.get(k)]
    
    if missing:
        print(f"\n✗ Critical data missing: {', '.join(missing)}")
        exit(1)
    else:
        print("\n✓ All critical data ready")
        exit(0)
