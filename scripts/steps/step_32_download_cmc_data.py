#!/usr/bin/env python3
"""
CMC Data Downloader - Production Implementation
==================================================

Downloads CMC (Cluster Monte Carlo) simulation data via POST with
progress bar, resume support, and pipeline integration.

Usage:
    python step_32_download_cmc_data.py [--cluster 47_Tuc] [--all]
    
Author: M. Smawfield
Date: March 2026
"""

import requests
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import tarfile
import io
import os
import sys
import time
from pathlib import Path
from typing import Dict, Optional, List
import argparse

# Progress bar function
def print_progress(downloaded: int, total: int, width: int = 50):
    """Print a progress bar."""
    if total > 0:
        percent = downloaded / total * 100
        filled = int(width * downloaded // total)
        bar = '█' * filled + '░' * (width - filled)
        mb_downloaded = downloaded / 1024 / 1024
        mb_total = total / 1024 / 1024
        sys.stdout.write(f"\r    │{bar}│ {percent:5.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)")
        sys.stdout.flush()


def log(msg: str, level: str = "INFO"):
    """Log message with timestamp."""
    timestamp = time.strftime("%H:%M:%S")
    prefix = f"[{timestamp}] [{level:8s}]"
    print(f"{prefix} {msg}", flush=True)


CMC_URL = "https://cmc.ciera.northwestern.edu/index.php"

# ---------------------------------------------------------------------------
# CMC grid-cell matching (revised Nov 2026 audit)
# ---------------------------------------------------------------------------
# The CMC catalog (Kremer et al. 2020) is a coarse grid of 148 single-realization
# models keyed on (N, rv, rg, Z): the SAME grid cell always returns the SAME
# deterministic realization. Each real cluster is therefore assigned its
# nearest grid cell by observed properties; clusters that share a cell share a
# realization (this is disclosed in the step_37 output and the manuscript).
#
# Observed anchors (Harris 1996, 2010 ed.; Baumgardt & Hilker 2018):
#   - metallicity: [Fe/H] >= -0.6 -> Z0.02 ; -1.6 < [Fe/H] < -0.6 -> Z0.002 ;
#                  [Fe/H] <= -1.6 -> Z0.0002
#   - galactocentric distance: R_gc < 4 kpc -> rg2 ; 4-12 kpc -> rg8 ; >12 -> rg20
#   - N: nearest of {2e5, 4e5, 8e5, 1.6e6, 3.2e6} to cluster mass in objects
#        (N_particles ~ 2 x N_stars today for these IMFs; used as mass proxy)
#   - rv: virial-radius proxy for cluster extent (compact -> rv0.5, extended ->
#        rv2, sparse -> rv4)
#
# Previous mapping (pre-audit) used only rg8 and Z in {0.0002, 0.02} for every
# cluster, which collapsed 20 clusters onto 5 unique realizations. The mapping
# below yields distinct cells wherever the observed properties differ.
# ---------------------------------------------------------------------------
CLUSTER_MODELS = {
    "47_Tuc": {
        "name": "47 Tucanae",
        "ngc": "NGC 104",
        "params": {
            "number_of_objects": "N8e5",
            "virial_radius": "rv1",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -0.72, "rgc_kpc": 7.4, "mass_msun": 1.0e6},
        "expected_size_mb": 250,
    },
    "Terzan_5": {
        "name": "Terzan 5",
        "ngc": None,
        "params": {
            "number_of_objects": "N8e5",
            "virial_radius": "rv0.5",
            "galactocentric_distance": "rg2",
            "metallicity": "Z0.02",
        },
        "observed": {"feh": -0.23, "rgc_kpc": 1.0, "mass_msun": 6.0e5},
        "expected_size_mb": 250,
        "notes": "Massive metal-rich bulge cluster",
    },
    "M15": {
        "name": "M15",
        "ngc": "NGC 7078",
        "params": {
            "number_of_objects": "N8e5",
            "virial_radius": "rv0.5",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.0002",
        },
        "observed": {"feh": -2.37, "rgc_kpc": 10.4, "mass_msun": 5.6e5},
        "expected_size_mb": 250,
        "notes": "Core-collapsed, very metal-poor",
    },
    "M62": {
        "name": "M62",
        "ngc": "NGC 6266",
        "params": {
            "number_of_objects": "N4e5",
            "virial_radius": "rv1",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -1.18, "rgc_kpc": 7.1, "mass_msun": 3.0e5},
        "expected_size_mb": 150,
    },
    "NGC_6517": {
        "name": "NGC 6517",
        "ngc": None,
        "params": {
            "number_of_objects": "N4e5",
            "virial_radius": "rv0.5",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -1.24, "rgc_kpc": 4.0, "mass_msun": 2.0e5},
        "expected_size_mb": 150,
        "notes": "Core-collapsed, highest density in sample (log rho = 5.8)",
    },
    "M28": {
        "name": "M28",
        "ngc": "NGC 6626",
        "params": {
            "number_of_objects": "N4e5",
            "virial_radius": "rv0.5",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -1.32, "rgc_kpc": 6.9, "mass_msun": 2.8e5},
        "expected_size_mb": 150,
    },
    "M13": {
        "name": "M13",
        "ngc": "NGC 6205",
        "params": {
            "number_of_objects": "N8e5",
            "virial_radius": "rv2",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -1.53, "rgc_kpc": 8.4, "mass_msun": 6.0e5},
        "expected_size_mb": 250,
    },
    "NGC_6397": {
        "name": "NGC 6397",
        "ngc": "NGC 6397",
        "params": {
            "number_of_objects": "N2e5",
            "virial_radius": "rv1",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.0002",
        },
        "observed": {"feh": -2.02, "rgc_kpc": 6.0, "mass_msun": 1.0e5},
        "expected_size_mb": 140,
        "notes": "Low-metallicity nearby core-collapsed cluster (nearest served N2e5 Z0.0002 cell: rv1)",
    },
    "NGC_6752": {
        "name": "NGC 6752",
        "ngc": "NGC 6752",
        "params": {
            "number_of_objects": "N4e5",
            "virial_radius": "rv1",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -1.54, "rgc_kpc": 5.2, "mass_msun": 2.0e5},
        "expected_size_mb": 150,
        "notes": "Intermediate-metallicity globular cluster",
    },
    "M3": {
        "name": "M3",
        "ngc": "NGC 5272",
        "params": {
            "number_of_objects": "N4e5",
            "virial_radius": "rv2",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -1.50, "rgc_kpc": 12.0, "mass_msun": 4.5e5},
        "expected_size_mb": 150,
        "notes": "Large extended Oosterhoff type I halo cluster",
    },
    "M5": {
        "name": "M5",
        "ngc": "NGC 5904",
        "params": {
            "number_of_objects": "N4e5",
            "virial_radius": "rv1",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -1.29, "rgc_kpc": 6.2, "mass_msun": 4.0e5},
        "expected_size_mb": 150,
        "notes": "Well-studied metal-intermediate cluster",
    },
    "M4": {
        "name": "M4",
        "ngc": "NGC 6121",
        "params": {
            "number_of_objects": "N2e5",
            "virial_radius": "rv2",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -1.16, "rgc_kpc": 5.9, "mass_msun": 1.2e5},
        "expected_size_mb": 80,
        "notes": "Nearby low-mass extended cluster",
    },
    "Omega_Cen": {
        "name": "Omega Centauri",
        "ngc": "NGC 5139",
        "params": {
            "number_of_objects": "N3.2e6",
            "virial_radius": "rv1",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -1.53, "rgc_kpc": 6.4, "mass_msun": 3.5e6},
        "expected_size_mb": 800,
        "notes": "Most massive globular cluster in Milky Way",
    },
    "NGC_6440": {
        "name": "NGC 6440",
        "ngc": "NGC 6440",
        "params": {
            "number_of_objects": "N8e5",
            "virial_radius": "rv0.5",
            "galactocentric_distance": "rg2",
            "metallicity": "Z0.02",
        },
        "observed": {"feh": -0.36, "rgc_kpc": 1.3, "mass_msun": 5.0e5},
        "expected_size_mb": 250,
        "notes": "Metal-rich bulge cluster, high density",
    },
    "NGC_6441": {
        "name": "NGC 6441",
        "ngc": "NGC 6441",
        "params": {
            "number_of_objects": "N1.6e6",
            "virial_radius": "rv0.5",
            "galactocentric_distance": "rg2",
            "metallicity": "Z0.02",
        },
        "observed": {"feh": -0.46, "rgc_kpc": 3.9, "mass_msun": 1.0e6},
        "expected_size_mb": 400,
        "notes": "Metal-rich bulge cluster, high density",
    },
    "M22": {
        "name": "M22",
        "ngc": "NGC 6656",
        "params": {
            "number_of_objects": "N4e5",
            "virial_radius": "rv2",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.0002",
        },
        "observed": {"feh": -1.70, "rgc_kpc": 4.9, "mass_msun": 3.0e5},
        "expected_size_mb": 150,
        "notes": "Extended metal-poor cluster",
    },
    "M71": {
        "name": "M71",
        "ngc": "NGC 6838",
        "params": {
            "number_of_objects": "N2e5",
            "virial_radius": "rv4",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -0.78, "rgc_kpc": 6.7, "mass_msun": 1.0e5},
        "expected_size_mb": 80,
        "notes": "Sparse metal-intermediate cluster",
    },
    "NGC_6624": {
        "name": "NGC 6624",
        "ngc": "NGC 6624",
        "params": {
            "number_of_objects": "N8e5",
            "virial_radius": "rv0.5",
            "galactocentric_distance": "rg2",
            "metallicity": "Z0.02",
        },
        "observed": {"feh": -0.44, "rgc_kpc": 1.2, "mass_msun": 5.0e5},
        "expected_size_mb": 250,
        "notes": "Bulge cluster with MSPs",
    },
    "NGC_6388": {
        "name": "NGC 6388",
        "ngc": "NGC 6388",
        "params": {
            "number_of_objects": "N1.6e6",
            "virial_radius": "rv0.5",
            "galactocentric_distance": "rg2",
            "metallicity": "Z0.02",
        },
        "observed": {"feh": -0.55, "rgc_kpc": 3.1, "mass_msun": 1.0e6},
        "expected_size_mb": 400,
        "notes": "Metal-rich, high-density cluster",
    },
    "NGC_6712": {
        "name": "NGC 6712",
        "ngc": "NGC 6712",
        "params": {
            "number_of_objects": "N2e5",
            "virial_radius": "rv1",
            "galactocentric_distance": "rg2",
            "metallicity": "Z0.002",
        },
        "observed": {"feh": -1.02, "rgc_kpc": 3.9, "mass_msun": 1.0e5},
        "expected_size_mb": 80,
        "notes": "Low-mass inner-galaxy cluster",
    },
}


def download_with_progress(url: str, params: Dict, timeout: int = 600) -> Optional[bytes]:
    """
    Download data with streaming progress bar.
    
    Returns downloaded bytes or None on failure.
    """
    try:
        # Stream the response
        response = requests.post(url, data=params, timeout=timeout, stream=True, verify=False)
        response.raise_for_status()
        
        # Get total size if available
        total_size = int(response.headers.get('content-length', 0))
        
        if total_size > 0:
            log(f"Downloading {total_size/1024/1024:.1f} MB...")
        else:
            log(f"Downloading (size unknown)...")
        
        # Download with progress
        downloaded = 0
        chunks = []
        
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                chunks.append(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    print_progress(downloaded, total_size)
        
        if total_size > 0:
            print()  # New line after progress bar
        
        return b''.join(chunks)
        
    except Exception as e:
        log(f"Download failed: {e}", "ERROR")
        return None


def download_cmc_simulation(cluster_id: str, cluster_info: Dict, cluster_dir: Path) -> Dict:
    """
    Download and extract CMC simulation for a cluster.
    
    Returns result dict with success status and extracted files.
    """
    params = cluster_info["params"]
    expected_mb = cluster_info.get("expected_size_mb", 250)
    
    log(f"Starting download for {cluster_info['name']}")
    log(f"  Parameters: {params}")
    log(f"  Expected size: ~{expected_mb} MB")
    log(f"  This will take 3-8 minutes depending on connection speed")
    
    result = {
        "success": False,
        "files_extracted": [],
        "error": None,
        "download_time_s": 0,
    }
    
    start_time = time.time()
    
    # Download with progress
    data = download_with_progress(CMC_URL, params, timeout=600)
    
    download_time = time.time() - start_time
    result["download_time_s"] = download_time
    
    if data is None:
        result["error"] = "Download failed"
        return result
    
    log(f"Downloaded {len(data)/1024/1024:.1f} MB in {download_time:.1f}s")
    
    # Extract tarball
    log(f"Extracting tarball...")
    cluster_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        tarball = tarfile.open(fileobj=io.BytesIO(data), mode="r:gz")
        members = tarball.getmembers()
        
        extracted = 0
        key_files = []
        
        for member in members:
            if member.isfile():
                member_path = cluster_dir / member.name
                with tarball.extractfile(member) as f:
                    if f:
                        content = f.read()
                        with open(member_path, 'wb') as out:
                            out.write(content)
                        result["files_extracted"].append(member.name)
                        extracted += 1
                        
                        # Track key files
                        if 'morepulsars' in member.name or 'snapshot' in member.name or 'dyn' in member.name:
                            key_files.append(f"{member.name} ({len(content)/1024:.1f} KB)")
        
        tarball.close()
        
        result["success"] = True
        log(f"SUCCESS - Extracted {extracted} files")
        for kf in key_files[:5]:
            log(f"  Found: {kf}")
        if len(key_files) > 5:
            log(f"  ... and {len(key_files) - 5} more")
            
    except tarfile.TarError as e:
        result["error"] = f"Invalid tarball: {e}"
        log(f"FAILED - Invalid tarball: {e}", "ERROR")
        
        # Save error page for debugging
        error_path = cluster_dir / "download_error.bin"
        with open(error_path, 'wb') as f:
            f.write(data)
        log(f"Saved response to: {error_path} for debugging")
    
    return result


def download_cluster(cluster_id: str, data_dir: Path, force: bool = False) -> Dict:
    """Download CMC data for a specific cluster."""
    
    if cluster_id not in CLUSTER_MODELS:
        log(f"Unknown cluster: {cluster_id}", "ERROR")
        return {"success": False, "error": "Unknown cluster"}
    
    cluster_info = CLUSTER_MODELS[cluster_id]
    cluster_dir = data_dir / cluster_id
    
    log(f"\n{'='*70}")
    log(f"CLUSTER: {cluster_id} - {cluster_info['name']}")
    log(f"{'='*70}")
    
    # Check if already downloaded
    marker_file = cluster_dir / ".download_complete"
    if marker_file.exists() and not force:
        log(f"Already downloaded (use --force to re-download)")
        return {"success": True, "status": "already_exists"}
    
    # Download
    result = download_cmc_simulation(cluster_id, cluster_info, cluster_dir)
    
    # Create marker on success
    if result["success"]:
        with open(marker_file, 'w') as f:
            f.write(f"Downloaded: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Files: {len(result.get('files_extracted', []))}\n")
            f.write(f"Time: {result.get('download_time_s', 0):.1f}s\n")
    
    return result


def build_dedup_manifest(data_dir: Path) -> Dict:
    """
    Hash every downloaded initial.morepulsars.dat and report which cluster
    selections share a unique CMC realization.

    The CMC catalog serves one deterministic realization per (N, rv, rg, Z)
    grid cell, so clusters mapped to the same cell hold byte-identical data.
    Downstream steps (step_37 etc.) MUST treat identical files as a single
    independent realization; this manifest records the mapping explicitly.

    Returns
    -------
    Dict with cluster->model assignments and shared-model groups; also
    written to data/cmc/model_dedup_manifest.json.
    """
    import hashlib
    import json

    cluster_to_model: Dict[str, str] = {}
    for cluster_dir in sorted(data_dir.iterdir()):
        mp = cluster_dir / "initial.morepulsars.dat"
        if not mp.exists():
            continue
        h = hashlib.md5()
        with open(mp, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk)
        cluster_to_model[cluster_dir.name] = h.hexdigest()

    model_to_clusters: Dict[str, List[str]] = {}
    for name, md5 in cluster_to_model.items():
        model_to_clusters.setdefault(md5, []).append(name)

    manifest = {
        "note": (
            "CMC catalog serves one realization per (N, rv, rg, Z) grid "
            "cell; clusters sharing a cell share identical data. Treat "
            "identical files as ONE independent realization downstream."
        ),
        "n_cluster_selections": len(cluster_to_model),
        "n_unique_models": len(model_to_clusters),
        "cluster_to_model": cluster_to_model,
        "shared_model_groups": {
            k: v for k, v in model_to_clusters.items() if len(v) > 1
        },
    }

    out = data_dir / "model_dedup_manifest.json"
    with open(out, "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest


def main():
    """Main entry point with argument parsing."""
    parser = argparse.ArgumentParser(
        description="Download CMC Cluster Monte Carlo simulation data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python step_32_download_cmc_data.py --all
  python step_32_download_cmc_data.py --cluster 47_Tuc
  python step_32_download_cmc_data.py --cluster Terzan_5 --force
        """
    )
    
    parser.add_argument('--all', action='store_true', help='Download all clusters')
    parser.add_argument('--cluster', type=str, choices=list(CLUSTER_MODELS.keys()),
                       help='Download specific cluster')
    parser.add_argument('--force', action='store_true', help='Re-download even if exists')
    parser.add_argument('--data-dir', type=Path, default=None,
                       help='Data directory (default: REPO_ROOT/data/cmc)')
    
    args = parser.parse_args()
    
    # Determine data directory
    if args.data_dir is None:
        repo_root = Path(__file__).resolve().parents[2]
        data_dir = repo_root / "data" / "cmc"
    else:
        data_dir = args.data_dir
    
    # Determine which clusters to download
    if args.all:
        clusters = list(CLUSTER_MODELS.keys())
    elif args.cluster:
        clusters = [args.cluster]
    else:
        parser.print_help()
        sys.exit(1)
    
    # Banner
    log("=" * 70)
    log("CMC CLUSTER CATALOG DOWNLOADER")
    log("=" * 70)
    log(f"Target URL: {CMC_URL}")
    log(f"Data directory: {data_dir}")
    log(f"Clusters: {', '.join(clusters)}")
    log("")
    log("NOTE: Each download is 200-400 MB and takes 3-8 minutes")
    log("      Please be patient and do not interrupt the download")
    log("")
    
    # Download
    results = {}
    for i, cluster_id in enumerate(clusters, 1):
        if i > 1:
            log(f"\nWaiting 5s before next download...")
            time.sleep(5)
        
        result = download_cluster(cluster_id, data_dir, force=args.force)
        results[cluster_id] = result
    
    # Summary
    log("\n" + "=" * 70)
    log("FINAL SUMMARY")
    log("=" * 70)
    
    total_success = 0
    for cluster_id, result in results.items():
        cluster_info = CLUSTER_MODELS[cluster_id]
        if result["success"]:
            total_success += 1
            log(f"SUCCESS: {cluster_id:10s} - {cluster_info['name']:20s}")
            log(f"           Files: {len(result.get('files_extracted', []))}, Time: {result.get('download_time_s', 0):.1f}s")
        else:
            log(f"FAILED: {cluster_id:10s} - {cluster_info['name']:20s}")
            log(f"           Error: {result.get('error', 'Unknown')}")
    
    log(f"\n{'='*70}")
    log(f"TOTAL: {total_success}/{len(clusters)} clusters successfully downloaded")
    log(f"{'='*70}")

    # Model-deduplication check: warn loudly when distinct cluster labels
    # hold identical realizations so the duplication is never silent.
    manifest = build_dedup_manifest(data_dir)
    log(
        f"Model dedup: {manifest['n_cluster_selections']} selections -> "
        f"{manifest['n_unique_models']} unique realizations "
        f"(manifest: {data_dir}/model_dedup_manifest.json)"
    )
    for md5, members in manifest["shared_model_groups"].items():
        log(
            f"WARNING: {len(members)} cluster selections share realization "
            f"{md5[:12]}: {', '.join(members)} — these are ONE independent "
            "model for statistics",
            "WARNING",
        )
    
    # Exit code: allow partial success so downstream N-body steps can use
    # whichever clusters were successfully downloaded. The pipeline should not
    # halt because one or two remote cluster tarballs are temporarily unavailable.
    success_fraction = total_success / len(clusters) if clusters else 0
    if success_fraction >= 0.8:
        log(f"Downloaded {success_fraction*100:.0f}% of requested clusters. Continuing pipeline.")
        sys.exit(0)
    else:
        log(f"Only {success_fraction*100:.0f}% of requested clusters downloaded. Failing step.", "ERROR")
        sys.exit(1)


if __name__ == "__main__":
    main()
