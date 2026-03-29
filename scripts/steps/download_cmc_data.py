#!/usr/bin/env python3
"""
CMC Data Downloader - Direct Download Implementation
======================================================

Downloads CMC (Cluster Monte Carlo) simulation data directly via POST 
to https://cmc.ciera.northwestern.edu/index.php

Author: M. Smawfield
Date: March 2026
"""

import requests
import tarfile
import io
from pathlib import Path
from typing import Dict, Optional, List
import time

CMC_URL = "https://cmc.ciera.northwestern.edu/index.php"

# Cluster to CMC parameter mapping (best-fit models)
CLUSTER_MODELS = {
    "47_Tuc": {
        "name": "47 Tucanae",
        "params": {
            "number_of_objects": "N8e5",
            "virial_radius": "rv1",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.02",
        },
    },
    "Terzan_5": {
        "name": "Terzan 5",
        "params": {
            "number_of_objects": "N1.6e6",
            "virial_radius": "rv0.5",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.02",
        },
    },
    "M15": {
        "name": "M15",
        "params": {
            "number_of_objects": "N8e5",
            "virial_radius": "rv0.5",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.0002",
        },
    },
    "M62": {
        "name": "M62",
        "params": {
            "number_of_objects": "N8e5",
            "virial_radius": "rv1",
            "galactocentric_distance": "rg8",
            "metallicity": "Z0.0002",
        },
    },
}


def download_cmc_simulation(params: Dict[str, str], cluster_dir: Path, timeout: int = 120) -> Dict:
    """Download CMC simulation via POST request."""
    result = {"success": False, "files_extracted": [], "error": None}
    
    try:
        print(f"    POST to {CMC_URL}")
        print(f"    Params: {params}")
        
        response = requests.post(CMC_URL, data=params, timeout=timeout)
        response.raise_for_status()
        
        print(f"    Content-Length: {len(response.content)} bytes")
        
        cluster_dir.mkdir(parents=True, exist_ok=True)
        
        # Try to extract as tarball
        try:
            tarball = tarfile.open(fileobj=io.BytesIO(response.content), mode="r:gz")
            for member in tarball.getmembers():
                if member.isfile():
                    member_path = cluster_dir / member.name
                    with tarball.extractfile(member) as f:
                        if f:
                            with open(member_path, 'wb') as out:
                                out.write(f.read())
                            result["files_extracted"].append(member.name)
                            print(f"      Extracted: {member.name}")
            tarball.close()
            result["success"] = True
        except tarfile.TarError as e:
            result["error"] = f"Invalid tarball: {e}"
            error_path = cluster_dir / "download_error.html"
            with open(error_path, 'wb') as f:
                f.write(response.content)
            print(f"    Saved error page to: {error_path}")
            
    except Exception as e:
        result["error"] = str(e)
    
    return result


def download_cluster(cluster_id: str, data_dir: Path) -> Dict:
    """Download CMC data for a specific cluster."""
    if cluster_id not in CLUSTER_MODELS:
        return {"success": False, "error": "Unknown cluster"}
    
    cluster_info = CLUSTER_MODELS[cluster_id]
    cluster_dir = data_dir / cluster_id
    
    print(f"\n  Downloading CMC for {cluster_info['name']}")
    
    marker_file = cluster_dir / ".download_complete"
    if marker_file.exists():
        print(f"    Already downloaded")
        return {"success": True}
    
    result = download_cmc_simulation(cluster_info["params"], cluster_dir)
    
    if result["success"]:
        with open(marker_file, 'w') as f:
            f.write(f"Downloaded: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    return result


def download_all_clusters(data_dir: Path) -> Dict:
    """Download all clusters."""
    print("=" * 70)
    print("CMC Cluster Catalog - Direct Download")
    print("=" * 70)
    
    results = {}
    for cluster_id in CLUSTER_MODELS:
        result = download_cluster(cluster_id, data_dir)
        results[cluster_id] = result
        time.sleep(2)
    
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    
    for cluster_id, result in results.items():
        status = "✓" if result["success"] else "✗"
        print(f"  {status} {cluster_id}")
        if result["success"] and "files_extracted" in result:
            print(f"     Files: {len(result['files_extracted'])}")
    
    return results


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    data_dir = repo_root / "data" / "cmc"
    download_all_clusters(data_dir)
