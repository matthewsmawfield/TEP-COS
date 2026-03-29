#!/usr/bin/env python3
"""
CMC Data Parser - Parse CMC Cluster Catalog Files
====================================================

Parses CMC (Cluster Monte Carlo) output files:
- initial.morepulsars.dat (neutron star properties)
- output.window.snapshot.h5 (HDF5 stellar snapshots)
- initial.dyn.dat (dynamical evolution)

Author: M. Smawfield
Date: March 2026
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List
import warnings

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False
    warnings.warn("h5py not installed. Cannot read CMC HDF5 snapshots.")


class CMCParser:
    """Parser for CMC Cluster Monte Carlo output files."""
    
    def __init__(self, cluster_dir: Path):
        """
        Initialize parser for a specific cluster.
        
        Parameters
        ----------
        cluster_dir : Path
            Directory containing CMC output files
        """
        self.cluster_dir = Path(cluster_dir)
        self.cluster_name = self.cluster_dir.name
        
        # File paths
        self.morepulsars_file = self.cluster_dir / "initial.morepulsars.dat"
        self.snapshot_file = self.cluster_dir / "output.window.snapshot.h5"
        self.dyn_file = self.cluster_dir / "initial.dyn.dat"
        
    def parse_morepulsars(self) -> Optional[pd.DataFrame]:
        """
        Parse initial.morepulsars.dat file.
        
        Returns DataFrame with neutron star properties:
        - id0: ID number
        - m0: Mass [M_sun]
        - B0: Magnetic field [G]
        - P0: Spin period [sec]
        - r: Distance from cluster center [pc]
        - vr: Radial velocity [km/s]
        - vt: Tangential velocity [km/s]
        - binflag: Binary flag (1 if in binary)
        - a: Semi-major axis [AU] (if binary)
        - ecc: Eccentricity (if binary)
        - bacc0: Mass accreted [M_sun]
        - tacc0: Time spent accreting [Myr]
        """
        if not self.morepulsars_file.exists():
            return None
        
        try:
            # Read whitespace-delimited file
            df = pd.read_csv(
                self.morepulsars_file,
                delim_whitespace=True,
                comment='#',
                low_memory=False
            )
            
            # Filter to neutron stars (startype0 == 13)
            if 'startype0' in df.columns:
                df = df[df['startype0'] == 13].copy()
            
            if len(df) == 0:
                return None
            
            # Compute derived quantities
            df = self._compute_derived(df)
            
            return df
            
        except Exception as e:
            print(f"Error parsing {self.morepulsars_file}: {e}")
            return None
    
    def _compute_derived(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute derived quantities for pulsars."""
        
        # Total velocity
        if 'vr' in df.columns and 'vt' in df.columns:
            df['v_total'] = np.sqrt(df['vr']**2 + df['vt']**2)
        
        # Total position (3D)
        if 'r' in df.columns:
            df['r_3d'] = df['r']  # Simplified - CMC gives 1D radial distance
        
        # Compute line-of-sight acceleration (simplified)
        # In full implementation, would use cluster potential
        G = 4.302e-3  # pc (km/s)^2 / M_sun
        
        if 'r' in df.columns:
            # Estimate enclosed mass (simplified)
            r = df['r'].values
            
            # Velocity dispersion as proxy for mass
            if 'v_total' in df.columns:
                sigma = df['v_total'].values
                M_enc = sigma**2 * r / G
                
                # Gravitational acceleration (km/s per pc)
                a_grav = G * M_enc / r**2
                
                # Convert to m/s^2
                pc_to_km = 3.086e13
                df['a_grav_ms2'] = a_grav * 1e3 / pc_to_km
                
                # Period derivative contribution
                c = 3e8  # m/s
                df['pdot_fraction'] = df['a_grav_ms2'] / c
                df['log_pdot_contrib'] = np.log10(df['pdot_fraction'].abs() + 1e-25)
        
        return df
    
    def parse_snapshot(self, snapshot_key: Optional[str] = None) -> Optional[pd.DataFrame]:
        """
        Parse HDF5 snapshot file.
        
        Parameters
        ----------
        snapshot_key : str, optional
            Specific snapshot to read. If None, reads latest.
        
        Returns DataFrame with stellar population.
        """
        if not HAS_H5PY:
            warnings.warn("h5py not available. Cannot parse HDF5 snapshots.")
            return None
        
        if not self.snapshot_file.exists():
            return None
        
        try:
            with h5py.File(self.snapshot_file, 'r') as f:
                # List available snapshots
                keys = list(f.keys())
                
                if not keys:
                    return None
                
                # Select snapshot
                if snapshot_key is None:
                    # Use last snapshot (most evolved)
                    snapshot_key = sorted(keys)[-1]
                
                if snapshot_key not in f:
                    available = ', '.join(keys[:5])
                    print(f"Snapshot {snapshot_key} not found. Available: {available}...")
                    return None
                
                # Read data
                snap = f[snapshot_key]
                
                # Convert to DataFrame
                data = {}
                for col in snap.dtype.names:
                    data[col] = snap[col][:]
                
                df = pd.DataFrame(data)
                
                # Filter to neutron stars (startype == 13)
                if 'startype' in df.columns:
                    df = df[df['startype'] == 13].copy()
                
                return df
                
        except Exception as e:
            print(f"Error reading snapshot {self.snapshot_file}: {e}")
            return None
    
    def get_cluster_properties(self) -> Dict:
        """Extract cluster global properties from dyn file."""
        
        props = {
            'total_mass': None,
            'core_radius': None,
            'half_mass_radius': None,
            'velocity_dispersion': None,
            'central_density': None,
        }
        
        if not self.dyn_file.exists():
            return props
        
        try:
            # Parse dyn file
            # Format: time, mass, rc, rh, etc.
            df = pd.read_csv(
                self.dyn_file,
                delim_whitespace=True,
                comment='#',
                header=None,
                low_memory=False
            )
            
            if len(df) > 0:
                # Use last timestep (final state)
                last = df.iloc[-1]
                
                if len(last) >= 3:
                    props['total_mass'] = float(last[1]) if pd.notna(last[1]) else None
                    props['core_radius'] = float(last[2]) if pd.notna(last[2]) else None
                
                if len(last) >= 5:
                    props['half_mass_radius'] = float(last[4]) if pd.notna(last[4]) else None
            
            return props
            
        except Exception as e:
            print(f"Error parsing {self.dyn_file}: {e}")
            return props
    
    def get_all_pulsars(self) -> Optional[pd.DataFrame]:
        """
        Get all pulsars from both morepulsars and snapshot.
        
        Merges data from both sources for comprehensive pulsar sample.
        """
        # Try morepulsars first (most detailed for NS)
        df_more = self.parse_morepulsars()
        
        if df_more is not None and len(df_more) > 0:
            return df_more
        
        # Fallback to snapshot
        df_snap = self.parse_snapshot()
        
        return df_snap


def load_all_cmc_clusters(data_dir: Path) -> Dict[str, CMCParser]:
    """Load CMC data for all available clusters."""
    
    clusters = {}
    
    if not data_dir.exists():
        return clusters
    
    for cluster_dir in data_dir.iterdir():
        if cluster_dir.is_dir():
            parser = CMCParser(cluster_dir)
            
            # Check if has data
            if parser.morepulsars_file.exists() or parser.snapshot_file.exists():
                clusters[cluster_dir.name] = parser
    
    return clusters


if __name__ == "__main__":
    # Test parsing
    repo_root = Path(__file__).resolve().parents[2]
    data_dir = repo_root / "data" / "cmc"
    
    clusters = load_all_cmc_clusters(data_dir)
    
    print("CMC Data Parser Test")
    print("=" * 60)
    
    for name, parser in clusters.items():
        print(f"\nCluster: {name}")
        print("-" * 40)
        
        # Try to parse
        pulsars = parser.get_all_pulsars()
        
        if pulsars is not None:
            print(f"  Found {len(pulsars)} pulsars")
            print(f"  Columns: {', '.join(pulsars.columns[:5])}...")
        else:
            print("  No pulsar data found")
        
        # Cluster properties
        props = parser.get_cluster_properties()
        if props['total_mass']:
            print(f"  Mass: {props['total_mass']:.2e} M_sun")
        if props['core_radius']:
            print(f"  Core radius: {props['core_radius']:.3f} pc")
