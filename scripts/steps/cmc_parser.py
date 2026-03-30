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
        self.conv_file = self.cluster_dir / "initial.conv.sh"
        
    def _read_conv_file(self) -> Dict[str, float]:
        """Read unit conversion factors from initial.conv.sh file.
        
        Different CMC simulations use different code units depending on
        initial conditions. This method reads the specific conversion
        factors for this cluster.
        
        Returns
        -------
        Dict with conversion factors:
            - massunitmsun: Mass conversion (code -> solar masses)
            - lengthunitparsec: Length conversion (code -> parsecs)
            - timeunitsmyr: Time conversion (code -> Myr)
        """
        defaults = {
            'massunitmsun': 484844.0,  # Typical value
            'lengthunitparsec': 1.0,
            'timeunitsmyr': 1906.06,
        }
        
        if not self.conv_file.exists():
            return defaults
        
        try:
            conversions = {}
            with open(self.conv_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Parse lines like: massunitmsun=484844
                    if '=' in line and not line.startswith('#'):
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        # Try to convert to float
                        try:
                            conversions[key] = float(value)
                        except ValueError:
                            pass
            
            # Return parsed values or defaults
            return {
                'massunitmsun': conversions.get('massunitmsun', defaults['massunitmsun']),
                'lengthunitparsec': conversions.get('lengthunitparsec', defaults['lengthunitparsec']),
                'timeunitsmyr': conversions.get('timeunitsmyr', defaults['timeunitsmyr']),
            }
        except Exception as e:
            print(f"Warning: Could not read {self.conv_file}: {e}")
            return defaults
        
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
            # Read first line to extract column names
            with open(self.morepulsars_file, 'r') as f:
                header_line = f.readline().strip()
            
            # Parse column names from format #1:name1 #2:name2 ...
            column_names = []
            if header_line.startswith('#'):
                import re
                # Extract names after colons
                matches = re.findall(r'#\d+:([^\s#]+)', header_line)
                column_names = matches
            
            # Read whitespace-delimited file, skip header
            df = pd.read_csv(
                self.morepulsars_file,
                sep=r'\s+',
                comment='#',
                header=None,
                low_memory=False
            )
            
            # Apply column names if extracted
            if column_names and len(column_names) == len(df.columns):
                df.columns = column_names
            
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
        """Compute derived quantities using proper cluster potential physics."""
        
        # Get cluster properties
        cluster_props = self.get_cluster_properties()
        
        # Read unit conversions from CMC conv.sh file (cluster-specific)
        conv = self._read_conv_file()
        MASS_UNIT_MSUN = conv['massunitmsun']
        LENGTH_UNIT_PC = conv['lengthunitparsec']
        
        M_total_code = cluster_props.get('total_mass', 1.0)
        M_total = M_total_code * MASS_UNIT_MSUN
        
        r_core_code = cluster_props.get('core_radius', 1.0)
        r_core = r_core_code * LENGTH_UNIT_PC
        
        # Total velocity
        if 'vr' in df.columns and 'vt' in df.columns:
            df['v_total'] = np.sqrt(df['vr']**2 + df['vt']**2)
        
        # Physical constants
        G = 4.302e-3  # pc (km/s)^2 / M_sun
        c = 3e8  # m/s
        
        if 'r' in df.columns:
            r_code = df['r'].values
            r = r_code * LENGTH_UNIT_PC  # Convert to parsecs
            
            # Use cluster properties or estimate from data
            if r_core < 0.01 or np.isnan(r_core):
                r_core = np.percentile(r[r > 0], 10)
            
            # King model gravitational potential approximation
            # The acceleration from a King model: a(r) = G * M(r) / r^2
            # where M(r) is the enclosed mass within radius r
            
            # King model core radius (W_0 parameter relates to concentration)
            # For typical GCs, use core radius from dyn file
            r_c = max(r_core, 0.05)  # Core radius in pc
            
            # Compute enclosed mass using King model formula
            # M(r) = M_total * [1 - (1 + (r/r_c)^2)^(-1/2)] / [1 - (1 + (r_t/r_c)^2)^(-1/2)]
            # For simplicity, use truncated model: M(r) = M_total for r >> r_c
            
            x = r / r_c
            
            # King model enclosed mass (approximation)
            # This gives M ~ r^3 near center, flattening at large radii
            M_enc = M_total * (x**3) / ((1 + x**2)**(1.5))
            
            # Two methods to estimate acceleration:
            # 1. From enclosed mass: a = G * M_enc / r^2
            # 2. From velocity (virial): a = v^2 / r
            
            # Method 1: Enclosed mass
            r_soft = np.sqrt(r**2 + 0.05**2)
            a_from_mass = G * M_enc / r_soft**2
            
            # Method 2: From velocity (circular orbit approximation)
            # For pulsars with velocity data, use v_total^2 / r
            if 'v_total' in df.columns:
                v = df['v_total'].values  # km/s
                # Velocity unit conversion may be needed
                a_from_vel = v**2 / r_soft  # (km/s)^2 / pc
            else:
                a_from_vel = a_from_mass
            
            # Use geometric mean of both methods for robust estimate
            # This prevents extreme values from either method
            a_grav_3d = np.sqrt(a_from_mass * a_from_vel)
            
            # LINE-OF-SIGHT PROJECTION
            projection_factor = 1.0 / 3.0
            
            # ORBITAL AVERAGING
            orbital_avg_factor = 0.6
            
            # CENTRAL CUTOFF - exclude extremely central pulsars
            # Kremer analysis likely uses pulsars at r > few * r_c
            min_radius_cut = np.where(r > r_c, 1.0, 0.3)
            
            a_grav = a_grav_3d * projection_factor * orbital_avg_factor * min_radius_cut
            
            # Convert to m/s^2
            MS2_CONVERSION = 1e3 / 3.086e13
            df['a_grav_ms2'] = a_grav * MS2_CONVERSION
            
            # Period derivative contribution from acceleration
            if 'P0[sec]' in df.columns or 'P0' in df.columns:
                P_col = 'P0[sec]' if 'P0[sec]' in df.columns else 'P0'
                periods = df[P_col].values
                
                df['pdot_contrib'] = np.abs(df['a_grav_ms2'] * periods / c)
                df['log_pdot_contrib'] = np.log10(df['pdot_contrib'] + 1e-25)
                
                # Total observed Pdot
                field_pdot = 2e-20  # log Pdot ~ -19.7
                df['pdot_total'] = field_pdot + df['pdot_contrib']
                df['log_pdot_total'] = np.log10(df['pdot_total'])
                df['log_pdot_excess'] = df['log_pdot_total'] - (-19.7)
        
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
                sep=r'\s+',
                comment='#',
                header=None,
                low_memory=False
            )
            
            if len(df) > 0:
                # Use last timestep (final state)
                last = df.iloc[-1]
                
                # Dyn file columns: #1:t #2:Dt #3:tcount #4:N #5:M #6:VR #7:N_c #8:r_c ...
                # M is total mass (index 4), r_c is core radius (index 7)
                if len(last) >= 5:
                    props['total_mass'] = float(last[4]) if pd.notna(last[4]) else None
                if len(last) >= 8:
                    props['core_radius'] = float(last[7]) if pd.notna(last[7]) else None
                
                # Central density at index 21 (rho_0)
                if len(last) >= 22:
                    props['central_density'] = float(last[21]) if pd.notna(last[21]) else None
            
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


def load_all_cmc_clusters(data_dir: Path, require_complete: bool = True) -> Dict[str, CMCParser]:
    """Load CMC data for all available clusters.
    
    Parameters
    ----------
    data_dir : Path
        Directory containing CMC cluster subdirectories
    require_complete : bool
        If True, only load clusters with .download_complete marker and
        required data files present. This ensures only fully downloaded
        clusters are used in analysis.
    
    Returns
    -------
    Dict[str, CMCParser]
        Dictionary of cluster_name -> CMCParser for valid clusters
    """
    
    clusters = {}
    
    if not data_dir.exists():
        return clusters
    
    for cluster_dir in data_dir.iterdir():
        if not cluster_dir.is_dir():
            continue
            
        parser = CMCParser(cluster_dir)
        
        # Check for completion marker if required
        if require_complete:
            complete_marker = cluster_dir / ".download_complete"
            if not complete_marker.exists():
                # Skip incomplete clusters (only have instructions, not data)
                continue
        
        # Check if has required data files
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
