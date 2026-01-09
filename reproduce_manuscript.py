#!/usr/bin/env python3
"""
TEP-COS Manuscript Reproduction Script
======================================

This script orchestrates the full reproduction of the results presented in
"The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars"
(Smawfield 2026).

It runs the analysis steps in the order presented in the manuscript.
"""

import subprocess
import sys
import time
from pathlib import Path

# Configuration
PROJECT_ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts" / "steps"
LOG_DIR = PROJECT_ROOT / "logs" / "reproduction"

# Ensure log directory exists
LOG_DIR.mkdir(parents=True, exist_ok=True)

def run_step(script_name, description):
    """Run a single analysis step."""
    script_path = SCRIPTS_DIR / script_name
    log_path = LOG_DIR / f"{script_path.stem}.log"
    
    print(f"\n{'='*80}")
    print(f"RUNNING: {description}")
    print(f"SCRIPT:  {script_name}")
    print(f"LOG:     {log_path}")
    print(f"{'='*80}\n")
    
    start_time = time.time()
    
    if not script_path.exists():
        print(f"ERROR: Script not found: {script_path}")
        return False
        
    with open(log_path, "w") as log_file:
        try:
            # Run with unbuffered output to capture progress
            process = subprocess.Popen(
                [sys.executable, str(script_path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )
            
            # Stream output to both console and log
            for line in process.stdout:
                print(line, end="")
                log_file.write(line)
                
            process.wait()
            
            duration = time.time() - start_time
            
            if process.returncode == 0:
                print(f"\n✓ SUCCESS ({duration:.1f}s)")
                return True
            else:
                print(f"\n✗ FAILED (Exit Code: {process.returncode})")
                return False
                
        except Exception as e:
            print(f"\n✗ EXCEPTION: {e}")
            return False

def main():
    print("TEP-COS MANUSCRIPT REPRODUCTION PIPELINE")
    print("========================================")
    print(f"Root: {PROJECT_ROOT}")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Guide: See README.md for details")
    
    # --- HELPER STEPS ---
    # Convert any raw data formats if needed
    run_step("convert_cds_to_rdb.py", "Data Formatting: Convert CDS light curves to RDB")
    
    # --- SECTION 3: PULSAR TIMING ---
    steps_pulsar = [
        ("step_5_27_hybrid_maximum_analysis.py", "Sec 3.2: Construct Maximal Pulsar Sample"),
        ("step_5_10_pulsar_population_controls.py", "Sec 3.3: Population Controls"),
        ("step_5_31_per_cluster_controlled_residuals.py", "Sec 3.4: Density Scaling Test"),
        ("step_5_13_cluster_acceleration_figure.py", "Sec 3.5: Newtonian Baseline Figure (Fig 4.5)"),
        ("step_5_32_density_scaling_figure.py", "Sec 3.5: Density Scaling Figure (Fig 4.6)"),
        ("step_5_11_binary_pulsar_analysis.py", "Sec 3.7: Binary vs Isolated (Cluster)"),
        ("step_5_11_binary_spatial_figure.py", "Sec 3.7: Binary Spatial Figure"),
        ("step_5_12_field_binary_analysis.py", "Sec 3.8: Field Binary Control"),
        ("step_5_9_freire_gcpsr_radial_analysis.py", "Sec 3.10: Radial Diagnostics"),
    ]
    
    # --- SECTION 4: GRAVITATIONAL LENSING ---
    steps_lensing = [
        ("step_3_0_cosmograil_temporal_shear.py", "Sec 4.2: Temporal Shear Analysis (Main)"),
        ("step_3_0_temporal_shear_figure.py", "Sec 4.2: Temporal Shear Figure"),
        ("step_3_2_cosmograil_validation.py", "Sec 4.3: Validation & Injection-Recovery"),
        ("step_3_2_microlensing_figure.py", "Sec 4.3: Microlensing Comparison Figure"),
        ("step_3_10_instrumental_consistency.py", "Sec 4.4: Instrumental Consistency"),
        ("step_3_16_j1004_analysis.py", "Sec 4.5: High-z Cluster Lens (SDSS J1004)"),
        ("step_4_0_lensing_summary_figure.py", "Sec 4: Lensing Summary Figure"),
    ]
    
    # --- SECTION 5: SYNTHESIS ---
    steps_synthesis = [
        ("step_5_40_tep_summary_figure.py", "Sec 5: TEP Cosmology Summary Figure"),
    ]
    
    # --- APPENDIX ---
    steps_appendix = [
        ("step_7_0_sn_ia_stretch_test.py", "App A: SN Ia Stretch vs Host Dispersion"),
    ]
    
    all_steps = steps_pulsar + steps_lensing + steps_appendix + steps_synthesis
    
    failures = []
    
    for script, desc in all_steps:
        success = run_step(script, desc)
        if not success:
            failures.append(script)
    
    print("\n" + "="*80)
    print("REPRODUCTION COMPLETE")
    print("="*80)
    
    if failures:
        print(f"✗ The following {len(failures)} steps FAILED:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print("✓ All analysis steps completed successfully.")
        sys.exit(0)

if __name__ == "__main__":
    main()
