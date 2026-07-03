#!/usr/bin/env python3
"""
TEP-COS Analysis Pipeline Master Script (run_pipeline.py)
===========================================================
Orchestrates the full analysis pipeline for Paper 2:
"The Temporal Equivalence Principle: Cosmological Tests"

This script serves as the central controller for the TEP-COS analysis.
It executes the scientific workflow in a strictly ordered sequence, ensuring
data integrity and dependency management between steps.

Workflow Steps:
1.  **Data Ingestion**: Downloads raw data (Freire GCpsr, ATNF psrcat),
    reconstructs catalogs, and prepares pulsar samples.
2.  **Population Controls**: Matches GC and Field MSPs by period and magnetic field proxy.
3.  **Density Scaling**: Tests suppressed density scaling against Newtonian predictions.
4.  **Statistical Validation**: Covariance-aware tests, LOOCV, bootstrap, permutation tests.
5.  **Binary Analysis**: Compares binary vs isolated MSPs in GCs and field.
6.  **Sensitivity Analysis**: Tests robustness to rho_intra assumption.
7.  **Power Analysis**: Validates statistical power of differential tests.
8.  **Monte Carlo Validation**: Validates Type I error, power, and bias.

Usage:
    python run_pipeline.py
    python run_pipeline.py --skip-validation  # Skip long validation steps

Author: Matthew Lukin Smawfield
Date: March 2026
"""

import argparse
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

# Ensure project root is in path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(PROJECT_ROOT))

from scripts.steps.step_00_data_acquisition import check_data_health, ensure_data
from scripts.utils.logger import TEPLogger, print_status, print_table, set_step_logger

steps_dir = PROJECT_ROOT / "scripts" / "steps"


def check_and_acquire_data(pipeline_logger, skip_prompt=False):
    """Check data availability and acquire missing data."""
    print_status("=" * 80, "TITLE")
    print_status("DATA ACQUISITION PHASE", "TITLE")
    print_status("=" * 80, "TITLE")

    # Run data acquisition
    results = ensure_data("all", verbose=True)

    # Check for critical failures
    critical = ["freire_gcpsr", "atnf_psrcat"]
    missing_critical = [k for k in critical if not results.get(k)]

    if missing_critical:
        print_status(f"✗ CRITICAL DATA MISSING: {', '.join(missing_critical)}", "ERROR")
        print_status("Pipeline cannot proceed without pulsar data.", "ERROR")
        return False

    # Warn about non-critical missing data
    if not results.get("pantheon_plus"):
        print_status(
            "⚠ Pantheon+ data missing - supernova steps will be skipped", "WARNING"
        )

    print_status("✓ Data acquisition phase complete", "SUCCESS")
    return True


def run_step(script_with_args, description, logs_dir, pipeline_logger):
    """Execute a single pipeline step with proper logging."""
    # Parse script name and arguments
    parts = script_with_args.split()
    script_name = parts[0]
    script_args = parts[1:] if len(parts) > 1 else []

    script_path = steps_dir / script_name
    log_path = logs_dir / f"{script_path.stem}.log"

    print_status(f">>> {description}", "TITLE")
    print_status(f"Script: {script_name}", "INFO")
    print_status(f"Log: {log_path}", "INFO")

    start_time = time.time()

    if not script_path.exists():
        print_status(f"ERROR: Script not found: {script_path}", "ERROR")
        return False

    # Setup step-specific logger (but keep pipeline_logger as the active one)
    step_logger = TEPLogger(f"step_{script_path.stem}", log_file_path=log_path)

    try:
        # Setup environment with proper PYTHONPATH
        env = os.environ.copy()
        env["PYTHONPATH"] = (
            str(PROJECT_ROOT / "scripts") + os.pathsep + env.get("PYTHONPATH", "")
        )

        # Run with unbuffered output
        cmd = [sys.executable, str(script_path)] + script_args
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=str(PROJECT_ROOT),
            env=env,
        )

        # Stream output to console, step log, AND master log
        try:
            for line in process.stdout:
                try:
                    print(line, end="")
                except BrokenPipeError:
                    # Stdout pipe closed (e.g., user pressed Ctrl+C or pager quit)
                    pass
                step_logger.info(line.rstrip())
                pipeline_logger.info(line.rstrip())  # Also log to master
        except BrokenPipeError:
            # Handle case where process.stdout itself is closed
            pass

        process.wait()
        duration = time.time() - start_time

        if process.returncode == 0:
            print_status(f"✓ SUCCESS ({duration:.1f}s)", "SUCCESS")
            return True
        else:
            print_status(f"✗ FAILED (Exit Code: {process.returncode})", "ERROR")
            return False

    except Exception as e:
        print_status(f"✗ EXCEPTION: {e}", "ERROR")
        step_logger.error(traceback.format_exc())
        pipeline_logger.error(traceback.format_exc())
        return False


def run_pipeline():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument(
        "--skip-validation",
        action="store_true",
        help="Skip long validation steps (rho_sensitivity, power_analysis, monte_carlo)",
    )
    ap.add_argument(
        "--skip-figures", action="store_true", help="Skip figure generation steps"
    )
    ap.add_argument(
        "--parallel",
        action="store_true",
        help="Enable parallel processing using all available CPU cores (M4 Pro optimized)",
    )
    args = ap.parse_args()

    # Setup Global Logger
    logs_dir = PROJECT_ROOT / "logs"
    logs_dir.mkdir(exist_ok=True)

    pipeline_logger = TEPLogger(
        "pipeline_master", log_file_path=logs_dir / "pipeline_master.log"
    )
    set_step_logger(pipeline_logger)

    print_status("TEP-COS ANALYSIS PIPELINE INITIATED", "TITLE")
    print_status(f"Project Root: {PROJECT_ROOT}", "INFO")
    print_status(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}", "INFO")
    if args.parallel:
        import multiprocessing

        n_cores = multiprocessing.cpu_count()
        print_status(f"Parallel Mode: ENABLED ({n_cores} CPU cores)", "SUCCESS")
    print_status("Starting execution sequence...", "INFO")

    # Data acquisition phase
    if not check_and_acquire_data(pipeline_logger):
        print_status("\n✗ Pipeline aborted: Critical data unavailable", "ERROR")
        sys.exit(1)

    start_time = time.time()
    step_times = {}
    failures = []

    # Define all steps
    steps_data = [
        (
            "step_05_freire_gcpsr_radial_analysis.py",
            "Data Ingestion: Freire GC Pulsar Radial Analysis",
        ),
        (
            "step_02_pulsar_population_controls.py",
            "Data Ingestion: Population Controls",
        ),
        (
            "step_03_field_control_purity_audit.py",
            "Data Ingestion: Field Control Purity Audit",
        ),
    ]

    steps_core_pulsar = [
        ("step_06_hybrid_maximum_analysis.py", "Core: Hybrid Maximum Analysis"),
        (
            "step_07_per_cluster_controlled_residuals.py",
            "Core: Per-Cluster Controlled Residuals",
        ),
        ("step_04_pm_corrected_controls.py", "Core: PM-Corrected Field Controls"),
        ("step_11_full_density_scaling.py", "Core: Full Density Scaling Simulation"),
        (
            "step_12_hierarchical_density_scaling.py",
            "Core: Hierarchical Mixed-Effects Density Scaling",
        ),
        (
            "step_13_covariance_validation.py",
            "Core: Covariance-Aware Statistical Validation",
        ),
        (
            "step_08_signed_pdot_analysis.py",
            "Core: Signed P-dot Analysis",
        ),
        (
            "step_09_pdot_over_p_analysis.py",
            "Core: Pdot/P Parallel Analysis",
        ),
        (
            "step_10_signed_pdot_cmc_comparison.py",
            "Core: Signed Pdot CMC Comparison",
        ),
    ]

    steps_binary = [
        (
            "step_15_binary_pulsar_analysis.py",
            "Binary: GC Binary vs Isolated Analysis",
        ),
        (
            "step_16_binary_screening_model.py",
            "Binary: Temporal Topology Competition Model",
        ),
        ("step_17_curvature_verification.py", "Binary: Curvature Verification"),
        ("step_18_binary_model_validation.py", "Binary: Model Validation"),
        ("step_20_field_binary_analysis.py", "Binary: Field Binary Control"),
        (
            "step_21_integrated_binary_control.py",
            "Binary: Integrated Binary Control Test",
        ),
    ]

    steps_literature = [
        (
            "step_14_cmc_literature_comparison.py",
            "Literature: CMC Literature Consensus",
        ),
    ]

    steps_validation = [
        (
            "step_22_outlier_exclusion_sensitivity.py",
            "Validation: Outlier Exclusion Sensitivity",
        ),
        ("step_23_shklovskii_sensitivity.py", "Validation: Shklovskii Sensitivity"),
        ("step_25_rho_sensitivity.py", "Validation: Rho_intra Sensitivity Analysis"),
        ("step_26_power_analysis.py", "Validation: Power Analysis"),
        ("step_27_monte_carlo_validation.py", "Validation: Monte Carlo Validation"),
        (
            "step_28_sensitivity_cmc_report.py",
            "Validation: Sensitivity & CMC Comparison Report",
        ),
        (
            "step_24_equal_cluster_weighting.py",
            "Validation: Equal Cluster Weighting",
        ),
        ("step_50_injection_recovery.py", "Validation: Injection-Recovery Test"),
        ("step_51_matching_leakage_audit.py", "Validation: Matching Leakage Audit"),
        ("step_52_signed_observable_battery.py", "Validation: Signed Observable Battery"),
        ("step_54_environment_axis_scan.py", "Validation: Environment Axis Scan"),
        ("step_55_cluster_bootstrap.py", "Validation: Cluster Bootstrap"),
    ]

    steps_nbody_pushback = [
        ("step_32_download_cmc_data.py --all", "N-Body: Download CMC Data (All Clusters)"),
        ("step_37_cmc_gold_standard_analysis.py", "N-Body: CMC Gold Standard Test"),
        (
            "step_29_pulsar_dynamical_calibration.py",
            "N-Body Pushback: Dynamical Calibration",
        ),
        ("step_30_sensitivity_analysis.py", "N-Body Pushback: Sensitivity Analysis"),
        ("step_31_cmc_real_comparison.py", "N-Body Pushback: CMC Real Comparison"),
        (
            "step_33_theoretical_uncertainty.py",
            "N-Body Pushback: Theoretical Uncertainty",
        ),
        ("step_34_bayesian_posterior.py", "N-Body Pushback: Bayesian Posterior"),
        # NOTE: step_56 removed - insufficient MSPs with radial positions for meaningful analysis
        ("step_35_core_collapse_test.py", "N-Body Pushback: Core Collapse Test"),
        (
            "step_36_systematic_ceiling.py",
            "N-Body Pushback: Systematic Ceiling Analysis",
        ),
        (
            "step_41_exotic_physics_quantification.py",
            "N-Body Pushback: Exotic Physics Quantification",
        ),
        (
            "step_42_gamma_parameter_free_derivation.py",
            "Theory: Parameter-Free Γ Derivation",
        ),
        (
            "step_43_exact_gamma_derivation_proof.py",
            "Theory: Exact Gamma Derivation Proof",
        ),
        (
            "step_44_kappa_msp_prior.py",
            "Theory: κ_MSP Cross-Domain Prior Export",
        ),
        (
            "step_38_cmc_binary_forensic_audit.py",
            "N-Body Pushback: CMC Binary Forensic Audit",
        ),
        (
            "step_39_cmc_period_sensitivity.py",
            "N-Body Pushback: CMC Period Sensitivity",
        ),
        (
            "step_40_cmc_uncertainty_stack.py",
            "N-Body Pushback: CMC Uncertainty Stack",
        ),
        (
            "step_49_pta_mock_observation.py",
            "N-Body Pushback: PTA Mock Observation",
        ),
        (
            "step_19_derived_shielding_validation.py",
            "N-Body Pushback: Derived Shielding Validation",
        ),
    ]

    steps_figures = [
        ("step_45_binary_spatial_figure.py", "Figure: Binary Spatial Distribution"),
        ("step_47_density_scaling_figure.py", "Figure: Density Scaling"),
        ("step_46_cluster_acceleration_figure.py", "Figure: Cluster Acceleration"),
        ("step_48_tep_summary_figure.py", "Figure: TEP Summary"),
    ]

    steps_appendix = [
        ("step_62_sn_ia_stretch_test.py --fast", "Appendix: SN Ia σ-mB Correlation"),
        ("step_63_sn_ia_robustness.py", "Appendix: SN Ia Robustness Validation"),
        ("step_64_tep_vs_mass_step.py", "Appendix: TEP vs Mass Step Discrimination"),
        ("step_65_sn_ia_audit.py", "Appendix: SN Ia Deep Audit"),
        ("step_59_manga_spatially_resolved.py", "Appendix: MaNGA Age Gradients (Archived)"),
        ("step_60_manga_test_e_age_discrepancy.py", "Appendix: MaNGA LW vs MW Age Discrepancy (Archived)"),
        ("step_61_sdss_test_cc_manganese_clock.py", "Appendix: SDSS Mn Clock Test"),
    ]

    # Build execution list based on arguments
    all_steps = steps_data.copy()
    all_steps.extend(steps_core_pulsar)
    all_steps.extend(steps_binary)
    all_steps.extend(steps_literature)

    if not args.skip_validation:
        all_steps.extend(steps_validation)
        all_steps.extend(steps_nbody_pushback)

    if not args.skip_figures:
        all_steps.extend(steps_figures)

    all_steps.extend(steps_appendix)

    # Execute steps
    for script, desc in all_steps:
        t0 = time.time()
        success = run_step(script, desc, logs_dir, pipeline_logger)
        step_times[desc] = time.time() - t0

        if not success:
            failures.append((script, desc))
            print_status(f"Step failed: {desc}", "ERROR")
            print_status("Pipeline halted due to step failure.", "ERROR")
            sys.exit(1)

    # Final Summary
    total_time = time.time() - start_time

    print_status("=" * 80, "TITLE")
    print_status("PIPELINE EXECUTION SUMMARY", "TITLE")
    print_status("=" * 80, "TITLE")

    # Execution Times Table
    headers = ["Step", "Duration (s)", "Status"]
    rows = []
    for step_name, duration in step_times.items():
        status = "COMPLETED"
        for fail_script, fail_desc in failures:
            if fail_desc == step_name:
                status = "FAILED"
                break
        rows.append([step_name[:40], f"{duration:.2f}", status])

    if rows:
        print_table(headers, rows, title="Execution Timing")

    # Report failures
    if failures:
        print_status(f"\n✗ The following {len(failures)} steps FAILED:", "ERROR")
        for script, desc in failures:
            print_status(f"  - {desc} ({script})", "ERROR")

    print_status(f"\nTotal Execution Time: {total_time:.2f} seconds", "SUCCESS")
    print_status(f"Results Directory: {PROJECT_ROOT}/results/outputs/", "INFO")
    print_status(f"Logs Directory: {logs_dir}/", "INFO")

    if failures:
        print_status("\nPipeline finished with ERRORS.", "WARNING")
        sys.exit(1)
    else:
        print_status("\n✓ Pipeline finished successfully.", "SUCCESS")
        sys.exit(0)


if __name__ == "__main__":
    run_pipeline()
