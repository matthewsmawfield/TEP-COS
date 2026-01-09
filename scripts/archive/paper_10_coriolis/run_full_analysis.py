#!/usr/bin/env python3
"""TEP-COS Pipeline Orchestrator

Runs the Cosmic Coriolis pipeline end-to-end:
- Step 1.0: Acquire MaNGA dapall + optional MAPS subset
- Step 2.0: Compute per-galaxy asymmetry proxy and test CMB-axis dipole

Author: Matthew Lukin Smawfield
License: CC-BY-4.0
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List

PROJECT_ROOT = Path(__file__).resolve().parent
STEP_DIR = PROJECT_ROOT / "scripts" / "steps"

STEP_1_0 = STEP_DIR / "step_1_0_data_acquisition.py"
STEP_2_0 = STEP_DIR / "step_2_0_cosmic_coriolis_analysis.py"
STEP_2_1 = STEP_DIR / "step_2_1_distance_residual_dipole.py"


def run_step(label: str, cmd: List[str]) -> None:
    print("\n" + "=" * 80)
    print(f"{datetime.now():%Y-%m-%d %H:%M:%S} | {label}")
    print("=" * 80)
    subprocess.run(cmd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Step 1.0 + Step 2.0 Cosmic Coriolis pipeline")

    parser.add_argument("--skip-step1", action="store_true")
    parser.add_argument("--skip-step2", action="store_true")
    parser.add_argument("--skip-step2-1", action="store_true")
    parser.add_argument("--python", default=sys.executable)

    parser.add_argument("--daptype", default="HYB10-MILESHC-MASTARSSP")
    parser.add_argument("--drpver", default="v3_1_1")
    parser.add_argument("--dapver", default="3.1.0")
    parser.add_argument("--sample-n", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)

    parser.add_argument(
        "--stratify-xcmb",
        action="store_true",
        help="When selecting from dapall in Step 1, stratify sample to balance x_CMB bins (improves dipole sensitivity at fixed N).",
    )
    parser.add_argument("--xcmb-bins", type=int, default=10)

    parser.add_argument(
        "--no-size-estimate",
        action="store_true",
        help="Disable HTTP HEAD size estimates for Step 1 downloads.",
    )
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=2.0,
        help="Progress log interval (seconds) for Step 1 downloads.",
    )
    parser.add_argument(
        "--chunk-mb",
        type=float,
        default=1.0,
        help="Download chunk size (MiB) for Step 1 downloads.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing downloaded files in Step 1.",
    )

    parser.add_argument(
        "--download-drpall",
        action="store_true",
        help="Download drpall catalog in Step 1 (required for Step 2.1 distance-residual dipole).",
    )

    parser.add_argument(
        "--download-method",
        choices=["http", "rsync"],
        default="http",
        help="Step 1 download method. SDSS recommends rsync (preferably via dtn.sdss.org) for bulk transfers.",
    )
    parser.add_argument(
        "--rsync-host",
        default="dtn.sdss.org",
        help="Rsync host used by Step 1 when --download-method rsync.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=2,
        help="Max concurrent downloads in Step 1 (keep small to avoid SDSS throttling).",
    )
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-backoff", type=float, default=2.0)

    parser.add_argument("--velocity-source", choices=["stellar", "gas"], default="stellar")
    parser.add_argument("--gas-line", default="Ha-6564")

    parser.add_argument("--rre-min", type=float, default=0.8)
    parser.add_argument("--rre-max", type=float, default=1.2)
    parser.add_argument("--rre-sys-max", type=float, default=0.1)
    parser.add_argument("--delta-phi-deg", type=float, default=20.0)

    parser.add_argument("--tf-q0", type=float, default=0.2)
    parser.add_argument("--tf-ba-min", type=float, default=0.2)
    parser.add_argument("--tf-ba-max", type=float, default=0.85)

    parser.add_argument(
        "--output-tag",
        type=str,
        default=None,
        help="Optional tag appended to Step 2 outputs (CSV/JSON/MD and figures) to prevent overwrites across reruns.",
    )
    parser.add_argument(
        "--min-galaxies",
        type=int,
        default=20,
        help="Minimum galaxies required to run the Step 2 dipole fit/permutation test (below this, Step 2 writes tables/reports and exits cleanly).",
    )
    parser.add_argument("--n-perm", type=int, default=2000)
    parser.add_argument("--n-axis-rand", type=int, default=2000, help="Number of random axes for look-elsewhere control")
    parser.add_argument("--use-norm", action="store_true")

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not STEP_1_0.exists():
        raise FileNotFoundError(f"Missing: {STEP_1_0}")
    if not STEP_2_0.exists():
        raise FileNotFoundError(f"Missing: {STEP_2_0}")
    if not args.skip_step2_1 and not STEP_2_1.exists():
        raise FileNotFoundError(f"Missing: {STEP_2_1}")

    if not args.skip_step1:
        print("\n" + "-" * 80)
        print("Pipeline preflight")
        print("-" * 80)
        print(f"DRPVER={args.drpver} | DAPVER={args.dapver} | DAPTYPE={args.daptype}")
        print(f"Sample size (MAPS): N={args.sample_n} | Seed={args.seed}")
        print(
            "Step 1 downloads to: data/dapall/ and data/maps/<daptype>/... (see Step 1 logs for full plan)"
        )
        print("-" * 80)

    if not args.skip_step1:
        run_step(
            "STEP 1.0 | Data Acquisition",
            [
                args.python,
                str(STEP_1_0),
                "--download-dapall",
                "--select-from-dapall",
                "--download-maps",
                "--daptype",
                args.daptype,
                "--drpver",
                args.drpver,
                "--dapver",
                args.dapver,
                "--sample-n",
                str(args.sample_n),
                "--seed",
                str(args.seed),
                "--cmb-ra-deg",
                "168.0",
                "--cmb-dec-deg",
                "-7.0",
                "--xcmb-bins",
                str(args.xcmb_bins),
                "--progress-interval",
                str(args.progress_interval),
                "--chunk-mb",
                str(args.chunk_mb),
                "--download-method",
                args.download_method,
                "--rsync-host",
                args.rsync_host,
                "--max-workers",
                str(args.max_workers),
                "--retries",
                str(args.retries),
                "--retry-backoff",
                str(args.retry_backoff),
            ]
            + (["--stratify-xcmb"] if args.stratify_xcmb else [])
            + (["--download-drpall"] if args.download_drpall else [])
            + (["--no-size-estimate"] if args.no_size_estimate else [])
            + (["--overwrite"] if args.overwrite else []),
        )

    if args.skip_step2:
        print("Skipping Step 2.0 (per --skip-step2)")
        return

    run_step(
        "STEP 2.0 | Cosmic Coriolis Analysis",
        [
            args.python,
            str(STEP_2_0),
            "--daptype",
            args.daptype,
            "--velocity-source",
            args.velocity_source,
            "--gas-line",
            args.gas_line,
            "--rre-min",
            str(args.rre_min),
            "--rre-max",
            str(args.rre_max),
            "--rre-sys-max",
            str(args.rre_sys_max),
            "--delta-phi-deg",
            str(args.delta_phi_deg),
            "--min-galaxies",
            str(args.min_galaxies),
            "--n-perm",
            str(args.n_perm),
        ]
        + (["--use-norm"] if args.use_norm else [])
        + (["--output-tag", args.output_tag] if args.output_tag else []),
    )

    if args.skip_step2_1:
        print("Skipping Step 2.1 (per --skip-step2-1)")
        print("\nPipeline complete")
        return

    run_step(
        "STEP 2.1 | Distance-Residual Dipole",
        [
            args.python,
            str(STEP_2_1),
            "--daptype",
            args.daptype,
            "--drpver",
            args.drpver,
            "--velocity-source",
            "stellar",
            "--rre-min",
            str(args.rre_min),
            "--rre-max",
            str(args.rre_max),
            "--rre-sys-max",
            str(args.rre_sys_max),
            "--delta-phi-deg",
            str(args.delta_phi_deg),
            "--tf-q0",
            str(args.tf_q0),
            "--tf-ba-min",
            str(args.tf_ba_min),
            "--tf-ba-max",
            str(args.tf_ba_max),
            "--n-perm",
            str(args.n_perm),
            "--seed",
            str(args.seed),
            "--min-galaxies",
            str(args.min_galaxies),
            "--n-axis-rand",
            str(args.n_axis_rand),
        ]
        + (["--output-tag", args.output_tag] if args.output_tag else []),
    )

    print("\nPipeline complete")


if __name__ == "__main__":
    main()
