#!/usr/bin/env python3
"""
TEP-COS Manuscript Reproduction Script (Legacy)
================================================

DEPRECATED: This script is kept for backward compatibility.
Please use run_pipeline.py instead, which provides:
  - Better logging (master log + per-step logs)
  - Command-line arguments (--skip-validation, --only-core, etc.)
  - Execution timing summary
  - All new validation steps

This script now simply calls run_pipeline.py for compatibility.
"""

import sys
import subprocess
from pathlib import Path

if __name__ == "__main__":
    print("=" * 80)
    print("DEPRECATED: reproduce_manuscript.py")
    print("Please use: python run_pipeline.py")
    print("=" * 80)
    print()
    
    # Forward to scripts/run_pipeline.py
    pipeline_script = Path(__file__).parent / "scripts" / "run_pipeline.py"
    result = subprocess.run([sys.executable, str(pipeline_script)] + sys.argv[1:])
    sys.exit(result.returncode)
