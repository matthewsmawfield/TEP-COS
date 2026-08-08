# TEP-COS Pipeline Scripts

## Overview

The pipeline is orchestrated by `run_pipeline.py`, which executes analysis steps in strict dependency order. Each step writes outputs to `results/outputs/` (JSON, CSV, Markdown) and figures to `results/figures/` (PDF, PNG).

## Running the Pipeline

```bash
# Full pipeline
python3 scripts/run_pipeline.py

# Skip long validation steps
python3 scripts/run_pipeline.py --skip-validation

# Skip figure generation
python3 scripts/run_pipeline.py --skip-figures

# Enable parallel processing
python3 scripts/run_pipeline.py --parallel
```

## Step Structure

Steps are numbered `step_00` through `step_65` and grouped into phases:

| Phase | Steps | Description |
|-------|-------|-------------|
| Data Ingestion | 00–05 | Download ATNF/Freire catalogs, CMC data; build population controls |
| Core Pulsar | 06–13 | Hybrid sample, per-cluster residuals, density scaling, covariance validation |
| Binary | 15–21 | Binary vs isolated analysis, screening model, field controls |
| Literature | 14 | CMC literature consensus meta-analysis |
| Validation | 22–28, 50–55 | Outlier exclusion, Shklovskii, rho sensitivity, power, Monte Carlo, injection-recovery, bootstrap |
| N-Body Pushback | 29–44 | Dynamical calibration, CMC gold standard, uncertainty stack, PTA mock, kappa prior |
| Figures | 45–48 | Binary spatial, density scaling, acceleration, TEP summary |
| Appendix | 59–65 | MaNGA, SDSS, SN Ia stretch/robustness, mass-step discrimination |

## Building the Manuscript

```bash
# Build markdown from HTML components
python3 scripts/build_manuscript.py

# Build static site (HTML + markdown)
cd site && npm run build
```

## Key Outputs

| Output | Step | Key Numbers |
|--------|------|-------------|
| `step_06_hybrid_maximum_analysis.json` | 06 | 199 GC + 351 field, 0.63 dex raw, 0.40 dex controlled |
| `step_12_hierarchical_density_results.json` | 12 | Γ = 0.39 ± 0.08 |
| `step_13_covariance_validation.json` | 13 | 5.6σ covariance-aware significance |
| `step_14_cmc_literature.json` | 14 | Γ_N = 0.75, 4.1σ rejection |
| `step_40_cmc_uncertainty_stack.json` | 40 | 12.7σ nominal, 3.6σ conservative |
| `step_37_cmc_gold_standard.json` | 37 | 20 clusters, 18,813 MSPs |
