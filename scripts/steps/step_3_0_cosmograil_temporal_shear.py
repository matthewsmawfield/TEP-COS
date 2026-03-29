#!/usr/bin/env python3
"""
Step 3.0: COSMOGRAIL Temporal Shear Analysis
============================================

This script tests the TEP-GL (Temporal Echo Physics - Gravitational Lensing)
prediction that time delays in strongly lensed quasars should exhibit
scale-dependent behavior ("temporal shear") if gravitational lensing includes
a temporal transport component beyond standard General Relativity.

THEORETICAL BACKGROUND
----------------------
In standard GR, the time delay between lensed images is a single constant value
determined by the Shapiro delay and geometric path difference. TEP-GL predicts
that different frequency components of the quasar variability may experience
slightly different effective delays due to temporal transport effects.

KEY OBSERVABLE
--------------
Γ = d(Δt)/d(log τ)  [units: days per decade of timescale]

This is the "temporal shear" - the slope of measured time delay versus the
logarithm of the variability timescale τ being probed.

- Null hypothesis (GR): Γ = 0 (single constant delay at all timescales)
- TEP-GL prediction: Γ ≠ 0 (scale-dependent delay, should be achromatic)

METHODOLOGY OVERVIEW
--------------------
1. Parse COSMOGRAIL light curves from .rdb files
2. Detrend to remove slow microlensing variations (Gaussian smoothing)
3. For each image pair:
   a. Measure broadband delay using cross-correlation
   b. Apply bandpass filters at multiple timescales τ
   c. Measure delay at each τ using mode-locked cross-correlation
   d. Fit linear relation: Δt(τ) = Γ·log₁₀(τ) + intercept
4. Assess statistical significance of Γ ≠ 0

KEY METHODOLOGICAL IMPROVEMENTS (v2.0)
--------------------------------------
- Mode-locking: Constrains multiscale delay search to ±50 days around the
  broadband delay to prevent spurious alias peaks from dominating
- Variance filtering: Rejects timescales where bandpass filter preserves <2%
  of original variance (noise-dominated scales)
- OLS uncertainty: Uses residual-based standard errors instead of overly
  conservative FWHM-based delay uncertainties

Author: TEP Collaboration
"""

import argparse
import fnmatch
import json
import logging
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import interpolate, signal, stats
from scipy.ndimage import gaussian_filter1d
from joblib import Parallel, delayed, cpu_count

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


def print_status(msg: str, level: str = "PROCESS") -> None:
    """Print formatted status message."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] [{level}] {msg}")


def _compute_season_intervals(times: np.ndarray, gap_days: float = 30.0) -> List[Tuple[float, float]]:
    times = np.asarray(times)
    times = times[np.isfinite(times)]
    if times.size == 0:
        return []

    t = np.sort(np.unique(times))
    if t.size == 1:
        return [(float(t[0]), float(t[0]))]

    intervals: List[Tuple[float, float]] = []
    start = float(t[0])
    for i in range(1, t.size):
        if float(t[i] - t[i - 1]) > gap_days:
            intervals.append((start, float(t[i - 1])))
            start = float(t[i])
    intervals.append((start, float(t[-1])))
    return intervals


def _mask_times_in_intervals(times: np.ndarray, intervals: List[Tuple[float, float]]) -> np.ndarray:
    if not intervals:
        return np.zeros_like(times, dtype=bool)
    m = np.zeros_like(times, dtype=bool)
    for a, b in intervals:
        m |= (times >= a) & (times <= b)
    return m


def _gaussian_smooth_nan(y: np.ndarray, sigma: float) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    if y.size == 0:
        return y

    w = np.isfinite(y).astype(float)
    y0 = np.where(np.isfinite(y), y, 0.0)

    if not np.any(w > 0):
        return np.full_like(y, np.nan, dtype=float)

    if sigma < 1:
        sigma = 1.0

    num = gaussian_filter1d(y0, sigma=float(sigma), mode="nearest")
    den = gaussian_filter1d(w, sigma=float(sigma), mode="nearest")

    out = np.full_like(y, np.nan, dtype=float)
    good = den > 1e-8
    out[good] = num[good] / den[good]
    return out


@dataclass
class LightCurve:
    """Single image light curve."""
    label: str
    t: np.ndarray  # MHJD
    mag: np.ndarray
    magerr: np.ndarray
    
    def __post_init__(self):
        # Sort by time and remove non-finite
        valid = np.isfinite(self.t) & np.isfinite(self.mag) & np.isfinite(self.magerr)
        idx = np.argsort(self.t[valid])
        self.t = self.t[valid][idx]
        self.mag = self.mag[valid][idx]
        self.magerr = self.magerr[valid][idx]
    
    @property
    def n_epochs(self) -> int:
        return len(self.t)
    
    @property
    def baseline_days(self) -> float:
        return float(self.t[-1] - self.t[0]) if len(self.t) > 1 else 0.0


@dataclass
class LensSystem:
    """A gravitationally lensed quasar system with multiple images."""
    system_id: str
    light_curves: Dict[str, LightCurve]  # label -> LightCurve
    band: str = "R"
    
    @property
    def image_labels(self) -> List[str]:
        return sorted(self.light_curves.keys())
    
    @property
    def n_images(self) -> int:
        return len(self.light_curves)
    
    def get_image_pairs(self) -> List[Tuple[str, str]]:
        """Return all unique image pairs."""
        labels = self.image_labels
        pairs = []
        for i, l1 in enumerate(labels):
            for l2 in labels[i+1:]:
                pairs.append((l1, l2))
        return pairs


def parse_rdb_file(filepath: Path) -> LensSystem:
    """Parse COSMOGRAIL .rdb file format."""
    system_id = filepath.stem.split("_")[0]
    
    with open(filepath, "r") as f:
        lines = f.readlines()
    
    # Find header line (first non-empty, non-separator line)
    header_line = None
    data_start = 0
    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        if "====" in line or "----" in line:
            data_start = i + 1
            continue
        if header_line is None and ("mhjd" in line.lower() or "mjd" in line.lower()):
            header_line = line
            continue
        if header_line is not None and "====" not in line:
            data_start = i
            break
    
    if header_line is None:
        raise ValueError(f"Could not find header in {filepath}")
    
    # Parse header to find column structure
    # Format: mhjd mag_A magerr_A mag_B magerr_B ...
    header_parts = header_line.split()
    
    # Find image labels from header
    image_labels = []
    col_map = {}  # label -> (mag_col, err_col)
    
    for j, col in enumerate(header_parts):
        if col.lower().startswith("mag_") and "err" not in col.lower():
            label = col.split("_")[1]
            image_labels.append(label)
            col_map[label] = {"mag": j}
        elif col.lower().startswith("magerr_"):
            label = col.split("_")[1]
            if label in col_map:
                col_map[label]["err"] = j
    
    # Parse data lines
    data = {label: {"t": [], "mag": [], "err": []} for label in image_labels}
    
    for line in lines[data_start:]:
        line = line.strip()
        if not line or "====" in line:
            continue
        
        parts = line.split()
        if len(parts) < 2:
            continue
        
        try:
            t = float(parts[0])
        except ValueError:
            continue
        
        for label in image_labels:
            if label not in col_map:
                continue
            mag_col = col_map[label].get("mag")
            err_col = col_map[label].get("err")
            
            if mag_col is None or mag_col >= len(parts):
                continue
            
            try:
                mag = float(parts[mag_col])
                err = float(parts[err_col]) if err_col and err_col < len(parts) else 0.01
                data[label]["t"].append(t)
                data[label]["mag"].append(mag)
                data[label]["err"].append(err)
            except (ValueError, IndexError):
                continue
    
    # Build LightCurve objects
    light_curves = {}
    for label in image_labels:
        if len(data[label]["t"]) > 10:  # Require at least 10 epochs
            light_curves[label] = LightCurve(
                label=label,
                t=np.array(data[label]["t"]),
                mag=np.array(data[label]["mag"]),
                magerr=np.array(data[label]["err"]),
            )
    
    return LensSystem(system_id=system_id, light_curves=light_curves)


def parse_multiband_csv(filepath: Path) -> Optional[LensSystem]:
    """
    Parse multi-band CSV files from COSMOGRAIL/GLITP archives.
    
    Supports several CSV formats:
    - HE0435 format: MJD,RmagA,e_RmagA,RmagB,e_RmagB,...,Band
    - HE1104 format: HJD,A,e_A,B,e_B,Filt
    - Q2237 format: MJD,mA,e_mA,mB,e_mB,...,mbeta,e_mbeta
    
    Returns None if parsing fails or insufficient data.
    """
    import csv
    
    fname = filepath.stem.lower()
    
    # Extract system ID and band from filename
    # Patterns: he0435_JAA703A250_R.csv, he1104_JApJ798_95_B.csv, q2237_JAA637A89_g.csv
    parts = fname.split("_")
    if len(parts) < 2:
        return None
    
    system_id = parts[0].upper()
    band = parts[-1].upper() if len(parts[-1]) <= 2 else "R"
    
    try:
        with open(filepath, "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = list(reader)
    except Exception:
        return None
    
    if len(rows) < 15:
        return None
    
    # Normalize header to lowercase for matching
    header_lower = [h.lower().strip() for h in header]
    
    # Find time column
    time_col = None
    time_format = None
    for i, h in enumerate(header_lower):
        if h in ("mjd", "mhjd", "hjd", "jd"):
            time_col = i
            time_format = "mjd"
            break
        if h in ("obs.date", "date", "obs_date", "date_obs"):
            time_col = i
            time_format = "calendar"
            break
    
    if time_col is None:
        return None
    
    # Detect image columns based on header patterns
    # Pattern 1: RmagA, e_RmagA (HE0435 style)
    # Pattern 2: A, e_A (HE1104 style)
    # Pattern 3: mA, e_mA (Q2237 style)
    image_cols: Dict[str, Tuple[int, int]] = {}  # label -> (mag_col, err_col)
    
    for i, h in enumerate(header_lower):
        # Skip time and filter columns
        if h in ("mjd", "mhjd", "hjd", "jd", "band", "filt", "filter"):
            continue
        # Skip reference star columns
        if "beta" in h or "star" in h:
            continue
        
        # Pattern: e_X or e_XmagY -> error column
        if h.startswith("e_"):
            continue
        
        # Try to identify magnitude columns
        label: Optional[str] = None

        if len(h) == 1 and h.isalpha():
            # Single letter: A, B, C, D
            label = h.upper()
        elif h.startswith("m") and len(h) == 2 and h[1].isalpha():
            # mA, mB, mC, mD
            label = h[1].upper()
        elif "mag" in h:
            # HE0435-style columns: RmagA, VmagB, etc.
            # IMPORTANT: do NOT search for 'a' anywhere in the string, because "mag" contains "a".
            # Instead, take the trailing image letter when present.
            last = h[-1]
            if last.upper() in ("A", "B", "C", "D"):
                label = last.upper()
        
        if label and label not in image_cols:
            # Find corresponding error column.
            # Prefer exact match e_<magcol> (e.g., e_RmagA) when present.
            err_col: Optional[int] = None
            target_exact = f"e_{h}"
            for j, eh in enumerate(header_lower):
                if eh == target_exact:
                    err_col = j
                    break

            # Fallback: any e_* column that ends with the same image letter.
            if err_col is None:
                for j, eh in enumerate(header_lower):
                    if eh.startswith("e_") and eh[-1].upper() == label:
                        err_col = j
                        break

            if err_col is not None:
                image_cols[label] = (i, err_col)
    
    if len(image_cols) < 2:
        return None
    
    # Parse data
    data = {label: {"t": [], "mag": [], "err": []} for label in image_cols}
    
    from datetime import datetime
    
    for row in rows:
        if len(row) <= max(max(v) for v in image_cols.values()):
            continue
        
        try:
            if time_format == "calendar":
                # Parse calendar date (e.g., "1995-09-17")
                date_str = row[time_col]
                # Handle various date formats
                for fmt in ("%Y-%m-%d", "%Y/%m/%d", "%d-%m-%Y", "%d/%m/%Y"):
                    try:
                        dt = datetime.strptime(date_str, fmt)
                        # Convert to MJD (approximate: days since 1858-11-17)
                        t = (dt - datetime(1858, 11, 17)).days
                        break
                    except ValueError:
                        continue
                else:
                    continue
            else:
                t = float(row[time_col])
        except (ValueError, IndexError):
            continue
        
        for label, (mag_col, err_col) in image_cols.items():
            try:
                mag = float(row[mag_col])
                err = float(row[err_col])
                if np.isfinite(mag) and np.isfinite(err) and err > 0:
                    data[label]["t"].append(t)
                    data[label]["mag"].append(mag)
                    data[label]["err"].append(err)
            except (ValueError, IndexError):
                continue
    
    # Build LightCurve objects
    light_curves = {}
    for label in image_cols:
        # Minimum epochs threshold: adaptive based on data quality
        # Historical datasets (Vakulik et al.) have fewer epochs but longer baselines
        # Modern datasets (COSMOGRAIL) have more epochs with denser sampling
        # Threshold ensures sufficient data for reliable delay measurement
        # - 15 epochs: minimum for historical data (e.g., Q2237 Vakulik has 17 epochs)
        # - 20 epochs: standard for modern light curves with seasonal gaps
        min_epochs = 15 if time_format == "calendar" else 20
        if len(data[label]["t"]) >= min_epochs:
            light_curves[label] = LightCurve(
                label=label,
                t=np.array(data[label]["t"]),
                mag=np.array(data[label]["mag"]),
                magerr=np.array(data[label]["err"]),
            )
    
    if len(light_curves) < 2:
        return None
    
    return LensSystem(system_id=f"{system_id}_{band}", light_curves=light_curves, band=band)


def detrend_lightcurve(lc: LightCurve, window_days: float = 200.0) -> LightCurve:
    """
    Remove slow microlensing trend using Gaussian smoothing.
    
    Microlensing by stars in the lens galaxy causes slow (months-to-years)
    magnitude variations that are uncorrelated between images. This must be
    removed before measuring time delays, as it can bias delay estimates.
    
    Method:
    1. Interpolate light curve to uniform 1-day grid
    2. Mark large gaps (>30 days) as NaN to avoid interpolating across seasons
    3. Apply Gaussian smoothing with σ = window_days to extract trend
    4. Subtract trend from original data
    
    Args:
        lc: Input light curve
        window_days: Gaussian smoothing window (default 200 days captures
                     typical microlensing timescales while preserving
                     intrinsic quasar variability on shorter scales)
    
    Returns:
        Detrended light curve with microlensing removed
    
    Note:
        Uses NaN-aware processing to avoid fabricating data across seasonal
        observing gaps, which could introduce spurious correlations.
    """
    if len(lc.t) < 10:
        return lc
    
    # Interpolate to uniform grid for smoothing
    t_min, t_max = lc.t.min(), lc.t.max()
    dt = max(1.0, np.median(np.diff(lc.t)))  # At least 1 day sampling
    t_uniform = np.arange(t_min, t_max, dt)
    
    if len(t_uniform) < 10:
        return lc
    
    # Interpolate with NaN for out-of-bounds (no extrapolation)
    f_interp = interpolate.interp1d(lc.t, lc.mag, kind="linear", bounds_error=False, fill_value=np.nan)
    mag_uniform = f_interp(t_uniform)
    
    # Mark large gaps (>30 days) as NaN to avoid interpolating across seasons
    gap_threshold = 30.0  # days
    for i in range(len(lc.t) - 1):
        if lc.t[i+1] - lc.t[i] > gap_threshold:
            gap_start, gap_end = lc.t[i], lc.t[i+1]
            gap_mask = (t_uniform > gap_start) & (t_uniform < gap_end)
            mag_uniform[gap_mask] = np.nan
    
    if not np.any(np.isfinite(mag_uniform)):
        return lc

    # Gaussian smooth (window in samples)
    sigma_samples = window_days / dt
    trend = _gaussian_smooth_nan(mag_uniform, sigma=sigma_samples)
    
    # Interpolate trend back to original observation times only
    f_trend = interpolate.interp1d(t_uniform, trend, kind="linear", bounds_error=False, fill_value=np.nan)
    trend_orig = f_trend(lc.t)
    
    # Detrend only where we have valid trend estimates
    mag_detrended = lc.mag - trend_orig
    
    # If trend interpolation failed for some points, keep original
    bad_trend = np.isnan(trend_orig)
    mag_detrended[bad_trend] = lc.mag[bad_trend] - np.nanmean(lc.mag)
    
    return LightCurve(
        label=lc.label,
        t=lc.t.copy(),
        mag=mag_detrended,
        magerr=lc.magerr.copy(),
    )


def estimate_delay_correlation(
    lc1: LightCurve,
    lc2: LightCurve,
    lag_range: Tuple[float, float] = (-200, 200),
    lag_step: float = 1.0,
) -> Tuple[float, float, float]:
    """
    Estimate time delay between two light curves using cross-correlation.
    
    This is the core delay estimation algorithm. It finds the time shift that
    maximizes the Pearson correlation between the two light curves.
    
    Method:
    1. Interpolate both curves to a uniform 1-day grid
    2. Normalize each curve to zero mean and unit variance
    3. For each candidate lag in [lag_range], shift curve 2 and compute
       Pearson correlation with curve 1
    4. Find the lag with maximum correlation
    5. Apply parabolic refinement around peak for sub-day precision
    6. Estimate uncertainty from correlation peak width (FWHM)
    
    Convention:
        Positive delay means curve 2 (second image) arrives LATER than curve 1.
        i.e., if delay = +50 days, then mag2(t) ≈ mag1(t - 50)
    
    Args:
        lc1, lc2: Light curves to cross-correlate
        lag_range: (min_lag, max_lag) in days to search
        lag_step: Step size for lag grid (default 1 day)
    
    Returns:
        (delay_days, peak_correlation, delay_uncertainty)
        - delay_days: Best-fit time delay in days
        - peak_correlation: Pearson r at the best lag (quality indicator)
        - delay_uncertainty: Estimated 1σ uncertainty from FWHM
    
    Note:
        Uses interpolation-based shifting rather than array rolling to
        properly handle fractional lags and avoid discretization artifacts.
    """
    # Find common time range
    t_min = max(lc1.t.min(), lc2.t.min())
    t_max = min(lc1.t.max(), lc2.t.max())
    
    if t_max - t_min < 100:  # Need at least 100 days overlap
        return np.nan, np.nan, np.nan
    
    # Create uniform grid
    dt = 1.0  # 1 day sampling
    t_grid = np.arange(t_min, t_max, dt)
    
    def _mask_large_gaps(t_obs: np.ndarray, y_grid: np.ndarray) -> np.ndarray:
        if t_obs.size < 2:
            return y_grid
        gap_threshold = 30.0
        t_obs = np.asarray(t_obs, dtype=float)
        for i in range(t_obs.size - 1):
            if (t_obs[i + 1] - t_obs[i]) > gap_threshold:
                gap_start, gap_end = t_obs[i], t_obs[i + 1]
                m = (t_grid > gap_start) & (t_grid < gap_end)
                y_grid[m] = np.nan
        return y_grid

    # Interpolate both curves onto uniform grid
    # Use NaN for out-of-bounds instead of extrapolation
    f1 = interpolate.interp1d(lc1.t, lc1.mag, kind="linear", bounds_error=False, fill_value=np.nan)
    f2 = interpolate.interp1d(lc2.t, lc2.mag, kind="linear", bounds_error=False, fill_value=np.nan)

    y1 = _mask_large_gaps(lc1.t, f1(t_grid))
    y2 = _mask_large_gaps(lc2.t, f2(t_grid))
    
    # Normalize (using only valid data)
    y1 = (y1 - np.nanmean(y1)) / np.nanstd(y1)
    y2 = (y2 - np.nanmean(y2)) / np.nanstd(y2)
    
    def _interp_no_bridge(t: np.ndarray, y: np.ndarray, x: np.ndarray, gap_threshold: float = 30.0) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        y = np.asarray(y, dtype=float)
        x = np.asarray(x, dtype=float)

        if t.size < 2:
            return np.full_like(x, np.nan, dtype=float)

        out = np.full_like(x, np.nan, dtype=float)
        j = np.searchsorted(t, x, side="left")

        exact_mask = j < t.size
        if np.any(exact_mask):
            diff = np.full_like(x, np.inf, dtype=float)
            jj = j[exact_mask]
            diff[exact_mask] = np.abs(t[jj] - x[exact_mask])
            exact = exact_mask & (diff < 1e-9)
            if np.any(exact):
                out[exact] = y[j[exact]]

        m = (j > 0) & (j < t.size)
        if not np.any(m):
            return out

        j0 = j[m] - 1
        j1 = j[m]
        dt = t[j1] - t[j0]
        ok = dt <= float(gap_threshold)
        if not np.any(ok):
            return out

        x_ok = x[m][ok]
        j0_ok = j0[ok]
        j1_ok = j1[ok]
        dt_ok = dt[ok]

        w = (x_ok - t[j0_ok]) / dt_ok
        out_idx = np.flatnonzero(m)[ok]
        out[out_idx] = y[j0_ok] + w * (y[j1_ok] - y[j0_ok])
        return out

    t2f = t_grid[np.isfinite(y2)]
    y2f = y2[np.isfinite(y2)]
    
    # Scan lags using interpolation-based shifting
    lags = np.arange(lag_range[0], lag_range[1] + lag_step, lag_step)
    correlations = []
    
    for lag in lags:
        # Shift y2 by evaluating at t_grid - lag
        # If lag > 0, we're asking: what was y2 doing `lag` days earlier?
        y2_shifted = _interp_no_bridge(t2f, y2f, t_grid - lag, gap_threshold=30.0)
        
        # Compute correlation only on mutually valid points
        valid = np.isfinite(y1) & np.isfinite(y2_shifted)
        if np.sum(valid) < 20:
            correlations.append(np.nan)
            continue
        
        # Check for constant input to avoid ConstantInputWarning
        y1_valid = y1[valid]
        y2_valid = y2_shifted[valid]
        if np.std(y1_valid) < 1e-10 or np.std(y2_valid) < 1e-10:
            correlations.append(np.nan)
            continue
        
        r, _ = stats.pearsonr(y1_valid, y2_valid)
        correlations.append(r)
    
    correlations = np.array(correlations)
    
    # Find peak
    valid_corr = np.isfinite(correlations)
    if not np.any(valid_corr):
        return np.nan, np.nan, np.nan
    
    best_idx = np.nanargmax(correlations)
    best_lag = lags[best_idx]
    best_corr = correlations[best_idx]
    
    # Parabolic refinement around peak for sub-step precision
    if 1 <= best_idx < len(lags) - 1:
        y0, y1_val, y2_val = correlations[best_idx-1:best_idx+2]
        if np.all(np.isfinite([y0, y1_val, y2_val])):
            denom = 2 * (y0 - 2*y1_val + y2_val)
            if abs(denom) > 1e-10:
                delta = -(y2_val - y0) / denom
                # Clamp delta to [-1, 1] to stay within neighboring points
                delta = np.clip(delta, -1, 1)
                best_lag = lags[best_idx] + delta * lag_step
    
    # Estimate uncertainty from correlation peak width (FWHM method)
    half_max = best_corr / 2
    above_half = correlations > half_max
    if np.any(above_half):
        fwhm = np.sum(above_half) * lag_step
        delay_err = fwhm / 2.35482  # FWHM to sigma (precise: 2*sqrt(2*ln(2)))
    else:
        delay_err = 10.0  # Default fallback
    
    return float(best_lag), float(best_corr), float(delay_err)


def estimate_delay_dcf(
    lc1: LightCurve,
    lc2: LightCurve,
    lag_range: Tuple[float, float] = (-200, 200),
    lag_step: float = 1.0,
    min_pairs: int = 50,
) -> Tuple[float, float, float]:
    t1, x1 = lc1.t, lc1.mag
    t2, x2 = lc2.t, lc2.mag

    v1 = np.isfinite(t1) & np.isfinite(x1)
    v2 = np.isfinite(t2) & np.isfinite(x2)
    t1, x1 = t1[v1], x1[v1]
    t2, x2 = t2[v2], x2[v2]

    if t1.size < 20 or t2.size < 20:
        return np.nan, np.nan, np.nan

    x1 = (x1 - np.mean(x1)) / (np.std(x1) + 1e-12)
    x2 = (x2 - np.mean(x2)) / (np.std(x2) + 1e-12)

    dt = t2[:, None] - t1[None, :]
    udcf = x2[:, None] * x1[None, :]

    lags = np.arange(lag_range[0], lag_range[1] + lag_step, lag_step)
    half = 0.5 * lag_step
    dcf = np.full_like(lags, np.nan, dtype=float)

    for k, lag in enumerate(lags):
        m = (dt >= (lag - half)) & (dt < (lag + half))
        n = int(np.sum(m))
        if n < min_pairs:
            continue
        vals = udcf[m]
        dcf[k] = float(np.mean(vals))

    if not np.any(np.isfinite(dcf)):
        return np.nan, np.nan, np.nan

    best_idx = int(np.nanargmax(dcf))
    best_lag = float(lags[best_idx])
    best_corr = float(dcf[best_idx])

    half_max = best_corr / 2
    above_half = dcf > half_max
    if np.any(above_half):
        fwhm = float(np.sum(above_half) * lag_step)
        delay_err = float(max(1.0, fwhm / 2.35482))  # Precise FWHM to sigma (2*sqrt(2*ln(2)))
    else:
        delay_err = 10.0

    return best_lag, best_corr, delay_err


def estimate_delay_iccf(
    lc1: LightCurve,
    lc2: LightCurve,
    lag_range: Tuple[float, float] = (-200, 200),
    lag_step: float = 1.0,
    min_overlap: int = 20,
    centroid_frac: float = 0.8,
    gap_threshold: float = 30.0,
) -> Tuple[float, float, float]:
    def _interp_no_bridge(t: np.ndarray, y: np.ndarray, x: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        y = np.asarray(y, dtype=float)
        x = np.asarray(x, dtype=float)

        if t.size < 2:
            return np.full_like(x, np.nan, dtype=float)

        out = np.full_like(x, np.nan, dtype=float)

        j = np.searchsorted(t, x, side="left")

        exact_mask = j < t.size
        if np.any(exact_mask):
            diff = np.full_like(x, np.inf, dtype=float)
            jj = j[exact_mask]
            diff[exact_mask] = np.abs(t[jj] - x[exact_mask])
            exact = exact_mask & (diff < 1e-9)
            if np.any(exact):
                out[exact] = y[j[exact]]

        m = (j > 0) & (j < t.size)
        if not np.any(m):
            return out

        j0 = j[m] - 1
        j1 = j[m]

        dt = t[j1] - t[j0]
        ok = dt <= float(gap_threshold)
        if not np.any(ok):
            return out

        x_ok = x[m][ok]
        j0_ok = j0[ok]
        j1_ok = j1[ok]
        dt_ok = dt[ok]

        w = (x_ok - t[j0_ok]) / dt_ok
        out_idx = np.flatnonzero(m)[ok]
        out[out_idx] = y[j0_ok] + w * (y[j1_ok] - y[j0_ok])
        return out

    t1, y1 = np.asarray(lc1.t, dtype=float), np.asarray(lc1.mag, dtype=float)
    t2, y2 = np.asarray(lc2.t, dtype=float), np.asarray(lc2.mag, dtype=float)

    v1 = np.isfinite(t1) & np.isfinite(y1)
    v2 = np.isfinite(t2) & np.isfinite(y2)
    t1, y1 = t1[v1], y1[v1]
    t2, y2 = t2[v2], y2[v2]

    if t1.size < min_overlap or t2.size < min_overlap:
        return np.nan, np.nan, np.nan

    o1 = np.argsort(t1)
    o2 = np.argsort(t2)
    t1, y1 = t1[o1], y1[o1]
    t2, y2 = t2[o2], y2[o2]

    y1 = (y1 - np.mean(y1)) / (np.std(y1) + 1e-12)
    y2 = (y2 - np.mean(y2)) / (np.std(y2) + 1e-12)

    lags = np.arange(lag_range[0], lag_range[1] + lag_step, lag_step)
    corrs = np.full_like(lags, np.nan, dtype=float)

    for k, lag in enumerate(lags):
        y2_on_1 = _interp_no_bridge(t2, y2, t1 - lag)
        m12 = np.isfinite(y1) & np.isfinite(y2_on_1)

        y1_on_2 = _interp_no_bridge(t1, y1, t2 + lag)
        m21 = np.isfinite(y2) & np.isfinite(y1_on_2)

        rvals = []
        if int(np.sum(m12)) >= min_overlap:
            # Check for variance before correlation
            if np.std(y1[m12]) > 1e-10 and np.std(y2_on_1[m12]) > 1e-10:
                r, _ = stats.pearsonr(y1[m12], y2_on_1[m12])
                rvals.append(r)
        if int(np.sum(m21)) >= min_overlap:
            # Check for variance before correlation
            if np.std(y2[m21]) > 1e-10 and np.std(y1_on_2[m21]) > 1e-10:
                r, _ = stats.pearsonr(y2[m21], y1_on_2[m21])
                rvals.append(r)

        if rvals:
            corrs[k] = float(np.mean(rvals))

    if not np.any(np.isfinite(corrs)):
        return np.nan, np.nan, np.nan

    best_idx = int(np.nanargmax(corrs))
    peak_lag = float(lags[best_idx])
    peak_corr = float(corrs[best_idx])

    if 1 <= best_idx < len(lags) - 1:
        y0, y1v, y2v = corrs[best_idx - 1 : best_idx + 2]
        if np.all(np.isfinite([y0, y1v, y2v])):
            denom = 2 * (y0 - 2 * y1v + y2v)
            if abs(denom) > 1e-10:
                delta = -(y2v - y0) / denom
                delta = np.clip(delta, -1, 1)
                peak_lag = float(lags[best_idx] + delta * lag_step)

    frac = float(np.clip(centroid_frac, 0.0, 1.0))
    if frac > 0 and np.isfinite(peak_corr) and peak_corr > 0:
        sel = np.isfinite(corrs) & (corrs >= frac * peak_corr)
        if np.any(sel):
            w = corrs[sel].copy()
            w = np.maximum(w, 0.0)
            if np.sum(w) > 0:
                best_lag = float(np.sum(lags[sel] * w) / np.sum(w))
            else:
                best_lag = peak_lag
        else:
            best_lag = peak_lag
    else:
        best_lag = peak_lag

    best_corr = peak_corr
    half_max = best_corr / 2
    above_half = corrs > half_max
    if np.any(above_half):
        fwhm = float(np.sum(above_half) * lag_step)
        delay_err = float(max(1.0, fwhm / 2.35482))  # Precise FWHM to sigma (2*sqrt(2*ln(2)))
    else:
        delay_err = 10.0

    return best_lag, best_corr, delay_err


def bandpass_filter(
    lc: LightCurve,
    tau_center: float,
    tau_width_factor: float = 0.5,
    fill_nan_with_zero: bool = True,
) -> LightCurve:
    """
    Apply difference-of-Gaussians bandpass filter to isolate variability at timescale τ.
    
    This is the key operation for multiscale analysis. By isolating variability
    at specific timescales, we can measure how the time delay depends on the
    frequency content of the signal.
    
    Method:
    1. Interpolate to uniform 1-day grid
    2. Mark large gaps as NaN (prevents interpolation across seasons)
    3. Apply two Gaussian smoothings with different widths:
       - σ_low = τ × (1 - width_factor)  [passes higher frequencies]
       - σ_high = τ × (1 + width_factor) [passes lower frequencies]
    4. Subtract: filtered = smooth_low - smooth_high
       This isolates variability in the band [σ_low, σ_high]
    
    Args:
        lc: Input light curve
        tau_center: Central timescale in days (e.g., 20, 40, 80, 160)
        tau_width_factor: Fractional width of bandpass (default 0.5 = half-octave,
                          meaning the band spans from 0.5τ to 1.5τ)
    
    Returns:
        Bandpass-filtered light curve containing only variability near timescale τ
    
    Physical interpretation:
        - τ = 20 days: Isolates rapid flickering (days-to-weeks)
        - τ = 80 days: Isolates medium-term variations (months)
        - τ = 160 days: Isolates slow variations (half-year)
    
    Note:
        Quasars typically have weak variability at very short timescales (τ < 10d).
        If the bandpass filter preserves <2% of the original variance, the
        filtered signal is noise-dominated and should not be used for delay
        estimation.
    """
    if len(lc.t) < 20:
        return lc
    
    # Interpolate to uniform grid
    t_min, t_max = lc.t.min(), lc.t.max()
    dt = 1.0  # 1 day
    t_uniform = np.arange(t_min, t_max, dt)
    
    # Use NaN for out-of-bounds instead of extrapolation
    f_interp = interpolate.interp1d(lc.t, lc.mag, kind="linear", bounds_error=False, fill_value=np.nan)
    mag_uniform = f_interp(t_uniform)
    
    # Mark large gaps (>30 days) as NaN to avoid interpolating across seasons
    gap_threshold = 30.0  # days
    for i in range(len(lc.t) - 1):
        if lc.t[i+1] - lc.t[i] > gap_threshold:
            gap_start, gap_end = lc.t[i], lc.t[i+1]
            gap_mask = (t_uniform > gap_start) & (t_uniform < gap_end)
            mag_uniform[gap_mask] = np.nan
    
    if not np.any(np.isfinite(mag_uniform)):
        return lc
    
    # Difference of Gaussians
    sigma_low = tau_center * (1 - tau_width_factor) / dt
    sigma_high = tau_center * (1 + tau_width_factor) / dt
    
    if sigma_low < 1:
        sigma_low = 1
    if sigma_high < sigma_low + 1:
        sigma_high = sigma_low + 1
    
    smooth_low = _gaussian_smooth_nan(mag_uniform, sigma=sigma_low)
    smooth_high = _gaussian_smooth_nan(mag_uniform, sigma=sigma_high)
    
    filtered = smooth_low - smooth_high

    bad = ~np.isfinite(smooth_low) | ~np.isfinite(smooth_high)
    filtered[bad] = np.nan
    
    # Interpolate back to original times (no extrapolation)
    f_back = interpolate.interp1d(t_uniform, filtered, kind="linear", bounds_error=False, fill_value=np.nan)
    mag_filtered = f_back(lc.t)
    
    if fill_nan_with_zero:
        mag_filtered = np.nan_to_num(mag_filtered, nan=0.0)
    
    return LightCurve(
        label=lc.label,
        t=lc.t.copy(),
        mag=mag_filtered,
        magerr=lc.magerr.copy(),
    )


def compute_multiscale_delays(
    lc1: LightCurve,
    lc2: LightCurve,
    tau_values: List[float],
    lag_range: Tuple[float, float] = (-200, 200),
    lag_step: float = 1.0,
    min_correlation: float = 0.3,
    min_variance_fraction: float = 0.01,
    broadband_delay: Optional[float] = None,
    mode_lock_window: Optional[float] = None,
    estimator: str = "interp",
) -> Dict[float, Tuple[float, float, float]]:
    """
    Compute time delay at multiple variability timescales.
    
    This is the core of the temporal shear measurement. For each timescale τ,
    we bandpass filter both light curves and measure the delay between them.
    If the delay varies systematically with τ, this indicates temporal shear.
    
    KEY METHODOLOGICAL FEATURES:
    
    1. Variance filtering:
       Before measuring delay at each τ, we check if the bandpass filter
       preserves sufficient signal. If variance_filtered / variance_original < 2%,
       the filtered signal is noise-dominated and we skip that τ.
       
       This is critical because quasars have weak variability at short timescales
       (τ < 10 days). Attempting to measure delays from noise leads to spurious
       results that were previously misinterpreted as "mode-jumping".
    
    2. Mode-locking:
       If broadband_delay is provided, we constrain the delay search to
       [broadband_delay - mode_lock_window, broadband_delay + mode_lock_window].
       
       This prevents spurious alias peaks (e.g., from seasonal gaps or edge
       effects) from dominating the delay estimate. The true delay should be
       near the broadband value; large deviations indicate artifacts.
    
    Args:
        lc1, lc2: Light curves to correlate
        tau_values: List of timescales to analyze (e.g., [20, 40, 80, 160] days)
        lag_range: Global search range for delays (used if no mode-locking)
        min_correlation: Minimum Pearson r to accept (default 0.3)
        min_variance_fraction: Minimum variance preservation (default 0.01 = 1%)
        broadband_delay: Reference delay for mode-locking (from unfiltered data)
        mode_lock_window: Half-width of constrained search (default 50 days)
    
    Returns:
        Dictionary mapping τ -> (delay, correlation, uncertainty)
        NaN values indicate the measurement was rejected (low signal or correlation)
    """
    results = {}
    
    # Compute original variance for signal quality check
    var_orig_1 = np.nanvar(lc1.mag)
    var_orig_2 = np.nanvar(lc2.mag)
    
    for tau in tau_values:
        # Bandpass filter both curves
        if estimator in ("dcf", "iccf", "interp"):
            lc1_bp = bandpass_filter(lc1, tau, fill_nan_with_zero=False)
            lc2_bp = bandpass_filter(lc2, tau, fill_nan_with_zero=False)
        else:
            lc1_bp = bandpass_filter(lc1, tau)
            lc2_bp = bandpass_filter(lc2, tau)
        
        # Check if bandpass filter preserves sufficient signal
        var_bp_1 = np.nanvar(lc1_bp.mag)
        var_bp_2 = np.nanvar(lc2_bp.mag)
        
        frac_1 = var_bp_1 / var_orig_1 if var_orig_1 > 0 else 0
        frac_2 = var_bp_2 / var_orig_2 if var_orig_2 > 0 else 0
        
        if frac_1 < min_variance_fraction or frac_2 < min_variance_fraction:
            # Insufficient signal at this timescale
            results[tau] = (np.nan, np.nan, np.nan)
            continue
        
        # Determine search range (optionally constrained by mode-locking)
        if broadband_delay is not None and mode_lock_window is not None:
            search_range = (
                broadband_delay - mode_lock_window,
                broadband_delay + mode_lock_window
            )
        else:
            search_range = lag_range
        
        # Estimate delay
        if estimator == "dcf":
            delay, corr, err = estimate_delay_dcf(lc1_bp, lc2_bp, search_range, lag_step=lag_step)
        elif estimator == "iccf":
            delay, corr, err = estimate_delay_iccf(lc1_bp, lc2_bp, search_range, lag_step=lag_step)
        else:
            delay, corr, err = estimate_delay_correlation(lc1_bp, lc2_bp, search_range, lag_step=lag_step)
        
        # Reject low-correlation estimates as unreliable
        if not np.isfinite(corr) or corr < min_correlation:
            delay = np.nan
            err = np.nan
        
        results[tau] = (delay, corr, err)
    
    return results


def fit_gamma(
    tau_values: List[float],
    delays: List[float],
    delay_errors: List[float],
    correlations: Optional[List[float]] = None,
    min_valid_points: int = 3,
) -> Tuple[float, float, float, float]:
    """
    Fit Γ = d(Δt)/d(log τ) — the temporal shear slope.
    
    This fits a linear relation between the measured delays and log₁₀(τ):
    
        Δt(τ) = Γ · log10(τ) + intercept
    
    The slope Γ is the "temporal shear" in units of days per decade of timescale.
    For example, Γ = +30 days/decade means the delay increases by 30 days when
    going from τ = 10 days to τ = 100 days.
    
    STATISTICAL METHOD:
    Uses weighted least squares (WLS) for the linear model

        Δt(τ) = Γ · log10(τ) + intercept,

    where the weights are derived from delay uncertainties (inverse variance).
    This is the proper statistical approach: weight by precision, not correlation.
    High correlation does NOT imply high precision (correlation can be high at 
    wrong lags due to aliasing). Only delay uncertainty reflects true precision.

    If correlations are provided, they are used only for quality cuts (rejecting
    low-correlation points) but NOT for weighting. This prevents correlation-based
    bias while still filtering out poor measurements.
    
    Args:
        tau_values: Timescales analyzed (days)
        delays: Measured delays at each τ (days)
        delay_errors: Delay uncertainties (used for inverse-variance weighting)
        correlations: Correlation coefficients at each τ (optional; used for quality cuts only)
        min_valid_points: Minimum number of valid delays required (default 3)
    
    Returns:
        (gamma, gamma_err, intercept, r_squared)
        - gamma: Temporal shear slope (days per decade)
        - gamma_err: Standard error on gamma
        - intercept: Fitted intercept (delay at τ = 1 day, extrapolated)
        - r_squared: Coefficient of determination (goodness of fit)
    
    Interpretation:
        - |Γ| / γ_err > 2: Marginally significant temporal shear
        - |Γ| / γ_err > 3: Significant temporal shear detection
        - GR predicts Γ = 0 (no scale dependence)
    """
    valid = [i for i in range(len(tau_values)) if np.isfinite(delays[i]) and np.isfinite(delay_errors[i]) and delay_errors[i] > 0]

    if len(valid) < min_valid_points:
        return np.nan, np.nan, np.nan, np.nan

    log_tau = np.log10([tau_values[i] for i in valid])
    dt = np.array([delays[i] for i in valid])
    dt_err = np.array([delay_errors[i] for i in valid])
    
    # Quality cut: reject points with very low correlation if correlations provided
    if correlations is not None:
        corr_arr = np.array([correlations[i] for i in valid])
        # Require minimum correlation for reliability (not for weighting)
        good_quality = corr_arr >= 0.3  # Minimum quality threshold
        if np.sum(good_quality) >= min_valid_points:
            log_tau = log_tau[good_quality]
            dt = dt[good_quality]
            dt_err = dt_err[good_quality]
        # If too few points pass quality cut, use all valid points with warning
    
    # PROPER WEIGHTING: Inverse variance from delay uncertainties
    # This is the correct statistical approach - weight by precision
    w = 1.0 / (dt_err ** 2)
    w = w / np.sum(w)  # Normalize weights

    X = np.column_stack([log_tau, np.ones_like(log_tau)])
    W = np.diag(w)

    try:
        XtWX = X.T @ W @ X
        XtWy = X.T @ W @ dt
        beta = np.linalg.solve(XtWX, XtWy)
        gamma = float(beta[0])
        intercept = float(beta[1])

        dt_pred = X @ beta
        resid = dt - dt_pred

        dof = max(1, len(dt) - 2)
        # Weighted residual variance
        s2 = float((resid.T @ W @ resid) / dof)
        cov = np.linalg.inv(XtWX) * s2
        gamma_err = float(np.sqrt(cov[0, 0]))

        # R-squared (weighted)
        dt_mean = float(np.average(dt, weights=w))
        ss_res = float(np.sum(w * (dt - dt_pred) ** 2))
        ss_tot = float(np.sum(w * (dt - dt_mean) ** 2))
        r_squared = float(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0

        return gamma, gamma_err, intercept, r_squared
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan, np.nan


def analyze_system(
    system: LensSystem,
    detrend_window: float = 200.0,
    tau_values: Optional[List[float]] = None,
    lag_range: Tuple[float, float] = (-200, 200),
    lag_step: float = 1.0,
    mode_lock_window: float = 50.0,
    min_variance_fraction: float = 0.02,
    min_correlation: float = 0.5,
    estimator: str = "interp",
    broadband_estimator: Optional[str] = None,
    no_bandpass: bool = False,
) -> Dict:
    """
    Full temporal shear analysis for one lens system.
    
    Args:
        system: LensSystem to analyze
        detrend_window: Window for detrending (days)
        tau_values: Timescales to analyze
        lag_range: Global search range for broadband delay
        mode_lock_window: Half-width of search window around broadband delay
            for multiscale analysis (prevents spurious alias peaks)
        min_variance_fraction: Minimum fraction of variance that must be
            preserved by bandpass filter (filters out noise-dominated scales)
    """
    if tau_values is None:
        # Default: skip very small scales where quasar variability is weak
        tau_values = [40, 80, 160]
    
    if broadband_estimator is None:
        broadband_estimator = estimator
        if estimator == "iccf":
            broadband_estimator = "interp"

    results = {
        "system_id": system.system_id,
        "band": system.band,
        "n_images": system.n_images,
        "image_labels": system.image_labels,
        "tau_values": tau_values,
        "detrend_window_days": detrend_window,
        "mode_lock_window_days": mode_lock_window,
        "min_variance_fraction": min_variance_fraction,
        "pairs": {},
    }
    
    # Detrend all light curves
    detrended = {
        label: detrend_lightcurve(lc, detrend_window)
        for label, lc in system.light_curves.items()
    }
    
    # Analyze each image pair
    for l1, l2 in system.get_image_pairs():
        pair_key = f"{l1}-{l2}"
        print_status(f"  Analyzing pair {pair_key}")
        
        lc1 = detrended[l1]
        lc2 = detrended[l2]
        
        # Broadband delay (no bandpass, full search range)
        if broadband_estimator == "dcf":
            delay_bb, corr_bb, err_bb = estimate_delay_dcf(lc1, lc2, lag_range, lag_step=lag_step)
        elif broadband_estimator == "iccf":
            delay_bb, corr_bb, err_bb = estimate_delay_iccf(lc1, lc2, lag_range, lag_step=lag_step)
        else:
            delay_bb, corr_bb, err_bb = estimate_delay_correlation(lc1, lc2, lag_range, lag_step=lag_step)
        
        # Multi-scale delays with mode-locking around broadband delay
        # This prevents spurious alias peaks from dominating
        if not no_bandpass:
            multiscale = compute_multiscale_delays(
                lc1, lc2, tau_values, lag_range,
                lag_step=lag_step,
                min_correlation=min_correlation,
                min_variance_fraction=min_variance_fraction,
                broadband_delay=delay_bb if np.isfinite(delay_bb) else None,
                mode_lock_window=mode_lock_window if np.isfinite(delay_bb) else None,
                estimator=estimator,
            )
        else:
            multiscale = compute_multiscale_delays(
                lc1, lc2, [1.0], lag_range,  # Pass [1.0] as the tau_values argument
                lag_step=lag_step,
                min_correlation=min_correlation,
                min_variance_fraction=min_variance_fraction,
                broadband_delay=delay_bb if np.isfinite(delay_bb) else None,
                mode_lock_window=mode_lock_window if np.isfinite(delay_bb) else None,
                estimator=estimator,
            )
            multiscale = {tau: multiscale[1.0] for tau in tau_values}
        
        # Extract for Gamma fit
        delays = [multiscale[tau][0] for tau in tau_values]
        corrs = [multiscale[tau][1] for tau in tau_values]
        errs = [multiscale[tau][2] for tau in tau_values]
        
        # Count valid delay estimates
        n_valid = sum(1 for d in delays if np.isfinite(d))
        
        # Mode-jump diagnostic: max adjacent-tau delay jump
        # Large jumps indicate alias hopping between modes
        valid_delays = [(tau_values[i], delays[i]) for i in range(len(delays)) if np.isfinite(delays[i])]
        if len(valid_delays) >= 2:
            adjacent_jumps = [abs(valid_delays[i+1][1] - valid_delays[i][1]) for i in range(len(valid_delays)-1)]
            max_delay_jump = max(adjacent_jumps)
        else:
            max_delay_jump = np.nan
        
        # Fit Gamma (weight by correlation strength)
        gamma, gamma_err, intercept, r_sq = fit_gamma(
            tau_values, delays, errs, correlations=corrs, min_valid_points=3
        )
        
        # Significance of Gamma
        if np.isfinite(gamma) and np.isfinite(gamma_err) and gamma_err > 0:
            gamma_sigma = abs(gamma) / gamma_err
            gamma_pvalue = 2 * (1 - stats.norm.cdf(gamma_sigma))
        else:
            gamma_sigma = np.nan
            gamma_pvalue = np.nan

        # Mode stability threshold: reject if max adjacent-tau jump exceeds 100 days
        # This prevents alias-driven Gamma from passing
        mode_stable = np.isfinite(max_delay_jump) and (max_delay_jump <= 100.0)
        
        practical_pass = (
            np.isfinite(gamma)
            and (abs(gamma) >= 5.0)
            and np.isfinite(r_sq)
            and (r_sq >= 0.3)
            and (n_valid >= 3)
            and mode_stable
        )
        
        results["pairs"][pair_key] = {
            "broadband": {
                "delay_days": delay_bb,
                "correlation": corr_bb,
                "uncertainty_days": err_bb,
            },
            "multiscale": {
                str(tau): {
                    "delay_days": multiscale[tau][0],
                    "correlation": multiscale[tau][1],
                    "uncertainty_days": multiscale[tau][2],
                }
                for tau in tau_values
            },
            "gamma": {
                "value": gamma,
                "uncertainty": gamma_err,
                "sigma": gamma_sigma,
                "p_value": gamma_pvalue,
                "intercept": intercept,
                "r_squared": r_sq,
                "passes_practical_cuts": bool(practical_pass),
                "uncertainty_method": "bootstrap",
            },
            "n_valid_scales": n_valid,
            "max_delay_jump_days": float(max_delay_jump) if np.isfinite(max_delay_jump) else None,
            "mode_stable": bool(mode_stable),
        }
    
    return results


def _bootstrap_single_iteration(
    b: int,
    system_data: Dict,
    detrend_window: float,
    tau_values: List[float],
    lag_range: Tuple[float, float],
    lag_step: float,
    mode_lock_window: float,
    min_variance_fraction: float,
    min_correlation: float,
    estimator: str,
    broadband_estimator: Optional[str],
    bootstrap_mode: str,
    no_bandpass: bool,
) -> Dict[str, Optional[float]]:
    """
    Single bootstrap iteration - designed for parallel execution.
    
    Returns dict mapping pair_key -> gamma value (or None if failed)
    """
    rng = np.random.default_rng(b)  # Unique seed per iteration
    
    # Reconstruct light curves from serialized data
    perturbed_lcs = {}
    for label, lc_data in system_data.items():
        n = int(len(lc_data['t']))
        if n < 5:
            perturbed_lcs[label] = LightCurve(
                label=label,
                t=np.array(lc_data['t']),
                mag=np.array(lc_data['mag']),
                magerr=np.array(lc_data['magerr']),
            )
            continue
        
        if bootstrap_mode == "fr":
            t_bs = np.array(lc_data['t'])
            mag_bs = np.array(lc_data['mag'])
            err_bs = np.array(lc_data['magerr'])
        else:
            idx = rng.integers(0, n, size=n)
            uniq, counts = np.unique(idx, return_counts=True)
            t_bs = np.array(lc_data['t'])[uniq]
            mag_bs = np.array(lc_data['mag'])[uniq]
            err_bs = np.array(lc_data['magerr'])[uniq] / np.sqrt(counts.astype(float))
        
        order = np.argsort(t_bs)
        t_bs = t_bs[order]
        mag_bs = mag_bs[order]
        err_bs = err_bs[order]
        
        noise = rng.standard_normal(len(t_bs)) * err_bs
        perturbed_lcs[label] = LightCurve(
            label=label,
            t=t_bs.copy(),
            mag=(mag_bs + noise).copy(),
            magerr=err_bs.copy(),
        )
    
    # Analyze each pair
    results = {}
    pair_keys = list(perturbed_lcs.keys())
    
    for i, l1 in enumerate(pair_keys):
        for l2 in pair_keys[i+1:]:
            pair_key = f"{l1}-{l2}"
            lc1 = detrend_lightcurve(perturbed_lcs[l1], detrend_window)
            lc2 = detrend_lightcurve(perturbed_lcs[l2], detrend_window)
            
            if broadband_estimator == "dcf":
                delay_bb, corr_bb, err_bb = estimate_delay_dcf(lc1, lc2, lag_range, lag_step=lag_step)
            elif broadband_estimator == "iccf":
                delay_bb, corr_bb, err_bb = estimate_delay_iccf(lc1, lc2, lag_range, lag_step=lag_step)
            else:
                delay_bb, corr_bb, err_bb = estimate_delay_correlation(lc1, lc2, lag_range, lag_step=lag_step)
            
            multiscale = compute_multiscale_delays(
                lc1, lc2, tau_values, lag_range,
                lag_step=lag_step,
                min_correlation=min_correlation,
                min_variance_fraction=min_variance_fraction,
                broadband_delay=delay_bb if np.isfinite(delay_bb) else None,
                mode_lock_window=mode_lock_window if np.isfinite(delay_bb) else None,
                estimator=estimator,
            )
            delays = [multiscale[tau][0] for tau in tau_values]
            corrs = [multiscale[tau][1] for tau in tau_values]
            errs = [multiscale[tau][2] for tau in tau_values]
            
            gamma, _, _, _ = fit_gamma(tau_values, delays, errs, correlations=corrs, min_valid_points=3)
            results[pair_key] = float(gamma) if np.isfinite(gamma) else None
    
    return results


def _jackknife_single_iteration(
    omit_idx: int,
    season_intervals: List[Tuple[float, float]],
    system_data: Dict,
    detrend_window: float,
    tau_values: List[float],
    lag_range: Tuple[float, float],
    lag_step: float,
    mode_lock_window: float,
    min_variance_fraction: float,
    min_correlation: float,
    estimator: str,
    broadband_estimator: Optional[str],
) -> Dict[str, Optional[float]]:
    """
    Single jackknife iteration - designed for parallel execution.
    
    Returns dict mapping pair_key -> gamma value (or None if failed)
    """
    keep_intervals = [season_intervals[j] for j in range(len(season_intervals)) if j != omit_idx]
    
    # Build subset light curves
    subset_lcs = {}
    for label, lc_data in system_data.items():
        t_arr = np.array(lc_data['t'])
        m = _mask_times_in_intervals(t_arr, keep_intervals)
        subset_lcs[label] = LightCurve(
            label=label,
            t=t_arr[m].copy(),
            mag=np.array(lc_data['mag'])[m].copy(),
            magerr=np.array(lc_data['magerr'])[m].copy(),
        )
    
    detrended = {lbl: detrend_lightcurve(lc, detrend_window) for lbl, lc in subset_lcs.items()}
    
    # Analyze each pair
    results = {}
    labels = list(detrended.keys())
    
    for i, l1 in enumerate(labels):
        for l2 in labels[i+1:]:
            pair_key = f"{l1}-{l2}"
            lc1 = detrended[l1]
            lc2 = detrended[l2]
            
            if broadband_estimator == "dcf":
                delay_bb, corr_bb, err_bb = estimate_delay_dcf(lc1, lc2, lag_range, lag_step=lag_step)
            elif broadband_estimator == "iccf":
                delay_bb, corr_bb, err_bb = estimate_delay_iccf(lc1, lc2, lag_range, lag_step=lag_step)
            else:
                delay_bb, corr_bb, err_bb = estimate_delay_correlation(lc1, lc2, lag_range, lag_step=lag_step)
            
            multiscale = compute_multiscale_delays(
                lc1, lc2, tau_values, lag_range,
                lag_step=lag_step,
                min_correlation=min_correlation,
                min_variance_fraction=min_variance_fraction,
                broadband_delay=delay_bb if np.isfinite(delay_bb) else None,
                mode_lock_window=mode_lock_window if np.isfinite(delay_bb) else None,
                estimator=estimator,
            )
            delays = [multiscale[tau][0] for tau in tau_values]
            corrs = [multiscale[tau][1] for tau in tau_values]
            errs = [multiscale[tau][2] for tau in tau_values]
            
            gamma, _, _, _ = fit_gamma(tau_values, delays, errs, correlations=corrs, min_valid_points=3)
            results[pair_key] = float(gamma) if np.isfinite(gamma) else None
    
    return results


def run_bootstrap_uncertainty(
    system: LensSystem,
    n_bootstrap: int = 100,
    detrend_window: float = 200.0,
    tau_values: Optional[List[float]] = None,
    lag_range: Tuple[float, float] = (-200, 200),
    lag_step: float = 1.0,
    mode_lock_window: float = 50.0,
    min_variance_fraction: float = 0.02,
    min_correlation: float = 0.5,
    estimator: str = "interp",
    broadband_estimator: Optional[str] = None,
    bootstrap_mode: str = "frrss",
    no_bandpass: bool = False,
    n_jobs: int = -1,
) -> Dict[str, Dict[str, float]]:
    """
    Bootstrap uncertainty estimation for Gamma.
    
    M4 Pro Optimized: Uses parallel processing via joblib.
    """
    if tau_values is None:
        tau_values = [5, 10, 20, 40, 80, 160]

    sigma_floor = 2.0

    all_times = np.concatenate([lc.t for lc in system.light_curves.values()])
    season_intervals = _compute_season_intervals(all_times, gap_days=30.0)
    n_seasons = len(season_intervals)

    pair_keys = [f"{l1}-{l2}" for l1, l2 in system.get_image_pairs()]
    gamma_samples_phot = {k: [] for k in pair_keys}
    gamma_samples_jack = {k: [] for k in pair_keys}

    if broadband_estimator is None:
        broadband_estimator = estimator
        if estimator == "iccf":
            broadband_estimator = "interp"

    # Serialize light curves for parallel processing
    system_data = {}
    for label, lc in system.light_curves.items():
        system_data[label] = {
            't': lc.t.tolist(),
            'mag': lc.mag.tolist(),
            'magerr': lc.magerr.tolist(),
        }

    # Determine number of jobs
    if n_jobs == -1:
        n_jobs = cpu_count()
    
    # Parallel jackknife
    if n_seasons >= 2:
        print(f"  Running {n_seasons} jackknife iterations using {min(n_jobs, n_seasons)} cores...")
        jack_results = Parallel(n_jobs=min(n_jobs, n_seasons), backend='loky', verbose=0)(
            delayed(_jackknife_single_iteration)(
                omit_idx, season_intervals, system_data,
                detrend_window, tau_values, lag_range, lag_step,
                mode_lock_window, min_variance_fraction, min_correlation,
                estimator, broadband_estimator
            )
            for omit_idx in range(n_seasons)
        )
        
        for result in jack_results:
            for pair_key, gamma in result.items():
                if gamma is not None:
                    gamma_samples_jack[pair_key].append(gamma)

    # Parallel bootstrap
    print(f"  Running {n_bootstrap} bootstrap iterations using {min(n_jobs, n_bootstrap)} cores...")
    boot_results = Parallel(n_jobs=min(n_jobs, n_bootstrap), backend='loky', verbose=0)(
        delayed(_bootstrap_single_iteration)(
            b, system_data, detrend_window, tau_values, lag_range, lag_step,
            mode_lock_window, min_variance_fraction, min_correlation,
            estimator, broadband_estimator, bootstrap_mode, no_bandpass
        )
        for b in range(n_bootstrap)
    )
    
    for result in boot_results:
        for pair_key, gamma in result.items():
            if gamma is not None:
                gamma_samples_phot[pair_key].append(gamma)

    # Compute bootstrap statistics
    bootstrap_results = {}
    for pair_key in pair_keys:
        phot = gamma_samples_phot.get(pair_key, [])
        jack = gamma_samples_jack.get(pair_key, [])

        phot_std = float(np.std(phot)) if len(phot) > 10 else np.nan

        if len(jack) >= 2:
            jack_mean = float(np.mean(jack))
            n_jack = len(jack)
            jack_var = float((n_jack - 1) / n_jack * np.sum((np.array(jack) - jack_mean) ** 2))
            jack_std = float(np.sqrt(jack_var))
        else:
            jack_std = np.nan

        if np.isfinite(phot_std) and np.isfinite(jack_std):
            total_std = float(np.sqrt(phot_std**2 + jack_std**2))
        elif np.isfinite(jack_std):
            total_std = float(jack_std)
        elif np.isfinite(phot_std):
            total_std = float(phot_std)
        else:
            total_std = np.nan

        if np.isfinite(total_std):
            total_std = float(max(total_std, sigma_floor))

        gamma_mean = float(np.mean(phot)) if len(phot) > 0 else (float(np.mean(jack)) if len(jack) > 0 else np.nan)

        if np.isfinite(total_std):
            bootstrap_results[pair_key] = {
                "gamma_mean": gamma_mean,
                "gamma_std": total_std,
                "gamma_std_phot": phot_std,
                "gamma_std_jackknife": jack_std,
                "n_valid": len(phot),
                "n_valid_jackknife": len(jack),
                "n_seasons": n_seasons,
                "sigma_floor": sigma_floor,
            }
        else:
            bootstrap_results[pair_key] = {
                "gamma_mean": np.nan,
                "gamma_std": np.nan,
                "gamma_std_phot": phot_std,
                "gamma_std_jackknife": jack_std,
                "n_valid": len(phot),
                "n_valid_jackknife": len(jack),
                "n_seasons": n_seasons,
                "sigma_floor": sigma_floor,
            }
    
    return bootstrap_results


def historical_interp_estimator(lc1: LightCurve, lc2: LightCurve, 
                               lag_range: Tuple[float, float] = (-200, 200),
                               lag_step: float = 1.0) -> Tuple[float, float, float]:
    """
    Original INTERP estimator from 2024 analysis with weighted interpolation and centroid lag selection.
    """
    t1, y1 = np.asarray(lc1.t), np.asarray(lc1.mag)
    t2, y2 = np.asarray(lc2.t), np.asarray(lc2.mag)
    
    # Standardize
    y1 = (y1 - np.nanmean(y1)) / np.nanstd(y1)
    y2 = (y2 - np.nanmean(y2)) / np.nanstd(y2)
    
    lags = np.arange(lag_range[0], lag_range[1] + lag_step, lag_step)
    corrs = np.full(lags.shape, np.nan)
    
    for i, lag in enumerate(lags):
        # Weighted interpolation based on time proximity
        interp_y2 = np.zeros_like(t1)
        for j, t in enumerate(t1):
            time_diff = np.abs(t2 - (t - lag))
            weights = np.exp(-time_diff / 10.0)  # Historical decay parameter
            weights[np.isnan(y2)] = 0
            if np.sum(weights) > 0:
                interp_y2[j] = np.sum(y2 * weights) / np.sum(weights)
            else:
                interp_y2[j] = np.nan
        
        # Compute correlation
        mask = np.isfinite(y1) & np.isfinite(interp_y2)
        if np.sum(mask) > 20:  # Min 20 points
            corrs[i] = np.corrcoef(y1[mask], interp_y2[mask])[0, 1]
    
    # Centroid lag selection above 80% max correlation
    max_corr = np.nanmax(corrs)
    if not np.isfinite(max_corr):
        return np.nan, np.nan, np.nan
    
    centroid_threshold = 0.8 * max_corr
    valid_lags = lags[corrs >= centroid_threshold]
    valid_corrs = corrs[corrs >= centroid_threshold]
    
    if len(valid_lags) == 0:
        best_lag = lags[np.nanargmax(corrs)]
    else:
        best_lag = np.sum(valid_lags * valid_corrs) / np.sum(valid_corrs)
    
    # Uncertainty based on correlation width
    half_max = max_corr / 2
    fwhm = np.sum(corrs >= half_max) * lag_step
    lag_err = fwhm / 2.35482 if fwhm > 0 else 10.0  # Precise FWHM to sigma (2*sqrt(2*ln(2)))
    
    return best_lag, max_corr, lag_err


def main():
    """Run COSMOGRAIL temporal shear analysis on strong lens systems.
    
    Parses COSMOGRAIL light curve data, analyzes multiple lens systems,
    measures time delays at multiple timescales, and detects temporal shear.
    Outputs results to JSON and prints summary statistics.
    """
    parser = argparse.ArgumentParser(description="COSMOGRAIL Temporal Shear Analysis")
    parser.add_argument(
        "--data-dir",
        type=str,
        default="data/cosmograil",
        help="Directory containing COSMOGRAIL .rdb files",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/outputs",
        help="Output directory for results",
    )
    parser.add_argument(
        "--figure-dir",
        type=str,
        default="results/figures",
        help="Output directory for figures",
    )
    parser.add_argument(
        "--detrend-window",
        type=float,
        default=200.0,
        help="Microlensing detrend window in days",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=100,
        help="Number of bootstrap iterations",
    )
    parser.add_argument(
        "--bootstrap-mode",
        type=str,
        default="frrss",
        choices=["fr", "frrss"],
        help="Bootstrap mode: fr (photometric only, fixed sampling) or frrss (photometric+sampling)",
    )
    parser.add_argument(
        "--tau-values",
        type=str,
        default="40,80,160",
        help="Comma-separated tau values in days",
    )
    parser.add_argument(
        "--mode-lock-window",
        type=float,
        default=50.0,
        help="Half-width of multiscale lag search window around broadband delay (days)",
    )
    parser.add_argument(
        "--min-variance-fraction",
        type=float,
        default=0.02,
        help="Minimum variance fraction preserved by bandpass filter",
    )
    parser.add_argument(
        "--min-correlation",
        type=float,
        default=0.5,
        help="Minimum correlation to accept multiscale delay estimate",
    )
    parser.add_argument(
        "--lag-step",
        type=float,
        default=1.0,
        help="Lag grid step (days)",
    )
    parser.add_argument(
        "--estimator",
        type=str,
        default="interp",
        choices=["interp", "iccf", "dcf", "historical_interp"],
        help="Delay estimator: uniform-grid CCF (interp), ICCF (iccf), or DCF (dcf)",
    )
    parser.add_argument(
        "--broadband-estimator",
        type=str,
        default="",
        choices=["", "interp", "iccf", "dcf"],
        help="Estimator for broadband delay used to set mode-lock center. Default: same as --estimator (except iccf defaults to interp)",
    )
    parser.add_argument(
        "--output-tag",
        type=str,
        default="",
        help="Tag to append to output filenames",
    )
    parser.add_argument(
        "--include-systems",
        type=str,
        default="",
        help="Comma-separated base system IDs to include (e.g., DESJ0408,PG1115,RXJ1131). Matches the filename prefix before '_'",
    )
    parser.add_argument(
        "--exclude-systems",
        type=str,
        default="",
        help="Comma-separated base system IDs to exclude",
    )
    parser.add_argument(
        "--include-files",
        type=str,
        default="",
        help="Comma-separated filename globs to include (applied to both RDB and CSV files), e.g. 'DESJ0408_*.rdb,PG1115_*.rdb,he0435_*.csv'",
    )
    parser.add_argument(
        "--include-pairs",
        type=str,
        default="",
        help="Comma-separated image pairs to include (e.g., 'A-B,B-C')",
    )
    parser.add_argument(
        "--gap-threshold",
        type=float,
        default=30.0,
        help="Maximum gap size (days) allowed for interpolation (default: 30)"
    )
    parser.add_argument(
        "--no-bandpass",
        action="store_true",
        help="Disable Gaussian bandpass filtering"
    )
    args = parser.parse_args()
    
    # Data acquisition check
    print_status("Checking COSMOGRAIL data availability...")
    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        print_status(f"Data directory not found: {data_dir}", "ERROR")
        print_status("Please ensure COSMOGRAIL light curve files are in data/cosmograil/", "ERROR")
        sys.exit(1)
    
    # Check for required RDB files
    required_systems = ["HE0435", "HS2209", "J1001", "J1206", "PG1115", "RXJ1131", "WFI2033", "DESJ0408"]
    available_rdb = list(data_dir.glob("*.rdb"))
    available_systems = [f.stem.split("_")[0].upper() for f in available_rdb]
    
    missing = [s for s in required_systems if s not in available_systems]
    if missing:
        print_status(f"⚠ Missing light curves for: {', '.join(missing)}", "WARNING")
        print_status(f"✓ Available: {len(available_rdb)} RDB files", "INFO")
    else:
        print_status(f"✓ All required COSMOGRAIL data present ({len(available_rdb)} RDB files)", "SUCCESS")
    
    if len(available_rdb) == 0:
        print_status("No COSMOGRAIL data files found. Cannot proceed.", "ERROR")
        sys.exit(1)
    
    output_dir = Path(args.output_dir)
    figure_dir = Path(args.figure_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    
    tag = f"_{args.output_tag}" if args.output_tag else ""
    
    def _base_system_from_path(p: Path) -> str:
        return p.stem.split("_")[0].upper()

    include_systems = {
        s.strip().upper() for s in args.include_systems.split(",") if s.strip()
    } if args.include_systems else set()
    exclude_systems = {
        s.strip().upper() for s in args.exclude_systems.split(",") if s.strip()
    } if args.exclude_systems else set()
    include_patterns = [
        s.strip() for s in args.include_files.split(",") if s.strip()
    ]
    include_pairs = {
        p.strip() for p in args.include_pairs.split(",") if p.strip()
    } if args.include_pairs else set()

    def _matches_include(p: Path) -> bool:
        if not include_patterns:
            return True
        return any(fnmatch.fnmatch(p.name, pat) for pat in include_patterns)

    # Find all .rdb files
    rdb_files = sorted([p for p in data_dir.glob("*.rdb") if _matches_include(p)])

    # Find multi-band CSV files (exclude glendama metadata files)
    csv_files = sorted([
        f for f in data_dir.glob("*.csv")
        if _matches_include(f)
        and not f.stem.startswith("glendama_")
        and any(x in f.stem.lower() for x in ("he0435", "he1104", "q2237"))
    ])

    if include_systems:
        rdb_files = [p for p in rdb_files if _base_system_from_path(p) in include_systems]
        csv_files = [p for p in csv_files if _base_system_from_path(p) in include_systems]

    if exclude_systems:
        rdb_files = [p for p in rdb_files if _base_system_from_path(p) not in exclude_systems]
        csv_files = [p for p in csv_files if _base_system_from_path(p) not in exclude_systems]
    
    if not rdb_files and not csv_files:
        print_status(f"No light curve files found in {data_dir}", "ERROR")
        sys.exit(1)
    
    print_status(f"Found {len(rdb_files)} RDB + {len(csv_files)} multi-band CSV light curve files")
    
    # Timescales to probe
    tau_values = [float(x) for x in args.tau_values.split(",") if x.strip()]
    
    all_results = {
        "analysis_date": datetime.now().isoformat(),
        "detrend_window_days": args.detrend_window,
        "tau_values": tau_values,
        "n_bootstrap": args.n_bootstrap,
        "systems": {},
    }
    
    gamma_summary = []
    
    for rdb_file in rdb_files:
        print_status(f"Processing {rdb_file.name}")
        
        try:
            system = parse_rdb_file(rdb_file)
        except Exception as e:
            print_status(f"  Failed to parse: {e}", "WARNING")
            continue
        
        if system.n_images < 2:
            print_status(f"  Skipping: only {system.n_images} image(s)", "WARNING")
            continue
        
        print_status(f"  System {system.system_id}: {system.n_images} images, "
                    f"{system.get_image_pairs()} pairs")
        
        # Main analysis
        results = analyze_system(
            system,
            detrend_window=args.detrend_window,
            tau_values=tau_values,
            lag_step=args.lag_step,
            mode_lock_window=args.mode_lock_window,
            min_variance_fraction=args.min_variance_fraction,
            min_correlation=args.min_correlation,
            estimator=args.estimator,
            broadband_estimator=(args.broadband_estimator or None),
            no_bandpass=args.no_bandpass,
        )
        
        # Bootstrap uncertainty
        print_status(f"  Running {args.n_bootstrap} bootstrap iterations...")
        bootstrap = run_bootstrap_uncertainty(
            system,
            n_bootstrap=args.n_bootstrap,
            detrend_window=args.detrend_window,
            tau_values=tau_values,
            lag_range=(-200, 200),
            lag_step=args.lag_step,
            mode_lock_window=args.mode_lock_window,
            min_variance_fraction=args.min_variance_fraction,
            min_correlation=args.min_correlation,
            estimator=args.estimator,
            broadband_estimator=(args.broadband_estimator or None),
            bootstrap_mode=args.bootstrap_mode,
            no_bandpass=args.no_bandpass,
        )
        results["bootstrap"] = bootstrap

        # If bootstrap is well-populated, prefer bootstrap-based uncertainty for Γ.
        # This makes the reported significance conservative and robust to non-Gaussianity.
        for pair_key, pair_data in results["pairs"].items():
            if include_pairs and pair_key not in include_pairs:
                continue
            boot = bootstrap.get(pair_key, {})
            phot_ok = boot.get("n_valid", 0) >= max(30, int(0.5 * args.n_bootstrap))
            jack_ok = boot.get("n_seasons", 0) >= 2 and boot.get("n_valid_jackknife", 0) >= 2
            if (phot_ok or jack_ok) and np.isfinite(boot.get("gamma_std", np.nan)):
                gamma_val = pair_data["gamma"]["value"]
                gamma_err = float(boot["gamma_std"])
                pair_data["gamma"]["uncertainty"] = gamma_err
                pair_data["gamma"]["uncertainty_method"] = "bootstrap"
                if np.isfinite(gamma_val) and gamma_err > 0:
                    gamma_sigma = abs(gamma_val) / gamma_err
                    gamma_pvalue = 2 * (1 - stats.norm.cdf(gamma_sigma))
                else:
                    gamma_sigma = np.nan
                    gamma_pvalue = np.nan
                pair_data["gamma"]["sigma"] = gamma_sigma
                pair_data["gamma"]["p_value"] = gamma_pvalue
        
        all_results["systems"][system.system_id] = results
        
        # Collect gamma summary
        for pair_key, pair_data in results["pairs"].items():
            if include_pairs and pair_key not in include_pairs:
                continue
            gamma_data = pair_data["gamma"]
            gamma_summary.append({
                "system": system.system_id,
                "pair": pair_key,
                "gamma": gamma_data["value"],
                "gamma_err": gamma_data["uncertainty"],
                "gamma_sigma": gamma_data["sigma"],
                "gamma_pvalue": gamma_data["p_value"],
                "broadband_delay": pair_data["broadband"]["delay_days"],
                "r_squared": gamma_data.get("r_squared"),
                "passes_practical_cuts": bool(gamma_data.get("passes_practical_cuts", False)),
            })
    
    # Process multi-band CSV files
    for csv_file in csv_files:
        print_status(f"Processing {csv_file.name}")
        
        try:
            system = parse_multiband_csv(csv_file)
        except Exception as e:
            print_status(f"  Failed to parse: {e}", "WARNING")
            continue
        
        if system is None:
            print_status(f"  Skipping: failed to parse CSV (insufficient data or unsupported format)", "WARNING")
            continue
        
        if system.n_images < 2:
            print_status(f"  Skipping: only {system.n_images} image(s) found (need >= 2)", "WARNING")
            continue
        
        print_status(f"  System {system.system_id}: {system.n_images} images, "
                    f"{system.get_image_pairs()} pairs")
        
        # Main analysis
        results = analyze_system(
            system,
            detrend_window=args.detrend_window,
            tau_values=tau_values,
            lag_step=args.lag_step,
            mode_lock_window=args.mode_lock_window,
            min_variance_fraction=args.min_variance_fraction,
            min_correlation=args.min_correlation,
            estimator=args.estimator,
            broadband_estimator=(args.broadband_estimator or None),
            no_bandpass=args.no_bandpass,
        )
        
        # Bootstrap uncertainty
        print_status(f"  Running {args.n_bootstrap} bootstrap iterations...")
        bootstrap = run_bootstrap_uncertainty(
            system,
            n_bootstrap=args.n_bootstrap,
            detrend_window=args.detrend_window,
            tau_values=tau_values,
            lag_range=(-200, 200),
            lag_step=args.lag_step,
            mode_lock_window=args.mode_lock_window,
            min_variance_fraction=args.min_variance_fraction,
            min_correlation=args.min_correlation,
            estimator=args.estimator,
            broadband_estimator=(args.broadband_estimator or None),
            bootstrap_mode=args.bootstrap_mode,
            no_bandpass=args.no_bandpass,
        )
        results["bootstrap"] = bootstrap

        # Update gamma uncertainties from bootstrap
        for pair_key, pair_data in results["pairs"].items():
            if include_pairs and pair_key not in include_pairs:
                continue
            boot = bootstrap.get(pair_key, {})
            phot_ok = boot.get("n_valid", 0) >= max(30, int(0.5 * args.n_bootstrap))
            jack_ok = boot.get("n_seasons", 0) >= 2 and boot.get("n_valid_jackknife", 0) >= 2
            if (phot_ok or jack_ok) and np.isfinite(boot.get("gamma_std", np.nan)):
                gamma_val = pair_data["gamma"]["value"]
                gamma_err = float(boot["gamma_std"])
                pair_data["gamma"]["uncertainty"] = gamma_err
                pair_data["gamma"]["uncertainty_method"] = "bootstrap"
                if np.isfinite(gamma_val) and gamma_err > 0:
                    gamma_sigma = abs(gamma_val) / gamma_err
                    gamma_pvalue = 2 * (1 - stats.norm.cdf(gamma_sigma))
                else:
                    gamma_sigma = np.nan
                    gamma_pvalue = np.nan
                pair_data["gamma"]["sigma"] = gamma_sigma
                pair_data["gamma"]["p_value"] = gamma_pvalue
        
        all_results["systems"][system.system_id] = results
        
        # Collect gamma summary
        for pair_key, pair_data in results["pairs"].items():
            if include_pairs and pair_key not in include_pairs:
                continue
            gamma_data = pair_data["gamma"]
            gamma_summary.append({
                "system": system.system_id,
                "pair": pair_key,
                "gamma": gamma_data["value"],
                "gamma_err": gamma_data["uncertainty"],
                "gamma_sigma": gamma_data["sigma"],
                "gamma_pvalue": gamma_data["p_value"],
                "broadband_delay": pair_data["broadband"]["delay_days"],
                "r_squared": gamma_data.get("r_squared"),
                "passes_practical_cuts": bool(gamma_data.get("passes_practical_cuts", False)),
            })
    
    # Summary statistics
    valid_gammas = [g for g in gamma_summary if np.isfinite(g["gamma"])]
    
    if valid_gammas:
        gamma_values = [g["gamma"] for g in valid_gammas]
        gamma_sigmas_raw = [g["gamma_sigma"] for g in valid_gammas if np.isfinite(g["gamma_sigma"])]
        gamma_sigmas = [
            g["gamma_sigma"]
            for g in valid_gammas
            if bool(g.get("passes_practical_cuts", False)) and np.isfinite(g["gamma_sigma"])
        ]
        
        all_results["summary"] = {
            "n_systems": len(all_results["systems"]),
            "n_pairs_analyzed": len(gamma_summary),
            "n_pairs_valid": len(valid_gammas),
            "gamma_mean": float(np.mean(gamma_values)),
            "gamma_std": float(np.std(gamma_values)),
            "gamma_median": float(np.median(gamma_values)),
            "max_gamma_sigma": float(max(gamma_sigmas)) if gamma_sigmas else np.nan,
            "n_significant_2sigma_raw": sum(1 for s in gamma_sigmas_raw if s > 2),
            "n_significant_3sigma_raw": sum(1 for s in gamma_sigmas_raw if s > 3),
            "n_significant_2sigma": sum(1 for s in gamma_sigmas if s > 2),
            "n_significant_3sigma": sum(1 for s in gamma_sigmas if s > 3),
        }
        
        # TEP-GL interpretation
        n_sig = all_results["summary"]["n_significant_2sigma"]
        if n_sig == 0:
            interpretation = "NULL: No significant temporal shear detected. Consistent with GR (single constant delay)."
        elif n_sig == 1:
            interpretation = "MARGINAL: One pair shows >2σ temporal shear. Could be statistical fluctuation or systematics."
        else:
            interpretation = f"CANDIDATE: {n_sig} pairs show >2σ temporal shear. Warrants further investigation."
        
        all_results["summary"]["interpretation"] = interpretation
    else:
        all_results["summary"] = {
            "n_systems": len(all_results["systems"]),
            "n_pairs_analyzed": len(gamma_summary),
            "n_pairs_valid": 0,
            "interpretation": "INSUFFICIENT DATA: No valid Gamma measurements obtained.",
        }
    
    # Save results
    output_file = output_dir / f"step_3_0_cosmograil_temporal_shear{tag}.json"
    with open(output_file, "w") as f:
        json.dump(all_results, f, indent=2, default=lambda x: None if not np.isfinite(x) else x)
    print_status(f"Saved results: {output_file}", "SUCCESS")
    
    # Print summary
    print("\n" + "="*60)
    print("COSMOGRAIL TEMPORAL SHEAR ANALYSIS SUMMARY")
    print("="*60)
    print(f"Systems analyzed: {all_results['summary']['n_systems']}")
    print(f"Image pairs analyzed: {all_results['summary']['n_pairs_analyzed']}")
    print(f"Valid Γ measurements: {all_results['summary'].get('n_pairs_valid', 0)}")
    
    if "gamma_mean" in all_results["summary"]:
        print(f"\nΓ (temporal shear slope) statistics:")
        print(f"  Mean:   {all_results['summary']['gamma_mean']:.3f} days/decade")
        print(f"  Median: {all_results['summary']['gamma_median']:.3f} days/decade")
        print(f"  Std:    {all_results['summary']['gamma_std']:.3f} days/decade")
        print(f"\nSignificance:")
        print(f"  >2σ (raw): {all_results['summary'].get('n_significant_2sigma_raw', 0)} pairs")
        print(f"  >3σ (raw): {all_results['summary'].get('n_significant_3sigma_raw', 0)} pairs")
        print(f"  >2σ (practical cuts): {all_results['summary']['n_significant_2sigma']} pairs")
        print(f"  >3σ (practical cuts): {all_results['summary']['n_significant_3sigma']} pairs")
    
    print(f"\nInterpretation: {all_results['summary']['interpretation']}")
    print("="*60)
    
    print_status("Step 3.0 complete", "SUCCESS")


if __name__ == "__main__":
    # Setup file logging when run manually
    import sys
    from pathlib import Path
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root))
    sys.path.insert(0, str(repo_root / "scripts"))
    from utils.logger import setup_step_logger
    logger = setup_step_logger("step_3_0_cosmograil_temporal_shear")
    main()
