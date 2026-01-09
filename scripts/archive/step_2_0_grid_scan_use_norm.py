#!/usr/bin/env python3

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def _parse_float_list(s: str) -> list[float]:
    out: list[float] = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        out.append(float(part))
    return out


def _parse_windows(s: str) -> list[tuple[float, float]]:
    out: list[tuple[float, float]] = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        lo_s, hi_s = part.split("-")
        out.append((float(lo_s), float(hi_s)))
    return out


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-perm", type=int, default=2000)
    parser.add_argument("--n-axis-rand", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--delta-phi-deg", type=str, default="10,15,20,30")
    parser.add_argument("--rre-windows", type=str, default="0.7-1.3,0.8-1.2,0.9-1.1")
    parser.add_argument("--gas-line", type=str, default="Ha-6564")
    parser.add_argument("--out-csv", type=str, default="results/outputs/step_2_0_grid_scan_use_norm_summary.csv")

    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    step_script = repo_root / "scripts" / "steps" / "step_2_0_cosmic_coriolis_analysis.py"
    if not step_script.exists():
        raise FileNotFoundError(step_script)

    dphis = _parse_float_list(args.delta_phi_deg)
    windows = _parse_windows(args.rre_windows)

    out_csv = (repo_root / args.out_csv).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []

    tracers: list[tuple[str, list[str]]] = [
        ("stellar", ["--velocity-source", "stellar"]),
        ("gas", ["--velocity-source", "gas", "--gas-line", str(args.gas_line)]),
    ]

    for tracer_name, tracer_args in tracers:
        for dphi in dphis:
            for rre_min, rre_max in windows:
                tag = f"norm_{tracer_name}_p{int(round(dphi))}_r{int(round(100*rre_min)):02d}{int(round(100*rre_max)):02d}"

                cmd = [
                    sys.executable,
                    str(step_script),
                    "--use-norm",
                    "--output-tag",
                    tag,
                    "--delta-phi-deg",
                    str(float(dphi)),
                    "--rre-min",
                    str(float(rre_min)),
                    "--rre-max",
                    str(float(rre_max)),
                    "--n-perm",
                    str(int(args.n_perm)),
                    "--n-axis-rand",
                    str(int(args.n_axis_rand)),
                    "--seed",
                    str(int(args.seed)),
                ] + tracer_args

                print(f"[grid] {tracer_name} dphi={dphi} rre=[{rre_min},{rre_max}] tag={tag}")
                subprocess.run(cmd, cwd=str(repo_root), check=True)

                js = repo_root / "results" / "outputs" / f"step_2_0_cosmic_coriolis_summary_{tag}.json"
                if not js.exists():
                    raise FileNotFoundError(js)

                s = _read_json(js)

                def _p(d: dict | None) -> float:
                    if not isinstance(d, dict):
                        return float("nan")
                    v = d.get("p_value")
                    return float(v) if v is not None else float("nan")

                axis_rand = s.get("axis_randomization") if isinstance(s, dict) else None
                axis_rand_legacy = _p(axis_rand.get("legacy") if isinstance(axis_rand, dict) else None)
                axis_rand_axis = _p(axis_rand.get("axis") if isinstance(axis_rand, dict) else None)
                axis_rand_wedge = _p(axis_rand.get("wedge") if isinstance(axis_rand, dict) else None)

                perm = s.get("permutation")
                perm_axis = s.get("permutation_axis")
                perm_wedge = s.get("permutation_wedge")

                fit = s.get("fit")
                fit_axis = s.get("fit_axis")
                fit_wedge = s.get("fit_wedge")

                rows.append(
                    {
                        "tracer": tracer_name,
                        "gas_line": s.get("gas_line"),
                        "use_norm": bool(s.get("use_norm")),
                        "delta_phi_deg": float(s.get("delta_phi_deg")),
                        "rre_min": float(s.get("rre_min")),
                        "rre_max": float(s.get("rre_max")),
                        "n_galaxies": int(s.get("n_galaxies")),
                        "n_wedge": int(fit_wedge.get("n")) if isinstance(fit_wedge, dict) and fit_wedge.get("n") is not None else 0,
                        "a_legacy": float(fit.get("a")) if isinstance(fit, dict) else float("nan"),
                        "p_pair_legacy": float(perm.get("p_value_pair")) if isinstance(perm, dict) else float("nan"),
                        "a_axis": float(fit_axis.get("a")) if isinstance(fit_axis, dict) else float("nan"),
                        "p_pair_axis": float(perm_axis.get("p_value_pair")) if isinstance(perm_axis, dict) else float("nan"),
                        "a_wedge": float(fit_wedge.get("a")) if isinstance(fit_wedge, dict) else float("nan"),
                        "p_pair_wedge": float(perm_wedge.get("p_value_pair")) if isinstance(perm_wedge, dict) else float("nan"),
                        "p_axis_rand_legacy": axis_rand_legacy,
                        "p_axis_rand_axis": axis_rand_axis,
                        "p_axis_rand_wedge": axis_rand_wedge,
                        "output_tag": tag,
                    }
                )

    fieldnames = list(rows[0].keys()) if rows else []
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote grid summary: {out_csv}")


if __name__ == "__main__":
    main()
