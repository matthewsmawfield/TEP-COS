import json
import numpy as np
from scipy import stats

# Load Observed Data
with open('results/outputs/step_5_31_per_cluster_controlled_residuals.json', 'r') as f:
    obs_data = json.load(f)

# Load Newtonian Data
with open('results/outputs/step_5_32_full_density_scaling.json', 'r') as f:
    newton_data = json.load(f)

# 1. Observed Slope
obs_clusters = obs_data['clusters']
# We need to map cluster names to densities.
# Use densities from Newtonian simulation which has the mapping
rho_map = {c['name']: c['rho_c_log'] for c in newton_data['clusters']}

obs_x = []
obs_y = []

print("Cluster | Log(rho) | Residual")
for name, data in obs_clusters.items():
    if name in rho_map:
        rho = rho_map[name]
        res = data['controlled_residual']
        obs_x.append(rho)
        obs_y.append(res)
        print(f"{name:<20} | {rho:.2f} | {res:+.3f}")

slope_obs, intercept_obs, r_obs, p_obs, err_obs = stats.linregress(obs_x, obs_y)

print("-" * 30)
print(f"Observed Slope: {slope_obs:.3f} ± {err_obs:.3f}")
print(f"Observed R:     {r_obs:.3f} (p={p_obs:.4f})")
print("-" * 30)
print(f"Newtonian Slope: {newton_data['slope']:.3f}")
print("-" * 30)
