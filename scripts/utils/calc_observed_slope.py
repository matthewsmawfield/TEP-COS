
import json
import numpy as np
from scipy import stats

def calculate_observed_slope():
    # Load Observed Residuals
    with open('results/outputs/step_5_31_per_cluster_controlled_residuals.json', 'r') as f:
        obs_data = json.load(f)
    
    # Load Cluster Params (for density)
    with open('results/outputs/step_5_32_full_density_scaling.json', 'r') as f:
        sim_data = json.load(f)
        
    # Map name to density
    rho_map = {c['name']: c['rho_c_log'] for c in sim_data['clusters']}
    
    # Collect X (rho) and Y (residual)
    x = []
    y = []
    
    print(f"{'Cluster':<20} {'Rho':<6} {'Res':<6}")
    print("-" * 35)
    
    for name, data in obs_data['clusters'].items():
        # Match names (approximate)
        rho = None
        for k, v in rho_map.items():
            if name in k or k in name:
                rho = v
                break
        
        if rho is not None:
            res = data['controlled_residual']
            x.append(rho)
            y.append(res)
            print(f"{name:<20} {rho:<6.2f} {res:+.3f}")
            
    # Linear Regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    
    print("-" * 35)
    print(f"Observed Slope: {slope:.4f}")
    print(f"Observed R:     {r_value:.4f}")
    print(f"Observed P:     {p_value:.4e}")
    print(f"N Clusters:     {len(x)}")

if __name__ == "__main__":
    calculate_observed_slope()
