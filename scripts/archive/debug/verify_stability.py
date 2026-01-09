import json
import numpy as np
from scipy import stats

def verify_stability():
    """
    Verify the stability of reported detections by excluding small-scale smoothing (tau=5, 10)
    to ensure the signal is not driven by mode-jumping artifacts like Q0957.
    """
    with open('results/outputs/step_3_0_cosmograil_temporal_shear_v3_expanded.json', 'r') as f:
        data = json.load(f)
    
    # Systems to check (The claimed detections)
    targets = [
        ('DESJ0408', 'A-D', -333),
        ('DESJ0408', 'B-D', -129),
        ('PG1115', 'B-C', -207),
        ('PG1115', 'A-B', 156),
        ('J1206', 'A-B', -103)
    ]
    
    print(f"{'System':<10} {'Pair':<5} {'Gamma_Full':<10} {'Gamma_Robust':<10} {'Status'}")
    print("-" * 60)
    
    for sys_id, pair_id, expected_gamma in targets:
        pair_data = data['systems'][sys_id]['pairs'][pair_id]
        ms = pair_data['multiscale']
        
        # Collect all points
        taus_all = []
        delays_all = []
        
        # Collect robust points (tau >= 20)
        taus_rob = []
        delays_rob = []
        
        for t_str, res in ms.items():
            t = int(t_str)
            d = res.get('delay_days')
            
            if d is not None and np.isfinite(d):
                taus_all.append(t)
                delays_all.append(d)
                
                if t >= 20:
                    taus_rob.append(t)
                    delays_rob.append(d)
        
        # Fit Full
        if len(taus_all) >= 2:
            res_full = stats.linregress(np.log10(taus_all), delays_all)
            g_full = res_full.slope
        else:
            g_full = np.nan
            
        # Fit Robust
        if len(taus_rob) >= 2:
            res_rob = stats.linregress(np.log10(taus_rob), delays_rob)
            g_rob = res_rob.slope
        else:
            g_rob = np.nan
            
        # Check consistency
        # We look for sign flip or massive change in magnitude
        consistent = False
        if np.isfinite(g_full) and np.isfinite(g_rob):
            # Sign consistency
            if np.sign(g_full) == np.sign(g_rob):
                consistent = True
            
            # Magnitude consistency (factor of 2?)
            # Or just check if robust is still "large"
            if abs(g_rob) < 20: # If it collapses to near-zero
                consistent = False
                
        status = "STABLE" if consistent else "UNSTABLE"
        
        print(f"{sys_id:<10} {pair_id:<5} {g_full:<10.1f} {g_rob:<10.1f} {status}")

if __name__ == "__main__":
    verify_stability()
