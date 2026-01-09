import json
import numpy as np

def check_quality():
    with open('results/outputs/step_3_0_cosmograil_temporal_shear_v3_expanded.json', 'r') as f:
        data = json.load(f)
        
    print(f"{'System':<10} {'Pair':<5} {'Gamma':<8} {'MaxJump':<8} {'Linearity(R2)':<10} {'Verdict'}")
    print("-" * 65)
    
    detections = []
    
    for sys_id, sys_data in data['systems'].items():
        for pair_id, pair_data in sys_data['pairs'].items():
            gamma_val = pair_data['gamma']['value']
            if gamma_val is None or not np.isfinite(gamma_val):
                continue
                
            sigma = pair_data['gamma']['sigma']
            is_detection = sigma > 3.0
            
            ms = pair_data['multiscale']
            sorted_taus = sorted([int(k) for k in ms.keys()])
            delays = []
            for t in sorted_taus:
                d = ms[str(t)].get('delay_days')
                if d is not None and np.isfinite(d):
                    delays.append(d)
            
            if len(delays) < 2:
                continue
                
            # Calc max jump between consecutive scales
            delays = np.array(delays)
            jumps = np.abs(np.diff(delays))
            max_jump = np.max(jumps) if len(jumps) > 0 else 0
            
            r2 = pair_data['gamma']['r_squared']
            
            # Verdict
            verdict = "OK"
            if max_jump > 50:
                verdict = "JUMP"
            
            if is_detection:
                row = f"{sys_id:<10} {pair_id:<5} {gamma_val:<8.1f} {max_jump:<8.1f} {r2:<10.3f} {verdict}"
                detections.append(row)
            
    print("\nSIGNIFICANT DETECTIONS (>3 sigma):")
    for row in detections:
        print(row)

if __name__ == "__main__":
    check_quality()
