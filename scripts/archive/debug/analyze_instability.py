import json
import numpy as np

def analyze_instability():
    with open('results/outputs/step_3_0_cosmograil_temporal_shear_v3_expanded.json', 'r') as f:
        data = json.load(f)
    
    # Systems to check
    targets = [
        ('DESJ0408', 'A-D'),
        ('DESJ0408', 'B-D'),
        ('PG1115', 'B-C'),
        ('PG1115', 'A-B'),
        ('J1206', 'A-B')
    ]
    
    print(f"{'System':<10} {'Pair':<5} {'Gamma':<8} {'Tau vs Delay'}")
    print("-" * 80)
    
    for sys_id, pair_id in targets:
        pair_data = data['systems'][sys_id]['pairs'][pair_id]
        gamma = pair_data['gamma']['value']
        
        ms = pair_data['multiscale']
        taus = sorted([int(k) for k in ms.keys()])
        
        print(f"{sys_id:<10} {pair_id:<5} {gamma:<8.1f}")
        for t in taus:
            d = ms[str(t)].get('delay_days')
            unc = ms[str(t)].get('uncertainty_days')
            if d is not None:
                print(f"  Tau={t:<3}: Delay={d:>8.2f} +/- {unc:>6.2f}")
            else:
                print(f"  Tau={t:<3}: Delay=    None")
        print("")

if __name__ == "__main__":
    analyze_instability()
