import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from astropy.cosmology import Planck18 as cosmo

def run_microlensing_scaling():
    """
    Calculates the microlensing optical depth proxy (Kappa_eff) for the lens sample
    and tests its correlation with source redshift.
    
    Kappa_eff ~ (Sigma_star / Sigma_crit)
    Sigma_crit = (c^2 / 4pi G) * (Ds / (Dd * Dds))
    So Kappa_eff propto (Dd * Dds) / Ds  (assuming Sigma_star is roughly constant or uncorrelated)
    
    Actually, optical depth tau ~ integral(rho) dl.
    For a singular isothermal sphere (SIS), tau ~ 1/2 everywhere inside Einstein radius?
    But for microlensing specifically, it depends on the surface mass density of stars.
    
    We use the geometric efficiency factor D_d * D_ds / D_s as a proxy for the 'lensing efficiency'
    distance weighting, though the stellar density Sigma_star is the dominant physical variable.
    
    However, the counter-argument is that higher z_s means longer path, but microlensing is local to the lens.
    The optical depth for microlensing is proportional to the surface mass density in stars at the image position.
    
    Let's check if the geometric factor correlates with z_s.
    """
    print("--- Step 3.18: Microlensing Redshift Scaling ---")
    
    # System: (z_l, z_s)
    systems = [
        {"name": "Q2237",    "zl": 0.039, "zs": 1.695},
        {"name": "RXJ1131",  "zl": 0.295, "zs": 0.658},
        {"name": "HS2209",   "zl": 0.280, "zs": 1.070},
        {"name": "PG1115",   "zl": 0.311, "zs": 1.722},
        {"name": "J1001",    "zl": 0.415, "zs": 1.838},
        {"name": "HE0435",   "zl": 0.454, "zs": 1.693},
        {"name": "DESJ0408", "zl": 0.597, "zs": 2.375},
        {"name": "WFI2033",  "zl": 0.661, "zs": 1.662},
        {"name": "HE1104",   "zl": 0.729, "zs": 2.319},
        {"name": "J1206",    "zl": 0.745, "zs": 1.789}
    ]
    
    z_s = []
    kappa_proxy = []
    
    for sys in systems:
        zl = sys['zl']
        zs = sys['zs']
        
        # Angular diameter distances (Mpc)
        Dd = cosmo.angular_diameter_distance(zl).value
        Ds = cosmo.angular_diameter_distance(zs).value
        Dds = cosmo.angular_diameter_distance_z1z2(zl, zs).value
        
        # Sigma_crit inverse propto Dd * Dds / Ds
        # Efficiency ~ Dd * Dds / Ds
        eff = (Dd * Dds) / Ds
        
        z_s.append(zs)
        kappa_proxy.append(eff)
        
    z_s = np.array(z_s)
    kappa_proxy = np.array(kappa_proxy)
    
    # Correlation
    r, p = stats.pearsonr(z_s, kappa_proxy)
    print(f"Correlation (z_s vs Efficiency): r={r:.3f}, p={p:.3e}")
    
    # Plot
    plt.figure(figsize=(6, 5))
    plt.scatter(z_s, kappa_proxy, c='black', alpha=0.7)
    
    # Linear fit for display
    m, b = np.polyfit(z_s, kappa_proxy, 1)
    plt.plot(z_s, m*z_s + b, 'r--', alpha=0.5, label=f'r={r:.2f}')
    
    for i, txt in enumerate(systems):
        plt.annotate(txt['name'], (z_s[i], kappa_proxy[i]), fontsize=8, alpha=0.7, xytext=(3,3), textcoords='offset points')
        
    plt.xlabel('Source Redshift $z_s$')
    plt.ylabel(r'Lensing Efficiency Proxy ($D_d D_{ds} / D_s$)')
    plt.title('Microlensing Geometric Scaling Check')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    os.makedirs('site/figures', exist_ok=True)
    plt.savefig('site/figures/microlensing_redshift_scaling.png', dpi=300)
    print("Saved site/figures/microlensing_redshift_scaling.png")

if __name__ == "__main__":
    import os
    run_microlensing_scaling()
