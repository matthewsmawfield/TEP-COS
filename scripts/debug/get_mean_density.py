import pandas as pd
import numpy as np

# Load data same as step_5_33
csv_path = "results/outputs/step_5_10_pulsar_population_controls.csv"
df = pd.read_csv(csv_path)
gc_df = df[df['environment'] == 'globular_cluster'].copy()

CLUSTER_DENSITIES = {
    "Terzan 5": 5.50, "47 Tuc (NGC 104)": 4.88, "NGC 6517": 5.80,
    "M28 (NGC 6626)": 4.52, "M62 (NGC 6266)": 5.16, "M13 (NGC 6205)": 3.79,
    "M15 (NGC 7078)": 5.05, "M5 (NGC 5904)": 3.53, "Terzan 1": 5.00,
    "NGC 6752": 4.30, "M2 (NGC 7089)": 4.15, "Omega Centauri (NGC 5139)": 3.12,
    "M53 (NGC 5024)": 2.96, "M3 (NGC 5272)": 3.68, "M71 (NGC 6838)": 2.29,
    "NGC 6397": 5.68, "NGC 1851": 5.09, "NGC 6522": 5.50,
    "NGC 6544": 5.20, "NGC 6624": 5.60, "NGC 6760": 3.80,
    "M22 (NGC 6656)": 2.97, "M80 (NGC 6093)": 4.79, "M92 (NGC 6341)": 4.30,
    "NGC 6712": 3.70, "NGC 6652": 4.50, "M14 (NGC 6402)": 3.44,
    "NGC 6539": 3.30, "M4 (NGC 6121)": 2.85
}

gc_df['log_rho_c'] = gc_df['cluster'].map(CLUSTER_DENSITIES)
gc_df = gc_df.dropna(subset=['log_rho_c', 'logPdot_abs', 'logP', 'log_b_proxy'])

counts = gc_df['cluster'].value_counts()
valid_clusters = counts[counts >= 3].index
gc_df = gc_df[gc_df['cluster'].isin(valid_clusters)]

mean_log_rho = gc_df['log_rho_c'].mean()
print(f"Mean log_rho_c: {mean_log_rho}")
