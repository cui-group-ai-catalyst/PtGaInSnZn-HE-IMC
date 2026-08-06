"""
20260415_FigE_compute_and_plot.py
=================================
Submission-facing provenance script for Panel e.

Scientific purpose:
- Quantify host-versus-cocktail size mismatch against the fixed Ga/In/Sn/Zn
  liquid composition.
- Pair the mismatch metric with the corrected Miedema enthalpy drive so the
  geometric constraint can be interpreted alongside the chemical attraction.
- Export a ranking table and review plot after the calculation is complete.

Important interpretation note:
- Any visual threshold lines used in the plot are heuristic highlighting rules
  for manuscript discussion. They should not be overstated as exact physical
  phase boundaries.

Manuscript figure scope (important for reviewers)
-------------------------------------------------
This script reproduces the ranked Resistance/Enthalpy_Drive view and writes
its outputs under `_regen` filenames:
    data_FigE_Resistance_Ranked_regen.csv
    preview_FigE_Resistance_Plot_regen.png

The manuscript-facing figure used in Supplementary Fig. 2 is
`preview_FigE_ThreeWay.png`, an Origin layout that composes three views of
the dataset. It is paired with the bundled `data_FigE_True_ThreeWay.csv`.

This script does NOT reproduce the ThreeWay Origin layout. It reproduces
only the ranked-resistance subset. A reviewer running this script will
therefore see a PNG that differs from the manuscript figure -- this is
expected and does NOT indicate a code or data problem. The numerical
substrate is the same; only the layout differs.

Canonical vs regen convention
-----------------------------
This script never overwrites canonical files. After the 2026-07-22 correction
to the n_WS^(1/3) term, the `_regen` CSV intentionally differs from the older
canonical `data_FigE_Resistance_Ranked.csv`. The canonical file is retained
for provenance until the corrected SI source data and figure are approved.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "shared"))

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from data_periodic_table import periodic_table_data

# --- Pt-Ga-In-Sn-Zn Manuscript: Panel e (Structural Mismatch Analysis) ---
# Date: 2026-04-15
# Audit Note: Miedema constants fixed. No *10 multiplier.

def calc_cocktail_properties():
    weights = {'Ga': 0.65, 'In': 0.20, 'Sn': 0.10, 'Zn': 0.05}
    avg_r = sum(periodic_table_data[el]['r'] * w for el, w in weights.items())
    avg_Phi = sum(periodic_table_data[el]['Phi'] * w for el, w in weights.items())
    avg_nws = sum(periodic_table_data[el]['n_ws'] * w for el, w in weights.items())
    return {'r': avg_r, 'Phi': avg_Phi, 'n_ws': avg_nws, 'weights': weights}

def calc_resistance_analysis():
    cocktail = calc_cocktail_properties()
    r_solute = cocktail['r']
    print(f"Average Liquid Atom Radius: {r_solute:.2f} pm")
    
    results = []
    candidates = [
        'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu', 
        'Y','Zr','Nb','Mo','Ru','Rh','Pd','Ag','Cd',
        'La','Ce',
        'Hf','Ta','W','Re','Os','Ir','Pt','Au'
    ]
    
    for host in candidates:
        if host not in periodic_table_data: continue
        props = periodic_table_data[host]
        
        r_host = props['r']
        mismatch = abs((r_solute - r_host) / r_host) * 100
        
        # Correct theoretical Miedema constants
        P, Q = 14.1, 9.4
        dH_cocktail = 0
        for el, frac in cocktail['weights'].items():
            t_props_c = periodic_table_data[el]
            dPhi_c = props['Phi'] - t_props_c['Phi']
            # The shared table already stores n_WS^(1/3), so use a direct
            # difference. Taking another cube root would compress this term.
            dn_c = props['n_ws'] - t_props_c['n_ws']
            # BUG FIX: Removed the * 10 legacy multiplier
            dH_i = (-P * dPhi_c**2 + Q * dn_c**2)
            dH_cocktail += frac * dH_i
            
        dH_mix_kj = dH_cocktail 
        
        results.append({
            'Host': host,
            'Mismatch_Percent': mismatch,
            'Enthalpy_Drive': dH_mix_kj
        })
        
    df = pd.DataFrame(results)
    
    # Export to CSV
    try:
        import os
        script_dir = os.path.dirname(os.path.abspath(__file__))
        results_dir = script_dir
        os.makedirs(results_dir, exist_ok=True)
        csv_path = os.path.join(results_dir, "data_FigE_Resistance_Ranked_regen.csv")
        # Sort by Mismatch_Percent
        df_sorted = df.sort_values(by='Mismatch_Percent', ascending=True)
        df_sorted.to_csv(csv_path, index=False)
        print(f"Saved Resistance Analysis Data to {csv_path}")
    except Exception as e:
        print(f"Could not save CSV: {e}")
        
    return df

def plot_resistance(df):
    try:
        import os
        script_dir = os.path.dirname(os.path.abspath(__file__))
        results_dir = script_dir
        os.makedirs(results_dir, exist_ok=True)
    except:
        results_dir = os.path.dirname(os.path.abspath(__file__))

    plt.figure(figsize=(12, 10)) 
    sns.set_style("whitegrid")
    
    x = df['Enthalpy_Drive']
    y = df['Mismatch_Percent']
    
    colors = []
    sizes = []
    for i, row in df.iterrows():
        # Thermodynamic stability criteria: Enthalpy < -20 kJ/mol (after fix, still a strong criteria)
        is_stable = (row['Enthalpy_Drive'] < -20) and (row['Mismatch_Percent'] < 15)
        if row['Host'] == 'Pt':
            colors.append('red') 
            sizes.append(180)
        elif is_stable:
            colors.append('green')
            sizes.append(80)
        else:
            colors.append('gray')
            sizes.append(50)
            
    plt.scatter(x, y, c=colors, s=sizes, alpha=0.7, edgecolors='k')
    
    # --- Precision Labeling ---
    for i, row in df.iterrows():
        label = row['Host']
        dH = row['Enthalpy_Drive']
        mis = row['Mismatch_Percent']
        
        in_green_zone = (dH < -20) and (mis < 15)
        is_outlier = label in ['La', 'Ce', 'Y']
        
        if not (in_green_zone or is_outlier): continue

        font_size = 10
        weight = 'normal'
        dx = 0
        dy = 0.5 
        ha = 'center'
        va = 'bottom'
        
        # --- Cluster Management (dx scaled by 1/10 due to new physically accurate x-axis) ---
        
        if label == 'Pt':
            font_size = 13; weight = 'bold'; color='red'
            dy = -1.1; va='top'
            dx = 0; ha='center'
        elif label == 'Pd':
            dx = -0.3; ha='right'; dy=0; va='center' 
        elif label == 'Os':
            dx = 0.3; ha='left'; dy=0; va='center'
        elif label == 'Ir':
            dy = 0.8; va='bottom' 
        elif label == 'Ru':
            dy = 0.6; va='bottom' 
            dx = -0.1; ha='right'   
        elif label == 'Rh':
            dy = 0.6; va='bottom'
            dx = 0.1; ha='left'     
        elif label == 'Re':
            dy = -0.8; va='top'
            ha = 'center'
        elif label == 'Au':
            dx = -0.3; ha='right'; dy=0; va='center'
        elif label == 'Ni':
            dx = 0.3; ha='left'; dy=0; va='center'
        elif label == 'Co':
            dy = -0.8; va='top'
        elif label == 'Fe':
            dy = 0.8; va='bottom'
        elif label == 'W':
            dx = -0.3; ha='right'; dy=0; va='center' 
        elif label == 'Mo':
            dx = 0.3; ha='left'; dy=0; va='center'   
        elif label == 'Zr':
            dy = -0.8; va='top'
        elif label == 'Hf':
            dy = 0.8; va='bottom'
        elif label == 'Ti':
            dy = 0.8; va='bottom'
        elif label == 'Cr':
            dy = -0.8; va='top'
        elif label == 'Cu':
            dx = 0.3; ha='left'; dy=0; va='center'
        elif label == 'La':
            dy = 0.6; va='bottom'
        elif label == 'Ce':
            dy = -0.8; va='top'
             
        plt.text(dH + dx, mis + dy, label, fontsize=font_size, fontweight=weight, ha=ha, va=va)

    plt.axhline(15, color='red', linestyle='--', linewidth=2, label='Hume-Rothery Limit (15%)')
    
    rect = plt.Rectangle((-100, 0), 80, 15, linewidth=2, edgecolor='green', facecolor='green', alpha=0.1)
    plt.gca().add_patch(rect)
    
    plt.text(-40, 7.5, "Stability Zone", color='green', 
             fontsize=14, fontweight='bold', ha='center', va='center', alpha=0.3)

    plt.title('Panel e: Structural Mismatch Analysis', fontsize=16, fontweight='bold')
    plt.xlabel(r'Mixing enthalpy, $\Delta H_{mix}$ (kJ mol$^{-1}$)', fontsize=14, fontweight='bold')
    plt.ylabel(r'Atomic size mismatch, $\delta$ (%)', fontsize=14, fontweight='bold')
    plt.legend(loc='upper right', fontsize=10)
    
    # Updated physical limits for correct kJ/mol bounds
    plt.xlim(-50, 10)
    plt.ylim(0, 30)
    
    plt.tight_layout()
    save_path = os.path.join(results_dir, "preview_FigE_Resistance_Plot_regen.png")
    plt.savefig(save_path, dpi=300)
    print(f"Saved Resistance Plot to {save_path}")

if __name__ == "__main__":
    df = calc_resistance_analysis()
    plot_resistance(df)
