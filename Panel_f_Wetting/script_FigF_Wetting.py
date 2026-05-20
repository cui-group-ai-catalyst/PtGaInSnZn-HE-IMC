"""
script_FigF_Wetting.py
======================
Panel f — solid-liquid interfacial free energy (γ_SL) across 28 host
candidates wetted by the fixed Ga/In/Sn/Zn liquid cocktail.

Method (matches manuscript Methods eq. 1 and eq. 4):
    γ_SL = (sqrt(γ_S) − sqrt(γ_L))^2 + ΔH_mix × f_surf
        γ_S        = 1.15 × γ_L,host                (Miedema rule)
        ΔH_mix     = Σ_X y_X · ΔH(host, X)          (cocktail mixing enthalpy)
        ΔH(M, X)   = −P (Δφ*)^2 + Q (Δn_WS^(1/3))^2 (simplified Miedema, P=14.1, Q=9.4)
        f_surf     = 1000 / (1.091 · N_A^(1/3) · V_m^(2/3))   (J m⁻² per kJ mol⁻¹)

Element parameters: shared/data_periodic_table.py (φ*, n_WS^(1/3), γ_L).
Molar volumes (cm³/mol): tabulated locally from de Boer 1988.

Convention note: the 'n_ws' field in shared/data_periodic_table.py stores
n_WS^(1/3) directly (Pt = 1.78 = n_WS^(1/3); raw n_WS ≈ 5.6 d.u.). The
Miedema formula calls for Δ(n_WS^(1/3)) so the field is used as-is below.

Outputs (rewritten on each run):
    data_FigF_Wetting_Ranked_regen.csv   28-host γ_SL ranking
    preview_FigF_Wetting_regen.png       ΔH_mix vs γ_SL scatter

Reproducibility: pure NumPy/Pandas; no UMA dependency; deterministic.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "shared"))

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from data_periodic_table import periodic_table_data

# --- Pillar 2: The Enablers (Wetting / Interface Energy) ---

def calc_cocktail_properties():
    weights = {'Ga': 0.65, 'In': 0.20, 'Sn': 0.10, 'Zn': 0.05}
    avg_gamma_L = sum(periodic_table_data[el]['gamma_L'] * w for el, w in weights.items())
    avg_Phi = sum(periodic_table_data[el]['Phi'] * w for el, w in weights.items())
    avg_nws = sum(periodic_table_data[el]['n_ws'] * w for el, w in weights.items())
    return {'gamma_L': avg_gamma_L, 'Phi': avg_Phi, 'n_ws': avg_nws, 'weights': weights}

# Molar volumes (cm3/mol) for f_surf calculation — de Boer 1988 Table 4.1
VM_DATA = {
    'Pt': 9.09, 'Ir': 8.52, 'Rh': 8.28, 'Pd': 8.56, 'Ru': 8.17, 'Os': 8.42,
    'Re': 8.86, 'Au': 10.21, 'W': 9.47, 'Mo': 9.38, 'Ni': 6.59, 'Co': 6.67,
    'Fe': 7.09, 'Cu': 7.11, 'Cr': 7.23, 'Mn': 7.35, 'V': 8.32, 'Ti': 10.64,
    'Sc': 15.00, 'Ag': 10.27, 'Cd': 13.10, 'Y': 19.88, 'Zr': 14.02, 'Nb': 10.83,
    'Hf': 13.44, 'Ta': 10.85, 'La': 22.39, 'Ce': 20.69,
}

def calc_wetting_analysis():
    """
    Corrected Miedema interfacial energy model (de Boer 1988).
    
    gamma_SL = gamma_SL_geom + gamma_SL_chem  [units: J/m2]
    
    gamma_SL_geom = (sqrt(gamma_S) - sqrt(gamma_L_cocktail))^2
        gamma_S = 1.15 * gamma_L_host  (Miedema empirical)
        
    gamma_SL_chem = dH_mix * f_surf
        dH_mix = weighted sum of binary dH(M, L)  [kJ/mol], CORRECTED (no *10)
        f_surf = 1000 / (1.091 * N_A^(1/3) * Vm^(2/3))  [J/m2 per kJ/mol]
        Vm = element-specific molar volume from VM_DATA  [m3/mol]
        
    Reference: de Boer et al., Cohesion in Metals, North-Holland, 1988.
    """
    cocktail = calc_cocktail_properties()
    P, Q = 14.1, 9.4  # Correct Miedema constants (no *10)
    N_A = 6.022e23
    
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
        
        gamma_S = props['gamma_L'] * 1.15
        
        # --- Corrected dH_mix: no *10 ---
        dH_cocktail = 0
        for el, frac in cocktail['weights'].items():
            lp = periodic_table_data[el]
            dPhi_c = props['Phi'] - lp['Phi']
            # NOTE: shared/data_periodic_table.py stores n_WS^(1/3) directly in the
            # 'n_ws' field (Pt = 1.78 = n_WS^(1/3), not the raw electron density
            # ~5.6 d.u.). The Miedema formula calls for Δ(n_WS^(1/3)), so use the
            # tabulated value directly — do NOT take the cube root again.
            dn_c = props['n_ws'] - lp['n_ws']
            dH_i = -P * dPhi_c**2 + Q * dn_c**2   # Miedema: -P (Δφ*)^2 + Q (Δn_WS^(1/3))^2
            dH_cocktail += frac * dH_i
        
        dH_mix_kj = dH_cocktail  # kJ/mol, correct scale
        
        # --- Corrected f_surf: element-specific via Vm (de Boer 1988) ---
        Vm_cm3 = VM_DATA.get(host, 9.0)  # cm3/mol, fallback 9.0
        Vm_m3 = Vm_cm3 * 1e-6            # m3/mol
        f_surf = 1000.0 / (1.091 * (N_A**(1/3)) * (Vm_m3**(2/3)))
        
        interaction_term = dH_mix_kj * f_surf   # J/m2, CORRECTED
        Gamma_SL_geom = (np.sqrt(gamma_S) - np.sqrt(cocktail['gamma_L']))**2
        gamma_sl = Gamma_SL_geom + interaction_term
        
        status = 'Wets' if gamma_sl < 0 else 'Does Not Wet'
        results.append({
            'Host': host,
            'Delta_H_mix': round(dH_mix_kj, 4),
            'f_surf': round(f_surf, 5),
            'Gamma_SL_geom': round(Gamma_SL_geom, 6),
            'Gamma_SL_chem': round(interaction_term, 6),
            'Gamma_SL': round(gamma_sl, 6),
            'Status': status
        })
        
    df = pd.DataFrame(results)
    
    # Export to CSV
    try:
        import os
        script_dir = os.path.dirname(os.path.abspath(__file__))
        results_dir = script_dir
        os.makedirs(results_dir, exist_ok=True)
        csv_path = os.path.join(results_dir, "data_FigF_Wetting_Ranked_regen.csv")
        # Sort by Gamma_SL for ranking
        df_sorted = df.sort_values(by='Gamma_SL', ascending=True)
        df_sorted.to_csv(csv_path, index=False)
        print(f"Saved Wetting Analysis Data to {csv_path}")
    except Exception as e:
        print(f"Could not save CSV: {e}")
        
    return df

def plot_wetting(df):
    try:
        import os
        script_dir = os.path.dirname(os.path.abspath(__file__))
        results_dir = script_dir
        os.makedirs(results_dir, exist_ok=True)
    except:
        results_dir = os.path.dirname(os.path.abspath(__file__))

    plt.figure(figsize=(11, 9))
    sns.set_style("whitegrid")
    
    x = df['Delta_H_mix']
    y = df['Gamma_SL']
    
    colors = ['red' if el == 'Pt' else 'blue' for el in df['Host']]
    sizes = [150 if el == 'Pt' else 60 for el in df['Host']]
    
    plt.scatter(x, y, c=colors, s=sizes, alpha=0.7, edgecolors='k')
    
    # --- Refined Labeling ---
    labels_to_show = [
        'Pt', 'Pd', 'Ni', 'Au', 
        'W', 'Ta', 'Re', 'Mo', 
        'La', 'Ce', 'Co',       
        'Ru', 'Os', 'Ir', 'Rh'  
    ]
    
    for i, row in df.iterrows():
        label = row['Host']
        if label not in labels_to_show: continue
        
        x_pos = row['Delta_H_mix']
        y_pos = row['Gamma_SL']
        
        font_size = 9 
        weight = 'normal'
        dx = 0
        dy = 0
        ha = 'center'
        va = 'center'
        
        if label == 'Pt':
            font_size = 12
            weight = 'bold'
            dy = -0.2 # Below
            
        # 1. Ru / Os / Ir / Rh (Horizontal Spread)
        elif label == 'Ru':
            ha = 'right'
            dx = -0.12
        elif label == 'Os':
            ha = 'right'
            dx = -0.12
        elif label == 'Rh':
            ha = 'left'
            dx = 0.12
        elif label == 'Ir':
            ha = 'left'
            dx = 0.12
        elif label == 'Re':
            ha = 'left'
            dx = 0.12

        # 2. La / Ce (Merge)
        elif label == 'La':
            label = 'La/Ce' # Merge label
            ha = 'left'
            dx = 0.1
            dy = 0.05
        elif label == 'Ce':
            continue # Skip Ce, already covered by La/Ce
            
        # 3. Co (Near La/Ce)
        elif label == 'Co':
            ha = 'right'
            dx = -0.1
            dy = -0.05
            
        # 4. W / Ta / Mo
        elif label == 'W':
            ha = 'left'
            dx = 0.1
        elif label == 'Ta':
            ha = 'left'
            dx = 0.1
        elif label == 'Mo':
            ha = 'right'
            dx = -0.1
            
        # 5. Ni / Pd / Au
        elif label == 'Ni':
            ha = 'right'
            dx = -0.15
        elif label == 'Pd':
            ha = 'right'
            dx = -0.15
        elif label == 'Au':
            ha = 'left'
            dx = 0.15

        plt.text(x_pos + dx, y_pos + dy, label, fontsize=font_size, fontweight=weight, ha=ha, va=va)

    plt.axhline(0, color='green', linestyle='--', linewidth=2, label='Perfect Wetting Limit')
    plt.axvline(0, color='black', linestyle='-', linewidth=1)
    
    x_max = max(abs(x.min()), abs(x.max())) * 1.15
    plt.xlim(-x_max, x_max)
    plt.ylim(bottom=min(-0.5, y.min() * 1.1), top=max(1.5, y.max()*1.1))

    plt.title('Pillar 2: The Enabler - Wetting Analysis (Sparse & Merged Labels)', fontsize=14, fontweight='bold')
    plt.xlabel(r'Mixing enthalpy, $\Delta H_{mix}$ (kJ mol$^{-1}$)', fontsize=12, fontweight='bold')
    plt.ylabel(r'Solid-liquid interfacial free energy, $\gamma_{SL}$ (J m$^{-2}$)', fontsize=12, fontweight='bold')
    plt.legend()
    
    plt.text(-x_max*0.8, -0.2, "Quadrant III\nSpontaneous Wetting\n(Pt Zone)", color='green', 
             fontsize=10, fontweight='bold', bbox=dict(facecolor='white', alpha=0.8, edgecolor='green'))
    
    plt.text(x_max*0.5, 1.0, "Quadrant I\nNon-Wetting\n(Repulsive)", color='red', 
             fontsize=10, fontweight='bold', bbox=dict(facecolor='white', alpha=0.8, edgecolor='red'))
    
    plt.tight_layout()
    save_path = os.path.join(results_dir, "preview_FigF_Wetting_regen.png")
    plt.savefig(save_path, dpi=300)
    print(f"Saved Wetting Plot to {save_path}")

if __name__ == "__main__":
    df = calc_wetting_analysis()
    plot_wetting(df)
