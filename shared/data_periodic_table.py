
# Data Source: Miedema Parameters & Atomic Properties (Hardcoded for Reliability)
# Elements: 4th, 5th, 6th Period Transition Metals + Post-Transition Metals
# Format: {Symbol: {Phi (V), n_ws (density units), V_m (cm^3/mol), Tm (K), r (pm), gamma_L (J/m2), Chi (Pauling)}}

periodic_table_data = {
    # --- Period 4 ---
    'Sc': {'Phi': 3.25, 'n_ws': 1.27, 'V': 15.0, 'Tm': 1814, 'r': 162, 'gamma_L': 1.00, 'Chi': 1.36},
    'Ti': {'Phi': 3.65, 'n_ws': 1.47, 'V': 10.6, 'Tm': 1941, 'r': 147, 'gamma_L': 1.65, 'Chi': 1.54},
    'V':  {'Phi': 4.25, 'n_ws': 1.64, 'V': 8.36, 'Tm': 2183, 'r': 134, 'gamma_L': 1.85, 'Chi': 1.63},
    'Cr': {'Phi': 4.65, 'n_ws': 1.73, 'V': 7.23, 'Tm': 2180, 'r': 128, 'gamma_L': 1.70, 'Chi': 1.66},
    'Mn': {'Phi': 4.45, 'n_ws': 1.61, 'V': 7.39, 'Tm': 1519, 'r': 127, 'gamma_L': 1.10, 'Chi': 1.55},
    'Fe': {'Phi': 4.93, 'n_ws': 1.77, 'V': 7.10, 'Tm': 1811, 'r': 126, 'gamma_L': 1.87, 'Chi': 1.83},
    'Co': {'Phi': 5.10, 'n_ws': 1.75, 'V': 6.70, 'Tm': 1768, 'r': 125, 'gamma_L': 1.88, 'Chi': 1.88},
    'Ni': {'Phi': 5.20, 'n_ws': 1.75, 'V': 6.60, 'Tm': 1728, 'r': 124, 'gamma_L': 1.77, 'Chi': 1.91},
    'Cu': {'Phi': 4.45, 'n_ws': 1.47, 'V': 7.12, 'Tm': 1358, 'r': 128, 'gamma_L': 1.30, 'Chi': 1.90},
    'Zn': {'Phi': 4.10, 'n_ws': 1.32, 'V': 9.17, 'Tm': 693,  'r': 134, 'gamma_L': 0.78, 'Chi': 1.65}, # Target
    'Ga': {'Phi': 4.10, 'n_ws': 1.34, 'V': 11.8, 'Tm': 303,  'r': 135, 'gamma_L': 0.72, 'Chi': 1.81}, # Target
    'Ge': {'Phi': 4.55, 'n_ws': 1.37, 'V': 13.6, 'Tm': 1211, 'r': 122, 'gamma_L': 0.60, 'Chi': 2.01},

    # --- Period 5 ---
    'Y':  {'Phi': 3.20, 'n_ws': 1.11, 'V': 19.9, 'Tm': 1799, 'r': 180, 'gamma_L': 0.90, 'Chi': 1.22},
    'Zr': {'Phi': 3.40, 'n_ws': 1.39, 'V': 14.0, 'Tm': 2128, 'r': 160, 'gamma_L': 1.46, 'Chi': 1.33},
    'Nb': {'Phi': 4.00, 'n_ws': 1.62, 'V': 10.8, 'Tm': 2750, 'r': 146, 'gamma_L': 1.90, 'Chi': 1.6},
    'Mo': {'Phi': 4.65, 'n_ws': 1.77, 'V': 9.40, 'Tm': 2896, 'r': 139, 'gamma_L': 2.25, 'Chi': 2.16},
    'Tc': {'Phi': 5.30, 'n_ws': 1.90, 'V': 8.50, 'Tm': 2430, 'r': 136, 'gamma_L': 2.00, 'Chi': 1.9}, 
    'Ru': {'Phi': 5.40, 'n_ws': 1.83, 'V': 8.30, 'Tm': 2607, 'r': 134, 'gamma_L': 2.20, 'Chi': 2.2},
    'Rh': {'Phi': 5.40, 'n_ws': 1.76, 'V': 8.30, 'Tm': 2237, 'r': 134, 'gamma_L': 2.00, 'Chi': 2.28},
    'Pd': {'Phi': 5.45, 'n_ws': 1.67, 'V': 8.90, 'Tm': 1828, 'r': 137, 'gamma_L': 1.50, 'Chi': 2.20},
    'Ag': {'Phi': 4.35, 'n_ws': 1.39, 'V': 10.3, 'Tm': 1235, 'r': 144, 'gamma_L': 0.92, 'Chi': 1.93},
    'Cd': {'Phi': 4.05, 'n_ws': 1.16, 'V': 13.1, 'Tm': 594,  'r': 151, 'gamma_L': 0.60, 'Chi': 1.69},
    'In': {'Phi': 3.90, 'n_ws': 1.17, 'V': 15.7, 'Tm': 430,  'r': 167, 'gamma_L': 0.56, 'Chi': 1.78}, # Target
    'Sn': {'Phi': 4.15, 'n_ws': 1.25, 'V': 16.3, 'Tm': 505,  'r': 140, 'gamma_L': 0.56, 'Chi': 1.96}, # Target
    'Sb': {'Phi': 4.40, 'n_ws': 1.26, 'V': 18.2, 'Tm': 904,  'r': 140, 'gamma_L': 0.38, 'Chi': 2.05},

    # --- Period 6 ---
    'Hf': {'Phi': 3.55, 'n_ws': 1.43, 'V': 13.6, 'Tm': 2506, 'r': 159, 'gamma_L': 1.60, 'Chi': 1.3},
    'Ta': {'Phi': 4.05, 'n_ws': 1.63, 'V': 10.9, 'Tm': 3290, 'r': 146, 'gamma_L': 2.10, 'Chi': 1.5},
    'W':  {'Phi': 4.80, 'n_ws': 1.84, 'V': 9.60, 'Tm': 3695, 'r': 139, 'gamma_L': 2.50, 'Chi': 2.36},
    'Re': {'Phi': 5.40, 'n_ws': 1.96, 'V': 8.85, 'Tm': 3459, 'r': 137, 'gamma_L': 2.70, 'Chi': 1.9},
    'Os': {'Phi': 5.40, 'n_ws': 1.98, 'V': 8.43, 'Tm': 3306, 'r': 135, 'gamma_L': 2.50, 'Chi': 2.2},
    'Ir': {'Phi': 5.55, 'n_ws': 1.83, 'V': 8.54, 'Tm': 2719, 'r': 136, 'gamma_L': 2.20, 'Chi': 2.20},
    'Pt': {'Phi': 5.65, 'n_ws': 1.78, 'V': 9.10, 'Tm': 2041, 'r': 139, 'gamma_L': 1.80, 'Chi': 2.28}, # Host
    'Au': {'Phi': 5.15, 'n_ws': 1.57, 'V': 10.2, 'Tm': 1337, 'r': 144, 'gamma_L': 1.14, 'Chi': 2.54},
    'Hg': {'Phi': 4.20, 'n_ws': 1.13, 'V': 14.8, 'Tm': 234,  'r': 151, 'gamma_L': 0.48, 'Chi': 2.00},
    'Tl': {'Phi': 3.90, 'n_ws': 1.12, 'V': 17.2, 'Tm': 577,  'r': 170, 'gamma_L': 0.46, 'Chi': 1.62},
    'Pb': {'Phi': 4.10, 'n_ws': 1.15, 'V': 18.3, 'Tm': 600,  'r': 175, 'gamma_L': 0.46, 'Chi': 2.33},
    'Bi': {'Phi': 4.15, 'n_ws': 1.05, 'V': 21.3, 'Tm': 544,  'r': 150, 'gamma_L': 0.37, 'Chi': 2.02},
    
    # --- Lanthanides (Rare Earths - Key Controls) ---
    'La': {'Phi': 3.05, 'n_ws': 1.09, 'V': 22.5, 'Tm': 1193, 'r': 187, 'gamma_L': 0.72, 'Chi': 1.10},
    'Ce': {'Phi': 3.05, 'n_ws': 1.13, 'V': 21.0, 'Tm': 1068, 'r': 181, 'gamma_L': 0.72, 'Chi': 1.12},

    # --- Others for Reference (Al, Mg) ---
    'Al': {'Phi': 4.20, 'n_ws': 1.47, 'V': 10.0, 'Tm': 933,  'r': 143, 'gamma_L': 0.91, 'Chi': 1.61},
    'Mg': {'Phi': 3.45, 'n_ws': 1.17, 'V': 14.0, 'Tm': 923,  'r': 160, 'gamma_L': 0.56, 'Chi': 1.31}
}
