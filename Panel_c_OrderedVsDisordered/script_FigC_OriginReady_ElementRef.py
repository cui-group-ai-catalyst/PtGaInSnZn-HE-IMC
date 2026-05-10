"""
script_FigC_OriginReady_ElementRef.py
=====================================
Panel c — element-referenced formation-enthalpy preview.

Formation enthalpy is reported on the stable-element-reference (SER) axis:
    ΔH_f = E_alloy - Σ x_i · E_i^elem
where each E_i^elem is the UMA-s-1p1 single-point energy of the stable
elemental reference structure of Pt, Ga, In, Sn, Zn (see SI methods).
Ordered (L1₂) and disordered (random, 30 configs) supercells are at
identical overall composition, so changing the reference state shifts both
curves uniformly; the ordering gap (≈ 16.04 kJ mol⁻¹ atom⁻¹) is therefore
reference-state invariant.

Inputs
------
- data_FigC_Long.csv
  (1 ordered + 30 disordered rows; column ZeroK_ElementRef_Hf_kJ_mol holds
  the SER-referenced ΔH_f.)

Outputs (written to the panel directory)
----------------------------------------
data_FigC_OriginReady_Wide_regen.csv
    Wide format: Ordered_ElementRef_Hf_kJ_mol, Disordered_ElementRef_Hf_kJ_mol.
data_FigC_OriginReady_Long_regen.csv
    Long format: Type, Config_ID, Seed, ZeroK_ElementRef_Hf_kJ_mol.
data_FigC_Summary_regen.csv
    Summary statistics: ordered / disordered_mean / disordered_std / gap.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SRC = SCRIPT_DIR / "data_FigC_Long.csv"


def main() -> None:
    if not SRC.exists():
        raise FileNotFoundError(f"Missing source CSV: {SRC}")

    df = pd.read_csv(SRC)

    ordered = df[df["Type"] == "Ordered_L12_Equimolar"].copy()
    disord = df[df["Type"] == "Disordered_Random_Equimolar"].copy()

    if len(ordered) != 1:
        raise ValueError(f"Expected 1 ordered row, got {len(ordered)}")
    if len(disord) != 30:
        raise ValueError(f"Expected 30 disordered rows, got {len(disord)}")

    e_ord = float(ordered["ZeroK_ElementRef_Hf_kJ_mol"].iloc[0])
    e_dis_arr = disord["ZeroK_ElementRef_Hf_kJ_mol"].to_numpy(dtype=float)
    e_dis_mean = float(np.mean(e_dis_arr))
    e_dis_std = float(np.std(e_dis_arr, ddof=1))  # sample std
    gap = e_dis_mean - e_ord

    # --- Wide (同旧 Origin_Ready 格式，仅保留 ElementRef 列) ---
    n_rows = max(1, len(e_dis_arr))
    wide = pd.DataFrame({
        "Ordered_ElementRef_Hf_kJ_mol": [e_ord] + [np.nan] * (n_rows - 1),
        "Disordered_ElementRef_Hf_kJ_mol": list(e_dis_arr),
    })
    wide_path = RESULTS_DIR / "data_FigC_OriginReady_Wide_regen.csv"
    wide.to_csv(wide_path, index=False, float_format="%.6f")
    print(f"[Panel c] Origin-ready wide -> {wide_path}")

    # --- Long（带 Type/Seed，便于 Origin 分组） ---
    long_rows = [{
        "Type": "Ordered_L12_Equimolar",
        "Config_ID": 0,
        "Seed": 0,
        "ZeroK_ElementRef_Hf_kJ_mol": e_ord,
    }]
    for _, r in disord.iterrows():
        long_rows.append({
            "Type": "Disordered_Random_Equimolar",
            "Config_ID": int(r["Config_ID"]),
            "Seed": int(r["Seed"]),
            "ZeroK_ElementRef_Hf_kJ_mol": float(r["ZeroK_ElementRef_Hf_kJ_mol"]),
        })
    long_df = pd.DataFrame(long_rows)
    long_path = RESULTS_DIR / "data_FigC_OriginReady_Long_regen.csv"
    long_df.to_csv(long_path, index=False, float_format="%.6f")
    print(f"[Panel c] Origin-ready long -> {long_path}")

    # --- Summary ---
    summary = pd.DataFrame([
        {"Metric": "Ordered_ElementRef_Hf_kJ_mol",       "Value": round(e_ord, 4)},
        {"Metric": "Disordered_Mean_ElementRef_Hf_kJ_mol", "Value": round(e_dis_mean, 4)},
        {"Metric": "Disordered_Std_ElementRef_Hf_kJ_mol",  "Value": round(e_dis_std, 4)},
        {"Metric": "Disordered_N",                        "Value": int(len(e_dis_arr))},
        {"Metric": "Ordering_Gap_kJ_mol",                  "Value": round(gap, 4)},
        {"Metric": "Ordering_Gap_sign_note",
         "Value": "positive = ordered is lower (more stable)"},
    ])
    sum_path = RESULTS_DIR / "data_FigC_Summary_regen.csv"
    summary.to_csv(sum_path, index=False)
    print(f"[Panel c] summary           -> {sum_path}")

    print("\n-- Panel c key numbers (element-reference axis) --")
    print(f"  Ordered             = {e_ord:.4f} kJ/mol")
    print(f"  Disordered mean     = {e_dis_mean:.4f} +/- {e_dis_std:.4f} kJ/mol (N=30)")
    print(f"  Ordering gap (dis-ord) = {gap:+.4f} kJ/mol")
    print("  Note: gap unchanged vs historical -4.356 display axis (same composition pair)")


if __name__ == "__main__":
    main()
