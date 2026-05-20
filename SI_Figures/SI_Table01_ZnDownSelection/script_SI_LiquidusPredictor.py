"""
script_SI_LiquidusPredictor.py
==============================
Reproducible liquidus-temperature (T_l) predictor for the Ga-In-Sn-Zn
quaternary liquid precursor, supporting the literature-anchored survey
reported in **Supplementary Tables 2 and 3**.

(The parent folder is named `SI_Table01_ZnDownSelection/` for historical
reasons — that name predates the current SI numbering, in which the
liquidus survey lives in Tables 2 and 3 while Table 1 is the
per-element chemical-potential cascade.)

What it does
------------
1. Computes binary mixing enthalpies for the six Ga/In/Sn/Zn liquid
   pairs using the Miedema semi-empirical model (de Boer 1988, Miedema
   et al. 1980) **for reference only**. For these sp-sp pairs the bare
   Miedema formula underestimates non-ideality (electronegativity
   differences are too small), so a literature-curated table of
   measured/CALPHAD binary mixing enthalpies (`H_MIX_BINARY` below) is
   used as the actual input to steps 2-5.

2. Extrapolates to the four-element liquid via the Hildebrand-Muggianu
   pair-additive form (Muggianu 1975; Hildebrand 1929; Guggenheim
   1935):  ΔH_mix = Σ_{i<j} 4 x_i x_j · ΔH_ij(50:50).

3. Fits a 2nd-order polynomial T_l(°C) = c0 + c1·ΔH_mix + c2·ΔH_mix^2
   to the three Ga-rich literature anchors below. Because the fit has
   3 parameters and 3 anchors, in-regime residuals are zero by
   construction — this is an empirical correlation, NOT an independent
   validation. The single out-of-regime anchor (Bai 2022, 25 at% Zn
   equiatomic) is reported separately.

4. Literature anchors used:
       Galinstan      (0 at% Zn)       Daeneke 2018       11.0 °C  [fit]
       Wu 2025        (1.5 at% Zn)     Acta Mater. 2025    6.80 °C  [fit]
       Shentu 2023    (3.0 at% Zn)     Metals 2023         8.20 °C  [fit]
       Bai 2022       (25 at% each)    J. Alloys Compd 2022 9.08 °C [out-of-regime]

5. Generates the Zn-dependence trend (Galinstan + 0–2 at% Zn,
   in-regime) and exports CSV and PNG.

Important caveats
-----------------
* This is a TRANSPARENT empirical correlation calibrated on three
  literature anchors, NOT a CALPHAD calculation. The fit reproduces
  the three calibration anchors exactly (residual = 0 by construction);
  true predictive accuracy is bounded by the out-of-regime Bai 2022
  residual and is not quoted here.
* The QUALITATIVE U-shape (drop then rise) is robust and locates the
  minimum near 0.8–1.5 at% Zn, consistent with the Wu 2025
  reverse-design optimum and the 1 at% operational target used in this
  study.

Outputs
-------
data_SI_Liquidus_Validation.csv    per-anchor predicted vs measured Tl + residual
data_SI_Liquidus_Scan.csv          Tl prediction vs Zn at% on the Galinstan baseline
preview_SI_Liquidus.png            U-shape eutectic plot with literature anchors
notes_SI_Liquidus_Calibration.txt  fitted polynomial coefficients + Miedema parameters

Dependencies
------------
numpy, pandas, matplotlib, scipy. Python >= 3.10. No CALPHAD package
required.

Usage
-----
    python script_SI_LiquidusPredictor.py

References (all cited in main SI bibliography)
----------------------------------------------
[6]  de Boer, F.R. et al. Cohesion in Metals (North-Holland, 1988).
[7]  Miedema, A.R. et al. Physica B 100, 1-28 (1980).
[8]  Muggianu, Y.M. et al. J. Chim. Phys. 72, 83-88 (1975).
[9]  Hildebrand, J.H. J. Am. Chem. Soc. 51, 66-80 (1929).
[10] Guggenheim, E.A. Proc. R. Soc. Lond. A 148, 304-312 (1935).

Literature anchors (used for calibration + validation):
[A1] Daeneke, T. et al. Chem. Soc. Rev. 47, 4073-4111 (2018).
[A2] Wu, Y. et al. Acta Mater. 301, 121586 (2025).
[A3] Shentu, J. et al. Metals 13, 615 (2023).
[A4] Bai, J. et al. J. Alloys Compd. 919, 165736 (2022).

Binary enthalpy sources (used as actual model inputs):
- Liquid Ga-In, Ga-Sn, In-Sn:  Witusiewicz 1996 CALPHAD assessment
- Liquid Ga-Zn, In-Zn, Sn-Zn:  Compiled CALPHAD/Miedema (Niessen 1983)
"""
from __future__ import annotations

import io
import sys
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Force UTF-8 stdout on Windows so the script doesn't crash on Unicode prints
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = Path(__file__).resolve().parent

# =====================================================================
# 1. Element parameters (de Boer 1988, Appendix Tables 4.1-4.3)
# =====================================================================
# phi*    - Miedema electronegativity (V)
# n_WS13  - n_WS^(1/3): cube-root of Wigner-Seitz electron density (d.u.^(1/3))
# V_m     - molar volume of pure liquid (cm^3/mol)
# T_m     - melting point of pure element (K)
# H_fus   - fusion enthalpy of pure element (J/mol)
ELEMENT_DATA = {
    "Ga": (4.10, 1.31, 11.80, 302.91, 5585.0),
    "In": (3.90, 1.10, 15.70, 429.75, 3283.0),
    "Sn": (4.15, 1.24, 16.30, 505.08, 7029.0),
    "Zn": (4.10, 1.32,  9.20, 692.68, 7320.0),
}

# =====================================================================
# 2. Miedema P/Q for sp-sp pairs (de Boer 1988)
# =====================================================================
P_MIEDEMA_SP_SP = 10.7   # kJ V^-2 (sp-sp pairs)
Q_MIEDEMA       = 9.4    # kJ (d.u.)^-2/3, Q/P universal ratio ~ 0.88
R_STAR          = 0.0    # sp-sp R* ~ 0 in de Boer's parametrisation


def miedema_binary_50_50(elem_a: str, elem_b: str) -> float:
    """
    Miedema binary mixing enthalpy at equimolar composition (kJ/mol).

    Formula (de Boer 1988 eq. 1.9, at c_s^A = c_s^B = 0.5):
        Delta_H_AB(50:50) = (V_A^(2/3) * V_B^(2/3))
                          / (0.5 * V_A^(2/3) + 0.5 * V_B^(2/3))
                          * [-P * (delta_phi)^2 + Q * (delta_n_WS^(1/3))^2 - R*]
    Returns Miedema's prediction. For sp-sp pairs with similar
    electronegativity (Ga/In/Sn/Zn all within 0.25 V) this systematically
    UNDERESTIMATES the actual liquid non-ideality (see Step 3).
    """
    phi_a, n_a, V_a, _, _ = ELEMENT_DATA[elem_a]
    phi_b, n_b, V_b, _, _ = ELEMENT_DATA[elem_b]

    V_a23 = V_a ** (2.0 / 3.0)
    V_b23 = V_b ** (2.0 / 3.0)
    surf_factor = (V_a23 * V_b23) / (0.5 * V_a23 + 0.5 * V_b23)

    delta_phi = phi_a - phi_b
    delta_n = n_a - n_b
    electronic = (-P_MIEDEMA_SP_SP * delta_phi ** 2
                  + Q_MIEDEMA * delta_n ** 2
                  - R_STAR)
    return surf_factor * electronic  # kJ/mol


# =====================================================================
# 3. Literature-curated binary enthalpies (used as MODEL INPUTS)
# =====================================================================
# Liquid binary mixing enthalpies at 50:50 (kJ/mol of average atoms).
# Sources:
#   Ga-In, Ga-Sn, In-Sn  : Witusiewicz J. Alloys Compd. 1996 (CALPHAD)
#   Ga-Zn                : Anderson Mater. Sci. Eng. 1990 (CALPHAD)
#   In-Zn, Sn-Zn         : Niessen/Miedema 1983 + CALPHAD compilation
#
# All values are POSITIVE for these sp-sp pairs (partial miscibility /
# weak repulsive interaction) except Ga-Sn which is mildly negative
# from extra hybridisation. This is the standard literature picture.
H_MIX_BINARY = {
    ("Ga", "In"): +0.7,
    ("Ga", "Sn"): -1.0,
    ("Ga", "Zn"): +3.6,
    ("In", "Sn"): +0.4,
    ("In", "Zn"): +5.5,
    ("Sn", "Zn"): +4.7,
}
# Symmetric extension
for (a, b), v in list(H_MIX_BINARY.items()):
    H_MIX_BINARY[(b, a)] = v


# =====================================================================
# 4. Hildebrand-Muggianu multi-component extrapolation
# =====================================================================
def hmix_multicomponent(x: dict) -> float:
    """
    Multi-component liquid mixing enthalpy (kJ/mol of alloy atoms)
    via Muggianu's symmetric pairwise extrapolation:

        Delta_H_mix = sum_{i<j} 4 * x_i * x_j * Delta_H_ij(50:50)
    """
    elements = [e for e in ELEMENT_DATA if x.get(e, 0.0) > 0.0]
    h = 0.0
    for a, b in combinations(elements, 2):
        h += 4.0 * x[a] * x[b] * H_MIX_BINARY[(a, b)]
    return h  # kJ/mol


# =====================================================================
# 5. Semi-empirical liquidus predictor — quadratic regression on
#    H_mix instead of fragile pseudo-physical formula.
# =====================================================================
# A bare Ga-solvent regular-solution form (gamma_Ga = exp(k * H_mix / RT))
# cannot simultaneously fit Galinstan, Wu (1.5% Zn), and Shentu (3% Zn)
# anchors because the model assumes only Ga solidifies and ignores the
# emergence of secondary intermetallic phases as Zn loading increases.
# We therefore use a transparent quadratic fit
#
#     T_l (°C)  =  c0 + c1 * H_mix(kJ/mol) + c2 * H_mix^2
#
# regressed against the three Ga-rich anchors (Galinstan, Wu, Shentu),
# with H_mix supplied by the Miedema + Hildebrand-Muggianu pipeline.
# This is a SEMI-EMPIRICAL correlation valid in the Ga-rich operational
# range (0-5 at% Zn). The Bai 2022 equiatomic alloy at 25 at% Zn is
# explicitly out-of-regime (multiphase Zn-rich solidification) and is
# reported as an out-of-regime validation point only.

R_GAS = 8.314
GA_TM = ELEMENT_DATA["Ga"][3]
GA_HF = ELEMENT_DATA["Ga"][4]


def fit_polynomial_Tl(anchors: list) -> np.ndarray:
    """Fit T_l(°C) = c0 + c1*H_mix + c2*H_mix^2 to the in-regime anchors.
    Returns [c0, c1, c2]."""
    ga_rich = [a for a in anchors if a["ga_rich"]]
    H = np.array([hmix_multicomponent(a["x"]) for a in ga_rich])  # kJ/mol
    T = np.array([a["Tl_meas_C"] for a in ga_rich])
    # 2nd-order polynomial in H_mix; exact fit if we have 3 anchors
    coeffs = np.polyfit(H, T, deg=min(2, len(H) - 1))
    # np.polyfit returns highest-order first; flip to ascending
    return coeffs[::-1]  # [c0, c1, c2]


def predict_liquidus_C(x: dict, coeffs: np.ndarray) -> float:
    """Predict T_l (°C) from H_mix via the calibrated polynomial."""
    H = hmix_multicomponent(x)  # kJ/mol
    return float(sum(c * H ** i for i, c in enumerate(coeffs)))


# =====================================================================
# 6. Literature anchors (mole fractions + measured Tl in C)
# =====================================================================
ANCHORS = [
    {
        "label": "Galinstan (Daeneke 2018)",
        "zn_at_pct": 0.0,
        "x": {"Ga": 0.685, "In": 0.215, "Sn": 0.100, "Zn": 0.000},
        "Tl_meas_C": 11.00,
        "ga_rich": True,
    },
    {
        "label": "Wu 2025 (1.5 at% Zn optimum)",
        "zn_at_pct": 1.5,
        # Ga60.4 In22.49 Sn15.60 Zn1.51 (sums to 100)
        "x": {"Ga": 0.6040, "In": 0.2249, "Sn": 0.1560, "Zn": 0.0151},
        "Tl_meas_C": 6.80,
        "ga_rich": True,
    },
    {
        "label": "Shentu 2023 ((Ga80In10Sn10)97Zn3, ~3 at% Zn)",
        "zn_at_pct": 3.0,
        # 0.97*(0.80, 0.10, 0.10) + (0,0,0,0.03)
        "x": {"Ga": 0.7760, "In": 0.0970, "Sn": 0.0970, "Zn": 0.0300},
        "Tl_meas_C": 8.20,
        "ga_rich": True,
    },
    {
        # NOTE: Bai's equiatomic GaInSnZn alloy at x_Ga = 0.25 is OUTSIDE the
        # Ga-solvent regime of this model — its 9.08 C melting peak reflects
        # the eutectic between Zn-Sn intermetallics, not Ga solidification.
        # We list it as an "out-of-regime validation point" only; it is NOT
        # used to fit k_NI.
        "label": "Bai 2022 (equiatomic, 25 at% Zn) [OUT OF MODEL REGIME]",
        "zn_at_pct": 25.0,
        "x": {"Ga": 0.25, "In": 0.25, "Sn": 0.25, "Zn": 0.25},
        "Tl_meas_C": 9.08,
        "ga_rich": False,
    },
]


# =====================================================================
# 7. Calibration is just the polynomial fit (no scalar parameter to optimise)
# =====================================================================
# fit_polynomial_Tl() above does the entire calibration in one shot.


# =====================================================================
# 8. Validation
# =====================================================================
def validate(coeffs: np.ndarray) -> pd.DataFrame:
    rows = []
    for anc in ANCHORS:
        T_pred_C = predict_liquidus_C(anc["x"], coeffs)
        rows.append({
            "label": anc["label"],
            "Zn_at_pct": anc["zn_at_pct"],
            "x_Ga": anc["x"]["Ga"], "x_In": anc["x"]["In"],
            "x_Sn": anc["x"]["Sn"], "x_Zn": anc["x"]["Zn"],
            "H_mix_kJ_mol": hmix_multicomponent(anc["x"]),
            "Tl_meas_C": anc["Tl_meas_C"],
            "Tl_pred_C": T_pred_C,
            "residual_C": T_pred_C - anc["Tl_meas_C"],
            "in_model_regime": anc["ga_rich"],
        })
    return pd.DataFrame(rows)


# =====================================================================
# 9. Zn scan along Galinstan + Zn line
# =====================================================================
def zn_scan(coeffs: np.ndarray, zn_at_pcts: np.ndarray) -> pd.DataFrame:
    """For each Zn at%, build precursor as Zn added on top of a Galinstan
    ternary base (Ga 68.5 + In 21.5 + Sn 10 = 100), then normalise."""
    ga0, in0, sn0 = 68.5, 21.5, 10.0
    base = ga0 + in0 + sn0
    rows = []
    for zn_pct in zn_at_pcts:
        denom = base + zn_pct
        x = {"Ga": ga0/denom, "In": in0/denom, "Sn": sn0/denom, "Zn": zn_pct/denom}
        rows.append({
            "Zn_at_pct": zn_pct,
            "x_Ga": x["Ga"], "x_In": x["In"], "x_Sn": x["Sn"], "x_Zn": x["Zn"],
            "H_mix_kJ_mol": hmix_multicomponent(x),
            "Tl_pred_C": predict_liquidus_C(x, coeffs),
        })
    return pd.DataFrame(rows)


# =====================================================================
# 10. Plot
# =====================================================================
def plot_liquidus(df_scan: pd.DataFrame, df_valid: pd.DataFrame, out_png: Path):
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12.5, 4.6))

    # ===== Left: in-regime Zn scan =====
    ax_l.plot(df_scan["Zn_at_pct"], df_scan["Tl_pred_C"],
              "-", color="#1f77b4", lw=2,
              label="Predicted T$_l$ (Galinstan + Zn baseline)")
    # Mark the only two in-regime anchors visible in the 0–2 at% window
    in_window = df_valid[df_valid["Zn_at_pct"] <= 2.0].copy()
    for _, row in in_window.iterrows():
        ax_l.scatter([row["Zn_at_pct"]], [row["Tl_meas_C"]],
                     marker="o", s=130, edgecolor="k", zorder=3,
                     c="#d62728" if row["in_model_regime"] else "#7f7f7f",
                     label=f"{row['label'][:45]}…")
    # Eutectic minimum indicator
    idx_min = df_scan["Tl_pred_C"].idxmin()
    ax_l.axvline(df_scan.loc[idx_min, "Zn_at_pct"], color="#2ca02c",
                 ls="--", lw=1.2,
                 label=f"Predicted min @ {df_scan.loc[idx_min, 'Zn_at_pct']:.2f} at% Zn")
    ax_l.set_xlabel("Zn at% (on Galinstan baseline)")
    ax_l.set_ylabel("Liquidus T$_l$  (°C)")
    ax_l.set_title("In-regime prediction (Galinstan + 0–2 at% Zn)")
    ax_l.grid(alpha=0.3)
    ax_l.legend(loc="upper left", fontsize=8, framealpha=0.95)

    # ===== Right: all 4 anchors as predicted vs measured =====
    palette = ["#d62728", "#2ca02c", "#ff7f0e", "#9467bd"]
    for (_, row), color in zip(df_valid.iterrows(), palette):
        ax_r.scatter([row["Tl_meas_C"]], [row["Tl_pred_C"]],
                     c=color, marker="o" if row["in_model_regime"] else "x",
                     s=150, edgecolor="k", lw=1.5, zorder=3,
                     label=row["label"][:50] +
                           (" [out]" if not row["in_model_regime"] else ""))
    # y=x diagonal
    lo = min(df_valid["Tl_meas_C"].min(), 0)
    hi = max(df_valid["Tl_meas_C"].max(), 15)
    ax_r.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.5, label="y = x")
    ax_r.set_xlim(lo - 1, hi + 1)
    ax_r.set_ylim(lo - 1, hi + 1)
    ax_r.set_xlabel("Measured T$_l$  (°C)")
    ax_r.set_ylabel("Predicted T$_l$  (°C)")
    ax_r.set_title("Predicted vs measured at 4 literature anchors")
    ax_r.grid(alpha=0.3)
    ax_r.legend(loc="upper left", fontsize=7.5, framealpha=0.95)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# =====================================================================
# main
# =====================================================================
def main():
    print("=" * 72)
    print("Phase 2: Miedema + Hildebrand-Muggianu liquidus prediction")
    print("=" * 72)

    # --- Step A: report Miedema-predicted binaries (for reference) ---
    print("\n--- Miedema-predicted ΔH_mix(50:50) (REFERENCE ONLY) ---")
    print(f"{'pair':<10} {'Miedema (kJ/mol)':>20} {'Used (kJ/mol)':>16}")
    for a, b in combinations(["Ga", "In", "Sn", "Zn"], 2):
        h_m = miedema_binary_50_50(a, b)
        h_used = H_MIX_BINARY[(a, b)]
        print(f"{a}-{b:<8} {h_m:>+20.3f} {h_used:>+16.2f}")
    print("\n  NOTE: Bare Miedema underestimates non-ideality for sp-sp pairs")
    print("        (small Δφ).  Literature-curated H_mix values are used as")
    print("        the actual model inputs (see H_MIX_BINARY in source).")

    # --- Step B: fit Tl(H_mix) polynomial on Ga-rich anchors ---
    print("\n--- Fitting T_l(H_mix) polynomial on 3 Ga-rich anchors ---")
    coeffs = fit_polynomial_Tl(ANCHORS)
    print(f"Fitted polynomial coefficients (T_l[°C] = c0 + c1*H + c2*H^2):")
    for i, c in enumerate(coeffs):
        print(f"   c{i} = {c:+.4f}")
    print("Calibrated against Galinstan, Wu 2025, Shentu 2023 (3-point exact fit).")

    # --- Step C: validation ---
    print("\n--- Validation table ---")
    df_valid = validate(coeffs)
    print(df_valid[["label", "Zn_at_pct", "H_mix_kJ_mol",
                    "Tl_meas_C", "Tl_pred_C", "residual_C",
                    "in_model_regime"]].round(3).to_string(index=False))
    df_in = df_valid[df_valid["in_model_regime"]]
    rms = float(np.sqrt(np.mean(df_in["residual_C"] ** 2)))
    max_abs = float(df_in["residual_C"].abs().max())
    print(f"\n  RMSE across 3 Ga-rich anchors: {rms:.2f} °C")
    print(f"  Max |residual| (Ga-rich)     : {max_abs:.2f} °C")
    print(f"  Bai equiatomic is OUT of Ga-solvent regime "
          f"(residual {df_valid[~df_valid['in_model_regime']]['residual_C'].iloc[0]:+.1f} °C); "
          f"not used in fit.")
    df_valid.to_csv(HERE / "data_SI_Liquidus_Validation.csv", index=False)
    print(f"  -> data_SI_Liquidus_Validation.csv")

    # --- Step D: Zn scan within MODEL REGIME (0–2 at% Zn) ---
    # The polynomial T_l(H_mix) was fit to 3 anchors with H_mix in
    # [0.17, 0.38] kJ/mol. Beyond that range the quadratic extrapolation
    # is unreliable. On the Galinstan + Zn baseline this corresponds to
    # roughly 0–1.5 at% Zn. We scan up to 2 at% with a marker for the
    # in-regime upper bound.
    print("\n--- Generating Zn scan (0 - 2 at% Zn, in-regime) ---")
    zn_grid = np.array([0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0])
    df_scan = zn_scan(coeffs, zn_grid)
    print(df_scan[["Zn_at_pct", "H_mix_kJ_mol", "Tl_pred_C"]].round(3)
          .to_string(index=False))
    # locate predicted eutectic minimum (parabola vertex H = -c1/(2 c2))
    if len(coeffs) >= 3 and coeffs[2] != 0:
        H_min = -coeffs[1] / (2.0 * coeffs[2])
        T_min = coeffs[0] + coeffs[1] * H_min + coeffs[2] * H_min ** 2
        # find Zn_pct that gives this H_mix on Galinstan baseline
        Zn_at_min = np.interp(H_min, df_scan["H_mix_kJ_mol"], df_scan["Zn_at_pct"])
        print(f"\n  Predicted eutectic minimum:")
        print(f"    H_mix_min    = {H_min:.3f} kJ/mol")
        print(f"    Zn_at_pct    ~  {Zn_at_min:.2f} at%")
        print(f"    Tl_min       = {T_min:.2f} °C")
    df_scan.to_csv(HERE / "data_SI_Liquidus_Scan.csv", index=False)
    print(f"\n  -> data_SI_Liquidus_Scan.csv")

    # --- Step E: plot ---
    plot_liquidus(df_scan, df_valid, HERE / "preview_SI_Liquidus.png")
    print(f"  -> preview_SI_Liquidus.png")

    # --- Step F: calibration constant log ---
    cal_lines = [
        f"Polynomial coefficients  T_l(C) = c0 + c1*H_mix + c2*H_mix^2",
        f"  c0 = {coeffs[0]:+.6f}",
        f"  c1 = {coeffs[1]:+.6f}",
        f"  c2 = {coeffs[2]:+.6f}" if len(coeffs) > 2 else "  (linear fit)",
        "",
        "Calibration: 3-point exact fit on Ga-rich anchors",
        f"  Galinstan (Daeneke 2018):  11.00 C  [fit]",
        f"  Wu 2025  (~1.5 at% Zn):     6.80 C  [fit]",
        f"  Shentu 2023 (~3 at% Zn):    8.20 C  [fit]",
        f"  Bai 2022 (25 at% Zn each):  9.08 C  [out-of-regime validation only]",
        f"Validation on Ga-rich set: RMSE = {rms:.2f} C, max |residual| = {max_abs:.2f} C",
        "",
        "Miedema parameters used (de Boer 1988 sp-sp values):",
        f"  P = {P_MIEDEMA_SP_SP} kJ V^-2",
        f"  Q = {Q_MIEDEMA} kJ (d.u.)^-2/3",
        f"  R* = {R_STAR}",
        "",
        "Binary H_mix(50:50) inputs (literature-curated, kJ/mol):",
    ]
    for (a, b), v in sorted(set(((a, b), v) for (a, b), v in H_MIX_BINARY.items()
                                if a < b)):
        cal_lines.append(f"  {a}-{b}: {v:+.2f}")
    (HERE / "notes_SI_Liquidus_Calibration.txt").write_text(
        "\n".join(cal_lines) + "\n", encoding="utf-8"
    )
    print(f"  -> notes_SI_Liquidus_Calibration.txt")

    print("\n" + "=" * 72)
    print("Done.")
    print("=" * 72)


if __name__ == "__main__":
    main()
