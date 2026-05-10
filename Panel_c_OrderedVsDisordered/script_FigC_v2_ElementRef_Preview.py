"""
Panel c v2 — element-referenced preview plot.

Input : data_FigC_Long.csv (already in element reference)
Output: 20260430_FigC_v2_ElementRef_Preview.png
        Origin-ready CSV is unchanged (already correct).

This script generates a publication-quality preview using the *same*
element-referenced ΔH_f values that appear in:
  - data_FigC_Summary.csv
  - SI Note S1
  - Manuscript Version 3 main text

Key numbers (matching the source CSV, NOT the deprecated internal Pt3Ga
reference):
  ΔH_f^HEI  = -30.008 kJ mol-1 atom-1
  ΔH_f^HEA  = -13.964 +/- 2.262 kJ mol-1 atom-1 (N=30)
  Ordering gap = 16.044 kJ mol-1 atom-1
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

LONG = SCRIPT_DIR / "data_FigC_Long.csv"
SUMMARY = SCRIPT_DIR / "data_FigC_Summary.csv"


def main():
    df = pd.read_csv(LONG)
    summary = pd.read_csv(SUMMARY).set_index("Metric")["Value"]

    # Pull verified numbers from the summary CSV — single source of truth
    e_ord = float(summary["Ordered_ElementRef_Hf_kJ_mol"])
    e_dis_mean = float(summary["Disordered_Mean_ElementRef_Hf_kJ_mol"])
    e_dis_std = float(summary["Disordered_Std_ElementRef_Hf_kJ_mol"])
    gap = float(summary["Ordering_Gap_kJ_mol"])

    e_dis = df[df["Type"] == "Disordered_Random_Equimolar"]["ZeroK_ElementRef_Hf_kJ_mol"].to_numpy()

    fig, ax = plt.subplots(figsize=(4.0, 5.5))
    rng = np.random.default_rng(42)

    # HEI: single point on left
    x_ord = 0.0
    ax.scatter([x_ord], [e_ord], s=180, c="#C73E3A", marker="o", edgecolor="black",
               linewidth=1.2, label="HE-IMC (Ordered L1$_2$)", zorder=5)
    ax.annotate(f"{e_ord:.2f}", xy=(x_ord, e_ord), xytext=(-30, 0),
                textcoords="offset points", fontsize=11, color="#C73E3A",
                ha="right", va="center", fontweight="bold")

    # HEA: 30 jittered points on right
    x_dis = 1.0 + rng.uniform(-0.18, 0.18, size=len(e_dis))
    ax.scatter(x_dis, e_dis, s=40, c="#7B6FA8", alpha=0.65,
               edgecolor="#3B345E", linewidth=0.4, label="HEA (random, $N=30$)", zorder=4)

    # HEA box / mean band
    ax.hlines(e_dis_mean, 0.78, 1.22, color="#3B345E", linewidth=2.0, zorder=6)
    ax.fill_between([0.78, 1.22],
                    e_dis_mean - e_dis_std, e_dis_mean + e_dis_std,
                    color="#7B6FA8", alpha=0.18, zorder=2)
    ax.annotate(f"{e_dis_mean:.2f} ± {e_dis_std:.2f}",
                xy=(1.22, e_dis_mean), xytext=(8, 0),
                textcoords="offset points", fontsize=10.5, color="#3B345E",
                ha="left", va="center", fontweight="bold")

    # Ordering gap arrow
    ax.annotate("",
                xy=(0.42, e_ord), xytext=(0.42, e_dis_mean),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax.text(0.42, (e_ord + e_dis_mean) / 2,
            f"  ΔH$_{{ordering}}$\n  = {gap:.2f} kJ mol$^{{-1}}$",
            fontsize=10, ha="left", va="center")

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["HE-IMC\n(Ordered L1$_2$)", "HEA\n(Disordered FCC)"], fontsize=10)
    ax.set_xlim(-0.45, 1.65)
    ax.set_ylabel("Formation enthalpy, $\\Delta H_f$ (kJ mol$^{-1}$ atom$^{-1}$)",
                  fontsize=11)
    ax.set_title("Panel c v2 — element-referenced", fontsize=11)

    # y-axis range: leave headroom around values
    y_lo = min(e_ord, e_dis.min()) - 3
    y_hi = max(e_ord, e_dis.max()) + 3
    ax.set_ylim(y_lo, y_hi)
    ax.invert_yaxis()  # convention: more negative at top

    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend(loc="lower right", fontsize=9, frameon=True)

    out_png = RESULTS_DIR / "preview_FigC_v2_ElementRef_regen.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Saved Panel c v2 preview: {out_png}")
    print(f"Values used: HEI={e_ord:.4f}, HEA mean={e_dis_mean:.4f} +/- {e_dis_std:.4f},"
          f" gap={gap:.4f} kJ mol-1 atom-1")


if __name__ == "__main__":
    main()
