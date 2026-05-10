"""
script_D_panel_a_plot_v2.py
===========================
Render the v2 quantitative Panel a schematic with:
 - Three tiers (liquid top / Pt middle / HEI bottom)
 - Five element-resolved trajectories with DIFFERENT heights
 - Driving-force arrow + Delta_G_rxn annotation

v2 key visual changes:
 - Liquid tier: 4 elements at different heights (+0.76 ~ +1.83 kJ/mol atom)
 - HEI tier: 5 elements at different depths (-9.86 ~ -15.61 kJ/mol atom)
 - Driving force labels per element

Inputs: panel_a_v2_tier_summary.csv, delta_mu_v2_0K.csv, delta_G_rxn_v2_summary.csv
Outputs: panel_a_v2_schematic.png, .pdf

Author: 2026-04-23 v2
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch

ELEMENT_COLORS: dict[str, str] = {
    "Ga": "#E07A5F",
    "In": "#3D5A80",
    "Sn": "#81B29A",
    "Zn": "#F2CC8F",
    "Pt": "#2D2D2D",
}

# Horizontal fan positions
X_LIQ_FAN = {"Ga": 0.08, "In": 0.15, "Sn": 0.22, "Zn": 0.29}
X_PT = 0.50
X_HEI_FAN = {"Pt": 0.71, "Ga": 0.78, "In": 0.85, "Sn": 0.92, "Zn": 0.99}


def load(out: Path):
    tier = pd.read_csv(out / "panel_a_v2_tier_summary.csv")
    dmu = pd.read_csv(out / "delta_mu_v2_0K.csv")
    dG = pd.read_csv(out / "delta_G_rxn_v2_summary.csv").iloc[0]
    return tier, dmu, dG


def plot(tier: pd.DataFrame, dmu: pd.DataFrame, dG: pd.Series, out_path: Path):
    fig, ax = plt.subplots(figsize=(12.0, 7.5))

    hei = tier[tier["tier"] == "HEI_bottom"]
    liq = tier[tier["tier"] == "liquid_top"]
    hei_min = hei["mu_per_atom_kJmol"].min()
    hei_max = hei["mu_per_atom_kJmol"].max()
    liq_max = liq["mu_per_atom_kJmol"].max()

    # ?? Shaded tier bands ??
    ax.axhspan(liq["mu_per_atom_kJmol"].min() - 0.3, liq_max + 0.3,
               xmin=0.01, xmax=0.38, facecolor="#EAE6D9", alpha=0.5, zorder=1)
    ax.axhspan(-0.15, 0.15, xmin=0.43, xmax=0.57,
               facecolor="#DDDDDD", alpha=0.5, zorder=1)
    ax.axhspan(hei_min - 0.5, hei_max + 0.5,
               xmin=0.62, xmax=1.05, facecolor="#CFE0CF", alpha=0.5, zorder=1)

    # ?? Element trajectories ??
    beta_els = ["Ga", "In", "Sn", "Zn"]
    for el in beta_els:
        color = ELEMENT_COLORS[el]
        mu_s = liq[liq["element"] == el]["mu_per_atom_kJmol"].values[0]
        mu_e = hei[hei["element"] == el]["mu_per_atom_kJmol"].values[0]
        xs, xe = X_LIQ_FAN[el], X_HEI_FAN[el]
        ax.plot([xs, xe], [mu_s, mu_e], color=color, linestyle=(0, (6, 3)),
                linewidth=2.0, alpha=0.85, zorder=3)
        ax.scatter([xs], [mu_s], s=130, color=color, edgecolor="white",
                   linewidths=1.4, zorder=6)
        ax.scatter([xe], [mu_e], s=150, color=color, edgecolor="white",
                   linewidths=1.4, zorder=6)
        ax.annotate(el, xy=(xs, mu_s), xytext=(xs, mu_s + 0.7),
                    fontsize=12, fontweight="bold", color=color,
                    ha="center", va="bottom")
        dy_end = 0.7 if el in ("Ga", "Sn") else -0.7
        va_end = "bottom" if el in ("Ga", "Sn") else "top"
        ax.annotate(el, xy=(xe, mu_e), xytext=(xe, mu_e + dy_end),
                    fontsize=11, fontweight="bold", color=color,
                    ha="center", va=va_end)

    # Pt trajectory (middle -> HEI)
    mu_Pt_s = 0.0
    mu_Pt_e = hei[hei["element"] == "Pt"]["mu_per_atom_kJmol"].values[0]
    ax.plot([X_PT, X_HEI_FAN["Pt"]], [mu_Pt_s, mu_Pt_e],
            color=ELEMENT_COLORS["Pt"], linestyle=(0, (6, 3)),
            linewidth=2.0, alpha=0.85, zorder=3)
    ax.scatter([X_PT], [mu_Pt_s], s=160, color="#2D2D2D", edgecolor="white",
               linewidths=1.4, zorder=6)
    ax.scatter([X_HEI_FAN["Pt"]], [mu_Pt_e], s=160, color="#2D2D2D",
               edgecolor="white", linewidths=1.4, zorder=6)
    ax.annotate("Pt", xy=(X_PT, mu_Pt_s), xytext=(X_PT, mu_Pt_s - 1.0),
                fontsize=12.5, fontweight="bold", color="#2D2D2D",
                ha="center", va="top")
    ax.annotate("Pt", xy=(X_HEI_FAN["Pt"], mu_Pt_e),
                xytext=(X_HEI_FAN["Pt"], mu_Pt_e - 0.7),
                fontsize=11, fontweight="bold", color="#2D2D2D",
                ha="center", va="top")

    # ?? Tier titles ??
    ax.text(0.185, liq_max + 1.8,
            "Liquid cocktail\n"
            r"(Ga$_{0.65}$In$_{0.20}$Sn$_{0.10}$Zn$_{0.05}$)",
            ha="center", va="center", fontsize=11, fontweight="semibold")
    ax.text(X_PT, liq_max + 1.8,
            "Pt substrate\n(s, fcc, SER = 0)",
            ha="center", va="center", fontsize=11, fontweight="semibold")
    ax.text(0.86, liq_max + 1.8,
            "HEI product\n"
            r"Pt$_3$(Ga,In,Sn,Zn) L1$_2$",
            ha="center", va="center", fontsize=11, fontweight="semibold")

    # ?? Driving-force arrow ??
    dG_per_atom = float(dG["delta_G_rxn_per_atom_kJmol"])
    dG_per_fu = float(dG["delta_G_rxn_per_fu_kJmol"])
    arrow = FancyArrowPatch(
        posA=(1.06, liq_max), posB=(1.06, hei["mu_per_atom_kJmol"].mean()),
        arrowstyle="-|>", mutation_scale=24,
        color="#C8102E", linewidth=2.4, zorder=7)
    ax.add_patch(arrow)
    ax.text(1.08, (liq_max + hei["mu_per_atom_kJmol"].mean()) / 2,
            f"$\\Delta G_\\mathrm{{rxn}}$\n"
            f"= {dG_per_atom:.1f}\n  kJ/mol atom\n"
            f"= {dG_per_fu:.1f}\n  kJ/mol f.u.",
            ha="left", va="center", fontsize=9.5, color="#C8102E",
            fontweight="semibold")

    # ?? Driving-force summary box ??
    lines = [r"$F_i = \mu_i^\mathrm{start} - \mu_i^\mathrm{HEI}$"
             + "\n  (kJ/mol atom)\n"]
    for _, row in dmu.iterrows():
        lines.append(
            f"  {row['element']:<2s}  {row['driving_force_per_atom_kJmol']:+.2f}")
    box_txt = "\n".join(lines)
    ax.text(0.02, hei_min - 3.8, box_txt, fontsize=9.5,
            family="DejaVu Sans Mono",
            bbox=dict(facecolor="#F8F8F5", edgecolor="#BBB",
                      boxstyle="round,pad=0.5"),
            va="bottom", ha="left")

    # ?? Axis cosmetics ??
    ax.set_xlim(-0.03, 1.28)
    ax.set_ylim(hei_min - 4.5, liq_max + 3.8)
    ax.set_ylabel("Chemical potential per atom of f.u. (kJ/mol)", fontsize=11.5)
    ax.set_xticks([0.185, X_PT, 0.86])
    ax.set_xticklabels(["start\n(liquid reservoir)",
                        "start\n(Pt substrate)",
                        "end\n(HEI sink)"], fontsize=10.5)
    ax.axhline(0, color="#999", linewidth=0.7, linestyle=":", zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_title(
        r"Panel a  chemical-potential landscape (SER reference, 0 K enthalpic)"
        "\n"
        r"liquid reservoir $\to$ Pt substrate $\to$ HEI sink",
        fontsize=13, fontweight="bold", pad=12)

    plt.tight_layout()
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    out = Path(__file__).resolve().parent.parent / "outputs"
    tier, dmu, dG = load(out)
    out_path = out / "panel_a_v2_schematic"
    plot(tier, dmu, dG, out_path)
    meta = {
        "script": "script_D_panel_a_plot_v2.py",
        "inputs": ["panel_a_v2_tier_summary.csv", "delta_mu_v2_0K.csv",
                    "delta_G_rxn_v2_summary.csv"],
        "outputs": ["panel_a_v2_schematic.png", "panel_a_v2_schematic.pdf"],
    }
    (out / "script_D_v2_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Panel a v2 rendered to:\n  {out_path.with_suffix('.png')}\n  {out_path.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
