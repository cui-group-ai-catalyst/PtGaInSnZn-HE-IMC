# -*- coding: utf-8 -*-
"""
script_D_panel_a_plot_v3.py
===========================
Render the v3 quantitative Panel a schematic with:
 - Three tiers (liquid top / Pt middle / HEI bottom)
 - Five element-resolved trajectories at distinct heights
 - Driving-force arrow and dG_rxn annotation
 - CEF-fitted HEI data (UMA, Panel g consistent)

v3 changes vs v2:
 - HEI tier uses UMA-CEF values (more negative, more spread)
 - Liquid tier unchanged (SER + DH_fus + Miedema)
 - Layout rescaled for larger y-range

Inputs: panel_a_v3_tier_summary.csv, delta_mu_v3_0K.csv,
        delta_G_rxn_v3_summary.csv
Outputs: panel_a_v3_schematic.png, .pdf

Author: 2026-04-23 v3
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

X_LIQ_FAN = {"Ga": 0.08, "In": 0.15, "Sn": 0.22, "Zn": 0.29}
X_PT = 0.50
X_HEI_FAN = {"Pt": 0.71, "Ga": 0.78, "In": 0.85, "Sn": 0.92, "Zn": 0.99}


def load(out: Path):
    tier = pd.read_csv(out / "panel_a_v3_tier_summary.csv")
    dmu = pd.read_csv(out / "delta_mu_v3_0K.csv")
    dG = pd.read_csv(out / "delta_G_rxn_v3_summary.csv").iloc[0]
    return tier, dmu, dG


def plot(tier: pd.DataFrame, dmu: pd.DataFrame, dG: pd.Series, out_path: Path):
    fig, ax = plt.subplots(figsize=(12.0, 8.5))

    hei = tier[tier["tier"] == "HEI_bottom"]
    liq = tier[tier["tier"] == "liquid_top"]
    hei_min = hei["mu_per_atom_kJmol"].min()
    hei_max = hei["mu_per_atom_kJmol"].max()
    liq_max = liq["mu_per_atom_kJmol"].max()
    liq_min = liq["mu_per_atom_kJmol"].min()

    ax.axhspan(liq_min - 0.5, liq_max + 0.5,
               xmin=0.01, xmax=0.38, facecolor="#EAE6D9", alpha=0.5, zorder=1)
    ax.axhspan(-0.8, 0.8, xmin=0.43, xmax=0.57,
               facecolor="#DDDDDD", alpha=0.5, zorder=1)
    ax.axhspan(hei_min - 1.5, hei_max + 1.5,
               xmin=0.62, xmax=1.05, facecolor="#CFE0CF", alpha=0.5, zorder=1)

    beta_els = ["Ga", "In", "Sn", "Zn"]
    y_range = liq_max - hei_min
    label_offset = y_range * 0.03

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
        ax.annotate(el, xy=(xs, mu_s), xytext=(xs, mu_s + label_offset),
                    fontsize=12, fontweight="bold", color=color,
                    ha="center", va="bottom")
        dy_end = label_offset if el in ("Ga", "Sn") else -label_offset
        va_end = "bottom" if el in ("Ga", "Sn") else "top"
        ax.annotate(el, xy=(xe, mu_e), xytext=(xe, mu_e + dy_end),
                    fontsize=11, fontweight="bold", color=color,
                    ha="center", va=va_end)

    mu_Pt_s = 0.0
    mu_Pt_e = hei[hei["element"] == "Pt"]["mu_per_atom_kJmol"].values[0]
    ax.plot([X_PT, X_HEI_FAN["Pt"]], [mu_Pt_s, mu_Pt_e],
            color=ELEMENT_COLORS["Pt"], linestyle=(0, (6, 3)),
            linewidth=2.0, alpha=0.85, zorder=3)
    ax.scatter([X_PT], [mu_Pt_s], s=160, color="#2D2D2D", edgecolor="white",
               linewidths=1.4, zorder=6)
    ax.scatter([X_HEI_FAN["Pt"]], [mu_Pt_e], s=160, color="#2D2D2D",
               edgecolor="white", linewidths=1.4, zorder=6)
    ax.annotate("Pt", xy=(X_PT, mu_Pt_s), xytext=(X_PT, mu_Pt_s - label_offset),
                fontsize=12.5, fontweight="bold", color="#2D2D2D",
                ha="center", va="top")
    ax.annotate("Pt", xy=(X_HEI_FAN["Pt"], mu_Pt_e),
                xytext=(X_HEI_FAN["Pt"], mu_Pt_e - label_offset),
                fontsize=11, fontweight="bold", color="#2D2D2D",
                ha="center", va="top")

    title_y = liq_max + y_range * 0.12
    ax.text(0.185, title_y,
            "Liquid cocktail\n"
            r"(Ga$_{0.65}$In$_{0.20}$Sn$_{0.10}$Zn$_{0.05}$)",
            ha="center", va="center", fontsize=11, fontweight="semibold")
    ax.text(X_PT, title_y,
            "Pt substrate\n(s, fcc, SER = 0)",
            ha="center", va="center", fontsize=11, fontweight="semibold")
    ax.text(0.86, title_y,
            "HEI product\n"
            r"Pt$_3$(Ga,In,Sn,Zn) L1$_2$",
            ha="center", va="center", fontsize=11, fontweight="semibold")

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

    lines = [r"$F_i = \mu_i^\mathrm{start} - \mu_i^\mathrm{HEI}$"
             + "\n  (kJ/mol atom)\n"]
    for _, row in dmu.iterrows():
        lines.append(
            f"  {row['element']:<2s}  {row['driving_force_per_atom_kJmol']:+.2f}")
    box_txt = "\n".join(lines)
    ax.text(0.02, hei_min - y_range * 0.12, box_txt, fontsize=9.5,
            family="DejaVu Sans Mono",
            bbox=dict(facecolor="#F8F8F5", edgecolor="#BBB",
                      boxstyle="round,pad=0.5"),
            va="bottom", ha="left")

    ax.set_xlim(-0.03, 1.30)
    ax.set_ylim(hei_min - y_range * 0.18, liq_max + y_range * 0.22)
    ax.set_ylabel(r"Chemical potential per atom of f.u., $\mu_i$ (kJ$\cdot$mol$^{-1}$)",
                  fontsize=11.5)
    ax.set_xticks([0.185, X_PT, 0.86])
    ax.set_xticklabels(["start\n(liquid reservoir)",
                        "start\n(Pt substrate)",
                        "end\n(HEI sink)"], fontsize=10.5)
    ax.axhline(0, color="#999", linewidth=0.7, linestyle=":", zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_title(
        "Panel a  chemical-potential landscape (SER reference, 0 K enthalpic)\n"
        r"liquid reservoir $\to$ Pt substrate $\to$ HEI sink"
        "  |  CEF fitted to Panel g (165 UMA points)",
        fontsize=12.5, fontweight="bold", pad=12)

    plt.tight_layout()
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    out = Path(__file__).resolve().parent.parent / "outputs"
    tier, dmu, dG = load(out)
    out_path = out / "panel_a_v3_schematic"
    plot(tier, dmu, dG, out_path)
    meta = {
        "script": "script_D_panel_a_plot_v3.py",
        "inputs": ["panel_a_v3_tier_summary.csv", "delta_mu_v3_0K.csv",
                    "delta_G_rxn_v3_summary.csv"],
        "outputs": ["panel_a_v3_schematic.png", "panel_a_v3_schematic.pdf"],
        "v3_change": "UMA-CEF data; rescaled layout for larger y-range",
    }
    (out / "script_D_v3_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Panel a v3 rendered to:\n  {out_path.with_suffix('.png')}"
          f"\n  {out_path.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
