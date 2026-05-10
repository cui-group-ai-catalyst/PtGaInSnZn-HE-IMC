"""
script_D_panel_a_plot.py
========================
Render a quantitative schematic of Panel a using the numbers produced by
Scripts A-C. The chart encodes:

  - Three tiers on y-axis (chemical potential, kJ/mol atom of species)
      Top    : liquid cocktail (Ga, In, Sn, Zn), mu_i^L ~ 0
      Middle : Pt substrate (s, fcc), mu_Pt = 0 (pure reference)
      Bottom : HEI (Pt3(Ga,In,Sn,Zn) L1_2), mu_i^HEI ~ -12.6
  - Five element lines connecting starting point -> HEI endpoint, with
    slope proportional to delta_mu_i (all ~ +12.6 kJ/mol here).
  - A main "driving-force arrow" spanning avg mu_start -> avg mu_HEI
    labelled with Delta_G_rxn per atom.

Inputs:
  outputs/panel_a_tier_summary.csv      (from Script C)
  outputs/delta_mu_0K.csv               (from Script C)
  outputs/delta_G_rxn_summary.csv       (from Script C)

Outputs:
  outputs/panel_a_schematic.png
  outputs/panel_a_schematic.pdf

Author: 2026-04-23
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
    "Ga": "#E07A5F",   # warm terracotta
    "In": "#3D5A80",   # deep blue
    "Sn": "#81B29A",   # sage green
    "Zn": "#F2CC8F",   # mustard
    "Pt": "#2D2D2D",   # near-black
}

ELEMENT_ORDER: list[str] = ["Ga", "In", "Sn", "Zn", "Pt"]

# x-axis coordinates for the three tiers; beta elements are fanned
# horizontally near the liquid and HEI ends for visual separation.
X_LIQUID_FAN: dict[str, float] = {
    "Ga": 0.08,
    "In": 0.14,
    "Sn": 0.20,
    "Zn": 0.26,
}
X_PT: float = 0.50
X_HEI_FAN: dict[str, float] = {
    "Pt": 0.74,
    "Ga": 0.80,
    "In": 0.86,
    "Sn": 0.92,
    "Zn": 0.98,
}


def load_data(out_dir: Path):
    tier = pd.read_csv(out_dir / "panel_a_tier_summary.csv")
    dmu = pd.read_csv(out_dir / "delta_mu_0K.csv")
    dG = pd.read_csv(out_dir / "delta_G_rxn_summary.csv").iloc[0]
    return tier, dmu, dG


def plot(tier: pd.DataFrame, dmu: pd.DataFrame, dG: pd.Series, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(11.0, 7.0))

    hei_mus_all = tier[tier["tier"] == "HEI_bottom"]
    hei_mu_mean = hei_mus_all["mu_kJmol"].mean()
    hei_mu_min = hei_mus_all["mu_kJmol"].min()
    hei_mu_max = hei_mus_all["mu_kJmol"].max()

    # ---- tier shaded bands ----
    # Top band: liquid reservoir (narrow, around 0)
    ax.axhspan(-0.35, 0.15, xmin=0.02, xmax=0.33,
               facecolor="#EAE6D9", alpha=0.55, zorder=1)
    # Middle tick: Pt solid
    ax.axhspan(-0.08, 0.08, xmin=0.43, xmax=0.57,
               facecolor="#DDDDDD", alpha=0.55, zorder=1)
    # Bottom band: HEI sink
    ax.axhspan(hei_mu_min - 0.12, hei_mu_max + 0.12, xmin=0.66, xmax=1.02,
               facecolor="#CFE0CF", alpha=0.55, zorder=1)

    # ---- element trajectories (start -> HEI) ----
    for el in ELEMENT_ORDER:
        color = ELEMENT_COLORS[el]
        mu_hei = hei_mus_all[hei_mus_all["element"] == el]["mu_kJmol"].values[0]
        x_end = X_HEI_FAN[el]

        if el == "Pt":
            mu_start = tier[tier["tier"] == "Pt_middle"]["mu_kJmol"].values[0]
            x_start = X_PT
        else:
            mu_start = tier[(tier["tier"] == "liquid_top") & (tier["element"] == el)][
                "mu_kJmol"
            ].values[0]
            x_start = X_LIQUID_FAN[el]

        # trajectory line
        ax.plot(
            [x_start, x_end], [mu_start, mu_hei],
            color=color, linestyle=(0, (6, 3)), linewidth=1.8, alpha=0.85, zorder=3,
        )
        # start marker
        ax.scatter([x_start], [mu_start], s=120, color=color, edgecolor="white",
                   linewidths=1.4, zorder=6)
        # end marker
        ax.scatter([x_end], [mu_hei], s=140, color=color, edgecolor="white",
                   linewidths=1.4, zorder=6)

        # start label (above point for beta, below for Pt)
        if el == "Pt":
            ax.annotate(el, xy=(x_start, mu_start),
                        xytext=(x_start, mu_start - 0.9),
                        fontsize=11.5, fontweight="bold", color=color,
                        ha="center", va="top")
        else:
            ax.annotate(el, xy=(x_start, mu_start),
                        xytext=(x_start, mu_start + 0.55),
                        fontsize=11.5, fontweight="bold", color=color,
                        ha="center", va="bottom")
        # end label (staggered)
        dy_end = 0.55 if el in ("Ga", "Sn") else -0.75
        va_end = "bottom" if el in ("Ga", "Sn") else "top"
        ax.annotate(el, xy=(x_end, mu_hei),
                    xytext=(x_end, mu_hei + dy_end),
                    fontsize=11.0, fontweight="bold", color=color,
                    ha="center", va=va_end)

    # ---- tier titles ----
    ax.text(0.17, 3.7,
            "Liquid cocktail\n(Ga 0.65 / In 0.20 / Sn 0.10 / Zn 0.05)",
            ha="center", va="center", fontsize=11, fontweight="semibold")
    ax.text(X_PT, 3.7,
            "Pt substrate\n(s, fcc)",
            ha="center", va="center", fontsize=11, fontweight="semibold")
    ax.text(0.86, 3.7,
            "HEI product\n$\\mathrm{Pt_3(Ga,In,Sn,Zn)}$ (L1$_2$)",
            ha="center", va="center", fontsize=11, fontweight="semibold")

    # ---- driving-force arrow (DOWNWARD to indicate ?G < 0, i.e. spontaneous) ----
    dG_per_atom = float(dG["delta_G_rxn_per_atom_kJmol"])
    dG_per_fu = float(dG["delta_G_rxn_per_fu_kJmol"])

    arrow = FancyArrowPatch(
        posA=(1.04, 0.0),
        posB=(1.04, hei_mu_mean),
        arrowstyle="-|>",
        mutation_scale=24,
        color="#C8102E",
        linewidth=2.4,
        zorder=7,
    )
    ax.add_patch(arrow)
    ax.text(
        1.065, hei_mu_mean / 2,
        (f"$\\Delta G_\\mathrm{{rxn}}$\n"
         f"= {dG_per_atom:+.1f}\n  kJ/mol atom\n"
         f"= {dG_per_fu:+.1f}\n  kJ/mol f.u.\n\n"
         f"(spontaneous\n formation\n of HEI)"),
        ha="left", va="center", fontsize=9.5, color="#C8102E", fontweight="semibold",
    )

    # ---- axis cosmetics ----
    ax.set_xlim(-0.02, 1.25)
    ax.set_ylim(hei_mu_min - 3.5, 5.2)
    ax.set_ylabel(r"Chemical potential $\mu_i$ (kJ/mol atom of species)",
                  fontsize=11.5)
    ax.set_xticks([0.17, X_PT, 0.86])
    ax.set_xticklabels([
        "start\n(liquid reservoir)",
        "start\n(Pt substrate)",
        "end\n(HEI sink)",
    ], fontsize=10.5)
    ax.axhline(0, color="#999", linewidth=0.7, linestyle=":", zorder=0)
    ax.axhline(hei_mu_mean, color="#5A8F5A", linewidth=0.7,
               linestyle=":", alpha=0.6, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ---- title ----
    ax.set_title(
        "Panel a  chemical-potential landscape: liquid reservoir $\\to$ HEI sink\n"
        "(0 K enthalpic approximation; Ga$_{0.65}$In$_{0.20}$Sn$_{0.10}$Zn$_{0.05}$ cocktail)",
        fontsize=12.5, fontweight="bold", pad=12,
    )

    # ---- Delta_mu summary box (inside axes, lower-left) ----
    rows_txt = "\n".join(
        f"  {row['element']:<2s}  {row['delta_mu_kJmol']:+.2f}"
        for _, row in dmu.iterrows()
    )
    summary_txt = (
        r"$\Delta\mu_i \equiv \mu_i^{(0)} - \mu_i^\mathrm{HEI}$" + "\n" +
        "  (kJ/mol atom)\n" + rows_txt
    )
    ax.text(
        0.03, hei_mu_min - 3.3, summary_txt, fontsize=9.5,
        family="DejaVu Sans Mono",
        bbox=dict(facecolor="#F8F8F5", edgecolor="#BBBBBB", boxstyle="round,pad=0.5"),
        va="bottom", ha="left",
    )

    plt.tight_layout()
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    out_dir = Path(__file__).resolve().parent.parent / "outputs"
    tier, dmu, dG = load_data(out_dir)

    out_path = out_dir / "panel_a_schematic"
    plot(tier, dmu, dG, out_path)

    meta = {
        "script": "script_D_panel_a_plot.py",
        "inputs": [
            "panel_a_tier_summary.csv",
            "delta_mu_0K.csv",
            "delta_G_rxn_summary.csv",
        ],
        "outputs": [str(out_path.with_suffix(".png").name), str(out_path.with_suffix(".pdf").name)],
        "description": "Panel a quantitative schematic (three tiers + element trajectories)",
    }
    (out_dir / "script_D_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("Rendered Panel a schematic to:")
    print(f"  {out_path.with_suffix('.png')}")
    print(f"  {out_path.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
