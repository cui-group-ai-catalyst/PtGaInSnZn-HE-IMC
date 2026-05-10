"""
script_FigF_TripleConsensus_Replot.py
=====================================
Lightweight replot utility — reads `data_FigF_TripleConsensus_Data.csv`
and regenerates the host-on-x / |ΔH_f|-on-y bar chart with the Miedema
overlay. Does NOT load UMA / fairchem, so it runs in the standard Python
environment shipped with this package (no GPU / heavy ML dependency).

The companion `script_FigF_TripleConsensus.py` (with `--rerun-uma`)
produces `data_FigF_TripleConsensus_Data.csv`; subsequent style changes
should be made here, in the presentation layer.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR

DATA_CSV    = SCRIPT_DIR / "data_FigF_TripleConsensus_Data.csv"
SUMMARY_CSV = SCRIPT_DIR / "data_FigF_TripleConsensus_Summary.csv"
OUT_ORIGIN  = SCRIPT_DIR / "data_FigF_TripleConsensus_Origin_replot.csv"
OUT_PLOT    = SCRIPT_DIR / "preview_FigF_TripleConsensus_replot.png"

THREE_ABS = [
    "MP_DFT_abs_Hf_kJ_mol_atom",
    "UMA_abs_Hf_kJ_mol_atom",
    "CHGNet_abs_Hf_kJ_mol_atom",
]


def load_rho_dict() -> dict:
    if not SUMMARY_CSV.exists():
        raise FileNotFoundError(f"Missing summary CSV: {SUMMARY_CSV}")
    sdf = pd.read_csv(SUMMARY_CSV)
    return {row["Metric"]: float(row["Value"]) for _, row in sdf.iterrows()}


def save_origin_wide(df: pd.DataFrame) -> None:
    """Origin-ready wide table (3-method version, Miedema dropped from display).

    Columns kept per host, all sorted by consensus rank:
      - Host, Formula, c_host, Size_Pass
      - signed Hf per method (MP-DFT / UMA / CHGNet) + the method ranks
      - abs Hf per method (lets Origin rebuild the error bars if needed)
      - Bar_Mean_absHf / Bar_Std_absHf (ready-to-plot bar + error)
      - Consensus_Rank_Mean (row ordering key, 3-method mean here)
    """
    df = df.copy()
    df["Bar_Mean_absHf_kJ_mol_atom"] = df[THREE_ABS].mean(axis=1).round(4)
    df["Bar_Std_absHf_kJ_mol_atom"]  = df[THREE_ABS].std(axis=1, ddof=1).round(4)
    # Re-compute consensus over only the 3 retained methods
    df["Consensus_Rank_Mean"] = df[["MP_DFT_Rank", "UMA_Rank", "CHGNet_Rank"]].mean(axis=1).round(2)
    df = df.sort_values(["Consensus_Rank_Mean", "UMA_Rank"]).reset_index(drop=True)

    cols = [
        "Host", "Formula", "c_host", "Size_Pass",
        "MP_DFT_Hf_kJ_mol_atom",   "MP_DFT_Rank",
        "UMA_Hf_kJ_mol_atom",      "UMA_Rank",
        "CHGNet_Hf_kJ_mol_atom",   "CHGNet_Rank",
        "MP_DFT_abs_Hf_kJ_mol_atom",
        "UMA_abs_Hf_kJ_mol_atom",
        "CHGNet_abs_Hf_kJ_mol_atom",
        "Bar_Mean_absHf_kJ_mol_atom",
        "Bar_Std_absHf_kJ_mol_atom",
        "Consensus_Rank_Mean",
    ]
    df[cols].to_csv(OUT_ORIGIN, index=False, float_format="%.4f")


def plot_consensus(df: pd.DataFrame, rho: dict) -> None:
    """Clean bar chart: hosts along X (sorted by consensus rank), bars = mean of
    (|MP-DFT|, |UMA|, |CHGNet|), error bars = std of those three.
    Miedema and individual method markers intentionally omitted per user feedback
    (Miedema magnitudes diverge too strongly from the DFT/ML family).
    """
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
    })

    order = df.sort_values(["Consensus_Rank_Mean", "UMA_Rank"]).reset_index(drop=True)
    x = np.arange(len(order))

    bar_mean = order[THREE_ABS].mean(axis=1).to_numpy()
    bar_std  = order[THREE_ABS].std(axis=1, ddof=1).to_numpy()

    fig, ax = plt.subplots(figsize=(12.0, 5.8))

    colors  = ["#C62828" if h == "Pt" else "#455A64" for h in order["Host"]]
    hatches = ["" if sp else "///" for sp in order["Size_Pass"]]

    for xi, mean_i, std_i, color, hatch in zip(x, bar_mean, bar_std, colors, hatches):
        ax.bar(xi, mean_i, width=0.72, color=color, alpha=0.90,
               edgecolor="black", linewidth=0.8, hatch=hatch, zorder=2)
        ax.errorbar(xi, mean_i, yerr=std_i, fmt="none",
                    ecolor="black", elinewidth=1.1, capsize=4, capthick=1.1,
                    zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{h}\n{f}" for h, f in zip(order["Host"], order["Formula"])],
        fontsize=8.5,
    )
    ax.set_ylabel(r"$|\Delta H_f|$  (kJ mol$^{-1}$ atom$^{-1}$)",
                  fontsize=11, fontweight="bold")
    ax.set_title(
        "Panel f (2026-04-18) | Binary M-Ga formation enthalpy magnitude by host\n"
        "Bar = mean over {MP-DFT, UMA, CHGNet};  error bar = std across the three methods",
        fontsize=10.5, fontweight="bold", pad=10,
    )
    ax.grid(axis="y", linestyle=":", alpha=0.35)
    y_max = float((bar_mean + bar_std).max()) * 1.12
    ax.set_ylim(0, y_max)
    ax.set_xlim(-0.6, len(order) - 0.4)

    handles = [
        Patch(facecolor="#C62828", edgecolor="black", label="Pt (focus host)"),
        Patch(facecolor="#455A64", edgecolor="black", label="Other hosts"),
        Patch(facecolor="white",   edgecolor="black", hatch="///",
              label="Size_Pass = False"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8.5, frameon=True)

    info = (
        f"Spearman rho (size-pass only, N={int(df['Size_Pass'].sum())}):\n"
        f"  MP-DFT  vs UMA    = {rho['MP_DFT_vs_UMA_sizepass']:+.3f}\n"
        f"  MP-DFT  vs CHGNet = {rho.get('MP_DFT_vs_CHGNet_all', float('nan')):+.3f} (all-15)\n"
        f"  UMA     vs CHGNet = {rho.get('UMA_vs_CHGNet_all',    float('nan')):+.3f} (all-15)"
    )
    ax.text(0.012, 0.98, info, transform=ax.transAxes, fontsize=8.5,
            color="#333333", va="top", ha="left",
            bbox={"facecolor": "white", "alpha": 0.92, "edgecolor": "#9E9E9E"})

    plt.tight_layout()
    fig.savefig(OUT_PLOT, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    if not DATA_CSV.exists():
        raise FileNotFoundError(f"Missing shipped data CSV: {DATA_CSV}")

    df = pd.read_csv(DATA_CSV)
    rho = load_rho_dict()

    # Redefine consensus over the 3 retained methods (Miedema no longer part of display)
    df["Consensus_Rank_Mean"] = (
        df[["MP_DFT_Rank", "UMA_Rank", "CHGNet_Rank"]].mean(axis=1).round(2)
    )
    df = df.sort_values(["Consensus_Rank_Mean", "UMA_Rank"]).reset_index(drop=True)

    save_origin_wide(df)
    print(f"[Panel f replot] Origin wide -> {OUT_ORIGIN}")

    plot_consensus(df, rho)
    print(f"[Panel f replot] preview PNG -> {OUT_PLOT}")

    bar_mean = df[THREE_ABS].mean(axis=1)
    bar_std  = df[THREE_ABS].std(axis=1, ddof=1)
    out = pd.DataFrame({
        "Host":                df["Host"],
        "Formula":              df["Formula"],
        "Bar_Mean":           bar_mean.round(3),
        "Bar_Std":            bar_std.round(3),
        "Miedema":            df["Miedema_abs_dH_kJ_mol_atom"].round(3),
        "Consensus_Rank_Mean": df["Consensus_Rank_Mean"],
    })
    print("\n-- Bar values (sorted by consensus rank) --")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
