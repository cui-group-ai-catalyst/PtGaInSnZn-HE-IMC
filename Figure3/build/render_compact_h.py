"""Render the compact single-panel Figure 3h row-profile plot."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import fig3_paths

ROOT = fig3_paths.OUTPUT
SOURCE = fig3_paths.SOURCE_DATA / "h_row_summary.csv"
ROW_REPEAT_NM = 0.3984872634717192
INK = "#202124"
FAMILY_A = "#242628"
FAMILY_B = "#777D80"
NEUTRAL = "#B2B7B9"


def main() -> None:
    data = pd.read_csv(SOURCE)

    # Sixteen complete central rows; the partial edge row is retained upstream only.
    data = data[data.projected_row.between(-9, 6) & data.n_columns.ge(10)].copy()
    data = data.sort_values("projected_row").reset_index(drop=True)
    data["display_row"] = np.arange(1, len(data) + 1)
    data["row_family"] = np.where(data.row_class.eq(0), "A", "B")

    pairs = []
    indexed = data.set_index("projected_row")
    for low_row in range(-9, 6, 2):
        high_row = low_row + 1
        low = float(indexed.loc[low_row, "median_normalized_intensity"])
        high = float(indexed.loc[high_row, "median_normalized_intensity"])
        pairs.append(
            {
                "family_B_row": low_row,
                "family_A_row": high_row,
                "family_B_median": low,
                "family_A_median": high,
                "signed_fractional_contrast": (high - low) / ((high + low) / 2),
            }
        )
    pairs = pd.DataFrame(pairs)
    median_contrast = float(pairs.signed_fractional_contrast.median())

    data.to_csv(ROOT / "h_compact_source_data.csv", index=False)
    pairs.to_csv(ROOT / "h_compact_row_pairs.csv", index=False)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 6.0,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 2.2,
            "ytick.major.size": 2.2,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )

    width_in = 62 / 25.4
    height_in = 30 / 25.4
    fig, ax = plt.subplots(figsize=(width_in, height_in), facecolor="white")
    fig.subplots_adjust(left=0.19, right=0.96, bottom=0.26, top=0.78)

    x = data.display_row.to_numpy(float)
    y = data.median_normalized_intensity.to_numpy(float)
    ax.plot(x, y, color=NEUTRAL, lw=0.75, zorder=1)

    for family, color, marker, facecolor in (
        ("A", FAMILY_A, "o", FAMILY_A),
        ("B", FAMILY_B, "s", "white"),
    ):
        one = data[data.row_family.eq(family)]
        lower = one.median_normalized_intensity - one.q25_normalized_intensity
        upper = one.q75_normalized_intensity - one.median_normalized_intensity
        ax.errorbar(
            one.display_row,
            one.median_normalized_intensity,
            yerr=np.vstack([lower, upper]),
            fmt=marker,
            ms=3.1,
            color=color,
            markerfacecolor=facecolor,
            markeredgecolor=color,
            markeredgewidth=0.75,
            ecolor=color,
            elinewidth=0.75,
            capsize=1.4,
            capthick=0.65,
            zorder=3,
        )

    ax.axhline(1.0, color="#D5D9DB", lw=0.65, zorder=0)
    ax.set_xlim(0.4, 17.0)
    ax.set_ylim(0.82, 1.20)
    ax.set_xticks([1, 4, 7, 10, 13, 16])
    ax.set_yticks([0.8, 1.0, 1.2])
    ax.set_xlabel("Projected row", labelpad=2)
    ax.set_ylabel("Normalized phase", labelpad=2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(16.45, float(data.iloc[-1].median_normalized_intensity), "A", color=FAMILY_A,
            ha="left", va="center", fontsize=5.5, fontweight="bold")
    ax.text(15.45, float(data.iloc[-2].median_normalized_intensity), "B", color=FAMILY_B,
            ha="left", va="center", fontsize=5.5, fontweight="bold")

    # Same-family repeat: two adjacent projected-row spacings.
    y_bracket = 1.212
    ax.plot([2, 4], [y_bracket, y_bracket], color=INK, lw=0.65, clip_on=False)
    ax.plot([2, 2], [y_bracket - 0.006, y_bracket], color=INK, lw=0.65, clip_on=False)
    ax.plot([4, 4], [y_bracket - 0.006, y_bracket], color=INK, lw=0.65, clip_on=False)
    ax.text(3, y_bracket + 0.004, f"{ROW_REPEAT_NM:.2f} nm repeat", ha="center", va="bottom",
            fontsize=5.15, color=INK, clip_on=False)

    ax.text(
        0.99,
        1.115,
        f"8/8 pairs | median contrast {median_contrast * 100:.1f}%",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.15,
        color=INK,
        clip_on=False,
    )
    ax.text(-0.20, 1.12, "h", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.0, fontweight="bold", color=INK)

    stem = ROOT / "figure3h_compact_row_profile"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".tif"), dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)

    print(f"rows={len(data)}, columns={int(data.n_columns.sum())}, pairs={len(pairs)}")
    print(f"median_contrast={median_contrast:.6f}")


if __name__ == "__main__":
    main()
