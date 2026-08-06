"""Recalculate Fig. 3h from the exact Layer 8 ROI used for canonical g2.

The analysis uses lattice geometry and measured image intensities only. Existing
manual bright/dim labels are deliberately excluded from grouping and fitting.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import rankdata

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import fig3_paths


PX_NM = 0.01138848395
INK = "#202124"
CLASS0 = "#C83E4D"
CLASS1 = "#2A7F9E"
GRID = "#D6DADD"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def auc(values: np.ndarray, positive: np.ndarray) -> float:
    values = np.asarray(values, float)
    positive = np.asarray(positive, bool)
    n1 = int(positive.sum())
    n0 = int((~positive).sum())
    ranks = rankdata(values)
    return float((ranks[positive].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def derive_g2_roi(labels: pd.DataFrame, size: int = 320) -> dict[str, float | int]:
    ordered = labels[labels["spatial_region"].eq("ordered-intermetallic")]
    cx = float(ordered["x_px"].median())
    cy = float(ordered["y_px"].median())
    x0 = int(np.clip(round(cx - size / 2), 0, 1315 - size))
    y0 = int(np.clip(round(cy - size / 2), 0, 1315 - size))
    return {"cx_px": cx, "cy_px": cy, "x0_px": x0, "y0_px": y0, "size_px": size}


def register_columns(columns: pd.DataFrame, model_path: Path) -> pd.DataFrame:
    model = np.load(model_path)
    origin = np.asarray(model["origin"], float)
    basis = np.asarray(model["basis"], float)
    xy = columns[["x_px", "y_px"]].to_numpy(float)
    ij = np.rint(np.linalg.solve(basis, xy.T - origin[:, None]).T).astype(int)
    predicted = origin + ij @ basis.T
    result = columns.copy()
    result["lattice_i"] = ij[:, 0]
    result["lattice_j"] = ij[:, 1]
    result["projected_row"] = ij[:, 0] + ij[:, 1]
    result["row_class"] = result["projected_row"] & 1
    result["lattice_fit_residual_px"] = np.linalg.norm(xy - predicted, axis=1)
    return result


def integrate_columns(
    image: np.ndarray,
    columns: pd.DataFrame,
    aperture_radius: int,
    background_inner: int,
    background_outer: int,
    trend_radius: float,
) -> pd.DataFrame:
    half = background_outer
    yy, xx = np.mgrid[-half : half + 1, -half : half + 1]
    radius2 = xx * xx + yy * yy
    aperture = radius2 <= aperture_radius**2
    annulus = (radius2 >= background_inner**2) & (radius2 <= background_outer**2)
    integrated = []
    backgrounds = []
    for x, y in columns[["x_px", "y_px"]].to_numpy(float):
        xi, yi = int(round(x)), int(round(y))
        patch = image[yi - half : yi + half + 1, xi - half : xi + half + 1]
        background = float(np.median(patch[annulus]))
        integrated.append(float(np.sum(patch[aperture] - background)))
        backgrounds.append(background)

    result = columns.copy()
    result["local_background"] = backgrounds
    result["background_subtracted_integrated_intensity"] = integrated
    xy = result[["x_px", "y_px"]].to_numpy(float)
    values = result["background_subtracted_integrated_intensity"].to_numpy(float)
    trend = []
    for point in xy:
        local = np.sum((xy - point) ** 2, axis=1) <= trend_radius**2
        trend.append(float(np.median(values[local])))
    result["local_spatial_trend"] = trend
    result["normalized_phase_intensity"] = values / np.asarray(trend)
    return result


def summarize_rows(columns: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row_id, group in columns.groupby("projected_row", sort=True):
        values = group["normalized_phase_intensity"].to_numpy(float)
        rows.append(
            {
                "projected_row": int(row_id),
                "row_class": int(row_id) & 1,
                "n_columns": int(len(group)),
                "median_normalized_intensity": float(np.median(values)),
                "q25_normalized_intensity": float(np.percentile(values, 25)),
                "q75_normalized_intensity": float(np.percentile(values, 75)),
                "x_min_px": float(group["x_px"].min()),
                "x_max_px": float(group["x_px"].max()),
            }
        )
    return pd.DataFrame(rows)


def pair_rows(row_summary: pd.DataFrame) -> pd.DataFrame:
    by_row = row_summary.set_index("projected_row")
    rows = []
    for even_row in sorted(row_summary.loc[row_summary.row_class.eq(0), "projected_row"]):
        odd_row = int(even_row + 1)
        if odd_row not in by_row.index:
            continue
        i0 = float(by_row.loc[even_row, "median_normalized_intensity"])
        i1 = float(by_row.loc[odd_row, "median_normalized_intensity"])
        rows.append(
            {
                "class0_row": int(even_row),
                "class1_row": odd_row,
                "class0_median": i0,
                "class1_median": i1,
                "class0_to_class1_ratio": i0 / i1,
                "signed_fractional_contrast": (i0 - i1) / ((i0 + i1) / 2),
            }
        )
    return pd.DataFrame(rows)


def spatial_blocks(columns: pd.DataFrame, roi: dict[str, float | int]) -> pd.DataFrame:
    result = columns.copy()
    x_edges = np.linspace(float(roi["x0_px"]), float(roi["x0_px"]) + 320, 4)
    y_edges = np.linspace(float(roi["y0_px"]), float(roi["y0_px"]) + 320, 3)
    result["x_block"] = np.clip(np.digitize(result.x_px, x_edges[1:-1]), 0, 2)
    result["y_block"] = np.clip(np.digitize(result.y_px, y_edges[1:-1]), 0, 1)
    rows = []
    for (xb, yb), group in result.groupby(["x_block", "y_block"]):
        p0 = group.loc[group.row_class.eq(0), "normalized_phase_intensity"].to_numpy(float)
        p1 = group.loc[group.row_class.eq(1), "normalized_phase_intensity"].to_numpy(float)
        if min(len(p0), len(p1)) < 5:
            continue
        rows.append(
            {
                "x_block": int(xb),
                "y_block": int(yb),
                "n_class0": int(len(p0)),
                "n_class1": int(len(p1)),
                "median_class0": float(np.median(p0)),
                "median_class1": float(np.median(p1)),
                "class0_to_class1_ratio": float(np.median(p0) / np.median(p1)),
                "auc_class0_brighter": auc(np.r_[p0, p1], np.r_[np.ones(len(p0)), np.zeros(len(p1))]),
            }
        )
    return pd.DataFrame(rows)


def bootstrap_pair_median(pairs: pd.DataFrame, seed: int = 20260805) -> tuple[float, float, float]:
    values = pairs["signed_fractional_contrast"].to_numpy(float)
    rng = np.random.default_rng(seed)
    samples = np.median(rng.choice(values, size=(20000, len(values)), replace=True), axis=1)
    return float(np.median(values)), float(np.percentile(samples, 2.5)), float(np.percentile(samples, 97.5))


def sensitivity_analysis(
    image: np.ndarray,
    registered: pd.DataFrame,
    roi: dict[str, float | int],
) -> pd.DataFrame:
    rows = []
    for aperture in (3, 4, 5):
        for trend_radius in (70.0, 90.0, 110.0):
            for shrink in (0, 10, 20):
                x0 = int(roi["x0_px"]) + max(9, shrink)
                y0 = int(roi["y0_px"]) + max(9, shrink)
                x1 = int(roi["x0_px"]) + 319 - max(9, shrink)
                y1 = int(roi["y0_px"]) + 319 - max(9, shrink)
                subset = registered[
                    registered.x_px.between(x0, x1)
                    & registered.y_px.between(y0, y1)
                    & registered.lattice_fit_residual_px.le(5.0)
                ].copy()
                measured = integrate_columns(image, subset, aperture, 6, 9, trend_radius)
                row_summary = summarize_rows(measured)
                pairs = pair_rows(row_summary)
                p0 = measured.loc[measured.row_class.eq(0), "normalized_phase_intensity"].to_numpy(float)
                p1 = measured.loc[measured.row_class.eq(1), "normalized_phase_intensity"].to_numpy(float)
                rows.append(
                    {
                        "aperture_radius_px": aperture,
                        "trend_radius_px": trend_radius,
                        "roi_shrink_px": shrink,
                        "n_columns": int(len(measured)),
                        "n_row_pairs": int(len(pairs)),
                        "median_class0": float(np.median(p0)),
                        "median_class1": float(np.median(p1)),
                        "class0_to_class1_ratio": float(np.median(p0) / np.median(p1)),
                        "auc_class0_brighter": auc(np.r_[p0, p1], np.r_[np.ones(len(p0)), np.zeros(len(p1))]),
                        "row_pairs_class0_brighter": int((pairs.signed_fractional_contrast > 0).sum()),
                        "median_pair_contrast": float(pairs.signed_fractional_contrast.median()),
                    }
                )
    return pd.DataFrame(rows)


def render_figure(
    output: Path,
    image: np.ndarray,
    columns: pd.DataFrame,
    row_summary: pd.DataFrame,
    pairs: pd.DataFrame,
    roi: dict[str, float | int],
    ci: tuple[float, float, float],
) -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.5,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig = plt.figure(figsize=(6.9, 3.15), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.08, 0.92], left=0.055, right=0.985,
                          bottom=0.17, top=0.91, wspace=0.28)
    ax_image = fig.add_subplot(gs[0, 0])
    ax_rows = fig.add_subplot(gs[0, 1])

    x0, y0, size = int(roi["x0_px"]), int(roi["y0_px"]), int(roi["size_px"])
    crop = image[y0 : y0 + size, x0 : x0 + size]
    lo, hi = np.percentile(crop, [1, 99.7])
    ax_image.imshow(crop, cmap="gray", vmin=lo, vmax=hi, interpolation="nearest")
    for row_id, group in columns.groupby("projected_row"):
        local_x = group.x_px.to_numpy() - x0
        local_y = group.y_px.to_numpy() - y0
        order = np.argsort(local_x)
        color = CLASS0 if int(row_id) % 2 == 0 else CLASS1
        if len(group) >= 3:
            ax_image.plot(local_x[order], local_y[order], color=color, lw=0.75, alpha=0.55)
        ax_image.scatter(local_x, local_y, s=13, facecolors="none", edgecolors=color,
                         linewidths=0.6, alpha=0.9)
    bar_px = 1.0 / PX_NM
    ax_image.plot([18, 18 + bar_px], [size - 18, size - 18], color="white", lw=2.5)
    ax_image.text(18 + bar_px / 2, size - 25, "1 nm", color="white", ha="center", va="bottom",
                  fontsize=7, fontweight="bold")
    ax_image.set_title(r"Layer 8 g2 ROI: projected $[110]$ rows", pad=5)
    ax_image.set_xticks([])
    ax_image.set_yticks([])
    for spine in ax_image.spines.values():
        spine.set_visible(False)

    for row_class, color, label in ((0, CLASS0, "row class 0"), (1, CLASS1, "row class 1")):
        one = row_summary[row_summary.row_class.eq(row_class)]
        ax_rows.errorbar(
            one.projected_row,
            one.median_normalized_intensity,
            yerr=[
                one.median_normalized_intensity - one.q25_normalized_intensity,
                one.q75_normalized_intensity - one.median_normalized_intensity,
            ],
            fmt="o",
            color=color,
            ms=4.2,
            capsize=1.8,
            elinewidth=0.8,
            label=label,
            zorder=3,
        )
    ax_rows.plot(row_summary.projected_row, row_summary.median_normalized_intensity,
                 color="#8A9094", lw=0.75, zorder=1)
    ax_rows.axhline(1.0, color=GRID, lw=0.8, zorder=0)
    ax_rows.set_xlabel("Projected row index")
    ax_rows.set_ylabel("Normalized phase intensity")
    ax_rows.set_xticks(row_summary.projected_row.iloc[::2])
    ax_rows.tick_params(length=2.5)
    ax_rows.spines[["top", "right"]].set_visible(False)
    ax_rows.legend(frameon=False, fontsize=6.7, loc="lower left")
    median_contrast, low, high = ci
    ax_rows.text(
        0.98,
        0.97,
        f"class 0 > class 1: {(pairs.signed_fractional_contrast > 0).sum()}/{len(pairs)} pairs\n"
        f"median contrast = {median_contrast:.3f}\n"
        f"descriptive bootstrap interval [{low:.3f}, {high:.3f}]",
        transform=ax_rows.transAxes,
        ha="right",
        va="top",
        fontsize=6.7,
        bbox=dict(boxstyle="square,pad=0.3", fc="white", ec=GRID, lw=0.6, alpha=0.92),
    )
    ax_rows.set_title("Row-resolved real-space contrast", pad=5)

    ax_image.text(-0.07, 1.04, "h1", transform=ax_image.transAxes, fontsize=9.5,
                  fontweight="bold", color=INK)
    ax_rows.text(-0.13, 1.04, "h2", transform=ax_rows.transAxes, fontsize=9.5,
                 fontweight="bold", color=INK)
    for ext in ("png", "pdf", "svg", "tif"):
        kwargs = {"dpi": 600} if ext in {"png", "tif"} else {}
        if ext == "tif":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(output / f"figure3h_l8_projected_row_contrast.{ext}", **kwargs)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=fig3_paths.OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    labels_path = fig3_paths.DATA / "l8_three_region_intensity_labels.csv"
    columns_path = fig3_paths.DATA / "corrected_atom_columns_intensity.csv"
    model_path = fig3_paths.DATA / "complete_lattice_model.npz"
    image_path = fig3_paths.DATA / "gray8_plain_layer_08_of_16_scale2nm.tif"

    labels = pd.read_csv(labels_path)
    source_columns = pd.read_csv(columns_path)
    image = np.asarray(Image.open(image_path), float)
    roi = derive_g2_roi(labels)
    margin = 9
    x0, y0 = int(roi["x0_px"]), int(roi["y0_px"])
    in_roi = source_columns[
        source_columns.x_px.between(x0 + margin, x0 + 319 - margin)
        & source_columns.y_px.between(y0 + margin, y0 + 319 - margin)
    ].copy()
    registered = register_columns(in_roi, model_path)
    registered = registered[registered.lattice_fit_residual_px.le(5.0)].copy()
    measured = integrate_columns(image, registered, 4, 6, 9, 90.0)
    row_summary = summarize_rows(measured)
    pairs = pair_rows(row_summary)
    blocks = spatial_blocks(measured, roi)
    sensitivity = sensitivity_analysis(image, register_columns(in_roi, model_path), roi)
    ci = bootstrap_pair_median(pairs)

    measured.drop(columns=[c for c in ("manual_atom_class",) if c in measured]).to_csv(
        args.output / "h_column_source_data.csv", index=False
    )
    row_summary.to_csv(args.output / "h_row_summary.csv", index=False)
    pairs.to_csv(args.output / "h_adjacent_row_pairs.csv", index=False)
    blocks.to_csv(args.output / "h_spatial_block_audit.csv", index=False)
    sensitivity.to_csv(args.output / "h_sensitivity_audit.csv", index=False)

    p0 = measured.loc[measured.row_class.eq(0), "normalized_phase_intensity"].to_numpy(float)
    p1 = measured.loc[measured.row_class.eq(1), "normalized_phase_intensity"].to_numpy(float)
    result = {
        "panel_g_provenance": {
            "g1": "Layer 8, Pt-rich region",
            "g2": "Layer 8, ordered-intermetallic region",
            "g3": "Layer 16, liquid-metal-rich region",
        },
        "g2_roi": roi,
        "analysis_definition": {
            "lattice_registration": "fixed geometric origin and basis; nearest integer i,j; residual <= 5 px",
            "row_definition": "projected [110] rows are i+j = constant; row class is (i+j) mod 2",
            "intensity": "4 px circular aperture minus median 6-9 px annular background",
            "normalization": "divide by 90 px local median spatial trend using both row classes",
            "manual_atom_class_labels_used_in_grouping_or_intensity_analysis": False,
            "upstream_geometry_note": (
                "Measured column positions and the fixed lattice geometry are inherited from the "
                "audited reconstruction; only atom-class labels are excluded here."
            ),
        },
        "primary_results": {
            "n_columns": int(len(measured)),
            "n_contiguous_rows": int(row_summary.projected_row.nunique()),
            "n_adjacent_row_pairs": int(len(pairs)),
            "pairs_class0_brighter": int((pairs.signed_fractional_contrast > 0).sum()),
            "median_pair_signed_fractional_contrast": ci[0],
            "bootstrap_95_interval_technical": [ci[1], ci[2]],
            "spatial_blocks_class0_brighter": int((blocks.class0_to_class1_ratio > 1).sum()),
            "n_spatial_blocks": int(len(blocks)),
        },
        "secondary_descriptive_results": {
            "n_class0": int(len(p0)),
            "n_class1": int(len(p1)),
            "median_class0": float(np.median(p0)),
            "median_class1": float(np.median(p1)),
            "median_ratio_class0_to_class1": float(np.median(p0) / np.median(p1)),
            "auc_class0_brighter": auc(np.r_[p0, p1], np.r_[np.ones(len(p0)), np.zeros(len(p1))]),
        },
        "sensitivity": {
            "n_parameter_combinations": int(len(sensitivity)),
            "all_median_ratios_above_one": bool((sensitivity.class0_to_class1_ratio > 1).all()),
            "median_ratio_range": [
                float(sensitivity.class0_to_class1_ratio.min()),
                float(sensitivity.class0_to_class1_ratio.max()),
            ],
            "auc_range": [
                float(sensitivity.auc_class0_brighter.min()),
                float(sensitivity.auc_class0_brighter.max()),
            ],
            "minimum_fraction_of_row_pairs_with_same_direction": float(
                (sensitivity.row_pairs_class0_brighter / sensitivity.n_row_pairs).min()
            ),
        },
        "scientific_boundary": (
            "This analysis tests whether the g2 reciprocal-space modulation has a repeated, "
            "row-specific real-space intensity counterpart consistent with the [110] projected "
            "L12 motif. It does not identify individual elements or independently determine "
            "site occupancies."
        ),
        "input_sha256": {
            str(labels_path): sha256(labels_path),
            str(columns_path): sha256(columns_path),
            str(model_path): sha256(model_path),
            str(image_path): sha256(image_path),
        },
    }
    (args.output / "h_analysis_results.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    render_figure(args.output, image, measured, row_summary, pairs, roi, ci)
    print(json.dumps(result["primary_results"], indent=2))
    print(json.dumps(result["secondary_descriptive_results"], indent=2))


if __name__ == "__main__":
    main()
