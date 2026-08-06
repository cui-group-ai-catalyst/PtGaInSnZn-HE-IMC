"""Audit atom coordinates for all 16 registered layers before 3D rendering.

The corrected layer-8 two-sublattice map is used as the lateral template.
For each layer, template sites inside the current spherical radius are refined
to a nearby local intensity maximum.  DeepSeek's existing per-layer detections
are compared to the refined template but are never silently overwritten.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree


PT_COLOR = "#D85C41"
LM_COLOR = "#3E78B2"
ACCEPT_COLOR = "#18D8E8"
REVIEW_COLOR = "#FF9F1C"
REJECT_COLOR = "#9B9B9B"
EXTRA_COLOR = "#FF2E88"


def refine_site(image: np.ndarray, enhanced: np.ndarray, px: float, py: float, search_radius: int = 7) -> dict:
    height, width = image.shape
    ix, iy = int(round(px)), int(round(py))
    if ix < 12 or iy < 12 or ix >= width - 12 or iy >= height - 12:
        return {"x": px, "y": py, "shift": np.nan, "peak": np.nan, "prominence": np.nan, "snr": np.nan}

    search = enhanced[iy - search_radius : iy + search_radius + 1, ix - search_radius : ix + search_radius + 1]
    ry, rx = np.unravel_index(np.argmax(search), search.shape)
    peak_x = ix - search_radius + rx
    peak_y = iy - search_radius + ry

    radius = 3
    patch = enhanced[peak_y - radius : peak_y + radius + 1, peak_x - radius : peak_x + radius + 1]
    weights = np.clip(patch - np.percentile(patch, 20), 0, None)
    yy, xx = np.indices(patch.shape)
    total = float(np.sum(weights))
    if total > 0:
        refined_x = peak_x - radius + float(np.sum(xx * weights) / total)
        refined_y = peak_y - radius + float(np.sum(yy * weights) / total)
    else:
        refined_x, refined_y = float(peak_x), float(peak_y)

    local = enhanced[peak_y - 10 : peak_y + 11, peak_x - 10 : peak_x + 11]
    ly, lx = np.indices(local.shape)
    rr = np.hypot(lx - 10, ly - 10)
    annulus = local[(rr >= 6) & (rr <= 10)]
    background = float(np.median(annulus))
    noise = float(1.4826 * np.median(np.abs(annulus - background)))
    peak_value = float(enhanced[peak_y, peak_x])
    prominence = peak_value - background
    snr = prominence / max(noise, 1e-6)
    shift = float(np.hypot(refined_x - px, refined_y - py))
    return {"x": refined_x, "y": refined_y, "shift": shift, "peak": peak_value, "prominence": prominence, "snr": snr}


def status_from_metrics(shift: float, snr: float) -> str:
    if np.isfinite(shift) and np.isfinite(snr) and shift <= 5.5 and snr >= 2.0:
        return "accept"
    if np.isfinite(shift) and np.isfinite(snr) and shift <= 8.0 and snr >= 1.0:
        return "review"
    return "reject"


def draw_circle(ax: plt.Axes, cx: float, cy: float, radius_px: float) -> None:
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(cx + radius_px * np.cos(theta), cy + radius_px * np.sin(theta), color="red", lw=1.0, ls="--")
    ax.plot(cx, cy, marker="+", ms=7, mew=1.0, color="red")


def clean_axis(ax: plt.Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(args.workspace))
    from particle_geometry import CX, CY, PX_NM, R_NM

    stack_data = np.load(args.workspace / "output" / "registered_stack.npz")
    stack = np.asarray(stack_data["registered"], dtype=float)
    db = np.load(args.workspace / "output" / "atoms_db.npz", allow_pickle=True)
    deepseek_metrics = [np.asarray(layer) for layer in db["metrics"]]
    template = pd.read_csv(args.template)
    template_x = template["x_px"].to_numpy(float)
    template_y = template["y_px"].to_numpy(float)
    template_class = template["manual_atom_class"].to_numpy(str)
    template_ids = template["atom_id"].to_numpy(int)

    detail_rows = []
    summary_rows = []
    layer_plot_data = []

    for layer_index in range(16):
        image = stack[layer_index]
        enhanced = gaussian_filter(image, 1.1) - gaussian_filter(image, 9.0)
        radius_px = float(R_NM[layer_index] / PX_NM)
        template_radius = np.hypot(template_x - CX, template_y - CY)
        inside_template = template_radius <= radius_px

        deepseek = deepseek_metrics[layer_index]
        deepseek_radius = np.hypot(deepseek[:, 0] - CX, deepseek[:, 1] - CY)
        deepseek_inside = deepseek[deepseek_radius <= radius_px]

        layer_template_ids = template_ids[inside_template]
        layer_classes = template_class[inside_template]
        initial_x = template_x[inside_template]
        initial_y = template_y[inside_template]
        refined = [refine_site(image, enhanced, px, py) for px, py in zip(initial_x, initial_y)]
        refined_x = np.asarray([row["x"] for row in refined])
        refined_y = np.asarray([row["y"] for row in refined])
        shifts = np.asarray([row["shift"] for row in refined])
        peaks = np.asarray([row["peak"] for row in refined])
        prominences = np.asarray([row["prominence"] for row in refined])
        snr = np.asarray([row["snr"] for row in refined])
        statuses = np.asarray([status_from_metrics(shift, signal) for shift, signal in zip(shifts, snr)])

        # Remove duplicate assignments to the same local maximum from the accepted set.
        valid_for_duplicates = np.isfinite(refined_x) & np.isfinite(refined_y)
        duplicate = np.zeros(len(refined_x), dtype=bool)
        candidate_indices = np.where(valid_for_duplicates)[0]
        order = candidate_indices[np.argsort(-np.nan_to_num(snr[candidate_indices], nan=-np.inf))]
        kept_indices = []
        for idx in order:
            if kept_indices:
                distances = np.hypot(refined_x[kept_indices] - refined_x[idx], refined_y[kept_indices] - refined_y[idx])
                if np.min(distances) < 5.0:
                    duplicate[idx] = True
                    continue
            kept_indices.append(idx)
        statuses[duplicate] = "reject"

        refined_tree = cKDTree(np.c_[refined_x, refined_y]) if len(refined_x) else None
        if len(deepseek_inside) and refined_tree is not None:
            deepseek_to_template, _ = refined_tree.query(deepseek_inside[:, :2])
            deepseek_extra = deepseek_to_template > 6.0
        else:
            deepseek_to_template = np.zeros(len(deepseek_inside))
            deepseek_extra = np.ones(len(deepseek_inside), dtype=bool)

        if len(deepseek_inside):
            deepseek_tree = cKDTree(deepseek_inside[:, :2])
            template_to_deepseek, _ = deepseek_tree.query(np.c_[refined_x, refined_y])
        else:
            template_to_deepseek = np.full(len(refined_x), np.inf)
        missed_by_deepseek = template_to_deepseek > 6.0

        accepted = statuses == "accept"
        review = statuses == "review"
        rejected = statuses == "reject"

        if np.any(accepted):
            weight = np.clip(prominences[accepted], 0, None)
            if np.sum(weight) <= 0:
                weight = np.ones(np.sum(accepted))
            centroid_x = float(np.average(refined_x[accepted], weights=weight))
            centroid_y = float(np.average(refined_y[accepted], weights=weight))
        else:
            centroid_x = centroid_y = np.nan

        for idx in range(len(refined_x)):
            detail_rows.append(
                {
                    "layer": layer_index + 1,
                    "template_atom_id": int(layer_template_ids[idx]),
                    "manual_atom_class": layer_classes[idx],
                    "x_template_px": initial_x[idx],
                    "y_template_px": initial_y[idx],
                    "x_refined_px": refined_x[idx],
                    "y_refined_px": refined_y[idx],
                    "refinement_shift_px": shifts[idx],
                    "peak_bandpass": peaks[idx],
                    "local_prominence": prominences[idx],
                    "local_snr": snr[idx],
                    "status": statuses[idx],
                    "duplicate_peak": int(duplicate[idx]),
                    "nearest_deepseek_distance_px": template_to_deepseek[idx],
                    "missed_by_deepseek": int(missed_by_deepseek[idx]),
                    "radius_from_configured_center_nm": float(np.hypot(refined_x[idx] - CX, refined_y[idx] - CY) * PX_NM),
                }
            )

        summary_rows.append(
            {
                "layer": layer_index + 1,
                "configured_center_x_px": CX,
                "configured_center_y_px": CY,
                "configured_radius_nm": R_NM[layer_index],
                "template_sites_inside_radius": len(refined_x),
                "accepted": int(np.sum(accepted)),
                "review": int(np.sum(review)),
                "rejected": int(np.sum(rejected)),
                "deepseek_points_inside_radius": len(deepseek_inside),
                "deepseek_unmatched_extra": int(np.sum(deepseek_extra)),
                "template_sites_missed_by_deepseek": int(np.sum(missed_by_deepseek)),
                "median_refinement_shift_px": float(np.nanmedian(shifts)),
                "median_local_snr": float(np.nanmedian(snr)),
                "accepted_brightness_centroid_x_px": centroid_x,
                "accepted_brightness_centroid_y_px": centroid_y,
                "centroid_offset_from_config_px": float(np.hypot(centroid_x - CX, centroid_y - CY)) if np.isfinite(centroid_x) else np.nan,
            }
        )

        layer_plot_data.append(
            {
                "image": image,
                "radius_px": radius_px,
                "deepseek": deepseek_inside,
                "deepseek_extra": deepseek_extra,
                "x": refined_x,
                "y": refined_y,
                "classes": layer_classes,
                "accepted": accepted,
                "review": review,
                "rejected": rejected,
                "summary": summary_rows[-1],
            }
        )

    with (args.output / "all_layers_coordinate_review.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(detail_rows[0]))
        writer.writeheader()
        writer.writerows(detail_rows)
    with (args.output / "layer_coordinate_review_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    # One detailed review image per layer.
    for layer_index, data in enumerate(layer_plot_data):
        image = data["image"]
        lo, hi = np.percentile(image[280:1200, :1120], [2, 99.6])
        fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.1), facecolor="white")
        titles = ["configured particle geometry", "DeepSeek points currently used", "template-refined coordinates"]
        for ax, title in zip(axes, titles):
            ax.imshow(image, cmap="gray", vmin=lo, vmax=hi, interpolation="none")
            draw_circle(ax, CX, CY, data["radius_px"])
            ax.set_xlim(300, 1120)
            ax.set_ylim(1130, 300)
            ax.set_title(title, fontsize=9)
            clean_axis(ax)

        deepseek = data["deepseek"]
        extra = data["deepseek_extra"]
        axes[1].scatter(deepseek[~extra, 0], deepseek[~extra, 1], s=9, facecolors="none", edgecolors="#F5D547", linewidths=0.6, label="matched")
        axes[1].scatter(deepseek[extra, 0], deepseek[extra, 1], s=12, marker="x", color=EXTRA_COLOR, linewidths=0.7, label="unmatched extra")
        axes[1].legend(loc="lower left", fontsize=6, frameon=False)

        x, y = data["x"], data["y"]
        accepted, review, rejected = data["accepted"], data["review"], data["rejected"]
        classes = data["classes"]
        pt = classes == "Pt-like"
        axes[2].scatter(x[accepted & pt], y[accepted & pt], s=10, facecolors="none", edgecolors=PT_COLOR, linewidths=0.65, label="accepted Pt-like")
        axes[2].scatter(x[accepted & ~pt], y[accepted & ~pt], s=10, facecolors="none", edgecolors=LM_COLOR, linewidths=0.65, label="accepted LM-like")
        axes[2].scatter(x[review], y[review], s=12, facecolors="none", edgecolors=REVIEW_COLOR, linewidths=0.75, label="review")
        axes[2].scatter(x[rejected], y[rejected], s=9, marker="x", color=REJECT_COLOR, linewidths=0.6, label="reject")
        axes[2].legend(loc="lower left", fontsize=5.6, frameon=False, ncol=2)

        summary = data["summary"]
        fig.suptitle(
            f"Layer {layer_index + 1:02d} coordinate audit   R={summary['configured_radius_nm']:.2f} nm   "
            f"template={summary['template_sites_inside_radius']}   accept/review/reject="
            f"{summary['accepted']}/{summary['review']}/{summary['rejected']}   "
            f"DeepSeek extra={summary['deepseek_unmatched_extra']}",
            fontsize=10,
            y=0.985,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(args.output / f"layer_{layer_index + 1:02d}_coordinate_review.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Compact 4x4 overview of the refined coordinates.
    fig, axes = plt.subplots(4, 4, figsize=(12.0, 12.0), facecolor="black")
    for layer_index, (ax, data) in enumerate(zip(axes.flat, layer_plot_data)):
        image = data["image"]
        lo, hi = np.percentile(image[280:1200, :1120], [2, 99.6])
        ax.imshow(image, cmap="gray", vmin=lo, vmax=hi, interpolation="none")
        draw_circle(ax, CX, CY, data["radius_px"])
        x, y = data["x"], data["y"]
        accepted, review = data["accepted"], data["review"]
        classes = data["classes"]
        pt = classes == "Pt-like"
        ax.scatter(x[accepted & pt], y[accepted & pt], s=3.5, color=PT_COLOR, linewidths=0)
        ax.scatter(x[accepted & ~pt], y[accepted & ~pt], s=3.5, color=LM_COLOR, linewidths=0)
        ax.scatter(x[review], y[review], s=5, facecolors="none", edgecolors=REVIEW_COLOR, linewidths=0.45)
        ax.set_xlim(300, 1120)
        ax.set_ylim(1130, 300)
        ax.set_title(
            f"L{layer_index + 1}: {data['summary']['accepted']}/{data['summary']['review']}/{data['summary']['rejected']}",
            color="white",
            fontsize=7,
        )
        clean_axis(ax)
    fig.suptitle("16-layer coordinate review: accepted / review / rejected", color="white", fontsize=12, y=0.995)
    fig.tight_layout(pad=0.6)
    fig.savefig(args.output / "all_layers_coordinate_review_montage.png", dpi=300, bbox_inches="tight", facecolor="black")
    plt.close(fig)

    print(args.output / "all_layers_coordinate_review_montage.png")
    for row in summary_rows:
        print(
            f"L{row['layer']:02d} R={row['configured_radius_nm']:.2f}nm "
            f"template={row['template_sites_inside_radius']} "
            f"A/R/X={row['accepted']}/{row['review']}/{row['rejected']} "
            f"DS={row['deepseek_points_inside_radius']} extra={row['deepseek_unmatched_extra']} "
            f"miss={row['template_sites_missed_by_deepseek']} shift={row['median_refinement_shift_px']:.2f}px "
            f"centroid_offset={row['centroid_offset_from_config_px']:.1f}px"
        )


if __name__ == "__main__":
    main()
