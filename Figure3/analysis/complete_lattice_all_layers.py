"""Complete the full atom lattice in every spherical slice before 3D rendering.

Strong sites use the local layer intensity peak. Weak sites are retained at a
structure-predicted lattice position rather than discarded. Atom class is
inherited from manual labels where available and otherwise inferred from the
user-guided spatial region plus ordered-lattice parity.
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
from sklearn.cluster import KMeans

from audit_layer_coordinates import refine_site, status_from_metrics


PT_COLOR = "#D85C41"
LM_COLOR = "#3E78B2"
WEAK_COLOR = "#FFD166"


def fit_lattice(points: np.ndarray, center: tuple[float, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    tree = cKDTree(points)
    distance, index = tree.query(points, k=9)
    vectors = []
    for atom_index in range(len(points)):
        for dist, neighbor in zip(distance[atom_index, 1:], index[atom_index, 1:]):
            if 14 <= dist <= 26:
                vector = points[neighbor] - points[atom_index]
                if vector[1] < 0 or (abs(vector[1]) < 1e-9 and vector[0] < 0):
                    vector = -vector
                vectors.append(vector)
    vectors = np.asarray(vectors)
    angles = np.arctan2(vectors[:, 1], vectors[:, 0])
    features = np.c_[np.cos(2 * angles), np.sin(2 * angles)]
    clusters = KMeans(n_clusters=3, random_state=0, n_init=20).fit(features).labels_
    median_vectors = np.asarray([np.median(vectors[clusters == group], axis=0) for group in range(3)])
    # The two oblique directions have the largest vertical components; their
    # difference gives the third, near-horizontal nearest-neighbour direction.
    chosen = np.argsort(np.abs(median_vectors[:, 1]))[-2:]
    basis = median_vectors[chosen].T
    if np.linalg.det(basis) < 0:
        basis = basis[:, ::-1]

    origin = points[np.argmin(np.hypot(points[:, 0] - center[0], points[:, 1] - center[1]))].copy()
    for _ in range(8):
        ij = np.rint(np.linalg.solve(basis, points.T - origin[:, None]).T).astype(int)
        design = np.c_[np.ones(len(points)), ij]
        fit_x = np.linalg.lstsq(design, points[:, 0], rcond=None)[0]
        fit_y = np.linalg.lstsq(design, points[:, 1], rcond=None)[0]
        origin = np.asarray([fit_x[0], fit_y[0]])
        basis = np.asarray([[fit_x[1], fit_x[2]], [fit_y[1], fit_y[2]]])
    ij = np.rint(np.linalg.solve(basis, points.T - origin[:, None]).T).astype(int)
    predicted = origin + ij @ basis.T
    residual = np.linalg.norm(points - predicted, axis=1)
    return origin, basis, ij, residual


def lattice_nodes(origin: np.ndarray, basis: np.ndarray, cx: float, cy: float, radius_px: float) -> tuple[np.ndarray, np.ndarray]:
    span = int(np.ceil(radius_px / max(np.min(np.linalg.norm(basis, axis=0)), 1.0))) + 8
    ii, jj = np.meshgrid(np.arange(-span, span + 1), np.arange(-span, span + 1), indexing="ij")
    ij = np.c_[ii.ravel(), jj.ravel()]
    xy = origin + ij @ basis.T
    keep = np.hypot(xy[:, 0] - cx, xy[:, 1] - cy) <= radius_px
    return ij[keep], xy[keep]


def draw_circle(ax: plt.Axes, cx: float, cy: float, radius_px: float) -> None:
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(cx + radius_px * np.cos(theta), cy + radius_px * np.sin(theta), color="red", lw=0.9, ls="--")
    ax.plot(cx, cy, marker="+", ms=6, mew=0.9, color="red")


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

    template = pd.read_csv(args.template)
    template_xy = template[["x_px", "y_px"]].to_numpy(float)
    template_pt = (template["manual_atom_class"].to_numpy(str) == "Pt-like")
    template_region = template["spatial_region"].to_numpy(str)
    template_tree = cKDTree(template_xy)

    middle_mask = np.hypot(template_xy[:, 0] - CX, template_xy[:, 1] - CY) * PX_NM <= R_NM[7]
    origin, basis, template_ij_middle, residual = fit_lattice(template_xy[middle_mask], (CX, CY))
    all_template_ij = np.rint(np.linalg.solve(basis, template_xy.T - origin[:, None]).T).astype(int)
    all_template_pred = origin + all_template_ij @ basis.T
    all_template_residual = np.linalg.norm(template_xy - all_template_pred, axis=1)

    parity = (all_template_ij[:, 0] + all_template_ij[:, 1]) & 1
    accuracy_if_one_is_pt = np.mean((parity == 1) == template_pt)
    pt_parity = 1 if accuracy_if_one_is_pt >= 0.5 else 0
    parity_accuracy = max(accuracy_if_one_is_pt, 1.0 - accuracy_if_one_is_pt)

    # Keep the closest manual annotation for every occupied lattice key.
    manual_by_key = {}
    for index, key in enumerate(map(tuple, all_template_ij)):
        if key not in manual_by_key or all_template_residual[index] < all_template_residual[manual_by_key[key]]:
            manual_by_key[key] = index

    stack = np.asarray(np.load(args.workspace / "output" / "registered_stack.npz")["registered"], dtype=float)
    all_rows = []
    summary_rows = []
    plot_layers = []

    for layer_index in range(16):
        radius_px = float(R_NM[layer_index] / PX_NM)
        node_ij, ideal_xy = lattice_nodes(origin, basis, CX, CY, radius_px)
        image = stack[layer_index]
        enhanced = gaussian_filter(image, 1.1) - gaussian_filter(image, 9.0)
        refinement = [refine_site(image, enhanced, px, py) for px, py in ideal_xy]
        refined_xy = np.asarray([[row["x"], row["y"]] for row in refinement])
        shift = np.asarray([row["shift"] for row in refinement])
        snr = np.asarray([row["snr"] for row in refinement])
        prominence = np.asarray([row["prominence"] for row in refinement])
        status = np.asarray([status_from_metrics(one_shift, one_snr) for one_shift, one_snr in zip(shift, snr)])

        strong = status == "accept"
        weak_peak = status == "review"
        structure_only = status == "reject"

        if np.any(strong):
            translation = np.median(refined_xy[strong] - ideal_xy[strong], axis=0)
        else:
            translation = np.zeros(2)
        completed_xy = ideal_xy + translation
        completed_xy[strong | weak_peak] = refined_xy[strong | weak_peak]
        position_source = np.full(len(node_ij), "lattice_completed", dtype=object)
        position_source[weak_peak] = "weak_peak_refined"
        position_source[strong] = "strong_peak_refined"

        node_pt = np.empty(len(node_ij), dtype=bool)
        class_source = np.empty(len(node_ij), dtype=object)
        nearest_distance, nearest_index = template_tree.query(ideal_xy)
        for node_index, key in enumerate(map(tuple, node_ij)):
            if key in manual_by_key and all_template_residual[manual_by_key[key]] <= 5.0:
                manual_index = manual_by_key[key]
                node_pt[node_index] = template_pt[manual_index]
                class_source[node_index] = "manual_label"
                continue
            region = template_region[nearest_index[node_index]]
            if region == "Pt-rich":
                node_pt[node_index] = True
                class_source[node_index] = "Pt-rich_region"
            elif region == "liquid-metal":
                node_pt[node_index] = False
                class_source[node_index] = "liquid_region"
            else:
                node_pt[node_index] = ((node_ij[node_index, 0] + node_ij[node_index, 1]) & 1) == pt_parity
                class_source[node_index] = "ordered_lattice_parity"

        for node_index in range(len(node_ij)):
            all_rows.append(
                {
                    "layer": layer_index + 1,
                    "lattice_i": int(node_ij[node_index, 0]),
                    "lattice_j": int(node_ij[node_index, 1]),
                    "x_ideal_px": ideal_xy[node_index, 0],
                    "y_ideal_px": ideal_xy[node_index, 1],
                    "x_completed_px": completed_xy[node_index, 0],
                    "y_completed_px": completed_xy[node_index, 1],
                    "atom_class": "Pt-like" if node_pt[node_index] else "LM-like",
                    "class_source": class_source[node_index],
                    "position_source": position_source[node_index],
                    "local_snr": snr[node_index],
                    "local_prominence": prominence[node_index],
                    "refinement_shift_px": shift[node_index],
                    "configured_radius_nm": R_NM[layer_index],
                }
            )

        summary_rows.append(
            {
                "layer": layer_index + 1,
                "configured_radius_nm": R_NM[layer_index],
                "complete_lattice_nodes": len(node_ij),
                "strong_peak_refined": int(np.sum(strong)),
                "weak_peak_refined": int(np.sum(weak_peak)),
                "lattice_completed": int(np.sum(structure_only)),
                "Pt_like": int(np.sum(node_pt)),
                "LM_like": int(np.sum(~node_pt)),
                "layer_translation_x_px": translation[0],
                "layer_translation_y_px": translation[1],
                "median_snr": float(np.nanmedian(snr)),
            }
        )
        plot_layers.append(
            {
                "image": image,
                "radius_px": radius_px,
                "xy": completed_xy,
                "pt": node_pt,
                "strong": strong,
                "weak_peak": weak_peak,
                "structure_only": structure_only,
            }
        )

    with (args.output / "complete_lattice_coordinates_all_layers.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    with (args.output / "complete_lattice_layer_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    np.savez_compressed(
        args.output / "complete_lattice_model.npz",
        origin=origin,
        basis=basis,
        pt_parity=pt_parity,
        parity_accuracy=parity_accuracy,
        layer8_fit_residual=residual,
    )

    # Per-layer two-panel review: raw geometry and completed lattice.
    for layer_index, data in enumerate(plot_layers):
        image = data["image"]
        lo, hi = np.percentile(image[280:1200, :1120], [2, 99.6])
        fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.2), facecolor="white")
        for ax in axes:
            ax.imshow(image, cmap="gray", vmin=lo, vmax=hi, interpolation="none")
            draw_circle(ax, CX, CY, data["radius_px"])
            ax.set_xlim(300, 1120)
            ax.set_ylim(1130, 300)
            clean_axis(ax)
        axes[0].set_title("registered layer + configured boundary")
        axes[1].set_title("complete lattice: filled=peak, hollow=structure-completed")

        xy, pt = data["xy"], data["pt"]
        measured = data["strong"] | data["weak_peak"]
        completed = data["structure_only"]
        axes[1].scatter(xy[measured & pt, 0], xy[measured & pt, 1], s=9, color=PT_COLOR, linewidths=0, label="Pt-like peak")
        axes[1].scatter(xy[measured & ~pt, 0], xy[measured & ~pt, 1], s=9, color=LM_COLOR, linewidths=0, label="LM-like peak")
        axes[1].scatter(xy[completed & pt, 0], xy[completed & pt, 1], s=11, facecolors="none", edgecolors=PT_COLOR, linewidths=0.65, label="Pt-like completed")
        axes[1].scatter(xy[completed & ~pt, 0], xy[completed & ~pt, 1], s=11, facecolors="none", edgecolors=LM_COLOR, linewidths=0.65, label="LM-like completed")
        axes[1].legend(loc="lower left", fontsize=5.7, frameon=False, ncol=2)
        row = summary_rows[layer_index]
        fig.suptitle(
            f"Layer {layer_index + 1:02d} complete lattice   nodes={row['complete_lattice_nodes']}   "
            f"strong/weak/completed={row['strong_peak_refined']}/{row['weak_peak_refined']}/{row['lattice_completed']}",
            fontsize=10,
            y=0.98,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(args.output / f"layer_{layer_index + 1:02d}_complete_lattice.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # 4x4 discussion montage. All structural nodes are present; hollow nodes are inferred.
    fig, axes = plt.subplots(4, 4, figsize=(12.0, 12.0), facecolor="black")
    for layer_index, (ax, data) in enumerate(zip(axes.flat, plot_layers)):
        image = data["image"]
        lo, hi = np.percentile(image[280:1200, :1120], [2, 99.6])
        ax.imshow(image, cmap="gray", vmin=lo, vmax=hi, interpolation="none")
        draw_circle(ax, CX, CY, data["radius_px"])
        xy, pt = data["xy"], data["pt"]
        measured = data["strong"] | data["weak_peak"]
        completed = data["structure_only"]
        ax.scatter(xy[measured & pt, 0], xy[measured & pt, 1], s=3.4, color=PT_COLOR, linewidths=0)
        ax.scatter(xy[measured & ~pt, 0], xy[measured & ~pt, 1], s=3.4, color=LM_COLOR, linewidths=0)
        ax.scatter(xy[completed & pt, 0], xy[completed & pt, 1], s=4.5, facecolors="none", edgecolors=PT_COLOR, linewidths=0.35)
        ax.scatter(xy[completed & ~pt, 0], xy[completed & ~pt, 1], s=4.5, facecolors="none", edgecolors=LM_COLOR, linewidths=0.35)
        ax.set_xlim(300, 1120)
        ax.set_ylim(1130, 300)
        clean_axis(ax)
        row = summary_rows[layer_index]
        ax.set_title(
            f"L{layer_index + 1}: N={row['complete_lattice_nodes']}  peak={row['strong_peak_refined'] + row['weak_peak_refined']}  completed={row['lattice_completed']}",
            color="white",
            fontsize=7,
        )
    fig.suptitle("Complete 16-layer lattice: filled = image-supported, hollow = structure-completed", color="white", fontsize=12, y=0.995)
    fig.tight_layout(pad=0.6)
    fig.savefig(args.output / "complete_lattice_all_layers_montage.png", dpi=300, bbox_inches="tight", facecolor="black")
    plt.close(fig)

    print("origin", origin.tolist())
    print("basis", basis.tolist())
    print("fit_residual_median_px", float(np.median(residual)))
    print("fit_residual_p95_px", float(np.percentile(residual, 95)))
    print("parity_accuracy", float(parity_accuracy))
    for row in summary_rows:
        print(
            f"L{row['layer']:02d} N={row['complete_lattice_nodes']} "
            f"strong/weak/completed={row['strong_peak_refined']}/{row['weak_peak_refined']}/{row['lattice_completed']} "
            f"Pt/LM={row['Pt_like']}/{row['LM_like']}"
        )
    print(args.output / "complete_lattice_all_layers_montage.png")


if __name__ == "__main__":
    main()
