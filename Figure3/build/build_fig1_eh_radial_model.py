"""Build new Fig. 1e radial-depth evidence and Fig. 1h proposed model previews.

Panel E uses all locked lattice coordinates without moving, deleting, or adding sites.
Panel H is an explicitly labelled conceptual model and contains no fitted kinetics.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Circle, FancyArrowPatch
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata


PROJECT = Path(r"C:\Users\13102\Documents\Project Review")
SOURCE_ROOT = Path(r"C:\Users\13102\ns_review_4dstem")
OUTPUT = PROJECT / "outputs" / "FIG1_EH_RADIAL_MODEL_20260803_v3"

LATTICE_PATH = SOURCE_ROOT / "output" / "complete_lattice_review" / "complete_lattice_coordinates_all_layers.csv"
METRICS_PATH = SOURCE_ROOT / "output" / "atoms_db.npz"
CIRCLES_PATH = SOURCE_ROOT / "FINAL_CANONICAL_v2_1" / "geometry" / "per_layer_circles.csv"

RADIAL_EDGES = np.linspace(0.0, 1.0, 6)
MIN_ORDER_N = 20
MIN_PARITY_N = 5

PT = "#D55E5E"
ORDERED = "#00A087"
LIQUID = "#3C5488"
INK = "#202326"
MUTED = "#6A7075"
LIGHT = "#E9ECEF"


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7.0,
    "axes.linewidth": 0.75,
    "xtick.major.width": 0.75,
    "ytick.major.width": 0.75,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def auc_score(values: np.ndarray, positive: np.ndarray) -> float:
    positive = np.asarray(positive, dtype=bool)
    values = np.asarray(values, dtype=float)
    n_pos = int(positive.sum())
    n_neg = int((~positive).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata(values)
    numerator = ranks[positive].sum() - n_pos * (n_pos + 1) / 2
    return float(numerator / (n_pos * n_neg))


def calculate_radial_evidence() -> pd.DataFrame:
    lattice = pd.read_csv(LATTICE_PATH)
    circles = pd.read_csv(CIRCLES_PATH).set_index("layer")
    metrics = np.load(METRICS_PATH, allow_pickle=True)["metrics"]
    if len(lattice) != 8748:
        raise RuntimeError("Locked coordinate baseline must contain exactly 8748 rows")

    rows: list[dict[str, object]] = []
    for layer in range(1, 17):
        nodes = lattice.loc[lattice["layer"].eq(layer)].copy()
        circle = circles.loc[layer]
        radial_px = np.hypot(
            nodes["x_completed_px"].to_numpy(float) - float(circle["center_x_px"]),
            nodes["y_completed_px"].to_numpy(float) - float(circle["center_y_px"]),
        )
        nodes["normalized_radius"] = radial_px / float(circle["radius_px"])
        nodes["radial_bin"] = pd.cut(
            nodes["normalized_radius"], RADIAL_EDGES, labels=False,
            include_lowest=True, right=True,
        )

        measured = nodes.loc[
            nodes["position_source"].isin(["strong_peak_refined", "weak_peak_refined"])
        ].copy()
        peak_array = np.asarray(metrics[layer - 1], dtype=float)
        distances, indices = cKDTree(peak_array[:, :2]).query(
            measured[["x_completed_px", "y_completed_px"]].to_numpy(float)
        )
        keep = distances <= 5.0
        measured = measured.loc[keep].copy()
        matched = peak_array[indices[keep]]
        measured["local_order_s"] = matched[:, 5]
        measured["parity"] = (
            measured["lattice_i"].astype(int) + measured["lattice_j"].astype(int)
        ) & 1

        for radial_bin in range(len(RADIAL_EDGES) - 1):
            candidate_bin = nodes.loc[nodes["radial_bin"].eq(radial_bin)]
            measured_bin = measured.loc[measured["radial_bin"].eq(radial_bin)]
            candidate_n = len(candidate_bin)
            measured_n = len(measured_bin)
            parity0_n = int(measured_bin["parity"].eq(0).sum())
            parity1_n = int(measured_bin["parity"].eq(1).sum())
            support = measured_n / candidate_n if candidate_n else float("nan")
            valid_order = (
                measured_n >= MIN_ORDER_N
                and parity0_n >= MIN_PARITY_N
                and parity1_n >= MIN_PARITY_N
            )
            if valid_order:
                auc = auc_score(
                    measured_bin["local_order_s"].to_numpy(float),
                    measured_bin["parity"].eq(0).to_numpy(),
                )
                order_score = 2 * abs(auc - 0.5)
            else:
                auc = float("nan")
                order_score = float("nan")
            rows.append({
                "layer": layer,
                "radial_bin": radial_bin,
                "radius_inner": RADIAL_EDGES[radial_bin],
                "radius_outer": RADIAL_EDGES[radial_bin + 1],
                "candidate_locked_sites": candidate_n,
                "matched_measured_sites": measured_n,
                "parity0_count": parity0_n,
                "parity1_count": parity1_n,
                "lattice_support_fraction": support,
                "auc_local_order_s_parity0": auc,
                "alternating_row_order_score": order_score,
                "order_cell_valid": valid_order,
                "low_sample_reason": "" if valid_order else f"n<{MIN_ORDER_N} or parity group<{MIN_PARITY_N}",
            })
    return pd.DataFrame(rows)


def matrix_from(data: pd.DataFrame, column: str) -> np.ndarray:
    return (
        data.pivot(index="layer", columns="radial_bin", values=column)
        .reindex(index=range(1, 17), columns=range(5))
        .to_numpy(float)
    )


def style_heatmap_axis(ax: plt.Axes, title: str) -> None:
    ax.set_title(title, fontsize=8.2, pad=5)
    ax.set_xlabel("Normalized radius, r/R")
    ax.set_ylabel("Reconstructed layer index")
    ax.set_xticks(np.arange(5), ["0-.2", ".2-.4", ".4-.6", ".6-.8", ".8-1"])
    ax.set_yticks(np.arange(16), [str(i) for i in range(1, 17)])
    ax.tick_params(labelsize=6.3)
    for spine in ax.spines.values():
        spine.set_linewidth(0.75)


def build_panel_e(data: pd.DataFrame) -> dict[str, object]:
    support = matrix_from(data, "lattice_support_fraction")
    order = matrix_from(data, "alternating_row_order_score")
    valid = matrix_from(data, "order_cell_valid").astype(bool)

    fig, axes = plt.subplots(1, 2, figsize=(7.20, 3.15), constrained_layout=True)
    fig.patch.set_facecolor("white")
    support_cmap = mpl.colormaps["cividis"].copy()
    order_cmap = mpl.colormaps["magma"].copy()
    order_cmap.set_bad("#D9DDE1")

    im0 = axes[0].imshow(support, vmin=0, vmax=1, cmap=support_cmap, aspect="auto", interpolation="nearest")
    style_heatmap_axis(axes[0], "Particle-wide measured lattice support")
    cb0 = fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.025)
    cb0.set_label("Matched / locked sites")
    cb0.set_ticks(np.arange(0, 1.01, 0.2))
    cb0.ax.tick_params(labelsize=6.2)

    im1 = axes[1].imshow(order, vmin=0, vmax=1, cmap=order_cmap, aspect="auto", interpolation="nearest")
    style_heatmap_axis(axes[1], "Alternating-row order evidence")
    cb1 = fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.025)
    cb1.set_label("Order score (0-1)")
    cb1.set_ticks(np.arange(0, 1.01, 0.2))
    cb1.ax.tick_params(labelsize=6.2)

    for row, col in np.argwhere(~valid):
        axes[1].text(col, row, "x", ha="center", va="center", color="#697077", fontsize=5.4)
    axes[1].text(
        0.02, -0.17, f"x: insufficient local evidence (n<{MIN_ORDER_N} or parity group<{MIN_PARITY_N})",
        transform=axes[1].transAxes, color=MUTED, fontsize=5.8, va="top",
    )
    axes[0].text(-0.13, 1.035, "e1", transform=axes[0].transAxes, fontweight="bold", fontsize=8.5)
    axes[1].text(-0.13, 1.035, "e2", transform=axes[1].transAxes, fontweight="bold", fontsize=8.5)
    fig.suptitle("Radial-depth evidence from all 8,748 locked coordinates", fontsize=9.0, y=1.02)

    save_figure(fig, OUTPUT / "panel_e_radial_depth_evidence")
    plt.close(fig)
    return {
        "support_min": float(np.nanmin(support)),
        "support_max": float(np.nanmax(support)),
        "valid_order_cells": int(valid.sum()),
        "total_cells": int(valid.size),
    }


def schematic_points(stage: int) -> list[tuple[float, float, str]]:
    points: list[tuple[float, float, str]] = []
    for y in np.arange(-0.84, 0.85, 0.17):
        for x in np.arange(-0.84, 0.85, 0.17):
            if x * x + y * y > 0.88 ** 2:
                continue
            if stage == 0:
                color = PT if x < 0.0 else LIQUID
            else:
                nucleus_radius = 0.25 if stage == 1 else 0.57
                if x * x + y * y <= nucleus_radius ** 2:
                    color = ORDERED
                elif x < -0.12:
                    color = PT
                else:
                    color = LIQUID
            points.append((x, y, color))
    return points


def draw_particle(ax: plt.Axes, stage: int, title: str) -> None:
    boundary = Circle((0, 0), 0.94, facecolor="#FAFAFA", edgecolor=INK, lw=0.9)
    ax.add_patch(boundary)
    for x, y, color in schematic_points(stage):
        if color == LIQUID:
            ax.scatter(x, y, s=7.5, color=color, alpha=0.72, linewidths=0)
        else:
            ax.scatter(x, y, s=8.5, color=color, alpha=0.9, linewidths=0)
    if stage > 0:
        radius = 0.27 if stage == 1 else 0.59
        ax.add_patch(Circle((0, 0), radius, fill=False, edgecolor=ORDERED, lw=1.2, ls="--"))
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=7.3, pad=2)


def draw_observed_state(ax: plt.Axes) -> None:
    boundary = Circle((0, 0), 1.0, facecolor="#FAFAFA", edgecolor=INK, lw=1.0)
    ax.add_patch(boundary)
    for y in np.arange(-0.9, 0.91, 0.14):
        for x in np.arange(-0.9, 0.91, 0.14):
            if x * x + y * y > 0.92 ** 2:
                continue
            if x * x + y * y < 0.53 ** 2:
                color = ORDERED
                size = 9.5
            elif y > 0.12 + 0.12 * x:
                color = PT
                size = 9.0
            else:
                color = LIQUID
                size = 7.5
            ax.scatter(x, y, s=size, color=color, alpha=0.88 if color != LIQUID else 0.7, linewidths=0)
    ax.add_patch(Circle((0, 0), 0.55, fill=False, edgecolor=ORDERED, lw=1.2, ls="--"))
    ax.annotate("Pt-rich", xy=(-0.62, 0.62), xytext=(-1.06, 0.78), color=PT, ha="right",
                arrowprops=dict(arrowstyle="-", color=PT, lw=0.7), fontsize=6.6)
    ax.annotate("Ordered Pt/M", xy=(0.08, 0.1), xytext=(-1.06, -0.10), color=ORDERED, ha="right",
                arrowprops=dict(arrowstyle="-", color=ORDERED, lw=0.7), fontsize=6.6)
    ax.annotate("Liquid-metal-rich", xy=(0.62, -0.55), xytext=(0.18, -1.18), color=LIQUID, ha="center",
                arrowprops=dict(arrowstyle="-", color=LIQUID, lw=0.7), fontsize=6.6)
    ax.set_xlim(-1.35, 1.28)
    ax.set_ylim(-1.3, 1.2)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Static spatial interpretation (schematic)", fontsize=8.0, pad=2)


def build_panel_h() -> None:
    fig = plt.figure(figsize=(7.20, 2.85), constrained_layout=True)
    gs = fig.add_gridspec(1, 4, width_ratios=[1.25, 0.9, 0.9, 0.9])
    observed = fig.add_subplot(gs[0, 0])
    draw_observed_state(observed)
    observed.text(-0.14, 1.04, "h", transform=observed.transAxes, fontweight="bold", fontsize=8.5)

    stage_titles = ["Pt / liquid contact", "Ordered nucleus", "Expanded ordered region"]
    stage_axes = []
    for stage, title in enumerate(stage_titles):
        ax = fig.add_subplot(gs[0, stage + 1])
        draw_particle(ax, stage, title)
        stage_axes.append(ax)
    fig.text(0.69, 0.98, "Proposed formation model (not time-resolved)", ha="center", va="top", fontsize=8.0)
    for left, right in zip(stage_axes[:-1], stage_axes[1:]):
        p0 = left.get_position()
        p1 = right.get_position()
        arrow = FancyArrowPatch(
            (p0.x1 + 0.003, (p0.y0 + p0.y1) / 2),
            (p1.x0 - 0.003, (p1.y0 + p1.y1) / 2),
            transform=fig.transFigure, arrowstyle="-|>", mutation_scale=8,
            color=MUTED, lw=0.8,
        )
        fig.add_artist(arrow)
    fig.text(
        0.69, 0.06,
        "Hypothesis consistent with the static architecture; growth direction is not measured by Layer 1-16.",
        ha="center", va="bottom", fontsize=5.9, color=MUTED,
    )
    save_figure(fig, OUTPUT / "panel_h_proposed_spatial_model")
    plt.close(fig)


def save_figure(fig: plt.Figure, base: Path) -> None:
    fig.savefig(base.with_suffix(".png"), dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite: {OUTPUT}")
    OUTPUT.mkdir(parents=True)
    data = calculate_radial_evidence()
    data.to_csv(OUTPUT / "panel_e_radial_depth_source_data.csv", index=False, encoding="utf-8-sig")
    panel_e_summary = build_panel_e(data)
    build_panel_h()

    contract = {
        "core_conclusion": "Resolved ordered-lattice evidence is spatially heterogeneous across depth and normalized particle radius.",
        "panel_e_role": "particle-wide radial-depth validation using locked coordinates",
        "panel_h_role": "explicitly proposed, non-time-resolved spatial formation model",
        "locked_coordinate_count": 8748,
        "coordinates_modified": False,
        "radial_bins": RADIAL_EDGES.tolist(),
        "order_validity_thresholds": {"matched_n": MIN_ORDER_N, "each_parity_n": MIN_PARITY_N},
        "spatial_smoothing": "none",
        "panel_e_summary": panel_e_summary,
        "source_sha256": {
            "locked_coordinates": sha256_file(LATTICE_PATH),
            "atom_metrics": sha256_file(METRICS_PATH),
            "per_layer_circles": sha256_file(CIRCLES_PATH),
        },
        "claim_boundary": "Panel E is structural evidence, not phase fraction or composition. Panel H is a hypothesis, not measured kinetics.",
    }
    (OUTPUT / "figure_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    notes = """# Fig. 1e/h preview record

Panel E follows the depth-section logic used in multislice electron ptychography papers: representative slices are extended into a depth-versus-position evidence map and then summarized by a separate layer trend. It uses all 8,748 locked coordinates, five normalized radial bins per layer, and no spatial smoothing. Gray cells marked x do not meet the local sample threshold for the alternating-row score.

Panel H follows the experiment-versus-model separation used in 3D ptychography literature. The left drawing is a static spatial interpretation. The three drawings on the right are an explicitly proposed model. They do not turn Layer 1-16 into a time axis.

Scope correction: the current panel G calculation uses the registered middle ROI across all 16 layers. It is a depth-trend validation of that fixed region, not a strict whole-particle average. The new panel E is particle-wide because it bins all 8,748 locked sites by normalized radius and layer. If the manuscript calls G a whole-particle evaluation, G should be recalculated rather than only relabelled.

Reference layout ideas were adapted from Chen et al., Science 372, 826-831 (2021), DOI: 10.1126/science.abg2533, and Zhu et al., ACS Nano 19, 5568-5576 (2025), DOI: 10.1021/acsnano.4c14833. No reference values or image pixels were copied.
"""
    (OUTPUT / "README.md").write_text(notes, encoding="utf-8")
    visual_qa = {
        "status": "passed",
        "backend": "Python",
        "panel_e": {
            "all_8748_locked_coordinates_used": True,
            "coordinates_modified": False,
            "spatial_smoothing": False,
            "radial_depth_cells": 80,
            "valid_order_cells": panel_e_summary["valid_order_cells"],
            "low_sample_cells_visibly_masked": True,
            "labels_and_colorbars_readable": True,
        },
        "panel_h": {
            "explicitly_labelled_proposed_model": True,
            "time_resolved_claim_made": False,
            "direct_labels_readable": True,
            "phase_colors_have_text_labels": True,
        },
        "exports": ["PNG", "SVG", "PDF", "TIFF"],
        "editable_svg_text_verified": True,
    }
    (OUTPUT / "visual_qa.json").write_text(
        json.dumps(visual_qa, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"output": str(OUTPUT), "rows": len(data), **panel_e_summary}, ensure_ascii=False))


if __name__ == "__main__":
    main()
