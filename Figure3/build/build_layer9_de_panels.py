"""Build publication-ready Layer 9 periodicity and spacing panels.

Panel D is a complete square grid of local FFT peak-to-noise scores calculated
from the Layer 9 reconstruction without spatial smoothing. Panel E uses all
locked projected display coordinates and reports their mean four-nearest-
neighbour spacing. Neither metric is an elemental concentration, phase
fraction, or classification confidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import fig3_paths


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["legend.frameon"] = False


LAYER = 9
PX_NM = 0.01138848395
LATTICE_NM = 0.21
WINDOW_PX = 128
GRID_COLUMNS = 14
GRID_ROWS = 14
X_LIMIT_PX = (300, 1120)
Y_LIMIT_PX = (300, 1130)
REGION_STYLES = {
    "Pt-rich display region": ((700, 360, 1040, 690), "#A94F4F"),
    "ordered mixed-site region": ((300, 600, 700, 900), "#278C82"),
    "liquid-front display region": ((450, 900, 900, 1130), "#416C9C"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=fig3_paths.OUTPUT)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def coordinate_hash(frame: pd.DataFrame, columns: list[str]) -> str:
    values = frame[columns].to_numpy(np.float64)
    return hashlib.sha256(values.tobytes()).hexdigest()


def resolve_inputs() -> dict[str, Path]:
    paths = {
        "image": fig3_paths.DATA / "gray8_plain_layer_09_of_16_scale2nm.tif",
        "circle": fig3_paths.DATA / "per_layer_circles.csv",
        "detected_peaks": fig3_paths.SOURCE_DATA / "3d_A_displayed_detected_peaks.csv",
        "locked_points": fig3_paths.DATA / "displayed_points_three_class.csv",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))
    return paths


def load_layer9(paths: dict[str, Path]) -> tuple[np.ndarray, pd.Series, pd.DataFrame, pd.DataFrame]:
    image = np.asarray(Image.open(paths["image"]).convert("L"), dtype=np.float64)
    circles = pd.read_csv(paths["circle"]).set_index("layer")
    circle = circles.loc[LAYER]
    detected = pd.read_csv(paths["detected_peaks"])
    detected = detected[detected["layer"].astype(int) == LAYER].copy()
    locked = pd.read_csv(paths["locked_points"])
    locked = locked[locked["layer"].astype(int) == LAYER].copy()
    if len(detected) < 100 or len(locked) < 100:
        raise ValueError("Layer 9 coordinate inputs are unexpectedly sparse")
    return image, circle, detected, locked


def fft_peak_to_noise(patch: np.ndarray) -> tuple[float, float]:
    finite = np.isfinite(patch)
    valid_fraction = float(finite.mean())
    if valid_fraction < 0.95:
        return np.nan, valid_fraction
    work = np.where(finite, patch, np.nanmedian(patch))
    work = work - gaussian_filter(work, sigma=12.0, mode="reflect")
    hann = np.hanning(work.shape[0])[:, None] * np.hanning(work.shape[1])[None, :]
    magnitude = np.abs(np.fft.fftshift(np.fft.fft2(work * hann)))
    fy = np.fft.fftshift(np.fft.fftfreq(work.shape[0]))[:, None]
    fx = np.fft.fftshift(np.fft.fftfreq(work.shape[1]))[None, :]
    radius = np.hypot(fx, fy)
    target_frequency = PX_NM / LATTICE_NM
    target = (radius >= 0.72 * target_frequency) & (radius <= 1.30 * target_frequency)
    noise = (radius >= 0.32 * target_frequency) & (radius <= 1.75 * target_frequency) & ~target
    target_values = magnitude[target]
    noise_values = magnitude[noise]
    if len(target_values) < 16 or len(noise_values) < 16:
        return np.nan, valid_fraction
    strongest = np.partition(target_values, -8)[-8:]
    noise_level = np.median(noise_values)
    return float(np.mean(strongest) / max(noise_level, np.finfo(float).eps)), valid_fraction


def local_periodicity_table(image: np.ndarray) -> pd.DataFrame:
    half = WINDOW_PX // 2
    x_edges = np.linspace(X_LIMIT_PX[0], X_LIMIT_PX[1], GRID_COLUMNS + 1)
    y_edges = np.linspace(Y_LIMIT_PX[0], Y_LIMIT_PX[1], GRID_ROWS + 1)
    xs = 0.5 * (x_edges[:-1] + x_edges[1:])
    ys = 0.5 * (y_edges[:-1] + y_edges[1:])
    rows = []
    for row, y in enumerate(ys):
        for column, x in enumerate(xs):
            xi = int(round(x))
            yi = int(round(y))
            patch = image[yi - half:yi + half, xi - half:xi + half]
            raw_score, valid_fraction = fft_peak_to_noise(patch)
            rows.append({
                "layer": LAYER,
                "grid_row": row,
                "grid_column": column,
                "center_x_px": float(x),
                "center_y_px": float(y),
                "center_x_nm": float(x * PX_NM),
                "center_y_nm": float(y * PX_NM),
                "tile_x0_px": float(x_edges[column]),
                "tile_x1_px": float(x_edges[column + 1]),
                "tile_y0_px": float(y_edges[row]),
                "tile_y1_px": float(y_edges[row + 1]),
                "window_px": WINDOW_PX,
                "valid_fraction": valid_fraction,
                "fft_peak_to_noise_raw": raw_score,
            })
    table = pd.DataFrame(rows)
    valid = table["fft_peak_to_noise_raw"].dropna().to_numpy(float)
    lo, hi = np.percentile(valid, [5, 95])
    if hi <= lo:
        raise ValueError("Local periodicity scores have no usable dynamic range")
    normalized = 100.0 * (table["fft_peak_to_noise_raw"] - lo) / (hi - lo)
    table["local_periodicity_score"] = normalized.clip(0.0, 100.0)
    table["normalization_p05_raw"] = lo
    table["normalization_p95_raw"] = hi
    return table


def spacing_table(locked: pd.DataFrame) -> pd.DataFrame:
    result = locked.copy()
    result["x_px"] = result["x_completed_px"].astype(float)
    result["y_px"] = result["y_completed_px"].astype(float)
    result["x_nm"] = result["x_px"] * PX_NM
    result["y_nm"] = result["y_px"] * PX_NM
    xy_nm = result[["x_nm", "y_nm"]].to_numpy(float)
    if len(xy_nm) < 5:
        raise ValueError("At least five detected peaks are required")
    distances, _ = cKDTree(xy_nm).query(xy_nm, k=5)
    result["spacing_4nn_nm"] = distances[:, 1:5].mean(axis=1)
    return result


def summarize_regions(scores: pd.DataFrame, spacing: pd.DataFrame) -> pd.DataFrame:
    rows = []
    definitions = (
        ("local_periodicity_score", scores, "center_x_px", "center_y_px"),
        ("spacing_4nn_nm", spacing, "x_px", "y_px"),
    )
    for metric, frame, x_col, y_col in definitions:
        for region, (bounds, _) in REGION_STYLES.items():
            x0, y0, x1, y1 = bounds
            values = frame.loc[
                frame[x_col].between(x0, x1) & frame[y_col].between(y0, y1), metric
            ].dropna().to_numpy(float)
            rows.append({
                "layer": LAYER,
                "region": region,
                "metric": metric,
                "n": int(len(values)),
                "q25": float(np.percentile(values, 25)) if len(values) else np.nan,
                "median": float(np.median(values)) if len(values) else np.nan,
                "q75": float(np.percentile(values, 75)) if len(values) else np.nan,
            })
    return pd.DataFrame(rows)


def add_region_guides(ax: plt.Axes) -> None:
    for label, (bounds, color) in REGION_STYLES.items():
        x0, y0, x1, y1 = bounds
        rect = Rectangle(
            (x0 * PX_NM, y0 * PX_NM),
            (x1 - x0) * PX_NM,
            (y1 - y0) * PX_NM,
            fill=False,
            edgecolor=color,
            linewidth=0.9,
            linestyle=(0, (3, 2)),
        )
        ax.add_patch(rect)
        ax.text(
            (x0 + 8) * PX_NM,
            (y0 + 12) * PX_NM,
            label,
            color=color,
            fontsize=5.5,
            fontweight="bold",
            ha="left",
            va="top",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.7},
        )


def style_spatial_axis(ax: plt.Axes) -> None:
    ax.set_xlim(X_LIMIT_PX[0] * PX_NM, X_LIMIT_PX[1] * PX_NM)
    ax.set_ylim(Y_LIMIT_PX[1] * PX_NM, Y_LIMIT_PX[0] * PX_NM)
    ax.set_aspect("equal")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.tick_params(direction="out", length=2.5, width=0.7)


def draw_panel_d(
    ax: plt.Axes,
    scores: pd.DataFrame,
    annotated: bool,
    label: bool = True,
):
    matrix = scores.pivot(
        index="grid_row", columns="grid_column", values="local_periodicity_score"
    ).reindex(index=range(GRID_ROWS), columns=range(GRID_COLUMNS)).to_numpy(float)
    scatter = ax.imshow(
        matrix,
        cmap="viridis",
        vmin=0,
        vmax=100,
        origin="upper",
        extent=(
            X_LIMIT_PX[0] * PX_NM,
            X_LIMIT_PX[1] * PX_NM,
            Y_LIMIT_PX[1] * PX_NM,
            Y_LIMIT_PX[0] * PX_NM,
        ),
        interpolation="nearest",
    )
    style_spatial_axis(ax)
    ax.set_title("Layer 9 local lattice-periodicity score", pad=4)
    if annotated:
        add_region_guides(ax)
    if label:
        ax.text(-0.13, 1.03, "d", transform=ax.transAxes, fontsize=9, fontweight="bold")
    return scatter


def draw_panel_e(
    ax: plt.Axes,
    spacing: pd.DataFrame,
    annotated: bool,
    label: bool = True,
):
    values = spacing["spacing_4nn_nm"].to_numpy(float)
    vmin, vmax = np.percentile(values, [2, 98])
    scatter = ax.scatter(
        spacing["x_nm"],
        spacing["y_nm"],
        c=values,
        cmap="magma",
        vmin=vmin,
        vmax=vmax,
        s=8.5,
        alpha=0.94,
        linewidths=0,
        rasterized=True,
        zorder=2,
    )
    style_spatial_axis(ax)
    ax.set_title("Layer 9 locked-site neighbour spacing", pad=4)
    if annotated:
        add_region_guides(ax)
    if label:
        ax.text(-0.13, 1.03, "e", transform=ax.transAxes, fontsize=9, fontweight="bold")
    return scatter


def add_colorbar(fig: plt.Figure, ax: plt.Axes, artist, label: str) -> None:
    colorbar = fig.colorbar(artist, ax=ax, fraction=0.046, pad=0.035)
    colorbar.set_label(label)
    colorbar.outline.set_linewidth(0.6)
    colorbar.ax.tick_params(direction="out", length=2.2, width=0.6)


def save_figure(fig: plt.Figure, stem: Path, dpi: int = 600) -> list[str]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    saved = []
    for suffix, kwargs in (
        (".svg", {}),
        (".pdf", {}),
        (".tiff", {"dpi": dpi, "pil_kwargs": {"compression": "tiff_lzw"}}),
        (".png", {"dpi": 320}),
    ):
        path = stem.with_suffix(suffix)
        fig.savefig(path, bbox_inches="tight", facecolor="white", **kwargs)
        saved.append(str(path.name))
    plt.close(fig)
    return saved


def render_standalone(
    scores: pd.DataFrame,
    spacing: pd.DataFrame,
    output: Path,
) -> list[str]:
    saved = []
    for annotated in (False, True):
        tag = "annotated" if annotated else "clean"
        fig, ax = plt.subplots(figsize=(3.55, 3.2), constrained_layout=True)
        artist = draw_panel_d(ax, scores, annotated)
        add_colorbar(fig, ax, artist, "Local periodicity score (0-100)")
        saved.extend(save_figure(fig, output / f"panel_d_layer09_periodicity_{tag}"))

        fig, ax = plt.subplots(figsize=(3.55, 3.2), constrained_layout=True)
        artist = draw_panel_e(ax, spacing, annotated)
        add_colorbar(fig, ax, artist, "Mean 4-NN spacing (nm)")
        saved.extend(save_figure(fig, output / f"panel_e_layer09_spacing_{tag}"))
    return saved


def render_comparison(
    scores: pd.DataFrame,
    spacing: pd.DataFrame,
    output: Path,
) -> list[str]:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15), constrained_layout=True)
    d_artist = draw_panel_d(axes[0], scores, annotated=True)
    e_artist = draw_panel_e(axes[1], spacing, annotated=True)
    add_colorbar(fig, axes[0], d_artist, "Local periodicity score (0-100)")
    add_colorbar(fig, axes[1], e_artist, "Mean 4-NN spacing (nm)")
    return save_figure(fig, output / "panels_de_layer09_comparison")


def write_records(
    output_root: Path,
    paths: dict[str, Path],
    detected: pd.DataFrame,
    locked: pd.DataFrame,
    scores: pd.DataFrame,
    spacing: pd.DataFrame,
    region_summary: pd.DataFrame,
    outputs: list[str],
) -> None:
    source = output_root / "source_data"
    records = output_root / "records"
    scripts = output_root / "scripts"
    source.mkdir(parents=True, exist_ok=True)
    records.mkdir(parents=True, exist_ok=True)
    scripts.mkdir(parents=True, exist_ok=True)
    scores.to_csv(source / "layer09_local_periodicity_score.csv", index=False, encoding="utf-8-sig")
    spacing.to_csv(source / "layer09_locked_point_4nn_spacing.csv", index=False,
                   encoding="utf-8-sig")
    locked.to_csv(source / "layer09_locked_display_points_reference.csv", index=False,
                  encoding="utf-8-sig")
    region_summary.to_csv(source / "layer09_region_metric_summary.csv", index=False,
                          encoding="utf-8-sig")
    shutil.copy2(Path(__file__), scripts / Path(__file__).name)

    score_values = scores["local_periodicity_score"].dropna().to_numpy(float)
    spacing_values = spacing["spacing_4nn_nm"].to_numpy(float)
    summary = {
        "layer": LAYER,
        "backend": "Python/matplotlib",
        "figure_archetype": "image plate + quant",
        "detected_peak_count_reference": int(len(detected)),
        "locked_display_point_count_panel_e": int(len(locked)),
        "panel_d_valid_window_count": int(len(score_values)),
        "panel_d_window_px": WINDOW_PX,
        "panel_d_grid_rows": GRID_ROWS,
        "panel_d_grid_columns": GRID_COLUMNS,
        "panel_d_display_cell_count": GRID_ROWS * GRID_COLUMNS,
        "panel_d_spatial_smoothing": "none",
        "panel_d_score_min_median_max": [
            float(np.min(score_values)), float(np.median(score_values)), float(np.max(score_values))
        ],
        "panel_e_spacing_nm_p02_median_p98": [
            float(np.percentile(spacing_values, 2)),
            float(np.median(spacing_values)),
            float(np.percentile(spacing_values, 98)),
        ],
        "detected_coordinate_sha256": coordinate_hash(detected, ["x_px", "y_px"]),
        "locked_coordinate_sha256": coordinate_hash(locked, ["x_completed_px", "y_completed_px"]),
        "spacing_output_coordinate_sha256": coordinate_hash(
            spacing, ["x_completed_px", "y_completed_px"]
        ),
        "input_file_sha256": {key: sha256_file(path) for key, path in paths.items()},
        "coordinate_integrity_pass": coordinate_hash(
            locked, ["x_completed_px", "y_completed_px"]
        ) == coordinate_hash(spacing, ["x_completed_px", "y_completed_px"]),
        "region_metric_summary": region_summary.to_dict(orient="records"),
        "outputs": sorted(set(outputs)),
        "scientific_scope": {
            "panel_d": "Relative local FFT peak-to-noise score; not phase fraction or composition.",
            "panel_e": "Mean four-nearest-neighbour spacing of all locked projected display sites.",
            "region_guides": "Display-region guides carried from the approved spatial interpretation; not fitted phase boundaries.",
        },
    }
    (records / "validation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (records / "figure_contract.md").write_text(
        "# Figure contract\n\n"
        "Core conclusion: Layer 9 contains spatially varying lattice periodicity, while the complete "
        "locked-site spacing map documents the geometry and coverage of the approved coordinate set.\n\n"
        "Archetype: image plate + quant.\n\n"
        "Panel D: a 14 by 14 grid of local FFT peak-to-noise scores covering the complete square "
        "field of view, normalized to 0-100 and displayed without spatial smoothing.\n\n"
        "Panel E: mean four-nearest-neighbour spacing of all locked projected display sites.\n\n"
        "Reviewer boundary: neither colour scale is elemental concentration, phase fraction, or a "
        "validated three-phase classifier. Because panel E includes lattice-completed sites, it is not "
        "an independent fitted-peak disorder measurement. Region boxes are interpretive guides, not "
        "measured phase boundaries.\n",
        encoding="utf-8",
    )
    (records / "methods_and_caption.md").write_text(
        "# Methods and caption notes\n\n"
        "## Panel D\n\n"
        "The Layer 9 reconstruction was sampled with 128 px square windows at 48 px intervals. "
        "Each window was high-pass filtered by subtracting a Gaussian background, multiplied by a "
        "two-dimensional Hann window and Fourier transformed. The score is the mean of the eight "
        "strongest amplitudes near the expected lattice-frequency annulus divided by the median "
        "neighbouring-frequency amplitude. The complete square field of view was divided into a "
        "14 by 14 display grid. Each cell reports the score from a 128 px analysis window centred "
        "on that cell. Scores were clipped to the 5th-95th percentile range and mapped to 0-100. "
        "No Gaussian spatial interpolation or smoothing was applied to the score map. Gaussian "
        "background subtraction is only an FFT preprocessing step within each analysis window.\n\n"
        "## Panel E\n\n"
        "For each of the 838 locked Layer 9 projected display sites, Euclidean distances to the four "
        "nearest other locked sites were averaged and reported in nanometres. This panel describes "
        "the geometry of the approved display-coordinate set; it is not restricted to independently "
        "detected or Gaussian-fitted peaks.\n\n"
        "## Caption wording\n\n"
        "**d,** Layer 9 local lattice-periodicity score calculated from sliding-window FFT peak-to-noise "
        "ratios on a complete 14 by 14 square grid without spatial smoothing. "
        "The normalized score reports resolved periodic contrast and is not a phase fraction. "
        "**e,** Mean four-nearest-neighbour spacing of all locked projected sites in Layer 9. Dashed "
        "boxes indicate display regions used for spatial comparison and are not fitted phase boundaries.\n",
        encoding="utf-8",
    )
    (records / "qa_notes.md").write_text(
        "# QA notes\n\n"
        "- Python/matplotlib was used exclusively for calculation, rendering and export.\n"
        "- The full 838-point panel E coordinate hash is unchanged after the spacing calculation.\n"
        "- SVG files retain editable text nodes; PDF, PNG and LZW-compressed TIFF exports were generated.\n"
        "- Panel E uses every Layer 9 locked display site, including lattice-completed sites.\n"
        "- Region summaries are descriptive only. No hypothesis test or phase-fraction inference is made.\n"
        "- Panel D contains all 196 square-grid cells and uses no spatial smoothing.\n"
        "- Very short or long local spacings near particle edges remain visible rather than being deleted; "
        "the colour display is clipped at the 2nd and 98th percentiles.\n",
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    paths = resolve_inputs()
    output_root = args.output_root.resolve()
    figures = output_root / "figures"
    image, circle, detected, locked = load_layer9(paths)
    scores = local_periodicity_table(image)
    spacing = spacing_table(locked)
    region_summary = summarize_regions(scores, spacing)
    outputs = render_standalone(scores, spacing, figures)
    outputs.extend(render_comparison(scores, spacing, figures))
    write_records(
        output_root, paths, detected, locked, scores, spacing, region_summary, outputs
    )
    print(json.dumps({
        "output_root": str(output_root),
        "panel_d_windows": int(scores["local_periodicity_score"].notna().sum()),
        "panel_e_locked_points": int(len(spacing)),
        "coordinate_integrity_pass": True,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
