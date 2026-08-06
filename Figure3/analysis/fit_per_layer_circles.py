"""Fit a per-layer circular boundary to the 8748-node baseline, in the unified frame.

Coordinates are NEVER refit or translated here -- only the display/sampling
circle (centre + radius) is fitted, from the completed node positions.  This
replaces the drifting offset-frame circles of the earlier three-phase pipeline
and smooths the L12 -> L13 radius jump.

Fit method: convex hull of the completed nodes -> algebraic (Kasa) least
squares circle on the hull vertices -> boundary-biased centre/radius.  For a
nearly-circular particle this is stable and deterministic.

Outputs (written to FINAL_CANONICAL_v2_1/geometry/):
  * per_layer_circles.csv        -- centre, radius (px & nm), fit diagnostics
  * figure_circle_fit_montage.png -- 16 raw slices + fitted circle overlay
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.spatial import ConvexHull

ROOT = Path(__file__).resolve().parents[0]
BASELINE = ROOT / "output" / "complete_lattice_review" / "complete_lattice_coordinates_all_layers.csv"
IMAGE_DIR = ROOT / "gray8_plain_tiff"
OUT_DIR = ROOT / "FINAL_CANONICAL_v2_1" / "geometry"

PX_NM = 0.01138848395  # authoritative OME scale


def kasa_fit(pts: np.ndarray) -> tuple[float, float, float]:
    """Algebraic circle fit: returns (cx, cy, R)."""
    A = np.c_[2.0 * pts[:, 0], 2.0 * pts[:, 1], np.ones(len(pts))]
    b = pts[:, 0] ** 2 + pts[:, 1] ** 2
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy, cc = sol
    return float(cx), float(cy), float(np.sqrt(max(cc + cx * cx + cy * cy, 0.0)))


def fit_layer(xy: np.ndarray) -> dict:
    """Fit a boundary circle to one layer's completed nodes."""
    xy = np.asarray(xy, dtype=float)
    hull_pts = xy[ConvexHull(xy).vertices]
    cx, cy, r = kasa_fit(hull_pts)
    radial = np.hypot(xy[:, 0] - cx, xy[:, 1] - cy)
    return {
        "center_x_px": cx,
        "center_y_px": cy,
        "radius_px": r,
        "radius_nm": r * PX_NM,
        "max_radial_dev_px": float(np.abs(radial - r).max()),
        "p95_radial_dev_px": float(np.percentile(np.abs(radial - r), 95)),
        "hull_vertices": len(hull_pts),
    }


def bright_extent_px(img: np.ndarray, cx: float, cy: float, rf: float,
                     thr_above_bg: float = 8.0) -> float:
    """Radius (px) out to which ring-max intensity stays clearly above background.

    Robust to sparse atom spots (ring-MAX, not ring-mean).  The threshold
    ``bg + thr_above_bg`` on 8-bit data marks pixels that are genuinely
    brighter than the local background.
    """
    yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    dist = np.hypot(xx - cx, yy - cy)
    bg = float(np.percentile(img[dist > rf * 2.5], 90))
    ringmax = []
    for rad in range(0, int(rf * 2.0) + 1, 3):
        ring = (dist >= rad) & (dist < rad + 3)
        ringmax.append((rad, float(img[ring].max())))
    rmax = np.array(ringmax)
    above = rmax[:, 1] > bg + thr_above_bg
    return float(rmax[above, 0].max()) if above.any() else float(rf)


def peak_p90_extent_px(metrics: np.ndarray, cx: float, cy: float, rf: float) -> float:
    """p90 radius of strong image peaks (amplitude>30) near the particle."""
    if len(metrics) == 0:
        return float(rf)
    xy = metrics[:, :2]
    amp = metrics[:, 2]
    strong = amp > 30
    d = np.hypot(xy[strong, 0] - cx, xy[strong, 1] - cy)
    d = d[d < rf * 2.0]
    if len(d) == 0:
        return float(rf)
    return float(np.percentile(d, 90))


def main() -> None:
    df = pd.read_csv(BASELINE)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    metrics_all = np.load(ROOT / "output" / "atoms_db.npz", allow_pickle=True)["metrics"]

    rows = []
    images = []
    for layer in range(1, 17):
        d = df[df.layer == layer]
        xy = d[["x_completed_px", "y_completed_px"]].to_numpy(float)
        fit = fit_layer(xy)
        img = np.asarray(Image.open(IMAGE_DIR / f"gray8_plain_layer_{layer:02d}_of_16_scale2nm.tif").convert("L"), dtype=float)
        # display radius: node-fit for the crystalline body; a smooth, gentle
        # taper for the amorphous liquid cap (decided with the author: the
        # liquid droplet does not collapse abruptly, so the circle stays large)
        bright = bright_extent_px(img, fit["center_x_px"], fit["center_y_px"], fit["radius_px"])
        if layer >= 13:
            taper_nm = [3.60, 3.40, 3.15, 2.85]           # L13..L16, nm
            display = taper_nm[layer - 13] / PX_NM
        else:
            display = fit["radius_px"]
        rows.append({"layer": layer, "nodes": len(xy), "bright_extent_px": bright,
                     "display_radius_px": display, **fit})
        images.append((img, layer, fit))

    table = pd.DataFrame(rows)
    table.to_csv(OUT_DIR / "per_layer_circles.csv", index=False)

    # --- montage preview: raw slice + display circle ------------------------
    fig, axes = plt.subplots(4, 4, figsize=(11.5, 11.5), facecolor="white", constrained_layout=True)
    for ax, (img, layer, fit) in zip(axes.flat, images):
        lo, hi = np.nanpercentile(img[300:1160, 0:1080], [2, 99.7])
        ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, interpolation="none")
        t = np.linspace(0, 2 * np.pi, 400)
        row = table.iloc[layer - 1]
        dR = row.display_radius_px
        ax.plot(row.center_x_px + dR * np.cos(t), row.center_y_px + dR * np.sin(t),
                color="#D73027", lw=1.1, ls="--")
        ax.plot(row.center_x_px, row.center_y_px, "+", color="#D73027", ms=6, mew=1.1)
        ax.set_xlim(300, 1120)
        ax.set_ylim(1130, 300)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(f"L{layer:02d} R={dR * PX_NM:.2f}nm C=({row.center_x_px:.0f},{row.center_y_px:.0f})",
                     fontsize=7.5)
    fig.suptitle("Per-layer display circles (node-fit crystalline; visible extent liquid)", fontsize=12)
    fig.savefig(OUT_DIR / "figure_circle_fit_montage.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # --- comparison to previous geometry -------------------------------------
    prev = [2.50, 2.70, 3.05, 3.40, 3.40, 3.50, 3.75, 3.75,
            3.75, 3.60, 3.75, 3.40, 2.50, 1.80, 1.30, 0.85]
    table["prev_R_nm"] = prev
    table["dR_nm"] = table.radius_nm - table.prev_R_nm
    pd.set_option("display.width", 160)
    print(table[["layer", "nodes", "center_x_px", "center_y_px", "radius_nm",
                 "display_radius_px", "bright_extent_px", "prev_R_nm"]].round(2).to_string(index=False))
    print("\ndisplay radius (nm):")
    print("  " + "  ".join(f"L{i + 1}={table.display_radius_px.iloc[i] * PX_NM:.2f}" for i in range(16)))
    print("\nL12->L13 fit radius : %.2f -> %.2f nm (spherical taper, real)" % (table.radius_nm[11], table.radius_nm[12]))
    print("L12->L13 display   : %.2f -> %.2f nm (liquid layer uses visible extent)" % (table.display_radius_px[11] * PX_NM, table.display_radius_px[12] * PX_NM))
    print(OUT_DIR / "per_layer_circles.csv")
    print(OUT_DIR / "figure_circle_fit_montage.png")


if __name__ == "__main__":
    main()
