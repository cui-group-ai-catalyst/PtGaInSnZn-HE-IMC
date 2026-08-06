#!/usr/bin/env python
"""
Depth-resolved atom-column analysis of a 16-slice 4D-STEM / MEP reconstruction.

Reproduces the analysis style of:
  * ACS Nano 2025, 19, 5568  (4D MEP of PMW)   -- Fig 2c-f, Fig 3
  * Science 2021, 372, 826  (multislice ptychography) -- Fig 4
for a stack of 16 reconstructed phase-slice images (Pt3(GaInSnZn) L12-HEA NP).

Pipeline
--------
  0.  Calibration & masking (px/nm from OME metadata or 2 nm scale bar).
  1.  Optional sub-pixel rigid alignment (skimage phase cross-correlation).
  2.  Atom-column detection + 2D Gaussian fitting per slice -> (x, y, z, A).
  3.  Local normalization:  order proxy  s = (A - <A_NN>) / (A + <A_NN>)
  4.  Grid mapping + 3-D order-volume export (.npz and ParaView .vtk).
  5.  Figures:
        fig1_depth_sectioning  (Science Fig 4 A/D style)
        fig2_slices_crosssection (ACS Nano Fig 2 c-f style)
        fig3_order_fwhm        (ACS Nano Fig 3 style, pseudo-Voigt FWHM)
        fig4_3d_scatter        (atom cloud coloured by order proxy)

Notes / caveats
---------------
* The input TIFFs are 8-bit images; amplitudes are therefore a *relative
  contrast / phase proxy*, not quantitative phase in radians.  A z spacing and
  a true nm calibration must be supplied for absolute numbers.
* All spatial units are px unless an axis says otherwise; ``px_nm`` converts.

Example
-------
    python pipeline_atom_depth.py --input gray8_plain_tiff \
        --output output --dz-nm 0.5 --px-nm 0.01136
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.optimize import curve_fit, least_squares
from scipy.spatial import cKDTree
from scipy.special import erfc
from skimage.feature import peak_local_max
from skimage.registration import phase_cross_correlation

LAYER_RE = re.compile(r"^(?:gray8_plain_)?layer_(\d+)_of_16(?:_scale2nm)?\.tif$")
DEFAULT_PX_NM = 0.01138848395    # OME PhysicalSizeX from the original RGB TIFFs
# bottom-left annotation (scale bar + "2 nm" text)
ANNOTATION_BOX = (35, 1185, 235, 1285)   # x0, y0, x1, y1


# --------------------------------------------------------------------------- #
# 0. I/O and calibration
# --------------------------------------------------------------------------- #
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path,
                   help="folder with layer_01..16_of_16_scale2nm.tif")
    p.add_argument("--output", type=Path, default=Path("output"))
    p.add_argument("--px-nm", type=float, default=None)
    p.add_argument("--dz-nm", type=float, default=0.5,
                   help="slice spacing in nm (0.5 used by ACS Nano MEP paper).")
    p.add_argument("--min-distance", type=float, default=8.0,
                   help="peak detection min distance in px")
    p.add_argument("--no-align", action="store_true",
                   help="skip registration (frames are already aligned)")
    p.add_argument("--z-smooth", type=float, default=1.0,
                   help="sigma (slices) for optional z smoothing of maps")
    p.add_argument("--skip-fit", action="store_true",
                   help="reuse saved atoms_db.npz / columns_db.npz "
                        "(skip the slow 2D-Gaussian fitting)")
    return p.parse_args()


def list_layers(input_dir: Path) -> list[Path]:
    hits: list[tuple[int, Path]] = []
    for path in input_dir.iterdir():
        m = LAYER_RE.match(path.name)
        if m:
            hits.append((int(m.group(1)), path))
    hits.sort()
    return [p for _, p in hits]


def read_px_nm(path: Path, fallback: float = DEFAULT_PX_NM) -> float:
    """Pixel size from OME-XML PhysicalSizeX if present, else scale-bar px/nm."""
    try:
        with Image.open(path) as im:
            desc = str(im.tag_v2.get(270, ""))
        m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', desc)
        if m:
            return float(m.group(1))
    except Exception:
        pass
    return fallback


def annotation_mask(shape: tuple[int, int]) -> np.ndarray:
    mask = np.ones(shape, dtype=bool)
    x0, y0, x1, y1 = ANNOTATION_BOX
    mask[y0:y1, x0:x1] = False
    return mask


def load_stack(paths: list[Path]) -> np.ndarray:
    stack = [np.asarray(Image.open(p).convert("L"), dtype=np.float32) for p in paths]
    return np.stack(stack)


# --------------------------------------------------------------------------- #
# 1. Alignment
# --------------------------------------------------------------------------- #
def register_stack(stack: np.ndarray, ref: int = 7) -> np.ndarray:
    """Sub-pixel rigid registration to one reference slice (default #8)."""
    out = stack.copy()
    ref_img = stack[ref, 300:1200, :1050]
    shifts = []
    for i in range(stack.shape[0]):
        mov = stack[i, 300:1200, :1050]
        shift, _, _ = phase_cross_correlation(
            ref_img, mov, upsample_factor=20, normalization="phase"
        )
        shifts.append(shift.tolist())
        if np.linalg.norm(shift) > 1e-4:
            out[i] = ndi.shift(stack[i], shift=shift, order=1, mode="nearest",
                               prefilter=False)
    print("registration shifts (row, col) vs layer 8:")
    for i, s in enumerate(shifts):
        print(f"  layer {i+1:2d}: {s[0]:+6.2f}, {s[1]:+6.2f} px")
    return out


# --------------------------------------------------------------------------- #
# 2. Atom detection + 2D Gaussian fitting
# --------------------------------------------------------------------------- #
def detect_peaks(image: np.ndarray, mask: np.ndarray, min_distance: float,
                 thresh: float) -> np.ndarray:
    return peak_local_max(
        image, min_distance=int(round(min_distance)), threshold_abs=thresh,
        exclude_border=10,
    )


def fit_gaussian(image: np.ndarray, y0: float, x0: float,
                 half=6) -> np.ndarray | None:
    """Fit a 2D Gaussian on a (2*half+1) window; returns [A, x, y, sx, sy]."""
    yy, xx = np.mgrid[y0 - half : y0 + half + 1, x0 - half : x0 + half + 1]
    lo_y, hi_y = max(0, y0 - half), min(image.shape[0], y0 + half + 1)
    lo_x, hi_x = max(0, x0 - half), min(image.shape[1], x0 + half + 1)
    if (hi_y - lo_y) < 5 or (hi_x - lo_x) < 5:
        return None
    data = image[lo_y:hi_y, lo_x:hi_x].astype(np.float64)
    gy, gx = np.mgrid[lo_y:hi_y, lo_x:hi_x]
    amp = float(image[int(y0), int(x0)])

    def resid(p):
        A, bx, by, sx, sy = p
        model = A * np.exp(-0.5 * ((gx - bx) ** 2 / sx**2 + (gy - by) ** 2 / sy**2))
        return (model - data).ravel()

    p0 = [amp, float(x0), float(y0), 3.0, 3.0]
    lb = [amp * 0.2, x0 - 3, y0 - 3, 1.2, 1.2]
    ub = [amp * 1.5, x0 + 3, y0 + 3, 7.0, 7.0]
    try:
        r = least_squares(resid, p0, bounds=(lb, ub), max_nfev=150)
        A, bx, by, sx, sy = r.x
        return np.array([A, bx, by, sx, sy], dtype=np.float64)
    except Exception:
        return None


def fit_all_atoms(stack: np.ndarray, mask: np.ndarray, min_distance: float,
                  ) -> list[dict]:
    """Fit every slice; returns one dict per layer with atom arrays."""
    layers = []
    for z, img in enumerate(stack):
        valid = img[mask]
        bg, sd = float(np.median(valid)), float(np.std(valid))
        thresh = bg + 3.5 * sd
        peaks = detect_peaks(img, mask, min_distance, thresh)
        rows = []
        for (y0, x0) in peaks:
            r = fit_gaussian(img, y0, x0)
            if r is None:
                continue
            A, x, y, sx, sy = r
            rows.append((x, y, A, sx, sy))
        arr = np.array(rows, dtype=np.float64).reshape(-1, 5)
        layers.append({"z": z, "n": len(arr), "atoms": arr,
                       "background_mean": bg, "background_std": sd})
        print(f"layer {z+1:2d}: {len(arr):4d} atoms fitted  (thresh={thresh:.1f}, "
              f"bg={bg:.1f}+/-{sd:.1f})")
    return layers


# --------------------------------------------------------------------------- #
# 3. Local normalization (order proxy)
# --------------------------------------------------------------------------- #
def order_metric(layers: list[dict], nn_radius: float) -> np.ndarray:
    """Per-atom order proxy s = (A - <A_NN>) / (A + <A_NN>).

    s>0 : atom brighter than its nearest neighbours (e.g. ordered sublattice A)
    s<0 : atom dimmer than its neighbours (sublattice B / disordered site)
    """
    out = []
    for d in layers:
        atoms = d["atoms"]
        if len(atoms) == 0:
            out.append(np.zeros((0, 6)))
            continue
        pos = atoms[:, :2]
        A = atoms[:, 2]
        tree = cKDTree(pos)
        nb = tree.query_ball_point(pos, r=nn_radius)
        s = np.empty(len(atoms))
        for i, idxs in enumerate(nb):
            nn = [j for j in idxs if j != i]
            if not nn:
                s[i] = 0.0
                continue
            mA = float(np.mean(A[nn]))
            s[i] = (A[i] - mA) / (A[i] + mA + 1e-12)
        out.append(np.column_stack([atoms, s, np.full(len(atoms), d["z"])]))
    return out


def track_columns(metrics: list[np.ndarray], z_nm: np.ndarray,
                  tol_px: float = 6.0) -> dict:
    """Track atom columns through depth.

    A common 'template' of column positions is built from the sharpest slice
    (most atoms), then for every slice the nearest detected atom within
    ``tol_px`` is assigned to each template column.  Returns

      x, y            : template positions (px)
      A(z), s(z)      : amplitude / order proxy vs depth per column
      n(z)            : whether the column was detected at depth z
    """
    n = max(len(m) for m in metrics)
    best = int(np.argmax([len(m) for m in metrics]))
    template = metrics[best]
    order = np.argsort(template[:, 2])[::-1]   # by amplitude desc
    # de-duplicate within tol so template columns are well separated
    x0s, y0s = [], []
    for i in order:
        if any(abs(template[i, 0] - x) < tol_px and abs(template[i, 1] - y) < tol_px
               for x, y in zip(x0s, y0s)):
            continue
        x0s.append(template[i, 0])
        y0s.append(template[i, 1])
    x0s, y0s = np.array(x0s), np.array(y0s)
    ncol = len(x0s)
    A = np.full((ncol, len(metrics)), np.nan)
    S = np.full((ncol, len(metrics)), np.nan)
    for z, m in enumerate(metrics):
        if len(m) == 0:
            continue
        tree = cKDTree(m[:, :2])
        dist, idx = tree.query(np.c_[x0s, y0s])
        ok = dist < tol_px
        A[ok, z] = m[idx[ok], 2]
        S[ok, z] = m[idx[ok], 5]
    return {"x": x0s, "y": y0s, "A": A, "s": S, "z_nm": z_nm,
            "template_layer": best + 1}


# --------------------------------------------------------------------------- #
# 4. Grid mapping / volume
# --------------------------------------------------------------------------- #
def map_atoms_to_grid(metrics: list[np.ndarray], shape: tuple[int, int],
                      scale: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Splat per-atom order proxies onto a z-slice stack (Gaussian weighted).

    scale>1 downsamples (e.g. 2 -> shape/2).  Returns (volume, count_volume).
    """
    h, w = shape
    nz = len(metrics)
    out_h, out_w = int(h / scale), int(w / scale)
    vol = np.zeros((nz, out_h, out_w), dtype=np.float32)
    cnt = np.zeros((nz, out_h, out_w), dtype=np.float32)
    sigma = 1.2  # smooth width in downsampled grid px
    for z, m in enumerate(metrics):
        if len(m) == 0:
            continue
        x, y, s = m[:, 0] / scale, m[:, 1] / scale, m[:, 5]
        xi = np.clip(np.round(x).astype(int), 0, out_w - 1)
        yi = np.clip(np.round(y).astype(int), 0, out_h - 1)
        acc = np.zeros((out_h, out_w), np.float32)
        cc = np.zeros((out_h, out_w), np.float32)
        np.add.at(acc, (yi, xi), s)
        np.add.at(cc, (yi, xi), 1.0)
        acc = ndi.gaussian_filter(acc, sigma)
        cc = ndi.gaussian_filter(cc, sigma)
        ok = cc > 0.05
        vol[z, ok] = acc[ok] / cc[ok]
        cnt[z, ok] = 1.0
    return vol, cnt


def write_vtk(path: Path, vol: np.ndarray, spacing: tuple[float, float, float],
              origin=(0.0, 0.0, 0.0)) -> None:
    """Write a binary VTK structured-points file for ParaView volume rendering."""
    nz, ny, nx = vol.shape
    data = np.ascontiguousarray(vol, dtype="<f4")   # little-endian float32
    with path.open("wb") as f:
        f.write(b"# vtk DataFile Version 3.0\norder_metric\nBINARY\n"
                b"DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx} {ny} {nz}\n".encode())
        f.write(f"ORIGIN {origin[0]} {origin[1]} {origin[2]}\n".encode())
        f.write(f"SPACING {spacing[0]} {spacing[1]} {spacing[2]}\n".encode())
        f.write(f"POINT_DATA {nx * ny * nz}\n".encode())
        f.write(b"SCALARS order_metric float 1\nLOOKUP_TABLE default\n")
        f.write(data.tobytes())
    print("wrote", path)


# --------------------------------------------------------------------------- #
# 5. Pseudo-Voigt / Gaussian fits for depth profiles
# --------------------------------------------------------------------------- #
def pseudo_voigt(z, A, z0, sigma, gamma, c):
    """Pseudo-Voigt (equal-weighted Gaussian+Lorentzian) + baseline c."""
    return (c + A * (0.5 * np.exp(-0.5 * ((z - z0) / sigma) ** 2)
                     + 0.5 / (1 + ((z - z0) / gamma) ** 2)))


def errfunc_profile(z, A, z0, sigma, c):
    """Error-function (S-shaped) depth profile: A*erfc(-(z-z0)/sigma)+c."""
    return c + A * erfc((z0 - z) / (sigma * np.sqrt(2.0)))


def fwhm_pseudo_voigt(sigma, gamma):
    """FWHM of the pseudo-Voigt built from gauss sigma and lorentz gamma."""
    fg, fl = 2.0 * sigma * np.sqrt(2 * np.log(2.0)), 2.0 * gamma
    # standard pseudo-Voigt FWHM approximation (equal weights)
    return 0.5 * fg + 0.5 * fl


def fit_depth_profiles(metrics: list[np.ndarray], z_nm: np.ndarray,
                       rois: list[tuple[int, int, int, int]],
                       radius_px: float = 3.0) -> list[dict]:
    """Fit pseudo-Voigt to order-metric vs depth for each ROI (mean over ROI)."""
    results = []
    for (y0, x0, h, w) in rois:
        prof = []
        for z, m in enumerate(metrics):
            if len(m) == 0:
                prof.append(np.nan)
                continue
            sel = (m[:, 0] > x0 - w / 2) & (m[:, 0] < x0 + w / 2) & \
                  (m[:, 1] > y0 - h / 2) & (m[:, 1] < y0 + h / 2)
            if sel.sum() == 0:
                prof.append(np.nan)
            else:
                prof.append(float(np.nanmean(m[sel, 5])))
        prof = np.array(prof)
        ok = ~np.isnan(prof)
        try:
            p0 = [np.nanmax(prof) - np.nanmin(prof), z_nm[ok][np.argmax(prof[ok])],
                  1.5, 1.5, np.nanmin(prof)]
            popt, _ = curve_fit(pseudo_voigt, z_nm[ok], prof[ok], p0=p0,
                                maxfev=20000)
            A, z0, sigma, gamma, c = popt
            results.append({"roi": (y0, x0, h, w), "z0_nm": z0,
                            "fwhm_nm": fwhm_pseudo_voigt(sigma, gamma),
                            "sigma": sigma, "gamma": gamma, "A": A, "c": c,
                            "profile": prof, "z_nm": z_nm,
                            "popt": popt})
        except Exception as e:
            print("  fit failed for ROI", (y0, x0, h, w), e)
            results.append({"roi": (y0, x0, h, w), "profile": prof,
                            "z_nm": z_nm, "fwhm_nm": np.nan, "z0_nm": np.nan})
    return results


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def style_axes(ax, fontsize=7):
    ax.tick_params(labelsize=fontsize, width=0.6, length=3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def add_scalebar(ax, px_nm, y_size, nm=2.0, x0=20):
    length = nm / px_nm
    ax.plot([x0, x0 + length], [y_size - 30] * 2, color="w", lw=2,
            solid_capstyle="butt")
    ax.text(x0 + length / 2, y_size - 42, f"{nm:.0f} nm", color="w",
            ha="center", va="top", fontsize=6,
            bbox=dict(facecolor="black", edgecolor="none", pad=1))


def _draw_circle(ax, pmask: np.ndarray, scale: float) -> None:
    """Outline the round-particle mask on a downsampled image axis."""
    ys, xs = np.where(pmask)
    if len(ys) == 0:
        return
    # fit circle to mask boundary (already circular in this pipeline)
    cx, cy = xs.mean() / scale, ys.mean() / scale
    r = np.sqrt(((xs - cx * scale) ** 2 + (ys - cy * scale) ** 2).mean()) / scale
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(cx + r * np.cos(th), cy + r * np.sin(th), color="red", lw=0.9,
            ls="--")


def fig_depth_sectioning(stack, normalized, z_nm, px_nm, columns, out: Path,
                         n_cols=6):
    """Science Fig 4 A/D style: slices at depths + column depth profiles."""
    h, w = stack.shape[1:]
    selected = [3, 6, 8, 10, 13]
    vmin, vmax = np.nanpercentile(normalized, [1, 99.5])
    fig = plt.figure(figsize=(7.6, 6.4), facecolor="white")
    gs = fig.add_gridspec(2, 6, height_ratios=[1.6, 1.4], hspace=0.5, wspace=0.35)
    for i, idx in enumerate(selected):
        ax = fig.add_subplot(gs[0, i])
        ax.imshow(normalized[idx - 1], cmap="gray", vmin=vmin, vmax=vmax,
                  interpolation="none")
        ax.set_title(f"z = {z_nm[idx-1]:.1f} nm", fontsize=7)
        ax.set_axis_off()
        if i == 0:
            add_scalebar(ax, px_nm, h)
    # depth profiles of the same tracked columns through the stack
    ax = fig.add_subplot(gs[1, :])
    x0s, y0s = columns["x"], columns["y"]
    A = columns["A"]
    # pick well-separated bright columns in a central band
    in_band = (y0s > 500) & (y0s < 700) & (x0s > 400) & (x0s < 1000)
    cand = np.where(in_band)[0]
    med = np.nanmedian(A[cand], axis=1)
    cand = cand[np.argsort(med)[::-1]]
    keep = []
    for c in cand:
        if len(keep) >= n_cols:
            break
        if all(abs(x0s[c] - x0s[k]) > 80 for k in keep):
            keep.append(c)
    for c in keep:
        prof = A[c]
        ax.plot(z_nm, prof, "o-", ms=3, lw=0.8, alpha=0.85)
    ax.set_xlabel("z (nm)")
    ax.set_ylabel("atom amplitude (8-bit)")
    ax.set_title("depth profiles of tracked atom columns (Science Fig 4 style)")
    style_axes(ax)
    fig.suptitle("Depth sectioning of the 16-slice stack", fontsize=10, y=0.99)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved", out)


def fig_slices_crosssection(metrics, z_nm, px_nm, shape, order_vol, scale, out,
                            pmask: np.ndarray | None = None):
    """ACS Nano Fig 2 c-f style: maps, depth cross-section, 3-D cloud."""
    nz, ny, nx = order_vol.shape
    selected = [3, 6, 8, 10, 13]
    vm = np.nanpercentile(order_vol, [2, 98])
    fig = plt.figure(figsize=(7.6, 7.8), facecolor="white")
    gs = fig.add_gridspec(3, 6, height_ratios=[1.25, 1.25, 1.4], hspace=0.42,
                          wspace=0.4)
    # top: order-metric maps at selected depths (with round-particle outline)
    for i, idx in enumerate(selected):
        ax = fig.add_subplot(gs[0, i])
        ax.imshow(order_vol[idx - 1], cmap="coolwarm", vmin=vm[0], vmax=vm[1],
                  interpolation="none")
        if pmask is not None:
            _draw_circle(ax, pmask, scale)
        ax.set_title(f"z = {z_nm[idx-1]:.1f} nm", fontsize=7)
        ax.set_axis_off()
        if i == 0:
            add_scalebar(ax, px_nm * scale, ny)  # downsampled px size
    # middle: x-z cross-section along a horizontal line y0
    y_line = 700
    ax = fig.add_subplot(gs[1, :])
    yl = int(y_line / scale)
    xz = order_vol[:, yl, :]
    ax.imshow(xz, cmap="coolwarm", vmin=vm[0], vmax=vm[1], aspect="auto",
              origin="lower",
              extent=[0, nx * scale * px_nm, z_nm[0], z_nm[-1]],
              interpolation="none")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("z (nm)")
    ax.set_title(f"depth section along y = {y_line*px_nm:.2f} nm "
                 "(ACS Nano Fig 2f style)")
    style_axes(ax)
    # bottom: 3-D scatter of atom order proxies (subsampled, inside round particle)
    ax = fig.add_subplot(gs[2, :], projection="3d")
    xs, ys, zs, ss = [], [], [], []
    for z, m in enumerate(metrics):
        if len(m) == 0:
            continue
        if pmask is not None:
            xi = np.clip(m[:, 0].astype(int), 0, pmask.shape[1] - 1)
            yi = np.clip(m[:, 1].astype(int), 0, pmask.shape[0] - 1)
            keep = pmask[yi, xi]
            m = m[keep]
        if len(m) == 0:
            continue
        xs.append(m[:, 0]); ys.append(m[:, 1])
        zs.append(np.full(len(m), z_nm[z]))
        ss.append(m[:, 5])
    X = np.concatenate(xs); Y = np.concatenate(ys)
    Z = np.concatenate(zs); S = np.concatenate(ss)
    rng = np.random.default_rng(0)
    pick = rng.choice(len(X), size=min(8000, len(X)), replace=False)
    im = ax.scatter(X[pick] * px_nm, Y[pick] * px_nm, Z[pick], c=S[pick],
                    s=2, cmap="coolwarm", vmin=vm[0], vmax=vm[1], depthshade=False)
    ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)"); ax.set_zlabel("z (nm)")
    ax.set_title("atom cloud, coloured by local order proxy")
    fig.colorbar(im, ax=ax, shrink=0.5, pad=0.1, label="order proxy")
    fig.suptitle("Normalised contrast / order proxy vs depth (ACS Nano Fig 2 style)",
                 fontsize=10, y=0.99)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved", out)


def fit_column_envelopes(columns: dict, z_nm: np.ndarray,
                         max_cols: int = 400) -> list[dict]:
    """Fit pseudo-Voigt to the *amplitude vs depth* envelope of tracked columns.

    The 8-bit frames carry no measurable per-slice L12 checkerboard order
    contrast, so this measures the *depth envelope* of the reconstruction
    (focus / effective sectioning), which is the well-defined signal here.
    """
    A = columns["A"]
    rng = np.random.default_rng(1)
    ncol = A.shape[0]
    pick = rng.choice(ncol, size=min(max_cols, ncol), replace=False)
    fits = []
    for c in pick:
        prof = A[c]
        ok = ~np.isnan(prof)
        if ok.sum() < 6:
            continue
        z = z_nm[ok]; y = prof[ok]
        try:
            p0 = [y.max() - y.min(), z[np.argmax(y)], 1.5, 1.5, y.min()]
            popt, _ = curve_fit(pseudo_voigt, z, y, p0=p0, maxfev=20000)
            A_, z0, sigma, gamma, c_ = popt
            fwhm = fwhm_pseudo_voigt(sigma, gamma)
            if not (0.3 <= fwhm <= 20.0):
                continue
            resid = y - pseudo_voigt(z, *popt)
            ss = np.sum((y - y.mean()) ** 2)
            r2 = 1 - np.sum(resid**2) / max(ss, 1e-12)
            fits.append({"col": int(c), "z0_nm": z0, "fwhm_nm": fwhm,
                         "r2": r2, "popt": popt, "z": z, "y": y})
        except Exception:
            continue
    return fits


def fig_order_fwhm(fits, z_nm, out):
    """ACS Nano Fig 3 style: depth envelopes + pseudo-Voigt FWHM."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), facecolor="white")
    ax = axes[0]
    for f in fits[:: max(1, len(fits) // 12)][:12]:
        ax.plot(f["z"], f["y"], "o", ms=2.5, alpha=0.5)
        zf = np.linspace(z_nm[0], z_nm[-1], 100)
        ax.plot(zf, pseudo_voigt(zf, *f["popt"]), "-", lw=1)
    ax.set_xlabel("z (nm)")
    ax.set_ylabel("atom amplitude (8-bit)")
    ax.set_title("column depth envelopes + pseudo-Voigt fits")
    style_axes(ax)
    fwhm = np.array([f["fwhm_nm"] for f in fits])
    ax2 = axes[1]
    if len(fwhm):
        ax2.hist(fwhm, bins=16, color="#4daf4a", edgecolor="white")
        ax2.axvline(fwhm.mean(), color="k", ls="--", lw=1)
        ax2.text(0.97, 0.9, f"mean FWHM =\n{fwhm.mean():.2f} nm\n"
                 f"(n={len(fwhm)})", transform=ax2.transAxes, ha="right",
                 va="top", fontsize=7)
    ax2.set_xlabel("depth-envelope FWHM (nm)")
    ax2.set_ylabel("count")
    ax2.set_title("reconstruction depth extent (FWHM)")
    style_axes(ax2)
    fig.suptitle("Depth envelope of atom columns (pseudo-Voigt FWHM)",
                 fontsize=10, y=0.98)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved", out, f"({len(fwhm)} stable fits)")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    args = parse_args()
    input_dir = args.input.expanduser().resolve()
    out_dir = args.output.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    paths = list_layers(input_dir)
    if len(paths) != 16:
        print("warning: expected 16 layers, found", len(paths))
    stack = load_stack(paths)
    mask = annotation_mask(stack.shape[1:])
    px_nm = args.px_nm or read_px_nm(paths[0])
    dz_nm = args.dz_nm
    z_nm = np.arange(1, stack.shape[0] + 1) * dz_nm
    print(f"loaded {len(paths)} slices, shape {stack.shape[1:]}, "
          f"px_nm={px_nm:.5f}, dz_nm={dz_nm}")

    # 1. align
    if args.no_align:
        registered = stack
    else:
        registered = register_stack(stack)

    # 2. atom fitting
    nn_radius = 1.6 * args.min_distance * 2.0   # ~2x nearest-neighbour
    if args.skip_fit:
        db = np.load(out_dir / "atoms_db.npz", allow_pickle=True)
        metrics = [np.asarray(a) for a in db["metrics"]]
        n_atoms = sum(len(m) for m in metrics)
        print(f"loaded {n_atoms} atoms from atoms_db.npz (skip-fit)")
    else:
        layers = fit_all_atoms(registered, mask, args.min_distance)
        n_atoms = sum(d["n"] for d in layers)
        print(f"total atoms fitted: {n_atoms}")
        if n_atoms == 0:
            raise SystemExit("no atoms found -- check --min-distance / thresholds")
        # 3. order proxy (nearest neighbours in the same slice)
        metrics = order_metric(layers, nn_radius)

    # 4. grid + volume
    scale = 2
    order_vol, cnt = map_atoms_to_grid(metrics, stack.shape[1:], scale=scale)
    if args.z_smooth and args.z_smooth > 0:
        order_vol = ndi.gaussian_filter(order_vol, sigma=(args.z_smooth, 1.0, 1.0))

    # exports
    np.savez_compressed(
        out_dir / "atoms_db.npz",
        metrics=np.asarray(metrics, dtype=object),
        px_nm=px_nm, dz_nm=dz_nm, z_nm=z_nm,
        nn_radius=nn_radius,
    )
    np.savez_compressed(out_dir / "order_volume.npz",
                        order_vol=order_vol, scale=scale, px_nm=px_nm,
                        dz_nm=dz_nm, z_nm=z_nm)
    write_vtk(out_dir / "order_volume.vtk", order_vol,
              spacing=(px_nm * scale, px_nm * scale, dz_nm))

    # 5. figures
    columns = track_columns(metrics, z_nm, tol_px=0.5 * nn_radius)
    print(f"tracked {len(columns['x'])} columns across {stack.shape[0]} slices "
          f"(template from layer {columns['template_layer']})")
    np.savez_compressed(out_dir / "columns_db.npz",
                        x=columns["x"], y=columns["y"], A=columns["A"],
                        s=columns["s"], z_nm=z_nm)

    vals = registered[:, mask]
    norm = (registered - vals.mean()) / vals.std()
    fig_depth_sectioning(registered, norm, z_nm, px_nm, columns,
                         out_dir / "fig1_depth_sectioning.png")
    # round-particle mask (if available) keeps figures honest to the round NP
    pmask = None
    mask_path = out_dir / "particle_mask.npz"
    if mask_path.exists():
        pmask = np.load(mask_path)["mask"]
        print(f"using round-particle mask (fraction {pmask.mean():.2f})")
    fig_slices_crosssection(metrics, z_nm, px_nm, stack.shape[1:],
                            order_vol, scale,
                            out_dir / "fig2_slices_crosssection.png",
                            pmask=pmask)

    # Depth-envelope FWHM (ACS Nano Fig 3 style) from tracked columns
    env_fits = fit_column_envelopes(columns, z_nm)
    fig_order_fwhm(env_fits, z_nm, out_dir / "fig3_order_fwhm.png")
    fwhm_vals = [f["fwhm_nm"] for f in env_fits]
    print(f"depth-envelope FWHM: median={np.median(fwhm_vals):.2f} nm, "
          f"mean={np.mean(fwhm_vals):.2f} nm (n={len(fwhm_vals)})")

    meta = {
        "input": [p.name for p in paths],
        "px_nm": px_nm, "dz_nm": dz_nm,
        "n_atoms": n_atoms,
        "n_atoms_per_layer": ([d["n"] for d in layers]
                              if not args.skip_fit
                              else [len(m) for m in metrics]),
        "min_distance_px": args.min_distance,
        "nn_radius_px": nn_radius,
        "order_metric": "(A - mean(A_NN)) / (A + mean(A_NN)) "
                        "(relative contrast proxy)",
        "depth_fwhm_nm_median": float(np.median(fwhm_vals)),
        "depth_fwhm_nm_mean": float(np.mean(fwhm_vals)),
        "depth_fwhm_n_count": len(fwhm_vals),
        "caveat": "8-bit image-derived contrast; relative order proxy, "
                  "not quantitative chemical phase",
    }
    (out_dir / "pipeline_metadata.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(meta, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
