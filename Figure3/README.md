# Figure 3 — Atomic ordering and three-dimensional heterogeneity of the reaction intermediate

This folder reproduces the numerical analyses behind Figure 3 of the manuscript
*Data-Guided Liquid-Metal Synthesis of High-Entropy Intermetallics*. It covers
the automatic atom-column detection, local fast-Fourier-transform (FFT)
analysis, three-dimensional rendering, whole-field periodic-contrast mapping and
layer-/radius-resolved lattice-support quantification of a 4D-STEM multislice
electron ptychography reconstruction of a single reaction-intermediate
nanoparticle.

The layout mirrors the repository convention: `analysis/` holds the computation
scripts, `build/` holds the figure-rendering scripts, `source_data/` holds the
canonical manuscript-facing result tables, `data/` holds the input datasets, and
`static/` + `interactive_3d/` hold the published figures.

## Layout

```text
Figure3/
  fig3_paths.py            path resolution (relative to this folder; override with FIG3_ROOT / FIG3_OUTPUT)
  analysis/                atom-column detection, lattice locking, per-layer geometry, row analysis
  build/                   figure build/render scripts (h, i1, i2, i3, j, g1-g3)
  source_data/             canonical result tables used by the manuscript figures
  data/                    input datasets (aligned reconstruction stack, detected columns, FFT inputs)
  static/                  published panel previews
  interactive_3d/          white-background interactive 3-D package (HTML + sources + static views)
```

## What each panel computes

| Panel | Analysis | Reproducible from bundled data? |
|---|---|---|
| g1-g3 | Raw FFT intensity maps of representative Pt-rich / ordered / liquid-metal-rich regions | Input matrix bundled; canonical rendering in OriginPro |
| h | Row-resolved projected-column phase-intensity profile (206 columns, 16 complete rows) | Yes — `analysis/analyze_h_l8_projected_rows.py` + `build/render_compact_h.py` |
| i1 | 3-D rendering of detected columns across all 16 layers (incl. interactive HTML) | Yes — `build/render_white_layered_3d_readable_v2.py` |
| i2 | Layer 9 whole-field local FFT periodic-contrast map (14x14 grid of 128x128-px windows) | Yes (scores) — `build/build_layer9_de_panels.py`; canonical heatmap in OriginPro |
| i3 | Measured lattice-support fraction, 16 layers x 5 in-plane radial bins (r/R) | 80-cell table bundled; canonical heatmap in OriginPro |
| j | Centre-minus-edge lattice-support difference (r/R 0-0.4 vs 0.8-1.0) | Summary table bundled; canonical figure in OriginPro |

## How to reproduce the non-Origin panels

From the release root, with the `ptgainsnzn` conda environment active
(`scikit-image`, `pillow` and `plotly` are used by the Figure 3 scripts):

### h — row-resolved projected-column modulation (full chain)

```bash
# 1) recompute the row summary from the Layer 8 ROI (writes h_row_summary.csv)
python Figure3/analysis/analyze_h_l8_projected_rows.py --output _fig3_out

# 2) render the compact panel from that summary (writes h_compact_source_data.csv + figure)
FIG3_OUTPUT=_fig3_out python Figure3/build/render_compact_h.py
```

The regenerated `h_row_summary.csv` is byte-identical to the bundled
`source_data/h_row_summary.csv`, and the compact profile reproduces 16 rows,
206 columns, 8 adjacent row pairs and a median relative contrast of 0.147395.

### i1 — 3-D rendering and interactive HTML

```bash
FIG3_OUTPUT=_fig3_out python Figure3/build/render_white_layered_3d_readable_v2.py
```

Regenerates the static views (PNG/SVG/PDF), the interactive HTML and the
Origin-ready CSV exports under `_fig3_out/`. The underlying source tables are
byte-identical to the bundled `source_data/3d_A_*` and `3d_B_*` files; the
rendered PNG can differ cosmetically from the bundled preview because of
matplotlib version/font drift.

### i2 — Layer 9 local periodic-contrast scores

```bash
python Figure3/build/build_layer9_de_panels.py --output-root _fig3_out
```

Recomputes the 196 (14x14) local FFT peak-to-noise scores and locked-point
4-nearest-neighbour spacings. The regenerated score grid matches the bundled
`source_data/layer09_periodic_contrast_input.csv` (the Origin export matrix) to
within ~1e-13.

## Panels that require OriginPro

g1-g3, and the canonical heatmap/figure renders of i2, i3 and j, were assembled
in OriginPro. The scripts in `build/` that drive Origin
(`build_c1_c2_true_numeric_axes_v8.py`, `build_c3_true_numeric_axes_v8.py`,
`build_e_origin_radial_support.py`, `build_hj_origin2025b.py`) are included for
transparency; the underlying source matrices and result tables are bundled in
`source_data/` so the plotted numbers can be inspected without Origin.

## Inputs in `data/`

- `registered_stack.npz` — 16-layer sub-pixel-aligned ptychographic
  reconstruction stack (the input to the lattice-locking analysis).
- `atoms_db.npz`, `columns_db.npz` — automatic atom-column detection outputs
  (local maxima refined by 2-D Gaussian fits).
- `complete_lattice_coordinates_all_layers.csv` — the locked candidate lattice
  (8,748 sites across 16 layers).
- `complete_lattice_model.npz`, `per_layer_circles.csv` — lattice geometry and
  per-layer fitted circles.
- `l8_three_region_intensity_labels.csv`, `corrected_atom_columns_intensity.csv` —
  Layer 8 region labels and measured column intensities for the h analysis.
- `gray8_plain_layer_08/09_of_16_scale2nm.tif` — Layer 8 / Layer 9 grayscale
  reconstructions used by the h and i2 analyses.
- `displayed_points_three_class.csv` — locked display points for the i2
  spacing analysis.

## Scientific claim boundaries

- Layer index is a reconstruction coordinate, not a calibrated physical depth
  or reaction time.
- The h row classes (A/B) are geometric projected-row classes, not elemental
  identities; h is not a formal chemical-order parameter.
- The i2 0-100 score and the i3/j lattice-support fraction are relative
  structural metrics, not order-parameter, phase-fraction or composition maps.
- All values are spatially correlated technical measurements from a single
  reconstructed particle dataset; no inferential p value is reported.
