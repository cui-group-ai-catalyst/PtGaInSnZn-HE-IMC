# Release notes - v1.4.0 (proposed)

Status: working-tree candidate. Combines the v1.3.2 scientific-audit
corrections already present in the working tree with a new Figure 3
reproducibility package. Exact-version archive (DOI) pending; version number is
a proposal and can be renumbered before archiving.

## Added

### Figure 3 - atomic ordering and 3-D heterogeneity of the reaction intermediate

A complete reproducibility package for the Figure 3 structural analyses,
derived from a 4D-STEM multislice electron ptychography reconstruction of a
single reaction-intermediate nanoparticle:

- `Figure3/analysis/` - automated atom-column detection (local maxima refined
  by 2-D Gaussian fits), lattice locking (8,748 candidate sites across 16
  layers), per-layer circle geometry, and the Layer 8 row-resolved column
  analysis.
- `Figure3/build/` - figure build/render scripts for panels g1-g3, h, i1, i2,
  i3 and j.
- `Figure3/source_data/` - canonical manuscript-facing result tables
  (h row profile and row pairs, j centre-edge summary, i2 periodic-contrast
  matrix, i3 80-cell lattice-support audit, g1-g3 FFT input matrices).
- `Figure3/data/` - input datasets, including the 16-layer aligned
  reconstruction stack, detected-column database, locked lattice and Layer 8/9
  reconstructions.
- `Figure3/interactive_3d/` - the white-background interactive 3-D package
  (HTML, origin-ready exports, records, static views).
- `Figure3/Figure3_Code_Map.csv` and `Figure3/README.md` following the
  repository's Figure 1 conventions.

Panels h, i1 and the i2 score computation are fully reproducible from the
bundled data; verification showed the regenerated h row summary and i2 score
grid match the bundled canonical tables to within floating-point tolerance
(h: 16 rows / 206 columns / 8 pairs / median relative contrast 0.147395;
i2: 196-window 14x14 grid, max absolute difference ~1e-13; i3/j Layer 9
centre-edge fraction 0.969/0.696, difference +0.274). Panels g1-g3 and the
canonical i2/i3/j heatmap renders were assembled in OriginPro; the input
matrices and Origin-driving scripts are bundled for transparency.

### Bundled v1.3.2 scientific-audit corrections

The working tree already contained the v1.3.2 corrections; they remain part of
this release:

- Removed repeated cube-root errors in the SI Fig. 2 and Panel e Miedema
  electron-density mismatch terms.
- Gibbs-Duhem-consistent multicomponent regular-solution liquid mixing.
- Panel b `c_host` screening/sensitivity clarification.
- Panel f wetting uncertainty zones.
- Ordered-ensemble Panel c/d extensions and liquid-thermodynamics validator.
- Reviewer demo and P1 validation entry point.

## Scope

Figure 3 values are spatially correlated technical measurements from a single
reconstructed particle dataset. The i2 periodic-contrast score and i3/j
lattice-support fraction are relative structural metrics, not order-parameter,
phase-fraction or composition maps; the h row classes are geometric projected-row
classes, not elemental identities. Layer index is a reconstruction coordinate,
not a calibrated physical depth or reaction time.
