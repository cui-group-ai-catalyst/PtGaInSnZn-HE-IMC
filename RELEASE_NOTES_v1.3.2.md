# Release notes - v1.3.2 scientific-audit revision

Status: development working tree; exact-version archive pending.

## Corrected

- Removed a repeated cube root from the SI Fig. 2 Miedema electron-density
  mismatch term.
- Removed the same repeated cube-root error from the Panel e liquid Miedema
  interaction term and renamed the local field to `n_ws_one_third`.
- Replaced the pairwise-binary liquid partial-molar approximation with the
  Gibbs-Duhem-consistent multicomponent regular-solution expression. The
  regenerated reaction energy is -150.5481 kJ mol^-1 f.u.^-1; the change from
  the prior output is -0.0127 kJ mol^-1 f.u.^-1 and does not change its sign.
- Clarified that the Panel b `c_host=0.25` value is a historical screening
  anchor rather than A3B stoichiometry, and added a `c_host=0.75` sensitivity
  comparison.
- Replaced binary wetting interpretation with favourable, indeterminate, and
  unfavourable zones under the stated +/-0.20 J m^-2 model uncertainty.
- Corrected the documented sign convention for total reaction free energy.
- Distinguished per-alloy-atom CEF interaction parameters from per-beta-site
  parameters.
- Clarified that a fixed Pt3X composition manifold determines total free
  energy and beta-sublattice diffusion potentials, not five unique absolute
  elemental chemical potentials.
- Corrected the SI Fig. 4 description to fixed-lattice single-point energies
  with one representative occupancy per integer composition.

## Added

- The original 31-row Panel c UMA energy table, deterministic ordered/random
  structure generator, optional CIF export, full UMA rerun, and JSON comparison
  report for seeds 100-129.
- A 30-member symmetry-distinct L1_2 B-sublattice ordered ensemble, retaining
  the historical ordered anchor, plus measured UMA energies, CIFs, a 30+30
  distribution comparison, confidence interval, and machine-readable checks.
- An exhaustive 68-class ordered sensitivity calculation showing that the
  30-class sample reproduces the degeneracy-weighted all-class mean to within
  0.001 kJ mol^-1 atom^-1.
- A revised Panel d support path driven by the measured 30 ordered and original
  30 disordered energies, with ensemble SD bands, raw 0 K points, ideal-entropy
  bounds, key-temperature CSVs, and machine-readable validation.
- A liquid-thermodynamics validator covering the Miedema table convention,
  Omega symmetry, the Euler/Gibbs-Duhem identity, CSV regeneration, and the
  downstream reaction-energy identity.
- The archived SI Fig. 3 structure tree (15 M-Ga binaries plus 16 elemental
  references), SHA-256 source manifest, executable `--rerun-uma` path, and
  machine-readable comparison against the canonical UMA columns.
- CEF interpolation validation with leave-one-composition-out and Ga-count
  group holdout checks.
- Gauge-invariant beta-sublattice diffusion-potential outputs.
- A configuration-driven experimental extension layer, kept separate from the
  frozen manuscript reproduction path.
- Regression tests and release-verifier checks for the new evidence artifacts.
- A schema-checked `system_manifest.json` and machine-readable module contracts.
- A unified P1 validation entry point covering CEF, LOOCV, composition-family
  holdout, and three-backend reference comparison.
- Pairwise RMSE, Spearman, top-k overlap, and explicit ranking-reversal outputs
  for Materials Project DFT, UMA-s-1p1, and CHGNet.
- A deterministic two-panel SI validation figure and static HTML/PDF reviewer
  report with enforced non-transferability and non-synthesizability boundaries.
- Generic tests using a shuffled three-element manifold to exclude hard-coded
  four-element names and row-order dependence.
- A one-command reviewer demo with quick, full-UMA, and custom-manifest paths,
  isolated outputs, step logs, manifest hashing, machine-readable metric
  extraction, and explicit supported/excluded scientific conclusions.
- Reviewer-demo unit tests showing that metric extraction does not depend on
  the manuscript module IDs.
- A runnable synthetic three-element example that replaces element labels,
  count columns, module IDs, backend IDs, and energy columns through the same
  manifest-driven entry point without making a scientific transfer claim.

## Scope

These changes improve reproducibility, units, uncertainty reporting, and
claim boundaries. They do not establish transferability to new hosts,
prototypes, nonmetal bonding classes, or synthesis routes.

The 30+30 Panel c and ensemble-driven Panel d results are currently supporting
evidence. Historical canonical Panel c/d files remain in place for traceability;
the manuscript, captions, and canonical filenames must be promoted together.
