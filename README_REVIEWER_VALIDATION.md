# Reviewer-facing bounded validation attachment

This attachment supports the response to Reviewer 2, Comment 1. It packages
the configuration, source code, bundled tabular inputs, tests, and generated
evidence needed to inspect the specific methodological claims made in the
response letter.

## Quick run

From this attachment root:

```bash
python experimental_extensions/run_validation.py
```

On Windows, after activating a Python environment with the required packages,
`RUN_VALIDATION.bat` provides a double-click launcher and opens the generated
static HTML report when validation passes.

Minimal post-processing dependencies are listed in
`requirements_validation.txt`. The complete manuscript reproduction
environment, including the packages used for upstream UMA/CHGNet calculations,
is recorded in `environment.yml`. The bounded validation command consumes the
bundled CSV results and does not require a GPU or model-checkpoint download.

## What the command verifies

1. The schema, system manifest, module contracts, scientific-scope states, and
   input files are internally consistent.
2. The fixed L1_2-Pt3(Ga,In,Sn,Zn) dataset contains all 165 unique non-negative
   integer compositions on eight mixing-sublattice sites.
3. The six-parameter pairwise CEF is compared with an endmember-only ablation,
   non-endmember leave-one-composition-out validation, and Ga-count
   composition-family holdout.
4. Materials Project DFT, UMA-s-1p1, and CHGNet are compared on the bundled
   15-host M-Ga binary reference set using RMSE, MAE, Spearman correlation,
   top-five overlap, and explicit ranking reversals.
5. JSON/CSV results, a two-panel SI-ready PNG/PDF figure, and static HTML/PDF
   evidence reports are regenerated under
   `experimental_extensions/outputs/system_validation/`.

Run the extension-layer tests with:

```bash
python -m unittest discover -s experimental_extensions/tests -v
```

## Interpretation boundary

The successful run demonstrates software-level configurability for the bundled
inputs, internal interpolation of the supplied fixed-L1_2 energy manifold, and
rank consistency on the bundled binary reference set. It does not establish a
new thermodynamic theory or machine-learning algorithm, validate UMA on the
five-component landscape, prove transferability to a new host or prototype,
predict synthesizability, or establish applicability to N-, O-, S-, or
P-containing compounds.

For a new chemistry, the structures and site assignments, reference states,
energy method, competing phases, configurational sampling, interface/kinetic
assumptions, and independent validation must be supplied and re-evaluated.

## Main files

- `experimental_extensions/system_manifest.json`: bounded system definition
  and claim states.
- `experimental_extensions/system_manifest.schema.json`: machine-checkable
  manifest schema.
- `experimental_extensions/module_contracts.json`: module inputs, outputs,
  validation level, and excluded inferences.
- `experimental_extensions/run_validation.py`: unified validation entry point.
- `experimental_extensions/FORMULAS_AND_VALIDATION.md`: generic formulas and
  validation definitions.
- `experimental_extensions/outputs/system_validation/evidence_report.html`:
  browser-readable static evidence report.
- `SHA256SUMS.txt`: attachment file hashes generated at packaging time.
