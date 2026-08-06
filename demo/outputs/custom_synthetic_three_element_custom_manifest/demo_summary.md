# Reviewer demo summary

- Status: **passed**
- Mode: `quick`
- System: `synthetic_three_element_custom_manifest`
- Manifest: `demo/examples/three_element/system_manifest.json`
- Generated: `2026-07-24T04:33:48+00:00`

## What was executed

- `release_verifier`: passed (0.318 s)
- `bounded_validation`: passed (2.562 s)

## Bounded numerical evidence

### Manifold: `renamed_composition_module`

- Composition rows: 15
- Mixing elements: A, B, C
- Training R2: 1.0
- Non-endmember LOOCV RMSE: 1.4332917616497527e-16
- Group-holdout RMSE range: 0.0 to 3.8459253727671276e-16

### Energy backends: `renamed_backend_module`

- all_materials / reference_backend__candidate_backend: Spearman rho=0.9999999999999999, top-k overlap=2/2, rank reversals=0

## Supported conclusions

- The workflow accepts a schema-checked system manifest instead of relying on a fixed row order.
- Supplied composition-energy tables can be checked for uniqueness and manifold completeness.
- Pairwise CEF interpolation can be evaluated by training fit, non-endmember LOOCV, and composition-family holdout.
- Matched supplied energy backends can be compared by error, rank correlation, top-k overlap, and ranking reversals.
- The custom supplied tables completed the same contract, validation, visualization, and reporting path.

## Not established

- Accuracy or transferability for a new host, ordered prototype, or bonding class.
- Applicability to N-, O-, S-, or P-containing compounds without new models and independent validation.
- Universal accuracy of UMA, CHGNet, DFT, or any other energy backend.
- Global phase stability, kinetic accessibility, or synthesizability.

The HTML/PDF report under `bounded_validation/` is a static view of the
computed JSON/CSV evidence; it is not the computational engine.
