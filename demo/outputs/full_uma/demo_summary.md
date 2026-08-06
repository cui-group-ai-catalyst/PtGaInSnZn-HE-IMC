# Reviewer demo summary

- Status: **passed**
- Mode: `full-uma`
- System: `pt3_gainsnzn_bounded_validation`
- Manifest: `experimental_extensions/system_manifest.json`
- Generated: `2026-07-24T04:44:07+00:00`

## What was executed

- `release_verifier`: passed (0.212 s)
- `bounded_validation`: passed (2.247 s)
- `panel_c_uma_rerun`: passed (19.932 s)
- `si_fig3_uma_rerun`: passed (20.121 s)

## Bounded numerical evidence

### Manifold: `fixed_l12_manifold`

- Composition rows: 165
- Mixing elements: Ga, In, Sn, Zn
- Training R2: 0.9992773421034635
- Non-endmember LOOCV RMSE: 0.14225046715799528
- Group-holdout RMSE range: 0.07180603344096442 to 0.23196219085891873

### Energy backends: `uma_mp_dft_binary_rank`

- all_hosts / mp_dft__uma_s_1p1: Spearman rho=0.8678571428571429, top-k overlap=4/5, rank reversals=15
- all_hosts / mp_dft__chgnet: Spearman rho=0.95, top-k overlap=5/5, rank reversals=8
- all_hosts / uma_s_1p1__chgnet: Spearman rho=0.9035714285714285, top-k overlap=4/5, rank reversals=13
- size_compatible_hosts / mp_dft__uma_s_1p1: Spearman rho=0.9510489510489512, top-k overlap=4/5, rank reversals=5
- size_compatible_hosts / mp_dft__chgnet: Spearman rho=0.9370629370629372, top-k overlap=5/5, rank reversals=5
- size_compatible_hosts / uma_s_1p1__chgnet: Spearman rho=0.9720279720279721, top-k overlap=4/5, rank reversals=4

## Supported conclusions

- The workflow accepts a schema-checked system manifest instead of relying on a fixed row order.
- Supplied composition-energy tables can be checked for uniqueness and manifold completeness.
- Pairwise CEF interpolation can be evaluated by training fit, non-endmember LOOCV, and composition-family holdout.
- Matched supplied energy backends can be compared by error, rank correlation, top-k overlap, and ranking reversals.
- Bundled calculations and validation artifacts can be reproduced deterministically within the stated scope.

## Not established

- Accuracy or transferability for a new host, ordered prototype, or bonding class.
- Applicability to N-, O-, S-, or P-containing compounds without new models and independent validation.
- Universal accuracy of UMA, CHGNet, DFT, or any other energy backend.
- Global phase stability, kinetic accessibility, or synthesizability.

The HTML/PDF report under `bounded_validation/` is a static view of the
computed JSON/CSV evidence; it is not the computational engine.
