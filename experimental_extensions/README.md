# Experimental extension layer

This directory is intentionally separate from the manuscript reproduction
path. The scripts under the panel folders remain the source of the canonical
manuscript outputs. Nothing in this directory is evidence that the workflow is
scientifically transferable to a new host, prototype, or bonding class.

The extension layer demonstrates one narrower capability: an integer
composition manifold and its pairwise CEF representation can be described by
a configuration file instead of hard-coded Pt/Ga/In/Sn/Zn constants.

## Regression example

From the release root, run:

```bash
python experimental_extensions/run_manifold.py \
  experimental_extensions/pt3_gainsnzn_regression.json
```

The example reads the canonical 165-composition L1_2-Pt3(Ga,In,Sn,Zn)
landscape and writes only to:

```text
experimental_extensions/outputs/pt3_gainsnzn_regression/
```

Expected checks include:

- 165 unique integer compositions on 8 mixing sites;
- 6 independent pair-interaction parameters;
- training R2 approximately 0.999277;
- non-endmember leave-one-out RMSE approximately 0.14225 kJ mol-1 atom-1.

Run the regression test with:

```bash
python -m unittest experimental_extensions.tests.test_regression
```

## Unified bounded-validation entry point

The P1 interface is defined by `system_manifest.json` and
`system_manifest.schema.json`. From the release root, run:

```bash
python experimental_extensions/run_validation.py
```

This single command performs:

- manifest and module-contract checks;
- complete-manifold, CEF, LOOCV, and composition-family holdout validation;
- pairwise comparison of Materials Project DFT, UMA-s-1p1, and CHGNet using
  RMSE, Spearman correlation, top-k overlap, and explicit ranking reversals;
- generation of a two-panel SI-ready validation figure;
- generation of machine-readable JSON/CSV plus static HTML/PDF reports.

Outputs are isolated under:

```text
experimental_extensions/outputs/system_validation/
```

The generic equations and exact interpretation boundaries are documented in
`FORMULAS_AND_VALIDATION.md`. `module_contracts.json` states which inputs and
outputs each executable module accepts and which scientific inferences are
excluded.

## What can be configured

- host and prototype labels;
- mixing-sublattice elements and site count;
- CSV count-column mapping;
- energy column and units;
- optional complete-integer-manifold validation;
- optional composition-family holdout element.
- a system-level module list and bounded claim states;
- any matched set of three or more tabulated energy backends;
- subset filters and top-k size for backend comparison.

## What this does not validate

For a new material system, users must still provide and independently validate:

1. physically appropriate structures and site assignments;
2. consistent elemental or compound reference states;
3. an interatomic potential or electronic-structure method suitable for the
   new chemistry;
4. configurational sampling and competing phases;
5. interface and kinetic models appropriate to the bonding class;
6. comparison against independent calculations, databases, or experiments.

A successful run therefore proves software-level input substitution and CEF
interpolation only. It does not prove synthesizability or scientific transfer.
