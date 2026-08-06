# Reviewer demo: adaptability with bounded validation

## Purpose

This is a one-command evidence package for Reviewer 2 Comment 1. It shows
that the released analysis is controlled by explicit data and module
contracts, can accept a schema-compliant system manifest, and automatically
checks the supplied numerical data. The web and PDF files are reports made
from computed JSON/CSV results; they are not the computational engine.

The precise claim is **software and data-interface adaptability**. The demo
does not by itself establish scientific transferability to a new host,
prototype, energy landscape, or bonding class.

## Setup

From the repository root:

```bash
conda env create -f environment.yml
conda activate ptgainsnzn
```

For the CPU-only quick mode, the smaller validation dependency set is enough:

```bash
python -m pip install -r requirements_validation.txt
```

## One-command modes

### Quick: released tables, no UMA checkpoint

```bash
python demo/run_reviewer_demo.py --mode quick
```

This runs the release integrity checks and the bounded-validation pipeline.
It validates the 165-row fixed L1_2-Pt3(Ga,In,Sn,Zn) composition manifold,
fits the pairwise CEF representation, performs non-endmember leave-one-out
and composition-family holdout tests, and compares MP-DFT, UMA-s-1p1, and
CHGNet on the same 15-entry M-Ga reference set.

Default output:

```text
demo/outputs/quick/
  demo_summary.json
  demo_summary.md
  logs/
  bounded_validation/
    validation_results.json
    validation_evidence.png
    validation_evidence.pdf
    evidence_report.html
    evidence_report.pdf
    fixed_l12_manifold/*.csv
    uma_mp_dft_binary_rank/*.csv
```

### Full UMA: recompute the model-dependent energies

```bash
python demo/run_reviewer_demo.py --mode full-uma
```

This adds two gated-checkpoint calculations:

- Panel c: one ordered structure and 30 deterministic disordered structures
  generated with seeds 100-129;
- SI Fig. 3: 15 binary structures and 16 elemental-reference structures.

The command requires the external UMA-s-1p1 checkpoint. Set
`UMA_CHECKPOINT`, or pass `--checkpoint PATH`; use `--device cuda` only when
the installed scientific environment supports it. Checkpoint weights are not
redistributed by this repository.

### Custom manifest: replace supplied inputs through the same contracts

```bash
python demo/run_reviewer_demo.py \
  --mode quick \
  --manifest demo/examples/three_element/system_manifest.json
```

The bundled `demo/examples/three_element/` directory is an executable
input-substitution example and a template for user-supplied tables.

The custom manifest may point to:

- a composition-energy CSV plus a manifold-regression config; and/or
- a matched CSV containing two or more precomputed energy-backend columns.

The input files must follow `experimental_extensions/module_contracts.json`
and the manifest must retain the bounded claim states enforced by
`experimental_extensions/system_manifest.schema.json`. The code does not
silently invent structures, reference states, or energies for a new system.

This custom-input route demonstrates that element names, count-column names,
site count, energy column, unit, module IDs, backend columns, subset filters,
and top-k settings are replaceable. A generic three-element shuffled-row unit
test additionally checks that the manifold code is not tied to four elements
or a fixed CSV row order.

## How it works

1. The system manifest declares the modules, scientific scope, input paths,
   units, and non-negotiable claim boundaries.
2. Contract validation rejects missing fields, unsupported modules, unsafe
   transferability claims, incomplete manifolds, duplicate compositions, and
   incompatible backend tables.
3. The manifold module fits the pairwise CEF terms from the supplied energies
   and reports fit, LOOCV, and composition-family holdout errors.
4. The backend module compares matched supplied columns using RMSE, MAE,
   Spearman correlation, top-k overlap, and explicit rank reversals.
5. The runner records commands, runtimes, input-manifest hash, numerical
   metrics, artifacts, and claim boundaries in `demo_summary.json`.

## What a passing run supports

- The released Pt3(Ga,In,Sn,Zn) calculations are reproducible within the
  bundled numerical scope.
- The software accepts schema-compliant substitutions of labels, composition
  mappings, supplied energy tables, and comparison backends.
- Within the fixed supplied manifold, interpolation quality and sensitivity
  to held-out compositions are quantified rather than asserted verbally.
- Agreement and disagreement between supplied backends are visible, including
  ranking reversals; the code does not reduce the comparison to one favourable
  correlation coefficient.

## What a passing run does not support

- scientific transferability to a new host or ordered prototype;
- direct use for N-, O-, S-, or P-containing compounds;
- universal accuracy of UMA, CHGNet, DFT, or another potential;
- competing-phase stability, kinetic accessibility, or synthesizability.

Those claims require new structures, reference states, bonding-appropriate
thermodynamic and interface models, validation of the selected energy method,
competing-phase analysis, and independent calculations or experiments.

## Tests

```bash
python -m unittest discover -s experimental_extensions/tests
python -m unittest discover -s demo/tests
python scripts/verify_release.py
```

The quick mode should finish in several minutes after installation. Full UMA
runtime depends on checkpoint loading, hardware, and the selected device.
