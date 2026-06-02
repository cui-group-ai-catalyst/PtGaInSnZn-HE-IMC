# Reviewer demo

This demo is designed for editors and reviewers who need to confirm that the
release can be installed and that a small CPU-only workflow reproduces key
reference outputs. It does not require UMA, Hugging Face access, or GPU compute.
For a less code-centric walkthrough, start with
`../REVIEWER_REPRODUCTION_GUIDE.md`.

## Setup

Run from the repository root:

```bash
conda env create -f environment.yml
conda activate ptgainsnzn
```

Expected installation time:

- 10-25 minutes for the CPU scientific stack.
- 20-45 minutes for the full environment including ML packages.

## Demo commands

```bash
python Panel_a_BinaryHeatmap/regen_FigA_data.py
python Panel_b_MultiComponentScatter/verify_FigB_central.py
python scripts/verify_release.py
```

## Expected outputs

Panel a should write:

`Panel_a_BinaryHeatmap/data_FigA_v2_FamilyOrdered_Origin_regen.csv`

Expected key value:

`DeltaH_mix(Pt-Ga) = -32.06 kJ/mol`

Panel b should write:

`Panel_b_MultiComponentScatter/data_FigB_central_anchors_regen.csv`

Expected key values:

- Pt rank: `1`
- Pt `dH_mix_kJmol`: `-24.440`
- Pt `dG_mix_500K_kJmol`: `-29.839`

`scripts/verify_release.py` should report all checks as passed for the release
files and key numerical anchors.

## Expected demo runtime

After installation, the demo should run in less than 5 minutes on a normal
desktop computer. On a typical workstation, the CPU computations themselves
should complete in less than 1 minute.

## What this demo proves

This demo verifies the fully self-contained analytical path:

- Miedema binary enthalpy calculation for Figure 1a.
- Hildebrand-Muggianu central cocktail ranking for Figure 1b.
- Presence and integrity of key reference CSVs used by UMA/post-processing
  panels.

It does not test gated UMA checkpoint access. UMA-dependent reproduction is
documented separately in `UMA_Checkpoint_Setup.md`.
