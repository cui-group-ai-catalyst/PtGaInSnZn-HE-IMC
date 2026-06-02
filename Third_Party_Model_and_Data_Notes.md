# Third-Party Models, Data and Parameters

This document declares every external model, dataset, and parameter table
that the scripts in this package consume but do not redistribute. Each entry
states the source, the version pinned in `environment.yml` (where
applicable), the upstream license, and the attribution requirement expected
by the upstream maintainer.

If you redistribute or modify this package, you remain responsible for
complying with these third-party terms.

---

## 1. UMA-s-1p1 (Universal Models for Atoms, Meta FAIR)

- **Used in**: Panels c, d, e and Supplementary Figs. 3 (UMA column), 4.
- **Source**: `facebook/UMA` Hugging Face repository,
  https://huggingface.co/facebook/UMA
- **File**: `uma-s-1p1.pt`, MD5 `36a2f071350be0ee4c15e7ebdd16dde1`, snapshot
  `38529caa`.
- **Loaded via**: `fairchem-core==2.13.0` (`FAIRChemCalculator`).
- **License**: governed by the terms posted on the upstream Hugging Face
  page; refer to that page for the current license text and acceptable-use
  policy. Acknowledgement of the use of the UMA model is mandatory.
- **Citation (required)**:
  Wood, B. M. et al. *UMA: A Family of Universal Models for Atoms.*
  arXiv:2506.23971 (2025).
- **Provisioning**: not redistributed in this repository. See
  `UMA_Checkpoint_Setup.md` for the offline-cache procedure.

## 2. CHGNet v0.4.0

- **Used in**: Supplementary Fig. 3 (CHGNet column, three-method consensus
  on Pt-Ga binary intermetallics, only).
- **Source**: `chgnet==0.4.0` on PyPI; underlying weights bundled with
  the package and downloaded automatically on first use.
- **License**: BSD-3-Clause (refer to the upstream repository for the
  authoritative copy).
- **Citation (required)**:
  Deng, B. et al. *CHGNet as a pretrained universal neural network potential
  for charge-informed atomistic modelling.* Nat. Mach. Intell. 5, 1031-1041
  (2023).
- **Provisioning**: installed via `environment.yml` (`pip` section).

## 3. Materials Project DFT (PAW-PBE) reference energies

- **Used in**: Supplementary Fig. 3 (MP-DFT column) for cross-validating
  Pt-Ga binary formation enthalpies; also used implicitly as the reference
  standard against which UMA single-point energies are benchmarked.
- **Source**: Materials Project (https://materialsproject.org), accessed
  through the public REST API (`mp_api.client`).
- **Data license**: CC-BY-4.0; redistribution requires attribution.
- **Citations (required)**:
  - Jain, A. et al. *Commentary: The Materials Project: A materials genome
    approach to accelerating materials innovation.* APL Mater. 1, 011002
    (2013).
  - Cite the specific MP entries by their `mp-XXXXX` IDs as listed in the
    bundled panel CSVs and notes where applicable.
- **Provisioning**: queried at runtime via `mp_api`. Reviewers without an
  API key can use the bundled CSV snapshots in `*_data*` files to reproduce
  downstream analysis without re-querying the REST endpoint.

## 4. Miedema atomic-parameter table

- **Used in**: Panels a, b, f, Supplementary Fig. 2.
- **File**: `shared/data_periodic_table.py`.
- **Underlying sources**:
  - Takeuchi, A. & Inoue, A. *Classification of bulk metallic glasses by
    atomic size difference, heat of mixing and period of constituent
    elements and its application to characterization of the main alloying
    element.* Mater. Trans. 46, 2817-2829 (2005).
  - de Boer, F. R., Boom, R., Mattens, W. C. M., Miedema, A. R. & Niessen,
    A. K. *Cohesion in Metals: Transition Metal Alloys.* North-Holland,
    Amsterdam (1988).
  - Miedema, A. R., de Chatel, P. F. & de Boer, F. R. *Cohesion in alloys —
    fundamentals of a semi-empirical model.* Physica B 100, 1-28 (1980).
- **License of the parameter values**: the numerical parameters published
  in the cited references are facts not subject to copyright in most
  jurisdictions; the python representation in this package is released
  under the MIT license of this repository.
- **Provisioning**: bundled in this repository.

## 5. SGTE pure-element thermodynamic data

- **Used in**: Panel d (Gibbs curves vs T) for SER references at finite
  temperature.
- **Source**:
  Dinsdale, A. T. *SGTE data for pure elements.* Calphad 15, 317-425
  (1991).
- **Provisioning**: relevant SER values are hard-coded in the panel
  scripts (no external file dependency).

## 6. Demeter (Athena/Artemis/Hephaestus)

- **Used in**: EXAFS fitting (described in the Methods of the manuscript;
  not part of the python pipeline in this repository).
- **Citation (when reusing manuscript-reported EXAFS results)**:
  Ravel, B. & Newville, M. *ATHENA, ARTEMIS, HEPHAESTUS: data analysis for
  X-ray absorption spectroscopy using IFEFFIT.* J. Synchrotron Rad. 12,
  537-541 (2005).

---

## What this repository does NOT redistribute

- The UMA-s-1p1 checkpoint binary.
- Any Materials Project DFT bulk dataset; only mp-IDs and queried scalar
  values are included for the Pt-Ga validation panel.
- Origin (.opju) figure templates and Word document drafts.
- Internal review memos.

If a reviewer needs material from this list for evaluation, please contact
the corresponding author.
