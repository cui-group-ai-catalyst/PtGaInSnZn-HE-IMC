# Element reference structures

These 5 CIF files are the conventional unit cells used as the **stable element
reference (SER)** when computing formation enthalpies in this study:

| Element | Structure | Source |
|---|---|---|
| Pt | fcc | Materials Project entry `mp-126` |
| Ga | alpha-Ga | Materials Project entry `mp-142` |
| In | bct (I4/mmm), a=3.25, c=4.95 Angstrom | Generated conventional cell |
| Sn | alpha-Sn (Fd-3m), a=6.4892 Angstrom | Generated conventional diamond |
| Zn | hcp, a=2.66, c=4.93696 Angstrom | Generated via ASE bulk constructor |

`element_reference_manifest.csv` lists CIF provenance.
`../UMA_Element_Reference_Energies.csv` lists the UMA-s-1p1 single-point
energies (eV/atom) used in the manuscript.

## Why these are bundled

The Nature Software Submission Checklist requires that reviewers can run the
software on their own data. The formation enthalpy values reported in
Figure 1c/d/e and Supplementary Figs. 3/4 depend on these 5 element reference
energies. By bundling the actual CIF inputs, reviewers can:

1. Run a different ML potential (CHGNet, MACE, M3GNet, DFT, etc.) on **exactly
   the same** element structures used in the manuscript, and obtain a fair
   apples-to-apples reference set.
2. Substitute their own reference set into
   `scripts/recompute_element_referenced_Hf.py` to re-derive
   `ElementRef_Hf_kJ_mol` from the bundled `Energy_eV_atom` values without
   touching the manuscript canonical results.

See README section "Using your own data" for the substitution workflow.
