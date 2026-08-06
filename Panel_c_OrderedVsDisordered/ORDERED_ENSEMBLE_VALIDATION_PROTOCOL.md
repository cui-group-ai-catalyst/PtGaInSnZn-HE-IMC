# Ordered B-sublattice ensemble validation protocol

Date frozen: 2026-07-27

## Research question

At fixed Pt24Ga2In2Sn2Zn2 composition and fixed 2 x 2 x 2 fcc cell, does the
energy advantage previously reported for one L1_2-ordered representative remain
when multiple distinct Ga/In/Sn/Zn arrangements are sampled on the eight B sites?

## Falsifiable hypothesis

Thirty symmetry-distinct B-sublattice arrangements, with Pt fixed on all 24 A
sites, will remain lower in UMA-s-1p1 single-point energy than the existing 30
full-lattice random configurations. The conclusion fails if the two energy
distributions overlap materially or if the mean ordering gap is not positive.

## Data contract

- Composition: Pt24Ga2In2Sn2Zn2 for every structure.
- Parent cell: ASE cubic fcc Pt conventional cell, a = 3.903 Angstrom.
- Supercell: diagonal 2 x 2 x 2, 32 atoms.
- Ordered population: Pt fixed on A sites; Ga/In/Sn/Zn restricted to B sites.
- Raw ordered assignment space: 8!/(2!^4) = 2520 assignments.
- Symmetry treatment: spglib operations of the Pt24Ga8 parent identify 68
  equivalence classes. The canonical historical ordered structure is retained,
  and 29 additional classes are selected without replacement using selection
  seed 20260727.
- Disordered comparator: the existing 30 full-lattice random occupations,
  seeds 100-129, from data_FigC_Raw_UMA_Energies.csv.
- Model: frozen UMA-s-1p1 checkpoint, OC20 head; no system-specific training.
- Energy protocol: fixed-cell, fixed-coordinate single-point energy.
- Element reference: the same five UMA elemental references in
  ../shared/UMA_Element_Reference_Energies.csv.
- Provenance mode: measured computational output from the locally installed
  checkpoint, not simulated or illustrative data.

## Metrics fixed before execution

1. Ordered and disordered mean formation enthalpy.
2. Sample standard deviation and standard error for each set.
3. Difference of means, disordered minus ordered.
4. Welch 95% confidence interval for the difference of means.
5. Non-overlap margin: minimum disordered energy minus maximum ordered energy.
6. Standardized mean difference (pooled-SD Cohen d).

No structure may be removed on the basis of its calculated energy. A structure
may be rejected only before scoring for invalid composition, A/B-site violation,
duplicate symmetry class, or duplicate occupancy fingerprint.

## Interpretation boundary

This is a robustness test inside a fixed L1_2/A1-like, fixed-composition,
fixed-cell manifold. It does not establish the global ground state, compare all
competing phases, include vibrational free energies, or demonstrate synthesis
kinetics. Symmetry-distinct motifs are sampled uniformly for robustness; the
result is not a degeneracy-weighted canonical thermodynamic ensemble.

