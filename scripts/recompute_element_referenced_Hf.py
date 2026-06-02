"""Recompute element-referenced formation enthalpy with a user-supplied
reference energy set.

Use case
--------
A reviewer who wants to check the manuscript Pt3(Ga,In,Sn,Zn) formation
enthalpies using a different reference set -- e.g. DFT-PBE energies, a
different ML potential such as CHGNet or MACE, or a custom hand-curated
reference -- can supply their own 5 element energies and re-derive the
Element-referenced formation enthalpy column for the bundled 165-prototype
B-sublattice landscape (Supplementary Fig. 4 / Panel g).

This script does NOT touch the manuscript canonical CSV. It writes a new
side-by-side comparison CSV under shared/.

Usage
-----
1. Default mode -- reproduces the bundled column using the bundled UMA
   reference energies (sanity check that the math is the same):

       python scripts/recompute_element_referenced_Hf.py

2. User mode -- supply your own 5 element energies in a CSV file:

       python scripts/recompute_element_referenced_Hf.py \\
           --user-refs my_element_energies.csv

   The user CSV must contain at least these two columns:

       Element,E_eV_atom
       Pt,-5.4500
       Ga,-2.5300
       In,-2.1700
       Sn,-3.4500
       Zn,-0.7200

Output
------
Writes shared/ElementRef_Hf_USER_vs_BUNDLED.csv with both the bundled
ElementRef_Hf_kJ_mol column and the recomputed user_ElementRef_Hf_kJ_mol
column for all 165 prototypes, plus a per-row difference. Prints a short
summary to stdout.

Formula
-------
For each prototype on the Pt3X8 32-atom 2x2x2 supercell (24 Pt on A sublattice
plus 8 B-sublattice atoms among Ga/In/Sn/Zn):

    ref_mix = (24 * E_Pt + n_Ga * E_Ga + n_In * E_In
                          + n_Sn * E_Sn + n_Zn * E_Zn) / 32
    DeltaH_f [kJ/mol/atom] = (Energy_eV_atom - ref_mix) * 96.485

Only the user's 5 element energies are substituted; the per-prototype
Energy_eV_atom (raw alloy energy from UMA-s-1p1 at manuscript prep time)
is kept as ground truth, since substituting that would require a fresh
UMA / DFT / other-potential single-point run on each of the 165 structures.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

EV_TO_KJMOL = 96.485
ELEMENTS = ["Pt", "Ga", "In", "Sn", "Zn"]
N_PT_A_SUBLATTICE = 24
N_TOTAL_ATOMS = 32


def read_refs(path: Path, value_col: str) -> dict[str, float]:
    refs: dict[str, float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            elem = row["Element"]
            if elem in ELEMENTS:
                refs[elem] = float(row[value_col])
    missing = [e for e in ELEMENTS if e not in refs]
    if missing:
        raise SystemExit(
            f"Missing reference energies for {missing} in {path}. "
            f"Required columns: Element, {value_col}."
        )
    return refs


def main() -> int:
    release_root = Path(__file__).resolve().parents[1]
    bundled_refs_path = release_root / "shared" / "UMA_Element_Reference_Energies.csv"
    landscape_path = (
        release_root
        / "SI_Figures"
        / "SI_Fig04_165CompositionLandscape"
        / "data_FigG_165_ElementReferenced_Hf.csv"
    )
    output_path = release_root / "shared" / "ElementRef_Hf_USER_vs_BUNDLED.csv"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--user-refs",
        type=Path,
        default=None,
        help=(
            "Path to a user-supplied CSV with columns Element, E_eV_atom. "
            "If omitted, the bundled UMA reference energies are used (sanity check)."
        ),
    )
    args = parser.parse_args()

    bundled_refs = read_refs(bundled_refs_path, value_col="UMA_E_eV_atom")
    if args.user_refs is None:
        user_refs = bundled_refs
        mode = "BUNDLED (sanity check)"
    else:
        user_refs = read_refs(args.user_refs, value_col="E_eV_atom")
        mode = f"USER ({args.user_refs.name})"

    rows_out = []
    max_abs_diff = 0.0
    sum_abs_diff = 0.0
    n = 0

    with landscape_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            n_ga = int(row["Ga_count"])
            n_in = int(row["In_count"])
            n_sn = int(row["Sn_count"])
            n_zn = int(row["Zn_count"])
            e_alloy = float(row["Energy_eV_atom"])

            ref_mix = (
                N_PT_A_SUBLATTICE * user_refs["Pt"]
                + n_ga * user_refs["Ga"]
                + n_in * user_refs["In"]
                + n_sn * user_refs["Sn"]
                + n_zn * user_refs["Zn"]
            ) / N_TOTAL_ATOMS

            user_hf = (e_alloy - ref_mix) * EV_TO_KJMOL
            bundled_hf = float(row["ElementRef_Hf_kJ_mol"])
            diff = user_hf - bundled_hf

            rows_out.append(
                {
                    "Rank": row["Rank"],
                    "Composition": row["Composition"],
                    "Ga_count": n_ga,
                    "In_count": n_in,
                    "Sn_count": n_sn,
                    "Zn_count": n_zn,
                    "Energy_eV_atom": f"{e_alloy:.6f}",
                    "Bundled_ElementRef_Hf_kJ_mol": f"{bundled_hf:.6f}",
                    "User_ElementRef_Hf_kJ_mol": f"{user_hf:.6f}",
                    "Diff_kJ_mol": f"{diff:.6f}",
                }
            )

            max_abs_diff = max(max_abs_diff, abs(diff))
            sum_abs_diff += abs(diff)
            n += 1

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows_out[0].keys()))
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Mode:                {mode}")
    print(f"Element refs used:   {user_refs}")
    print(f"Bundled refs:        {bundled_refs}")
    print(f"Prototypes processed: {n}")
    print(f"Max |User - Bundled|: {max_abs_diff:.6f} kJ/mol/atom")
    print(f"Mean |User - Bundled|: {sum_abs_diff / n:.6f} kJ/mol/atom")
    print(f"Output:               {output_path}")

    if args.user_refs is None and max_abs_diff > 1e-3:
        print(
            "WARNING: Bundled-vs-bundled sanity check produced "
            f"non-trivial diff {max_abs_diff:.6f} kJ/mol/atom -- "
            "this indicates the bundled ElementRef_Hf column does not "
            "match the formula applied here. Investigate before trusting "
            "USER results."
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
