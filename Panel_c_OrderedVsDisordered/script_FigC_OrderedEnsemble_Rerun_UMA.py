"""Generate and score 30 symmetry-distinct ordered B-sublattice structures.

The historical Panel-c ordered structure is retained as configuration 0. Twenty-
nine additional Pt24Ga2In2Sn2Zn2 structures are selected deterministically from
different symmetry classes of the 2 x 2 x 2 L1_2 B-sublattice assignment space.
Canonical manuscript files are never overwritten.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = SCRIPT_DIR / "data_FigC_OrderedEnsemble_Raw_UMA_Energies_regen.csv"
DEFAULT_MANIFEST = SCRIPT_DIR / "data_FigC_OrderedEnsemble_StructureManifest_regen.csv"
DEFAULT_REPORT = SCRIPT_DIR / "validation_FigC_OrderedEnsemble_UMA_regen.json"

LATTICE_CONSTANT_ANGSTROM = 3.903
N_ORDERED = 30
SELECTION_SEED = 20260727
EXPECTED_COMPOSITION = Counter({"Pt": 24, "Ga": 2, "In": 2, "Sn": 2, "Zn": 2})
EXPECTED_RAW_ASSIGNMENTS = 2520
EXPECTED_SYMMETRY_CLASSES = 68
EXPECTED_CHECKPOINT_MD5 = "36a2f071350be0ee4c15e7ebdd16dde1"
DEFAULT_SNAPSHOT = "38529caa2c51a9a8a0d71f0b56b79ac33bc9eceb"


def build_parent_supercell():
    from ase.build import bulk, make_supercell

    parent = bulk("Pt", "fcc", a=LATTICE_CONSTANT_ANGSTROM, cubic=True)
    parent[0].symbol = "Ga"
    atoms = make_supercell(parent, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    b_indices = tuple(index for index, atom in enumerate(atoms) if atom.symbol == "Ga")
    if len(atoms) != 32 or len(b_indices) != 8:
        raise ValueError(f"Unexpected parent: N={len(atoms)}, B-sites={len(b_indices)}")
    return atoms, b_indices


def enumerate_assignments():
    """Yield every raw Ga2In2Sn2Zn2 assignment on eight labelled B sites."""

    positions = set(range(8))
    for ga_sites in itertools.combinations(range(8), 2):
        remaining_1 = positions - set(ga_sites)
        for in_sites in itertools.combinations(sorted(remaining_1), 2):
            remaining_2 = remaining_1 - set(in_sites)
            for sn_sites in itertools.combinations(sorted(remaining_2), 2):
                zn_sites = remaining_2 - set(sn_sites)
                assignment = [""] * 8
                for sites, symbol in (
                    (ga_sites, "Ga"),
                    (in_sites, "In"),
                    (sn_sites, "Sn"),
                    (zn_sites, "Zn"),
                ):
                    for site in sites:
                        assignment[site] = symbol
                yield tuple(assignment)


def b_sublattice_permutations(parent, b_indices, symprec: float = 1.0e-5):
    import spglib

    lattice = np.asarray(parent.cell.array, dtype=float)
    fractional = np.asarray(parent.get_scaled_positions(), dtype=float)
    numbers = np.asarray(parent.numbers, dtype=int)
    symmetry = spglib.get_symmetry((lattice, fractional, numbers), symprec=symprec)
    if symmetry is None:
        raise RuntimeError("spglib did not return symmetry operations")

    b_fractional = fractional[list(b_indices)]
    permutations: set[tuple[int, ...]] = set()
    for rotation, translation in zip(symmetry["rotations"], symmetry["translations"]):
        permutation: list[int] = []
        for position in b_fractional:
            target = (rotation @ position + translation) % 1.0
            delta = b_fractional - target
            delta -= np.rint(delta)
            distances = np.linalg.norm(delta @ lattice, axis=1)
            mapped = int(np.argmin(distances))
            if float(distances[mapped]) > symprec:
                raise RuntimeError("Failed to map a B site under a parent symmetry operation")
            permutation.append(mapped)
        if len(set(permutation)) != 8:
            raise RuntimeError("A symmetry operation did not produce a B-site permutation")
        permutations.add(tuple(permutation))
    return sorted(permutations), int(len(symmetry["rotations"]))


def transform_assignment(assignment, permutation):
    transformed = [""] * 8
    for old_index, new_index in enumerate(permutation):
        transformed[new_index] = assignment[old_index]
    return tuple(transformed)


def canonical_assignment(assignment, permutations):
    return min(transform_assignment(assignment, permutation) for permutation in permutations)


def historical_assignment(b_indices):
    shuffled = list(b_indices)
    random.Random(0).shuffle(shuffled)
    by_atom_index = {index: "Ga" for index in b_indices}
    for index, symbol in zip(shuffled[:6], ("In", "In", "Sn", "Sn", "Zn", "Zn")):
        by_atom_index[index] = symbol
    return tuple(by_atom_index[index] for index in b_indices)


def select_ordered_assignments(
    n_ordered: int = N_ORDERED,
    selection_seed: int = SELECTION_SEED,
):
    parent, b_indices = build_parent_supercell()
    permutations, n_spglib_operations = b_sublattice_permutations(parent, b_indices)
    classes: dict[tuple[str, ...], list[tuple[str, ...]]] = defaultdict(list)
    for assignment in enumerate_assignments():
        classes[canonical_assignment(assignment, permutations)].append(assignment)

    n_raw = sum(len(members) for members in classes.values())
    if n_raw != EXPECTED_RAW_ASSIGNMENTS:
        raise AssertionError(f"Expected {EXPECTED_RAW_ASSIGNMENTS} raw assignments, got {n_raw}")
    if len(classes) != EXPECTED_SYMMETRY_CLASSES:
        raise AssertionError(
            f"Expected {EXPECTED_SYMMETRY_CLASSES} symmetry classes, got {len(classes)}"
        )

    historical = historical_assignment(b_indices)
    historical_class = canonical_assignment(historical, permutations)
    if not 1 <= n_ordered <= len(classes):
        raise ValueError(f"n_ordered must be between 1 and {len(classes)}, got {n_ordered}")

    remaining_classes = sorted(key for key in classes if key != historical_class)
    if n_ordered == len(classes):
        selected_classes = remaining_classes
        selection_mode = "all_symmetry_classes"
    else:
        selected_classes = random.Random(selection_seed).sample(
            remaining_classes, n_ordered - 1
        )
        selection_mode = "seeded_symmetry_class_sample"

    selected = [(historical_class, historical)]
    for class_key in selected_classes:
        selected.append((class_key, min(classes[class_key])))

    rows = []
    for config_id, (class_key, assignment) in enumerate(selected):
        atoms = parent.copy()
        for site_index, symbol in zip(b_indices, assignment):
            atoms[site_index].symbol = symbol
        validate_ordered_structure(atoms, b_indices)
        fingerprint = "|".join(assignment)
        rows.append(
            {
                "Type": "Ordered_L12_BSublattice",
                "Config_ID": config_id,
                "Selection_Seed": 0 if config_id == 0 else selection_seed,
                "Selection_Mode": selection_mode,
                "Is_Historical_Anchor": config_id == 0,
                "Assignment_B0_to_B7": fingerprint,
                "Assignment_SHA256": hashlib.sha256(fingerprint.encode("ascii")).hexdigest(),
                "Canonical_Class": "|".join(class_key),
                "Class_Degeneracy": len(classes[class_key]),
                "atoms": atoms,
            }
        )

    if len({row["Canonical_Class"] for row in rows}) != n_ordered:
        raise AssertionError("Selected ordered configurations are not symmetry-distinct")
    return rows, {
        "raw_assignment_count": n_raw,
        "symmetry_class_count": len(classes),
        "spglib_operation_count": n_spglib_operations,
        "unique_b_sublattice_permutations": len(permutations),
    }


def validate_ordered_structure(atoms, b_indices) -> None:
    observed = Counter(atoms.get_chemical_symbols())
    if len(atoms) != 32 or observed != EXPECTED_COMPOSITION:
        raise ValueError(f"Unexpected ordered structure: N={len(atoms)}, composition={observed}")
    b_set = set(b_indices)
    if any(atom.symbol != "Pt" for index, atom in enumerate(atoms) if index not in b_set):
        raise ValueError("An ordered structure places a non-Pt species on the A sublattice")
    if any(atoms[index].symbol == "Pt" for index in b_indices):
        raise ValueError("An ordered structure places Pt on the B sublattice")


def export_structures(rows, output_dir: Path) -> None:
    from ase.io import write

    output_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        write(output_dir / f"ordered_config_{int(row['Config_ID']):02d}.cif", row["atoms"])


def resolve_checkpoint(explicit: Path | None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit.expanduser())
    if os.environ.get("UMA_CHECKPOINT"):
        candidates.append(Path(os.environ["UMA_CHECKPOINT"]).expanduser())
    candidates.append(
        Path.home()
        / ".cache"
        / "fairchem"
        / "models--facebook--UMA"
        / "snapshots"
        / DEFAULT_SNAPSHOT
        / "checkpoints"
        / "uma-s-1p1.pt"
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    checked = "\n  - ".join(str(path) for path in candidates)
    raise FileNotFoundError("UMA checkpoint not found. Checked:\n  - " + checked)


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_calculator(checkpoint: Path, device: str):
    from fairchem.core.calculate.ase_calculator import FAIRChemCalculator

    return FAIRChemCalculator.from_model_checkpoint(
        str(checkpoint), task_name="oc20", device=device
    )


def score_structures(rows, calculator) -> None:
    for row in rows:
        atoms = row["atoms"]
        atoms.calc = calculator
        row["Energy_eV_atom"] = float(atoms.get_potential_energy() / len(atoms))
        print(
            f"ordered config={int(row['Config_ID']):02d} "
            f"degeneracy={int(row['Class_Degeneracy']):02d} "
            f"E={float(row['Energy_eV_atom']):+.9f} eV/atom"
        )


def write_manifest(rows, output: Path) -> None:
    fields = [
        "Type",
        "Config_ID",
        "Selection_Seed",
        "Selection_Mode",
        "Is_Historical_Anchor",
        "Assignment_B0_to_B7",
        "Assignment_SHA256",
        "Canonical_Class",
        "Class_Degeneracy",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_energies(rows, output: Path) -> None:
    fields = [
        "Type",
        "Config_ID",
        "Selection_Seed",
        "Selection_Mode",
        "Is_Historical_Anchor",
        "Assignment_B0_to_B7",
        "Assignment_SHA256",
        "Canonical_Class",
        "Class_Degeneracy",
        "Energy_eV_atom",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            output_row = {field: row[field] for field in fields}
            output_row["Energy_eV_atom"] = f"{float(row['Energy_eV_atom']):.9f}"
            writer.writerow(output_row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--device", default="auto", choices=("auto", "cpu", "cuda"))
    parser.add_argument("--n-ordered", type=int, default=N_ORDERED)
    parser.add_argument("--selection-seed", type=int, default=SELECTION_SEED)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--export-cifs", type=Path)
    parser.add_argument("--structures-only", action="store_true")
    parser.add_argument("--skip-checkpoint-md5", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows, enumeration = select_ordered_assignments(
        n_ordered=args.n_ordered,
        selection_seed=args.selection_seed,
    )
    write_manifest(rows, args.manifest)
    if args.export_cifs:
        export_structures(rows, args.export_cifs)
    if args.structures_only:
        if not args.export_cifs:
            raise ValueError("--structures-only requires --export-cifs")
        print(json.dumps({"n_structures": len(rows), **enumeration}, indent=2))
        return 0

    checkpoint = resolve_checkpoint(args.checkpoint)
    checkpoint_md5 = EXPECTED_CHECKPOINT_MD5 if args.skip_checkpoint_md5 else file_md5(checkpoint)
    if checkpoint_md5 != EXPECTED_CHECKPOINT_MD5:
        raise ValueError(
            f"Checkpoint MD5 mismatch: {checkpoint_md5}, expected {EXPECTED_CHECKPOINT_MD5}"
        )

    if args.device == "auto":
        import torch

        device = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        device = args.device
    calculator = load_calculator(checkpoint, device)
    score_structures(rows, calculator)
    write_energies(rows, args.output)

    report = {
        "status": "passed",
        "provenance_mode": "measured_computational",
        "model": "UMA-s-1p1",
        "task_name": "oc20",
        "checkpoint_file": checkpoint.name,
        "checkpoint_md5": checkpoint_md5,
        "device": device,
        "selection_seed": args.selection_seed,
        "selection_mode": rows[0]["Selection_Mode"],
        "n_ordered": len(rows),
        "historical_anchor_energy_eV_atom": float(rows[0]["Energy_eV_atom"]),
        **enumeration,
        "energy_file": args.output.name,
        "structure_manifest": args.manifest.name,
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
