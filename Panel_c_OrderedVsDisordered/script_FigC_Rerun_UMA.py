"""Rebuild the 31 fixed-composition Panel-c UMA single-point energies.

The ordered structure and the 30 disordered structures are generated
deterministically from the manuscript lattice constant and seeds 100-129.
Canonical manuscript files are never overwritten; regenerated files use the
``_regen`` suffix unless an explicit output path is supplied.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
from collections import Counter
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
CANONICAL_RAW = SCRIPT_DIR / "data_FigC_Raw_UMA_Energies.csv"
DEFAULT_OUTPUT = SCRIPT_DIR / "data_FigC_Raw_UMA_Energies_regen.csv"
DEFAULT_REPORT = SCRIPT_DIR / "validation_FigC_UMA_Rerun.json"

LATTICE_CONSTANT_ANGSTROM = 3.903
FIRST_DISORDERED_SEED = 100
N_DISORDERED = 30
EXPECTED_COMPOSITION = Counter({"Pt": 24, "Ga": 2, "In": 2, "Sn": 2, "Zn": 2})
DEFAULT_SNAPSHOT = "38529caa2c51a9a8a0d71f0b56b79ac33bc9eceb"


def build_ordered_equimolar():
    from ase.build import bulk, make_supercell

    primitive = bulk("Pt", "fcc", a=LATTICE_CONSTANT_ANGSTROM, cubic=True)
    primitive[0].symbol = "Ga"
    atoms = make_supercell(primitive, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    b_indices = [index for index, atom in enumerate(atoms) if atom.symbol == "Ga"]
    if len(b_indices) != 8:
        raise ValueError(f"Expected 8 B-sublattice sites, found {len(b_indices)}")

    rng = random.Random(0)
    rng.shuffle(b_indices)
    for index, symbol in zip(b_indices[:6], ["In", "In", "Sn", "Sn", "Zn", "Zn"]):
        atoms[index].symbol = symbol
    _validate_structure(atoms)
    return atoms


def build_disordered_equimolar(seed: int):
    from ase.build import bulk, make_supercell

    primitive = bulk("Pt", "fcc", a=LATTICE_CONSTANT_ANGSTROM, cubic=True)
    atoms = make_supercell(primitive, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    species = ["Pt"] * 24 + ["Ga"] * 2 + ["In"] * 2 + ["Sn"] * 2 + ["Zn"] * 2
    rng = random.Random(seed)
    rng.shuffle(species)
    for atom, symbol in zip(atoms, species):
        atom.symbol = symbol
    _validate_structure(atoms)
    return atoms


def _validate_structure(atoms) -> None:
    observed = Counter(atoms.get_chemical_symbols())
    if len(atoms) != 32 or observed != EXPECTED_COMPOSITION:
        raise ValueError(f"Unexpected Panel-c structure: N={len(atoms)}, composition={observed}")


def generate_structures() -> list[tuple[str, int, str, object]]:
    structures = [("Ordered_L12", 0, "N/A", build_ordered_equimolar())]
    fingerprints: set[tuple[str, ...]] = set()
    for config_id, seed in enumerate(
        range(FIRST_DISORDERED_SEED, FIRST_DISORDERED_SEED + N_DISORDERED), start=1
    ):
        atoms = build_disordered_equimolar(seed)
        fingerprint = tuple(atoms.get_chemical_symbols())
        if fingerprint in fingerprints:
            raise ValueError(f"Duplicate disordered occupation generated for seed {seed}")
        fingerprints.add(fingerprint)
        structures.append(("Disordered_Random", config_id, str(seed), atoms))
    return structures


def export_structures(structures, output_dir: Path) -> None:
    from ase.io import write

    output_dir.mkdir(parents=True, exist_ok=True)
    for row_type, _, seed, atoms in structures:
        name = "ordered_l12.cif" if row_type == "Ordered_L12" else f"disordered_seed_{seed}.cif"
        write(output_dir / name, atoms)


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
    raise FileNotFoundError(
        "UMA checkpoint not found. Pass --checkpoint or set UMA_CHECKPOINT. Checked:\n  - "
        + checked
    )


def load_calculator(checkpoint: Path, device: str):
    from fairchem.core.calculate.ase_calculator import FAIRChemCalculator

    return FAIRChemCalculator.from_model_checkpoint(
        str(checkpoint), task_name="oc20", device=device
    )


def score_structures(structures, calculator) -> list[dict[str, str | int | float]]:
    rows: list[dict[str, str | int | float]] = []
    for row_type, config_id, seed, atoms in structures:
        atoms.calc = calculator
        energy = float(atoms.get_potential_energy() / len(atoms))
        rows.append(
            {
                "Type": row_type,
                "Config_ID": config_id,
                "Seed": seed,
                "Energy_eV_atom": energy,
            }
        )
        print(f"{row_type:20s} config={config_id:02d} seed={seed:>3s} E={energy:+.9f} eV/atom")
    return rows


def write_rows(rows, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Type", "Config_ID", "Seed", "Energy_eV_atom"])
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "Energy_eV_atom": f"{float(row['Energy_eV_atom']):.9f}"})


def validate_against_canonical(rows, tolerance: float) -> dict:
    with CANONICAL_RAW.open(newline="", encoding="utf-8") as handle:
        canonical = list(csv.DictReader(handle))
    if len(rows) != len(canonical):
        raise AssertionError(f"Row-count mismatch: regenerated={len(rows)}, canonical={len(canonical)}")

    differences = []
    for regenerated, expected in zip(rows, canonical):
        for key in ("Type", "Config_ID", "Seed"):
            if str(regenerated[key]) != expected[key]:
                raise AssertionError(
                    f"Panel-c identity mismatch at config {expected['Config_ID']}: "
                    f"{key}={regenerated[key]!r}, expected {expected[key]!r}"
                )
        differences.append(abs(float(regenerated["Energy_eV_atom"]) - float(expected["Energy_eV_atom"])))

    max_difference = max(differences)
    return {
        "status": "passed" if max_difference <= tolerance else "failed",
        "canonical_file": CANONICAL_RAW.name,
        "n_structures": len(rows),
        "ordered_structures": 1,
        "disordered_structures": N_DISORDERED,
        "disordered_seed_range": [FIRST_DISORDERED_SEED, FIRST_DISORDERED_SEED + N_DISORDERED - 1],
        "max_abs_energy_difference_eV_atom": max_difference,
        "tolerance_eV_atom": tolerance,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, help="Path to uma-s-1p1.pt")
    parser.add_argument("--device", default=os.environ.get("UMA_DEVICE", "cpu"))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--export-cifs", type=Path, help="Optional directory for the 31 generated CIFs")
    parser.add_argument(
        "--structures-only",
        action="store_true",
        help="Generate/export structures without loading UMA; requires --export-cifs",
    )
    parser.add_argument("--tolerance", type=float, default=1.0e-5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    structures = generate_structures()
    if args.export_cifs:
        export_structures(structures, args.export_cifs)
        print(f"Exported {len(structures)} structures to {args.export_cifs}")
    if args.structures_only:
        if not args.export_cifs:
            raise ValueError("--structures-only requires --export-cifs")
        return 0

    checkpoint = resolve_checkpoint(args.checkpoint)
    calculator = load_calculator(checkpoint, args.device)
    rows = score_structures(structures, calculator)
    write_rows(rows, args.output)
    report = validate_against_canonical(rows, args.tolerance)
    report.update(
        {
            "model": "UMA-s-1p1",
            "checkpoint_file": checkpoint.name,
            "device": args.device,
            "regenerated_file": args.output.name,
        }
    )
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if report["status"] != "passed":
        raise AssertionError("Panel-c UMA regeneration exceeded the configured tolerance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
