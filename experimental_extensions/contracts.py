"""Validation helpers for bounded experimental-extension manifests."""
from __future__ import annotations

import json
from pathlib import Path


SUPPORTED_MODULE_KINDS = {"manifold_regression", "energy_backend_comparison"}

SAFE_CLAIM_BOUNDARY = {
    "software_configurability": "demonstrated_for_bundled_inputs",
    "within_manifold_interpolation": "validated_internal",
    "binary_energy_rank_consistency": "observed_reference_set",
    "new_host_or_prototype_transferability": "not_evaluated",
    "nonmetal_compound_transferability": "not_evaluated",
    "synthesizability": "not_predicted",
}


def resolve_relative(owner_path: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (owner_path.parent / path).resolve()


def read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def validate_manifest(manifest: dict) -> None:
    required = {
        "schema_version", "system_id", "status", "scientific_scope",
        "output_dir", "contract_catalog", "claim_boundary", "modules",
    }
    missing = sorted(required - set(manifest))
    if missing:
        raise ValueError(f"Manifest is missing required fields: {missing}")
    if manifest["schema_version"] != 1:
        raise ValueError("Only system manifest schema_version=1 is supported")
    if manifest["status"] != "bounded-validation-demo":
        raise ValueError("Manifest status must remain 'bounded-validation-demo'")
    if manifest["claim_boundary"] != SAFE_CLAIM_BOUNDARY:
        raise ValueError(
            "claim_boundary must use the bounded validation states defined for schema v1"
        )
    modules = manifest["modules"]
    if not isinstance(modules, list) or not modules:
        raise ValueError("Manifest modules must be a non-empty list")
    seen: set[str] = set()
    for module in modules:
        module_id = module.get("id")
        kind = module.get("kind")
        if not module_id or module_id in seen:
            raise ValueError(f"Module ids must be non-empty and unique: {module_id!r}")
        seen.add(module_id)
        if kind not in SUPPORTED_MODULE_KINDS:
            raise ValueError(f"Unsupported module kind: {kind!r}")
        if not isinstance(module.get("enabled"), bool):
            raise ValueError(f"Module {module_id!r} must define boolean enabled")
        if kind == "manifold_regression" and "config_path" not in module:
            raise ValueError(f"Module {module_id!r} requires config_path")
        if kind == "energy_backend_comparison" and "config" not in module:
            raise ValueError(f"Module {module_id!r} requires config")


def load_manifest(path: Path) -> dict:
    manifest = read_json(path)
    validate_manifest(manifest)
    return manifest


def load_contract_catalog(path: Path, manifest: dict) -> dict:
    catalog = read_json(resolve_relative(path, manifest["contract_catalog"]))
    if catalog.get("schema_version") != 1:
        raise ValueError("Only module-contract schema_version=1 is supported")
    contracts = catalog.get("contracts")
    if not isinstance(contracts, list):
        raise ValueError("Contract catalog must contain a contracts list")
    available = {item.get("kind") for item in contracts}
    required = {item["kind"] for item in manifest["modules"] if item["enabled"]}
    missing = sorted(required - available)
    if missing:
        raise ValueError(f"Missing module contracts: {missing}")
    return catalog
