from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import originpro as op

from origin_sciplot.origin_backend.export_utils import export_graph


SOURCE_DIR = Path(
    r"C:\Users\13102\ns_review_4dstem\DEEPSEEK_FINAL_v2_1"
    r"\FIG1C_FFT_RAW_ORIGIN_v1\origin_exports_magma_axes_FINAL_v7\C3_L16"
)
OUTPUT_DIR = Path(
    r"C:\Users\13102\ns_review_4dstem\DEEPSEEK_FINAL_v2_1"
    r"\FIG1C_FFT_RAW_ORIGIN_v1\C3_true_numeric_axes_FINAL_v8"
)


def array_hash(values: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(values, dtype=np.float64)
    return hashlib.sha256(contiguous.tobytes()).hexdigest()


def axis_state(layer, axis: str) -> dict[str, float | int | str]:
    return {
        "from": float(layer.get_float(f"{axis}.from")),
        "to": float(layer.get_float(f"{axis}.to")),
        "increment": float(layer.get_float(f"{axis}.inc")),
        "first_tick": float(layer.get_float(f"{axis}.firstTick")),
        "ticks": int(layer.get_int(f"{axis}.ticks")),
        "minor_ticks": int(layer.get_int(f"{axis}.minorTicks")),
        "major_length": float(layer.get_float(f"{axis}.ticklength")),
        "minor_length": float(layer.get_float(f"{axis}.mticklength")),
        "label_type": int(layer.get_int(f"{axis}.label.type")),
        "ticks_by_data": op.get_lt_str(f"layer.{axis}.ticksbydata$"),
        "minor_ticks_by_data": op.get_lt_str(f"layer.{axis}.mticksbydata$"),
    }


if OUTPUT_DIR.exists():
    raise FileExistsError(f"Refusing to overwrite existing directory: {OUTPUT_DIR}")
OUTPUT_DIR.mkdir(parents=True)
base_copy = OUTPUT_DIR / "base_v7_copy.opju"
shutil.copy2(SOURCE_DIR / "result.opju", base_copy)
shutil.copy2(SOURCE_DIR / "input_copy.csv", OUTPUT_DIR / "input_copy.csv")
shutil.copy2(SOURCE_DIR / "render-plan.json", OUTPUT_DIR / "source_render-plan.json")
shutil.copy2(
    SOURCE_DIR / "origin_verify_report.json",
    OUTPUT_DIR / "source_v7_origin_verify_report.json",
)

frame = pd.read_csv(SOURCE_DIR / "input_copy.csv")
qx = np.asarray([float(value) for value in frame.columns[1:]], dtype=float)
qy = frame.iloc[:, 0].to_numpy(dtype=float)
if not np.all(np.diff(qx) > 0):
    raise ValueError("qx coordinates must be strictly increasing")
if not np.all(np.diff(qy) < 0):
    raise ValueError("qy coordinates must be strictly decreasing")

major_values = (-10, -5, 0, 5, 10)
minor_values = tuple(
    value
    for value in range(-12, 13)
    if value not in major_values and qx.min() < value < qx.max()
)
major_text = " ".join(str(value) for value in major_values)
minor_text = " ".join(str(value) for value in minor_values)

op.set_show(False)
try:
    if not op.open(str(base_copy), readonly=False, asksave=False):
        raise RuntimeError("Origin could not open the saved C3 v7 project")

    graph = op.find_graph(0)
    matrix = op.find_sheet("m", 0)
    if graph is None or matrix is None:
        raise RuntimeError("C3 v7 graph or matrix sheet is missing")
    layer = graph[0]

    values_before = matrix.to_np2d().copy()
    hash_before = array_hash(values_before)
    xymap_before = tuple(float(value) for value in matrix.xymap)

    # Rows in the source table run from positive qy to negative qy.
    matrix.xymap = (float(qx[0]), float(qx[-1]), float(qy[0]), float(qy[-1]))
    layer.axis("x").set_limits(float(qx[0]), float(qx[-1]), 5.0)
    layer.axis("y").set_limits(float(qy[-1]), float(qy[0]), 5.0)

    for axis in ("x", "y"):
        layer.set_int(f"{axis}.label.type", 1)
        layer.set_int(f"{axis}.label.numFormat", 1)
        layer.set_int(f"{axis}.label.decimalPlaces", 0)
        layer.set_float(f"{axis}.label.rotate", 0.0)
        layer.set_int(f"{axis}.ticks", 10)  # major out (2) + minor out (8)
        layer.set_int(f"{axis}.minorTicks", 4)
        layer.set_float(f"{axis}.ticklength", 5.5)
        layer.set_float(f"{axis}.mticklength", 3.0)
        layer.set_float(f"{axis}.tickthickness", 1.8)
        layer.set_float(f"{axis}.mtickthickness", 1.2)
        layer.set_float(f"{axis}.firstTick", -10.0)
        layer.set_int(f"{axis}.showLabels", 1)
        layer.set_int(f"{axis}.label.show", 1)

    if not layer.obj.LT_execute(
        f'layer.x.ticksbydata$="{major_text}";'
        f'layer.y.ticksbydata$="{major_text}";'
        f'layer.x.mticksbydata$="{minor_text}";'
        f'layer.y.mticksbydata$="{minor_text}";'
        "doc -uw;"
    ):
        raise RuntimeError("Origin could not apply exact numeric tick positions")

    values_after = matrix.to_np2d().copy()
    hash_after = array_hash(values_after)
    if hash_after != hash_before or not np.array_equal(values_after, values_before, equal_nan=True):
        raise RuntimeError("Matrix intensity values changed while setting coordinate metadata")

    graph.activate()
    op.set_show(True)
    op.lt_exec("doc -uw;")

    output_opju = OUTPUT_DIR / "result.opju"
    if not op.save(str(output_opju)):
        raise RuntimeError("Origin could not save the C3 v8 project")
    exports = export_graph(
        op,
        graph,
        OUTPUT_DIR / "result.png",
        OUTPUT_DIR / "result.pdf",
        OUTPUT_DIR / "result.tif",
    )

    report = {
        "source_project": str(SOURCE_DIR / "result.opju"),
        "source_table": str(SOURCE_DIR / "input_copy.csv"),
        "output_project": str(output_opju),
        "matrix_shape": list(values_after.shape),
        "matrix_hash_before": hash_before,
        "matrix_hash_after": hash_after,
        "matrix_values_modified": False,
        "xymap_before": list(xymap_before),
        "xymap_after": [float(value) for value in matrix.xymap],
        "major_tick_values": list(major_values),
        "minor_tick_values": list(minor_values),
        "tick_direction": "out",
        "x_axis": axis_state(layer, "x"),
        "y_axis": axis_state(layer, "y"),
        "exports": exports,
    }
    (OUTPUT_DIR / "numeric_axis_readback.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
finally:
    op.exit()
