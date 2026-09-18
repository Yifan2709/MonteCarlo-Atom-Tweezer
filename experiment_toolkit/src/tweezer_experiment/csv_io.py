"""Control-table CSV I/O: user import (F2) and export tables (F5).

Canonical physical control table columns (units in the header row):
    time_us, aod_x_um, aod_y_um, aod_vx_um_per_us, aod_vy_um_per_us,
    aod_depth_uK, slm_factor
The header line ``# units: ...`` documents them; readers require exactly
these names so a mistyped column is an error, not a silent zero.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .errors import (ToolkitError, E_MISSING_FIELD, E_NOT_FINITE, E_TIME_ORDER,
                     E_BAD_VALUE, E_BAD_TYPE, E_LENGTH_MISMATCH)
from . import units as U

COLUMNS = ("time_us", "aod_x_um", "aod_y_um", "aod_vx_um_per_us",
           "aod_vy_um_per_us", "aod_depth_uK", "slm_factor")
REQUIRED_INPUT = ("time_us", "aod_x_um", "aod_y_um", "aod_depth_uK", "slm_factor")


def read_control_csv(path: str | Path, *, velocity_columns: bool = False,
                     interpolation: str = "linear") -> dict:
    """Read a user control CSV into SI arrays.

    ``velocity_columns=False`` derives velocities from the position columns
    with central differences (one-sided at the ends) — the declared rule,
    applied uniformly. Times must be strictly increasing and finite.
    """
    path = Path(path)
    if not path.is_file():
        raise ToolkitError(E_MISSING_FIELD, f"control CSV not found: {path}", str(path))
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise ToolkitError(E_BAD_TYPE, "control CSV is empty", str(path))
    header = [h.strip() for h in lines[0].split(",")]
    for col in REQUIRED_INPUT:
        if col not in header:
            raise ToolkitError(E_MISSING_FIELD, f"control CSV missing column '{col}'", str(path))
    unknown = [h for h in header if h not in COLUMNS]
    if unknown:
        raise ToolkitError(E_BAD_VALUE, f"unknown columns {unknown}; allowed {list(COLUMNS)}", str(path))
    rows = []
    for i, ln in enumerate(lines[1:], start=2):
        parts = [p.strip() for p in ln.split(",")]
        if len(parts) != len(header):
            raise ToolkitError(E_LENGTH_MISMATCH,
                               f"row {i}: {len(parts)} fields, header has {len(header)}", str(path))
        try:
            rows.append([float(p) for p in parts])
        except ValueError as exc:
            raise ToolkitError(E_BAD_TYPE, f"row {i}: non-numeric field ({exc})", str(path)) from exc
    data = np.asarray(rows, dtype=float)
    if not np.all(np.isfinite(data)):
        raise ToolkitError(E_NOT_FINITE, "non-finite value in control CSV", str(path))
    col = {name: data[:, header.index(name)].copy() for name in header}
    t_us = col["time_us"]
    if np.any(np.diff(t_us) <= 0):
        bad = int(np.argmax(np.diff(t_us) <= 0))
        raise ToolkitError(E_TIME_ORDER,
                           f"time_us must strictly increase (violation at row {bad + 2})", str(path))
    if np.any(col["aod_depth_uK"] < 0):
        raise ToolkitError(E_BAD_VALUE, "aod_depth_uK must be >= 0", str(path))
    if np.any((col["slm_factor"] < 0) | (col["slm_factor"] > 1)):
        raise ToolkitError(E_BAD_VALUE, "slm_factor must be in [0, 1]", str(path))
    vx_us, vy_us = col.get("aod_vx_um_per_us"), col.get("aod_vy_um_per_us")
    has_velocity = vx_us is not None and vy_us is not None
    if velocity_columns and not has_velocity:
        raise ToolkitError(E_MISSING_FIELD,
                           "velocity_columns=True but aod_v{x,y}_um_per_us columns absent", str(path))
    if not has_velocity:
        x, y = col["aod_x_um"], col["aod_y_um"]
        t = t_us
        vx = np.gradient(x, t)
        vy = np.gradient(y, t)
    else:
        vx, vy = vx_us, vy_us
    return {
        "time_s": U.us_to_s(t_us),
        "x_m": U.um_to_m(col["aod_x_um"]),
        "y_m": U.um_to_m(col["aod_y_um"]),
        "vx_m_per_s": U.um_to_m(vx) / U.us_to_s(1.0),
        "vy_m_per_s": U.um_to_m(vy) / U.us_to_s(1.0),
        "depth_j": U.temperature_uK_to_joule(col["aod_depth_uK"]),
        "slm_factor": col["slm_factor"],
        "interpolation": interpolation,
        "velocity_rule": "provided" if has_velocity else "central_differences",
        "source": str(path),
    }


def write_control_csv(path: str | Path, table: dict, metadata: dict | None = None) -> None:
    """Write the canonical physical control table (F5 layer 1).

    Numbers are written with full double precision (repr) so that reading the
    table back reproduces the compiled controls bit-exactly; read-back
    differences therefore measure device-side changes, not serialization.
    """
    n = len(table["time_s"])
    rows = ["# tweezer-experiment physical control table (v1)",
            f"# interpolation_between_rows: {table.get('interpolation', 'linear')}",
            f"# velocity_rule: {table.get('velocity_rule', 'provided')}",
            f"# source: {table.get('source', 'compiled-protocol')}"]
    for key, value in (metadata or {}).items():
        rows.append(f"# {key}: {value}")
    rows.append(",".join(COLUMNS))
    body = np.column_stack([
        U.s_to_us(np.asarray(table["time_s"])),
        U.m_to_um(np.asarray(table["x_m"])),
        U.m_to_um(np.asarray(table["y_m"])),
        np.asarray(table["vx_m_per_s"]) * 1.0,       # 1 m/s == 1 um/us exactly
        np.asarray(table["vy_m_per_s"]) * 1.0,
        U.joule_to_temperature_uK(np.asarray(table["depth_j"])),
        np.asarray(table["slm_factor"])])
    for row in body:
        rows.append(",".join(repr(float(v)) for v in row))
    Path(path).write_text("\n".join(rows) + "\n", encoding="utf-8")
    assert n == len(body)
