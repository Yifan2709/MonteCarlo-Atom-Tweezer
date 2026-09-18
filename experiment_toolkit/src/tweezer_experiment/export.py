"""Waveform export, hardware mapping and read-back verification (F5).

Two layers:
1. Physical control table (always): time, AOD x/y + velocities, AOD depth,
   SLM factor, with units and provenance in the header comments. Boundary
   time duplicates and switch events are preserved, never smoothed.
2. Calibrated neutral instruction table (only with a HardwareCalibration):
   resample at the hardware clock, map position -> RF frequency and depth ->
   amplitude through monotonic calibration tables, quantize, and check
   range/slew limits. Violations are reported and no compliant-marked
   instruction file is produced.

Read-back: the exported artifacts are re-imported, the device-side control
(quantized, clock-limited) is reconstructed and the simulation is rerun on
the SAME initial pool; the report quantifies what the export chain changed.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from .errors import ToolkitError, E_CAL_RANGE, E_MISSING_FIELD, E_BAD_VALUE
from . import units as U
from .csv_io import write_control_csv, read_control_csv
from .schemas import ExperimentBundle, HardwareCalibration
from .simulation import simulate


def _table_from_compiled(compiled) -> dict:
    c = compiled.controls
    return {"time_s": c[:, 0], "x_m": c[:, 1], "y_m": c[:, 2],
            "vx_m_per_s": c[:, 3], "vy_m_per_s": c[:, 4],
            "depth_j": c[:, 5], "slm_factor": c[:, 6],
            "interpolation": "linear", "velocity_rule": "provided",
            "source": "compiled-protocol"}


def _map_forward(channel, values):
    """Monotone-table mapping with explicit extrapolation rejection."""
    lo, hi = channel.valid_input
    if np.any(values < lo - 1e-9) or np.any(values > hi + 1e-9):
        bad = int(np.argmax((values < lo - 1e-9) | (values > hi + 1e-9)))
        raise ToolkitError(E_CAL_RANGE,
                           f"input {values[bad]:.6f} {channel.input_unit} outside calibrated range "
                           f"[{lo}, {hi}] {channel.input_unit} (first offending sample #{bad})",
                           "calibration")
    return np.interp(values, channel.input_axis, channel.output_si)


def _map_inverse(channel, values_si):
    return np.interp(values_si, channel.output_si, channel.input_axis)


def _quantize(values, lo, hi, bits):
    if bits is None:
        return values, 0.0
    levels = 2 ** bits - 1
    q = np.round((np.asarray(values) - lo) / (hi - lo) * levels) / levels * (hi - lo) + lo
    return q, float(np.max(np.abs(q - values)))


def export_controls(bundle: ExperimentBundle, out_dir: str | Path, *,
                    protocol=None, calibration: HardwareCalibration | None = None,
                    verify: bool = True, meta: dict | None = None) -> dict:
    from .protocol_compile import compile_protocol
    from .provenance import write_json
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    compiled = compile_protocol(bundle, protocol)
    table = _table_from_compiled(compiled)
    physical_path = out / "physical_controls.csv"
    write_control_csv(physical_path, table, meta or {"units": "SI in header names"})

    manifest = {
        "physical_table": {"path": str(physical_path),
                           "rows": int(len(table["time_s"])),
                           "duplicate_boundary_times": "preserved (segment left/right limits)",
                           "switch_events": [dict(e) for e in compiled.events],
                           "audit": compiled.audit},
        "device_instructions": None,
        "readback": None,
    }

    cal = calibration if calibration is not None else bundle.calibration
    if cal is not None:
        clock_s = 1.0 / cal.clock_hz
        duration = float(table["time_s"][-1])
        t_ck = np.arange(0.0, duration + 0.5 * clock_s, clock_s)
        t_ck[-1] = min(t_ck[-1], duration)
        x_um = U.m_to_um(np.interp(t_ck, table["time_s"], table["x_m"]))
        y_um = U.m_to_um(np.interp(t_ck, table["time_s"], table["y_m"]))
        depth_uK = U.joule_to_temperature_uK(np.interp(t_ck, table["time_s"], table["depth_j"]))
        slm = np.interp(t_ck, table["time_s"], table["slm_factor"])

        moves_y = bool(np.ptp(y_um) > 1e-9)
        chans = cal.channels
        if "aod_x" not in chans or "aod_depth" not in chans:
            raise ToolkitError(E_MISSING_FIELD,
                               "calibration needs channels 'aod_x' and 'aod_depth'", "calibration.channels")
        if moves_y and "aod_y" not in chans:
            raise ToolkitError(E_MISSING_FIELD,
                               "protocol moves y but calibration has no aod_y channel",
                               "calibration.channels.aod_y")
        rf_x = _map_forward(chans["aod_x"], x_um)
        amp = _map_forward(chans["aod_depth"], depth_uK)
        rf_y = _map_forward(chans["aod_y"], y_um) if moves_y else None

        checks = []
        limit = cal.limits
        if "max_rf_hz" in limit:
            checks.append({"check": "max_rf_hz", "value": float(np.max(rf_x)),
                           "limit": limit["max_rf_hz"],
                           "pass": bool(np.max(rf_x) <= limit["max_rf_hz"])})
            if moves_y:
                checks.append({"check": "max_rf_hz (y)", "value": float(np.max(rf_y)),
                               "limit": limit["max_rf_hz"],
                               "pass": bool(np.max(rf_y) <= limit["max_rf_hz"])})
        if "max_rf_step_hz" in limit:
            step = float(np.max(np.abs(np.diff(rf_x))))
            checks.append({"check": "max_rf_step_hz", "value": step,
                           "limit": limit["max_rf_step_hz"], "pass": bool(step <= limit["max_rf_step_hz"])})
        if "max_amplitude" in limit:
            checks.append({"check": "max_amplitude", "value": float(np.max(amp)),
                           "limit": limit["max_amplitude"], "pass": bool(np.max(amp) <= limit["max_amplitude"])})

        q_x, err_x = _quantize(rf_x, chans["aod_x"].output_si[0], chans["aod_x"].output_si[-1], cal.quantization_bits)
        q_a, err_a = _quantize(amp, chans["aod_depth"].output_si[0], chans["aod_depth"].output_si[-1], cal.quantization_bits)
        q_y, err_y = ((None, 0.0) if rf_y is None else
                      _quantize(rf_y, chans["aod_y"].output_si[0], chans["aod_y"].output_si[-1], cal.quantization_bits))

        all_pass = all(c["pass"] for c in checks)
        if all_pass:
            cols = ["time_s", "rf_x_hz", "amplitude", "slm_digital"]
            rows = [",".join(cols)]
            for i in range(len(t_ck)):
                row = [f"{t_ck[i]:.9g}", f"{q_x[i]:.9g}", f"{q_a[i]:.9g}",
                       "1" if slm[i] >= 0.5 else "0"]
                if q_y is not None:
                    row.insert(2, f"{q_y[i]:.9g}")
                    if i == 0:
                        rows[0] = ",".join(["time_s", "rf_x_hz", "rf_y_hz", "amplitude", "slm_digital"])
                rows.append(",".join(row))
            rf_path = out / "rf_instructions.csv"
            rf_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
            manifest["device_instructions"] = {
                "path": str(rf_path), "clock_hz": cal.clock_hz,
                "quantization_bits": cal.quantization_bits,
                "calibration_version": cal.version, "calibration_source": cal.source_tag,
                "format": "neutral RF table (synthetic calibration example; not a vendor protocol)",
                "quantization_max_abs_error": {"rf_x_hz": err_x, "amplitude": err_a,
                                               "rf_y_hz": err_y if q_y is not None else None},
                "limit_checks": checks, "limits_ok": True}
        else:
            manifest["device_instructions"] = {
                "produced": False,
                "reason": "hardware limit checks failed; no compliant instruction file written",
                "limit_checks": checks, "limits_ok": False}

    if verify:
        manifest["readback"] = _readback_verify(bundle, compiled, out, cal)
    from .provenance import write_json as _wj
    _wj(out / "export_manifest.json", manifest)
    return manifest


def _reconstructed_from_rf(rf_path: Path, cal: HardwareCalibration, compiled, dt_s: float):
    """Device-side control: invert calibrated/quantized RF back to physical."""
    data = np.genfromtxt(rf_path, delimiter=",", names=True)
    t_ck = np.asarray(data["time_s"], dtype=float)
    x_um = _map_inverse(cal.channels["aod_x"], np.asarray(data["rf_x_hz"], dtype=float))
    depth_uK = _map_inverse(cal.channels["aod_depth"], np.asarray(data["amplitude"], dtype=float))
    names = data.dtype.names
    y_um = (_map_inverse(cal.channels["aod_y"], np.asarray(data["rf_y_hz"], dtype=float))
            if "rf_y_hz" in names else np.zeros_like(x_um))
    slm = np.asarray(data["slm_digital"], dtype=float)
    t_grid = compiled.controls[:, 0]
    x_m = U.um_to_m(np.interp(t_grid, t_ck, x_um))
    y_m = U.um_to_m(np.interp(t_grid, t_ck, y_um))
    depth_j = U.temperature_uK_to_joule(np.interp(t_grid, t_ck, depth_uK))
    factor = np.interp(t_grid, t_ck, slm)
    vx = np.gradient(x_m, t_grid)
    vy = np.gradient(y_m, t_grid)
    new_controls = np.column_stack([t_grid, x_m, y_m, vx, vy, depth_j, factor])
    return replace(compiled, controls=new_controls)


def _dedupe_boundary_times(t, *cols):
    """Toolkit tables intentionally repeat boundary times (segment left/right
    limits). Collapse duplicates keeping the RIGHT-limit row (the incoming
    segment's control governs from that instant)."""
    t = np.asarray(t, dtype=float)
    keep = np.concatenate([[True], np.diff(t) > 0])
    return (t[keep],) + tuple(np.asarray(c)[keep] for c in cols)


def _reconstructed_from_csv(path: Path, compiled):
    raw_table = _read_table_with_boundary_rule(path)
    t_grid = compiled.controls[:, 0]
    if len(raw_table["time_s"]) == len(t_grid) and np.allclose(
            raw_table["time_s"], t_grid, rtol=0, atol=1e-15):
        # Same row count and identical times (our own full-precision export):
        # rebuild row-by-row so duplicated boundary nodes keep BOTH the left
        # and right limit values — bit-exact round trip.
        new_controls = np.column_stack([
            t_grid, raw_table["x_m"], raw_table["y_m"], raw_table["vx_m_per_s"],
            raw_table["vy_m_per_s"], raw_table["depth_j"], raw_table["slm_factor"]])
        return replace(compiled, controls=new_controls)
    t_ck, x_m, y_m, vx, vy, depth_j, factor = _dedupe_boundary_times(
        raw_table["time_s"], raw_table["x_m"], raw_table["y_m"], raw_table["vx_m_per_s"],
        raw_table["vy_m_per_s"], raw_table["depth_j"], raw_table["slm_factor"])
    table = {"time_s": t_ck, "x_m": x_m, "y_m": y_m, "vx_m_per_s": vx,
             "vy_m_per_s": vy, "depth_j": depth_j, "slm_factor": factor}
    t_grid = compiled.controls[:, 0]
    cols = [np.interp(t_grid, table["time_s"], table[key])
            for key in ("x_m", "y_m", "vx_m_per_s", "vy_m_per_s", "depth_j", "slm_factor")]
    new_controls = np.column_stack([t_grid, *cols])
    return replace(compiled, controls=new_controls)


def _read_table_with_boundary_rule(path: Path) -> dict:
    """Read a toolkit-exported physical table, accepting its documented
    duplicate boundary times (unlike the strict user-import reader)."""
    import csv as _csv
    lines = [ln for ln in Path(path).read_text(encoding="utf-8").splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]
    header = [h.strip() for h in lines[0].split(",")]
    rows = np.asarray([[float(v) for v in ln.split(",")] for ln in lines[1:]], dtype=float)
    col = {name: rows[:, i] for i, name in enumerate(header)}
    return {"time_s": col["time_us"] * 1e-6, "x_m": col["aod_x_um"] * 1e-6,
            "y_m": col["aod_y_um"] * 1e-6,
            "vx_m_per_s": col["aod_vx_um_per_us"] * 1.0,   # um/us == m/s
            "vy_m_per_s": col["aod_vy_um_per_us"] * 1.0,
            "depth_j": col["aod_depth_uK"] * 1.380649e-29,
            "slm_factor": col["slm_factor"]}


def _readback_verify(bundle: ExperimentBundle, compiled, out: Path, cal) -> dict:
    base = simulate(bundle, compiled=compiled)
    report = {"reference": {"final_capture": float(base.survival[-1]),
                            "curve": base.survival.tolist(),
                            "seed": bundle.numerics.seed, "shots": base.shots}}
    # layer 1: CSV round trip at full precision
    try:
        rb1 = _reconstructed_from_csv(out / "physical_controls.csv", compiled)
        r1 = simulate(bundle, compiled=rb1)
        scale = np.maximum(np.abs(compiled.controls[:, 1:]), 1e-30)
        control_rel_change = float(np.max(np.abs(rb1.controls[:, 1:] - compiled.controls[:, 1:])
                                          / scale))
        report["physical_csv_roundtrip"] = {
            "controls_max_relative_change": control_rel_change,
            "controls_exact": bool(control_rel_change < 1e-12),
            "max_abs_dS": float(np.max(np.abs(r1.survival - base.survival))),
            "final_capture": float(r1.survival[-1]),
            "note": "unit-conversion round trips can differ at 1 ULP; with stochastic "
                    "noise such residuals diverge chaotically, so dS also contains MC-level "
                    "chaos amplification — control fidelity is the exactness criterion"}
    except (ToolkitError, ValueError, OSError) as exc:
        report["physical_csv_roundtrip"] = {"error": str(exc)}
    # layer 2: device-side reconstruction from quantized clocked RF
    rf_path = out / "rf_instructions.csv"
    if cal is not None and rf_path.is_file():
        try:
            rb2 = _reconstructed_from_rf(rf_path, cal, compiled, compiled.dt_s)
            r2 = simulate(bundle, compiled=rb2)
            common = min(len(r2.survival), len(base.survival))
            report["rf_clock_quantized_roundtrip"] = {
                "max_abs_dS": float(np.max(np.abs(r2.survival[:common] - base.survival[:common]))),
                "final_capture": float(r2.survival[-1]),
                "note": "zero-order-hold RF steps resampled onto the dt grid by linear "
                        "interpolation; judges stay at the original segment boundaries"}
        except (ToolkitError, ValueError, OSError) as exc:
            report["rf_clock_quantized_roundtrip"] = {"error": str(exc)}
    return report
