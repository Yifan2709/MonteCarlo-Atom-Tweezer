"""Protocol compilation: specs -> control grids, checkpoints, audits (F2).

Compile rules (documented contract):
- Every segment must be divisible by the mechanical dt (E_NOT_DIVISIBLE
  otherwise, with the nearest valid duration reported).
- Segment grids are concatenated with **duplicated boundary times**: the left
  node belongs to the previous segment's last step and the right node starts
  the next segment, so an SLM switch changes the potential without moving
  the atom. Boundary discontinuities are collected as *events*, never
  smoothed.
- Checkpoints are exactly the segments carrying ``observe`` (plus the
  optional final hold): their count follows the protocol definition and is
  never tied to a 2*rounds+1 layout.
- Overshoot inside transfer segments and cubic CSV tables is measured on a
  dense grid and reported (disclosure, not rejection).
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from .errors import (ToolkitError, E_NOT_DIVISIBLE, E_BAD_VALUE, E_CONFLICT)
from .schemas import ProtocolSpec, Segment
from . import units as U
from .waveforms import transfer_controls
from .csv_io import read_control_csv

MAX_INSTANCES = 100_000


@dataclass(frozen=True)
class CompiledProtocol:
    controls: np.ndarray          # (N, 7): t, x, y, vx, vy, aod_depth_j, slm_factor
    segment_bounds: tuple         # ((i0, i1), ...) inclusive node indices per instance
    segment_names: tuple          # per instance, e.g. "pickup#3"
    checkpoints: tuple            # ((label, instance_idx, "aod"|"slm"), ...)
    duration_s: float
    dt_s: float
    events: tuple                 # ((t_s, field, left, right), ...) boundary jumps
    audit: dict

    @property
    def n_checkpoints(self) -> int:
        return len(self.checkpoints)


def apply_override(spec: ProtocolSpec, target: str, value) -> ProtocolSpec:
    """Return a new ProtocolSpec with one searchable field overridden.

    Targets: segment.<id>.duration_us | .ramp_fraction | .depth_to_uK | .path.
    Used by constrained search; device/preparation are never touched here.
    """
    _, seg_id, fieldname = target.split(".")
    if seg_id not in spec.segments:
        raise ToolkitError(E_BAD_VALUE, f"unknown segment '{seg_id}'", target)
    seg = spec.segments[seg_id]
    if fieldname == "duration_us":
        if seg.kind not in ("static", "transfer"):
            raise ToolkitError(E_BAD_VALUE, "duration_us applies to static/transfer segments", target)
        seg = replace(seg, duration_s=U.us_to_s(float(value)))
    elif seg.kind != "transfer":
        raise ToolkitError(E_BAD_VALUE, f"{fieldname} applies to transfer segments", target)
    if fieldname == "ramp_fraction":
        tr = dict(seg.transfer)
        tr["ramp_fraction"] = float(value)
        seg = replace(seg, transfer=tr)
    elif fieldname == "depth_to_uK":
        tr = dict(seg.transfer)
        tr["depth_to_j"] = U.temperature_uK_to_joule(float(value))
        seg = replace(seg, transfer=tr)
    elif fieldname == "path":
        tr = dict(seg.transfer)
        tr["family"] = str(value)
        seg = replace(seg, transfer=tr)
    segments = dict(spec.segments)
    segments[seg_id] = seg
    return replace(spec, segments=segments)


def _expand(spec: ProtocolSpec) -> list[tuple[str, Segment]]:
    out: list[tuple[str, Segment]] = []

    def walk(items, counts: dict[str, int]):
        for item in items:
            if item.segment:
                counts[item.segment] = counts.get(item.segment, 0) + 1
                n = counts[item.segment]
                label = item.segment if n == 1 else f"{item.segment}#{n}"
                out.append((label, spec.segments[item.segment]))
            else:
                for _ in range(item.repeat):
                    walk(item.items, counts)
        return

    walk(spec.sequence, {})
    if len(out) > MAX_INSTANCES:
        raise ToolkitError(E_BAD_VALUE, f"protocol expands to {len(out)} segments "
                                        f"(max {MAX_INSTANCES})", "protocol.sequence")
    return out


def _grid(duration_s: float, dt_s: float, label: str) -> np.ndarray:
    steps = int(round(duration_s / dt_s))
    if steps < 1 or not np.isclose(steps * dt_s, duration_s, rtol=0, atol=1e-15):
        nearest = steps * dt_s if steps >= 1 else dt_s
        raise ToolkitError(E_NOT_DIVISIBLE,
                           f"segment '{label}': duration {duration_s * 1e6:.6f} us is not a "
                           f"multiple of dt {dt_s * 1e6:.6f} us (nearest {nearest * 1e6:.6f} us)",
                           f"segment.{label}.duration_us")
    return np.linspace(0.0, duration_s, steps + 1)


def compile_protocol(bundle, protocol: ProtocolSpec | None = None) -> CompiledProtocol:
    spec = protocol if protocol is not None else bundle.protocol
    dt_s = bundle.numerics.dt_s
    instances = _expand(spec)
    controls_blocks: list[np.ndarray] = []
    names: list[str] = []
    bounds: list[tuple[int, int]] = []
    checkpoints: list[tuple[str, int, str]] = []
    events: list[tuple] = []
    audit: dict = {"overshoot": [], "velocity_jumps": [], "accel_notes": []}
    offset = 0
    start_time = 0.0

    for idx, (label, seg) in enumerate(instances):
        if seg.kind == "static":
            local = _grid(seg.duration_s, dt_s, label) if seg.duration_s is not None else None
            st = seg.static
            block = np.column_stack([
                local + start_time,
                np.full_like(local, st["center_m"][0]),
                np.full_like(local, st["center_m"][1]),
                np.zeros_like(local), np.zeros_like(local),
                np.full_like(local, st["depth_j"]),
                np.full_like(local, st["slm_factor"])])
            duration = seg.duration_s
        elif seg.kind == "transfer":
            tr = seg.transfer
            block_local = transfer_controls(tr["family"], tr["ramp_fraction"], seg.duration_s,
                                            tr["from_m"], tr["to_m"], tr["depth_from_j"],
                                            tr["depth_to_j"], tr.get("order", "depth_then_move"))
            local = _grid(seg.duration_s, dt_s, label)
            x, y, vx, vy, depth = block_local(local)
            block = np.column_stack([local + start_time, x, y, vx, vy, depth,
                                     np.full_like(local, tr["slm_factor"])])
            duration = seg.duration_s
            _audit_transfer(audit, label, tr, local, x, y, depth)
        else:  # csv
            table = read_control_csv(seg.csv["file"], velocity_columns=seg.csv["velocity_columns"],
                                     interpolation=seg.csv["interpolation"])
            duration = float(table["time_s"][-1] - table["time_s"][0])
            if duration <= 0:
                raise ToolkitError(E_BAD_VALUE, f"csv segment '{label}' covers zero time",
                                   f"segment.{label}.file")
            local = _grid(duration, dt_s, label)
            t_abs = table["time_s"] - table["time_s"][0]
            cols = [_interp(t_abs, table[key], local, seg.csv["interpolation"], key)
                    for key in ("x_m", "y_m", "vx_m_per_s", "vy_m_per_s", "depth_j", "slm_factor")]
            block = np.column_stack([local + start_time, *cols])
            if seg.csv["interpolation"] == "cubic":
                dense = np.linspace(0, duration, 2001)
                _audit_overshoot(audit, label, "aod_x_um",
                                 U.m_to_um(_interp(t_abs, table["x_m"], dense, "cubic", "x_m")),
                                 U.m_to_um(table["x_m"]))
                _audit_overshoot(audit, label, "aod_depth_uK",
                                 U.joule_to_temperature_uK(_interp(t_abs, table["depth_j"], dense, "cubic", "depth_j")),
                                 U.joule_to_temperature_uK(table["depth_j"]))
        controls_blocks.append(block)
        names.append(label)
        bounds.append((offset, offset + len(block) - 1))
        offset += len(block)
        start_time += duration
        if seg.observe is not None:
            checkpoints.append((f"after_{label}", idx, seg.observe))

    # optional final hold: static SLM hold with the AOD parked wherever the
    # protocol left it, depth zero (off), observe per spec.
    if spec.final_hold_s is not None:
        last = controls_blocks[-1]
        park_x, park_y = last[-1, 1], last[-1, 2]
        local = _grid(spec.final_hold_s, dt_s, "final_hold")
        block = np.column_stack([local + start_time,
                                 np.full_like(local, park_x), np.full_like(local, park_y),
                                 np.zeros_like(local), np.zeros_like(local),
                                 np.zeros_like(local), np.ones_like(local)])
        controls_blocks.append(block)
        names.append("final_hold")
        bounds.append((offset, offset + len(block) - 1))
        offset += len(block)
        start_time += spec.final_hold_s
        checkpoints.append(("after_final_hold", len(names) - 1, spec.final_observe or "slm"))

    controls = np.vstack(controls_blocks)

    # boundary events: compare right-start of each instance with left-end of
    # the previous one (positions, depth, slm factor, commanded velocities;
    # velocities are m/s == um/us)
    for j in range(1, len(bounds)):
        left = controls[bounds[j - 1][1]]
        right = controls[bounds[j][0]]
        t = right[0]
        pairs = {"aod_x_um": (left[1], right[1]), "aod_y_um": (left[2], right[2]),
                 "aod_vx_um_per_us": (left[3], right[3]), "aod_vy_um_per_us": (left[4], right[4]),
                 "aod_depth_uK": (left[5], right[5]), "slm_factor": (left[6], right[6])}
        for fieldname, (lv, rv) in pairs.items():
            if abs(rv - lv) > 1e-9 * max(1.0, abs(lv), abs(rv)):
                events.append({"t_s": float(t), "field": fieldname,
                               "left": float(lv), "right": float(rv),
                               "between": (names[j - 1], names[j])})
    return CompiledProtocol(controls, tuple(bounds), tuple(names), tuple(checkpoints),
                            float(start_time), dt_s, tuple(events), audit)


def _interp(t_src, v_src, t_new, rule: str, name: str) -> np.ndarray:
    if rule == "linear":
        return np.interp(t_new, t_src, v_src)
    if rule == "hold":
        idx = np.searchsorted(t_src, t_new, side="right") - 1
        return v_src[np.clip(idx, 0, len(t_src) - 1)]
    # cubic: monotone-free natural cubic (PchipFactory would smooth extrema);
    # use numpy-free scipy? No scipy dependency: fall back to local cubic via
    # np.polyfit on nearest 4 points would be slow; instead require linear for
    # phase 1 unless scipy present.
    try:
        from scipy.interpolate import CubicSpline
    except ImportError as exc:
        raise ToolkitError(E_BAD_VALUE,
                           "interpolation 'cubic' requires scipy; use linear or hold", "interpolation") from exc
    return CubicSpline(t_src, v_src)(t_new)


def _audit_transfer(audit: dict, label: str, tr: dict, local, x, y, depth) -> None:
    for name, arr, ends in (("aod_x_um", x, (tr["from_m"][0], tr["to_m"][0])),
                            ("aod_y_um", y, (tr["from_m"][1], tr["to_m"][1])),
                            ("aod_depth_uK", depth, (tr["depth_from_j"], tr["depth_to_j"]))):
        lo, hi = min(ends), max(ends)
        scale = max(hi - lo, 1e-30)
        over = max(arr.max() - hi, lo - arr.min()) / scale
        if over > 1e-9:
            audit["overshoot"].append({"segment": label, "field": name,
                                       "fraction_of_span": float(over)})
    if tr["family"] == "paper_manual":
        audit["accel_notes"].append(
            {"segment": label, "note": "paper_manual has acceleration jumps at move endpoints"})


def _audit_overshoot(audit: dict, label: str, fieldname: str, dense, nodes) -> None:
    lo, hi = nodes.min(), nodes.max()
    scale = max(hi - lo, 1e-30)
    over = max(dense.max() - hi, lo - dense.min()) / scale
    if over > 1e-9:
        audit["overshoot"].append({"segment": label, "field": fieldname,
                                   "fraction_of_span": float(over),
                                   "between": "cubic-interior beyond node extremes"})
