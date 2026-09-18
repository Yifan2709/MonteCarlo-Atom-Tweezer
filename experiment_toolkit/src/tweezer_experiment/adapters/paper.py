"""Optional paper adapters (A16). Never imported by the core package.

The research repository's Fig. 6d assets are explicit user-supplied
resources: this adapter loads the digitized ML trajectory nodes and the
vector-extracted reference markers from paths the caller provides. Remove
the files and the toolkit keeps working; nothing in the core reads them.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from ..errors import ToolkitError, E_MISSING_FIELD
from ..schemas import ProtocolSpec, Segment, SequenceItem
from .. import units as U


def load_ml_nodes(path: str | Path) -> dict:
    """Paper Fig. 6c ML trajectory: grid_us / pos_um / depth_mK arrays."""
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in ("grid_us", "pos_um", "depth_mK"):
        if key not in raw:
            raise ToolkitError(E_MISSING_FIELD, f"ML node file missing '{key}'", str(path))
    return {k: np.asarray(v, dtype=float) for k, v in
            {k: raw[k] for k in ("grid_us", "pos_um", "depth_mK")}.items()}


def write_ml_control_csv(nodes: dict, path: str | Path, *, aod_depth_uK: float) -> None:
    """Materialize ML nodes as a physical control CSV (cubic interpolation
    happens at compile time via the csv segment's interpolation rule)."""
    t = nodes["grid_us"]
    pos, depth_mK = nodes["pos_um"], nodes["depth_mK"]
    vx = np.gradient(pos, t)
    lines = ["# tweezer-experiment paper-ML replay table (from user-supplied nodes)",
             "# velocity_rule: provided (central differences of digitized nodes)",
             ",".join(["time_us", "aod_x_um", "aod_y_um", "aod_vx_um_per_us",
                       "aod_vy_um_per_us", "aod_depth_uK", "slm_factor"])]
    for i in range(len(t)):
        lines.append(",".join([f"{t[i]:.10g}", f"{pos[i]:.10g}", "0",
                               f"{vx[i]:.10g}", "0",
                               f"{depth_mK[i] * 1000.0:.10g}", "1"]))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def reverse_control_csv(src: str | Path, dst: str | Path) -> None:
    """Write the exact time reverse of a control CSV (x_d(t) = x_p(T-t),
    d_d(t) = d_p(T-t), v negated). Used for the paper's ML drop-off, which
    is the ML replay run backward — not an analytic family."""
    src, dst = Path(src), Path(dst)
    lines = [ln for ln in src.read_text(encoding="utf-8").splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]
    header = [h.strip() for h in lines[0].split(",")]
    rows = [[float(v) for v in ln.split(",")] for ln in lines[1:]]
    out = ["# time-reversed ML replay (drop-off): same increasing time axis,",
           "# position/depth from the reversed sequence, velocities negated"]
    out.append(",".join(header))
    n = len(rows)
    for i in range(n):
        current = dict(zip(header, rows[i]))           # time from row i
        mirror = dict(zip(header, rows[n - 1 - i]))    # values from row n-1-i
        current["aod_x_um"] = mirror["aod_x_um"]
        current["aod_y_um"] = mirror["aod_y_um"]
        current["aod_depth_uK"] = mirror["aod_depth_uK"]
        current["slm_factor"] = mirror.get("slm_factor", 1.0)
        current["aod_vx_um_per_us"] = -mirror["aod_vx_um_per_us"]
        current["aod_vy_um_per_us"] = -mirror["aod_vy_um_per_us"]
        out.append(",".join(repr(current[h]) for h in header))
    dst.write_text("\n".join(out) + "\n", encoding="utf-8")


def matched_roundtrip_spec(*, pickup_duration_us: float, wait_us: float = 100.0,
                           distance_um: float = 2.4, aod_depth_uK: float = 280.0,
                           repeats: int = 30, ml_csv: str | Path | None = None,
                           ramp_fraction: float = 0.48,
                           family: str = "paper_manual") -> ProtocolSpec:
    """The Fig. 6d short-transfer protocol as a generic ProtocolSpec.

    pickup (observe aod) -> wait with SLM off -> drop-off (observe slm);
    ``ml_csv`` replaces the analytic pickup with the digitized ML replay AND
    the drop-off with its exact time reverse (the paper replays the ML
    trajectory backward; it is not an analytic family).
    """
    pickup_seg = (Segment("pickup", "csv", None, None, None,
                          {"file": str(ml_csv), "interpolation": "cubic",
                           "velocity_columns": True}, "aod")
                  if ml_csv is not None else
                  Segment("pickup", "transfer", U.us_to_s(pickup_duration_us), None,
                          {"from_m": (0.0, 0.0), "to_m": (U.um_to_m(distance_um), 0.0),
                           "depth_from_j": 0.0, "depth_to_j": U.temperature_uK_to_joule(aod_depth_uK),
                           "family": family, "ramp_fraction": ramp_fraction,
                           "slm_factor": 1.0}, None, "aod"))
    if ml_csv is not None:
        dropoff = Segment("dropoff", "csv", None, None, None,
                          {"file": str(ml_csv) + ".reversed.csv", "interpolation": "cubic",
                           "velocity_columns": True}, "slm")
    else:
        dropoff = Segment("dropoff", "transfer", U.us_to_s(pickup_duration_us), None,
                          {"from_m": (U.um_to_m(distance_um), 0.0), "to_m": (0.0, 0.0),
                           "depth_from_j": U.temperature_uK_to_joule(aod_depth_uK),
                           "depth_to_j": 0.0,
                           "family": family, "ramp_fraction": ramp_fraction,
                           "order": "move_then_depth",
                           "slm_factor": 1.0}, None, "slm")
    wait = Segment("wait", "static", U.us_to_s(wait_us),
                   {"center_m": (U.um_to_m(distance_um), 0.0), "depth_j": U.temperature_uK_to_joule(aod_depth_uK),
                    "slm_factor": 0.0}, None, None, None)
    segments = {"pickup": pickup_seg, "wait": wait, "dropoff": dropoff}
    sequence = (SequenceItem(repeat=repeats, items=(SequenceItem(segment="pickup"),
                                                    SequenceItem(segment="wait"),
                                                    SequenceItem(segment="dropoff"))),)
    return ProtocolSpec(segments, sequence, None, None)


def load_reference_markers(path: str | Path) -> dict:
    """Vector-extracted Fig. 6d experimental markers (n, survival, errors)."""
    rows = list(csv.DictReader(Path(path).open(encoding="utf-8")))
    out: dict[str, dict] = {}
    for row in rows:
        w = row["waveform"]
        out.setdefault(w, {"n": [], "survival": [], "error_lower": [], "error_upper": []})
        for key in ("n", "survival", "error_lower", "error_upper"):
            out[w][key].append(float(row[key]))
    return {w: {k: np.asarray(v) for k, v in d.items()} for w, d in out.items()}
