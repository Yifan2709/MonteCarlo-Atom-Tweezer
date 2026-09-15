"""Two explicit protocols; Fig. 6d comparison never silently adds long moves.

Segment boundaries have separate left/right controls, so an SLM switch changes
the potential, not the atom's position or velocity. Finite pulse edges are exact
grid boundaries; mechanical motion continues during each pulse.
"""
from dataclasses import dataclass
import json
import numpy as np
from level2_joint_transfer.level2_simulation import build_waveform, physics_from_config
from level2_joint_transfer.waveforms import MLCubicWaveform
from level2c_pickup_survival.preflight2c import pickup_spec
from level3_3d_transport_lensing.transport_profiles import get_profile
from level5_finite_pulse_irb.dd_finite_pulses import dd_pulse_events
from level5_finite_pulse_irb.quantum_propagators import su2_propagator


@dataclass
class ContinuousProtocol:
    controls: np.ndarray  # t,x,y,vx,vy,AOD depth,SLM factor,Omega,phi
    segments: np.ndarray  # first,last,jitter stream,judge (0/1=AOD/2=SLM/3=extra AOD)
    names: list
    pulses: list
    duration_s: float
    distance_per_round_m: float
    arm: str
    dd_mode: str
    checkpoint_times_s: tuple


def waveform_for(key, cfg, ml_path):
    physics = physics_from_config(cfg.base)
    if key.startswith("manual_"):
        duration = float(key.removeprefix("manual_").removesuffix("us"))
        return build_waveform(pickup_spec(duration, cfg), physics)
    if key != "ml_400us":
        raise ValueError(f"Unknown waveform {key}")
    raw = json.loads(ml_path.read_text(encoding="utf-8"))
    t = np.asarray(raw["grid_us"]) * 1e-6
    # Preserve the actual stage-1 replay: raw digitized offsets/end depths and
    # ordinary CubicSpline, not endpoint normalization or monotone interpolation.
    return MLCubicWaveform(t[-1], t / t[-1], np.asarray(raw["pos_um"]) * 1e-6,
                          np.asarray(raw["depth_mK"]) * 1.380649e-26,
                          physics["aod_initial_center_m"], physics["aod_initial_depth_j"])


def build_protocol(waveform, wait_s, dt_s, *, arm="matched", long_distance_m=375e-6,
                   long_duration_s=1450e-6, remote_hold_s=54e-6,
                   omega=2 * np.pi * 24611, dd=True):
    if arm not in ("matched", "extended") or dt_s <= 0 or wait_s <= 0:
        raise ValueError("Invalid protocol arm or time step")
    T = waveform.duration_s
    end = float(waveform.center(T))
    depth = float(waveform.depth(T))
    profile = get_profile("adiabatic_sine")
    def short(t, reverse=False):
        s = T - t if reverse else t
        x = np.asarray(waveform.center(s))
        v = np.asarray(waveform.center_velocity(s)) * (-1 if reverse else 1)
        return x, np.zeros_like(t), v, np.zeros_like(t), waveform.depth(s)
    def hold(t, far=False):
        offset = long_distance_m / np.sqrt(2) if far else 0.0
        return (np.full_like(t, end + offset), np.full_like(t, offset),
                np.zeros_like(t), np.zeros_like(t), np.full_like(t, depth))
    def long(t, reverse=False):
        s = t / long_duration_s
        q = profile.q(1 - s if reverse else s)
        v = profile.q_dot(1 - s if reverse else s) * (-1 if reverse else 1)
        offset = long_distance_m / np.sqrt(2) * q
        v = long_distance_m / np.sqrt(2) / long_duration_s * v
        return end + offset, offset, v, v, np.full_like(t, depth)
    specs = [("pickup", T, lambda t: short(t), 1, 0, 0),
             ("pickup_wait", wait_s, hold, 0, 1, 1)]
    windows = []
    if arm == "extended":
        start = T + wait_s
        windows = [(start, start + long_duration_s),
                   (start + long_duration_s + remote_hold_s,
                    start + 2 * long_duration_s + remote_hold_s)]
        specs += [("outbound", long_duration_s, long, 0, 3, 3),
                  ("remote_hold", remote_hold_s, lambda t: hold(t, True), 0, 4, 3),
                  ("inbound", long_duration_s, lambda t: long(t, True), 0, 5, 3)]
    specs += [("dropoff", T, lambda t: short(t, True), 1, 2, 2)]
    duration = sum(s[1] for s in specs)
    mode = ("xy4" if arm == "matched" else "xy4_per_long_move") if dd else "none"
    pulses = dd_pulse_events(mode, 0, duration, omega, long_move_windows=windows)
    edges = np.array([x for p in pulses for x in (p.start, p.start + p.duration)])
    all_controls, segments, names = [], [], []
    start, offset = 0.0, 0
    for name, length, fn, slm, stream, judge in specs:
        steps = int(round(length / dt_s))
        if steps < 1 or not np.isclose(steps * dt_s, length, rtol=0, atol=1e-14):
            raise ValueError(f"Segment {name} must be divisible by dt")
        local = np.linspace(0, length, steps + 1)
        local = np.unique(np.r_[local, edges[(edges > start) & (edges < start + length)] - start])
        x, y, vx, vy, d = fn(local)
        t = local + start
        mid = (t[:-1] + t[1:]) / 2
        om, ph = np.zeros_like(t), np.zeros_like(t)
        for p in pulses:
            active = (mid >= p.start) & (mid < p.start + p.duration)
            om[:-1][active] = p.nominal_area / p.duration
            ph[:-1][active] = p.axis_phase
        all_controls.append(np.column_stack((t, x, y, vx, vy, d,
                                             np.full_like(t, slm), om, ph)))
        segments.append((offset, offset + len(t) - 1, stream, judge))
        names.append(name)
        offset += len(t)
        start += length
    return ContinuousProtocol(np.vstack(all_controls), np.array(segments, dtype=np.int64),
                              names, pulses, duration,
                              2 * abs(end - float(waveform.center(0))) +
                              (2 * long_distance_m if arm == "extended" else 0),
                              arm, mode, (T + wait_s, duration))


def ideal_checkpoint_unitaries(protocol, rounds):
    """Ideal finite-drive target at every observation, including partial XY4 blocks."""
    def until(t):
        u = np.eye(2, dtype=complex)
        for p in protocol.pulses:
            dt = min(max(t - p.start, 0.0), p.duration)
            if dt:
                elements = su2_propagator(0.0, p.nominal_area / p.duration, p.axis_phase, dt)
                u = np.asarray(elements).reshape(2, 2) @ u
        return u
    partial, full = [until(t) for t in protocol.checkpoint_times_s]
    out, previous = [np.eye(2, dtype=complex)], np.eye(2, dtype=complex)
    for _ in range(rounds):
        out.extend((partial @ previous, full @ previous))
        previous = full @ previous
    return np.asarray(out)
